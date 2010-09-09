#include <algorithm>
#include <iostream>
#include <vector>

// ROOT includes
#include "TAxis.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "TGraphErrors.h"

// CMS includes
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagHistograms.h"
#include "RecoBTag/PerformanceMeasurements/interface/PMHistograms.h"
#include "RecoBTag/PerformanceMeasurements/interface/TH1Store.h"

using namespace std;

bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event)
{
    // if the jsonVec is empty, then no JSON file was provided so all
    // events should pass
    if (jsonVec.empty())
        return true;

    bool (*funcPtr)(edm::LuminosityBlockRange const &,
                    edm::LuminosityBlockID const &) = &edm::contains;
                    edm::LuminosityBlockID lumiID (event.id().run(), 
                    event.id().luminosityBlock());

    std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
        std::find_if(jsonVec.begin(), jsonVec.end(),
                     boost::bind(funcPtr, _1, lumiID) );

    return jsonVec.end() != iter;
}

int main (int argc, char* argv[])
{
    ////////////////////////////////////////////////
    // // Command Line Options // //
    ////////////////////////////////////////////////
    optutl::CommandLineParser parser ("Plots Jet Pt");
    parser.stringValue ("outputFile") = "PM_results"; // .root added automatically

    parser.addOption ("Tagger",   optutl::CommandLineParser::kString,
                      "b tagger alias (TCHE,TCHP,JP,SSV, SSVHP, SSVHE)",
                      "TCHE");
    parser.addOption ("TaggerCut",   optutl::CommandLineParser::kDouble,
                      "b tag discriminator cut",
                      3.3);
    parser.addOption ("AwayTagger",   optutl::CommandLineParser::kString,
                      "b tagger alias (TCHE,TCHP,JP,SSV, SSVHP, SSVHE)",
                      "TCHP");
    parser.addOption ("AwayTaggerCut",   optutl::CommandLineParser::kDouble,
                      "b tag discriminator cut",
                      1.19);
    parser.addOption ("JetCollection",   optutl::CommandLineParser::kString,
                      "jet collection to use selectedPatJets,selectedPatJetsAK5PF, selectedPatJetsAK5Track)",
                      "selectedPatJets");
    parser.addOption ("TriggerName",   optutl::CommandLineParser::kString,
                      "Trigger to be used HLT_Mu5,HLT_BTagMu_Jet10U)",
                      "HLT_BTagMu_Jet10U");

    parser.addOption ("config", optutl::CommandLineParser::kString,
                      "python configuration file with JSON.",
                      "CONFIG");

    parser.addOption("hlt_path", optutl::CommandLineParser::kString,
                     "HLT path: HLT, REDIGI, PAT, etc.",
                     "HLT");

    // Parse the command line arguments
    parser.parseArguments (argc, argv);

    // Python Config file
    std::vector<edm::LuminosityBlockRange> lumiVector;
    if (parser.stringValue("config").size() &&
        "CONFIG" != parser.stringValue("config"))
    {
        cout << "Config: " << parser.stringValue("config") << endl;

        PythonProcessDesc builder(parser.stringValue("config"));

        const edm::ParameterSet &inputs =
             builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("inputs");

        if ( inputs.exists("lumisToProcess") )
        {
            const std::vector<edm::LuminosityBlockRange> &jsonVector =
                inputs.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess");

            lumiVector.resize(jsonVector.size());
            std::copy(jsonVector.begin(), jsonVector.end(), lumiVector.begin());
        }
    }
    

    //___________________
    // FWLite setup
    AutoLibraryLoader::enable();

    fwlite::ChainEvent events( parser.stringVector ("inputFiles") );

    // for now
    std::map< std::string, std::string > btagMap;
    btagMap["TCHE"] = "trackCountingHighEffBJetTags";
    btagMap["TCHP"] = "trackCountingHighPurBJetTags";
    btagMap["JP"] = "jetProbabilityBJetTags";
    btagMap["SSV"] = "simpleSecondaryVertexBJetTags";
    btagMap["SSVHP"] = "simpleSecondaryVertexHighPurBJetTags";
    btagMap["SSVHE"] = "simpleSecondaryVertexHighEffBJetTags";
    
    std::string Tagger =  btagMap.find( parser.stringValue( "Tagger" ) )->second;
    std::string awayTagger = btagMap.find(parser.stringValue( "AwayTagger" ) )->second; 
    double btag_cut_Tagger = parser.doubleValue( "TaggerCut" );
    double btag_cut_awayTagger = parser.doubleValue( "AwayTaggerCut" );
    std::string jetCollection = parser.stringValue( "JetCollection" ) ; 
    std::string triggername = parser.stringValue( "TriggerName" ) ; 

    int outputEvery = parser.integerValue ( "outputEvery" );
    int maxNevents = parser.integerValue ( "maxevents" );
    
    cout << "outputEvery = " << outputEvery << endl;

    //__________________
    // Histogram manager
    TH1Store *hstore = new TH1Store();
    PMHistograms histos( hstore );
    histos.Add();

    string hlt_path = parser.stringValue("hlt_path");
    if (!hlt_path.size())
        cout << "[warning] HLT Trigger won't be used." << endl;

    // Loop over events
    int nentry = 0;
    for(events.toBegin(); ! events.atEnd(); ++events)
    {
        if (!jsonContainsEvent (lumiVector, events) )
            continue;

        // Event is allowed. Process it normally
        if (maxNevents && nentry > maxNevents) break;

        ++nentry;

        if ( outputEvery!=0 && (nentry % outputEvery) == 0 ) 
        {
            cout << "Processing Event: " << nentry << endl;
        }
        if ( outputEvery == 0 )
            cout << "Processing Event: " << nentry << endl;

        // Process Trigger only in case HLT Path is given
        //
        if (hlt_path.size())
        {
            fwlite::Handle<edm::TriggerResults> triggerResults ;
            triggerResults.getByLabel(events, "TriggerResults","", hlt_path.c_str()) ;

            // Test if TriggerResults were successfully extracted
            //
            if (!triggerResults.isValid())
            {
                cerr << "[warning] Failed to extract trigger results." << endl;

                // Skip event
                //
                continue;
            }

            // Find Trigger with given Name
            //
            typedef vector<string> Triggers;

            const Triggers &triggerNames =
                events.triggerNames(*triggerResults).triggerNames();

            Triggers::const_iterator trigger = find(triggerNames.begin(),
                                                    triggerNames.end(),
                                                    triggername);

            if (triggerNames.end() == trigger)
                // Trigger failed
                continue;
        }

        // load object collections
        fwlite::Handle< vector< pat::Jet > > jetHandle;

        if (jetCollection == "selectedPatJetsAK5Track")
            jetHandle.getByLabel ( events, "selectedPatJetsAK5Track");
        else if (jetCollection == "selectedPatJetsAK5PF")
            jetHandle.getByLabel ( events, "selectedPatJetsAK5PF");
        else
            jetHandle.getByLabel ( events, "selectedPatJets");

        if (!jetHandle.isValid())
        {
            cerr << "Failed to extract jets." << endl;

            break;
        }

        fwlite::Handle< vector< pat::Muon > > muonHandle;
        muonHandle.getByLabel ( events, "selectedPatMuonsForPtRel");

        if (!muonHandle.isValid())
        {
            cerr << "Failed to extract muons." << endl;

            break;
        }

        TLorentzVector p4Jet;
        TLorentzVector p4MuJet;
        TLorentzVector p4Muon;
        TLorentzVector p4AwayJet;

        // Loop over jets
        for (vector< pat::Jet >::const_iterator jetIter = jetHandle->begin();
             jetIter != jetHandle->end(); ++jetIter)
        {
            bool TaggedJet = false;
            hstore->hist("jet_pt")->Fill (jetIter->pt()); // just for testing

            // get MC flavor of jet
            int JetFlavor = abs( jetIter->partonFlavour() );

            p4Jet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );
            bool hasLepton = false;

            //int tmptotmuon = 0;
            // loop over muons
            ////////////////////////////////
            double mu_highest_pt = 0;
            double ptrel = 0;
            double ptreltmp = 0;

            // Loop over muons
            for (vector< pat::Muon >::const_iterator muonIter = muonHandle->begin();
                 muonIter != muonHandle->end(); ++muonIter)
            {
                hstore->hist("muon_pt")->Fill (muonIter->innerTrack()->pt());

                reco::TrackRef trackMuon = muonIter->innerTrack();
                math::XYZTLorentzVector trackMuonP4(trackMuon->px(),
                                                    trackMuon->py(),
                                                    trackMuon->pz(),
                                                    sqrt(trackMuon->p() * trackMuon->p() + 0.1057*0.1057));
                
                // find a muon in a jet
                //check deltar
                double deltaR  = ROOT::Math::VectorUtil::DeltaR(jetIter->p4().Vect(),
                                                                trackMuonP4.Vect() );
                TVector3 tmpvecOrg(jetIter->p4().Vect().X(),jetIter->p4().Vect().Y(),  jetIter->p4().Vect().Z());
                TVector3 tmpvec = tmpvecOrg;

                TVector3 leptonvec(trackMuon->px(), trackMuon->py(),trackMuon->pz());
                tmpvec += leptonvec;

                ptreltmp = leptonvec.Perp(tmpvec);

                // muon in jet cuts
                if ( deltaR >= 0.4 || deltaR < 0.01 || ptreltmp <= -1.0 ) continue;

                hasLepton = true;

                if ( muonIter->pt() > mu_highest_pt )
                {
                    mu_highest_pt = muonIter->pt();
                    p4Muon.SetPtEtaPhiE(muonIter->pt(),
                                        muonIter->eta(),
                                        muonIter->phi(),
                                        muonIter->energy());
                    // recalculate pTrel
                    tmpvec = tmpvecOrg;
                    leptonvec.SetXYZ(muonIter->px(), muonIter->py(),muonIter->pz());
                    tmpvec += leptonvec;
                    ptrel = leptonvec.Perp(tmpvec);  // maximum
                }
                hstore->hist("jet_pTrel")->Fill (ptrel);
                hstore->hist("jet_deltaR")->Fill (deltaR);

                string suffix;
                switch(JetFlavor)
                {
                    case 5: suffix = "_b";
                            break;

                    case 4: suffix = "_c";
                            break;

                    case 21: // fall through
                    case 3 : // fall through
                    case 2 : // fall through
                    case 1 : // fall through
                             suffix = "_udsg";
                             break;
                }

                if (suffix.length())
                {
                    hstore->hist("jet_deltaR" + suffix)->Fill(deltaR);
                    hstore->hist("jet_deltaR" + suffix)->Fill(ptrel);
                }
            }//muons

            if (!hasLepton)
                continue;

            p4MuJet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );

            double btag   = jetIter -> bDiscriminator( Tagger );

            if (btag > btag_cut_Tagger ) TaggedJet = true;

            histos.FillHistos("n", p4MuJet, ptrel, JetFlavor, TaggedJet );

            // find away jet
            ////////////////////////////
            bool AwayTaggedJet = false;
            bool AwayMuonJet = false;
            TLorentzVector p4AwayMuon;
            TLorentzVector p4AwayTagged;
        
            for (vector< pat::Jet >::const_iterator awayjetIter = jetHandle->begin();
                 awayjetIter != jetHandle->end(); ++awayjetIter)
            {
                p4AwayJet.SetPtEtaPhiE(awayjetIter->pt(), awayjetIter->eta(), awayjetIter->phi(), awayjetIter->energy() );

                // skip muon in jet
                if ( p4AwayJet == p4MuJet ) continue;

                hstore->hist("awayjet_pt")->Fill (awayjetIter->pt()); // just for testing

                if ( !AwayTaggedJet )
                {
                    btag  = awayjetIter -> bDiscriminator( awayTagger );
                    //std::cout <<" bDiscriminator trackCountingHighEffBJetTags = " <<jetBDiscr_track_count_high_eff << std::endl;
                    if ( btag > btag_cut_awayTagger) //loose operating point
                    {
                        AwayTaggedJet = true;
                        p4AwayTagged = p4AwayJet;
                    }

                }

                if ( !AwayMuonJet )
                {
                    mu_highest_pt = 0;
                    for (vector< pat::Muon >::const_iterator muonIter = muonHandle->begin();
                         muonIter != muonHandle->end(); ++muonIter)
                    {
                        hstore->hist("muon_pt")->Fill (muonIter->innerTrack()->pt());

                        // find a muon in a jet
                        //check deltar
                        reco::TrackRef trackMuon = muonIter->innerTrack();
                        math::XYZTLorentzVector trackMuonP4(trackMuon->px(),
                                                            trackMuon->py(),
                                                            trackMuon->pz(),
                                                            sqrt(trackMuon->p() * trackMuon->p() + 0.1057*0.1057));
                        double deltaR  = ROOT::Math::VectorUtil::DeltaR(jetIter->p4().Vect(), trackMuonP4.Vect() );
                        TVector3 tmpvecOrg(awayjetIter->p4().Vect().X(),awayjetIter->p4().Vect().Y(),  awayjetIter->p4().Vect().Z());
                        TVector3 tmpvec;
                        tmpvec = tmpvecOrg;
                        TVector3 leptonvec(trackMuon->px(), trackMuon->py(),trackMuon->pz());
                        tmpvec += leptonvec;
                        double awayptrel = leptonvec.Perp(tmpvec);
                        // muon in jet cuts
                        if ( (deltaR >= 0.4) || (ptreltmp <= -1.0) ) continue;
                        // now we have a good muon in a jet
                        AwayMuonJet = true;
                        // pick the leading muon inside the jet
                        if ( trackMuon->pt() > mu_highest_pt )
                        {
                            mu_highest_pt = muonIter->pt();
                            p4AwayMuon.SetPtEtaPhiE(trackMuon->pt(),
                                                    trackMuon->eta(),
                                                    trackMuon->phi(),
                                                    sqrt(trackMuon->p() * trackMuon->p() + 0.1057*0.1057));
                            // recalculate pTrel
                            tmpvec = tmpvecOrg;
                            leptonvec.SetXYZ(trackMuon->px(), trackMuon->py(),trackMuon->pz());
                            tmpvec += leptonvec;
                            awayptrel = leptonvec.Perp(tmpvec);
                        }
                    }
                } // away muon jet
            } // close away jet loop

            if (!AwayTaggedJet)
                continue;

            histos.FillHistos("p", p4MuJet, ptrel, JetFlavor, TaggedJet );
            hstore->hist("deltaRnearjet")->Fill ( p4MuJet.DeltaR( p4AwayTagged ) );
            hstore->hist("deltaPhi")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

            string suffix;
            switch(JetFlavor)
            {
                case 5: suffix = "_b";
                        break;

                case 4: suffix = "_c";
                        break;

                case 21: suffix = "_g";
                         break;

                case 3 : // fall through
                case 2 : // fall through
                case 1 : // fall through
                         suffix = "_l";
                         break;
            }

            if (suffix.length())
            {
                hstore->hist("deltaPhi" + suffix)->Fill(p4AwayTagged.Phi() - p4MuJet.Phi());
            }
        }//jets
    }//events

    std::cout << " done. " << std::endl;
    
    std::string m_outputName  = parser.stringValue  ("outputFile");
    if (optutl::CommandLineParser::kStringVector ==
        parser.hasOption("inputFiles"))
    {
        hstore->write (m_outputName,
                       parser.argVec(),
                       parser.stringVector ("inputFiles"));
    }
    else
    {
        hstore->write (m_outputName,
                       parser.argVec());
    }
    
    
    return 0;
    
}

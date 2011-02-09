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
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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
using reco::Vertex;

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
        // Event is allowed. Process it normally
        //
        if (maxNevents && nentry > maxNevents) break;

        // Test if event is among the ones, specified in JSON file
        //
        if (!jsonContainsEvent (lumiVector, events) )
            continue;

        ++nentry;

        if ( outputEvery!=0 && (nentry % outputEvery) == 0 ) 
            cout << "Processing Event: " << nentry << endl;

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
        } // end processing Trigger

        // load object collections
        //
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

        typedef vector<Vertex> PVCollection;
        fwlite::Handle<PVCollection> primaryVertices;
        primaryVertices.getByLabel(events, "offlinePrimaryVertices");

        if (!primaryVertices.isValid() ||
            !primaryVertices->size())
        {
            cerr << "Failed to extract primary vertices or collection is empty." << endl;

            break;
        }

        const Vertex &primaryVertex = primaryVertices->at(0);

        // Create All Muons plots
        //
        {
            int muons = 0;
            for (vector<pat::Muon>::const_iterator muon = muonHandle->begin();
                 muonHandle->end() != muon;
                 ++muon)
            {
                const double dz = muon->vz() - primaryVertex.z();

                if (1 < muon->numberOfMatches() &&
                    abs(dz) < 2)
                {
                    hstore->hist("all_mu_pt")->Fill(muon->pt());
                    hstore->hist("all_mu_phi")->Fill(muon->phi());
                    hstore->hist("all_mu_eta")->Fill(muon->eta());
                    hstore->hist("all_mu_charge")->Fill(muon->charge());
                    hstore->hist("all_mu_dz")->Fill(dz);
                    hstore->hist("all_mu_ip")->Fill(muon->dB());

                    ++muons;
                }
            }

            hstore->hist("all_mu_number")->Fill(muons);
        }

        hstore->hist("all_jet_number")->Fill(jetHandle->size());
        // Loop over jets
        //
        for (vector< pat::Jet >::const_iterator jetIter = jetHandle->begin();
             jetIter != jetHandle->end(); ++jetIter)
        {
            TLorentzVector p4Jet;
            TLorentzVector p4MuJet;
            TLorentzVector p4Muon;
            TLorentzVector p4AwayJet;

            bool TaggedJet = false;
            hstore->hist("jet_pt")->Fill (jetIter->pt()); // just for testing

            // All jet plots
            //
            hstore->hist("all_jet_pt")->Fill(jetIter->pt());
            hstore->hist("all_jet_eta")->Fill(jetIter->eta());
            hstore->hist("all_jet_phi")->Fill(jetIter->phi());

            // get MC flavor of jet
            //
            int JetFlavor = abs( jetIter->partonFlavour() );

            p4Jet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );
            bool hasLepton = false;

            //int tmptotmuon = 0;
            // loop over muons
            ////////////////////////////////
            double mu_highest_pt = 0;
            double ptrel = 0;
            double ptreltmp = 0;
            double mu_ip = 0;
            double mu_dz = 0;
            double mu_charge = 0;
            int muonsInJet = 0;

            // Attempt to find muon inside the jet
            //
            for (vector< pat::Muon >::const_iterator muonIter = muonHandle->begin();
                 muonIter != muonHandle->end(); ++muonIter)
            {
                // {samvel} loop over jets -> loop over muons and plot Pt of
                //          each muon. What is the expected result???
                //
                hstore->hist("muon_pt")->Fill (muonIter->innerTrack()->pt());

                const double dz = muonIter->vz() - primaryVertex.z();
                if (1 >= muonIter->numberOfMatches() ||
                    abs(dz) >= 2 )

                    continue;

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
                //
                if (deltaR >= 0.4 || deltaR < 0.01 ||
                    ptreltmp <= -1.0)

                    continue;

                hasLepton = true;

                ++muonsInJet;

                // It may turn out that another muon was found so far. Then
                // pick the one with highest pT.
                //
                if ( muonIter->pt() > mu_highest_pt )
                {
                    mu_highest_pt = muonIter->pt();
                    p4Muon.SetPtEtaPhiE(muonIter->pt(),
                                        muonIter->eta(),
                                        muonIter->phi(),
                                        muonIter->energy());
                    // recalculate pTrel
                    // {samvel} Why??? pTrel is already available at: ptreltmp
                    //
                    tmpvec = tmpvecOrg;
                    leptonvec.SetXYZ(muonIter->px(), muonIter->py(),muonIter->pz());
                    tmpvec += leptonvec;
                    ptrel = leptonvec.Perp(tmpvec);  // maximum

                    mu_ip = muonIter->dB();
                    mu_dz = abs(muonIter->vz() - primaryVertex.z());
                    mu_charge = muonIter->charge();
                }

                hstore->hist("jet_pTrel")->Fill (ptrel);
                hstore->hist("jet_deltaR")->Fill (deltaR);

                // b-jet
                //
                if (5 == JetFlavor)
                {
                    hstore->hist("jet_pTrel_b")->Fill(ptrel);
                    hstore->hist("jet_deltaR_b")->Fill(deltaR);
                }

                // c-jet
                //
                if (4 == JetFlavor)
                {
                    hstore->hist("jet_pTrel_c")->Fill(ptrel);
                    hstore->hist("jet_deltaR_c")->Fill(deltaR);

                    hstore->hist("jet_pTrel_cl")->Fill(ptrel);
                    hstore->hist("jet_deltaR_cl")->Fill(deltaR);
                }

                // uds-jet
                //
                if (4 > JetFlavor && 0 < JetFlavor)
                {
                    hstore->hist("jet_pTrel_l")->Fill(ptrel);
                    hstore->hist("jet_deltaR_l")->Fill(deltaR);

                    hstore->hist("jet_pTrel_cl")->Fill(ptrel);
                    hstore->hist("jet_deltaR_cl")->Fill(deltaR);

                    hstore->hist("jet_pTrel_uds")->Fill(ptrel);
                    hstore->hist("jet_deltaR_uds")->Fill(deltaR);
                }

                // gluon
                //
                if (21 == JetFlavor)
                {
                    hstore->hist("jet_pTrel_g")->Fill(ptrel);
                    hstore->hist("jet_deltaR_g")->Fill(deltaR);

                    hstore->hist("jet_pTrel_cl")->Fill(ptrel);
                    hstore->hist("jet_deltaR_cl")->Fill(deltaR);

                    hstore->hist("jet_pTrel_l")->Fill(ptrel);
                    hstore->hist("jet_deltaR_l")->Fill(deltaR);
                }
            }//muons

            // Proceed if Jet contains lepton otherwise go to next jet.
            //
            if (!hasLepton)
                continue;

            p4MuJet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );

            double btag   = jetIter -> bDiscriminator( Tagger );

            if (btag > btag_cut_Tagger )
                TaggedJet = true;

            hstore->hist("mu_injet_pt")->Fill(p4Muon.Pt());
            hstore->hist("mu_injet_phi")->Fill(p4Muon.Phi());
            hstore->hist("mu_injet_eta")->Fill(p4Muon.Eta());
            hstore->hist("mu_injet_ip")->Fill(mu_ip);
            hstore->hist("mu_injet_dz")->Fill(mu_dz);
            hstore->hist("mu_injet_charge")->Fill(mu_charge);
            hstore->hist("mu_injet")->Fill(muonsInJet);

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
                //
                if ( p4AwayJet == p4MuJet ) continue;

                // {samvel} what is this? For each jet we plot awayjet. So, for
                // 5 jets with 2 jet with muon 4^2 = 16 entries will be added.
                //
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
            } // close away jet loop

            // Skip if away tagged jet is not found
            if (!AwayTaggedJet)
                continue;

            histos.FillHistos("p", p4MuJet, ptrel, JetFlavor, TaggedJet );
            hstore->hist("deltaRnearjet")->Fill ( p4MuJet.DeltaR( p4AwayTagged ) );
            hstore->hist("deltaPhi")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

            if (5 == JetFlavor)
                hstore->hist("deltaPhi_b")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

            // c-jet
            //
            if (4 == JetFlavor)
            {
                hstore->hist("deltaPhi_c")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

                hstore->hist("deltaPhi_cl")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );
            }

            // uds-jet
            //
            if (4 > JetFlavor && 0 < JetFlavor)
            {
                hstore->hist("deltaPhi_l")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

                hstore->hist("deltaPhi_cl")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

                hstore->hist("deltaPhi_uds")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );
            }

            // gluon
            //
            if (21 == JetFlavor)
            {
                hstore->hist("deltaPhi_g")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

                hstore->hist("deltaPhi_l")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );

                hstore->hist("deltaPhi_cl")->Fill (p4AwayTagged.Phi() - p4MuJet.Phi() );
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

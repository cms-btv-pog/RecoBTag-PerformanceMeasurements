
#include "RecoBTag/PerformanceMeasurements/interface/PerformanceAnalyzerWithElectrons.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// reco track and vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


// HepPDT // for simtracks
#include "HepPDT/ParticleID.hh"

//#include "SimGeneral/HepPDT/interface/HepPDTable.h"
//#include "SimGeneral/HepPDT/interface/HepParticleData.h"

// Root
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "TGraphErrors.h"

#include <map>
#include <cmath>

using namespace edm;
using namespace reco;

//
// constructors and destructor
//
PerformanceAnalyzerWithElectrons::PerformanceAnalyzerWithElectrons(const ParameterSet& iConfig)
{

    fbadeventscounter = 0;

    analyzer_ = "PerformanceAnalyzerWithElectrons"; // name of this analyzer

    if (fdebug) std::cout << "[PerformanceAnalyzerWithElectrons] call constructor" << std::endl;

    //now do whatever initialization is needed
    //simG4_=iConfig.getParameter<edm::InputTag>( "simG4" );

    // open output file to store results
    outputFile_  = iConfig.getUntrackedParameter<std::string>("outputFile");
    // create or update output root file
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
    // verbose std output
    fdebug = iConfig.getUntrackedParameter<bool>("debug", false);
    // default tree
    ftree = new TTree("summary","summary");
    ftree->AutoSave();


    fS8evt = new BTagEvent();
    ftree->Branch("s8.","BTagEvent",&fS8evt,64000,1); //system8 branch


    // get list of tracks
    recoTrackList_ = iConfig.getUntrackedParameter<std::string >("TrackCollection");
    // get list of PV
    //recoVtxList_ = iConfig.getUntrackedParameter<std::string >("PrimaryVertexCollection");

    ElectronCollectionTags_ = iConfig.getParameter<std::string>("Electrons");

    CaloJetCollectionTags_ = iConfig.getParameter<std::string>("Jets");

    GenJetCollectionTags_ = iConfig.getParameter<std::string>("GenJets");

    SimTrkCollectionTags_ = iConfig.getParameter<std::string>("SimTracks");

    bTagTrackEventIPTagInfos_ = iConfig.getParameter<std::string>("bTagTrackEventIPtagInfos");

    //jetFlavourIdentifier_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters"));
    //jetFlavourIdentifier2_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters2"));


    StoreTrackProba_ = iConfig.getParameter<bool>("StoreTrackProba");

    const edm::ParameterSet& jetCuts = iConfig.getParameter<edm::ParameterSet>("jetcuts");
    // jet cuts
    MinJetPt_ = jetCuts.getParameter<double>("MinPt");
    MaxJetEta_ = jetCuts.getParameter<double>("MaxEta");
    MinEmFraction_ = jetCuts.getParameter<double>("MinEmFraction");
    MaxEmFraction_ = jetCuts.getParameter<double>("MaxEmFraction");
    MinDeltaR_ = jetCuts.getParameter<double>("MinDeltaR");
    MinPtRel_  = jetCuts.getParameter<double>("MinPtRel");
    MaxRatio_ = jetCuts.getParameter<double>("MaxRatio");

    const edm::ParameterSet& electronCuts = iConfig.getParameter<edm::ParameterSet>("electroncuts");
    // electron cuts
    MinElectronPt_ = electronCuts.getParameter<double>("MinElectronPt");
    MaxElectronEta_ = electronCuts.getParameter<double>("MaxElectronEta");
    MaxElectronChi2_ = electronCuts.getParameter<double>("MaxElectronChi2");
    MinElectronNHits_ = electronCuts.getParameter<int>("MinElectronNHits");

    //
    // ue jetcorrections or not, and in cae which label to use
    //
    useJetCorr_ = iConfig.getParameter<bool>("useJetCorrections");
    if (useJetCorr_ == true)
    {
        jetCorrLabel_ =  iConfig.getParameter<std::string>("jetCorrectionsLabel");
        std::cout<<" Use JetCorrections with Label "<<jetCorrLabel_ <<std::endl;
    }
    else
    {
        jetCorrLabel_  = "Fake";
        std::cout<<" Do NOT use JetCorrections."<<std::endl;
    }



    // get list of taggers

    // Flavour identification
    flavourMatchOptionf = iConfig.getParameter<std::string>( "flavourMatchOption" );
    if (flavourMatchOptionf == "fastMC")
    {
        flavourSourcef = iConfig.getParameter<edm::InputTag>("flavourSource");
    }
    else if (flavourMatchOptionf == "genParticle")
    {
        flavourSourcef = iConfig.getParameter<edm::InputTag> ("flavourSource");
    }


    //
    // get operating points
    std::vector<edm::ParameterSet> config = iConfig.getUntrackedParameter<std::vector<edm::ParameterSet > >("OperatingPointsList");
    if (fdebug) std::cout << " get operating points, total list: " << config.size() << std::endl;

    for (std::vector<edm::ParameterSet>::const_iterator it = config.begin(); it != config.end() ; ++it)
    {
        std::string aalias = (*it).getUntrackedParameter<std::string> ("alias");
		edm::InputTag atag = (*it).getUntrackedParameter<edm::InputTag> ("collection");
        double min = (*it).getUntrackedParameter<double> ("MinimumDiscriminator");
        double max = (*it).getUntrackedParameter<double> ("MaximumDiscriminator");
		std::vector<edm::ParameterSet> avec = (*it).getUntrackedParameter<std::vector<edm::ParameterSet> >("OperatingPoints");
		std::map< std::string, double > mapOP;
        for (std::vector<edm::ParameterSet>::const_iterator imapOP = avec.begin(); imapOP != avec.end(); ++imapOP)
        {
			 mapOP[(*imapOP).getUntrackedParameter<std::string> ("name")] =
                (*imapOP).getUntrackedParameter<double> ("cut");
        }

		WorkingPoint tmpwp((*it).getUntrackedParameter<edm::InputTag> ("collection"),
                           aalias,
                           min,
                           max,
                           mapOP );

		wp.push_back(tmpwp);
		if (fdebug) (wp.end()-1)->print();
		//wp_map[alias] = tmpwp;

		TaggerPerformances_[aalias].Set(aalias);
        TaggerPerformances_[aalias].SetMinDiscriminator(min);
        TaggerPerformances_[aalias].SetMaxDiscriminator(max);

		//double acut = (*it).getUntrackedParameter<double> ("cut");
		//WorkingPoint tmp(atag, aname, acut);
        //wp.push_back( tmp );
        //if (fdebug) (wp.end()-1)->print();
    }


    fnselectors= wp.size();

    // get tagger for away jet
    fAwayJetTagger = iConfig.getParameter<std::string>("AwayJetTagger");

    // write performance plots?
    fWritePerformancePlots = iConfig.getParameter< bool > ("WritePerformancePlots");

    // include weights?
    fWeightHistograms = iConfig.getParameter< bool > ("WeightHistograms");
    fStoreWeightsInNtuple = iConfig.getParameter< bool > ("StoreWeightsInNtuple");
    fStorePtHat = iConfig.getParameter< bool > ("StorePtHat");

    topdir = rootFile_->mkdir("Histograms");
    topdir->cd();
    topdir->mkdir("ptrel");
    topdir->cd();
    topdir->mkdir("MCtruth");
    topdir->cd();
    topdir->mkdir("electron_in_jet");
    topdir->cd();
	if (fdebug) std::cout<< " ROOT directories created." << std::endl;

    histcounterf = new TH1I("histcounterf","counter",5,0,5);
    histcounterf->SetCanExtend(TH1::kAllAxes);
    rootFile_->cd();

    // initialize histograms
    EffHistos     = new BTagHistograms();
    PtrelHistos   = new BTagHistograms();
    ElecjetHistos   = new BTagHistograms();
    AwayjetHistos = new BTagHistograms();
    TaggedElecjetHistos   = new BTagHistograms();
    TaggedAwayjetHistos = new BTagHistograms();
    ElecjetHistos_mc   = new BTagHistograms();
    AwayjetHistos_mc = new BTagHistograms();
    TaggedElecjetHistos_mc   = new BTagHistograms();
    TaggedAwayjetHistos_mc = new BTagHistograms();

    EffHistos->Init("efficiencies");
    PtrelHistos->Init("ptrel");
    ElecjetHistos->Init("n");
    AwayjetHistos->Init("p");
    ElecjetHistos_mc->Init("n","b");
    AwayjetHistos_mc->Init("p","b");
    ElecjetHistos_mc->Init("n","cl");
    AwayjetHistos_mc->Init("p","cl");
    ElecjetHistos_mc->Init("n","c");
    AwayjetHistos_mc->Init("p","c");
    ElecjetHistos_mc->Init("n","l");
    AwayjetHistos_mc->Init("p","l");
    for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it!=wp.end(); ++it)
    {
		std::map<std::string, double > list_cuts = (*it).list();

		for(std::map<std::string, double >::const_iterator icut = list_cuts.begin(); icut != list_cuts.end(); ++icut) {
			std::string aalias = icut->first;
			EffHistos->Init("efficiencies",aalias);
			PtrelHistos->Init("ptrel",aalias);

			TaggedElecjetHistos->Init("ntag",aalias);
			TaggedAwayjetHistos->Init("ptag",aalias);
			TaggedElecjetHistos_mc->Init("ntag","b",aalias);
			TaggedAwayjetHistos_mc->Init("ptag","b",aalias);
			TaggedElecjetHistos_mc->Init("ntag","cl",aalias);
			TaggedAwayjetHistos_mc->Init("ptag","cl",aalias);
			TaggedElecjetHistos_mc->Init("ntag","c",aalias);
			TaggedAwayjetHistos_mc->Init("ptag","c",aalias);
			TaggedElecjetHistos_mc->Init("ntag","l",aalias);
			TaggedAwayjetHistos_mc->Init("ptag","l",aalias);
		}
    }

	if (fdebug) std::cout << "Histograms initialized" << std::endl;

	/*
    fperformanceTC2trk.Set("TC2trk");
    fperformanceTC3trk.Set("TC3trk");
    fperformanceMTC2trk.Set("MTC2trk");
    fperformanceMTC3trk.Set("MTC3trk");
    fperformanceTP.Set("TP");
    fperformanceJBP.Set("JBP");
    fperformanceSMT.Set("SMT");
    fperformanceSSV.Set("SSV");
    fperformanceCSV.Set("CSV");

    fperformanceTC2trk.SetMinDiscriminator(-1);
    fperformanceTC2trk.SetMaxDiscriminator(15);
    fperformanceTC3trk.SetMinDiscriminator(-1);
    fperformanceTC3trk.SetMaxDiscriminator(15);
    fperformanceMTC2trk.SetMinDiscriminator(-1);
    fperformanceMTC2trk.SetMaxDiscriminator(15);
    fperformanceMTC3trk.SetMinDiscriminator(-1);
    fperformanceMTC3trk.SetMaxDiscriminator(15);
    fperformanceTP.SetMinDiscriminator(0);
    fperformanceTP.SetMaxDiscriminator(1);
    fperformanceJBP.SetMinDiscriminator(0);
    fperformanceJBP.SetMaxDiscriminator(1);
    fperformanceSMT.SetMinDiscriminator(0);
    fperformanceSMT.SetMaxDiscriminator(1);
    fperformanceSSV.SetMinDiscriminator(0);
    fperformanceSSV.SetMaxDiscriminator(10);
    fperformanceCSV.SetMinDiscriminator(0);
    fperformanceCSV.SetMaxDiscriminator(1);
	*/

    feventcounter = 0;

}


PerformanceAnalyzerWithElectrons::~PerformanceAnalyzerWithElectrons()
{
    rootFile_->cd();
    ftree->Write();

    topdir->cd();
    topdir->cd("MCtruth");
    EffHistos->Save();
	std::vector< TGraph* > gVector;
    if (fWritePerformancePlots)
    {


		for (std::map<std::string, S8bPerformance>::const_iterator iperf = TaggerPerformances_.begin(); iperf!= TaggerPerformances_.end(); ++iperf )
    {

        S8bPerformance Perf = iperf->second;

        Perf.Eval();


        TGraphErrors *gTb = Perf.EfficiencyGraph("b");
        TGraphErrors *gTc = Perf.EfficiencyGraph("c");
        TGraphErrors *gTl = Perf.EfficiencyGraph("udsg");
        gVector.push_back( gTb );
        gVector.push_back( gTc );
        gVector.push_back( gTl );

		TGraph *discTl = Perf.DiscriminatorGraph("udsg");
        discTl->Sort();
        gVector.push_back( discTl );

	}

    }

	if (fWritePerformancePlots)
    {
		for (std::vector< TGraph* >::const_iterator iv = gVector.begin(); iv != gVector.end(); ++iv )
		{
			(*iv)->Write();
		}
	}

    topdir->cd("ptrel");
    PtrelHistos->Save();
    topdir->cd("electron_in_jet");
    ElecjetHistos->Save();
    AwayjetHistos->Save();
    TaggedElecjetHistos->Save();
    TaggedAwayjetHistos->Save();
    topdir->cd("MCtruth");
    ElecjetHistos_mc->Save();
    AwayjetHistos_mc->Save();
    TaggedElecjetHistos_mc->Save();
    TaggedAwayjetHistos_mc->Save();

	topdir->cd();
    topdir->cd("Histograms");
    histcounterf->Write();

    topdir->Write();

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

    rootFile_->Close();

    delete fS8evt;


    //delete ftree;
    //delete rootFile_;
}




//
// member functions
//
void PerformanceAnalyzerWithElectrons::beginJob(edm::EventSetup const& iSetup)
{
    std::cout << analyzer_ << " begin Job" << std::endl;
    rootFile_->cd();
}


void PerformanceAnalyzerWithElectrons::endJob()
{
    std::cout << analyzer_ << " Total events processed: " << feventcounter << std::endl;
    rootFile_->cd();
}



//______________________________________________________________________________________________________________________
void PerformanceAnalyzerWithElectrons::FillPerformance(reco::CaloJet jet, int JetFlavor, const edm::Event& event)
{

    std::map< std::string, bool > aMap;
    int ith_tagged = -1;

    std::map<std::string, bool> mymap;

    for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it != wp.end(); ++it)
    {

		edm::Handle<reco::JetTagCollection > jetTags;
		event.getByLabel((*it).inputTag(),jetTags);
		std::string alias = (*it).alias();

		std::string moduleLabel = (jetTags).provenance()->moduleLabel();
		if (mymap.find(moduleLabel) != mymap.end()) continue;
		mymap[moduleLabel] = true;

		ith_tagged = PFTools::TaggedJet(jet,jetTags);

		if (ith_tagged == -1) continue;

		TaggerPerformances_[alias].Add( (*jetTags)[ith_tagged].second, JetFlavor );

		/*
        edm::Handle<reco::JetTagCollection > jetTags;

        event.getByLabel((*it).inputTag(),jetTags);
        std::string moduleLabel = (jetTags).provenance()->moduleLabel();
        if (mymap.find(moduleLabel) != mymap.end()) continue;

        mymap[moduleLabel] = true;


        ith_tagged = PFTools::TaggedJet(jet,jetTags);


        if (ith_tagged == -1) continue;


        // *********************************
        // Track Counting taggers
        // *********************************
        if ( moduleLabel == "trackCountingHighEffBJetTags" )
        {

            fperformanceTC2trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
        else if ( moduleLabel == "trackCountingHighPurBJetTags" )
        {

            fperformanceTC3trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }

        // *********************************
        // Modified Track Counting taggers
        // *********************************
        if ( moduleLabel == "modifiedtrackCountingHighEffBJetTags" )
        {

            fperformanceMTC2trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
        else if ( moduleLabel == "modifiedtrackCountingHighPurBJetTags" )
        {

            fperformanceMTC3trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }

        // *********************************
        // Jet Probability taggers
        // *********************************
        else if ( moduleLabel == "jetProbabilityBJetTags" )
        {

            fperformanceTP.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
        // jetBProbabilityBJetTags
        else if ( moduleLabel == "jetBProbabilityBJetTags" )
        {

            fperformanceJBP.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }

        // *********************************
        // Soft Electron Tagger
        // *********************************
        else if ( moduleLabel == "softElectronBJetTags" )
        {

            fperformanceSMT.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }

        // *********************************
        // Secondary Vertex Tagger
        // *********************************
        else if ( moduleLabel == "simpleSecondaryVertexBJetTags" )
        {

            fperformanceSSV.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
        else if ( moduleLabel == "combinedSecondaryVertexBJetTags" )
        {

            fperformanceCSV.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }

		*/

    }


}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzerWithElectrons::FillEff(TLorentzVector p4ElecJet, int JetFlavor, std::map<std::string, bool> aMap, double weight)
{

    std::string flavor = "";
    for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap )
    {

        if (imap->second)
        {

            std::string tag = "_"+imap->first;
            flavor = "";

            EffHistos->Fill1d("jet_pt"+tag,p4ElecJet.Pt(), weight );
            EffHistos->Fill1d("jet_eta"+tag,p4ElecJet.Eta(), weight );

            if ( JetFlavor == 5 )
            {
                flavor = "_b";
            }
            else if ( JetFlavor == 4 )
            {
                flavor = "_c";
            }
            else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 )
            {
                flavor = "_udsg";
            }

            if (flavor!="")
            {
                EffHistos->Fill1d("jet_pt"+flavor+tag,p4ElecJet.Pt(), weight );
                EffHistos->Fill1d("jet_eta"+flavor+tag,p4ElecJet.Eta(),weight );
            }
        }
    }

    EffHistos->Fill1d("jet_pt",p4ElecJet.Pt(),weight );
    EffHistos->Fill1d("jet_eta",p4ElecJet.Eta(),weight );

    if ( JetFlavor == 5 )
    {
        flavor = "_b";
    }
    else if ( JetFlavor == 4 )
    {
        flavor = "_c";
    }
    else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 )
    {
        flavor = "_udsg";
    }

    if (flavor!="")
    {
        EffHistos->Fill1d("jet_pt"+flavor,p4ElecJet.Pt(), weight );
        EffHistos->Fill1d("jet_eta"+flavor,p4ElecJet.Eta(),weight );
    }

}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzerWithElectrons::FillPtrel(double ptrel, int JetFlavor, std::map<std::string, bool> aMap, double weight)
{

    std::string flavor = "";

    for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap )
    {

        if (imap->second)
        {

            std::string tag = "_"+imap->first;
            flavor = "";

            PtrelHistos->Fill1d("jet_ptrel"+tag,ptrel,weight);

            if ( JetFlavor == 5 )
            {
                flavor = "_b";
            }
            else if ( JetFlavor == 4 )
            {
                flavor = "_c";
            }
            else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 )
            {
                flavor = "_udsg";
            }

            if (flavor!="") PtrelHistos->Fill1d("jet_ptrel"+flavor+tag,ptrel,weight);

        }
    }

    PtrelHistos->Fill1d("jet_ptrel",ptrel,weight);

    if ( JetFlavor == 5 )
    {
        flavor = "_b";
    }
    else if ( JetFlavor == 4 )
    {
        flavor = "_c";
    }
    else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 )
    {
        flavor = "_udsg";
    }

    if (flavor!="") PtrelHistos->Fill1d("jet_ptrel"+flavor,ptrel,weight);

}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzerWithElectrons::FillHistos(std::string type, TLorentzVector p4ElecJet, double ptrel,
                                     int JetFlavor, std::map<std::string, bool> aMap, double weight)
{
    if ( type == "n")
    {
        ElecjetHistos->Fill2d(type+"_pT",p4ElecJet.Pt(),ptrel,weight);
        ElecjetHistos->Fill2d(type+"_eta",TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        if ( JetFlavor == 5 )
        {
            std::string flavor = "b";
            ElecjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            ElecjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
        {
            std::string flavor = "cl";
            ElecjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            ElecjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
        if ( JetFlavor == 4 )
        {
            std::string flavor = "c";
            ElecjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            ElecjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
        {
            std::string flavor = "l";
            ElecjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            ElecjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
    }
    else if ( type == "p")
    {

        AwayjetHistos->Fill2d(type+"_pT",p4ElecJet.Pt(),ptrel,weight);
        AwayjetHistos->Fill2d(type+"_eta",TMath::Abs(p4ElecJet.Eta()),ptrel, weight);
        if ( JetFlavor == 5 )
        {
            std::string flavor = "b";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
        {
            std::string flavor = "cl";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
        if ( JetFlavor == 4 )
        {
            std::string flavor = "c";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
        {
            std::string flavor = "l";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4ElecJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
        }
    }

    for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap )
    {

        if ( imap->second )
        {

            if ( type == "n")
            {

                TaggedElecjetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                TaggedElecjetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);

                if ( JetFlavor == 5 )
                {
                    std::string flavor = "b_";
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
                {
                    std::string flavor = "cl_";
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
                if ( JetFlavor == 4 )
                {
                    std::string flavor = "c_";
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
                {
                    std::string flavor = "l_";
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedElecjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
            }
            else if ( type == "p")
            {


                TaggedAwayjetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                TaggedAwayjetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);

                if ( JetFlavor == 5 )
                {
                    std::string flavor = "b_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
                {
                    std::string flavor = "cl_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
                if ( JetFlavor == 4 )
                {
                    std::string flavor = "c_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
                {
                    std::string flavor = "l_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4ElecJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4ElecJet.Eta()),ptrel,weight);
                }
            }
        }
    }

}

reco::JetFlavour PerformanceAnalyzerWithElectrons::getMatchedParton(const reco::CaloJet &jet)
{
    reco::JetFlavour jetFlavour;

    if (flavourMatchOptionf == "fastMC")
    {
        //
        // disabled
        //

        //        jetFlavour.underlyingParton4Vec(jet.p4());
        //const edm::RefToBase<reco::Jet> & caloRefTB = jetTag.jet();
        //const reco::CaloJetRef & caloRef = jet.castTo<reco::CaloJetRef>();
        //jetFlavour.flavour(flavoursMapf[caloRef]);

    }
    else if (flavourMatchOptionf == "genParticle")
    {

        //
        // try and use directly the map
        //

        //    RefToBase<Jet> testJet(jet);

        for ( JetFlavourMatchingCollection::const_iterator j  = theJetPartonMapf->begin();
                j != theJetPartonMapf->end();
                j ++ )
        {
            RefToBase<Jet> aJet  = (*j).first;
            //      const JetFlavour aFlav = (*j).second;
            if ( fabs(aJet->phi() - jet.phi()) < 1.e-5 && fabs(aJet->eta() - jet.eta())< 1.e-5 )
            {
                // matched
                jetFlavour = reco::JetFlavour (aJet->p4(), math::XYZPoint(0,0,0), (*j).second.getFlavour());
            }
        }

        return jetFlavour;


    }

    return jetFlavour;
}


// ------------ method called to produce the data  ------------
void
PerformanceAnalyzerWithElectrons::analyze(const Event& iEvent, const EventSetup& iSetup)
{

    // count
    histcounterf->Fill("Processed", 1);
    bool NjetsCut = false;
    bool NelecCut = false;
    bool Nelec_in_jetCut = false;
    bool Nelec_in_jet_away_taggedCut = false;

    // Trakcs
    Handle<reco::TrackCollection> recTrks;
    iEvent.getByLabel(recoTrackList_, recTrks);

    // initialize flavour identifiers
    edm::Handle<JetFlavourMatchingCollection> jetMC;

    if (flavourMatchOptionf == "fastMC")
    {
        iEvent.getByLabel(flavourSourcef, jetMC);
        for (JetFlavourMatchingCollection::const_iterator iter =
                    jetMC->begin(); iter != jetMC->end(); iter++)
            flavoursMapf.insert(std::pair <const edm::RefToBase<reco::Jet>, unsigned int>((*iter).first, (unsigned int)((*iter).second).getFlavour()));
        //      flavoursMapf[(*iter).first]= (unsigned int)((*iter).second).getFlavour();
    }
    else if (flavourMatchOptionf == "genParticle")
    {
        iEvent.getByLabel (flavourSourcef, theJetPartonMapf);
    }


    // generator candidates
    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);

    // Electrons
    Handle<reco::GsfElectronCollection> electronsColl;
    iEvent.getByLabel(ElectronCollectionTags_, electronsColl);

    // Calo Jets
    Handle<reco::CaloJetCollection> jetsColl;
    iEvent.getByLabel(CaloJetCollectionTags_, jetsColl);

    // Get the bTagTrackEventIPTagInfo collection
    Handle<std::vector<reco::TrackIPTagInfo> > bTagTrackEventIPTagInfos;

    if ( !bTagTrackEventIPTagInfos_.empty() )
        iEvent.getByLabel(bTagTrackEventIPTagInfos_, bTagTrackEventIPTagInfos);


    // initialize jet corrector
    const JetCorrector *acorrector = 0;
    if (useJetCorr_ == true)
    {
        acorrector = JetCorrector::getJetCorrector(jetCorrLabel_,iSetup);
    }


    Handle<reco::GenJetCollection> genjetsColl;
    iEvent.getByLabel(GenJetCollectionTags_, genjetsColl);


    const reco::CaloJetCollection recoJets =   *(jetsColl.product());
    const reco::GenJetCollection  genJets  =   *(genjetsColl.product());
    const reco::GsfElectronCollection    recoElectrons =  *(electronsColl.product());

    //const reco::VertexCollection recoPV = *(recVtxs.product());

    // MC

    //Handle<SimTrackContainer> simTrks;
    //iEvent.getByLabel( simG4_, simTrks);

    Handle<std::vector<reco::TrackIPTagInfo> > tagInfo;
    iEvent.getByLabel("impactParameterTagInfos", tagInfo);
    // WEIGHTS
    double weight = 1.;

    Handle< double> weightHandle;

    if (fWeightHistograms || fStoreWeightsInNtuple)
    {

        iEvent.getByLabel ("csaweightproducer","weight", weightHandle);

    }

    if (fWeightHistograms) weight = *weightHandle;

    if (fStoreWeightsInNtuple) fS8evt->evt_weight =  (*weightHandle);
    else fS8evt->evt_weight = 1.;


    if (fStorePtHat)
    {
        edm::Handle<int> genProcessID;
        iEvent.getByLabel( "genEventProcID", genProcessID );
        double processID = *genProcessID;
        edm::Handle<double> genEventScale;
        iEvent.getByLabel( "genEventScale", genEventScale );
        double pthat = *genEventScale;

        if (processID == 11 || processID == 12 || processID == 13
                || processID == 28 || processID == 68 || processID == 53) fS8evt->ptHat = pthat;
        else fS8evt->ptHat = -99;
    }


    fS8evt->Reset();

    fS8evt->event = iEvent.id().event();
    fS8evt->run = iEvent.id().run();
    fS8evt->evt_weight = weight ;
    //fS8evt->njets = recoJets.size();
    //fS8evt->nelectrons = recoElectrons.size();
    int total_nelectrons = 0;
    fS8evt->ngenjets = genJets.size();
    //fS8evt->nvertices = recoPV.size();

    CaloJetCollection::const_iterator jet;
    CaloJetCollection::const_iterator awayjet;
    //GenJetCollection::const_iterator genjet;
    reco::GsfElectronCollection::const_iterator electron;

    //reco::JetTagCollection::iterator btagite;
    //std::vector< BTagLeptonEvent > LeptonEvtVector;

    ///////////////////////////////
    // begin loop over jets
    //////////////////////////////
    TLorentzVector p4Jet;
    TLorentzVector p4ElecJet;
    TLorentzVector p4OppJet;
    TLorentzVector p4Electron;
    int ijet = 0;

	if (fdebug) std::cout << " begin loop over jets" << std::endl;

    for ( jet = recoJets.begin(); jet != recoJets.end(); ++jet )
    {

        // get jet corrections
        double jetcorrection = 1.;
        if (useJetCorr_ == true)
        {
            jetcorrection =  acorrector->correction(*jet);
        }
        // Jet quality cuts
        if ( (jet->pt() * jetcorrection ) <= MinJetPt_ || std::abs( jet->eta() ) >= MaxJetEta_ || jet->emEnergyFraction() < MinEmFraction_ || jet->emEnergyFraction() > MaxEmFraction_ ) continue;

        if ( !NjetsCut )
        {
            histcounterf->Fill("jets", 2 );
            NjetsCut = true;
        }

        // get MC flavor of jet
        //int JetFlavor = jetFlavourIdentifier_.identifyBasedOnPartons(*jet).flavour();
        int JetFlavor = abs(getMatchedParton(*jet).getFlavour());
        p4Jet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
        p4Jet = jetcorrection * p4Jet;
        int hasLepton = 0;
        int tmptotelectron = 0;
        BTagLeptonEvent leptonEvent;

        /////////////////////////////////
        // begin loop over electrons
        ////////////////////////////////
        double elec_highest_pt = 0;
        double ptrel = 0;
        double ptreltmp = 0;

        for ( electron = recoElectrons.begin(); electron != recoElectrons.end(); ++electron)
        {
            //      if (electron->track().isNull() || electron->standAloneElectron().isNull() || electron->combinedElectron().isNull()) continue;
            //if (electron->isGlobalElectron() == false) continue;
            //      std::cout <<" DEREF " <<electron->track().isValid()<<std::endl;
            //      std::cout <<" is the ref valid??? "<<electron->track()<<std::endl;
            GsfTrackRef mt = electron->gsfTrack();
            GsfTrack electronTrk = *(mt.get());
            //TrackingParticleRef TrueHitsTrk;
            int nhit = electronTrk.numberOfValidHits();//electronTrk.recHitsSize();

            // electron cuts
            double normChi2 = electronTrk.normalizedChi2();//(*(electron->combinedElectron())).chi2() / (*(electron->combinedElectron())).ndof();// use global fit
            bool fails = nhit <= MinElectronNHits_;
            fails = fails || electron->pt() <= MinElectronPt_;
            fails = fails || normChi2 >= MaxElectronChi2_;
            if ( fails ) continue;

            if ( !NelecCut )
            {
                histcounterf->Fill("N electrons", 3 );
                NelecCut = true;
            }

            // delta R(electron,jet)
            double deltaR  = ROOT::Math::VectorUtil::DeltaR(jet->p4().Vect(), electron->p4().Vect() );
            TVector3 tmpvecOrg(jet->p4().Vect().X(),jet->p4().Vect().Y(),  jet->p4().Vect().Z());
            TVector3 tmpvec;
            // use calibrated jet
            tmpvec = jetcorrection * tmpvecOrg;
            // find pTrel
            //TVector3 leptonvec(electronTrk.momentum().X(), electronTrk.momentum().Y(),electronTrk.momentum().Z());
            TVector3 leptonvec(electron->px(), electron->py(),electron->pz());
            //tmpvec += leptonvec;
            ptreltmp = leptonvec.Perp(tmpvec);
            double ratio = leptonvec.Mag() / tmpvecOrg.Mag();
            // electron in jet cuts
            if ( (deltaR >= MinDeltaR_ ) || (ptreltmp <= MinPtRel_ ) || (ratio >= MaxRatio_) ) continue;
            // now we have a good electron in a jet

            if ( !Nelec_in_jetCut )
            {
                histcounterf->Fill("N electron-in-jet", 4 );
                Nelec_in_jetCut = true;
            }

            total_nelectrons++;
            hasLepton = 1;
            tmptotelectron++;
            // pick the leading electron inside the jet
            if ( electron->pt() > elec_highest_pt )
            {
                elec_highest_pt = electron->pt();
                p4Electron.SetPtEtaPhiE(electron->pt(),
                                    electron->eta(),
                                    electron->phi(),
                                    electron->energy());
                // recalculate pTrel
                tmpvec = jetcorrection * tmpvecOrg;
                leptonvec.SetXYZ(electron->px(), electron->py(),electron->pz());
                //tmpvec += leptonvec;
                ptrel = leptonvec.Perp(tmpvec);  // maximum
            }

            // collect electron data
            leptonEvent.pdgid.push_back( 11 ); // electron only for the moment
            //  std::cout << "Electron energy " << electron->energy() <<std::endl;
            leptonEvent.e.push_back( electron->energy());
            leptonEvent.pt.push_back( electronTrk.pt());
            leptonEvent.eta.push_back( electronTrk.eta());
            leptonEvent.phi.push_back( electronTrk.phi());
            leptonEvent.charge.push_back( electronTrk.charge());

            leptonEvent.trkchi2.push_back( electronTrk.chi2());
            leptonEvent.trkndof.push_back( electronTrk.ndof());
            leptonEvent.chi2.push_back(    electronTrk.chi2() );
            leptonEvent.ndof.push_back(    electronTrk.ndof() );
            //leptonEvent.SArechits.push_back(  electronSA.numberOfValidHits() );
            leptonEvent.trkrechits.push_back( electronTrk.numberOfValidHits() );
            leptonEvent.d0.push_back(         electronTrk.d0());
            leptonEvent.d0sigma.push_back(    electronTrk.d0Error());

            leptonEvent.jet_ptrel.push_back( ptreltmp);
            leptonEvent.jet_deltaR.push_back( deltaR);


            // find a sim track
            if (flavourMatchOptionf == "genParticle")
            {

                for (size_t i = 0; i < genParticles->size(); ++ i)
                {
                    const GenParticle & p = (*genParticles)[i];

                    if ( abs(p.pdgId()) == 11 && p.status()==2 )
                    {

                        leptonEvent.mc_pt.push_back(           p.pt() );
                        leptonEvent.mc_phi.push_back(          p.phi() );
                        leptonEvent.mc_eta.push_back(          p.eta() );
                        leptonEvent.mc_e.push_back(            p.energy() );
                        leptonEvent.mc_charge.push_back(       p.charge() );
                        leptonEvent.mc_pdgid.push_back(        p.pdgId() );
                        leptonEvent.mc_mother_pdgid.push_back( p.mother()->pdgId() );

                    }

                }

            }

        } //close loop over electrons


        if ( hasLepton == 1 )
        {
            p4ElecJet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
            p4ElecJet = jetcorrection * p4ElecJet;
        }


        /////////////////////////////
        // find away jet
        ////////////////////////////
        bool AwayTaggedJet = false;
        bool AwayElectronJet = false;
        TLorentzVector p4AwayElectron;

        for ( awayjet = recoJets.begin(); awayjet != recoJets.end(); ++awayjet )
        {
            if ( hasLepton == 0 ) continue;

            TLorentzVector p4AwayJet;
            p4AwayJet.SetPtEtaPhiE(awayjet->pt(), awayjet->eta(), awayjet->phi(), awayjet->energy() );
            double jetcorrectionAway_ = 1.;
            if (useJetCorr_ == true)
            {
                jetcorrectionAway_ =   acorrector->correction(*awayjet);
            }
            p4AwayJet = p4AwayJet *jetcorrectionAway_;

            // Jet quality cuts
            if ( (awayjet->pt() * jetcorrectionAway_ )  <= MinJetPt_ || std::abs( awayjet->eta() ) >= MaxJetEta_ || awayjet->emEnergyFraction() < MinEmFraction_ || awayjet->emEnergyFraction() > MaxEmFraction_) continue;

            // skip electron in jet
            if ( p4AwayJet == p4ElecJet ) continue;

            // now we have an away jet

            // find an away tagged jet
			if (fdebug) std::cout << " find an away tagged jet" << std::endl;

            if ( !AwayTaggedJet )
            {

                std::map< std::string, bool > aBmap = PFTools::GetBTaggingMap(wp, *awayjet, iEvent);
                for (std::map<std::string,bool>::const_iterator imap = aBmap.begin(); imap != aBmap.end(); ++imap )
                {
                    if ( imap->first == fAwayJetTagger && imap->second )
                    {

                        AwayTaggedJet = true;

                        if ( ! Nelec_in_jet_away_taggedCut )
                        {
                            histcounterf->Fill("N elec-in-jet-away-tagged-jet", 5 );
                            Nelec_in_jet_away_taggedCut = true;
                        }
                    }

                }

            }



            // find an away electron in jet
            if ( !AwayElectronJet )
            {
                elec_highest_pt = 0;
                for ( electron = recoElectrons.begin(); electron != recoElectrons.end(); ++electron)
                {
                    //	  if (electron->track().isNull() || electron->standAloneElectron().isNull()|| electron->combinedElectron().isNull()) continue;
                    //if (electron->isGlobalElectron() == false) continue;
                    GsfTrack electronTrk = *(electron->gsfTrack().get());
                    //TrackingParticleRef TrueHitsTrk;
                    //Track electronSA = *electron->standAloneElectron();
                    int nhit = electronTrk.numberOfValidHits();

                    // electron cuts
                    double normChi2 = electronTrk.normalizedChi2();
                    if ( (nhit <= MinElectronNHits_ ) || (electron->pt()<= MinElectronPt_) || (normChi2 >= MaxElectronChi2_ ) ) continue;

                    // delta R(electron,jet)
                    double deltaR  = ROOT::Math::VectorUtil::DeltaR(awayjet->p4().Vect(), electron->p4().Vect() );
                    TVector3 tmpvecOrg(awayjet->p4().Vect().X(),awayjet->p4().Vect().Y(),  awayjet->p4().Vect().Z());
                    TVector3 tmpvec;
                    // use calibrated jet
                    tmpvec = jetcorrection * tmpvecOrg;
                    // find pTrel
                    //TVector3 leptonvec(electronTrk.momentum().X(), electronTrk.momentum().Y(),electronTrk.momentum().Z());
                    TVector3 leptonvec(electron->px(), electron->py(),electron->pz());
                    //tmpvec += leptonvec;
                    double awayptrel = leptonvec.Perp(tmpvec);
                    double awayratio = leptonvec.Mag() / tmpvecOrg.Mag();
                    // electron in jet cuts
                    if ( (deltaR >= MinDeltaR_ ) || (awayptrel <= MinPtRel_ ) || (awayratio >= MaxRatio_) ) continue;

                    // now we have a good electron in a jet
                    AwayElectronJet = true;
                    // pick the leading electron inside the jet
                    if ( electron->pt() > elec_highest_pt )
                    {
                        elec_highest_pt = electron->pt();
                        p4AwayElectron.SetPtEtaPhiE(electron->pt(),
                                                electron->eta(),
                                                electron->phi(),
                                                electron->energy());
                        // recalculate pTrel
                        tmpvec = jetcorrection * tmpvecOrg;
                        leptonvec.SetXYZ(electron->px(), electron->py(),electron->pz());
                        //tmpvec += leptonvec;
                        awayptrel = leptonvec.Perp(tmpvec);
                    }
                }
            }

        } // close away jet loop
		if (fdebug) std::cout << " get b-tagging to fill efficiency and performance plots" << std::endl;
        std::map<std::string, bool> thebtaggingmap = PFTools::GetBTaggingMap(wp, *jet, iEvent, ptrel);
        FillEff(p4Jet, JetFlavor, thebtaggingmap, weight );
		if (fdebug) std::cout << " done efficiency plots" << std::endl;
        FillPerformance(*jet, JetFlavor, iEvent );
		if (fdebug) std::cout << " done performance plots" << std::endl;

        if ( hasLepton == 1 )
        {
            p4ElecJet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
            p4ElecJet = jetcorrection * p4ElecJet;

            FillPtrel(ptrel, JetFlavor, thebtaggingmap, weight );
            FillHistos("n",p4ElecJet, ptrel, JetFlavor, thebtaggingmap, weight );

            if (AwayTaggedJet) FillHistos("p",p4ElecJet, ptrel, JetFlavor, thebtaggingmap, weight);
            //if (AwayElectronJet) FillHistos("q",p4ElecJet, ptrel, JetFlavor, this->GetBTaggingMap(*jet,jetTags_testManyByType));

        }

        fS8evt->lepton.push_back( leptonEvent );
        fS8evt->jet_hasLepton.push_back( hasLepton );
        fS8evt->jet_flavour.push_back(JetFlavor);
        fS8evt->jet_e.push_back(jet->energy());
        fS8evt->jet_pt.push_back(jet->pt());
        fS8evt->jet_eta.push_back(jet->eta());
        fS8evt->jet_phi.push_back(jet->phi());
        fS8evt->jet_et.push_back(jet->et());
        // get jet correction
        fS8evt->jetcorrection.push_back( jetcorrection );

        // find generated jet
        reco::GenJet genjet = PFTools::GetGenJet(*jet, genJets);
        if ( genjet.p() > 0 )
        {
            //fS8evt->genjet_p.push_back(genjet.p());
            fS8evt->genjet_pt.push_back(genjet.pt());
            fS8evt->genjet_eta.push_back(genjet.eta());
            fS8evt->genjet_phi.push_back(genjet.phi());
            fS8evt->genjet_e.push_back(genjet.energy());
            //fS8evt->genjet_vx.push_back(genjet.vx());
            //fS8evt->genjet_vy.push_back(genjet.vy());
            //fS8evt->genjet_vz.push_back(genjet.vz());
        }
        else
        {
            //fS8evt->genjet_p.push_back(-100000);
            fS8evt->genjet_pt.push_back(-100000);
            fS8evt->genjet_eta.push_back(-100000);
            fS8evt->genjet_phi.push_back(-100000);
            fS8evt->genjet_e.push_back(-100000);
            //fS8evt->genjet_vx.push_back(-100000);
            //fS8evt->genjet_vy.push_back(-100000);
            //fS8evt->genjet_vz.push_back(-100000);
        }

        //*************************************
        // TrackEvents
        // Running only on reco/aod samples
        // TrackCategories only in reco samples
        //*************************************

        if ( bTagTrackEventIPTagInfos.isValid() )
        {
            // Associate the jet and jettag
            // TODO : we need to use a more standard of jet matching
            int jetIndex = PFTools::TaggedJet(*jet, bTagTrackEventIPTagInfos);
            if (jetIndex < 0) continue;

            // Get a vector of reference to the selected tracks in each jet
            TrackRefVector tracks( (*bTagTrackEventIPTagInfos)[jetIndex].selectedTracks() );

            std::vector<reco::btag::TrackIPData> const & ipdata = (*tagInfo)[jetIndex].impactParameterData();

            // Create a new BTagTrackEvent
            BTagTrackEvent trackEvent;

            for ( size_t index=0; index < tracks.size(); index++ )
            {
                TrackRef track = tracks[index];

                // collect reco track data (including their categories)
                trackEvent.pt.push_back( track->pt() );
                trackEvent.eta.push_back( track->eta() );
                trackEvent.phi.push_back( track->phi() );
                trackEvent.charge.push_back( track->charge() );
                trackEvent.trkchi2.push_back( track->chi2() );
                trackEvent.trkndof.push_back( track->ndof() );
                trackEvent.trkrechits.push_back( track->numberOfValidHits() );
                trackEvent.d0.push_back( track->d0() );
                trackEvent.d0sigma.push_back( track->d0Error() );

                // Get the IP information.
                trackEvent.ip2d.push_back( ipdata[index].ip2d.value() );
                trackEvent.ip2dSigma.push_back( ipdata[index].ip2d.error() );
                trackEvent.ip3d.push_back( ipdata[index].ip3d.value() );
                trackEvent.ip3dSigma.push_back( ipdata[index].ip3d.error() );
                trackEvent.dta.push_back( ipdata[index].distanceToJetAxis.value() );

                // delta R(electron,jet)
                trackEvent.jet_deltaR.push_back(
                    ROOT::Math::VectorUtil::DeltaR(jet->p4().Vect(), track->momentum())
                );

                // ptrel calculation per track
                TVector3 trackvec(track->px(), track->py(), track->pz());
                TVector3 jetvec(jet->p4().Vect().X(),jet->p4().Vect().Y(),  jet->p4().Vect().Z());
                jetvec = jetcorrection * jetvec;
                trackEvent.jet_ptrel.push_back(
                    trackvec.Perp(jetvec)
                );
            }

            // Add the track event into BTagEvent.
            fS8evt->tracks.push_back( trackEvent );
        }

        // b tagging
        int ith_tagged = -1;
        //int isbtagged = 0;

        bool gotTCHE     = false;
        bool gotTCHEneg  = false;
        bool gotTCHP     = false;
        bool gotTCHPneg  = false;
        bool gotJP       = false;
        bool gotJPneg    = false;
        bool gotJPpos    = false;
        //bool gotSMT      = false;
        bool gotMTCHE     = false;
        bool gotMTCHP     = false;

        // use calibrated jet
        TVector3 tmpvec, tmpvecOrg(jet->p4().Vect().X(), jet->p4().Vect().Y(), jet->p4().Vect().Z());
        tmpvec = jetcorrection * tmpvecOrg;


        std::map<std::string, bool> mymap;

		if (fdebug) std::cout << " btag operating points" << std::endl;

        for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it != wp.end(); ++it)
        {
            //    for (size_t k=0; k<jetTags_testManyByType.size(); k++)
            //    {


            edm::Handle<reco::JetTagCollection > jetTags;
            //	  std::cout <<" Asking for "<< (*it).inputTag()<<std::endl;
            iEvent.getByLabel((*it).inputTag(),jetTags);
            std::string moduleLabel = (jetTags).provenance()->moduleLabel();
            if (mymap.find(moduleLabel) != mymap.end()) continue;
            mymap[moduleLabel] = true;

            ith_tagged = PFTools::TaggedJet(*jet,jetTags);

            if (ith_tagged == -1) continue;


            //*********************************
            // Track Counting taggers
            //*********************************

            // Get a vector of reference to the selected tracks in each jet
            TrackRefVector tracks((*tagInfo)[ith_tagged].selectedTracks());

            if ( moduleLabel == "trackCountingHighEffBJetTags" )
            {
                //			  std::vector< Measurement1D  > trackIP = (*tagInfo)[ith_tagged].impactParameters(0);
                std::vector< Measurement1D  > trackIP;

                std::vector<btag::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<btag::TrackIPData>::const_iterator itipdata = ipdata.begin();
                        itipdata != ipdata.end(); itipdata++)
                {
                    trackIP.push_back((*itipdata).ip3d );
                }


                if ((trackIP).size()>=2)
                {
                    float iptrack1 = (trackIP)[0].significance();
                    fS8evt->btag_TrkCounting_disc3D_1trk.push_back( iptrack1 );
                }
                fS8evt->btag_TrkCounting_disc3D_2trk.push_back( (*jetTags)[ith_tagged].second ); // 2nd trk, 3D
                gotTCHE = true;
            }
            else if ( moduleLabel == "trackCountingHighPurBJetTags" && tracks.size()>2)
            {

                fS8evt->btag_TrkCounting_disc3D_3trk.push_back( (*jetTags)[ith_tagged].second ); // 3rd trk, 3D
                gotTCHP = true;
            }
            else if ( moduleLabel == "modifiedtrackCountingHighEffBJetTags" )
            {
                //			  std::vector< Measurement1D  > trackIP = (*tagInfo)[ith_tagged].impactParameters(0);
                std::vector< Measurement1D  > trackIP;
                std::vector<btag::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<btag::TrackIPData>::const_iterator itipdata = ipdata.begin();
                        itipdata != ipdata.end(); itipdata++)
                {
                    trackIP.push_back((*itipdata).ip3d );
                }


                if ((trackIP).size()>=2)
                {
                    float iptrack1 = (trackIP)[0].significance();
                    fS8evt->btag_ModTrkCounting_disc3D_1trk.push_back( iptrack1 );
                }
                fS8evt->btag_ModTrkCounting_disc3D_2trk.push_back( (*jetTags)[ith_tagged].second ); // 2nd trk, 3D
                gotMTCHE = true;
            }
            else if ( moduleLabel == "modifiedtrackCountingHighPurBJetTags" )
            {
                fS8evt->btag_ModTrkCounting_disc3D_3trk.push_back( (*jetTags)[ith_tagged].second ); // 3rd trk, 3D
                gotMTCHP = true;
            }
            else if ( moduleLabel == "negativeTrackCounting2ndTrck" )
            {
                //			    std::vector< Measurement1D  > trackIP = (*tagInfo)[ith_tagged].impactParameters(0);
                std::vector< Measurement1D  > trackIP;
                std::vector<btag::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<btag::TrackIPData>::const_iterator itipdata = ipdata.begin();
                        itipdata != ipdata.end(); itipdata++)
                {
                    trackIP.push_back((*itipdata).ip3d );
                }


                if ((trackIP).size()>=2)
                {
                    float iptrack1 = (trackIP)[(trackIP).size()-1].significance();
                    fS8evt->btag_NegTag_disc3D_1trk.push_back( iptrack1 );
                }
                //std::cout << "discri neg tc " << (*jetTags)[ith_tagged].second << std::endl;
                fS8evt->btag_NegTag_disc3D_2trk.push_back( (*jetTags)[ith_tagged].second ); // 2nd trk, 3D
                gotTCHEneg = true;
            }
            else if ( moduleLabel == "negativeTrackCounting3rdTrck" )
            {
                //	 		    std::vector< Measurement1D  > trackIP = (*tagInfo)[ith_tagged].impactParameters(0);
                std::vector< Measurement1D  > trackIP;


                std::vector<btag::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<btag::TrackIPData>::const_iterator itipdata = ipdata.begin();
                        itipdata != ipdata.end(); itipdata++)
                {
                    trackIP.push_back((*itipdata).ip3d );
                }

                fS8evt->btag_NegTag_disc3D_3trk.push_back( (*jetTags)[ith_tagged].second ); // 3rd trk, 3D
                gotTCHPneg = true;
            }

            //*********************************
            // Jet Probability taggers
            //*********************************
            else if ( moduleLabel == "jetProbabilityBJetTags" )
            {

                fS8evt->btag_JetProb_disc3D.push_back( (*jetTags)[ith_tagged].second);

                gotJP = true;

                std::string moduleLabel = (jetTags).provenance()->moduleLabel();

                int NtrksInJet = (*tagInfo)[ith_tagged].tracks().size();
                fS8evt->jet_ntrks.push_back( NtrksInJet );
                //
                // get the probabilities
                //

                fS8evt->trackProbaVector_Size.push_back((*tagInfo)[ith_tagged].probabilities(0).size());
                if (StoreTrackProba_)
                {
                    int i=0;
                    std::vector< float > track_proba = (*tagInfo)[ith_tagged].probabilities(0) ;
                    std::vector< float > probabilities;
                    for (std::vector<float>::const_iterator it = track_proba.begin(); it!=track_proba.end(); ++it, i++)
                    {

                        double delta  = -2.;
                        delta = ROOT::Math::VectorUtil::DeltaR( (*(*jetTags)[ith_tagged].first).p4().Vect(), (*(*tagInfo)[ith_tagged].tracks()[i]).momentum());
                        if (delta <0.3) probabilities.push_back((*it));


                    }
                    fS8evt->jet_Tracks_Probability.push_back(probabilities);
                }

            }
            else if ( moduleLabel == "jetProbabilityJetTagsNegativeOnly" )
            {

                fS8evt->btag_negJetProb_disc3D.push_back( (*jetTags)[ith_tagged].second);

                gotJPneg = true;

            }
            else if ( moduleLabel == "jetProbabilityJetTagsPositiveOnly" )
            {

                fS8evt->btag_posJetProb_disc3D.push_back( (*jetTags)[ith_tagged].second);

                gotJPpos = true;

            }
            //*********************************
            // SoftLeptons Taggers
            //*********************************
            else if ( moduleLabel == "softElectronBJetTags" )
            {

                fS8evt->btag_SoftElectron_disc.push_back( (*jetTags)[ith_tagged].second);

                //gotSMT = true;

            }

        }


        //FillHistos("n",p4ElecJet, ptrel, JetFlavor,isbTaggedJet);

        if (!gotTCHE)    fS8evt->btag_TrkCounting_disc3D_2trk.push_back( -9999. );
        if (!gotTCHEneg) fS8evt->btag_NegTag_disc3D_2trk.push_back( -9999. );

        if (!gotTCHP)    fS8evt->btag_TrkCounting_disc3D_3trk.push_back( -9999. );
        if (!gotTCHPneg) fS8evt->btag_NegTag_disc3D_3trk.push_back( -9999. );

        if (!gotJP)      fS8evt->btag_JetProb_disc3D.push_back( -9999. );
        if (!gotJPneg)   fS8evt->btag_negJetProb_disc3D.push_back( -9999. );
        if (!gotJPpos)   fS8evt->btag_posJetProb_disc3D.push_back( -9999. );
        if (!gotMTCHE)   fS8evt->btag_ModTrkCounting_disc3D_2trk.push_back( -9999. );
        if (!gotMTCHP)   fS8evt->btag_ModTrkCounting_disc3D_2trk.push_back( -9999. );

        ijet++;
    } //end loop over reco jets

    fS8evt->njets = fS8evt->jet_pt.size();
    fS8evt->nelectrons = total_nelectrons;

// fill tree

    ftree->Fill();

    feventcounter++;

}

//define this as a plug-in
DEFINE_FWK_MODULE(PerformanceAnalyzerWithElectrons);

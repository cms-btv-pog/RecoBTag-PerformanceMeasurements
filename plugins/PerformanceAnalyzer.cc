
#include "RecoBTag/PerformanceMeasurements/interface/PerformanceAnalyzer.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// reco track and vertex
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

//#include "PhysicsTools/Utilities/interface/deltaR.h"
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

using namespace edm;
using namespace reco;


void WorkingPoint::print() const {
  std::cout <<" Working point "<<name()<<" Input Tag "<<inputTag()<<" Cut "<< cut()<<std::endl;
}



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PerformanceAnalyzer::PerformanceAnalyzer(const ParameterSet& iConfig)
{

    fbadeventscounter = 0;

    analyzer_ = "PerformanceAnalyzer"; // name of this analyzer

    //now do whatever initialization is needed
    //simG4_=iConfig.getParameter<edm::InputTag>( "simG4" );

    // open output file to store results
    outputFile_  = iConfig.getUntrackedParameter<std::string>("outputFile");
    // create or update output root file
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
    // verbose std output
    fverbose= iConfig.getUntrackedParameter<bool>("verbose", false);
    // default tree
    ftree = new TTree("summary","summary");
    ftree->AutoSave();


    fS8evt = new BTagEvent();
    ftree->Branch("s8.","BTagEvent",&fS8evt,64000,1); //system8 branch


    // get list of tracks
    recoTrackList_ = iConfig.getUntrackedParameter<std::string >("TrackCollection");
    // get list of PV
    //recoVtxList_ = iConfig.getUntrackedParameter<std::string >("PrimaryVertexCollection");

    MuonCollectionTags_ = iConfig.getParameter<std::string>("Muons");

    CaloJetCollectionTags_ = iConfig.getParameter<std::string>("Jets");

    GenJetCollectionTags_ = iConfig.getParameter<std::string>("GenJets");

    SimTrkCollectionTags_ = iConfig.getParameter<std::string>("SimTracks");

    bTagTrackEventIPTagInfos_ = iConfig.getParameter<std::string>("bTagTrackEventIPtagInfos");

    //jetFlavourIdentifier_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters"));
    //jetFlavourIdentifier2_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters2"));


    StoreTrackProba_ = iConfig.getParameter<bool>("StoreTrackProba");

    // jet cuts
    MinJetPt_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinPt");
    MaxJetEta_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MaxEta");
    MinDeltaR_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinDeltaR");
    MinPtRel_  = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinPtRel");
    // muon cuts
    MinMuonPt_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MinMuonPt");
    MaxMuonEta_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MaxMuonEta");
    MaxMuonChi2_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MaxMuonChi2");
    MinMuonNHits_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<int>("MinNHits");

    //
    // ue jetcorrections or not, and in cae which label to use
    //
    useJetCorr_ = iConfig.getParameter<bool>("useJetCorrections");
    if (useJetCorr_ == true){
		jetCorrLabel_ =  iConfig.getParameter<std::string>("jetCorrectionsLabel");
		std::cout<<" Use JetCorrections with Label "<<jetCorrLabel_ <<std::endl;
    }else{
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

    std::vector<edm::ParameterSet> config = iConfig.getUntrackedParameter<std::vector<edm::ParameterSet > >("bTagCutList");
    for (std::vector<edm::ParameterSet>::const_iterator it = config.begin(); it != config.end() ; ++it){
      wp.push_back(WorkingPoint(
				(*it).getUntrackedParameter<edm::InputTag> ("collection"),
				(*it).getUntrackedParameter<std::string> ("name"),
				(*it).getUntrackedParameter<double> ("cut")
				)
		   );
     		   (wp.end()-1)->print();
 
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
    topdir->mkdir("muon_in_jet");
    topdir->cd();
    histcounterf = new TH1I("histcounterf","counter",5,0,5);
    histcounterf->SetBit(TH1::kCanRebin);
    rootFile_->cd();
    
    // initialize histograms
    EffHistos     = new BTagHistograms();
    PtrelHistos   = new BTagHistograms();
    MujetHistos   = new BTagHistograms();
    AwayjetHistos = new BTagHistograms();
    TaggedMujetHistos   = new BTagHistograms();
    TaggedAwayjetHistos = new BTagHistograms();
    MujetHistos_mc   = new BTagHistograms();
    AwayjetHistos_mc = new BTagHistograms();
    TaggedMujetHistos_mc   = new BTagHistograms();
    TaggedAwayjetHistos_mc = new BTagHistograms();

    EffHistos->Init("efficiencies");
    PtrelHistos->Init("ptrel");
    MujetHistos->Init("n");
    AwayjetHistos->Init("p");
    MujetHistos_mc->Init("n","b");
    AwayjetHistos_mc->Init("p","b");
    MujetHistos_mc->Init("n","cl");
    AwayjetHistos_mc->Init("p","cl");
    MujetHistos_mc->Init("n","c");
    AwayjetHistos_mc->Init("p","c");
    MujetHistos_mc->Init("n","l");
    AwayjetHistos_mc->Init("p","l");
    for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it!=wp.end(); ++it){
      
      EffHistos->Init("efficiencies",(*it).name());
      PtrelHistos->Init("ptrel",(*it).name());
      
      TaggedMujetHistos->Init("ntag",(*it).name());
      TaggedAwayjetHistos->Init("ptag",(*it).name());
      TaggedMujetHistos_mc->Init("ntag","b",(*it).name());
      TaggedAwayjetHistos_mc->Init("ptag","b",(*it).name());
      TaggedMujetHistos_mc->Init("ntag","cl",(*it).name());
      TaggedAwayjetHistos_mc->Init("ptag","cl",(*it).name());
      TaggedMujetHistos_mc->Init("ntag","c",(*it).name());
      TaggedAwayjetHistos_mc->Init("ptag","c",(*it).name());
      TaggedMujetHistos_mc->Init("ntag","l",(*it).name());
      TaggedAwayjetHistos_mc->Init("ptag","l",(*it).name());

    }

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
    

    feventcounter = 0;

}


PerformanceAnalyzer::~PerformanceAnalyzer()
{
    rootFile_->cd();
    ftree->Write();

    topdir->cd();
    topdir->cd("MCtruth");
    EffHistos->Save();
    if (fWritePerformancePlots)
    {

        fperformanceTC2trk.Eval();
        fperformanceTC3trk.Eval();
        fperformanceMTC2trk.Eval();
        fperformanceMTC3trk.Eval();
        fperformanceTP.Eval();
        fperformanceJBP.Eval();
		fperformanceSMT.Eval();
		fperformanceSSV.Eval();
		fperformanceCSV.Eval();

        TGraphErrors *gTC2_b = new TGraphErrors(fperformanceTC2trk.GetN(),
                                                fperformanceTC2trk.GetArray("b").GetArray(),fperformanceTC2trk.GetArray("b").GetArray(),
                                                fperformanceTC2trk.GetArray("bErr").GetArray(),fperformanceTC2trk.GetArray("bErr").GetArray());

        TGraphErrors *gTC2_c = new TGraphErrors(fperformanceTC2trk.GetN(),
                                                fperformanceTC2trk.GetArray("b").GetArray(),fperformanceTC2trk.GetArray("c").GetArray(),
                                                fperformanceTC2trk.GetArray("bErr").GetArray(),fperformanceTC2trk.GetArray("cErr").GetArray());

        TGraphErrors *gTC2_udsg = new TGraphErrors(fperformanceTC2trk.GetN(),
												   fperformanceTC2trk.GetArray("b").GetArray(),fperformanceTC2trk.GetArray("udsg").GetArray(),
												   fperformanceTC2trk.GetArray("bErr").GetArray(),fperformanceTC2trk.GetArray("udsgErr").GetArray());
        TGraph *dTC2_udsg = new TGraph(fperformanceTC2trk.GetN(),
                                       fperformanceTC2trk.GetArray("udsg").GetArray(),
                                       fperformanceTC2trk.GetArray("discriminator").GetArray());

        TGraphErrors *gTC3_b = new TGraphErrors(fperformanceTC3trk.GetN(),
                                                fperformanceTC3trk.GetArray("b").GetArray(),fperformanceTC3trk.GetArray("b").GetArray(),
                                                fperformanceTC3trk.GetArray("bErr").GetArray(),fperformanceTC3trk.GetArray("bErr").GetArray());

        TGraphErrors *gTC3_c = new TGraphErrors(fperformanceTC3trk.GetN(),
                                                fperformanceTC3trk.GetArray("b").GetArray(),fperformanceTC3trk.GetArray("c").GetArray(),
                                                fperformanceTC3trk.GetArray("bErr").GetArray(),fperformanceTC3trk.GetArray("cErr").GetArray());

        TGraphErrors *gTC3_udsg = new TGraphErrors(fperformanceTC3trk.GetN(),
												   fperformanceTC3trk.GetArray("b").GetArray(),fperformanceTC3trk.GetArray("udsg").GetArray(),
												   fperformanceTC3trk.GetArray("bErr").GetArray(),fperformanceTC3trk.GetArray("udsgErr").GetArray());
        TGraph *dTC3_udsg = new TGraph(fperformanceTC3trk.GetN(),
                                       fperformanceTC3trk.GetArray("udsg").GetArray(),
                                       fperformanceTC3trk.GetArray("discriminator").GetArray());

        TGraphErrors *gTP_b = new TGraphErrors(fperformanceTP.GetN(),
                                               fperformanceTP.GetArray("b").GetArray(),fperformanceTP.GetArray("b").GetArray(),
                                               fperformanceTP.GetArray("bErr").GetArray(),fperformanceTP.GetArray("bErr").GetArray());

        TGraphErrors *gTP_c = new TGraphErrors(fperformanceTP.GetN(),
                                               fperformanceTP.GetArray("b").GetArray(),fperformanceTP.GetArray("c").GetArray(),
                                               fperformanceTP.GetArray("bErr").GetArray(),fperformanceTP.GetArray("cErr").GetArray());

        TGraphErrors *gTP_udsg = new TGraphErrors(fperformanceTP.GetN(),
												  fperformanceTP.GetArray("b").GetArray(),fperformanceTP.GetArray("udsg").GetArray(),
												  fperformanceTP.GetArray("bErr").GetArray(),fperformanceTP.GetArray("udsgErr").GetArray());
        TGraph *dTP_udsg = new TGraph(fperformanceTP.GetN(),
                                      fperformanceTP.GetArray("udsg").GetArray(),
                                      fperformanceTP.GetArray("discriminator").GetArray());

		TGraphErrors *gJBP_b = new TGraphErrors(fperformanceJBP.GetN(),
												fperformanceJBP.GetArray("b").GetArray(),fperformanceJBP.GetArray("b").GetArray(),
												fperformanceJBP.GetArray("bErr").GetArray(),fperformanceJBP.GetArray("bErr").GetArray());

        TGraphErrors *gJBP_c = new TGraphErrors(fperformanceJBP.GetN(),
												fperformanceJBP.GetArray("b").GetArray(),fperformanceJBP.GetArray("c").GetArray(),
												fperformanceJBP.GetArray("bErr").GetArray(),fperformanceJBP.GetArray("cErr").GetArray());

        TGraphErrors *gJBP_udsg = new TGraphErrors(fperformanceJBP.GetN(),
												   fperformanceJBP.GetArray("b").GetArray(),fperformanceJBP.GetArray("udsg").GetArray(),
												   fperformanceJBP.GetArray("bErr").GetArray(),fperformanceJBP.GetArray("udsgErr").GetArray());
        TGraph *dJBP_udsg = new TGraph(fperformanceJBP.GetN(),
									   fperformanceJBP.GetArray("udsg").GetArray(),
									   fperformanceJBP.GetArray("discriminator").GetArray());


		TGraphErrors *gSMT_b = new TGraphErrors(fperformanceSMT.GetN(),
												fperformanceSMT.GetArray("b").GetArray(),fperformanceSMT.GetArray("b").GetArray(),
												fperformanceSMT.GetArray("bErr").GetArray(),fperformanceSMT.GetArray("bErr").GetArray());

        TGraphErrors *gSMT_c = new TGraphErrors(fperformanceSMT.GetN(),
												fperformanceSMT.GetArray("b").GetArray(),fperformanceSMT.GetArray("c").GetArray(),
												fperformanceSMT.GetArray("bErr").GetArray(),fperformanceSMT.GetArray("cErr").GetArray());

        TGraphErrors *gSMT_udsg = new TGraphErrors(fperformanceSMT.GetN(),
												   fperformanceSMT.GetArray("b").GetArray(),fperformanceSMT.GetArray("udsg").GetArray(),
												   fperformanceSMT.GetArray("bErr").GetArray(),fperformanceSMT.GetArray("udsgErr").GetArray());
        TGraph *dSMT_udsg = new TGraph(fperformanceSMT.GetN(),
									   fperformanceSMT.GetArray("udsg").GetArray(),
									   fperformanceSMT.GetArray("discriminator").GetArray());

		TGraphErrors *gSSV_b = new TGraphErrors(fperformanceSSV.GetN(),
												fperformanceSSV.GetArray("b").GetArray(),fperformanceSSV.GetArray("b").GetArray(),
												fperformanceSSV.GetArray("bErr").GetArray(),fperformanceSSV.GetArray("bErr").GetArray());

        TGraphErrors *gSSV_c = new TGraphErrors(fperformanceSSV.GetN(),
												fperformanceSSV.GetArray("b").GetArray(),fperformanceSSV.GetArray("c").GetArray(),
												fperformanceSSV.GetArray("bErr").GetArray(),fperformanceSSV.GetArray("cErr").GetArray());

        TGraphErrors *gSSV_udsg = new TGraphErrors(fperformanceSSV.GetN(),
												   fperformanceSSV.GetArray("b").GetArray(),fperformanceSSV.GetArray("udsg").GetArray(),
												   fperformanceSSV.GetArray("bErr").GetArray(),fperformanceSSV.GetArray("udsgErr").GetArray());
        TGraph *dSSV_udsg = new TGraph(fperformanceSSV.GetN(),
									   fperformanceSSV.GetArray("udsg").GetArray(),
									   fperformanceSSV.GetArray("discriminator").GetArray());	

		TGraphErrors *gCSV_b = new TGraphErrors(fperformanceCSV.GetN(),
												fperformanceCSV.GetArray("b").GetArray(),fperformanceCSV.GetArray("b").GetArray(),
												fperformanceCSV.GetArray("bErr").GetArray(),fperformanceCSV.GetArray("bErr").GetArray());

        TGraphErrors *gCSV_c = new TGraphErrors(fperformanceCSV.GetN(),
												fperformanceCSV.GetArray("b").GetArray(),fperformanceCSV.GetArray("c").GetArray(),
												fperformanceCSV.GetArray("bErr").GetArray(),fperformanceCSV.GetArray("cErr").GetArray());

        TGraphErrors *gCSV_udsg = new TGraphErrors(fperformanceCSV.GetN(),
												   fperformanceCSV.GetArray("b").GetArray(),fperformanceCSV.GetArray("udsg").GetArray(),
												   fperformanceCSV.GetArray("bErr").GetArray(),fperformanceCSV.GetArray("udsgErr").GetArray());
        TGraph *dCSV_udsg = new TGraph(fperformanceCSV.GetN(),
									   fperformanceCSV.GetArray("udsg").GetArray(),
									   fperformanceCSV.GetArray("discriminator").GetArray());




        TGraphErrors *gMTC2_b = new TGraphErrors(fperformanceMTC2trk.GetN(),
												 fperformanceMTC2trk.GetArray("b").GetArray(),fperformanceMTC2trk.GetArray("b").GetArray(),
												 fperformanceMTC2trk.GetArray("bErr").GetArray(),fperformanceMTC2trk.GetArray("bErr").GetArray());

        TGraphErrors *gMTC2_c = new TGraphErrors(fperformanceMTC2trk.GetN(),
												 fperformanceMTC2trk.GetArray("b").GetArray(),fperformanceMTC2trk.GetArray("c").GetArray(),
												 fperformanceMTC2trk.GetArray("bErr").GetArray(),fperformanceMTC2trk.GetArray("cErr").GetArray());

        TGraphErrors *gMTC2_udsg = new TGraphErrors(fperformanceMTC2trk.GetN(),
													fperformanceMTC2trk.GetArray("b").GetArray(),fperformanceMTC2trk.GetArray("udsg").GetArray(),
													fperformanceMTC2trk.GetArray("bErr").GetArray(),fperformanceMTC2trk.GetArray("udsgErr").GetArray());
        TGraph *dMTC2_udsg = new TGraph(fperformanceMTC2trk.GetN(),
                                        fperformanceMTC2trk.GetArray("udsg").GetArray(),
                                        fperformanceMTC2trk.GetArray("discriminator").GetArray());

        TGraphErrors *gMTC3_b = new TGraphErrors(fperformanceMTC3trk.GetN(),
												 fperformanceMTC3trk.GetArray("b").GetArray(),fperformanceMTC3trk.GetArray("b").GetArray(),
												 fperformanceMTC3trk.GetArray("bErr").GetArray(),fperformanceMTC3trk.GetArray("bErr").GetArray());

        TGraphErrors *gMTC3_c = new TGraphErrors(fperformanceMTC3trk.GetN(),
												 fperformanceMTC3trk.GetArray("b").GetArray(),fperformanceMTC3trk.GetArray("c").GetArray(),
												 fperformanceMTC3trk.GetArray("bErr").GetArray(),fperformanceMTC3trk.GetArray("cErr").GetArray());

        TGraphErrors *gMTC3_udsg = new TGraphErrors(fperformanceMTC3trk.GetN(),
													fperformanceMTC3trk.GetArray("b").GetArray(),fperformanceMTC3trk.GetArray("udsg").GetArray(),
													fperformanceMTC3trk.GetArray("bErr").GetArray(),fperformanceMTC3trk.GetArray("udsgErr").GetArray());
        TGraph *dMTC3_udsg = new TGraph(fperformanceMTC3trk.GetN(),
                                        fperformanceMTC3trk.GetArray("udsg").GetArray(),
                                        fperformanceMTC3trk.GetArray("discriminator").GetArray());


        gTC2_b->SetName("gTC2_b");
        gTC2_c->SetName("gTC2_c");
        gTC2_udsg->SetName("gTC2_udsg");
        gTC3_b->SetName("gTC3_b");
        gTC3_c->SetName("gTC3_c");
        gTC3_udsg->SetName("gTC3_udsg");
        gTP_b->SetName("gTP_b");
        gTP_c->SetName("gTP_c");
        gTP_udsg->SetName("gTP_udsg");
		gJBP_b->SetName("gJBP_b");
        gJBP_c->SetName("gJBP_c");
        gJBP_udsg->SetName("gJBP_udsg");
		gSMT_b->SetName("gSMT_b");
        gSMT_c->SetName("gSMT_c");
        gSMT_udsg->SetName("gSMT_udsg");
		gSSV_b->SetName("gSSV_b");
        gSSV_c->SetName("gSSV_c");
        gSSV_udsg->SetName("gSSV_udsg");
 		gCSV_b->SetName("gCSV_b");
        gCSV_c->SetName("gCSV_c");
        gCSV_udsg->SetName("gCSV_udsg");
        gMTC2_b->SetName("gMTC2_b");
        gMTC2_c->SetName("gMTC2_c");
        gMTC2_udsg->SetName("gMTC2_udsg");
        gMTC3_b->SetName("gMTC3_b");
        gMTC3_c->SetName("gMTC3_c");
        gMTC3_udsg->SetName("gMTC3_udsg");


        dTC2_udsg->SetName("discTC2_udsg");
        dTC3_udsg->SetName("discTC3_udsg");
        dTP_udsg->SetName("discTP_udsg");
	dJBP_udsg->SetName("discJBP_udsg");
	dSMT_udsg->SetName("discSMT_udsg");
	dSSV_udsg->SetName("discSSV_udsg");
	dCSV_udsg->SetName("discCSV_udsg");
        dMTC2_udsg->SetName("discMTC2_udsg");
        dMTC3_udsg->SetName("discMTC3_udsg");
        dTC2_udsg->SetTitle("TC2trk discriminator vs udsg-mistagging");
        dTC3_udsg->SetTitle("TC3trk discriminator vs udsg-mistagging");
        dTP_udsg->SetTitle("TP discriminator vs udsg-mistagging");
	dJBP_udsg->SetTitle("JBP discriminator vs udsg-mistagging");
	dSMT_udsg->SetTitle("SMT discriminator vs udsg-mistagging");
	dSSV_udsg->SetTitle("SSV discriminator vs udsg-mistagging");
	dCSV_udsg->SetTitle("CSV discriminator vs udsg-mistagging");
        dMTC2_udsg->SetTitle("MTC2trk discriminator vs udsg-mistagging");
        dMTC3_udsg->SetTitle("MTC3trk discriminator vs udsg-mistagging");

        gTC2_b->SetTitle("Jet b-efficiency");
        gTC2_c->SetTitle("Jet c-mistagging");
        gTC2_udsg->SetTitle("Jet udsg-mistagging");
        gTC3_b->SetTitle("Jet b-efficiency");
        gTC3_c->SetTitle("Jet c-mistagging");
        gTC3_udsg->SetTitle("Jet udsg-mistagging");
        gTP_b->SetTitle("Jet b-efficiency");
        gTP_c->SetTitle("Jet c-mistagging");
        gTP_udsg->SetTitle("Jet udsg-mistagging");
	gJBP_b->SetTitle("Jet b-efficiency");
        gJBP_c->SetTitle("Jet c-mistagging");
        gJBP_udsg->SetTitle("Jet udsg-mistagging");
	gSMT_b->SetTitle("Jet b-efficiency");
        gSMT_c->SetTitle("Jet c-mistagging");
        gSMT_udsg->SetTitle("Jet udsg-mistagging");
	gSSV_b->SetTitle("Jet b-efficiency");
        gSSV_c->SetTitle("Jet c-mistagging");
        gSSV_udsg->SetTitle("Jet udsg-mistagging");
	gCSV_b->SetTitle("Jet b-efficiency");
        gCSV_c->SetTitle("Jet c-mistagging");
        gCSV_udsg->SetTitle("Jet udsg-mistagging");
        gMTC2_b->SetTitle("Jet b-efficiency");
        gMTC2_c->SetTitle("Jet c-mistagging");
        gMTC2_udsg->SetTitle("Jet udsg-mistagging");
        gMTC3_b->SetTitle("Jet b-efficiency");
        gMTC3_c->SetTitle("Jet c-mistagging");
        gMTC3_udsg->SetTitle("Jet udsg-mistagging");

        gTC2_b->Write();
        gTC2_c->Write();
        gTC2_udsg->Write();
        gTC3_b->Write();
        gTC3_c->Write();
        gTC3_udsg->Write();
        gTP_b->Write();
        gTP_c->Write();
        gTP_udsg->Write();
	gJBP_b->Write();
        gJBP_c->Write();
        gJBP_udsg->Write();
	gSMT_b->Write();
        gSMT_c->Write();
        gSMT_udsg->Write();
	gSSV_b->Write();
        gSSV_c->Write();
        gSSV_udsg->Write();
	gCSV_b->Write();
        gCSV_c->Write();
        gCSV_udsg->Write();
        gMTC2_b->Write();
        gMTC2_c->Write();
        gMTC2_udsg->Write();
        gMTC3_b->Write();
        gMTC3_c->Write();
        gMTC3_udsg->Write();


        dTC2_udsg->Write();
        dTC3_udsg->Write();
        dTP_udsg->Write();
        dJBP_udsg->Write();
        dSMT_udsg->Write();
	dSSV_udsg->Write();
	dCSV_udsg->Write();
	dMTC2_udsg->Write();
        dMTC3_udsg->Write();

    }

    topdir->cd("ptrel");
    PtrelHistos->Save();
    topdir->cd("muon_in_jet");
    MujetHistos->Save();
    AwayjetHistos->Save();
    TaggedMujetHistos->Save();
    TaggedAwayjetHistos->Save();
    topdir->cd("MCtruth");
    MujetHistos_mc->Save();
    AwayjetHistos_mc->Save();
    TaggedMujetHistos_mc->Save();
    TaggedAwayjetHistos_mc->Save();

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
void PerformanceAnalyzer::beginJob(edm::EventSetup const& iSetup)
{
    std::cout << analyzer_ << " begin Job" << std::endl;
    rootFile_->cd();
}


void PerformanceAnalyzer::endJob()
{

    std::cout << analyzer_ << " Total events processed: " << feventcounter << std::endl;

    rootFile_->cd();

}

reco::GenJet PerformanceAnalyzer::GetGenJet(reco::CaloJet calojet, reco::GenJetCollection genJetColl)
{

    reco::GenJet matchedJet;

    double predelta = 99999.;
    for (GenJetCollection::const_iterator genjet = genJetColl.begin(); genjet != genJetColl.end(); ++genjet)
    {

        double delta  = ROOT::Math::VectorUtil::DeltaR(genjet->p4().Vect(), calojet.p4().Vect() );

        if ( delta < 0.2 && delta<predelta )
        {
            matchedJet = *genjet;
            predelta = delta;
        }
    }

    return matchedJet;
}


int PerformanceAnalyzer::TaggedJet(reco::CaloJet calojet,edm::Handle<reco::JetTagCollection > jetTags )
{

    double small = 1.e-5;
    int result = -1; // no tagged
    //int ith = -1;

    //std::cout << "calo jet: pz = " << calojet.pz() << " pt = " << calojet.pt() << std::endl;
    //for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
    //  edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];

    //get label and module names


    //    std::cout <<" ECCO " << jetTags.product()<< std::endl;


    for (size_t t = 0; t < jetTags->size(); ++t)
    {
        edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].first;
        if (jet_p.isNull())
        {
            //std::cout << "-----------> JetTag::jet() returned null reference" << std::endl;
            continue;
        }
        //std::cout << "[TaggedJet]  calojet pt = " << calojet.pt() << " tagged jet pt = " << jet_p->pt() << std::endl;
        if (DeltaR<reco::Candidate>()( calojet, *jet_p ) < small)
        {

            result = (int) t;

        }
    }

    return result;
}


int PerformanceAnalyzer::TaggedJet(reco::CaloJet const & caloJet, edm::Handle<std::vector<reco::TrackIPTagInfo> > const & trackIPTagInfos )
{
    double small = 1.e-5;
    int result = -1;

    for (size_t t = 0; t < trackIPTagInfos->size(); ++t)
    {
        edm::RefToBase<reco::Jet> jet_p = (*trackIPTagInfos)[t].jet();
        if (jet_p.isNull())
        {
            //std::cout << "-----------> TrackIPTagInfos::jet() returned null reference" << std::endl;
            continue;
        }
        //std::cout << "[TaggedJet]  calojet pt = " << calojet.pt() << " tagged jet pt = " << jet_p->pt() << std::endl;
        if (DeltaR<reco::Candidate>()( caloJet, * jet_p ) < small)
        {
            result = (int) t;
            break;
        }
    }

    return result;
}


SimTrack PerformanceAnalyzer::GetGenTrk(reco::Track atrack, const edm::SimTrackContainer *simTrkColl, const edm::SimVertexContainer *simVtcs)
{

    SimTrack matchedTrk;
    edm::SimVertexContainer mysimVtcs = *simVtcs;

    double predelta = 99999.;
    for (SimTrackContainer::const_iterator gentrk = simTrkColl->begin(); gentrk != simTrkColl->end(); ++gentrk)
    {
        // VEB
        // double delta  = ROOT::Math::VectorUtil::DeltaR( TVector3((*gentrk).momentum().x(),(*gentrk).momentum().y(),(*gentrk).momentum().z()) , atrack.momentum() );
        double delta = reco::deltaR( (*gentrk).momentum(), atrack.momentum() );

        int type = (*gentrk).type();
        if ( delta < 0.2 && delta<predelta && ((*gentrk).charge() == atrack.charge() ) &&
			 ( abs(type)==11 || abs(type)==13 || abs(type)==15 || abs(type)==211 || abs(type)==321 ) )
        {
            matchedTrk = *gentrk;
            predelta = delta;
            math::XYZTLorentzVectorD v = (mysimVtcs)[(*gentrk).vertIndex()].position();

            //std::cout << "gentrk: vx = " << v.x() << std::endl;
            //std::cout << "rectrk: vx = " << atrack.vx() << std::endl;

        }
    }

    return matchedTrk;
}

int
PerformanceAnalyzer::GetMotherId(const edm::SimVertexContainer *simVtxColl, const edm::SimTrackContainer *simTrkColl, SimTrack muonMC)
{

    edm::SimVertexContainer mysimVtxColl = *simVtxColl;

    edm::SimTrackContainer mysimTrkColl = *simTrkColl;

    // fill map of simtracks
    std::map<unsigned, unsigned> geantToIndex;
    for ( unsigned it=0; it< mysimTrkColl.size(); ++it )
    {
        geantToIndex[ mysimTrkColl[it].trackId() ] = it;
    }

    // The origin vertex
    int vertexId = muonMC.vertIndex();
    SimVertex vertex = mysimVtxColl[vertexId];

    // The mother track
    int motherId = -1;
    if ( vertex.parentIndex() )  // there is a parent to this vertex
    {

        // geant id of the mother
        unsigned motherGeandId =   vertex.parentIndex();
        std::map<unsigned, unsigned >::iterator association
			= geantToIndex.find( motherGeandId );
        if (association != geantToIndex.end() )
            motherId = association->second;
    }


    int motherType = motherId == -1 ? 0 : mysimTrkColl[motherId].type();

    return motherType;

}

//______________________________________________________________________________________________________________________
std::map< std::string, bool >
PerformanceAnalyzer::GetBTaggingMap(reco::CaloJet jet, const Event& event, double ptrel )
{

    std::map< std::string, bool > aMap;

    int ith_tagged = -1;


    for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it != wp.end(); ++it){
      
      //start loop over all jetTags
      
      //    for (size_t k=0; k<jetTags_testManyByType.size(); k++)
      //    {
      
      edm::Handle<reco::JetTagCollection > jetTags;
      event.getByLabel((*it).inputTag(),jetTags);
      
      
      ith_tagged = this->TaggedJet(jet,jetTags);
      
      if (ith_tagged == -1) continue;
      
      //*********************************
      // Track Counting taggers
      //*********************************
      
      
      if ((*jetTags)[ith_tagged].second > (*it).cut()) {
	aMap[(*it).name()] = true;
      }else{
	aMap[(*it).name()] = false;
      }
      
    }
    
    
    return aMap;
}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillPerformance(reco::CaloJet jet, int JetFlavor, const edm::Event& event)
{

    std::map< std::string, bool > aMap;
    int ith_tagged = -1;

    std::map<std::string, bool> mymap;

    for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it != wp.end(); ++it){
      //    for (size_t k=0; k<jetTags_testManyByType.size(); k++)
      //    {
      
      
      
      edm::Handle<reco::JetTagCollection > jetTags;
      //      std::cout <<" Asking for "<< (*it).inputTag()<<std::endl;
      event.getByLabel((*it).inputTag(),jetTags);
      std::string moduleLabel = (jetTags).provenance()->moduleLabel();
      if (mymap.find(moduleLabel) != mymap.end()) continue;
      
      mymap[moduleLabel] == true;
      
      //std::cout <<" ECCO "<< jetTags->size()<<std::endl;
      
      ith_tagged = this->TaggedJet(jet,jetTags);
      
      
      if (ith_tagged == -1) continue;
      
      
      //*********************************
      // Track Counting taggers
      //*********************************
      if ( moduleLabel == "trackCountingHighEffBJetTags" )
        {
	  
	  fperformanceTC2trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      else if ( moduleLabel == "trackCountingHighPurBJetTags" )
        {
	  
	  fperformanceTC3trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      
      //*********************************
      // Modified Track Counting taggers
      //*********************************
      if ( moduleLabel == "modifiedtrackCountingHighEffBJetTags" )
        {
	  
	  fperformanceMTC2trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      else if ( moduleLabel == "modifiedtrackCountingHighPurBJetTags" )
        {
	  
	  fperformanceMTC3trk.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      
      //*********************************
      // Jet Probability taggers
      //*********************************
      else if ( moduleLabel == "jetProbabilityBJetTags" )
        {
	  
	  fperformanceTP.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      // jetBProbabilityBJetTags
      else if ( moduleLabel == "jetBProbabilityBJetTags" )
        {
	  
	  fperformanceJBP.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      
      //*********************************
      // Soft Muon Tagger
      //*********************************
      else if ( moduleLabel == "softMuonBJetTags" )
        {
	  
	  fperformanceSMT.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      
      //*********************************
      // Secondary Vertex Tagger
      //*********************************
      else if ( moduleLabel == "simpleSecondaryVertexBJetTags" )
        {
	  
	  fperformanceSSV.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      else if ( moduleLabel == "combinedSecondaryVertexBJetTags" )
        {
	  
	  fperformanceCSV.Add( (*jetTags)[ith_tagged].second, JetFlavor );
        }
      
    }
    
    
}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillEff(TLorentzVector p4MuJet, int JetFlavor, std::map<std::string, bool> aMap, double weight)
{

    std::string flavor = "";
    for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap )
    {

        if (imap->second)
        {

            std::string tag = "_"+imap->first;
            flavor = "";

            EffHistos->Fill1d("jet_pt"+tag,p4MuJet.Pt(), weight );
            EffHistos->Fill1d("jet_eta"+tag,p4MuJet.Eta(), weight );

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
                EffHistos->Fill1d("jet_pt"+flavor+tag,p4MuJet.Pt(), weight );
                EffHistos->Fill1d("jet_eta"+flavor+tag,p4MuJet.Eta(),weight );
            }
        }
    }

    EffHistos->Fill1d("jet_pt",p4MuJet.Pt(),weight );
    EffHistos->Fill1d("jet_eta",p4MuJet.Eta(),weight );

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
        EffHistos->Fill1d("jet_pt"+flavor,p4MuJet.Pt(), weight );
        EffHistos->Fill1d("jet_eta"+flavor,p4MuJet.Eta(),weight );
    }

}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillPtrel(double ptrel, int JetFlavor, std::map<std::string, bool> aMap, double weight)
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
void PerformanceAnalyzer::FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
                                     int JetFlavor, std::map<std::string, bool> aMap, double weight)
{
    if ( type == "n")
    {
        MujetHistos->Fill2d(type+"_pT",p4MuJet.Pt(),ptrel,weight);
        MujetHistos->Fill2d(type+"_eta",TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        if ( JetFlavor == 5 )
        {
            std::string flavor = "b";
            MujetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            MujetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
        {
            std::string flavor = "cl";
            MujetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            MujetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
          if ( JetFlavor == 4 )
        {
            std::string flavor = "c";
            MujetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            MujetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
        {
            std::string flavor = "l";
            MujetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            MujetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
  }
    else if ( type == "p")
    {

        AwayjetHistos->Fill2d(type+"_pT",p4MuJet.Pt(),ptrel,weight);
        AwayjetHistos->Fill2d(type+"_eta",TMath::Abs(p4MuJet.Eta()),ptrel, weight);
        if ( JetFlavor == 5 )
        {
            std::string flavor = "b";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
        {
            std::string flavor = "cl";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
        if ( JetFlavor == 4 )
        {
            std::string flavor = "c";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
        if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
        {
            std::string flavor = "l";
            AwayjetHistos_mc->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
            AwayjetHistos_mc->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
        }
    }

    for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap )
    {

        if ( imap->second )
        {

            if ( type == "n")
            {

                TaggedMujetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4MuJet.Pt(),ptrel,weight);
                TaggedMujetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);

                if ( JetFlavor == 5 )
                {
                    std::string flavor = "b_";
                    TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
                {
                    std::string flavor = "cl_";
                    TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
                if ( JetFlavor == 4 )
                {
                    std::string flavor = "c_";
                    TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
                {
                    std::string flavor = "l_";
                    TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
            }
            else if ( type == "p")
            {


                TaggedAwayjetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4MuJet.Pt(),ptrel,weight);
                TaggedAwayjetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);

                if ( JetFlavor == 5 )
                {
                    std::string flavor = "b_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
                {
                    std::string flavor = "cl_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
                if ( JetFlavor == 4 )
                {
                    std::string flavor = "c_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
                if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
                {
                    std::string flavor = "l_";
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
                    TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
            }
        }
    }

}

reco::JetFlavour PerformanceAnalyzer::getMatchedParton(const reco::CaloJet &jet)
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

        /*    for( reco::CandMatchMap::const_iterator f  = theJetPartonMapf->begin();
			  f != theJetPartonMapf->end(); f++) {
			  const reco::Candidate *theJetInTheMatchMap = &*(f->key);
			  const reco::Candidate *theMatchedParton    = &*(f->val);
			  if(theJetInTheMatchMap->hasMasterClone ()) {
			  const reco::CaloJet* theMasterClone = dynamic_cast<const reco::CaloJet*>(theJetInTheMatchMap->masterClone().get());
			  //std::cout << " masterclone pt = " << theMasterClone->pt() << " calo jet pt = " << jet.pt() << std::endl;
			  // FIXME, compare pointers rather than values:
			  //if ( fabs( theMasterClone->pt() - jet.pt() ) < 1.e-5 ) {
			  if ( fabs(theMasterClone->phi() - jet.phi()) < 1.e-5 && fabs(theMasterClone->eta() - jet.eta())< 1.e-5 ){
			  //std::cout << " it matches! " << std::endl;
			  jetFlavour.flavour(abs(theMatchedParton->pdgId()));
			  jetFlavour.underlyingParton4Vec(theMatchedParton->p4());
			  return jetFlavour;
			  }
			  }
			  }
        */

    }

    return jetFlavour;
}


// ------------ method called to produce the data  ------------
void
PerformanceAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  
  // count
  histcounterf->Fill("Processed", 1);
  bool NjetsCut = false;
  bool NmuCut = false;
  bool Nmu_in_jetCut = false;
  bool Nmu_in_jet_away_taggedCut = false;

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
  
  //
  
  // Muons
  Handle<reco::MuonCollection> muonsColl;
  iEvent.getByLabel(MuonCollectionTags_, muonsColl);
  
  // Calo Jets
  Handle<reco::CaloJetCollection> jetsColl;
  iEvent.getByLabel(CaloJetCollectionTags_, jetsColl);
  
  // Get the bTagTrackEventIPTagInfo collection
  Handle<std::vector<reco::TrackIPTagInfo> > bTagTrackEventIPTagInfos;
  
  if ( !bTagTrackEventIPTagInfos_.empty() )
    iEvent.getByLabel(bTagTrackEventIPTagInfos_, bTagTrackEventIPTagInfos);
  
  
  // initialize jet corrector
  const JetCorrector *acorrector = 0;
  if (useJetCorr_ == true){
    acorrector = JetCorrector::getJetCorrector(jetCorrLabel_,iSetup);
  }
  
  
  Handle<reco::GenJetCollection> genjetsColl;
  iEvent.getByLabel(GenJetCollectionTags_, genjetsColl);
  
  
  const reco::CaloJetCollection recoJets =   *(jetsColl.product());
  const reco::GenJetCollection  genJets  =   *(genjetsColl.product());
  const reco::MuonCollection    recoMuons =  *(muonsColl.product());
  
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
  //fS8evt->nmuons = recoMuons.size();
  int total_nmuons = 0;
  fS8evt->ngenjets = genJets.size();
  //fS8evt->nvertices = recoPV.size();
  
  CaloJetCollection::const_iterator jet;
  CaloJetCollection::const_iterator awayjet;
  //GenJetCollection::const_iterator genjet;
  reco::MuonCollection::const_iterator muon;
  
  //reco::JetTagCollection::iterator btagite;
  //std::vector< BTagLeptonEvent > LeptonEvtVector;
  
  ///////////////////////////////
  // begin loop over jets
  //////////////////////////////
  TLorentzVector p4Jet;
  TLorentzVector p4MuJet;
  TLorentzVector p4OppJet;
  TLorentzVector p4Muon;
  int ijet = 0;
  
  for ( jet = recoJets.begin(); jet != recoJets.end(); ++jet )
    {
      
      // get jet corrections
      double jetcorrection = 1.;
      if (useJetCorr_ == true){
	jetcorrection =  acorrector->correction(*jet, iEvent, iSetup);
      }
      // Jet quality cuts
      if ( (jet->pt() * jetcorrection ) <= MinJetPt_ || std::abs( jet->eta() ) >= MaxJetEta_ ) continue;
      
      if ( !NjetsCut ) { histcounterf->Fill("jets", 2 ); NjetsCut = true; }
      
      // get MC flavor of jet
      //int JetFlavor = jetFlavourIdentifier_.identifyBasedOnPartons(*jet).flavour();
      int JetFlavor = abs(getMatchedParton(*jet).getFlavour());
      p4Jet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
      p4Jet = jetcorrection * p4Jet;
      int hasLepton = 0;
      int tmptotmuon = 0;
      BTagLeptonEvent leptonEvent;
      
      /////////////////////////////////
      // begin loop over muons
      ////////////////////////////////
      double mu_highest_pt = 0;
      double ptrel = 0;
      double ptreltmp = 0;
      
      for ( muon = recoMuons.begin(); muon != recoMuons.end(); ++muon)
        {
	  //      if (muon->track().isNull() || muon->standAloneMuon().isNull() || muon->combinedMuon().isNull()) continue;
	  if (muon->isGlobalMuon() == false) continue;
	  //      std::cout <<" DEREF " <<muon->track().isValid()<<std::endl;
	  //      std::cout <<" is the ref valid??? "<<muon->track()<<std::endl;
	  TrackRef mt = muon->track();
	  Track muonTrk = *muon->track();
	  //TrackingParticleRef TrueHitsTrk;
	  Track muonSA = *muon->standAloneMuon();
	  int nhit = muonTrk.numberOfValidHits();//muonTrk.recHitsSize();
	  // muon cuts
	  double normChi2 = (*(muon->combinedMuon())).normalizedChi2();//(*(muon->combinedMuon())).chi2() / (*(muon->combinedMuon())).ndof();// use global fit
	  if ( (nhit <= MinMuonNHits_ ) || (muon->pt()<= MinMuonPt_) || (normChi2 >= MaxMuonChi2_ ) ) continue;
	  
	  if ( !NmuCut ) { histcounterf->Fill("N muons", 3 ); NmuCut = true; } 

	  // delta R(muon,jet)
	  double deltaR  = ROOT::Math::VectorUtil::DeltaR(jet->p4().Vect(), muon->p4().Vect() );
	  TVector3 tmpvecOrg(jet->p4().Vect().X(),jet->p4().Vect().Y(),  jet->p4().Vect().Z());
	  TVector3 tmpvec;
	  // use calibrated jet
	  tmpvec = jetcorrection * tmpvecOrg;
	  // find pTrel
	  //TVector3 leptonvec(muonTrk.momentum().X(), muonTrk.momentum().Y(),muonTrk.momentum().Z());
	  TVector3 leptonvec(muon->px(), muon->py(),muon->pz());
	  tmpvec += leptonvec;
	  ptreltmp = leptonvec.Perp(tmpvec);
	  // muon in jet cuts
	  if ( (deltaR >= MinDeltaR_ ) || (ptreltmp <= MinPtRel_ ) ) continue;
	  // now we have a good muon in a jet

	  if ( !Nmu_in_jetCut ) { histcounterf->Fill("N muon-in-jet", 4 ); Nmu_in_jetCut = true; }

	  total_nmuons++;
	  hasLepton = 1;
	  tmptotmuon++;
	  // pick the leading muon inside the jet
	  if ( muon->pt() > mu_highest_pt )
            {
	      mu_highest_pt = muon->pt();
	      p4Muon.SetPtEtaPhiE(muon->pt(),
				  muon->eta(),
				  muon->phi(),
				  muon->energy());
	      // recalculate pTrel
	      tmpvec = jetcorrection * tmpvecOrg;
	      leptonvec.SetXYZ(muon->px(), muon->py(),muon->pz());
	      tmpvec += leptonvec;
	      ptrel = leptonvec.Perp(tmpvec);  // maximum
            }
	  
	  // collect muon data
	  leptonEvent.pdgid.push_back( 13 ); // muon only for the moment
	  std::cout << "Muon energy " << muon->energy() <<std::endl;
	  leptonEvent.e.push_back( muon->energy());
	  leptonEvent.pt.push_back( muonTrk.pt());
	  leptonEvent.eta.push_back( muonTrk.eta());
	  leptonEvent.phi.push_back( muonTrk.phi());
	  leptonEvent.charge.push_back( muonTrk.charge());
	  
	  leptonEvent.trkchi2.push_back( muonTrk.chi2());
	  leptonEvent.trkndof.push_back( muonTrk.ndof());
	  leptonEvent.chi2.push_back(    (*(muon->combinedMuon())).chi2() );
	  leptonEvent.ndof.push_back(    (*(muon->combinedMuon())).ndof() );
	  leptonEvent.SArechits.push_back(  muonSA.numberOfValidHits() );
	  leptonEvent.trkrechits.push_back( muonTrk.numberOfValidHits() );
	  leptonEvent.d0.push_back(         muonTrk.d0());
	  leptonEvent.d0sigma.push_back(    muonTrk.d0Error());
	  
	  leptonEvent.jet_ptrel.push_back( ptreltmp);
	  leptonEvent.jet_deltaR.push_back( deltaR);
	  
	  
	  // find a sim track
	  if (flavourMatchOptionf == "genParticle")
            {
	      
	      for (size_t i = 0; i < genParticles->size(); ++ i)
                {
		  const GenParticle & p = (*genParticles)[i];
		  
		  if ( abs(p.pdgId()) == 13 && p.status()==2 )
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
	  
        } //close loop over muons
      
      
      if ( hasLepton == 1 )
        {
	  p4MuJet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
	  p4MuJet = jetcorrection * p4MuJet;
        }
      
      
      /////////////////////////////
      // find away jet
      ////////////////////////////
      bool AwayTaggedJet = false;
      bool AwayMuonJet = false;
      TLorentzVector p4AwayMuon;
      
      for ( awayjet = recoJets.begin(); awayjet != recoJets.end(); ++awayjet )
        {
	  if ( hasLepton == 0 ) continue;
	  
	  TLorentzVector p4AwayJet;
	  p4AwayJet.SetPtEtaPhiE(awayjet->pt(), awayjet->eta(), awayjet->phi(), awayjet->energy() );
	  double jetcorrectionAway_ = 1.;
	  if (useJetCorr_ == true){
	    jetcorrectionAway_ =   acorrector->correction(*awayjet, iEvent, iSetup);
	  }
	  p4AwayJet = p4AwayJet *jetcorrectionAway_;
	  
	  // Jet quality cuts
	  if ( (awayjet->pt() * jetcorrectionAway_ )  <= MinJetPt_ || std::abs( awayjet->eta() ) >= MaxJetEta_ ) continue;
	  
	  // skip muon in jet
	  if ( p4AwayJet == p4MuJet ) continue;
	  
	  // now we have an away jet
	  
	  // find an away tagged jet
	  if ( !AwayTaggedJet )
            {
	      
	      std::map< std::string, bool > aBmap = this->GetBTaggingMap(*awayjet,iEvent);
	      for (std::map<std::string,bool>::const_iterator imap = aBmap.begin(); imap != aBmap.end(); ++imap )
                {
                    if ( imap->first == fAwayJetTagger && imap->second )
                    {

                        AwayTaggedJet = true;

			if ( ! Nmu_in_jet_away_taggedCut ) { histcounterf->Fill("N mu-in-jet-away-tagged-jet", 5 ); Nmu_in_jet_away_taggedCut = true; }
                    }

                }

            }

	  

            // find an away muon in jet
            if ( !AwayMuonJet )
            {
                mu_highest_pt = 0;
                for ( muon = recoMuons.begin(); muon != recoMuons.end(); ++muon)
                {
                    //	  if (muon->track().isNull() || muon->standAloneMuon().isNull()|| muon->combinedMuon().isNull()) continue;
                    if (muon->isGlobalMuon() == false) continue;
                    Track muonTrk = *muon->track();
                    //TrackingParticleRef TrueHitsTrk;
                    Track muonSA = *muon->standAloneMuon();
                    int nhit = muonTrk.numberOfValidHits();

                    // muon cuts
                    double normChi2 = (*(muon->combinedMuon())).normalizedChi2();
                    if ( (nhit <= MinMuonNHits_ ) || (muon->pt()<= MinMuonPt_) || (normChi2 >= MaxMuonChi2_ ) ) continue;

                    // delta R(muon,jet)
                    double deltaR  = ROOT::Math::VectorUtil::DeltaR(awayjet->p4().Vect(), muon->p4().Vect() );
                    TVector3 tmpvecOrg(awayjet->p4().Vect().X(),awayjet->p4().Vect().Y(),  awayjet->p4().Vect().Z());
                    TVector3 tmpvec;
                    // use calibrated jet
                    tmpvec = jetcorrection * tmpvecOrg;
                    // find pTrel
                    //TVector3 leptonvec(muonTrk.momentum().X(), muonTrk.momentum().Y(),muonTrk.momentum().Z());
                    TVector3 leptonvec(muon->px(), muon->py(),muon->pz());
                    tmpvec += leptonvec;
                    double awayptrel = leptonvec.Perp(tmpvec);
                    // muon in jet cuts
                    if ( (deltaR >= MinDeltaR_ ) || (awayptrel <= MinPtRel_ ) ) continue;

                    // now we have a good muon in a jet
                    AwayMuonJet = true;
                    // pick the leading muon inside the jet
                    if ( muon->pt() > mu_highest_pt )
                    {
                        mu_highest_pt = muon->pt();
                        p4AwayMuon.SetPtEtaPhiE(muon->pt(),
                                                muon->eta(),
                                                muon->phi(),
                                                muon->energy());
                        // recalculate pTrel
                        tmpvec = jetcorrection * tmpvecOrg;
                        leptonvec.SetXYZ(muon->px(), muon->py(),muon->pz());
                        tmpvec += leptonvec;
                        awayptrel = leptonvec.Perp(tmpvec);
                    }
                }
            }

        } // close away jet loop

        std::map< std::string, bool > thebtaggingmap = this->GetBTaggingMap(*jet,iEvent, ptrel);
        FillEff(p4Jet, JetFlavor, thebtaggingmap, weight );
        FillPerformance(*jet, JetFlavor, iEvent );

        if ( hasLepton == 1 )
        {
            p4MuJet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
            p4MuJet = jetcorrection * p4MuJet;

            FillPtrel(ptrel, JetFlavor, thebtaggingmap, weight );
            FillHistos("n",p4MuJet, ptrel, JetFlavor, thebtaggingmap, weight );

            if (AwayTaggedJet) FillHistos("p",p4MuJet, ptrel, JetFlavor, thebtaggingmap, weight);
            //if (AwayMuonJet) FillHistos("q",p4MuJet, ptrel, JetFlavor, this->GetBTaggingMap(*jet,jetTags_testManyByType));

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
        reco::GenJet genjet = this->GetGenJet(*jet,genJets);
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
            int jetIndex = TaggedJet(*jet, bTagTrackEventIPTagInfos);
            if (jetIndex < 0) continue;

            // Get a vector of reference to the selected tracks in each jet
            TrackRefVector tracks( (*bTagTrackEventIPTagInfos)[jetIndex].selectedTracks() );

            std::vector<reco::TrackIPTagInfo::TrackIPData> const & ipdata = (*tagInfo)[jetIndex].impactParameterData();

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

                // delta R(muon,jet)
                trackEvent.jet_deltaR.push_back(
                    ROOT::Math::VectorUtil::DeltaR(jet->p4().Vect(), track->momentum())
					);

                // ptrel calculation per track
                TVector3 trackvec(track->px(), track->py(), track->pz());
                TVector3 jetvec(jet->p4().Vect().X(),jet->p4().Vect().Y(),  jet->p4().Vect().Z());
                jetvec = jetcorrection * jetvec;
                trackEvent.jet_ptrel.push_back(
                    trackvec.Perp(jetvec + trackvec)
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
        bool gotSMT      = false;
        bool gotMTCHE     = false;
        bool gotMTCHP     = false;

        // use calibrated jet
        TVector3 tmpvec, tmpvecOrg(jet->p4().Vect().X(), jet->p4().Vect().Y(), jet->p4().Vect().Z());
        tmpvec = jetcorrection * tmpvecOrg;


	std::map<std::string, bool> mymap;
	
	for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it != wp.end(); ++it){
	  //    for (size_t k=0; k<jetTags_testManyByType.size(); k++)
	  //    {
	  

	  edm::Handle<reco::JetTagCollection > jetTags;
	  //	  std::cout <<" Asking for "<< (*it).inputTag()<<std::endl;
	  iEvent.getByLabel((*it).inputTag(),jetTags);
	  std::string moduleLabel = (jetTags).provenance()->moduleLabel();
	  if (mymap.find(moduleLabel) != mymap.end()) continue;
	  mymap[moduleLabel] == true;
	  
	  ith_tagged = this->TaggedJet(*jet,jetTags);
	  
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

                std::vector<TrackIPTagInfo::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
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
                std::vector<TrackIPTagInfo::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
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
                std::vector<TrackIPTagInfo::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
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


                std::vector<TrackIPTagInfo::TrackIPData>  ipdata =  (*tagInfo)[ith_tagged].impactParameterData();

                for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
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
            else if ( moduleLabel == "softMuonBJetTags" )
            {

                fS8evt->btag_SoftMuon_disc.push_back( (*jetTags)[ith_tagged].second);

                gotSMT = true;

            }

        }

	
	//FillHistos("n",p4MuJet, ptrel, JetFlavor,isbTaggedJet);
	
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
fS8evt->nmuons = total_nmuons;

// fill tree

ftree->Fill();

feventcounter++;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PerformanceAnalyzer);

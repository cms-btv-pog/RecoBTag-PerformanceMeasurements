
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
#include "DataFormats/BTauReco/interface/JetTagFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
//#include "DataFormats/HepMCCandidate/interfave/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/BTauReco/interface/JetTracksAssociation.h"
//#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
//#include "DataFormats/BTauReco/interface/TrackCountingTagInfoFwd.h"
//#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfo.h"
//#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfoFwd.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "PhysicsTools/Utilities/interface/DeltaR.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

//generator level + CLHEP
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/GenVertex.h"

// HepPDT // for simtracks
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

  //JetTrackAssociatorTags_ = iConfig.getParameter<std::string>("JetTracks");

  SimTrkCollectionTags_ = iConfig.getParameter<std::string>("SimTracks");
  
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
  
  // get list of taggers
  bTaggerList_ = iConfig.getUntrackedParameter<std::vector<std::string> >("bTaggerList");
  fnselectors= bTaggerList_.size();
  for(std::vector<std::string>::iterator objectName = bTaggerList_.begin(); objectName != bTaggerList_.end(); ++objectName) {
    std::cout << "*objectName  " << *objectName << std::endl;
	  if ( *objectName == "trackCountingHighEffJetTags" ) {
		  
		  moduleLabel_.push_back("trackCountingHighEffJetTags");
		 
	  }
	  if ( *objectName == "trackCountingHighPurJetTags" ) {
		  
		  moduleLabel_.push_back("trackCountingHighPurJetTags");
		 
	  }
	  if ( *objectName == "negativeTrackCounting2ndTrck" ) {
		  
		  moduleLabel_.push_back("negativeTrackCounting2ndTrck");
		  
	  }
	  if ( *objectName == "negativeTrackCounting3rdTrck" ) {
		  
		  moduleLabel_.push_back("negativeTrackCounting3rdTrck");

	  }
	  if ( *objectName == "jetProbabilityJetTagsNegativeOnly" ) {

		  moduleLabel_.push_back("jetProbabilityJetTagsNegativeOnly");

	  }
	  if ( *objectName == "jetProbabilityJetTagsPositiveOnly" ) {

		  moduleLabel_.push_back("jetProbabilityJetTagsPositiveOnly");

	  }
	  if ( *objectName == "jetProbabilityJetTags" ) {
		  
		  moduleLabel_.push_back("jetProbabilityJetTags");
		  
	  }
	  if ( *objectName == "softElectronJetTags" ) {

		  moduleLabel_.push_back("softElectronJetTags");
	  }
	  if ( *objectName == "softMuonJetTags" ) {

		  moduleLabel_.push_back("softMuonJetTags");
		  
	  }
	  
  }

  // Flavour identification
  flavourMatchOptionf = iConfig.getParameter<std::string>( "flavourMatchOption" );
  if (flavourMatchOptionf == "fastMC") {
	  flavourSourcef = iConfig.getParameter<edm::InputTag>("flavourSource");
  } else if (flavourMatchOptionf == "hepMC") {
	  jetFlavourIdentifier_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters"));
  } else if (flavourMatchOptionf == "genParticle") {
	  flavourSourcef = iConfig.getParameter<edm::InputTag> ("flavourSource");
  }
  
  // get operating points
  bTagCutList_ = iConfig.getUntrackedParameter<std::vector<double> >("bTagCutList");
  fOPMap["TCL"] = bTagCutList_[0];
  fOPMap["TCM"] = bTagCutList_[1];
  fOPMap["TCT"] = bTagCutList_[2];
  fOPMap["JPL"] = bTagCutList_[3];
  fOPMap["JPM"] = bTagCutList_[4];
  fOPMap["JPT"] = bTagCutList_[5];
  
  // get tagger for away jet
  fAwayJetTagger = iConfig.getParameter<std::string>("AwayJetTagger");

  // write performance plots?
  fWritePerformancePlots = iConfig.getParameter< bool > ("WritePerformancePlots");
 
  // include weights?
  fIncludeWeights = iConfig.getParameter< bool > ("IncludeWeights");

  
  topdir = rootFile_->mkdir("Histograms");
  topdir->cd();
  topdir->mkdir("ptrel");
  topdir->cd();
  topdir->mkdir("efficiencies");
  topdir->cd();
  topdir->mkdir("samples");
  topdir->cd();
  rootFile_->cd();
  
  // initialize histograms
  EffHistos     = new BTagHistograms();
  PtrelHistos   = new BTagHistograms();
  MujetHistos   = new BTagHistograms();
  AwayjetHistos = new BTagHistograms();
  TaggedMujetHistos   = new BTagHistograms();
  TaggedAwayjetHistos = new BTagHistograms();

  EffHistos->Init("efficiencies");
  PtrelHistos->Init("ptrel");
  MujetHistos->Init("n");
  AwayjetHistos->Init("p");
  MujetHistos->Init("n","b");
  AwayjetHistos->Init("p","b");
  MujetHistos->Init("n","cl");
  AwayjetHistos->Init("p","cl");
	  
  for ( std::map<std::string,float>::const_iterator imap = fOPMap.begin(); imap != fOPMap.end(); ++imap) {

	  EffHistos->Init("efficiencies",imap->first);
	  PtrelHistos->Init("ptrel",imap->first);
	  
	  TaggedMujetHistos->Init("ntag",imap->first);
	  TaggedAwayjetHistos->Init("ptag",imap->first);
	  TaggedMujetHistos->Init("ntag","b",imap->first);
	  TaggedAwayjetHistos->Init("ptag","b",imap->first);
	  TaggedMujetHistos->Init("ntag","cl",imap->first);
	  TaggedAwayjetHistos->Init("ptag","cl",imap->first);
	  
  }

  fperformanceTC2trk.Set("TC2trk");
  fperformanceTC3trk.Set("TC3trk");
  fperformanceTP.Set("TP");
  fperformanceTC2trk.SetMinDiscriminator(-1);
  fperformanceTC2trk.SetMaxDiscriminator(15);
  fperformanceTC3trk.SetMinDiscriminator(-1);
  fperformanceTC3trk.SetMaxDiscriminator(15);
  fperformanceTP.SetMinDiscriminator(0);
  fperformanceTP.SetMaxDiscriminator(1);

  //simUnit_= 1.0;  // starting with CMSSW_1_2_x ??
  //if ( (edm::getReleaseVersion()).find("CMSSW_1_1_",0)!=std::string::npos){
  //  simUnit_=0.1;  // for use in  CMSSW_1_1_1 tutorial
  //}
  //simUnit_= 1.;  // apparently not, still need this

  feventcounter = 0;
  
}


PerformanceAnalyzer::~PerformanceAnalyzer()
{
	rootFile_->cd();
	ftree->Write();

	topdir->cd();
	topdir->cd("efficiencies");
	EffHistos->Save();
	if (fWritePerformancePlots) {
		
		fperformanceTC2trk.Eval();
		fperformanceTC3trk.Eval();
		fperformanceTP.Eval();
		
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
		
		gTC2_b->SetName("gTC2_b");
		gTC2_c->SetName("gTC2_c");
		gTC2_udsg->SetName("gTC2_udsg");
		gTC3_b->SetName("gTC3_b");
		gTC3_c->SetName("gTC3_c");
		gTC3_udsg->SetName("gTC3_udsg");
		gTP_b->SetName("gTP_b");
		gTP_c->SetName("gTP_c");
		gTP_udsg->SetName("gTP_udsg");

		dTC2_udsg->SetName("discTC2_udsg");
		dTC3_udsg->SetName("discTC3_udsg");
		dTP_udsg->SetName("discTP_udsg");
		dTC2_udsg->SetTitle("TC2trk discriminator vs udsg-mistagging");
		dTC3_udsg->SetTitle("TC3trk discriminator vs udsg-mistagging");
		dTP_udsg->SetTitle("TP discriminator vs udsg-mistagging");
		
		gTC2_b->SetTitle("Jet b-efficiency");
		gTC2_c->SetTitle("Jet c-mistagging");
		gTC2_udsg->SetTitle("Jet udsg-mistagging");
		gTC3_b->SetTitle("Jet b-efficiency");
		gTC3_c->SetTitle("Jet c-mistagging");
		gTC3_udsg->SetTitle("Jet udsg-mistagging");
		gTP_b->SetTitle("Jet b-efficiency");
		gTP_c->SetTitle("Jet c-mistagging");
		gTP_udsg->SetTitle("Jet udsg-mistagging");
				
		gTC2_b->Write();
		gTC2_c->Write();
		gTC2_udsg->Write();
		gTC3_b->Write();
		gTC3_c->Write();
		gTC3_udsg->Write();
		gTP_b->Write();
		gTP_c->Write();
		gTP_udsg->Write();

		dTC2_udsg->Write();
		dTC3_udsg->Write();
		dTP_udsg->Write();
		
	}
	
	topdir->cd("ptrel");
	PtrelHistos->Save();
	topdir->cd("samples");
	MujetHistos->Save();
	AwayjetHistos->Save();
	TaggedMujetHistos->Save();
	TaggedAwayjetHistos->Save();
	
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
void PerformanceAnalyzer::beginJob(edm::EventSetup const& iSetup){
	
	std::cout << analyzer_ << " begin Job" << std::endl;

	rootFile_->cd();
  
  /*
  // track matching MC
  edm::ESHandle<TrackAssociatorBase> theChiAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2",theChiAssociator);
  associatorByChi2 = (TrackAssociatorBase *) theChiAssociator.product();
  
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
  associatorByHits = (TrackAssociatorByHits *) theHitsAssociator.product();
  */
}


void PerformanceAnalyzer::endJob() {

	std::cout << analyzer_ << " Total events processed: " << feventcounter << std::endl;
	
	rootFile_->cd();

}

reco::GenJet PerformanceAnalyzer::GetGenJet(reco::CaloJet calojet, reco::GenJetCollection genJetColl) {

  reco::GenJet matchedJet;

  double predelta = 99999.;
  for (GenJetCollection::const_iterator genjet = genJetColl.begin(); genjet != genJetColl.end(); ++genjet) {

    double delta  = ROOT::Math::VectorUtil::DeltaR(genjet->p4().Vect(), calojet.p4().Vect() );

    if ( delta < 0.2 && delta<predelta ) {
      matchedJet = *genjet;
      predelta = delta;
    }
  }
  
  return matchedJet;
}


int PerformanceAnalyzer::TaggedJet(reco::CaloJet calojet, edm::Handle<std::vector<reco::JetTag> > jetTags ) {
	
	double small = 1.e-5;
	int result = -1; // no tagged
	//int ith = -1;

	//std::cout << "calo jet: pz = " << calojet.pz() << " pt = " << calojet.pt() << std::endl;
	//for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
	//  edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];

	//get label and module names
	std::string moduleLabel = (jetTags).provenance()->moduleLabel();
	std::string processname = (jetTags).provenance()->processName();
	
	//std::cout << "[TaggedJet] moduleLabel = " << moduleLabel << std::endl;
	//std::cout << "[TaggedJet] processName = " << processname << std::endl;
	
	
	for (size_t t = 0; t < jetTags->size(); ++t) {
		edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].jet();
		if (jet_p.isNull()) {
			//std::cout << "-----------> JetTag::jet() returned null reference" << std::endl; 
			continue;
		}
		//std::cout << "[TaggedJet]  calojet pt = " << calojet.pt() << " tagged jet pt = " << jet_p->pt() << std::endl;
		if (DeltaR<reco::Candidate>()( calojet, *jet_p ) < small) {
			
			result = (int) t;
			
		}
	}
	/*

	
	for ( reco::JetTagCollection::const_iterator jetTag = taggedColl.begin(); jetTag != taggedColl.end(); ++jetTag ) {		
		double deltar  = ROOT::Math::VectorUtil::DeltaR( calojet.p4().Vect(), (*jetTag).jet()->p4().Vect() );
		double deltapt = std::abs( calojet.pt() - (*jetTag).jet()->pt() );
		//std::cout << "   deltar = " << deltar << "  deltapt = " << deltapt << std::endl;
		// check if calo jet is a tagged jet
		if ( deltar < 0.05 && deltapt < small ) {
			result = ith;
			//std::cout << "  tag jet: pz = " << (*jetTag).jet().pz() << " pt = " << (*jetTag).jet().pt() << std::endl;
		}
		ith++;
	}
	*/
	
	//if (result==-1) std::cout << " no jet tagged" << std::endl;
	
	return result;
}



SimTrack PerformanceAnalyzer::GetGenTrk(reco::Track atrack, const edm::SimTrackContainer *simTrkColl, const edm::SimVertexContainer *simVtcs) {

	  SimTrack matchedTrk;
	  edm::SimVertexContainer mysimVtcs = *simVtcs;

  double predelta = 99999.;
  for (SimTrackContainer::const_iterator gentrk = simTrkColl->begin(); gentrk != simTrkColl->end(); ++gentrk) {

	  double delta  = ROOT::Math::VectorUtil::DeltaR( TVector3((*gentrk).momentum().x(),(*gentrk).momentum().y(),(*gentrk).momentum().z()) , atrack.momentum() );
	  int type = (*gentrk).type();
	if ( delta < 0.2 && delta<predelta && ((*gentrk).charge() == atrack.charge() ) && 
		( abs(type)==11 || abs(type)==13 || abs(type)==15 || abs(type)==211 || abs(type)==321 ) ) {
		  matchedTrk = *gentrk;
		  predelta = delta;
		  HepLorentzVector v = (mysimVtcs)[(*gentrk).vertIndex()].position();
		  
		  //std::cout << "gentrk: vx = " << v.x() << std::endl;
		  //std::cout << "rectrk: vx = " << atrack.vx() << std::endl;
		  
	  }
  }
  
  return matchedTrk;
}

int
PerformanceAnalyzer::GetMotherId(const edm::SimVertexContainer *simVtxColl, const edm::SimTrackContainer *simTrkColl, SimTrack muonMC) {

  edm::SimVertexContainer mysimVtxColl = *simVtxColl;

  edm::SimTrackContainer mysimTrkColl = *simTrkColl;

  // fill map of simtracks
  std::map<unsigned, unsigned> geantToIndex;
  for( unsigned it=0; it< mysimTrkColl.size(); ++it ) {
    geantToIndex[ mysimTrkColl[it].trackId() ] = it;
  }

  // The origin vertex
  int vertexId = muonMC.vertIndex();
  SimVertex vertex = mysimVtxColl[vertexId];

  // The mother track 
  int motherId = -1;
  if( vertex.parentIndex() ) { // there is a parent to this vertex
 
    // geant id of the mother
    unsigned motherGeandId =   vertex.parentIndex(); 
    std::map<unsigned, unsigned >::iterator association  
      = geantToIndex.find( motherGeandId );
    if(association != geantToIndex.end() )
      motherId = association->second;
  }
  
  
  int motherType = motherId == -1 ? 0 : mysimTrkColl[motherId].type();
  
  return motherType;

}

//______________________________________________________________________________________________________________________
std::map< std::string, bool >
PerformanceAnalyzer::GetBTaggingMap(reco::CaloJet jet,std::vector<edm::Handle<std::vector<reco::JetTag> > > jetTags_testManyByType) {

	std::map< std::string, bool > aMap;
	
	int ith_tagged = -1;

	aMap["TCL"] = false;
	aMap["TCM"] = false;
	aMap["TCT"] = false;
	aMap["JPL"] = false;
	aMap["JPM"] = false;
	aMap["JPT"] = false;
	
	//start loop over all jetTags
	
	for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
		
			edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];

			
			ith_tagged = this->TaggedJet(jet,jetTags);

			if (ith_tagged == -1) continue;
		
			std::string moduleLabel = (jetTags).provenance()->moduleLabel();
			std::string processName = (jetTags).provenance()->processName();


			
			// Track Counting taggers
			//*********************************


			if ( moduleLabel == "trackCountingHighEffJetTags" ) {

				if ( (*jetTags)[ith_tagged].discriminator() > fOPMap["TCL"] ) aMap["TCL"] = true; // 2nd trk, 3D
				if ( (*jetTags)[ith_tagged].discriminator() > fOPMap["TCM"] ) aMap["TCM"] = true; // 2nd trk, 3D
				
			}
			else if ( moduleLabel == "trackCountingHighPurJetTags" ) {

				 if ( (*jetTags)[ith_tagged].discriminator() > fOPMap["TCT"] ) aMap["TCT"] = true; // 3rd trk, 3D
				 
			}
			//*********************************
			// Jet Probability taggers
			//*********************************
			else if ( moduleLabel == "jetProbabilityJetTags" ) {

				if ( (*jetTags)[ith_tagged].discriminator() > fOPMap["JPL"] ) aMap["JPL"] = true;
				if ( (*jetTags)[ith_tagged].discriminator() > fOPMap["JPM"] ) aMap["JPM"] = true;
				if ( (*jetTags)[ith_tagged].discriminator() > fOPMap["JPT"] ) aMap["JPT"] = true;
				
			}
			
	}
	
	return aMap;
}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillPerformance(reco::CaloJet jet, int JetFlavor, std::vector<edm::Handle<std::vector<reco::JetTag> > > jetTags_testManyByType) {

	std::map< std::string, bool > aMap;
	int ith_tagged = -1;
	
	//start loop over all jetTags
	for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
		
		edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];
			
		ith_tagged = this->TaggedJet(jet,jetTags);

		if (ith_tagged == -1) continue;
		
		std::string moduleLabel = (jetTags).provenance()->moduleLabel();
			
		//*********************************
		// Track Counting taggers
		//*********************************		
		if ( moduleLabel == "trackCountingHighEffJetTags" ) {
			
			fperformanceTC2trk.Add( (*jetTags)[ith_tagged].discriminator(), JetFlavor );			
		}
		else if ( moduleLabel == "trackCountingHighPurJetTags" ) {

			fperformanceTC3trk.Add( (*jetTags)[ith_tagged].discriminator(), JetFlavor );		
		}
		//*********************************
		// Jet Probability taggers
		//*********************************
		else if ( moduleLabel == "jetProbabilityJetTags" ) {

			fperformanceTP.Add( (*jetTags)[ith_tagged].discriminator(), JetFlavor );
		}
	}
}
	
//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillEff(TLorentzVector p4MuJet, int JetFlavor, std::map<std::string, bool> aMap, double weight) {
  
	std::string flavor = "";
	for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap ) {

		if (imap->second) {

			std::string tag = "_"+imap->first;
			flavor = "";

			EffHistos->Fill1d("jet_pt"+tag,p4MuJet.Pt(), weight );
			EffHistos->Fill1d("jet_eta"+tag,p4MuJet.Eta(), weight );
			
			if ( JetFlavor == 5 ) {
				flavor = "_b";
			}
			else if ( JetFlavor == 4 ) {
				flavor = "_c";
			}
			else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 ) {
				flavor = "_udsg";
			}
		
			if (flavor!="") {
			  EffHistos->Fill1d("jet_pt"+flavor+tag,p4MuJet.Pt(), weight );
			  EffHistos->Fill1d("jet_eta"+flavor+tag,p4MuJet.Eta(),weight );
			}			
		}
	}

	EffHistos->Fill1d("jet_pt",p4MuJet.Pt(),weight );
	EffHistos->Fill1d("jet_eta",p4MuJet.Eta(),weight );
			
	if ( JetFlavor == 5 ) {
		flavor = "_b";
	}
	else if ( JetFlavor == 4 ) {
		flavor = "_c";
	}
	else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 ) {
		flavor = "_udsg";
	}
	
	if (flavor!="") {
	  EffHistos->Fill1d("jet_pt"+flavor,p4MuJet.Pt(), weight );
	  EffHistos->Fill1d("jet_eta"+flavor,p4MuJet.Eta(),weight );
	}
	
}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillPtrel(double ptrel, int JetFlavor, std::map<std::string, bool> aMap, double weight) {

	std::string flavor = "";
	
	for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap ) {

		if (imap->second) {

			std::string tag = "_"+imap->first;
			flavor = "";

			PtrelHistos->Fill1d("jet_ptrel"+tag,ptrel,weight);
			
			if ( JetFlavor == 5 ) {
				flavor = "_b";
			}
			else if ( JetFlavor == 4 ) {
				flavor = "_c";
			}
			else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 ) {
				flavor = "_udsg";
			}
		
			if (flavor!="") PtrelHistos->Fill1d("jet_ptrel"+flavor+tag,ptrel,weight);

		}
	}

	PtrelHistos->Fill1d("jet_ptrel",ptrel,weight);
	
	if ( JetFlavor == 5 ) {
		flavor = "_b";
	}
	else if ( JetFlavor == 4 ) {
		flavor = "_c";
	}
	else if ( (JetFlavor >0 && JetFlavor<4) || JetFlavor==21 ) {
		flavor = "_udsg";
	}
	
	if (flavor!="") PtrelHistos->Fill1d("jet_ptrel"+flavor,ptrel,weight);
	
}

//______________________________________________________________________________________________________________________
void PerformanceAnalyzer::FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
				     int JetFlavor, std::map<std::string, bool> aMap, double weight)
{
  
	if ( type == "n") {
		
			
	  MujetHistos->Fill2d(type+"_pT",p4MuJet.Pt(),ptrel,weight);
	  MujetHistos->Fill2d(type+"_eta",TMath::Abs(p4MuJet.Eta()),ptrel,weight);
		if ( JetFlavor == 5 ) {
			std::string flavor = "b";
			MujetHistos->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
			MujetHistos->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
		}
		if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 ) {
			std::string flavor = "cl";
			MujetHistos->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
			MujetHistos->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
		}
	}
	else if ( type == "p") {
					
	  AwayjetHistos->Fill2d(type+"_pT",p4MuJet.Pt(),ptrel,weight);
	  AwayjetHistos->Fill2d(type+"_eta",TMath::Abs(p4MuJet.Eta()),ptrel, weight);
		if ( JetFlavor == 5 ) {
			std::string flavor = "b";
			AwayjetHistos->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
			AwayjetHistos->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);	
		}
		if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 ) {
			std::string flavor = "cl";
			AwayjetHistos->Fill2d(type+"_pT_"+flavor,p4MuJet.Pt(),ptrel,weight);
			AwayjetHistos->Fill2d(type+"_eta_"+flavor,TMath::Abs(p4MuJet.Eta()),ptrel,weight);
		}
	}
				
	for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap ) {

		if ( imap->second ) {
			
			if ( type == "n") {
				
			  TaggedMujetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4MuJet.Pt(),ptrel,weight);
			  TaggedMujetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
				
				if ( JetFlavor == 5 ) {
					std::string flavor = "b_";
					TaggedMujetHistos->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
					TaggedMujetHistos->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
				}
				if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 ) {
					std::string flavor = "cl_";
					TaggedMujetHistos->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
					TaggedMujetHistos->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
				}
			}
			else if ( type == "p") {

				
			  TaggedAwayjetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4MuJet.Pt(),ptrel,weight);
			  TaggedAwayjetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
			
				if ( JetFlavor == 5 ) {
					std::string flavor = "b_";
					TaggedAwayjetHistos->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
					TaggedAwayjetHistos->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
				}
				if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 ) {
					std::string flavor = "cl_";
					TaggedAwayjetHistos->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
					TaggedAwayjetHistos->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
				}
			}
		}
	}
	
}

JetFlavour PerformanceAnalyzer::getMatchedParton(const reco::CaloJet &jet)
{
  JetFlavour jetFlavour;
  
  if (flavourMatchOptionf == "fastMC") {

    jetFlavour.underlyingParton4Vec(jet.p4());
    //const edm::RefToBase<reco::Jet> & caloRefTB = jetTag.jet();
    //const reco::CaloJetRef & caloRef = jet.castTo<reco::CaloJetRef>(); 
    //jetFlavour.flavour(flavoursMapf[caloRef]); 

  } else if (flavourMatchOptionf == "hepMC") {

    jetFlavour = jetFlavourIdentifier_.identifyBasedOnPartons(jet);

  } else if (flavourMatchOptionf == "genParticle") {

    for( reco::CandMatchMap::const_iterator f  = theJetPartonMapf->begin();
	 f != theJetPartonMapf->end(); f++) {
      const reco::Candidate *theJetInTheMatchMap = &*(f->key);    
      const reco::Candidate *theMatchedParton    = &*(f->val);
      if(theJetInTheMatchMap->hasMasterClone ()) {
	const reco::CaloJet* theMasterClone = dynamic_cast<const reco::CaloJet*>(theJetInTheMatchMap->masterClone().get());
	//std::cout << " masterclone pt = " << theMasterClone->pt() << " calo jet pt = " << jet.pt() << std::endl;
	// FIXME, compare pointers rather than values:
	if ( ( theMasterClone->pt() - jet.pt() ) < 1.e-5 ) {
		//std::cout << " it matches! " << std::endl;
	  jetFlavour.flavour(abs(theMatchedParton->pdgId()));
	  jetFlavour.underlyingParton4Vec(theMatchedParton->p4());
	  return jetFlavour;
	}
      }
    }

  }

  return jetFlavour;
}

// ------------ method called to produce the data  ------------
void
PerformanceAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{

		
  // initialize flavour identifiers
  edm::Handle<JetFlavourMatchingCollection> jetMC;

  if (flavourMatchOptionf == "fastMC") {
    iEvent.getByLabel(flavourSourcef, jetMC);
    for(JetFlavourMatchingCollection::const_iterator iter =
	  jetMC->begin(); iter != jetMC->end(); iter++)
      flavoursMapf.insert(*iter);
  } else if (flavourMatchOptionf == "hepMC") {
    jetFlavourIdentifier_.readEvent(iEvent);
  } else if (flavourMatchOptionf == "genParticle") {
    iEvent.getByLabel (flavourSourcef, theJetPartonMapf);
  }

  //jetFlavourIdentifier_.readEvent(iEvent);
	//jetFlavourIdentifier2_.readEvent(iEvent);
	
  // generator candidates
  Handle<CandidateCollection> genParticles;
  iEvent.getByLabel("genParticleCandidates", genParticles);

    // Trakcs
	Handle<reco::TrackCollection> recTrks;
	iEvent.getByLabel(recoTrackList_, recTrks);
	
	// Primary vertices
	//Handle<reco::VertexCollection> recVtxs;
	//iEvent.getByLabel(recoVtxList_, recVtxs );
	
	//
	//Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
	//iEvent.getByLabel(JetTrackAssociatorTags_, jetTracksAssociation);

	// Muons
	Handle<reco::MuonCollection> muonsColl;
	iEvent.getByLabel(MuonCollectionTags_, muonsColl);

	// Calo Jets
	Handle<reco::CaloJetCollection> jetsColl;
	iEvent.getByLabel(CaloJetCollectionTags_, jetsColl);

	//Handle<reco::CaloJetCollection> jetsCorrColl;
	//iEvent.getByLabel(CorrCaloJetCollectionTags_, jetsCorrColl);

    // initialize jet corrector
	const JetCorrector *acorrector = JetCorrector::getJetCorrector("MCJetCorrectorIcone5",iSetup);
	
	Handle<reco::GenJetCollection> genjetsColl;
	iEvent.getByLabel(GenJetCollectionTags_, genjetsColl);
		
	Handle<edm::SimVertexContainer> simVtxColl;
	Handle<edm::SimTrackContainer> simtrkColl;
	const edm::SimVertexContainer *simVtcs = 0;
	const edm::SimTrackContainer *simTrks = 0;
	if (flavourMatchOptionf == "hepMC") {

	  iEvent.getByLabel( "g4SimHits", simVtxColl);
	  simVtcs = simVtxColl.product();
	
	  iEvent.getByLabel(SimTrkCollectionTags_, simtrkColl);
	  simTrks = simtrkColl.product();
	  
	}


	// truth tracks
	//edm::Handle<TrackingParticleCollection>  TPCollectionH ;
	//iEvent.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
	//const TrackingParticleCollection tPC = *(TPCollectionH.product());

	// Tag Jets	
	std::vector<edm::Handle<std::vector<reco::JetTag> > > jetTags_testManyByType ;

	edm::Handle<reco::JetTagCollection> tagHandle;
	for (size_t im = 0; im!= moduleLabel_.size(); ++im ) {
	  iEvent.getByLabel(moduleLabel_[im], tagHandle);
		jetTags_testManyByType.push_back(tagHandle);
	}
	
	
	const reco::CaloJetCollection recoJets =   *(jetsColl.product());
	const reco::GenJetCollection  genJets  =   *(genjetsColl.product());
	const reco::MuonCollection    recoMuons =  *(muonsColl.product());
	

	
	//const reco::VertexCollection recoPV = *(recVtxs.product());
	
	// MC
  
	//Handle<SimTrackContainer> simTrks;
	//iEvent.getByLabel( simG4_, simTrks);
	
	Handle<std::vector<reco::TrackIPTagInfo> > tagInfo;
	iEvent.getByLabel("impactParameterTagInfos", tagInfo);
		
	Handle<HepMCProduct> evtMC;
  	      
	try{
		iEvent.getByLabel("source",evtMC);
		if(fverbose){
			std::cout << "source HepMCProduct found"<< std::endl;
		}
	
	} catch(const Exception&) {
	
		if(fverbose){
			std::cout << "no HepMCProduct found"<< std::endl;
		}
	}
  	// WEIGHTS   
	double weight = 1.;
	if (fIncludeWeights){
	  Handle< double> weightHandle;
	  iEvent.getByLabel ("csaweightproducer", weightHandle);
	  weight = *weightHandle;
	  
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
	
	for( jet = recoJets.begin(); jet != recoJets.end(); ++jet ) {

		// get jet corrections
		double jetcorrection =  acorrector->correction(*jet, iEvent, iSetup);
		// Jet quality cuts
		if ( (jet->pt() * jetcorrection ) <= MinJetPt_ || std::abs( jet->eta() ) >= MaxJetEta_ ) continue;
		// get MC flavor of jet
		//int JetFlavor = jetFlavourIdentifier_.identifyBasedOnPartons(*jet).flavour();
		int JetFlavor = getMatchedParton(*jet).flavour();
		
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
		for ( muon = recoMuons.begin(); muon != recoMuons.end(); ++muon) {
			
			Track muonTrk = *muon->track();
			//TrackingParticleRef TrueHitsTrk;
			Track muonSA = *muon->standAloneMuon();
			int nhit = muonTrk.numberOfValidHits();//muonTrk.recHitsSize();
			
			// muon cuts
			double normChi2 = (*(muon->combinedMuon())).normalizedChi2();//(*(muon->combinedMuon())).chi2() / (*(muon->combinedMuon())).ndof();// use global fit
			if ( (nhit <= MinMuonNHits_ ) || (muon->pt()<= MinMuonPt_) || (normChi2 >= MaxMuonChi2_ ) ) continue;

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
			ptrel = leptonvec.Perp(tmpvec);
			// muon in jet cuts
			if ( (deltaR >= MinDeltaR_ ) || (ptrel <= MinPtRel_ ) ) continue;

			// now we have a good muon in a jet
			total_nmuons++;
			hasLepton = 1;
			tmptotmuon++;
			// pick the leading muon inside the jet
			if ( muon->pt() > mu_highest_pt ) {
				mu_highest_pt = muon->pt();
				p4Muon.SetPtEtaPhiE(muon->pt(),
									muon->eta(),
									muon->phi(),
									muon->energy());
				// recalculate pTrel
				tmpvec = jetcorrection * tmpvecOrg;
				leptonvec.SetXYZ(muon->px(), muon->py(),muon->pz());
				tmpvec += leptonvec;
				ptrel = leptonvec.Perp(tmpvec);
			}

			// collect muon data
			leptonEvent.pdgid.push_back( 13 ); // muon only for the moment
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

			leptonEvent.jet_ptrel.push_back( ptrel);
			leptonEvent.jet_deltaR.push_back( deltaR);

			// find a sim track
			if (flavourMatchOptionf == "hepMC" ) {
			  SimTrack genlepton = this->GetGenTrk(muonTrk, simTrks, simVtcs );
												 
			  leptonEvent.mc_pt.push_back(           genlepton.momentum().perp());
			  leptonEvent.mc_phi.push_back(          genlepton.momentum().phi());
			  leptonEvent.mc_eta.push_back(          genlepton.momentum().pseudoRapidity());
			  leptonEvent.mc_e.push_back(           genlepton.momentum().e());
			  leptonEvent.mc_charge.push_back(       genlepton.charge());
			  leptonEvent.mc_pdgid.push_back(        genlepton.type());
			  leptonEvent.mc_mother_pdgid.push_back( GetMotherId(simVtcs, simTrks, genlepton));
												 
			} else if (flavourMatchOptionf == "genParticle") {
			  
			  for(size_t i = 0; i < genParticles->size(); ++ i) {
			    const Candidate & p = (*genParticles)[i];

			    if ( abs(p.pdgId()) == 13 && p.status()==2 ) {

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
		
		if ( hasLepton == 1 ) {
		  p4MuJet.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy() );
		  p4MuJet = jetcorrection * p4MuJet;
		}
		  

		/////////////////////////////
		// find away jet
		////////////////////////////
		bool AwayTaggedJet = false;
		bool AwayMuonJet = false;
		TLorentzVector p4AwayMuon;
		
		for( awayjet = recoJets.begin(); awayjet != recoJets.end(); ++awayjet ) {

		  if ( hasLepton == 1 ) continue;

			TLorentzVector p4AwayJet;
			p4AwayJet.SetPtEtaPhiE(awayjet->pt(), awayjet->eta(), awayjet->phi(), awayjet->energy() );
			p4AwayJet = p4AwayJet * ( acorrector->correction(*awayjet, iEvent, iSetup) );

			// Jet quality cuts
			if ( (awayjet->pt() * ( acorrector->correction(*awayjet, iEvent, iSetup) ) ) <= MinJetPt_ || std::abs( awayjet->eta() ) >= MaxJetEta_ ) continue;
		
			// skip muon in jet
			if ( p4AwayJet == p4MuJet ) continue;

			// now we have an away jet

			// find an away tagged jet
			if ( !AwayTaggedJet ) {

				std::map< std::string, bool > aBmap = this->GetBTaggingMap(*awayjet,jetTags_testManyByType);
				for (std::map<std::string,bool>::const_iterator imap = aBmap.begin(); imap != aBmap.end(); ++imap ) {

					if ( imap->first == fAwayJetTagger && imap->second ) {

						AwayTaggedJet = true;
					}
				}
				
			}
			
			// find an away muon in jet
			if ( !AwayMuonJet ) {

				mu_highest_pt = 0;
				for ( muon = recoMuons.begin(); muon != recoMuons.end(); ++muon) {

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
					if ( muon->pt() > mu_highest_pt ) {
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

		std::map< std::string, bool > thebtaggingmap = this->GetBTaggingMap(*jet,jetTags_testManyByType);
		FillEff(p4Jet, JetFlavor, thebtaggingmap, weight );
		FillPerformance(*jet, JetFlavor, jetTags_testManyByType );
			
		if ( hasLepton == 1 ) {
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
		if ( genjet.p() !=0 ) {
			//fS8evt->genjet_p.push_back(genjet.p());
			fS8evt->genjet_pt.push_back(genjet.pt());
			fS8evt->genjet_eta.push_back(genjet.eta());
			fS8evt->genjet_phi.push_back(genjet.phi());
			fS8evt->genjet_e.push_back(genjet.energy());
			//fS8evt->genjet_vx.push_back(genjet.vx());
			//fS8evt->genjet_vy.push_back(genjet.vy());
			//fS8evt->genjet_vz.push_back(genjet.vz());
		} else {
			//fS8evt->genjet_p.push_back(-100000);
			fS8evt->genjet_pt.push_back(-100000);
			fS8evt->genjet_eta.push_back(-100000);
			fS8evt->genjet_phi.push_back(-100000);
			fS8evt->genjet_e.push_back(-100000);
			//fS8evt->genjet_vx.push_back(-100000);
			//fS8evt->genjet_vy.push_back(-100000);
			//fS8evt->genjet_vz.push_back(-100000);
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

		//start loop over all jetTags	
		for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
			edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];
			
			ith_tagged = this->TaggedJet(*jet,jetTags);

			if (ith_tagged == -1) continue;
		
			std::string moduleLabel = (jetTags).provenance()->moduleLabel();
			std::string processName = (jetTags).provenance()->processName();

			//*********************************
			// Track Counting taggers
			//*********************************

			if ( moduleLabel == "trackCountingHighEffJetTags" ) {
			  
			  std::vector< Measurement1D  > trackIP = (*tagInfo)[ith_tagged].impactParameters(0);
			  if((trackIP).size()>=2){
			    float iptrack1 = (trackIP)[0].significance();
			    fS8evt->btag_TrkCounting_disc3D_1trk.push_back( iptrack1 );
			  }
			  
			  fS8evt->btag_TrkCounting_disc3D_2trk.push_back( (*jetTags)[ith_tagged].discriminator() ); // 2nd trk, 3D
			  
			  gotTCHE = true;
			  
			}
			else if ( moduleLabel == "trackCountingHighPurJetTags" ) {

				fS8evt->btag_TrkCounting_disc3D_3trk.push_back( (*jetTags)[ith_tagged].discriminator() ); // 3rd trk, 3D
				gotTCHP = true;
			}
			else if ( moduleLabel == "negativeTrackCounting2ndTrck" ) {

			  std::vector< Measurement1D  > trackIP = (*tagInfo)[ith_tagged].impactParameters(0);
			  if((trackIP).size()>=2){
			    float iptrack1 = (trackIP)[(trackIP).size()-1].significance();
			    fS8evt->btag_NegTag_disc3D_1trk.push_back( iptrack1 );
			  }
			

			  std::cout << "discri neg tc " << (*jetTags)[ith_tagged].discriminator() << std::endl;
				fS8evt->btag_NegTag_disc3D_2trk.push_back( (*jetTags)[ith_tagged].discriminator() ); // 2nd trk, 3D
				
				gotTCHEneg = true;
			}
			else if ( moduleLabel == "negativeTrackCounting3rdTrck" ) {
				
				fS8evt->btag_NegTag_disc3D_3trk.push_back( (*jetTags)[ith_tagged].discriminator() ); // 3rd trk, 3D
				
				gotTCHPneg = true;
			}

			//*********************************
			// Jet Probability taggers
			//*********************************
			else if ( moduleLabel == "jetProbabilityJetTags" ) {

				fS8evt->btag_JetProb_disc3D.push_back( (*jetTags)[ith_tagged].discriminator());
				
				gotJP = true;

				std::string moduleLabel = (jetTags).provenance()->moduleLabel();
				
				int NtrksInJet = (*jetTags)[ith_tagged].tracks().size();
				fS8evt->jet_ntrks.push_back( NtrksInJet );
				fS8evt->trackProvaVector_Size = (*tagInfo)[ith_tagged].probabilities(0).size();
				if(StoreTrackProba_){
				  int i=0;
				  std::vector< float > track_proba = (*tagInfo)[ith_tagged].probabilities(0) ;
				  std::vector< float > probabilities;
				  for(std::vector<float>::const_iterator it = track_proba.begin(); it!=track_proba.end(); ++it, i++){
				    
				    double delta  = -2.; 
				    delta = ROOT::Math::VectorUtil::DeltaR( (*(*jetTags)[ith_tagged].jet()).p4().Vect(), (*(*jetTags)[ith_tagged].tracks()[i]).momentum());
				    if(delta <0.3) probabilities.push_back((*it));
				    
				    
				  }
				  fS8evt->jet_Tracks_Probability.push_back(probabilities);
				}
				
			}
			else if ( moduleLabel == "jetProbabilityJetTagsNegativeOnly" ) {
			  
			  fS8evt->btag_negJetProb_disc3D.push_back( (*jetTags)[ith_tagged].discriminator());
			  
			  gotJPneg = true;
			  
			}
			else if ( moduleLabel == "jetProbabilityJetTagsPositiveOnly" ) {
			  
			  fS8evt->btag_posJetProb_disc3D.push_back( (*jetTags)[ith_tagged].discriminator());
			  
			  gotJPpos = true;
			  
			}
			//*********************************
			// SoftLeptons Taggers 
			//*********************************
			else if ( moduleLabel == "softMuonJetTags" ) {

				fS8evt->btag_SoftMuon_disc.push_back( (*jetTags)[ith_tagged].discriminator());
				
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

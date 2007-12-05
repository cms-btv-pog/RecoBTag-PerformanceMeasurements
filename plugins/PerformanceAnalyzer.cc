
#include "RecoBTag/PerformanceMeasurements/interface/PerformanceAnalyzer.h"

#include "PluginManager/ModuleDef.h"
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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/BTauReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfoFwd.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfoFwd.h"

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
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"

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
  simG4_=iConfig.getParameter<edm::InputTag>( "simG4" );
  
  // open output file to store results
  outputFile_  = iConfig.getUntrackedParameter<std::string>("outputFile");
  // create or update output root file
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
  // verbose std output
  verbose_= iConfig.getUntrackedParameter<bool>("verbose", false);
  // default tree
  ftree = new TTree("summary","summary");
  ftree->AutoSave();

  
  fS8evt = new BTagEvent();
  //std::cout << "summary branch created" << std::endl;
  ftree->Branch("s8.","BTagEvent",&fS8evt,64000,1); //system8 branch

  
  // get list of tracks
  recoTrackList_ = iConfig.getUntrackedParameter<std::string >("TrackCollection");
  // get list of PV
  recoVtxList_ = iConfig.getUntrackedParameter<std::string >("PrimaryVertexCollection");

  //
  MuonCollectionTags_ = iConfig.getParameter<std::string>("Muons");
  
  CaloJetCollectionTags_ = iConfig.getParameter<std::string>("Jets");

  //CorrCaloJetCollectionTags_ = iConfig.getParameter<std::string>("CorrJets");

  GenJetCollectionTags_ = iConfig.getParameter<std::string>("GenJets");

  //JetTrackAssociatorTags_ = iConfig.getParameter<std::string>("JetTracks");

  SimTrkCollectionTags_ = iConfig.getParameter<std::string>("SimTracks");
  
  jetFlavourIdentifier_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters"));
  jetFlavourIdentifier2_ = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters2"));


  MinJetEt_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinEt");
  MaxJetEta_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MaxEta");
  MinDeltaR_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinDeltaR");
  MinPtRel_  = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinPtRel");
  MinMuonPt_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MinMuonPt");
  MaxMuonEta_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MaxMuonEta");
  MaxMuonChi2_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MaxMuonChi2");
  MinMuonNHits_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<int>("MinNHits");
						 
  // get list of taggers, just one for the moment
  bTaggerList_ = iConfig.getUntrackedParameter<std::vector<std::string> >("bTaggerList");
  fnselectors= bTaggerList_.size();
  for(std::vector<std::string>::iterator objectName = bTaggerList_.begin(); objectName != bTaggerList_.end(); ++objectName) {
	  if ( *objectName == "TrackCountingHighEff" ) {
		  
		  moduleLabel_.push_back("trackCountingHighEffJetTags");
		 
	  }
	  if ( *objectName == "TrackCountingHighPur" ) {
		  
		  moduleLabel_.push_back("trackCountingHighPurJetTags");
		 
	  }
	  if ( *objectName == "TrackProbability" ) {
		  
		  moduleLabel_.push_back("trackProbabilityJetTags");
		  
	  }
	  if ( *objectName == "SoftElectron" ) {

		  moduleLabel_.push_back("softElectronJetTags");
	  }
	  if ( *objectName == "SoftMuon" ) {

		  moduleLabel_.push_back("softMuonJetTags");
		  
	  }
	  
  }
	  
  
  
  //simUnit_= 1.0;  // starting with CMSSW_1_2_x ??
  if ( (edm::getReleaseVersion()).find("CMSSW_1_1_",0)!=std::string::npos){
    simUnit_=0.1;  // for use in  CMSSW_1_1_1 tutorial
  }
  simUnit_= 1.;  // apparently not, still need this

  feventcounter = 0;
  
}


PerformanceAnalyzer::~PerformanceAnalyzer()
{
	ftree->Write();
		
	
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
    delete fS8evt;
	
	delete ftree;
	delete rootFile_;
}



//
// member functions
//
void PerformanceAnalyzer::beginJob(edm::EventSetup const& iSetup){
	
	std::cout << analyzer_ << " begin job" << std::endl;

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
	int ith = 0;

	//std::cout << "calo jet: pz = " << calojet.pz() << " pt = " << calojet.pt() << std::endl;
	//for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
	//  edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];

	//get label and module names
	std::string moduleLabel = (jetTags).provenance()->moduleLabel();

	for (size_t t = 0; t < jetTags->size(); ++t) {
		edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].jet();
		if (jet_p.isNull()) {
			/*std::cout << "-----------> JetTag::jet() returned null reference" << std::endl; */
			continue;
		}
		if (DeltaR<reco::Candidate>()( calojet, *jet_p ) < small) {
			
			result = ith;
			
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



SimTrack PerformanceAnalyzer::GetGenTrk(reco::Track atrack, edm::SimTrackContainer simTrkColl, edm::SimVertexContainer simVtcs) {

	  SimTrack matchedTrk;

  double predelta = 99999.;
  for (SimTrackContainer::const_iterator gentrk = simTrkColl.begin(); gentrk != simTrkColl.end(); ++gentrk) {

	  double delta  = ROOT::Math::VectorUtil::DeltaR( TVector3((*gentrk).momentum().x(),(*gentrk).momentum().y(),(*gentrk).momentum().z()) , atrack.momentum() );
	  int type = (*gentrk).type();
	if ( delta < 0.2 && delta<predelta && ((*gentrk).charge() == atrack.charge() ) && 
		( abs(type)==11 || abs(type)==13 || abs(type)==15 || abs(type)==211 || abs(type)==321 ) ) {
		  matchedTrk = *gentrk;
		  predelta = delta;
		  HepLorentzVector v = (simVtcs)[(*gentrk).vertIndex()].position();
		  
		  //std::cout << "gentrk: vx = " << v.x() << std::endl;
		  //std::cout << "rectrk: vx = " << atrack.vx() << std::endl;
		  
	  }
  }
  
  return matchedTrk;
}

int
PerformanceAnalyzer::GetMotherId(edm::SimVertexContainer simVtxColl, edm::SimTrackContainer simTrkColl, SimTrack muonMC) {

  // fill map of simtracks
  std::map<unsigned, unsigned> geantToIndex;
  for( unsigned it=0; it<simTrkColl.size(); ++it ) {
    geantToIndex[ simTrkColl[it].trackId() ] = it;
  }

  // The origin vertex
  int vertexId = muonMC.vertIndex();
  SimVertex vertex = simVtxColl[vertexId];

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
  
  
  int motherType = motherId == -1 ? 0 : simTrkColl[motherId].type();
  
  return motherType;

}

// ------------ method called to produce the data  ------------
void
PerformanceAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{

	// initialize flavour identifiers
	jetFlavourIdentifier_.readEvent(iEvent);
	jetFlavourIdentifier2_.readEvent(iEvent);
	
    // Trakcs
	Handle<reco::TrackCollection> recTrks;
	iEvent.getByLabel(recoTrackList_, recTrks);
	
	// Primary vertices
	Handle<reco::VertexCollection> recVtxs;
	iEvent.getByLabel(recoVtxList_, recVtxs );
	
	//
	//Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
	//iEvent.getByLabel(JetTrackAssociatorTags_, jetTracksAssociation);

	// Muons
	Handle<reco::MuonCollection> muonsColl;
	iEvent.getByLabel(MuonCollectionTags_, muonsColl);

	Handle<edm::SimTrackContainer> simtrkColl;
	iEvent.getByLabel(SimTrkCollectionTags_, simtrkColl);

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
	iEvent.getByLabel( "g4SimHits", simVtxColl);
	const edm::SimVertexContainer simVtcs =    *(simVtxColl.product());


	// truth tracks
	//edm::Handle<TrackingParticleCollection>  TPCollectionH ;
	//iEvent.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
	//const TrackingParticleCollection tPC = *(TPCollectionH.product());

	// Tag Jets
	std::vector<edm::Handle<std::vector<reco::JetTag> > > jetTags_testManyByType ;
	iEvent.getManyByType(jetTags_testManyByType);
	// Define the handles for the specific algorithms
	//edm::Handle<reco::SoftLeptonTagInfoCollection> jetsInfoHandle_sl;
	//edm::Handle<reco::TrackProbabilityTagInfoCollection> jetsInfoHandleTP;
	//edm::Handle<reco::TrackCountingTagInfoCollection> jetsInfoHandleTC;
	
	//edm::Handle<reco::JetTagCollection> tagHandle;
	//	const reco::JetTagCollection & btagCollTC;
	//const reco::JetTagCollection & btagCollTP;

	// Track Counting
	//iEvent.getByLabel(moduleLabel_[0], tagHandle);
	//const reco::JetTagCollection &btagCollTC = *(tagHandle.product());
	//edm::Handle<reco::TrackCountingTagInfoCollection> tagInfoHandleTC;
	//iEvent.getByLabel(moduleLabel_[0], tagInfoHandleTC);
	//const reco::TrackCountingTagInfoCollection & TrkCountingInfo = *(tagInfoHandleTC.product());
	
	// Track Probability
	//iEvent.getByLabel(moduleLabel_[1], tagHandle);
	//const reco::JetTagCollection &btagCollTP = *(tagHandle.product());
	//edm::Handle<reco::TrackProbabilityTagInfoCollection> tagInfoHandleTP;
	//iEvent.getByLabel(moduleLabel_[1], tagInfoHandleTP);
	//const reco::TrackProbabilityTagInfoCollection & JetProbInfo = *(tagInfoHandleTP.product());
	
	const reco::CaloJetCollection recoJets =   *(jetsColl.product());
	const reco::GenJetCollection  genJets  =   *(genjetsColl.product());
	const reco::MuonCollection    recoMuons =  *(muonsColl.product());
	const edm::SimTrackContainer simTrks =    *(simtrkColl.product());
	const reco::VertexCollection recoPV = *(recVtxs.product());
	
	// MC
  
	//Handle<SimTrackContainer> simTrks;
	//iEvent.getByLabel( simG4_, simTrks);

	bool MC=false;
	Handle<HepMCProduct> evtMC;
  
	      
	try{
		iEvent.getByLabel("source",evtMC);
		if(verbose_){
			std::cout << "source HepMCProduct found"<< std::endl;
		}
		MC=true;
	} catch(const Exception&) {
		MC=false;
		if(verbose_){
			std::cout << "no HepMCProduct found"<< std::endl;
		}
	}
  
  	  
	fS8evt->Reset();
	
	fS8evt->event = iEvent.id().event();
	fS8evt->run = iEvent.id().run();

	
	// loop over jets
	//fS8evt->njets = recoJets.size();
	//fS8evt->nmuons = recoMuons.size();
	int total_nmuons = 0;
	fS8evt->ngenjets = genJets.size();
	fS8evt->nvertices = recoPV.size();
	
	CaloJetCollection::const_iterator jet;
	//GenJetCollection::const_iterator genjet;
	reco::MuonCollection::const_iterator muon;

	//reco::JetTagCollection::iterator btagite;
	//std::vector< BTagLeptonEvent > LeptonEvtVector;
	
	int ijet = 0;	
	for( jet = recoJets.begin(); jet != recoJets.end(); ++jet ) {

		// initial set of cuts on jets
		double jetcorrection =  acorrector->correction(*jet, iEvent, iSetup);
		if ( (jet->pt() * jetcorrection ) <= MinJetEt_ || std::abs( jet->eta() ) >= MaxJetEta_ ) continue;
		int hasLepton = 0;
		int tmptotmuon = 0;
		BTagLeptonEvent leptonEvent;
		
		for ( muon = recoMuons.begin(); muon != recoMuons.end(); ++muon) {

			Track muonTrk = *muon->track();
			TrackingParticleRef TrueHitsTrk;

			// muon cuts
			double normChi2 = (*(muon->combinedMuon())).chi2() / (*(muon->combinedMuon())).ndof();// use global fit
			if ( (muon->pt()<= MinMuonPt_) || (normChi2 >= MaxMuonChi2_ ) ) continue;

			// lepton in jet cuts			
			double deltaR  = ROOT::Math::VectorUtil::DeltaR(jet->p4().Vect(), muonTrk.momentum() );// use tracker muon
			TVector3 tmpvec(jet->p4().Vect().X(),jet->p4().Vect().Y(),  jet->p4().Vect().Z());
			// use calibrated jet
			tmpvec = jetcorrection * tmpvec;
			TVector3 leptonvec(muonTrk.momentum().X(), muonTrk.momentum().Y(),muonTrk.momentum().Z());
			tmpvec += leptonvec;
			double ptrel = leptonvec.Perp(tmpvec);
			
			if ( (deltaR >= MinDeltaR_ ) || (ptrel <= MinPtRel_ ) ) continue;

			// now we have a good lepton in a jet
			total_nmuons++;
			hasLepton = 1;
			tmptotmuon++;
			
			leptonEvent.pdgid.push_back( 13 ); // muon only for the moment
			leptonEvent.e.push_back( muon->energy());
			leptonEvent.pt.push_back( muonTrk.pt());
			leptonEvent.eta.push_back( muonTrk.eta());
			leptonEvent.phi.push_back( muonTrk.phi());
			leptonEvent.charge.push_back( muonTrk.charge());
			//leptonEvent.p.push_back( muonTrk.p());
			leptonEvent.trkchi2.push_back( muonTrk.chi2());
			leptonEvent.trkndof.push_back( muonTrk.ndof());
			leptonEvent.chi2.push_back(    (*(muon->combinedMuon())).chi2() );
			leptonEvent.ndof.push_back(    (*(muon->combinedMuon())).ndof() );
			Track muonSA = *muon->standAloneMuon();
			int nhit = muonSA.recHitsSize();
			leptonEvent.SArechits.push_back(  nhit );
			leptonEvent.trkrechits.push_back( muonTrk.recHitsSize());
			leptonEvent.d0.push_back(         muonTrk.d0());
			leptonEvent.d0sigma.push_back(    muonTrk.d0Error());
			//leptonEvent.vx.push_back(         muonTrk.vx());
			//leptonEvent.vy.push_back(         muonTrk.vy());
			//leptonEvent.vz .push_back(        muonTrk.vz());
			
			// find a sim track
			SimTrack genlepton = this->GetGenTrk(muonTrk, simTrks, simVtcs );
			
			//leptonEvent.mc_p.push_back(            genlepton.momentum().vect().mag());
			leptonEvent.mc_pt.push_back(           genlepton.momentum().perp());
			leptonEvent.mc_phi.push_back(          genlepton.momentum().phi());
			leptonEvent.mc_eta.push_back(          genlepton.momentum().pseudoRapidity());
			leptonEvent.mc_e.push_back(           genlepton.momentum().e());
			leptonEvent.mc_charge.push_back(       genlepton.charge());
			//leptonEvent.mc_vx.push_back(           genlepton.momentum().x());
			//leptonEvent.mc_vy.push_back(           genlepton.momentum().y());
			//leptonEvent.mc_vz.push_back(           genlepton.momentum().z());
			leptonEvent.mc_pdgid.push_back(        genlepton.type());
			leptonEvent.mc_mother_pdgid.push_back( GetMotherId(simVtcs, simTrks, genlepton));
	
			leptonEvent.jet_ptrel.push_back( ptrel);
			leptonEvent.jet_deltaR.push_back( deltaR);
			
			//LeptonEvtVector.push_back( leptonEvent ));
		}

		fS8evt->lepton.push_back( leptonEvent );
		
		fS8evt->jet_hasLepton.push_back( hasLepton );
		
		fS8evt->jet_flavour_alg.push_back(jetFlavourIdentifier_.identifyBasedOnPartons(*jet).flavour());	
		fS8evt->jet_flavour_phy.push_back(jetFlavourIdentifier2_.identifyBasedOnPartons(*jet).flavour());	
		//fS8evt->jet_p.push_back(jet->p());
		fS8evt->jet_e.push_back(jet->energy());
		fS8evt->jet_pt.push_back(jet->pt());
		fS8evt->jet_eta.push_back(jet->eta());
		fS8evt->jet_phi.push_back(jet->phi());
		fS8evt->jet_et.push_back(jet->et());
		//fS8evt->jet_vx.push_back(jet->vx());
		//fS8evt->jet_vy.push_back(jet->vy());
		//fS8evt->jet_vz.push_back(jet->vz());
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
		int isbtagged = 0;

		bool gotTCHE = false;
		bool gotTCHP = false;
		bool gotJP   = false;
		bool gotSMT  = false;
		
		for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
			edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];
			
			ith_tagged = this->TaggedJet(*jet,jetTags);

			if (ith_tagged == -1) continue;
		
			std::string moduleLabel = (jetTags).provenance()->moduleLabel();

			
			if ( moduleLabel == "trackCountingHighEffJetTags" ) {

				fS8evt->btag_TrkCounting_disc3D_2trk.push_back( (*jetTags)[ith_tagged].discriminator() ); // 2nd trk, 3D
				
				gotTCHE = true;

				int NtrksInJet = (*jetTags)[ith_tagged].tracks().size();
				fS8evt->jet_ntrks.push_back( NtrksInJet );
			}
			else if ( moduleLabel == "trackCountingHighPurJetTags" ) {

				fS8evt->btag_TrkCounting_disc3D_3trk.push_back( (*jetTags)[ith_tagged].discriminator() ); // 3rd trk, 3D

				gotTCHP = true;
			}
			else if ( moduleLabel == "jetProbabilityJetTags" ) {

				fS8evt->btag_JetProb_disc3D.push_back( (*jetTags)[ith_tagged].discriminator());
				
				gotJP = true;

			}
			else if ( moduleLabel == "softMuonJetTags" ) {

				gotSMT = true;

			}

			if (!gotTCHE) fS8evt->btag_TrkCounting_disc3D_2trk.push_back( -9999. );
			
			if (!gotTCHE) fS8evt->btag_TrkCounting_disc3D_3trk.push_back( -9999. );
			
			if (!gotJP)   fS8evt->btag_JetProb_disc3D.push_back( -9999. );
			
		}

		
		/*
		if (ith_tagged != -1 ) {
			fS8evt->btag_TrkCounting_disc3D_1trk.push_back( TrkCountingInfo[ith_tagged].discriminator(1,0) ); // 1st trk, 3D
			fS8evt->btag_TrkCounting_disc3D_2trk.push_back( TrkCountingInfo[ith_tagged].discriminator(2,0) ); // 2nd trk, 3D
			fS8evt->btag_TrkCounting_disc3D_3trk.push_back( TrkCountingInfo[ith_tagged].discriminator(3,0) ); // 3nd trk, 3D
			fS8evt->btag_TrkCounting_disc2D_1trk.push_back( TrkCountingInfo[ith_tagged].discriminator(1,1) ); // 1st trk, 2D
			fS8evt->btag_TrkCounting_disc2D_2trk.push_back( TrkCountingInfo[ith_tagged].discriminator(2,1) ); // 2nd trk, 2D
			fS8evt->btag_TrkCounting_disc2D_2trk.push_back( TrkCountingInfo[ith_tagged].discriminator(3,1) ); // 3nd trk, 2D
			isbtagged = 1;

			int NtrksInJet = TrkCountingInfo[ith_tagged].tracks().size();
			fS8evt->jet_ntrks.push_back( NtrksInJet );
			//negative tags
			int goodtrksInJet2D = 0;
			int goodtrksInJet3D = 0;
			//search all tracks with an IPs > -10
			for(int n=0; n < NtrksInJet ; ++n) {
				if( TrkCountingInfo[ith_tagged].significance(n,1)  > -99.)
					goodtrksInJet2D++;
				if( TrkCountingInfo[ith_tagged].significance(n,0)  > -99.)
					goodtrksInJet3D++;
			}
            //an example of how can we get the 3 tracks with the smallests IPs
			
			if( goodtrksInJet2D > 2 ) {
				fS8evt->btag_NegTag_disc2D_1trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet2D-3,1) ); //2D
				fS8evt->btag_NegTag_disc2D_2trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet2D-2,1) );
				fS8evt->btag_NegTag_disc2D_3trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet2D-1,1) );
			} else {
				fS8evt->btag_NegTag_disc2D_1trk.push_back( -9999. );
				fS8evt->btag_NegTag_disc2D_2trk.push_back( -9999. );
				fS8evt->btag_NegTag_disc2D_3trk.push_back( -9999. );
			}
			
			if( goodtrksInJet3D > 2 ) {
				fS8evt->btag_NegTag_disc3D_1trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-3,0) ); //3D	
				//std::cout << "test -3   " << TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-3,0) << std::endl;

				fS8evt->btag_NegTag_disc3D_2trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-2,0) );
				//std::cout << "test -2   " << TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-2,0) << std::endl;
				fS8evt->btag_NegTag_disc3D_3trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-1,0) );
				//std::cout << "test -1   " << TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-1,0) << std::endl;
				
			} 
			if(goodtrksInJet3D ==2) {
				fS8evt->btag_NegTag_disc3D_1trk.push_back( -9999. );

				fS8evt->btag_NegTag_disc3D_2trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-2,0) );
				//std::cout << "test -2   " << TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-2,0) << std::endl;
				fS8evt->btag_NegTag_disc3D_3trk.push_back(
					TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-1,0) );
				//std::cout << "test -1   " << TrkCountingInfo[ith_tagged].significance(goodtrksInJet3D-1,0) << std::endl;
			}
			if(goodtrksInJet3D<2) {
				fS8evt->btag_NegTag_disc3D_1trk.push_back( -9999. );
				fS8evt->btag_NegTag_disc3D_2trk.push_back( -9999. );
				fS8evt->btag_NegTag_disc3D_3trk.push_back( -9999. );
			}
		} else {
			fS8evt->btag_TrkCounting_disc2D_1trk.push_back( -9999. );
			fS8evt->btag_TrkCounting_disc2D_2trk.push_back( -9999. );
			fS8evt->btag_TrkCounting_disc2D_3trk.push_back( -9999. );
			fS8evt->btag_TrkCounting_disc3D_1trk.push_back( -9999. );
			fS8evt->btag_TrkCounting_disc3D_2trk.push_back( -9999. );
			fS8evt->btag_TrkCounting_disc3D_3trk.push_back( -9999. );
			
			fS8evt->jet_ntrks.push_back( 0 );
			//fS8evt->btag_NegTag_disc2D_1trk.push_back( -9999. );
			//fS8evt->btag_NegTag_disc2D_2trk.push_back( -9999. );
			//fS8evt->btag_NegTag_disc2D_3trk.push_back( -9999. );
			fS8evt->btag_NegTag_disc3D_1trk.push_back( -9999. );
			fS8evt->btag_NegTag_disc3D_2trk.push_back( -9999. );
			fS8evt->btag_NegTag_disc3D_3trk.push_back( -9999. );
		}
   
		ith_tagged = -1;
		ith_tagged = this->TaggedJet(*jet,btagCollTP);
		if (ith_tagged != -1 ) {
		  fS8evt->btag_JetProb_disc3D.push_back(btagCollTP[ith_tagged].discriminator());
		  //fS8evt->btag_JetProb_disc3D.push_back( JetProbInfo[ith_tagged].discriminator(2,0) ); // 3D
			//fS8evt->btag_JetProb_disc2D.push_back( JetProbInfo[ith_tagged].discriminator(2,1) ); // 2D
			isbtagged = 1;
		} else {
			fS8evt->btag_JetProb_disc3D.push_back( -9999. );
			//fS8evt->btag_JetProb_disc2D.push_back( -9999. );
		}
		
		//fS8evt->jet_isbtagged.push_back( isbtagged );
		
		*/
								   
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

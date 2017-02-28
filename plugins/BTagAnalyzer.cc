// -*- C++ -*-
//
// Package:    BTagAnalyzerT
// Class:      BTagAnalyzerT
//
/**\class BTagAnalyzer BTagAnalyzer.cc RecoBTag/PerformanceMeasurements/plugins/BTagAnalyzer.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Andrea Jeremy
//         Created:  Thu Dec 20 10:00:00 CEST 2012
//
//
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "RecoBTag/PerformanceMeasurements/interface/CategoryFinder.h"

// reco track and vertex
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "TGraphErrors.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/TrackHistory/interface/TrackCategories.h"
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CommonTools/Utils/interface/TMVAEvaluator.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// reconstruct IP
#include "TrackingTools/IPTools/interface/IPTools.h"

// Math clusters to TrackingParticles
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoBTag/SecondaryVertex/interface/V0Filter.h"
#include "RecoBTag/ImpactParameter/plugins/IPProducer.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "FWCore/Utilities/interface/RegexMatch.h"
#include <boost/regex.hpp>

#include "RecoBTag/PerformanceMeasurements/interface/JetInfoBranches.h"
#include "RecoBTag/PerformanceMeasurements/interface/EventInfoBranches.h"
#include "RecoBTag/PerformanceMeasurements/interface/BookHistograms.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Common/interface/Provenance.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "fastjet/contrib/Njettiness.hh"

#include "RecoBTag/SecondaryVertex/interface/CombinedSVSoftLeptonComputer.h"

//
// constants, enums and typedefs
//
typedef edm::View<pat::Jet> PatJetCollection;
typedef std::vector<edm::Ptr<pat::Jet> > PatJetPtrCollection;

//
// class declaration
//

struct orderByPt {
  const std::string mCorrLevel;
  orderByPt(const std::string& fCorrLevel) : mCorrLevel(fCorrLevel) {}
  bool operator ()(edm::Ptr<pat::Jet> const& a, edm::Ptr<pat::Jet> const& b) {
    if( mCorrLevel=="Uncorrected" )
      return a->correctedJet("Uncorrected").pt() > b->correctedJet("Uncorrected").pt();
    else
      return a->pt() > b->pt();
  }
};

class simPrimaryVertex {
  public:
    simPrimaryVertex(double x1,double y1,double z1):x(x1),y(y1),z(z1),ptsq(0),nGenTrk(0){};
    double x,y,z;
    HepMC::FourVector ptot;
    //HepLorentzVector ptot;
    double ptsq;
    int nGenTrk;
    std::vector<int> finalstateParticles;
    std::vector<int> simTrackIndex;
    std::vector<int> genVertex;
    const reco::Vertex *recVtx;
};

const reco::TrackBaseRef toTrackRef(const reco::TrackRef & trk) {return reco::TrackBaseRef(trk);}
const reco::TrackBaseRef toTrackRef(const edm::Ptr<reco::Candidate> & cnd)
{
  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(cnd.get());

  if(pcand) // MiniAOD case
    return reco::TrackBaseRef(); // return null reference since no tracks are stored in MiniAOD
  else
  {
    const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(cnd.get());

    if ( (std::abs(pfcand->pdgId()) == 11 || pfcand->pdgId() == 22) && pfcand->gsfTrackRef().isNonnull() && pfcand->gsfTrackRef().isAvailable() )
      return reco::TrackBaseRef(pfcand->gsfTrackRef());
    else if ( pfcand->trackRef().isNonnull() && pfcand->trackRef().isAvailable() )
      return reco::TrackBaseRef(pfcand->trackRef());
    else
      return reco::TrackBaseRef();
  }
}

const math::XYZPoint & position(const reco::Vertex & sv) {return sv.position();}
const math::XYZPoint & position(const reco::VertexCompositePtrCandidate & sv) {return sv.vertex();}
const double xError(const reco::Vertex & sv) {return sv.xError();}
const double xError(const reco::VertexCompositePtrCandidate & sv) {return sqrt(sv.vertexCovariance(0,0));}
const double yError(const reco::Vertex & sv) {return sv.yError();}
const double yError(const reco::VertexCompositePtrCandidate & sv) {return sqrt(sv.vertexCovariance(1,1));}
const double zError(const reco::Vertex & sv) {return sv.zError();}
const double zError(const reco::VertexCompositePtrCandidate & sv) {return sqrt(sv.vertexCovariance(2,2));}
const double chi2(const reco::Vertex & sv) {return sv.chi2();}
const double chi2(const reco::VertexCompositePtrCandidate & sv) {return sv.vertexChi2();}
const double ndof(const reco::Vertex & sv) {return sv.ndof();}
const double ndof(const reco::VertexCompositePtrCandidate & sv) {return sv.vertexNdof();}
const unsigned int vtxTracks(const reco::Vertex & sv) {return sv.nTracks();}
const unsigned int vtxTracks(const reco::VertexCompositePtrCandidate & sv) {return sv.numberOfSourceCandidatePtrs();}


template<typename IPTI,typename VTX>
class BTagAnalyzerT : public edm::EDAnalyzer
{
public:
  explicit BTagAnalyzerT(const edm::ParameterSet&);
  ~BTagAnalyzerT();
  typedef IPTI IPTagInfo;
  typedef typename IPTI::input_container Tracks;
  typedef typename IPTI::input_container::value_type TrackRef;
  typedef VTX Vertex;
  typedef reco::TemplatedSecondaryVertexTagInfo<IPTI,VTX> SVTagInfo;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillHelpers(const edm::Event&);

  const IPTagInfo * toIPTagInfo(const pat::Jet & jet, const std::string & tagInfos);
  const SVTagInfo * toSVTagInfo(const pat::Jet & jet, const std::string & tagInfos);

  const Tracks toAllTracks(const pat::Jet & jet, const std::string & tagInfos, const reco::JetTagInfo & jetTagInfo, const int & iJetColl);

  bool findCat(const reco::Track* ,CategoryFinder& );
  std::vector< float > getTrackProbabilies(std::vector< float > , const int );
  double calculProbability(std::vector< float > );

  float calculPtRel(const reco::Track& theMuon, const pat::Jet& theJet);

  const edm::Ptr<reco::Muon> matchMuon(const edm::Ptr<reco::Candidate>& theMuon, const edm::View<reco::Muon>& muons);

  void setTracksPVBase(const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);
  void setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);

  void setTracksSV(const TrackRef & trackRef, const SVTagInfo *, int & isFromSV, int & iSV, float & SVweight);
  void vertexKinematicsAndChange(const Vertex & vertex, reco::TrackKinematics & vertexKinematics, Int_t & charge);
  
  void etaRelToTauAxis(const Vertex & vertex, const fastjet::PseudoJet & tauAxis, std::vector<float> & tau_trackEtaRel);

  bool NameCompatible(const std::string& pattern, const std::string& name);

  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::HepMCProduct>& evtMC);

  void processTrig(const edm::Handle<edm::TriggerResults>&, const std::vector<std::string>&) ;

  void processJets(const edm::Handle<PatJetCollection>&, const edm::Handle<PatJetCollection>&,
		   const std::vector<edm::Handle<PatJetCollection> >&,
		   const edm::Event&, const edm::EventSetup&, const int) ;

  void recalcNsubjettiness(const pat::Jet & jet, float & tau1, float & tau2, std::vector<fastjet::PseudoJet> & currentAxes) const;

  int isFromGSP(const reco::Candidate* c);

  bool isHardProcess(const int status);

  void matchGroomedJets(const edm::Handle<PatJetCollection>& jets,
			const edm::Handle<PatJetCollection>& matchedJets,
			std::vector<int>& matchedIndices);

  // ----------member data ---------------------------
  std::string outputFile_;
  //std::vector< std::string > moduleLabel_;

  bool runFatJets_ ;
  bool runSubJets_ ;
  bool allowJetSkipping_ ;

  edm::EDGetTokenT<GenEventInfoProduct> src_;  // Generator/handronizer module label
  edm::EDGetTokenT<edm::View<reco::Muon>> muonCollectionName_;
  edm::EDGetTokenT<std::vector<pat::Muon>> patMuonCollectionName_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionName_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticleCollectionName_;
  edm::EDGetTokenT<edm::TriggerResults> triggerTable_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> putoken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> putokenmini;
  edm::EDGetTokenT<int> ttbartop;
  edm::EDGetTokenT<int> ttbartoptrig;
  edm::EDGetTokenT<int> metfilterIntoken;

  std::string   branchNamePrefix_;
  edm::EDGetTokenT<PatJetCollection> JetCollectionTag_;
  std::vector<edm::EDGetTokenT<PatJetCollection>> SubJetCollectionTags_;
  std::vector<std::string> SubJetLabels_;

  std::string jetPBJetTags_;
  std::string jetPNegBJetTags_;
  std::string jetPPosBJetTags_;

  std::string jetBPBJetTags_;
  std::string jetBPNegBJetTags_;
  std::string jetBPPosBJetTags_;

  std::string trackCHEBJetTags_;
  std::string trackCNegHEBJetTags_;

  std::string trackCHPBJetTags_;
  std::string trackCNegHPBJetTags_;

  std::string combinedSVBJetTags_;
  std::string combinedSVNegBJetTags_;
  std::string combinedSVPosBJetTags_;

	std::string deepCSVBJetTags_;
	std::string deepCSVNegBJetTags_;
	std::string deepCSVPosBJetTags_;

  std::string combinedIVFSVBJetTags_;
  std::string combinedIVFSVPosBJetTags_;
  std::string combinedIVFSVNegBJetTags_;

  std::string simpleSVHighEffBJetTags_;
  std::string simpleSVNegHighEffBJetTags_;
  std::string simpleSVHighPurBJetTags_;
  std::string simpleSVNegHighPurBJetTags_;

  std::string softPFMuonBJetTags_;
  std::string softPFMuonNegBJetTags_;
  std::string softPFMuonPosBJetTags_;

  std::string softPFElectronBJetTags_;
  std::string softPFElectronNegBJetTags_;
  std::string softPFElectronPosBJetTags_;

  std::string doubleSVBJetTags_;

  std::string cMVABJetTags_;
  std::string cMVAv2BJetTags_;
  std::string cMVAv2NegBJetTags_;
  std::string cMVAv2PosBJetTags_;

  std::string ipTagInfos_;
  std::string svTagInfos_;
  std::string svNegTagInfos_;
  std::string softPFMuonTagInfos_;
  std::string softPFElectronTagInfos_;
  std::string bdsvTagInfos_;

  edm::EDGetTokenT<reco::VertexCollection> primaryVertexColl_;
  edm::EDGetTokenT<reco::TrackCollection> tracksColl_;

  std::string   SVComputer_;
  std::string   SVComputerSubJets_;

  std::string   ipTagInfosCTag_;
  std::string   svTagInfosCTag_;
  std::string   svNegTagInfosCTag_;
  std::string   softPFMuonTagInfosCTag_;
  std::string   softPFElectronTagInfosCTag_;
  std::string   SLComputer_;
  std::string   CvsBCJetTags_;
  std::string   CvsBNegCJetTags_;
  std::string   CvsBPosCJetTags_;
  std::string   CvsLCJetTags_;
  std::string   CvsLNegCJetTags_;
  std::string   CvsLPosCJetTags_;

  bool useTrackHistory_;
  TFile*  rootFile_;
  double minJetPt_;
  double maxJetEta_;

  int selTagger_;
  bool isData_;
  bool useSelectedTracks_;
  bool fillsvTagInfo_;
  bool fillPU_;
  bool fillQuarks_;
  bool fillGenPruned_;
  bool produceJetTrackTree_;
  bool produceJetTrackTruthTree_;
  bool produceAllTrackTree_;
  bool producePtRelTemplate_;
  bool storeEventInfo_;
  bool storePatMuons_;
  bool storeTagVariables_;
  bool storeTagVariablesSubJets_;
  bool storeCSVTagVariables_;
  bool storeCSVTagVariablesSubJets_;

  bool storeCTagVariables_;
  bool doCTag_; 

  bool use_ttbar_filter_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > ttbarproducerGen_;
  edm::EDGetTokenT<edm::View<pat::Electron>> ttbarproducerEle_;
  edm::EDGetTokenT<edm::View<pat::Muon>> ttbarproducerMuon_;
  edm::EDGetTokenT<edm::View<pat::MET>> ttbarproducerMET_;
  edm::EDGetTokenT<GenEventInfoProduct> generator;
  edm::EDGetTokenT<edm::HepMCProduct> generatorhep;
  edm::EDGetTokenT<GenEventInfoProduct> generatorevt;
  edm::EDGetTokenT<LHEEventProduct> generatorlhe;
  
  // Matching clusters to TP token
  edm::EDGetTokenT<ClusterTPAssociation> clusterTPMapToken_;

  // trigger list
  std::vector<std::string> triggerPathNames_;

  //rho
  edm::EDGetTokenT<double> rhoTag_;

  TrackClassifier classifier_;

  edm::Service<TFileService> fs;

  CategoryFinder cat0;
  CategoryFinder cat1;
  CategoryFinder cat2;
  CategoryFinder cat3;
  CategoryFinder cat4;
  CategoryFinder cat5;
  CategoryFinder cat6;
  CategoryFinder cat7;
  CategoryFinder cat8;
  CategoryFinder cat9;

  ///////////////
  // Ntuple info

  TTree *smalltree;

  //// Event info
  EventInfoBranches EventInfo;

  //// Jet info
  std::vector<JetInfoBranches> JetInfo;
  std::map<std::string, SubJetInfoBranches> SubJetInfo;
  std::vector<BookHistograms*> Histos;

  const  reco::Vertex  *pv;

  bool PFJet80 ;
  std::vector<std::string> PFJet80TriggerPathNames_;

  const GenericMVAJetTagComputer *computer ;

  const GenericMVAJetTagComputer *slcomputer ;

  edm::View<reco::Muon> muons ;

  edm::ESHandle<TransientTrackBuilder> trackBuilder ;
  edm::Handle<reco::VertexCollection> primaryVertex ;

  int cap0, cap1, cap2, cap3, cap4, cap5, cap6, cap7, cap8;
  int can0, can1, can2, can3, can4, can5, can6, can7, can8;

  // Generator/hadronizer type (information stored bitwise)
  unsigned int hadronizerType_;

  // PF jet ID
  PFJetIDSelectionFunctor pfjetIDLoose_;
  PFJetIDSelectionFunctor pfjetIDTight_;

  // helper class for associating PF candidates to jets
  IPProducerHelpers::FromJetAndCands m_helper;
  std::vector<IPProducerHelpers::FromJetAndCands> m_subjetHelper;

  const double beta_;
  const double R0_;

  const double maxSVDeltaRToJet_;

  const edm::FileInPath weightFile_;

  // MVA evaluators
  std::unique_ptr<TMVAEvaluator> evaluator_SV_;
  
  // track V0 filter
  reco::V0Filter trackPairV0Filter;

  // track cuts
  double distJetAxis_;
  double decayLength_;
  double deltaR_;
  double distJetAxisSubJets_;
  double decayLengthSubJets_;
  double deltaRSubJets_;
};

template<typename IPTI,typename VTX>
BTagAnalyzerT<IPTI,VTX>::BTagAnalyzerT(const edm::ParameterSet& iConfig):
  classifier_(iConfig, consumesCollector()) ,
  pv(0),
  PFJet80(0),
  computer(0),
  slcomputer(0),
  cap0(0),
  cap1(0),
  cap2(0),
  cap3(0),
  cap4(0),
  cap5(0),
  cap6(0),
  cap7(0),
  cap8(0),
  can0(0),
  can1(0),
  can2(0),
  can3(0),
  can4(0),
  can5(0),
  can6(0),
  can7(0),
  can8(0),
  hadronizerType_(0),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),
  pfjetIDTight_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT ),
  m_helper(iConfig, consumesCollector(),"Jets"),
  beta_(iConfig.getParameter<double>("beta")),
  R0_(iConfig.getParameter<double>("R0")),
  maxSVDeltaRToJet_(iConfig.getParameter<double>("maxSVDeltaRToJet")),
  weightFile_(iConfig.getParameter<edm::FileInPath>("weightFile")),
  trackPairV0Filter(iConfig.getParameter<edm::ParameterSet>("trackPairV0Filter")),
  distJetAxis_(iConfig.getParameter<double>("distJetAxis")),
  decayLength_(iConfig.getParameter<double>("decayLength")),
  deltaR_(iConfig.getParameter<double>("deltaR")),
  distJetAxisSubJets_(iConfig.getParameter<double>("distJetAxisSubJets")),
  decayLengthSubJets_(iConfig.getParameter<double>("decayLengthSubJets")),
  deltaRSubJets_(iConfig.getParameter<double>("deltaRSubJets"))
{
  //now do what ever initialization you need
  std::string module_type  = iConfig.getParameter<std::string>("@module_type");
  std::string module_label = iConfig.getParameter<std::string>("@module_label");
  std::cout << "Constructing " << module_type << ":" << module_label << std::endl;

  // Parameters
  runFatJets_ = iConfig.getParameter<bool>("runFatJets");
  runSubJets_ = iConfig.getParameter<bool>("runSubJets");
  allowJetSkipping_ = iConfig.getParameter<bool>("allowJetSkipping");

  minJetPt_  = iConfig.getParameter<double>("MinPt");
  maxJetEta_ = iConfig.getParameter<double>("MaxEta");

  selTagger_ = iConfig.getParameter<int>("selTagger");

  useSelectedTracks_    = iConfig.getParameter<bool> ("useSelectedTracks");
  fillsvTagInfo_    = iConfig.getParameter<bool> ("fillsvTagInfo");
  fillPU_    = iConfig.getParameter<bool> ("fillPU");
  fillQuarks_    = iConfig.getParameter<bool> ("fillQuarks");
  fillGenPruned_    = iConfig.getParameter<bool> ("fillGenPruned");
  useTrackHistory_      = iConfig.getParameter<bool> ("useTrackHistory");
  produceJetTrackTree_  = iConfig.getParameter<bool> ("produceJetTrackTree");
  produceJetTrackTruthTree_  = iConfig.getParameter<bool> ("produceJetTrackTruthTree");
  produceAllTrackTree_  = iConfig.getParameter<bool> ("produceAllTrackTree");
  producePtRelTemplate_ = iConfig.getParameter<bool> ("producePtRelTemplate");
  storeEventInfo_ = iConfig.getParameter<bool>("storeEventInfo");
  storePatMuons_ = iConfig.getParameter<bool>("storePatMuons");
  storeTagVariables_ = iConfig.getParameter<bool>("storeTagVariables");
  storeTagVariablesSubJets_ = iConfig.getParameter<bool>("storeTagVariablesSubJets");
  storeCSVTagVariables_ = iConfig.getParameter<bool>("storeCSVTagVariables");
  storeCSVTagVariablesSubJets_ = iConfig.getParameter<bool>("storeCSVTagVariablesSubJets");

  storeCTagVariables_ = iConfig.getParameter<bool>("storeCTagVariables");
  doCTag_             = iConfig.getParameter<bool>("doCTag");

  use_ttbar_filter_ = iConfig.getParameter<bool> ("use_ttbar_filter");
  ttbarproducerGen_ = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("ttbarproducer")),
  ttbarproducerEle_ = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("ttbarproducer"));
  ttbarproducerMuon_ = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("ttbarproducer"));
  ttbarproducerMET_ = consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("ttbarproducer"));
  rhoTag_               = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));

  // Matching Clusters to TrackingParticles consume token
  if( produceJetTrackTruthTree_ )
    clusterTPMapToken_ = consumes<ClusterTPAssociation>(iConfig.getParameter<edm::InputTag>("clusterTPMap"));

  // Modules
  primaryVertexColl_   = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexColl"));
  if( produceAllTrackTree_ && storeEventInfo_ )
    tracksColl_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksColl"));

  src_ = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("src"));
  generator = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  generatorhep = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared"));
  generatorlhe = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""));
  generatorevt = consumes<GenEventInfoProduct>(edm::InputTag("generator",""));

  branchNamePrefix_ = iConfig.getParameter<std::string>("BranchNamePrefix");
  JetCollectionTag_ = consumes<PatJetCollection>(iConfig.getParameter<edm::InputTag>("Jets"));
  SubJetLabels_         = iConfig.getParameter<std::vector<std::string> >("SubJetLabels");
  for(const auto& tag: iConfig.getParameter<std::vector<edm::InputTag> >("SubJets")){
    SubJetCollectionTags_.push_back(consumes<PatJetCollection>(tag));
  }
  SubJetLabels_         = iConfig.getParameter<std::vector<std::string> >("SubJetLabels");

  if ( runFatJets_ )
  {
    // initialize MVA evaluators
    evaluator_SV_.reset( new TMVAEvaluator() );
    std::vector<std::string> variables({"z_ratio",
                                        "trackSipdSig_3","trackSipdSig_2","trackSipdSig_1","trackSipdSig_0",
                                        "trackSipdSig_1_0","trackSipdSig_0_0","trackSipdSig_1_1","trackSipdSig_0_1",
                                        "trackSip2dSigAboveCharm_0","trackSip2dSigAboveBottom_0","trackSip2dSigAboveBottom_1",
                                        "tau0_trackEtaRel_0","tau0_trackEtaRel_1","tau0_trackEtaRel_2",
                                        "tau1_trackEtaRel_0","tau1_trackEtaRel_1","tau1_trackEtaRel_2",
                                        "tau_vertexMass_0","tau_vertexEnergyRatio_0","tau_vertexDeltaR_0","tau_flightDistance2dSig_0",
                                        "tau_vertexMass_1","tau_vertexEnergyRatio_1","tau_flightDistance2dSig_1",
                                        "jetNTracks","nSV"});
    // book TMVA readers
    std::vector<std::string> spectators({"massPruned", "flavour", "nbHadrons", "ptPruned", "etaPruned"});
    
    bool useGBRForest = true;
    evaluator_SV_->initialize("Color:Silent:Error", "BDTG", weightFile_.fullPath(), variables, spectators, useGBRForest);
  }

  if( runSubJets_ && ( SubJetCollectionTags_.size()>0 || SubJetLabels_.size()>0 ) )
  {
    if( SubJetCollectionTags_.size()!=SubJetLabels_.size() )
    {
      edm::LogWarning("SubjetLabelMismatch") << "The number of subjet labels does not match the number of subjet collections. Subjets will not be stored.";
      runSubJets_ = false;
    }
  }

  trackCHEBJetTags_    = iConfig.getParameter<std::string>("trackCHEBJetTags");
  trackCNegHEBJetTags_ = iConfig.getParameter<std::string>("trackCNegHEBJetTags");

  trackCHPBJetTags_    = iConfig.getParameter<std::string>("trackCHPBJetTags");
  trackCNegHPBJetTags_ = iConfig.getParameter<std::string>("trackCNegHPBJetTags");

  jetPBJetTags_        = iConfig.getParameter<std::string>("jetPBJetTags");
  jetPNegBJetTags_     = iConfig.getParameter<std::string>("jetPNegBJetTags");
  jetPPosBJetTags_     = iConfig.getParameter<std::string>("jetPPosBJetTags");

  jetBPBJetTags_        = iConfig.getParameter<std::string>("jetBPBJetTags");
  jetBPNegBJetTags_     = iConfig.getParameter<std::string>("jetBPNegBJetTags");
  jetBPPosBJetTags_     = iConfig.getParameter<std::string>("jetBPPosBJetTags");

  simpleSVHighEffBJetTags_      = iConfig.getParameter<std::string>("simpleSVHighEffBJetTags");
  simpleSVNegHighEffBJetTags_   = iConfig.getParameter<std::string>("simpleSVNegHighEffBJetTags");
  simpleSVHighPurBJetTags_      = iConfig.getParameter<std::string>("simpleSVHighPurBJetTags");
  simpleSVNegHighPurBJetTags_   = iConfig.getParameter<std::string>("simpleSVNegHighPurBJetTags");

  combinedSVBJetTags_     = iConfig.getParameter<std::string>("combinedSVBJetTags");
  combinedSVNegBJetTags_  = iConfig.getParameter<std::string>("combinedSVNegBJetTags");
  combinedSVPosBJetTags_  = iConfig.getParameter<std::string>("combinedSVPosBJetTags");

	deepCSVBJetTags_    = iConfig.getParameter<std::string>("deepCSVBJetTags");
	deepCSVNegBJetTags_ = iConfig.getParameter<std::string>("deepCSVNegBJetTags");
	deepCSVPosBJetTags_ = iConfig.getParameter<std::string>("deepCSVPosBJetTags");

  combinedIVFSVBJetTags_      = iConfig.getParameter<std::string>("combinedIVFSVBJetTags");
  combinedIVFSVPosBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVPosBJetTags");
  combinedIVFSVNegBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVNegBJetTags");

  softPFMuonBJetTags_       = iConfig.getParameter<std::string>("softPFMuonBJetTags");
  softPFMuonNegBJetTags_    = iConfig.getParameter<std::string>("softPFMuonNegBJetTags");
  softPFMuonPosBJetTags_    = iConfig.getParameter<std::string>("softPFMuonPosBJetTags");

  softPFElectronBJetTags_       = iConfig.getParameter<std::string>("softPFElectronBJetTags");
  softPFElectronNegBJetTags_    = iConfig.getParameter<std::string>("softPFElectronNegBJetTags");
  softPFElectronPosBJetTags_    = iConfig.getParameter<std::string>("softPFElectronPosBJetTags");

  doubleSVBJetTags_ = iConfig.getParameter<std::string>("doubleSVBJetTags");

  cMVABJetTags_ = iConfig.getParameter<std::string>("cMVABJetTags");
  cMVAv2BJetTags_ = iConfig.getParameter<std::string>("cMVAv2BJetTags");
  cMVAv2NegBJetTags_ = iConfig.getParameter<std::string>("cMVAv2NegBJetTags");
  cMVAv2PosBJetTags_ = iConfig.getParameter<std::string>("cMVAv2PosBJetTags");

  ipTagInfos_              = iConfig.getParameter<std::string>("ipTagInfos");
  svTagInfos_              = iConfig.getParameter<std::string>("svTagInfos");
  svNegTagInfos_           = iConfig.getParameter<std::string>("svNegTagInfos");
  softPFMuonTagInfos_      = iConfig.getParameter<std::string>("softPFMuonTagInfos");
  softPFElectronTagInfos_  = iConfig.getParameter<std::string>("softPFElectronTagInfos");
  bdsvTagInfos_            = iConfig.getParameter<std::string>("bdsvTagInfos");

  muonCollectionName_       = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollectionName"));
  patMuonCollectionName_    = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("patMuonCollectionName"));
  genParticleCollectionName_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  prunedGenParticleCollectionName_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"));

  triggerTable_             = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerTable"));

  SVComputer_               = iConfig.getParameter<std::string>("svComputer");
  SVComputerSubJets_        = iConfig.getParameter<std::string>("svComputerSubJets");

  ipTagInfosCTag_              = iConfig.getParameter<std::string>("ipTagInfosCTag");
  svTagInfosCTag_              = iConfig.getParameter<std::string>("svTagInfosCTag");
  svNegTagInfosCTag_              = iConfig.getParameter<std::string>("svNegTagInfosCTag");
  softPFMuonTagInfosCTag_      = iConfig.getParameter<std::string>("softPFMuonTagInfosCTag");
  softPFElectronTagInfosCTag_  = iConfig.getParameter<std::string>("softPFElectronTagInfosCTag");
  SLComputer_               = iConfig.getParameter<std::string>("slComputer");
  CvsBCJetTags_             = iConfig.getParameter<std::string>("CvsBCJetTags");
  CvsBNegCJetTags_             = iConfig.getParameter<std::string>("CvsBNegCJetTags");
  CvsBPosCJetTags_             = iConfig.getParameter<std::string>("CvsBPosCJetTags");
  CvsLCJetTags_             = iConfig.getParameter<std::string>("CvsLCJetTags");
  CvsLNegCJetTags_             = iConfig.getParameter<std::string>("CvsLNegCJetTags");
  CvsLPosCJetTags_             = iConfig.getParameter<std::string>("CvsLPosCJetTags");

  triggerPathNames_        = iConfig.getParameter<std::vector<std::string> >("TriggerPathNames");
  PFJet80TriggerPathNames_ = iConfig.getParameter<std::vector<std::string> >("PFJet80TriggerPathNames");

  putoken = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
  putokenmini = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
  ttbartop = consumes<int>(edm::InputTag("ttbarselectionproducer:topChannel"));
  ttbartoptrig = consumes<int>(edm::InputTag("ttbarselectionproducer:topTrigger"));
  metfilterIntoken = consumes<int>(edm::InputTag("ttbarselectionproducer:topMETFilter"));

  ///////////////
  // TTree

  smalltree = fs->make<TTree>("ttree", "ttree");

  //--------------------------------------
  // event information
  //--------------------------------------
  if( storeEventInfo_ )
  {
    EventInfo.RegisterTree(smalltree);
    if ( use_ttbar_filter_ )    EventInfo.RegisterTTbarTree(smalltree);
    if ( produceJetTrackTree_ ) EventInfo.RegisterJetTrackTree(smalltree);
    if ( produceAllTrackTree_ ) EventInfo.RegisterAllTrackTree(smalltree);
    if ( storePatMuons_ )       EventInfo.RegisterPatMuonTree(smalltree);
  }

  //--------------------------------------
  // jet information
  //--------------------------------------
  JetInfo.reserve(1+(runSubJets_ ? SubJetLabels_.size() : 0));
  JetInfo[0].RegisterTree(smalltree,branchNamePrefix_);
  if ( runFatJets_ )          JetInfo[0].RegisterFatJetSpecificTree(smalltree,branchNamePrefix_,produceJetTrackTree_);
  if ( produceJetTrackTree_ ) JetInfo[0].RegisterJetTrackTree(smalltree,branchNamePrefix_);
  if ( produceJetTrackTruthTree_ ) JetInfo[0].RegisterJetTrackTruthTree(smalltree,branchNamePrefix_);
  if ( producePtRelTemplate_ ) JetInfo[0].RegisterJetTrackIncTree(smalltree,branchNamePrefix_);
  if ( fillsvTagInfo_ )       JetInfo[0].RegisterJetSVTree(smalltree,branchNamePrefix_);
  if ( storeTagVariables_)    JetInfo[0].RegisterTagVarTree(smalltree,branchNamePrefix_);
  if ( storeCSVTagVariables_) JetInfo[0].RegisterCSVTagVarTree(smalltree,branchNamePrefix_);
  if ( storeCTagVariables_) JetInfo[0].RegisterCTagVarTree(smalltree,branchNamePrefix_);
  if ( runSubJets_ ) {
    for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
    {
      if ( runFatJets_ ) SubJetInfo[SubJetLabels_[i]].RegisterTree(smalltree,branchNamePrefix_,SubJetLabels_[i]);
      JetInfo[1+i].RegisterTree(smalltree,SubJetLabels_[i]+"SubJetInfo");
      JetInfo[1+i].RegisterSubJetSpecificTree(smalltree,SubJetLabels_[i]+"SubJetInfo");
      if ( produceJetTrackTree_ )         JetInfo[1+i].RegisterJetTrackTree(smalltree,SubJetLabels_[i]+"SubJetInfo");
      if ( producePtRelTemplate_ )        JetInfo[1+i].RegisterJetTrackIncTree(smalltree,SubJetLabels_[i]+"SubJetInfo");
      if ( fillsvTagInfo_ )               JetInfo[1+i].RegisterJetSVTree(smalltree,SubJetLabels_[i]+"SubJetInfo");
      if ( storeTagVariablesSubJets_ )    JetInfo[1+i].RegisterTagVarTree(smalltree,SubJetLabels_[i]+"SubJetInfo");
      if ( storeCSVTagVariablesSubJets_ ) JetInfo[1+i].RegisterCSVTagVarTree(smalltree,SubJetLabels_[i]+"SubJetInfo");

      // set up subjet helper classes
      edm::ParameterSet iConfigCopy(iConfig);
      const auto subjettags = iConfig.getParameter<std::vector<edm::InputTag> >("SubJets");
      iConfigCopy.insert( true, "Jets",        edm::Entry("Jets",        subjettags[i],                        true) );
      iConfigCopy.insert( true, "maxDeltaR",   edm::Entry("maxDeltaR",   iConfig.getParameter<double>("subJetMaxDeltaR"), true) );
      iConfigCopy.insert( true, "explicitJTA", edm::Entry("explicitJTA", iConfig.getParameter<bool>("subJetExplicitJTA"), true) );

      m_subjetHelper.push_back(IPProducerHelpers::FromJetAndCands(iConfigCopy, consumesCollector(), "Jets"));
    }
  }

  //// Book Histograms
  Histos.resize(1+(runSubJets_ ? SubJetLabels_.size() : 0));
  Histos[0] = new BookHistograms(fs->mkdir( "HistJets" )) ;
  if (runSubJets_)
  {
    for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
    {
      Histos[1+i] = new BookHistograms(fs->mkdir( "Hist"+SubJetLabels_[i]+"SubJets" )) ;
    }
  }

  std::cout << module_type << ":" << module_label << " constructed" << std::endl;
}

template<typename IPTI,typename VTX>
BTagAnalyzerT<IPTI,VTX>::~BTagAnalyzerT()
{
}

static std::vector<std::size_t> sortedIndexes(const std::vector<reco::btag::TrackIPData >& values)
{
  std::multimap<float, std::size_t> sortedIdx;
  std::vector<std::size_t> result;

  for(size_t i = 0; i < values.size(); ++i)
    sortedIdx.insert( std::pair<float, std::size_t>(values[i].ip3d.significance() , i) );

  for(std::multimap<float, std::size_t>::reverse_iterator it = sortedIdx.rbegin(); it != sortedIdx.rend(); ++it)
    result.push_back(it->second);

  return result;
}


//
// member functions
//

// ------------ method called to for each event  ------------
template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  //------------------------------------------------------
  // Event information
  //------------------------------------------------------

  EventInfo.Run = iEvent.id().run();

  EventInfo.BX  = iEvent.bunchCrossing();

  isData_ = iEvent.isRealData();

  if ( !isData_ && EventInfo.Run > 0 ) EventInfo.Run = -EventInfo.Run;

  EventInfo.Evt  = iEvent.id().event();
  EventInfo.LumiBlock  = iEvent.luminosityBlock();

  // Tag Jets
  if ( useTrackHistory_ ) classifier_.newEvent(iEvent, iSetup);

  // fill the helper classes
  if( !useSelectedTracks_ )
    fillHelpers(iEvent);

  edm::Handle <PatJetCollection> jetsColl;
  iEvent.getByToken (JetCollectionTag_, jetsColl);

  std::vector<edm::Handle<PatJetCollection> > subjetColls;
  if (runSubJets_)
  {
    subjetColls.resize( SubJetLabels_.size() );
    for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
      iEvent.getByToken(SubJetCollectionTags_[i], subjetColls[i]);
  }

  //------------------------------------------------------
  // Determine hadronizer type (done only once per job)
  //------------------------------------------------------
  if( !isData_ && hadronizerType_ == 0 )
  {
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByToken(src_, genEvtInfoProduct);

    std::string moduleName = "";
    const edm::Provenance& prov = iEvent.getProvenance(genEvtInfoProduct.id());
    if( genEvtInfoProduct.isValid() )
      moduleName = edm::moduleName(prov);

    if( moduleName.find("Pythia8")!=std::string::npos )
      hadronizerType_ |= ( 1 << 1 ); // set the 2nd bit
    else // assuming Pythia6
      hadronizerType_ |= ( 1 << 0 ); // set the 1st bit
  }

  //------------------------------------------------------
  // MC informations
  //------------------------------------------------------
  EventInfo.pthat      = -1.;
  EventInfo.nPUtrue    = -1.;
  EventInfo.nPU        = 0;
  EventInfo.ncQuarks   = 0;
  EventInfo.nbQuarks   = 0;
  EventInfo.nBHadrons  = 0;
  EventInfo.nDHadrons  = 0;
  EventInfo.nDaughters = 0;
  EventInfo.nGenV0     = 0;
  EventInfo.nGenlep    = 0;
  EventInfo.nGenquark  = 0;
  EventInfo.nGenPruned = 0;
  EventInfo.nTrkAll    = 0;
  EventInfo.nPatMuon   = 0;
  EventInfo.mcweight   = 1.;

  bool AreBHadrons = false;


  //---------------------------- Start MC info ---------------------------------------//
  if ( !isData_ && storeEventInfo_ ) {
    // EventInfo.pthat
    edm::Handle<GenEventInfoProduct> geninfos;
    iEvent.getByToken( generator,geninfos );
    EventInfo.mcweight=geninfos->weight();
    if (geninfos->binningValues().size()>0) EventInfo.pthat = geninfos->binningValues()[0];
    // pileup

    edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
    edm::Handle<std::vector <PileupSummaryInfo> > PupInfotest;
    bool checkPUname = iEvent.getByToken(putoken, PupInfo);

    if (checkPUname){
      iEvent.getByToken(putoken, PupInfo);
    }
    else{
      iEvent.getByToken(putokenmini, PupInfo);
    }

    std::vector<PileupSummaryInfo>::const_iterator ipu;
    for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
      if(fillPU_){
	for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
	  EventInfo.PU_bunch[EventInfo.nPU]      =  ipu->getBunchCrossing();
	  EventInfo.PU_z[EventInfo.nPU]          = (ipu->getPU_zpositions())[i];
	  EventInfo.PU_sumpT_low[EventInfo.nPU]  = (ipu->getPU_sumpT_lowpT())[i];
	  EventInfo.PU_sumpT_high[EventInfo.nPU] = (ipu->getPU_sumpT_highpT())[i];
	  EventInfo.PU_ntrks_low[EventInfo.nPU]  = (ipu->getPU_ntrks_lowpT())[i];
	  EventInfo.PU_ntrks_high[EventInfo.nPU] = (ipu->getPU_ntrks_highpT())[i];
	  ++EventInfo.nPU;
	}
      }
      EventInfo.nPUtrue = ipu->getTrueNumInteractions();
      if(fillPU_){
	if(EventInfo.nPU==0) EventInfo.nPU = ipu->getPU_NumInteractions(); // needed in case getPU_zpositions() is empty
      }
    }

  //------------------------------------------------------
  // generated particles
  //------------------------------------------------------
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticleCollectionName_, genParticles);

    edm::Handle<reco::GenParticleCollection> prunedGenParticles;
    iEvent.getByToken(prunedGenParticleCollectionName_, prunedGenParticles);

    if(fillGenPruned_){
      //loop over pruned GenParticles to fill branches for MC hard process particles and muons
      for(size_t i = 0; i < prunedGenParticles->size(); ++i){
	const GenParticle & iGenPart = (*prunedGenParticles)[i];
	int status = iGenPart.status();
	int pdgid = iGenPart.pdgId();
	int numMothers = iGenPart.numberOfMothers();

	//fill all the branches
	EventInfo.GenPruned_pT[EventInfo.nGenPruned] = iGenPart.pt();
	EventInfo.GenPruned_eta[EventInfo.nGenPruned] = iGenPart.eta();
	EventInfo.GenPruned_phi[EventInfo.nGenPruned] = iGenPart.phi();
	EventInfo.GenPruned_mass[EventInfo.nGenPruned] = iGenPart.mass();
	EventInfo.GenPruned_status[EventInfo.nGenPruned] = status;
	EventInfo.GenPruned_pdgID[EventInfo.nGenPruned] = pdgid;
	// if no mothers, set mother index to -1 (just so it's not >=0)
	if (numMothers == 0)
	  EventInfo.GenPruned_mother[EventInfo.nGenPruned] = -1;
	else{
	  //something new to distinguish from the no mothers case
	  int idx = -100;
	  //loop over the pruned genparticle list to get the mother's index
	  for( reco::GenParticleCollection::const_iterator mit = prunedGenParticles->begin(); mit != prunedGenParticles->end(); ++mit ) {
	    if( iGenPart.mother(0)==&(*mit) ) {
	      idx = std::distance(prunedGenParticles->begin(),mit);
	      break;
	    }
	  }
	  EventInfo.GenPruned_mother[EventInfo.nGenPruned] = idx;
	}
	++EventInfo.nGenPruned;
      } //end loop over pruned genparticles
    }//fillGenPruned_

    EventInfo.GenPVz = -1000.;

    for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      int ID = abs(genIt.pdgId());
      unsigned int nDaughters = genIt.numberOfDaughters();

      if ( isHardProcess(genIt.status()) ) EventInfo.GenPVz = genIt.vz();

      // prompt parton
      if ( isHardProcess(genIt.status()) && genIt.p4().pt() > 10. ) {
        EventInfo.Genquark_pT[EventInfo.nGenquark]     = genIt.p4().pt();
        EventInfo.Genquark_eta[EventInfo.nGenquark]    = genIt.p4().eta();
        EventInfo.Genquark_phi[EventInfo.nGenquark]    = genIt.p4().phi();
        EventInfo.Genquark_pdgID[EventInfo.nGenquark]  = genIt.pdgId();
        if (genIt.numberOfMothers()>0) EventInfo.Genquark_mother[EventInfo.nGenquark] = genIt.mother()->pdgId();   // protection
        else EventInfo.Genquark_mother[EventInfo.nGenquark] = -99999;
        //  cout << " status " << genIt.status() << " pdgId " << genIt.pdgId()  << " pT " << genIt.p4().pt() << " mother " << moth->pdgId() << endl;
        ++EventInfo.nGenquark;
      }

      // b and c quarks from the end of parton showering and before hadronization
      if(fillQuarks_){
	if ( ID == 4 || ID == 5 ) {
	  if( nDaughters > 0 ) {
	    int nparton_daughters = 0;
	    for (unsigned int d=0; d<nDaughters; ++d) {
	      int daughterID = abs(genIt.daughter(d)->pdgId());
	      if( (daughterID == 1 || daughterID == 2 || daughterID == 3 ||
		   daughterID == 4 || daughterID == 5 || daughterID == 6 || daughterID == 21))
		nparton_daughters++;
	    }
	    if( nparton_daughters == 0 ) {
	      if ( ID == 5 ) {
		EventInfo.bQuark_pT[EventInfo.nbQuarks]  = genIt.p4().pt();
		EventInfo.bQuark_eta[EventInfo.nbQuarks] = genIt.p4().eta();
		EventInfo.bQuark_phi[EventInfo.nbQuarks] = genIt.p4().phi();
		EventInfo.bQuark_pdgID[EventInfo.nbQuarks] = genIt.pdgId();
		EventInfo.bQuark_status[EventInfo.nbQuarks] = genIt.status();
		EventInfo.bQuark_fromGSP[EventInfo.nbQuarks] = isFromGSP(&genIt);
		++EventInfo.nbQuarks;
	      }
	      if ( ID == 4 ) {
		EventInfo.cQuark_pT[EventInfo.ncQuarks]  = genIt.p4().pt();
		EventInfo.cQuark_eta[EventInfo.ncQuarks] = genIt.p4().eta();
		EventInfo.cQuark_phi[EventInfo.ncQuarks] = genIt.p4().phi();
		EventInfo.cQuark_pdgID[EventInfo.ncQuarks] = genIt.pdgId();
		EventInfo.cQuark_status[EventInfo.ncQuarks] = genIt.status();
		EventInfo.cQuark_fromGSP[EventInfo.ncQuarks] = isFromGSP(&genIt);
		++EventInfo.ncQuarks;
	      }
	    }
	  }
	}
      }//fillQuarks_

      if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) AreBHadrons = true;
//       // Primary b Hadrons
//       if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
//         //  cout << " pdgId " << genIt.pdgId()  << " pT " << genIt.p4().pt() << " mother " << mother->pdgId() << endl;
//         EventInfo.BHadron_pT[EventInfo.nBHadrons]    = genIt.p4().pt();
//         EventInfo.BHadron_eta[EventInfo.nBHadrons]   = genIt.p4().eta();
//         EventInfo.BHadron_phi[EventInfo.nBHadrons]   = genIt.p4().phi();
//         EventInfo.BHadron_mass[EventInfo.nBHadrons]  = genIt.mass();
//         EventInfo.BHadron_pdgID[EventInfo.nBHadrons] = genIt.pdgId();
//         EventInfo.BHadron_status[EventInfo.nBHadrons] = genIt.status();
//         EventInfo.BHadron_mother[EventInfo.nBHadrons] = genIt.mother()->pdgId();
//         // check if any of the daughters is also B hadron
//         int hasBHadronDaughter = 0;
//         for (unsigned int d=0; d<nDaughters; ++d) {
//           int daughterID = abs(genIt.daughter(d)->pdgId());
//           if ( (daughterID/100)%10 == 5 || (daughterID/1000)%10 == 5 ) { hasBHadronDaughter = 1; break; }
//         }
//         EventInfo.BHadron_hasBdaughter[EventInfo.nBHadrons] = hasBHadronDaughter;
//         ++EventInfo.nBHadrons;
//       }

      // Final c Hadrons
      if ( (ID/100)%10 == 4 || (ID/1000)%10 == 4 ) {
        // check if daughter is not a DHadron
	bool hasNoDHadronAsDaughter = true;
	for (unsigned int d=0; d<nDaughters; d++) {
	  const Candidate* daughter = genIt.daughter(d);
	  int IDdaughter = abs(daughter->pdgId());
	  bool hasDHadron = (IDdaughter/100)%10 == 4 || (IDdaughter/1000)%10 == 4;
	  if ( hasDHadron ) hasNoDHadronAsDaughter = false;
	}

	// variables
	if ( hasNoDHadronAsDaughter ) {
	  EventInfo.DHadron_pT[EventInfo.nDHadrons]    = genIt.p4().pt();
	  EventInfo.DHadron_eta[EventInfo.nDHadrons]   = genIt.p4().eta();
	  EventInfo.DHadron_phi[EventInfo.nDHadrons]   = genIt.p4().phi();
	  EventInfo.DHadron_mass[EventInfo.nDHadrons]  = genIt.mass();
	  EventInfo.DHadron_pdgID[EventInfo.nDHadrons] = genIt.pdgId();

	  int nCharged = 0;
	  // Loop on daughters
	  int numberChargedDaughters = 0;
	  for (unsigned int d=0; d<nDaughters; d++) {
	    const Candidate* daughter = genIt.daughter(d);
	    int IDdaughter = abs(daughter->pdgId());
	    EventInfo.DHadron_DaughtersPdgID[EventInfo.nDaughters] = IDdaughter;
	    // take vertex of first daughter for SV
	    if ( d == 0 ) {
	      EventInfo.DHadron_SVx[EventInfo.nDHadrons] = daughter->vx();
   	      EventInfo.DHadron_SVy[EventInfo.nDHadrons] = daughter->vy();
    	      EventInfo.DHadron_SVz[EventInfo.nDHadrons] = daughter->vz();
	    }
	    EventInfo.nDaughters++;
	    if ( daughter->charge() != 0 ) numberChargedDaughters++;
	    // charged track multiplicity
	    unsigned int nDaughters2 = daughter->numberOfDaughters();
	    if (nDaughters2 == 0 && daughter->charge() != 0
				 && daughter->p4().pt()>1.) nCharged++;
	    for (unsigned int d2=0; d2<nDaughters2; d2++) {
	      const Candidate* daughter2 = daughter->daughter(d2);
	      unsigned int nDaughters3 = daughter2->numberOfDaughters();
	      if (nDaughters3 == 0 && daughter2->charge() != 0
	        		   && daughter2->p4().pt()>1.) nCharged++;
	      for (unsigned int d3=0; d3<nDaughters3; d3++) {
		const Candidate* daughter3 = daughter2->daughter(d3);
		unsigned int nDaughters4 = daughter3->numberOfDaughters();
		if (nDaughters4 == 0 && daughter3->charge() != 0
				     && daughter3->p4().pt()>1.) nCharged++;
		for (unsigned int d4=0; d4<nDaughters4; d4++) {
		  const Candidate* daughter4 = daughter3->daughter(d4);
		  unsigned int nDaughters5 = daughter4->numberOfDaughters();
		  if (nDaughters5 == 0 && daughter4->charge() != 0
				       && daughter4->p4().pt()>1.) nCharged++;
		} // D -> X -> X -> X -> X
	      } // D -> X -> X -> X
	    } // D -> X -> X
	  } // loop on daughters from final c hadron

	  EventInfo.DHadron_nDaughters[EventInfo.nDHadrons] = nDaughters;
	  EventInfo.DHadron_nChargedDaughters[EventInfo.nDHadrons] = numberChargedDaughters;
	  EventInfo.DHadron_nCharged[EventInfo.nDHadrons] = nCharged;
// cout << " D hadron " << EventInfo.nDHadrons << "  id " << EventInfo.DHadron_pdgID[EventInfo.nDHadrons] << "  pT " << EventInfo.DHadron_pT[EventInfo.nDHadrons] << "  eta " << EventInfo.DHadron_eta[EventInfo.nDHadrons] << "  phi " << EventInfo.DHadron_phi[EventInfo.nDHadrons] << endl;
	  EventInfo.nDHadrons++;
	} // final c hadron
      } // c hadron

      // Leptons
      if ( (ID == 11 || ID == 13 || ID == 15) && genIt.p4().pt() > 3. ) {
       if (genIt.numberOfMothers()>0) {  // protection for herwig ttbar mc
        const Candidate * moth1 = genIt.mother();
        if ( moth1->pdgId() != genIt.pdgId() ) {
          EventInfo.Genlep_pT[EventInfo.nGenlep]    = genIt.p4().pt();
          EventInfo.Genlep_eta[EventInfo.nGenlep]   = genIt.p4().eta();
          EventInfo.Genlep_phi[EventInfo.nGenlep]   = genIt.p4().phi();
          EventInfo.Genlep_pdgID[EventInfo.nGenlep] = genIt.pdgId();
          EventInfo.Genlep_status[EventInfo.nGenlep] = genIt.status();
          int ID1 = abs( moth1->pdgId() );
          int ID2 = -9999;
          int ID3 = -9999;
          int ID4 = -9999;
          int isWZ = 0, istau = 0, isB = 0, isD = 0;
          if (moth1->numberOfMothers()>0) {   // protection for herwig ttbar mc
            const Candidate * moth2 = moth1->mother();
            ID2 = abs( moth2->pdgId() );
            if (moth2->numberOfMothers()>0) {
              const Candidate * moth3 = moth2->mother();
              ID3 = abs( moth3->pdgId() );
              if (moth3->numberOfMothers()>0) {
                const Candidate * moth4 = moth3->mother();
                ID4 = abs( moth4->pdgId() );
              }
            }
          }
          if ( ID1 == 15 ) istau = 1;
          if ( ID1 == 24 || ID2 == 24 || ID3 == 24 || ID4 == 24 ||
              ID1 == 25 || ID2 == 25 || ID3 == 25 || ID4 == 25 ) isWZ = 1;
          if ( ID1 == 411 || ID1 == 421 || ID1 == 431 ||
              ID1 ==4112 || ID1 ==4122 || ID1 ==4132 ||
              ID1 ==4212 || ID1 ==4222 || ID1 ==4232 || ID1 ==4332 ) isD = 1;
          if ( ID2 == 411 || ID2 == 421 || ID2 == 431 ||
              ID2 ==4112 || ID2 ==4122 || ID2 ==4132 ||
              ID2 ==4212 || ID2 ==4222 || ID2 ==4232 || ID2 ==4332 ) isD = 1;
          if ( ID1 == 511 || ID1 == 521 || ID1 == 531 || ID1 == 541 ||
              ID1 ==5112 || ID1 ==5122 || ID1 ==5132 ||
              ID1 ==5212 || ID1 ==5222 || ID1 ==5232 || ID1 ==5332 ) isB = 1;
          if ( ID2 == 511 || ID2 == 521 || ID2 == 531 || ID2 == 541 ||
              ID2 ==5112 || ID2 ==5122 || ID2 ==5132 ||
              ID2 ==5212 || ID2 ==5222 || ID2 ==5232 || ID2 ==5332 ) isB = 1;
          if ( ID3 == 511 || ID3 == 521 || ID3 == 531 || ID3 == 541 ||
              ID3 ==5112 || ID3 ==5122 || ID3 ==5132 ||
              ID3 ==5212 || ID3 ==5222 || ID3 ==5232 || ID3 ==5332 ) isB = 1;
          if ( ID4 == 511 || ID4 == 521 || ID4 == 531 || ID4 == 541 ||
              ID4 ==5112 || ID4 ==5122 || ID4 ==5132 ||
              ID4 ==5212 || ID4 ==5222 || ID4 ==5232 || ID4 ==5332 ) isB = 1;
          if ( isB+isD != 0 ) EventInfo.Genlep_mother[EventInfo.nGenlep] = 5*isB + 4*isD;
          else                EventInfo.Genlep_mother[EventInfo.nGenlep] = 10*istau + 100*isWZ;
          //  cout << " lepton " << EventInfo.nGenlep << " pdgID " << EventInfo.Genlep_pdgID[EventInfo.nGenlep]
          //       << " moth1 " << moth1->pdgId() << " moth2 " << moth2->pdgId()
          //       << " moth3 " << moth3->pdgId() << " moth4 " << moth4->pdgId()
          //       << " pT " << EventInfo.Genlep_pT[EventInfo.nGenlep]
          //       << " from " << EventInfo.Genlep_mother[EventInfo.nGenlep] << endl;
          ++EventInfo.nGenlep;
        }
       }
      }

      // V0: K0s and Lambda
      if ( ID == 310 || ID == 3122 ) {
	 int nCharged = 0;
	 float SVx = 0., SVy = 0., SVz = 0.;
	 for (unsigned int d=0; d<nDaughters; d++) {
	   const Candidate* daughter = genIt.daughter(d);
	   if ( daughter->charge() != 0 && daughter->pt() > 1. ) {
	     nCharged++;
   	     SVx = daughter->vx();
   	     SVy = daughter->vy();
   	     SVz = daughter->vz();
	   }
	 }
	 bool inAcceptance = false;
	 if ( nCharged > 0 ) {
	   float Radius = TMath::Sqrt( SVx*SVx + SVy*SVy );
	   float AbsEta = abs( genIt.p4().eta() );
	   if ( (AbsEta < 2.0 && Radius < 7.3) ||
	        (AbsEta < 2.5 && Radius < 4.4) ||
	        (AbsEta > 1.8 && AbsEta < 2.5 && abs(SVz) < 34.5 ) )
	     inAcceptance = true;
	 }
	 if ( inAcceptance ) {
	   EventInfo.GenV0_pT[EventInfo.nGenV0]    = genIt.p4().pt();
	   EventInfo.GenV0_eta[EventInfo.nGenV0]   = genIt.p4().eta();
	   EventInfo.GenV0_phi[EventInfo.nGenV0]   = genIt.p4().phi();
	   EventInfo.GenV0_pdgID[EventInfo.nGenV0] = genIt.pdgId();
	   EventInfo.GenV0_SVx[EventInfo.nGenV0] = SVx;
	   EventInfo.GenV0_SVy[EventInfo.nGenV0] = SVy;
   	   EventInfo.GenV0_SVz[EventInfo.nGenV0] = SVz;
	   EventInfo.GenV0_nCharged[EventInfo.nGenV0] = nCharged;
	   EventInfo.nGenV0++;
	 }
      }

    } // end loop on generated particles

    //------------------------------------------------------
    // generated particles: looking for b hadrons only
    //------------------------------------------------------
    if ( AreBHadrons ) {
    for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      int ID = abs(genIt.pdgId());
      unsigned int nDaughters = genIt.numberOfDaughters();

      // b Hadrons
      if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
        //  cout << " pdgId " << genIt.pdgId()  << " pT " << genIt.p4().pt() << " mother " << mother->pdgId() << endl;
        EventInfo.BHadron_pT[EventInfo.nBHadrons]    = genIt.p4().pt();
        EventInfo.BHadron_eta[EventInfo.nBHadrons]   = genIt.p4().eta();
        EventInfo.BHadron_phi[EventInfo.nBHadrons]   = genIt.p4().phi();
        EventInfo.BHadron_mass[EventInfo.nBHadrons]  = genIt.mass();
        EventInfo.BHadron_pdgID[EventInfo.nBHadrons] = genIt.pdgId();
        if (genIt.numberOfMothers()>0) EventInfo.BHadron_mother[EventInfo.nBHadrons] = genIt.mother()->pdgId();
        else EventInfo.BHadron_mother[EventInfo.nBHadrons] = -99999;
        // check if any of the daughters is also B hadron
        int hasBHadronDaughter = 0;
        for (unsigned int d=0; d<nDaughters; ++d) {
          int daughterID = abs(genIt.daughter(d)->pdgId());
          if ( (daughterID/100)%10 == 5 || (daughterID/1000)%10 == 5 ) { hasBHadronDaughter = 1; break; }
        }
        EventInfo.BHadron_hasBdaughter[EventInfo.nBHadrons] = hasBHadronDaughter;
	int idD1 = -1, idD2 = -1;
	int nCharged = 0;
	// take vertex of first daughter for SV
	EventInfo.BHadron_SVx[EventInfo.nBHadrons] = (nDaughters>0 ? genIt.daughter(0)->vx() : -9999);
	EventInfo.BHadron_SVy[EventInfo.nBHadrons] = (nDaughters>0 ? genIt.daughter(0)->vy() : -9999);
	EventInfo.BHadron_SVz[EventInfo.nBHadrons] = (nDaughters>0 ? genIt.daughter(0)->vz() : -9999);

	if ( !hasBHadronDaughter ) {
	  // Loop on daughters: B -> X
	  for (unsigned int d=0; d<nDaughters; d++) {
	    const Candidate* daughter = genIt.daughter(d);
	    int IDdaughter = abs(daughter->pdgId());
	    bool isDHadron = (IDdaughter/100)%10 == 4 || (IDdaughter/1000)%10 == 4;
	    bool hasNoDHadronAsDaughter = true;
	    bool isDHadron2, hasNoDHadronAsDaughter2, isDHadron3;
	    unsigned int nDaughters2 = daughter->numberOfDaughters();

            // B -> D
	    if ( isDHadron ) {
	      for (unsigned int d2=0; d2<nDaughters2; d2++) {
	        const Candidate* daughter2 = daughter->daughter(d2);
	        int IDdaughter2 = abs(daughter2->pdgId());
	        isDHadron2 = (IDdaughter2/100)%10 == 4 || (IDdaughter2/1000)%10 == 4;
	        if ( isDHadron2 ) hasNoDHadronAsDaughter = false;
	      }
	      if ( hasNoDHadronAsDaughter ) {
	        for (int c=0; c<EventInfo.nDHadrons; c++) {
	          if ( abs( daughter->p4().pt()
		                   - EventInfo.DHadron_pT[c] ) < 0.0001 ) {
		    if ( idD1 < 0 ) idD1 = c;
		    else if ( idD2 < 0 ) idD2 = c;
// cout << " D hadron from B " << c << "  id " << daughter->pdgId() << "  pT " << daughter->p4().pt()  << "  eta " << daughter->p4().eta() << "  phi " << daughter->p4().phi() << endl;
		  }
	        }
	      }
	    }

	    if ( !isDHadron || !hasNoDHadronAsDaughter ) {
	      if (nDaughters2 == 0 && daughter->charge() != 0
	  	    		   && daughter->p4().pt() > 1.) nCharged++;
	      // B -> X -> X
	      for (unsigned int d2=0; d2<nDaughters2; d2++) {
	        const Candidate* daughter2 = daughter->daughter(d2);
	        int IDdaughter2 = abs(daughter2->pdgId());
	        isDHadron2 = (IDdaughter2/100)%10 == 4 || (IDdaughter2/1000)%10 == 4;
	        hasNoDHadronAsDaughter2 = true;
	        unsigned int nDaughters3 = daughter2->numberOfDaughters();

	        // B -> D*(*) -> D
	        if ( isDHadron2 ) {
	          for (unsigned int d3=0; d3<nDaughters3; d3++) {
	            const Candidate* daughter3 = daughter2->daughter(d3);
	            int IDdaughter3 = abs(daughter3->pdgId());
	            isDHadron3 = (IDdaughter3/100)%10 == 4 || (IDdaughter3/1000)%10 == 4;
	            if ( isDHadron3 ) hasNoDHadronAsDaughter2 = false;
	          }
	          if ( hasNoDHadronAsDaughter2 ) {
	            for (int c=0; c<EventInfo.nDHadrons; c++) {
	              if ( abs( daughter2->p4().pt()
		                       - EventInfo.DHadron_pT[c] ) < 0.0001 ) {
		        if ( idD1 < 0 ) idD1 = c;
		        else if ( idD2 < 0 ) idD2 = c;
// cout << " D hadron from B->D*(*) " << c << "  id " << daughter2->pdgId() << "  pT " << daughter2->p4().pt() << "  eta " << daughter2->p4().eta() << "  phi " << daughter2->p4().phi() << endl;
		      }
	            }
	          }
	        }

	        if ( !isDHadron2 || !hasNoDHadronAsDaughter2 ) {
	          if (nDaughters3 == 0 && daughter2->charge() != 0
	        		       && daughter2->p4().pt() > 1.) nCharged++;
	          // B -> X -> X -> X
	          for (unsigned int d3=0; d3<nDaughters3; d3++) {
	            const Candidate* daughter3 = daughter2->daughter(d3);
	            int IDdaughter3 = abs(daughter3->pdgId());
	            isDHadron3 = (IDdaughter3/100)%10 == 4 || (IDdaughter3/1000)%10 == 4;
	            unsigned int nDaughters4 = daughter3->numberOfDaughters();

	            // B -> D** -> D* -> D
	            if ( isDHadron3 ) {
	              for (int c=0; c<EventInfo.nDHadrons; c++) {
	                if ( abs( daughter3->p4().pt()
		                         - EventInfo.DHadron_pT[c] ) < 0.0001 ) {
		          if ( idD1 < 0 ) idD1 = c;
		          else if ( idD2 < 0 ) idD2 = c;
// cout << " D hadron from B->D**->D* " << c << "  id " << daughter3->pdgId() << "  pT " << daughter3->p4().pt() << "  eta " << daughter3->p4().eta() << "  phi " << daughter3->p4().phi() << endl;
		        }
	              }
	            }
                    else {
	              if (nDaughters4 == 0 && daughter3->charge() != 0
	        	  	           && daughter3->p4().pt() > 1.) nCharged++;
	              // B -> X -> X -> X -> X
	              for (unsigned int d4=0; d4<nDaughters4; d4++) {
	                const Candidate* daughter4 = daughter3->daughter(d4);
	                unsigned int nDaughters5 = daughter4->numberOfDaughters();
	                if (nDaughters5 == 0 && daughter4->charge() != 0
	        	   	             && daughter4->p4().pt()>1.) nCharged++;
	              }
	            }
	          } // B -> X -> X -> X
	        }
	      } // B -> X -> X
	    }
	  } // loop on daughters from final b hadron
// cout << " B hadron " << EventInfo.nBHadrons << "  id " << EventInfo.BHadron_pdgID[EventInfo.nBHadrons] << "  pT " << EventInfo.BHadron_pT[EventInfo.nBHadrons] << "  eta " << EventInfo.BHadron_eta[EventInfo.nBHadrons]<< "  phi " << EventInfo.BHadron_phi[EventInfo.nBHadrons] << "  ->D " << idD1 << " " << idD2 << "  nCh " << nCharged << endl;
	} // final b hadron

    	EventInfo.BHadron_nCharged[EventInfo.nBHadrons] = nCharged;
    	EventInfo.BHadron_DHadron1[EventInfo.nBHadrons] = idD1;
    	EventInfo.BHadron_DHadron2[EventInfo.nBHadrons] = idD2;
 	EventInfo.nBHadrons++;
      } // b hadron
    } // end loop on generated particles
    } // end if ( AreBHadrons )
  }

  //------------------------------------------------------
  // simulated PV
  //------------------------------------------------------
  if ( !isData_ && useTrackHistory_ ) {
    edm::Handle<edm::HepMCProduct> theHepMCProduct;
    iEvent.getByToken(generatorhep,theHepMCProduct);
    std::vector<simPrimaryVertex> simpv;
    simpv=getSimPVs(theHepMCProduct);
    //       cout << "simpv.size() " << simpv.size() << endl;
  }
  //---------------------------- End MC info ---------------------------------------//
  //   std::cout << "EventInfo.Evt:" <<EventInfo.Evt << std::endl;
  //   std::cout << "EventInfo.pthat:" <<EventInfo.pthat << std::endl;


  //------------------------------------------------------
  // ttbar information
  //------------------------------------------------------
  if (use_ttbar_filter_) {

    edm::Handle<int> pIn;
    iEvent.getByToken(ttbartop, pIn);
    EventInfo.ttbar_chan=*pIn;

    edm::Handle<int> triggerIn;
    iEvent.getByToken(ttbartoptrig,triggerIn);
    EventInfo.ttbar_trigWord=*triggerIn;

    int gctr(0);
    edm::Handle<edm::View<reco::GenParticle> > selGen;
    iEvent.getByToken(ttbarproducerGen_,selGen);
    for (size_t i = 0; i < selGen->size(); ++i)
      {
	const auto g = selGen->ptrAt(i);
	EventInfo.ttbar_gpt[gctr] = g->pt();
	EventInfo.ttbar_geta[gctr] = g->eta();
	EventInfo.ttbar_gphi[gctr] = g->phi();	
	EventInfo.ttbar_gm[gctr] = g->mass();
	EventInfo.ttbar_gid[gctr] = g->pdgId();
	gctr++;
      }
    EventInfo.ttbar_ng=gctr;

    edm::Handle<int> metfilterIn;
    iEvent.getByToken(metfilterIntoken,metfilterIn);
    EventInfo.ttbar_metfilterWord=*metfilterIn;

    int lctr(0);
    edm::Handle<edm::View<pat::Electron> > selElectrons;
    iEvent.getByToken(ttbarproducerEle_,selElectrons);
    for (size_t i = 0; i < selElectrons->size(); ++i)
      {
	const auto l = selElectrons->ptrAt(i);
	EventInfo.ttbar_lpt[lctr]  = l->pt();
	EventInfo.ttbar_leta[lctr] = l->eta();
	EventInfo.ttbar_lphi[lctr] = l->phi();
	EventInfo.ttbar_lm[lctr]   = 0;
	EventInfo.ttbar_lch[lctr]  = l->charge();
	EventInfo.ttbar_lid[lctr]  = 11;
	EventInfo.ttbar_lgid[lctr] = l->genParticle() ? l->genParticle()->pdgId() : 0;
	lctr++;
      }

    edm::Handle<edm::View<pat::Muon> > selMuons;
    iEvent.getByToken(ttbarproducerMuon_,selMuons);
    for (size_t i = 0; i < selMuons->size(); ++i)
      {
	const auto l = selMuons->ptrAt(i);
	EventInfo.ttbar_lpt[lctr]  = l->pt();
	EventInfo.ttbar_leta[lctr] = l->eta();
	EventInfo.ttbar_lphi[lctr] = l->phi();
	EventInfo.ttbar_lm[lctr]   = 0;
	EventInfo.ttbar_lch[lctr]  = l->charge();
	EventInfo.ttbar_lid[lctr]  = 13;
	EventInfo.ttbar_lgid[lctr] = l->genParticle() ? l->genParticle()->pdgId() : 0;
	lctr++;
      }
    EventInfo.ttbar_nl=lctr;

    edm::Handle<edm::View<pat::MET> > selMETs;
    iEvent.getByToken(ttbarproducerMET_,selMETs);
    EventInfo.ttbar_metpt=selMETs->ptrAt(0)->pt();
    EventInfo.ttbar_metphi=selMETs->ptrAt(0)->phi();

    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoTag_,rhoH);
    EventInfo.ttbar_rho = *rhoH;

    //generator information
    EventInfo.ttbar_nw=0;
    if(!isData_)
      {
	edm::Handle<GenEventInfoProduct> evt;
	iEvent.getByToken(generatorevt, evt);
	if(evt.isValid())
	  {
	    EventInfo.ttbar_allmepartons   = evt->nMEPartons();
	    EventInfo.ttbar_matchmepartons = evt->nMEPartonsFiltered();
	    EventInfo.ttbar_w[0]           = evt->weight();
	    EventInfo.ttbar_nw++;
	  }
	edm::Handle<LHEEventProduct> evet;
	iEvent.getByToken(generatorlhe, evet);
	if(evet.isValid())
	  {
	    double asdd=evet->originalXWGTUP();
	    for(unsigned int i=0; i<evet->weights().size();i++){
	      double asdde=evet->weights()[i].wgt;
	      EventInfo.ttbar_w[EventInfo.ttbar_nw]=EventInfo.ttbar_w[0]*asdde/asdd;
 	      EventInfo.ttbar_nw++;
	    }
	  }
      }
  }

  //------------------------------------------------------
  // PAT Muons
  //------------------------------------------------------
  edm::Handle<std::vector<pat::Muon> >  patMuonsHandle;
  if( storePatMuons_ )
  {
    iEvent.getByToken(patMuonCollectionName_,patMuonsHandle);

    for( std::vector<pat::Muon>::const_iterator it = patMuonsHandle->begin(); it != patMuonsHandle->end(); ++it )
    {
      if( !it->isGlobalMuon() ) continue;

      EventInfo.PatMuon_isGlobal[EventInfo.nPatMuon] = 1;
      EventInfo.PatMuon_isPF[EventInfo.nPatMuon]     = it->isPFMuon();
      EventInfo.PatMuon_nTkHit[EventInfo.nPatMuon]   = it->innerTrack()->hitPattern().numberOfValidHits();
      EventInfo.PatMuon_nPixHit[EventInfo.nPatMuon]  = it->innerTrack()->hitPattern().numberOfValidPixelHits();
      EventInfo.PatMuon_nOutHit[EventInfo.nPatMuon]  = it->innerTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS);
      EventInfo.PatMuon_nMuHit[EventInfo.nPatMuon]   = it->outerTrack()->hitPattern().numberOfValidMuonHits();
      EventInfo.PatMuon_nMatched[EventInfo.nPatMuon] = it->numberOfMatches();
      EventInfo.PatMuon_chi2[EventInfo.nPatMuon]     = it->globalTrack()->normalizedChi2();
      EventInfo.PatMuon_chi2Tk[EventInfo.nPatMuon]   = it->innerTrack()->normalizedChi2();
      EventInfo.PatMuon_pt[EventInfo.nPatMuon]       = it->pt();
      EventInfo.PatMuon_eta[EventInfo.nPatMuon]      = it->eta();
      EventInfo.PatMuon_phi[EventInfo.nPatMuon]      = it->phi();
      EventInfo.PatMuon_vz[EventInfo.nPatMuon]       = it->vz();
      EventInfo.PatMuon_IP[EventInfo.nPatMuon]       = it->dB(pat::Muon::PV3D);
      EventInfo.PatMuon_IPsig[EventInfo.nPatMuon]    = (it->dB(pat::Muon::PV3D))/(it->edB(pat::Muon::PV3D));
      EventInfo.PatMuon_IP2D[EventInfo.nPatMuon]     = it->dB(pat::Muon::PV2D);
      EventInfo.PatMuon_IP2Dsig[EventInfo.nPatMuon]  = (it->dB(pat::Muon::PV2D))/(it->edB(pat::Muon::PV2D));

      ++EventInfo.nPatMuon;
    }
  }


  //------------------------------------------------------
  // Muons
  //------------------------------------------------------
  edm::Handle<edm::View<reco::Muon> >  muonsHandle;
  iEvent.getByToken(muonCollectionName_,muonsHandle);
  muons = *muonsHandle;

  //----------------------------------------
  // Transient track for IP calculation
  //----------------------------------------
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  //------------------
  // Primary vertex
  //------------------
  iEvent.getByToken(primaryVertexColl_,primaryVertex);
  //bool newvertex = false;

  bool pvFound = (primaryVertex->size() != 0);
  if ( pvFound ) {
    pv = &(*primaryVertex->begin());
  }
  else {
    reco::Vertex::Error e;
    e(0,0)=0.0015*0.0015;
    e(1,1)=0.0015*0.0015;
    e(2,2)=15.*15.;
    reco::Vertex::Point p(0,0,0);
    pv=  new reco::Vertex(p,e,1,1,1);
    //newvertex = true;
  }
  //   GlobalPoint Pv_point = GlobalPoint((*pv).x(), (*pv).y(), (*pv).z());
  EventInfo.PVz = (*primaryVertex)[0].z();
  EventInfo.PVez = (*primaryVertex)[0].zError();

  EventInfo.nPV=0;
  for (unsigned int i = 0; i< primaryVertex->size() ; ++i) {
    EventInfo.PV_x[EventInfo.nPV]      = (*primaryVertex)[i].x();
    EventInfo.PV_y[EventInfo.nPV]      = (*primaryVertex)[i].y();
    EventInfo.PV_z[EventInfo.nPV]      = (*primaryVertex)[i].z();
    EventInfo.PV_ex[EventInfo.nPV]     = (*primaryVertex)[i].xError();
    EventInfo.PV_ey[EventInfo.nPV]     = (*primaryVertex)[i].yError();
    EventInfo.PV_ez[EventInfo.nPV]     = (*primaryVertex)[i].zError();
    EventInfo.PV_chi2[EventInfo.nPV]   = (*primaryVertex)[i].normalizedChi2();
    EventInfo.PV_ndf[EventInfo.nPV]    = (*primaryVertex)[i].ndof();
    EventInfo.PV_isgood[EventInfo.nPV] = (*primaryVertex)[i].isValid();
    EventInfo.PV_isfake[EventInfo.nPV] = (*primaryVertex)[i].isFake();

    ++EventInfo.nPV;
  }

  //------------------------------------------------------
  // Trigger info
  //------------------------------------------------------

  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(triggerTable_, trigRes);

  EventInfo.nBitTrigger = int(triggerPathNames_.size()/32)+1;
  for(int i=0; i<EventInfo.nBitTrigger; ++i) EventInfo.BitTrigger[i] = 0;

  std::vector<std::string> triggerList;
  edm::Service<edm::service::TriggerNamesService> tns;
  bool foundNames = tns->getTrigPaths(*trigRes,triggerList);
  if ( !foundNames ) edm::LogError("TriggerNamesNotFound") << "Could not get trigger names!";
  if ( trigRes->size() != triggerList.size() ) edm::LogError("TriggerPathLengthMismatch") << "Length of names and paths not the same: "
    << triggerList.size() << "," << trigRes->size() ;

  processTrig(trigRes, triggerList);

  PFJet80 = false;
  for(std::vector<std::string>::const_iterator itTrigPathNames = PFJet80TriggerPathNames_.begin();
      itTrigPathNames != PFJet80TriggerPathNames_.end(); ++itTrigPathNames)
  {
    std::vector<std::string>::const_iterator it = std::find(triggerPathNames_.begin(), triggerPathNames_.end(), *itTrigPathNames);
    if( it != triggerPathNames_.end() )
    {
      int triggerIdx = ( it - triggerPathNames_.begin() );
      int bitIdx = int(triggerIdx/32);
      if ( EventInfo.BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) PFJet80 = true;
    }
  }

  //------------- added by Camille-----------------------------------------------------------//
  edm::ESHandle<JetTagComputer> computerHandle;
  iSetup.get<JetTagComputerRecord>().get( SVComputer_.c_str(), computerHandle );

  computer = dynamic_cast<const GenericMVAJetTagComputer*>( computerHandle.product() );
  //------------- end added-----------------------------------------------------------//
  if (doCTag_){
     edm::ESHandle<JetTagComputer> slcomputerHandle;
     iSetup.get<JetTagComputerRecord>().get( SLComputer_.c_str(), slcomputerHandle );

     slcomputer = dynamic_cast<const GenericMVAJetTagComputer*>( slcomputerHandle.product() );
  }
  //------------------------------------------------------
  // All tracks info
  //------------------------------------------------------
  edm::Handle<reco::TrackCollection> tracksHandle;
  if( produceAllTrackTree_ && storeEventInfo_ )
  {
    iEvent.getByToken(tracksColl_,tracksHandle);
    for(reco::TrackCollection::const_iterator trk = tracksHandle->begin(); trk != tracksHandle->end(); ++trk)
    {
      EventInfo.TrkAll_d0[EventInfo.nTrkAll]        = -trk->dxy((*primaryVertex)[0].position());
      EventInfo.TrkAll_dz[EventInfo.nTrkAll]        = trk->dz((*primaryVertex)[0].position());
      EventInfo.TrkAll_p[EventInfo.nTrkAll]         = trk->p();
      EventInfo.TrkAll_pt[EventInfo.nTrkAll]        = trk->pt();
      EventInfo.TrkAll_eta[EventInfo.nTrkAll]       = trk->eta();
      EventInfo.TrkAll_phi[EventInfo.nTrkAll]       = trk->phi();
      EventInfo.TrkAll_chi2[EventInfo.nTrkAll]      = trk->normalizedChi2();
      EventInfo.TrkAll_charge[EventInfo.nTrkAll]    = trk->charge();
      EventInfo.TrkAll_nHitAll[EventInfo.nTrkAll]   = trk->numberOfValidHits();
      EventInfo.TrkAll_nHitPixel[EventInfo.nTrkAll] = trk->hitPattern().numberOfValidPixelHits();
      EventInfo.TrkAll_nHitStrip[EventInfo.nTrkAll] = trk->hitPattern().numberOfValidStripHits();
      EventInfo.TrkAll_nHitTIB[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTIBHits();
      EventInfo.TrkAll_nHitTID[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTIDHits();
      EventInfo.TrkAll_nHitTOB[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTOBHits();
      EventInfo.TrkAll_nHitTEC[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTECHits();
      EventInfo.TrkAll_nHitPXB[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidPixelBarrelHits();
      EventInfo.TrkAll_nHitPXF[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidPixelEndcapHits();
      EventInfo.TrkAll_isHitL1[EventInfo.nTrkAll]   = trk->hitPattern().hasValidHitInFirstPixelBarrel();
      EventInfo.TrkAll_nSiLayers[EventInfo.nTrkAll] = trk->hitPattern().trackerLayersWithMeasurement();
      EventInfo.TrkAll_nPxLayers[EventInfo.nTrkAll] = trk->hitPattern().pixelLayersWithMeasurement();

      reco::TrackRef trkRef( tracksHandle, trk - tracksHandle->begin() );

      setTracksPVBase(trkRef, primaryVertex,
                      EventInfo.TrkAll_PV[EventInfo.nTrkAll],
                      EventInfo.TrkAll_PVweight[EventInfo.nTrkAll]);

      ++EventInfo.nTrkAll;
    }
  }
  //------------------------------------------------------

  //------------------------------------------------------
  // Jet info
  //------------------------------------------------------
  int iJetColl = 0 ;
  //// Do jets
  processJets(jetsColl, jetsColl, subjetColls, iEvent, iSetup, iJetColl); // the second 'jetsColl' is a dummy input here
  if (runSubJets_) {
    // for fat jets we might have a different jet tag computer
    // if so, grab the fat jet tag computer
    if(SVComputerSubJets_!=SVComputer_)
    {
      iSetup.get<JetTagComputerRecord>().get( SVComputerSubJets_.c_str(), computerHandle );

      computer = dynamic_cast<const GenericMVAJetTagComputer*>( computerHandle.product() );
    }
    for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
    {
      iJetColl = 1+i ;
      processJets(subjetColls[i], jetsColl, subjetColls, iEvent, iSetup, iJetColl); // 'subjetColls' is a dummy input here
    }
  }
  //------------------------------------------------------

  //// Fill TTree
  if ( EventInfo.BitTrigger > 0 || EventInfo.Run < 0 ) {
    smalltree->Fill();
  }

  return;
}

template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::processTrig(const edm::Handle<edm::TriggerResults>& trigRes, const std::vector<std::string>& triggerList)
{
  for (unsigned int i = 0; i < trigRes->size(); ++i) {

    if ( !trigRes->at(i).accept() ) continue;


    for (std::vector<std::string>::const_iterator itTrigPathNames = triggerPathNames_.begin();
        itTrigPathNames != triggerPathNames_.end(); ++itTrigPathNames)
    {
      int triggerIdx = ( itTrigPathNames - triggerPathNames_.begin() );
      int bitIdx = int(triggerIdx/32);
      if ( NameCompatible(*itTrigPathNames,triggerList[i]) ) EventInfo.BitTrigger[bitIdx] |= ( 1 << (triggerIdx - bitIdx*32) );
    }

  } //// Loop over trigger names

  return;
}

// This is needed to get a TrackingParticle --> Cluster match (instead of Cluster-->TP) (only needed in processJets)
using P = std::pair<OmniClusterRef, TrackingParticleRef>;
bool compare(const P& i, const P& j) {
  	return i.second.index() > j.second.index();
}

template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::processJets(const edm::Handle<PatJetCollection>& jetsColl, const edm::Handle<PatJetCollection>& jetsColl2,
					  const std::vector<edm::Handle<PatJetCollection> >& subjetColls,
					  const edm::Event& iEvent, const edm::EventSetup& iSetup, const int iJetColl)
{

  // matching hit-clusters to TrackingParticles
  edm::Handle<ClusterTPAssociation> clusterToTPMap;
  if( produceJetTrackTruthTree_ )
    iEvent.getByToken(clusterTPMapToken_, clusterToTPMap);

  int numjet = 0;
  //JetInfo[iJetColl].nMuon = 0;
  JetInfo[iJetColl].nPFElectron = 0;
  JetInfo[iJetColl].nPFMuon = 0;

  JetInfo[iJetColl].nJet = 0;
  JetInfo[iJetColl].nTrack = 0;
  JetInfo[iJetColl].nTrackTruth = 0;
  JetInfo[iJetColl].nTrkInc = 0;
  JetInfo[iJetColl].nSV = 0;
  JetInfo[iJetColl].nTrkTagVar = 0;
  JetInfo[iJetColl].nSVTagVar = 0;
  JetInfo[iJetColl].nTrkTagVarCSV = 0;
  JetInfo[iJetColl].nTrkEtaRelTagVarCSV = 0;
  JetInfo[iJetColl].nTrkCTagVar = 0;
  JetInfo[iJetColl].nTrkEtaRelCTagVar = 0;
  JetInfo[iJetColl].nLeptons = 0;  

  //Initialize new test variables for AK4 jets: to be cleaned up in the future
  JetInfo[iJetColl].Jet_trackSip2dSig_AboveBottom_0[JetInfo[iJetColl].nJet] = -19.;
  JetInfo[iJetColl].Jet_trackSip2dSig_AboveBottom_1[JetInfo[iJetColl].nJet] = -19.;

  if ( runFatJets_ && runSubJets_ && iJetColl == 0 )
  {
    for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
      SubJetInfo[SubJetLabels_[i]].nSubJet = 0;
  }

  bool storeTagVariables    = storeTagVariables_;
  bool storeCSVTagVariables = storeCSVTagVariables_;
  bool storeCTagVariables   = storeCTagVariables_;
  if ( runSubJets_ && iJetColl > 0 )
  {
    storeTagVariables    = storeTagVariablesSubJets_;
    storeCSVTagVariables = storeCSVTagVariablesSubJets_;
  }

  double _distJetAxis = distJetAxis_;
  double _decayLength = decayLength_;
  double _deltaR      = deltaR_;
  if ( runSubJets_ && iJetColl > 0 )
  {
    _distJetAxis = distJetAxisSubJets_;
    _decayLength = decayLengthSubJets_;
    _deltaR      = deltaRSubJets_;
  }
  //// Loop over the jets
  for ( PatJetCollection::const_iterator pjet = jetsColl->begin(); pjet != jetsColl->end(); ++pjet ) {

    double ptjet  = pjet->pt()  ;
    double etajet = pjet->eta() ;
    double phijet = pjet->phi() ;

    if ( allowJetSkipping_ && ( ptjet < minJetPt_ || std::fabs( etajet ) > maxJetEta_ ) ) continue;

    //// overlap removal with lepton from ttbar selection
    if (use_ttbar_filter_) {
      float minDRlj(9999.);
      TLorentzVector thejet;
      thejet.SetPtEtaPhiM(ptjet, etajet, phijet, 0.);      
      for(int il=0; il<EventInfo.ttbar_nl; il++)
	{
	  TLorentzVector theLepton;
	  theLepton.SetPtEtaPhiM(EventInfo.ttbar_lpt[il],
				 EventInfo.ttbar_leta[il],
				 EventInfo.ttbar_lphi[il],
				 EventInfo.ttbar_lm[il]);
	  float dR(thejet.DeltaR(theLepton));
	  if(dR>minDRlj) continue;
	  minDRlj=dR;
	}
      if(EventInfo.ttbar_chan>=0 && minDRlj<0.4) continue;
    }
    //// end of removal

    int flavour  =-1  ;
    if ( !isData_ ) {
      flavour = abs( pjet->partonFlavour() );
      if ( flavour >= 1 && flavour <= 3 ) flavour = 1;
    }

		int hflav = pjet->hadronFlavour();		
		int pflav = pjet->partonFlavour();
		int cflav = 0; //~correct flavour definition
		if(!isData_) {
			if(hflav != 0) {
				cflav = hflav;
			}
			else { //not a heavy jet
				if(std::abs(pflav) == 4 || std::abs(pflav) == 5) {
					cflav = 0;
				}
				else {
					cflav = pflav;
				}
			}
		}

    JetInfo[iJetColl].Jet_partonid[JetInfo[iJetColl].nJet]  = pjet->genParton() ? pjet->genParton()->pdgId() : 0;
    JetInfo[iJetColl].Jet_area[JetInfo[iJetColl].nJet]      = pjet->jetArea();
    JetInfo[iJetColl].Jet_flavour[JetInfo[iJetColl].nJet]   = cflav;
		JetInfo[iJetColl].Jet_flavourCleaned[JetInfo[iJetColl].nJet]  = cflav;
    JetInfo[iJetColl].Jet_partonFlavour[JetInfo[iJetColl].nJet]   = pjet->partonFlavour();
    JetInfo[iJetColl].Jet_hadronFlavour[JetInfo[iJetColl].nJet]   = pjet->hadronFlavour();
    JetInfo[iJetColl].Jet_nbHadrons[JetInfo[iJetColl].nJet] = pjet->jetFlavourInfo().getbHadrons().size();
    JetInfo[iJetColl].Jet_ncHadrons[JetInfo[iJetColl].nJet] = pjet->jetFlavourInfo().getcHadrons().size();
    JetInfo[iJetColl].Jet_eta[JetInfo[iJetColl].nJet]       = pjet->eta();
    JetInfo[iJetColl].Jet_phi[JetInfo[iJetColl].nJet]       = pjet->phi();
    JetInfo[iJetColl].Jet_pt[JetInfo[iJetColl].nJet]        = pjet->pt();
    JetInfo[iJetColl].Jet_mass[JetInfo[iJetColl].nJet]      = pjet->mass();
    JetInfo[iJetColl].Jet_genpt[JetInfo[iJetColl].nJet]     = ( pjet->genJet()!=0 ? pjet->genJet()->pt() : -1. );

    // generator-level jet cleaning
    // against prompt leptons
    double dRlj;
    double dRele = 1000., dRmuo = 1000., dRtau = 1000.;
    for (int k = 0; k < EventInfo.nGenlep; k++) {
      if ( EventInfo.Genlep_mother[k]==0 ) continue; //protection for QCD samples
      if ( EventInfo.Genlep_mother[k]%10 !=0 ) continue;
      int lepid = TMath::Abs(EventInfo.Genlep_pdgID[k]);
      dRlj = reco::deltaR( pjet->eta(), pjet->phi(), EventInfo.Genlep_eta[k], EventInfo.Genlep_phi[k] );
      if ( lepid == 11 && dRlj < dRele ) dRele = dRlj;
      if ( lepid == 13 && dRlj < dRmuo ) dRmuo = dRlj;
      if ( lepid == 15 && dRlj < dRtau ) dRtau = dRlj;
    }
    if ( dRele < 0.2 || dRmuo < 0.2 || dRtau < 0.3 ) JetInfo[iJetColl].Jet_flavourCleaned[JetInfo[iJetColl].nJet] = -999;

    // against pileup:
    if ( JetInfo[iJetColl].Jet_genpt[JetInfo[iJetColl].nJet] < 8. )  JetInfo[iJetColl].Jet_flavourCleaned[JetInfo[iJetColl].nJet] = -100;

    // available JEC sets
    unsigned int nJECSets = pjet->availableJECSets().size();
    // PF jet ID
    pat::strbitset retpf = pfjetIDLoose_.getBitTemplate();
    retpf.set(false);
    JetInfo[iJetColl].Jet_looseID[JetInfo[iJetColl].nJet] = ( ( nJECSets>0 && pjet->isPFJet() ) ? ( pfjetIDLoose_( *pjet, retpf ) ? 1 : 0 ) : 0 );
    retpf.set(false);
    JetInfo[iJetColl].Jet_tightID[JetInfo[iJetColl].nJet] = ( ( nJECSets>0 && pjet->isPFJet() ) ? ( pfjetIDTight_( *pjet, retpf ) ? 1 : 0 ) : 0 );

    JetInfo[iJetColl].Jet_jes[JetInfo[iJetColl].nJet]      = ( nJECSets>0 ? pjet->pt()/pjet->correctedJet("Uncorrected").pt() : 1. );
    JetInfo[iJetColl].Jet_residual[JetInfo[iJetColl].nJet] = ( nJECSets>0 ? pjet->pt()/pjet->correctedJet("L3Absolute").pt() : 1. );
		JetInfo[iJetColl].Jet_uncorrpt[JetInfo[iJetColl].nJet] = ( nJECSets>0 ? pjet->correctedJet("Uncorrected").pt() : pjet->pt());

    if( runSubJets_ && iJetColl > 0 )
    {
      int fatjetIdx=-1;
      bool fatJetFound = false;
      // loop over fat jets
      for( int fj = 0; fj < (int)jetsColl2->size(); ++fj )
      {
        const PatJetPtrCollection & subjets = jetsColl2->at(fj).subjets(SubJetLabels_[iJetColl-1]);
        // loop over subjets
        for( int sj = 0; sj < (int)subjets.size(); ++sj )
        {
          if( pjet->originalObjectRef() == subjets.at(sj)->originalObjectRef() )
          {
            fatjetIdx = fj;
            fatJetFound = true;
            break;
          }
        }
        if( fatJetFound )
          break;
      }
      JetInfo[iJetColl].Jet_FatJetIdx[JetInfo[iJetColl].nJet] = fatjetIdx;

      if( ptjet==0. ) // special treatment for pT=0 subjets
      {
        ++numjet;
        ++JetInfo[iJetColl].nJet;
        continue;
      }
    }

    if ( runFatJets_ && iJetColl == 0 )
    {
      // N-subjettiness
      JetInfo[iJetColl].Jet_tau1[JetInfo[iJetColl].nJet] = ( pjet->hasUserFloat("Njettiness:tau1") ? pjet->userFloat("Njettiness:tau1") : pjet->userFloat("NjettinessAK8:tau1") );
      JetInfo[iJetColl].Jet_tau2[JetInfo[iJetColl].nJet] = ( pjet->hasUserFloat("Njettiness:tau2") ? pjet->userFloat("Njettiness:tau2") : pjet->userFloat("NjettinessAK8:tau2") );
      // SoftDrop kinematics
      JetInfo[iJetColl].Jet_ptSoftDrop[JetInfo[iJetColl].nJet]    = ( pjet->hasUserFloat("SoftDrop:Pt")         ? pjet->userFloat("SoftDrop:Pt")         : 0. );
      JetInfo[iJetColl].Jet_etaSoftDrop[JetInfo[iJetColl].nJet]   = ( pjet->hasUserFloat("SoftDrop:Eta")        ? pjet->userFloat("SoftDrop:Eta")        : 0. );
      JetInfo[iJetColl].Jet_phiSoftDrop[JetInfo[iJetColl].nJet]   = ( pjet->hasUserFloat("SoftDrop:Phi")        ? pjet->userFloat("SoftDrop:Phi")        : 0. );
      JetInfo[iJetColl].Jet_massSoftDrop[JetInfo[iJetColl].nJet]  = ( pjet->hasUserFloat("SoftDrop:Mass")       ? pjet->userFloat("SoftDrop:Mass")       : pjet->userFloat("ak8PFJetsCHSSoftDropMass") );
      JetInfo[iJetColl].Jet_jecF0SoftDrop[JetInfo[iJetColl].nJet] = ( pjet->hasUserFloat("SoftDrop:jecFactor0") ? pjet->userFloat("SoftDrop:jecFactor0") : 0. );
      // Pruned kinematics
      JetInfo[iJetColl].Jet_ptPruned[JetInfo[iJetColl].nJet]    = ( pjet->hasUserFloat("Pruned:Pt")         ? pjet->userFloat("Pruned:Pt")         : 0. );
      JetInfo[iJetColl].Jet_etaPruned[JetInfo[iJetColl].nJet]   = ( pjet->hasUserFloat("Pruned:Eta")        ? pjet->userFloat("Pruned:Eta")        : 0. );
      JetInfo[iJetColl].Jet_phiPruned[JetInfo[iJetColl].nJet]   = ( pjet->hasUserFloat("Pruned:Phi")        ? pjet->userFloat("Pruned:Phi")        : 0. );
      JetInfo[iJetColl].Jet_massPruned[JetInfo[iJetColl].nJet]  = ( pjet->hasUserFloat("Pruned:Mass")       ? pjet->userFloat("Pruned:Mass")       : pjet->userFloat("ak8PFJetsCHSPrunedMass") );
      JetInfo[iJetColl].Jet_jecF0Pruned[JetInfo[iJetColl].nJet] = ( pjet->hasUserFloat("Pruned:jecFactor0") ? pjet->userFloat("Pruned:jecFactor0") : 0. );
    }
    if ( runFatJets_ && runSubJets_ && iJetColl == 0 )
    {
      for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
      {
        SubJetInfo[SubJetLabels_[i]].Jet_nFirstSJ[JetInfo[iJetColl].nJet] = SubJetInfo[SubJetLabels_[i]].nSubJet;

        PatJetPtrCollection subjets = pjet->subjets(SubJetLabels_[i]);

        // loop over subjets
        for( int sj = 0; sj < (int)subjets.size(); ++sj )
        {
          SubJetInfo[SubJetLabels_[i]].SubJetIdx[SubJetInfo[SubJetLabels_[i]].nSubJet] = -1;
          // loop over subjet collection
          for( PatJetCollection::const_iterator sjIt = subjetColls[i]->begin(); sjIt != subjetColls[i]->end(); ++sjIt )
          {
            if( sjIt->originalObjectRef() == subjets.at(sj)->originalObjectRef() )
            {
              SubJetInfo[SubJetLabels_[i]].SubJetIdx[SubJetInfo[SubJetLabels_[i]].nSubJet] = ( sjIt - subjetColls[i]->begin() );
              break;
            }
          }
          ++SubJetInfo[SubJetLabels_[i]].nSubJet;
        }

        SubJetInfo[SubJetLabels_[i]].Jet_nSubJets[JetInfo[iJetColl].nJet] = subjets.size();
        SubJetInfo[SubJetLabels_[i]].Jet_nLastSJ[JetInfo[iJetColl].nJet] = SubJetInfo[SubJetLabels_[i]].nSubJet;

        // sort subjets by uncorrected Pt
        std::sort(subjets.begin(), subjets.end(), orderByPt("Uncorrected"));

        int nsubjettracks = 0, nsharedsubjettracks = 0;

        // take two leading subjets
        if( subjets.size()>1 ) // protection for pathological cases of groomed fat jets with only one constituent which results in an undefined subjet 2 index
        {
          for(int sj=0; sj<2; ++sj)
          {
            int subjetIdx = (sj==0 ? 0 : 1); // subjet index
            int compSubjetIdx = (sj==0 ? 1 : 0); // companion subjet index
            int nTracks = ( subjets.at(subjetIdx)->hasTagInfo(ipTagInfos_.c_str()) ? toIPTagInfo(*(subjets.at(subjetIdx)),ipTagInfos_)->selectedTracks().size() : 0 );

            for(int t=0; t<nTracks; ++t)
            {
              if( reco::deltaR( toIPTagInfo(*(subjets.at(subjetIdx)),ipTagInfos_)->selectedTracks().at(t)->eta(), toIPTagInfo(*(subjets.at(subjetIdx)),ipTagInfos_)->selectedTracks().at(t)->phi(), subjets.at(subjetIdx)->eta(), subjets.at(subjetIdx)->phi() ) < deltaRSubJets_ )
              {
                ++nsubjettracks;
                if( reco::deltaR( toIPTagInfo(*(subjets.at(subjetIdx)),ipTagInfos_)->selectedTracks().at(t)->eta(), toIPTagInfo(*(subjets.at(subjetIdx)),ipTagInfos_)->selectedTracks().at(t)->phi(), subjets.at(compSubjetIdx)->eta(), subjets.at(compSubjetIdx)->phi() ) < deltaRSubJets_ )
                {
                  if(sj==0) ++nsharedsubjettracks;
                }
              }
            }
          }
        }

        SubJetInfo[SubJetLabels_[i]].Jet_nsubjettracks[JetInfo[iJetColl].nJet] = nsubjettracks-nsharedsubjettracks;
        SubJetInfo[SubJetLabels_[i]].Jet_nsharedsubjettracks[JetInfo[iJetColl].nJet] = nsharedsubjettracks;
      }
    }

    // Get all TagInfo pointers
    const IPTagInfo *ipTagInfo = toIPTagInfo(*pjet,ipTagInfos_);
    const SVTagInfo *svTagInfo = toSVTagInfo(*pjet,svTagInfos_);
    const SVTagInfo *svNegTagInfo = toSVTagInfo(*pjet,svNegTagInfos_);
    const reco::CandSoftLeptonTagInfo *softPFMuTagInfo = pjet->tagInfoCandSoftLepton(softPFMuonTagInfos_.c_str());
    const reco::CandSoftLeptonTagInfo *softPFElTagInfo = pjet->tagInfoCandSoftLepton(softPFElectronTagInfos_.c_str());
    const reco::BoostedDoubleSVTagInfo *bdsvTagInfo = nullptr;
    if ( runFatJets_ && iJetColl == 0 )
    {
      bdsvTagInfo = pjet->tagInfoBoostedDoubleSV(bdsvTagInfos_.c_str());
    }

    //Get all CTagInfo pointers
    const IPTagInfo *ipTagInfoCTag = toIPTagInfo(*pjet,ipTagInfosCTag_);
    const SVTagInfo *svTagInfoCTag = toSVTagInfo(*pjet,svTagInfosCTag_);
    const reco::CandSoftLeptonTagInfo *softPFMuTagInfoCTag = pjet->tagInfoCandSoftLepton(softPFMuonTagInfosCTag_.c_str());
    const reco::CandSoftLeptonTagInfo *softPFElTagInfoCTag = pjet->tagInfoCandSoftLepton(softPFElectronTagInfosCTag_.c_str()); 

    // Re-calculate N-subjettiness
    std::vector<fastjet::PseudoJet> currentAxes;
    if ( runFatJets_ && iJetColl == 0 )
    {
      float tau1 = JetInfo[iJetColl].Jet_tau1[JetInfo[iJetColl].nJet];
      float tau2 = JetInfo[iJetColl].Jet_tau2[JetInfo[iJetColl].nJet];

      // re-calculate N-subjettiness
      recalcNsubjettiness(*pjet,tau1,tau2,currentAxes);

      // store N-subjettiness axes
      if(currentAxes.size()>0)
      {
        JetInfo[iJetColl].Jet_tauAxis1_px[JetInfo[iJetColl].nJet] = currentAxes[0].px();
        JetInfo[iJetColl].Jet_tauAxis1_py[JetInfo[iJetColl].nJet] = currentAxes[0].py();
        JetInfo[iJetColl].Jet_tauAxis1_pz[JetInfo[iJetColl].nJet] = currentAxes[0].pz();
      }
      if(currentAxes.size()>1)
      {
        JetInfo[iJetColl].Jet_tauAxis2_px[JetInfo[iJetColl].nJet] = currentAxes[1].px();
        JetInfo[iJetColl].Jet_tauAxis2_py[JetInfo[iJetColl].nJet] = currentAxes[1].py();
        JetInfo[iJetColl].Jet_tauAxis2_pz[JetInfo[iJetColl].nJet] = currentAxes[1].pz();
      }
    }

    //*****************************************************************
    // Taggers
    //*****************************************************************

    // Loop on Selected Tracks
    int nseltracks = 0;
    std::map<std::string,int> nsharedtracks;
    std::map<std::string,PatJetPtrCollection> subjets;
    if ( runFatJets_ && runSubJets_ && iJetColl == 0 )
    {
      for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
      {
        nsharedtracks[SubJetLabels_[i]] = 0;
        subjets[SubJetLabels_[i]] = pjet->subjets(SubJetLabels_[i]);
        // sort subjets by uncorrected Pt
        std::sort(subjets[SubJetLabels_[i]].begin(), subjets[SubJetLabels_[i]].end(), orderByPt("Uncorrected"));
      }
    }

    reco::TrackKinematics allKinematics;
    
    JetInfo[iJetColl].Jet_trackSip2dSigAboveCharm_0[JetInfo[iJetColl].nJet]  = -19.;
    JetInfo[iJetColl].Jet_trackSip2dSigAboveCharm_1[JetInfo[iJetColl].nJet]  = -19.;
    JetInfo[iJetColl].Jet_trackSip2dSigAboveBottom_0[JetInfo[iJetColl].nJet] = -19.;
    JetInfo[iJetColl].Jet_trackSip2dSigAboveBottom_1[JetInfo[iJetColl].nJet] = -19.;

    edm::RefToBase<pat::Jet> patJetRef = jetsColl->refAt(pjet - jetsColl->begin());
    edm::RefToBase<reco::Jet> jetRef(patJetRef);
    reco::JetTagInfo jetTagInfo(jetRef);

    const Tracks & selectedTracks( ipTagInfo->selectedTracks() );
    const Tracks & tracks = ( useSelectedTracks_ ? Tracks() : toAllTracks(*pjet,ipTagInfos_,jetTagInfo,iJetColl) );

    int ntagtracks = 0;
    if (useSelectedTracks_) ntagtracks = selectedTracks.size();
    else ntagtracks = tracks.size();

    JetInfo[iJetColl].Jet_ntracks[JetInfo[iJetColl].nJet] = ntagtracks;

    JetInfo[iJetColl].Jet_nFirstTrack[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nTrack;
    JetInfo[iJetColl].Jet_nFirstTrackTruth[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nTrackTruth;
    JetInfo[iJetColl].Jet_nFirstTrkInc[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkInc;

    unsigned int trackSize = selectedTracks.size();
    if ( !useSelectedTracks_ ) trackSize = tracks.size();

    for (unsigned int itt=0; itt < trackSize; ++itt)
    {
      const TrackRef ptrackRef = ( useSelectedTracks_ ? selectedTracks[itt] : tracks[itt]);
      const reco::Track * ptrackPtr = reco::btag::toTrack(ptrackRef);
      const reco::Track & ptrack = *ptrackPtr;

      reco::TransientTrack transientTrack = trackBuilder->build(ptrackRef);
      GlobalVector direction(pjet->px(), pjet->py(), pjet->pz());

      //--------------------------------
      Double_t decayLength=-1;
      TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), *pv, direction, transientTrack.field());
      if (closest.isValid())
        decayLength =  (closest.globalPosition() - RecoVertex::convertPos(pv->position())).mag();
      else
        decayLength = -1;

      Double_t distJetAxis =  IPTools::jetTrackDistance(transientTrack, direction, *pv).second.value();

      JetInfo[iJetColl].Track_dist[JetInfo[iJetColl].nTrack]     = distJetAxis;
      JetInfo[iJetColl].Track_length[JetInfo[iJetColl].nTrack]   = decayLength;

      JetInfo[iJetColl].Track_dxy[JetInfo[iJetColl].nTrack]      = ptrack.dxy(pv->position());
      JetInfo[iJetColl].Track_dz[JetInfo[iJetColl].nTrack]       = ptrack.dz(pv->position());
      JetInfo[iJetColl].Track_dxyError[JetInfo[iJetColl].nTrack]      = ptrack.dxyError();
      JetInfo[iJetColl].Track_dzError[JetInfo[iJetColl].nTrack]       = ptrack.dzError();
       
	 {
	    TransverseImpactPointExtrapolator extrapolator(transientTrack.field());
	    TrajectoryStateOnSurface closestOnTransversePlaneState =
	      extrapolator.extrapolate(transientTrack.impactPointState(),RecoVertex::convertPos(pv->position()));
	    GlobalPoint impactPoint    = closestOnTransversePlaneState.globalPosition();
	    GlobalVector IPVec(impactPoint.x()-pv->x(),impactPoint.y()-pv->y(),0.);
	    double prod = IPVec.dot(direction);
	    int sign = (prod>=0) ? 1 : -1;
	    JetInfo[iJetColl].Track_sign2D[JetInfo[iJetColl].nTrack]      = sign;
	 }       
	 {
	    AnalyticalImpactPointExtrapolator extrapolator(transientTrack.field());
	    TrajectoryStateOnSurface closestIn3DSpaceState =
	      extrapolator.extrapolate(transientTrack.impactPointState(),RecoVertex::convertPos(pv->position()));
	    GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
	    GlobalVector IPVec(impactPoint.x()-pv->x(),impactPoint.y()-pv->y(),impactPoint.z()-pv->z());
	    double prod = IPVec.dot(direction);
	    int sign = (prod>=0) ? 1 : -1;
	    JetInfo[iJetColl].Track_sign3D[JetInfo[iJetColl].nTrack]      = sign;
	 }       

      float deltaR = reco::deltaR( ptrackRef->eta(), ptrackRef->phi(),
                                   JetInfo[iJetColl].Jet_eta[JetInfo[iJetColl].nJet], JetInfo[iJetColl].Jet_phi[JetInfo[iJetColl].nJet] );

      bool pass_cut_trk = false;
      if (std::fabs(distJetAxis) < _distJetAxis && decayLength < _decayLength) pass_cut_trk = true;

      if (std::fabs(distJetAxis) < _distJetAxis && decayLength < _decayLength
                                        && deltaR < _deltaR) nseltracks++;

      if ( runFatJets_ && runSubJets_ && iJetColl == 0 && pass_cut_trk )
      {
        for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
        {
          if ( subjets[SubJetLabels_[i]].size()<2 ) continue;

          float dR1 = reco::deltaR( ptrackRef->eta(), ptrackRef->phi(),
                                    subjets[SubJetLabels_[i]].at(0)->eta(), subjets[SubJetLabels_[i]].at(0)->phi() );

          float dR2 = reco::deltaR( ptrackRef->eta(), ptrackRef->phi(),
                                    subjets[SubJetLabels_[i]].at(1)->eta(), subjets[SubJetLabels_[i]].at(1)->phi() );

          if ( dR1 < deltaRSubJets_ && dR2 < deltaRSubJets_ ) nsharedtracks[SubJetLabels_[i]]++;
        }
      }

      // track selection
      if ( (useSelectedTracks_ && pass_cut_trk) ||  !useSelectedTracks_) {

        if ( useSelectedTracks_ ) {
          JetInfo[iJetColl].Track_IP2D[JetInfo[iJetColl].nTrack]     = ipTagInfo->impactParameterData()[itt].ip2d.value();
          JetInfo[iJetColl].Track_IP2Dsig[JetInfo[iJetColl].nTrack]  = ipTagInfo->impactParameterData()[itt].ip2d.significance();
          JetInfo[iJetColl].Track_IP[JetInfo[iJetColl].nTrack]       = ipTagInfo->impactParameterData()[itt].ip3d.value();
          JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack]    = ipTagInfo->impactParameterData()[itt].ip3d.significance();
          JetInfo[iJetColl].Track_IP2Derr[JetInfo[iJetColl].nTrack]  = ipTagInfo->impactParameterData()[itt].ip2d.error();
          JetInfo[iJetColl].Track_IPerr[JetInfo[iJetColl].nTrack]    = ipTagInfo->impactParameterData()[itt].ip3d.error();
          JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]    = ipTagInfo->probabilities(0)[itt];
        }
        else {
          Measurement1D ip2d    = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pv).second;
          Measurement1D ip3d    = IPTools::signedImpactParameter3D(transientTrack, direction, *pv).second;
 
          JetInfo[iJetColl].Track_IP2D[JetInfo[iJetColl].nTrack]     = (ip2d.value());
          JetInfo[iJetColl].Track_IP2Dsig[JetInfo[iJetColl].nTrack]  = (ip2d.significance());
          JetInfo[iJetColl].Track_IP[JetInfo[iJetColl].nTrack]       = (ip3d.value());
          JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack]    = (ip3d.significance());
          JetInfo[iJetColl].Track_IP2Derr[JetInfo[iJetColl].nTrack]  = (ip2d.error());
          JetInfo[iJetColl].Track_IPerr[JetInfo[iJetColl].nTrack]    = (ip3d.error());
          JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]    = -1000;
        }

        JetInfo[iJetColl].Track_p[JetInfo[iJetColl].nTrack]        = ptrack.p();
        JetInfo[iJetColl].Track_pt[JetInfo[iJetColl].nTrack]       = ptrack.pt();
        JetInfo[iJetColl].Track_eta[JetInfo[iJetColl].nTrack]      = ptrack.eta();
        JetInfo[iJetColl].Track_phi[JetInfo[iJetColl].nTrack]      = ptrack.phi();
        JetInfo[iJetColl].Track_chi2[JetInfo[iJetColl].nTrack]     = ptrack.normalizedChi2();
        JetInfo[iJetColl].Track_charge[JetInfo[iJetColl].nTrack]   = ptrack.charge();

        JetInfo[iJetColl].Track_nHitAll[JetInfo[iJetColl].nTrack]  = ptrack.numberOfValidHits();
        JetInfo[iJetColl].Track_nHitPixel[JetInfo[iJetColl].nTrack]= ptrack.hitPattern().numberOfValidPixelHits();
        JetInfo[iJetColl].Track_nHitStrip[JetInfo[iJetColl].nTrack]= ptrack.hitPattern().numberOfValidStripHits();
        JetInfo[iJetColl].Track_nHitTIB[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTIBHits();
        JetInfo[iJetColl].Track_nHitTID[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTIDHits();
        JetInfo[iJetColl].Track_nHitTOB[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTOBHits();
        JetInfo[iJetColl].Track_nHitTEC[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTECHits();
        JetInfo[iJetColl].Track_nHitPXB[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidPixelBarrelHits();
        JetInfo[iJetColl].Track_nHitPXF[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidPixelEndcapHits();
        JetInfo[iJetColl].Track_isHitL1[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().hasValidHitInFirstPixelBarrel();

        setTracksPV(ptrackRef, primaryVertex,
                    JetInfo[iJetColl].Track_PV[JetInfo[iJetColl].nTrack],
                    JetInfo[iJetColl].Track_PVweight[JetInfo[iJetColl].nTrack]);

        if(JetInfo[iJetColl].Track_PV[JetInfo[iJetColl].nTrack]==0 &&
           JetInfo[iJetColl].Track_PVweight[JetInfo[iJetColl].nTrack]>0.5) { allKinematics.add(ptrackRef); }


        if( pjet->hasTagInfo(svTagInfos_.c_str()) )
        {
          setTracksSV(ptrackRef, svTagInfo,
                      JetInfo[iJetColl].Track_isfromSV[JetInfo[iJetColl].nTrack],
                      JetInfo[iJetColl].Track_SV[JetInfo[iJetColl].nTrack],
                      JetInfo[iJetColl].Track_SVweight[JetInfo[iJetColl].nTrack]);
        }
        else
        {
          JetInfo[iJetColl].Track_isfromSV[JetInfo[iJetColl].nTrack] = 0;
          JetInfo[iJetColl].Track_SV[JetInfo[iJetColl].nTrack] = -1;
          JetInfo[iJetColl].Track_SVweight[JetInfo[iJetColl].nTrack] = 0.;
        }

        // check if the track is a V0 decay product candidate
        JetInfo[iJetColl].Track_isfromV0[JetInfo[iJetColl].nTrack] = 0;
        // apply the V0 filter
        std::vector<TrackRef> trackPairV0Test(2);
        trackPairV0Test[0] = ptrackRef;

        for (unsigned int jtt=0; jtt < trackSize; ++jtt)
        {
          if (itt == jtt) continue;

          const TrackRef pairTrackRef = ( useSelectedTracks_ ? selectedTracks[jtt] : tracks[jtt]);

          trackPairV0Test[1] = pairTrackRef;

          if (!trackPairV0Filter(trackPairV0Test))
          {
            JetInfo[iJetColl].Track_isfromV0[JetInfo[iJetColl].nTrack] = 1;
            break;
          }
        }

        // compute decay length and distance to tau axis
        if ( runFatJets_ && iJetColl == 0 )
        {
          reco::TransientTrack transientTrack = trackBuilder->build(ptrack);
          GlobalVector direction(pjet->px(), pjet->py(), pjet->pz());

          if (currentAxes.size() > 1)
          {
            if (reco::deltaR2(ptrackRef->momentum(),currentAxes[1]) < reco::deltaR2(ptrackRef->momentum(),currentAxes[0]))
              direction = GlobalVector(currentAxes[1].px(), currentAxes[1].py(), currentAxes[1].pz());
            else
              direction = GlobalVector(currentAxes[0].px(), currentAxes[0].py(), currentAxes[0].pz());
          }
          else if (currentAxes.size() > 0)
            direction = GlobalVector(currentAxes[0].px(), currentAxes[0].py(), currentAxes[0].pz());

          float decayLengthTau=-1;
          float distTauAxis=-1;

          TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), *pv, direction, transientTrack.field());
          if (closest.isValid())
            decayLengthTau =  (closest.globalPosition() - RecoVertex::convertPos(pv->position())).mag();

          distTauAxis = std::abs(IPTools::jetTrackDistance(transientTrack, direction, *pv).second.value());

          JetInfo[iJetColl].Track_distTau[JetInfo[iJetColl].nTrack]   = distTauAxis;
          JetInfo[iJetColl].Track_lengthTau[JetInfo[iJetColl].nTrack] = decayLengthTau;
        }


        JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] = 0;

        if ( useTrackHistory_ && !isData_ ) {
           TrackCategories::Flags theFlag ;
           if(useSelectedTracks_) theFlag  = classifier_.evaluate( toTrackRef(ipTagInfo->selectedTracks()[itt]) ).flags();
	       else                   theFlag  = classifier_.evaluate( toTrackRef(tracks[itt]) ).flags();

           if ( theFlag[TrackCategories::BWeakDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 1);
           if ( theFlag[TrackCategories::CWeakDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 2);
           if ( theFlag[TrackCategories::SignalEvent] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 3);
           if ( theFlag[TrackCategories::ConversionsProcess] )JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 4);
           if ( theFlag[TrackCategories::KsDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 5);
           if ( theFlag[TrackCategories::LambdaDecay] )       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 6);
           if ( theFlag[TrackCategories::HadronicProcess] )   JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 7);
           if ( theFlag[TrackCategories::Fake] ) 	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 8);
		   if ( theFlag[TrackCategories::SharedInnerHits] )   JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 9);
        
        
            //********************************************************************
		    //
			//	Match track to TrackingParticle and retrieve TrackTruth info
			//
			//********************************************************************
			if (produceJetTrackTruthTree_){
			
				std::pair<TrackingParticleRef, double> res;
				TrackingParticleRef tpr;
				double quality_tpr;
				res = classifier_.history().getMatchedTrackingParticle();
				tpr = res.first;
				quality_tpr = res.second;
	
				JetInfo[iJetColl].Track_TPAssociationQuality[JetInfo[iJetColl].nTrack] = quality_tpr;
				JetInfo[iJetColl].Track_idxMatchedTP[JetInfo[iJetColl].nTrack] = -99; // default in case no match was found
	
				// Match TP to hit-cluster (re-ordering according to TP rather than clusters and look for equal_range of a given tpr)
				auto clusterTPmap = clusterToTPMap->map();
				std::sort(clusterTPmap.begin(), clusterTPmap.end(), compare);
				auto clusterRange = std::equal_range(clusterTPmap.begin(), clusterTPmap.end(),std::make_pair(OmniClusterRef(), tpr), compare);
	
				if (quality_tpr != 0) {
	
					// to match Track to TP
					JetInfo[iJetColl].Track_idxMatchedTP[JetInfo[iJetColl].nTrack] = JetInfo[iJetColl].nTrackTruth;
					// to match TP to Track
					JetInfo[iJetColl].TrackTruth_idxMatchedTrack[JetInfo[iJetColl].nTrackTruth] = JetInfo[iJetColl].nTrack;
		
					JetInfo[iJetColl].TrackTruth_p[JetInfo[iJetColl].nTrackTruth] = tpr->p();
					JetInfo[iJetColl].TrackTruth_pt[JetInfo[iJetColl].nTrackTruth] = tpr->pt();
					JetInfo[iJetColl].TrackTruth_eta[JetInfo[iJetColl].nTrackTruth] = tpr->eta();
					JetInfo[iJetColl].TrackTruth_phi[JetInfo[iJetColl].nTrackTruth] = tpr->phi();
					JetInfo[iJetColl].TrackTruth_charge[JetInfo[iJetColl].nTrackTruth] = tpr->charge();
					JetInfo[iJetColl].TrackTruth_pdgid[JetInfo[iJetColl].nTrackTruth] = tpr->pdgId();
		
					// calculate dxy,dz,IP,...
					TrackingParticle::Point vertex_pv = pv->position();
					TrackingParticle::Point vertex_tpr = tpr->vertex();
					TrackingParticle::Vector momentum_tpr = tpr->momentum();
					float dxy_tpr = (-(vertex_tpr.x()-vertex_pv.x())*momentum_tpr.y()+(vertex_tpr.y()-vertex_pv.y())*momentum_tpr.x())/tpr->pt();
					float dz_tpr = (vertex_tpr.z()-vertex_pv.z()) - ((vertex_tpr.x()-vertex_pv.x())*momentum_tpr.x()+(vertex_tpr.y()-vertex_pv.y())*momentum_tpr.y())/sqrt(momentum_tpr.perp2()) * momentum_tpr.z()/sqrt(momentum_tpr.perp2());
					JetInfo[iJetColl].TrackTruth_dxy[JetInfo[iJetColl].nTrackTruth] = dxy_tpr;
					JetInfo[iJetColl].TrackTruth_dz[JetInfo[iJetColl].nTrackTruth] = dz_tpr;
		
					// calculate hits
					int n_pix_hits = 0;
					int n_strip_hits = 0;
					if( clusterRange.first != clusterRange.second ) {
						for( auto ip=clusterRange.first; ip != clusterRange.second; ++ip ) {
							const OmniClusterRef& cluster = ip->first;
							if (cluster.isPixel() && cluster.isValid()){ n_pix_hits+=1;}
							if (cluster.isStrip() && cluster.isValid()){ n_strip_hits+=1;}
						}
					}
					JetInfo[iJetColl].TrackTruth_nHitAll[JetInfo[iJetColl].nTrackTruth] = n_pix_hits+n_strip_hits;
					JetInfo[iJetColl].TrackTruth_nHitStrip[JetInfo[iJetColl].nTrackTruth] = n_strip_hits;
					JetInfo[iJetColl].TrackTruth_nHitPixel[JetInfo[iJetColl].nTrackTruth] = n_pix_hits;


					++JetInfo[iJetColl].nTrackTruth;
	
				}
        	}
        	// ************ end of track truth calculations ****************************
        }
        // ************ end of track history calculations ****************************

        JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = -1;

        TLorentzVector track4P, jet4P;
        track4P.SetPtEtaPhiM(ptrack.pt(), ptrack.eta(), ptrack.phi(), 0. );
        jet4P.SetPtEtaPhiM( pjet->pt(), pjet->eta(), pjet->phi(), pjet->energy() );


        if ( jet4P.DeltaR(track4P) < _deltaR && std::fabs(distJetAxis) < _distJetAxis && decayLength < _decayLength ) {

          if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) Histos[iJetColl]->TrackProbaNeg->Fill(-JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]);
          if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 && PFJet80 ) Histos[iJetColl]->TrackProbJet80->Fill(-JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]);

          if ( findCat( &ptrack, *&cat0 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 0;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat0->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat0->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat0->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat0->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat1 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack]  = 1;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat1->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat1->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat1->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat1->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat2 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 2;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat2->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat2->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat2->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat2->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat3 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 3;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat3->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat3->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat3->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat3->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat4 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 4;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat4->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat4->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat4->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat4->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat5 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 5;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat5->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat5->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat5->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat5->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat6 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 6;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat6->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat6->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat6->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat6->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat7 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 7;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat7->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat7->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat7->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat7->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat8 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 8;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat8->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat8->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat8->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat8->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          else if ( findCat( &ptrack, *&cat9 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 9;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat9->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat9->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat9->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat9->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }
        }

        ++JetInfo[iJetColl].nTrack;
      }
      if ( useSelectedTracks_ ) JetInfo[iJetColl].Jet_ntracks[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrack-JetInfo[iJetColl].Jet_nFirstTrack[JetInfo[iJetColl].nJet];

      if ( producePtRelTemplate_ ) {
        if ( ptrack.quality(reco::TrackBase::highPurity)
            // Remove the tracks that are pixel-less
            && ptrack.algo()!=(reco::TrackBase::pixelLessStep)
            && ptrack.algo()!=(reco::TrackBase::tobTecStep)
            && deltaR < 0.4
            && ptrack.pt() > 5.
            && ptrack.numberOfValidHits() >= 11
            && ptrack.hitPattern().numberOfValidPixelHits() >= 2
            && ptrack.normalizedChi2() < 10
            && ptrack.hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS) <= 2
            && ptrack.dz( pv->position() ) < 1. ) {

          if ( useSelectedTracks_ ) {
            JetInfo[iJetColl].TrkInc_IP[JetInfo[iJetColl].nTrkInc]    = ipTagInfo->impactParameterData()[itt].ip3d.value();
            JetInfo[iJetColl].TrkInc_IPsig[JetInfo[iJetColl].nTrkInc] = ipTagInfo->impactParameterData()[itt].ip3d.significance();
          }
          JetInfo[iJetColl].TrkInc_pt[JetInfo[iJetColl].nTrkInc]    = ptrack.pt();
          JetInfo[iJetColl].TrkInc_eta[JetInfo[iJetColl].nTrkInc]   = ptrack.eta();
          JetInfo[iJetColl].TrkInc_phi[JetInfo[iJetColl].nTrkInc]   = ptrack.phi();
          JetInfo[iJetColl].TrkInc_ptrel[JetInfo[iJetColl].nTrkInc] = calculPtRel( ptrack , *pjet);

          ++JetInfo[iJetColl].nTrkInc;
        }
      }

    } //// end loop on tracks



    JetInfo[iJetColl].Jet_nseltracks[JetInfo[iJetColl].nJet] = nseltracks;

    if ( runFatJets_ && runSubJets_ && iJetColl == 0 )
    {
      for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
        SubJetInfo[SubJetLabels_[i]].Jet_nsharedtracks[JetInfo[iJetColl].nJet] = nsharedtracks[SubJetLabels_[i]];
    }

    JetInfo[iJetColl].Jet_nLastTrack[JetInfo[iJetColl].nJet]   = JetInfo[iJetColl].nTrack;
    JetInfo[iJetColl].Jet_nLastTrackTruth[JetInfo[iJetColl].nJet]   = JetInfo[iJetColl].nTrackTruth;
    JetInfo[iJetColl].Jet_nLastTrkInc[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nTrkInc;

    math::XYZTLorentzVector allSum = allKinematics.weightedVectorSum() ; // allKinematics.vectorSum()

    // PFMuon information
    int nSM = (pjet->hasTagInfo(softPFMuonTagInfos_.c_str()) ? softPFMuTagInfo->leptons() : 0);
    JetInfo[iJetColl].Jet_nFirstSM[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nPFMuon;

    for (size_t leptIdx = 0; leptIdx < (size_t)nSM; ++leptIdx) {

      JetInfo[iJetColl].PFMuon_IdxJet[JetInfo[iJetColl].nPFMuon]    = JetInfo[iJetColl].nJet;
      JetInfo[iJetColl].PFMuon_pt[JetInfo[iJetColl].nPFMuon]        = softPFMuTagInfo->lepton(leptIdx)->pt();
      JetInfo[iJetColl].PFMuon_eta[JetInfo[iJetColl].nPFMuon]       = softPFMuTagInfo->lepton(leptIdx)->eta();
      JetInfo[iJetColl].PFMuon_phi[JetInfo[iJetColl].nPFMuon]       = softPFMuTagInfo->lepton(leptIdx)->phi();
      JetInfo[iJetColl].PFMuon_ptrel[JetInfo[iJetColl].nPFMuon]     = (softPFMuTagInfo->properties(leptIdx).ptRel);
      JetInfo[iJetColl].PFMuon_ratio[JetInfo[iJetColl].nPFMuon]     = (softPFMuTagInfo->properties(leptIdx).ratio);
      JetInfo[iJetColl].PFMuon_ratioRel[JetInfo[iJetColl].nPFMuon]  = (softPFMuTagInfo->properties(leptIdx).ratioRel);
      JetInfo[iJetColl].PFMuon_deltaR[JetInfo[iJetColl].nPFMuon]    = (softPFMuTagInfo->properties(leptIdx).deltaR);

      JetInfo[iJetColl].PFMuon_IP[JetInfo[iJetColl].nPFMuon]        = (softPFMuTagInfo->properties(leptIdx).sip3d);
      JetInfo[iJetColl].PFMuon_IP2D[JetInfo[iJetColl].nPFMuon]      = (softPFMuTagInfo->properties(leptIdx).sip2d);
      JetInfo[iJetColl].PFMuon_IPsig[JetInfo[iJetColl].nPFMuon]        = (softPFMuTagInfo->properties(leptIdx).sip3dsig);
      JetInfo[iJetColl].PFMuon_IP2Dsig[JetInfo[iJetColl].nPFMuon]      = (softPFMuTagInfo->properties(leptIdx).sip2dsig);

      JetInfo[iJetColl].PFMuon_nMuHit[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_nTkHit[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_nPixHit[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_nOutHit[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_nTkLwM[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_nPixLwM[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_nMatched[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_chi2[JetInfo[iJetColl].nPFMuon] = 99;
      JetInfo[iJetColl].PFMuon_chi2Tk[JetInfo[iJetColl].nPFMuon]= 99;
      JetInfo[iJetColl].PFMuon_isGlobal[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] = 0;
      JetInfo[iJetColl].PFMuon_dz[JetInfo[iJetColl].nPFMuon] = 99;
      JetInfo[iJetColl].PFMuon_GoodQuality[JetInfo[iJetColl].nPFMuon] = 0;
      
      const edm::Ptr<reco::Muon> muonPtr = matchMuon( softPFMuTagInfo->lepton(leptIdx), muons );
      if ( muonPtr.isNonnull() && muonPtr.isAvailable() && muonPtr->isGlobalMuon() ) {

        JetInfo[iJetColl].PFMuon_nMuHit[JetInfo[iJetColl].nPFMuon] = muonPtr->outerTrack()->hitPattern().numberOfValidMuonHits();
        JetInfo[iJetColl].PFMuon_nTkHit[JetInfo[iJetColl].nPFMuon] = muonPtr->innerTrack()->hitPattern().numberOfValidHits();
        JetInfo[iJetColl].PFMuon_nPixHit[JetInfo[iJetColl].nPFMuon] = muonPtr->innerTrack()->hitPattern().numberOfValidPixelHits();
        JetInfo[iJetColl].PFMuon_nOutHit[JetInfo[iJetColl].nPFMuon] = muonPtr->innerTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS);
        JetInfo[iJetColl].PFMuon_nTkLwM[JetInfo[iJetColl].nPFMuon] = muonPtr->innerTrack()->hitPattern().trackerLayersWithMeasurement();
        JetInfo[iJetColl].PFMuon_nPixLwM[JetInfo[iJetColl].nPFMuon] = muonPtr->innerTrack()->hitPattern().pixelLayersWithMeasurement();
        JetInfo[iJetColl].PFMuon_nMatched[JetInfo[iJetColl].nPFMuon] = muonPtr->numberOfMatchedStations();
        JetInfo[iJetColl].PFMuon_chi2[JetInfo[iJetColl].nPFMuon] = muonPtr->globalTrack()->normalizedChi2();
        JetInfo[iJetColl].PFMuon_chi2Tk[JetInfo[iJetColl].nPFMuon]= muonPtr->innerTrack()->normalizedChi2();
        JetInfo[iJetColl].PFMuon_isGlobal[JetInfo[iJetColl].nPFMuon] = 1;
        JetInfo[iJetColl].PFMuon_dz[JetInfo[iJetColl].nPFMuon] = muonPtr->muonBestTrack()->dz(pv->position());
        JetInfo[iJetColl].PFMuon_GoodQuality[JetInfo[iJetColl].nPFMuon] = 1;

        if (muonPtr->outerTrack()->hitPattern().numberOfValidMuonHits()>0 &&
            muonPtr->numberOfMatches()>1 && muonPtr->innerTrack()->hitPattern().numberOfValidHits()>10 &&
	    muonPtr->innerTrack()->hitPattern().numberOfValidPixelHits()>1 &&
	    muonPtr->innerTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS)<3 &&
	    muonPtr->globalTrack()->normalizedChi2()<10. && muonPtr->innerTrack()->normalizedChi2()<10.)
          JetInfo[iJetColl].PFMuon_GoodQuality[JetInfo[iJetColl].nPFMuon] = 2;

        if ( useTrackHistory_ && !isData_ ) {
          TrackCategories::Flags theFlagP = classifier_.evaluate( toTrackRef(muonPtr->globalTrack()) ).flags();
          if ( theFlagP[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 1));
          if ( theFlagP[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 2));
          if ( theFlagP[TrackCategories::SignalEvent] )           JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 3));
          if ( theFlagP[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 4));
          if ( theFlagP[TrackCategories::KsDecay] )            JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 5));
          if ( theFlagP[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 6));
          if ( theFlagP[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 7));
          if ( theFlagP[TrackCategories::Fake] )               JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 8));
          if ( theFlagP[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].PFMuon_hist[JetInfo[iJetColl].nPFMuon] += int(pow(10., -1 + 9));
        }
      }

      ++JetInfo[iJetColl].nPFMuon;
    }
    JetInfo[iJetColl].Jet_nLastSM[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nPFMuon;
    JetInfo[iJetColl].Jet_nSM[JetInfo[iJetColl].nJet] = nSM;

    // PFElectron information
    int nSE = (pjet->hasTagInfo(softPFElectronTagInfos_.c_str()) ? softPFElTagInfo->leptons() : 0);
    JetInfo[iJetColl].Jet_nFirstSE[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nPFElectron;

    for (size_t leptIdx = 0; leptIdx < (size_t)nSE; ++leptIdx) {

      JetInfo[iJetColl].PFElectron_IdxJet[JetInfo[iJetColl].nPFElectron]    = JetInfo[iJetColl].nJet;
      JetInfo[iJetColl].PFElectron_pt[JetInfo[iJetColl].nPFElectron]        = softPFElTagInfo->lepton(leptIdx)->pt();
      JetInfo[iJetColl].PFElectron_eta[JetInfo[iJetColl].nPFElectron]       = softPFElTagInfo->lepton(leptIdx)->eta();
      JetInfo[iJetColl].PFElectron_phi[JetInfo[iJetColl].nPFElectron]       = softPFElTagInfo->lepton(leptIdx)->phi();
      JetInfo[iJetColl].PFElectron_ptrel[JetInfo[iJetColl].nPFElectron]     = (softPFElTagInfo->properties(leptIdx).ptRel);
      JetInfo[iJetColl].PFElectron_ratio[JetInfo[iJetColl].nPFElectron]     = (softPFElTagInfo->properties(leptIdx).ratio);
      JetInfo[iJetColl].PFElectron_ratioRel[JetInfo[iJetColl].nPFElectron]  = (softPFElTagInfo->properties(leptIdx).ratioRel);
      JetInfo[iJetColl].PFElectron_deltaR[JetInfo[iJetColl].nPFElectron]    = (softPFElTagInfo->properties(leptIdx).deltaR);
      JetInfo[iJetColl].PFElectron_IP[JetInfo[iJetColl].nPFElectron]        = (softPFElTagInfo->properties(leptIdx).sip3d);
      JetInfo[iJetColl].PFElectron_IP2D[JetInfo[iJetColl].nPFElectron]      = (softPFElTagInfo->properties(leptIdx).sip2d);
      JetInfo[iJetColl].PFElectron_mva_e_pi[JetInfo[iJetColl].nPFElectron]  = (softPFElTagInfo->properties(leptIdx).elec_mva);

      ++JetInfo[iJetColl].nPFElectron;
    }
    JetInfo[iJetColl].Jet_nLastSE[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nPFElectron;
    JetInfo[iJetColl].Jet_nSE[JetInfo[iJetColl].nJet] = nSE;


    // b-tagger discriminants
    float Proba  = pjet->bDiscriminator(jetPBJetTags_.c_str());
    float ProbaN = pjet->bDiscriminator(jetPNegBJetTags_.c_str());
    float ProbaP = pjet->bDiscriminator(jetPPosBJetTags_.c_str());

    float Bprob  = pjet->bDiscriminator(jetBPBJetTags_.c_str());
    float BprobN = pjet->bDiscriminator(jetBPNegBJetTags_.c_str());
    float BprobP = pjet->bDiscriminator(jetBPPosBJetTags_.c_str());

    float CombinedSvtx  = pjet->bDiscriminator(combinedSVBJetTags_.c_str());
    float CombinedSvtxN = pjet->bDiscriminator(combinedSVNegBJetTags_.c_str());
    float CombinedSvtxP = pjet->bDiscriminator(combinedSVPosBJetTags_.c_str());

    float DeepCSVb   = (deepCSVBJetTags_.size()) ? pjet->bDiscriminator((deepCSVBJetTags_+":probb"   ).c_str()) : -10;
    float DeepCSVc   = (deepCSVBJetTags_.size()) ? pjet->bDiscriminator((deepCSVBJetTags_+":probc"   ).c_str()) : -10;
    float DeepCSVl   = (deepCSVBJetTags_.size()) ? pjet->bDiscriminator((deepCSVBJetTags_+":probudsg").c_str()) : -10;
    float DeepCSVbb  = (deepCSVBJetTags_.size()) ? pjet->bDiscriminator((deepCSVBJetTags_+":probbb"  ).c_str()) : -10;
    float DeepCSVcc  = (deepCSVBJetTags_.size()) ? pjet->bDiscriminator((deepCSVBJetTags_+":probcc"  ).c_str()) : -10;
    float DeepCSVbN  = (deepCSVNegBJetTags_.size()) ? pjet->bDiscriminator((deepCSVNegBJetTags_+":probb"   ).c_str()) : -10;
    float DeepCSVcN  = (deepCSVNegBJetTags_.size()) ? pjet->bDiscriminator((deepCSVNegBJetTags_+":probc"   ).c_str()) : -10;
    float DeepCSVlN  = (deepCSVNegBJetTags_.size()) ? pjet->bDiscriminator((deepCSVNegBJetTags_+":probudsg").c_str()) : -10;
    float DeepCSVbbN = (deepCSVNegBJetTags_.size()) ? pjet->bDiscriminator((deepCSVNegBJetTags_+":probbb"  ).c_str()) : -10;
    float DeepCSVccN = (deepCSVNegBJetTags_.size()) ? pjet->bDiscriminator((deepCSVNegBJetTags_+":probcc"  ).c_str()) : -10;
    float DeepCSVbP  = (deepCSVPosBJetTags_.size()) ? pjet->bDiscriminator((deepCSVPosBJetTags_+":probb"   ).c_str()) : -10;
    float DeepCSVcP  = (deepCSVPosBJetTags_.size()) ? pjet->bDiscriminator((deepCSVPosBJetTags_+":probc"   ).c_str()) : -10;
    float DeepCSVlP  = (deepCSVPosBJetTags_.size()) ? pjet->bDiscriminator((deepCSVPosBJetTags_+":probudsg").c_str()) : -10;
    float DeepCSVbbP = (deepCSVPosBJetTags_.size()) ? pjet->bDiscriminator((deepCSVPosBJetTags_+":probbb"  ).c_str()) : -10;
    float DeepCSVccP = (deepCSVPosBJetTags_.size()) ? pjet->bDiscriminator((deepCSVPosBJetTags_+":probcc"  ).c_str()) : -10;

    float CombinedIVF     = pjet->bDiscriminator(combinedIVFSVBJetTags_.c_str());
    float CombinedIVF_P   = pjet->bDiscriminator(combinedIVFSVPosBJetTags_.c_str());
    float CombinedIVF_N   = pjet->bDiscriminator(combinedIVFSVNegBJetTags_.c_str());

    float Svtx    = pjet->bDiscriminator(simpleSVHighEffBJetTags_.c_str());
    float SvtxN   = pjet->bDiscriminator(simpleSVNegHighEffBJetTags_.c_str());
    float SvtxHP  = pjet->bDiscriminator(simpleSVHighPurBJetTags_.c_str());
    float SvtxNHP = pjet->bDiscriminator(simpleSVNegHighPurBJetTags_.c_str());

    float SoftM  = pjet->bDiscriminator(softPFMuonBJetTags_.c_str());
    float SoftMN = pjet->bDiscriminator(softPFMuonNegBJetTags_.c_str());
    float SoftMP = pjet->bDiscriminator(softPFMuonPosBJetTags_.c_str());

    float SoftE  = pjet->bDiscriminator(softPFElectronBJetTags_.c_str());
    float SoftEN = pjet->bDiscriminator(softPFElectronNegBJetTags_.c_str());
    float SoftEP = pjet->bDiscriminator(softPFElectronPosBJetTags_.c_str());

    float DoubleSV = pjet->bDiscriminator(doubleSVBJetTags_.c_str());

    float cMVA = pjet->bDiscriminator(cMVABJetTags_.c_str());
    float cMVAv2 = pjet->bDiscriminator(cMVAv2BJetTags_.c_str());
    float cMVAv2Neg = pjet->bDiscriminator(cMVAv2NegBJetTags_.c_str());
    float cMVAv2Pos = pjet->bDiscriminator(cMVAv2PosBJetTags_.c_str());

    float CvsB = pjet->bDiscriminator(CvsBCJetTags_.c_str());
    float CvsBNeg = pjet->bDiscriminator(CvsBNegCJetTags_.c_str());
    float CvsBPos = pjet->bDiscriminator(CvsBPosCJetTags_.c_str());
    float CvsL = pjet->bDiscriminator(CvsLCJetTags_.c_str());
    float CvsLNeg = pjet->bDiscriminator(CvsLNegCJetTags_.c_str());
    float CvsLPos = pjet->bDiscriminator(CvsLPosCJetTags_.c_str());

    // Jet information
    JetInfo[iJetColl].Jet_DeepCSVBDisc[JetInfo[iJetColl].nJet]   = DeepCSVb + DeepCSVbb  ;
    JetInfo[iJetColl].Jet_DeepCSVBDiscN[JetInfo[iJetColl].nJet]  = DeepCSVbN + DeepCSVbbN;
    JetInfo[iJetColl].Jet_DeepCSVBDiscP[JetInfo[iJetColl].nJet]  = DeepCSVbP + DeepCSVbbP;
    JetInfo[iJetColl].Jet_DeepCSVCvsBDisc[JetInfo[iJetColl].nJet]   = (DeepCSVc  != -1) ? (DeepCSVc  + DeepCSVcc )/(1-(DeepCSVl )) : -1;
    JetInfo[iJetColl].Jet_DeepCSVCvsBDiscN[JetInfo[iJetColl].nJet]  = (DeepCSVcN != -1) ? (DeepCSVcN + DeepCSVccN)/(1-(DeepCSVlN)) : -1;
    JetInfo[iJetColl].Jet_DeepCSVCvsBDiscP[JetInfo[iJetColl].nJet]  = (DeepCSVcP != -1) ? (DeepCSVcP + DeepCSVccP)/(1-(DeepCSVlP)) : -1;
    JetInfo[iJetColl].Jet_DeepCSVCvsLDisc[JetInfo[iJetColl].nJet]   = (DeepCSVc  != -1) ? (DeepCSVc  + DeepCSVcc )/(1-(DeepCSVb + DeepCSVbb  )) : -1;
    JetInfo[iJetColl].Jet_DeepCSVCvsLDiscN[JetInfo[iJetColl].nJet]  = (DeepCSVcN != -1) ? (DeepCSVcN + DeepCSVccN)/(1-(DeepCSVbN + DeepCSVbbN)) : -1;
    JetInfo[iJetColl].Jet_DeepCSVCvsLDiscP[JetInfo[iJetColl].nJet]  = (DeepCSVcP != -1) ? (DeepCSVcP + DeepCSVccP)/(1-(DeepCSVbP + DeepCSVbbP)) : -1;

    JetInfo[iJetColl].Jet_DeepCSVb[JetInfo[iJetColl].nJet]   = DeepCSVb  ;
    JetInfo[iJetColl].Jet_DeepCSVc[JetInfo[iJetColl].nJet]   = DeepCSVc  ;
    JetInfo[iJetColl].Jet_DeepCSVl[JetInfo[iJetColl].nJet]   = DeepCSVl  ;
    JetInfo[iJetColl].Jet_DeepCSVbb[JetInfo[iJetColl].nJet]  = DeepCSVbb ;
    JetInfo[iJetColl].Jet_DeepCSVcc[JetInfo[iJetColl].nJet]  = DeepCSVcc ;
    JetInfo[iJetColl].Jet_DeepCSVbN[JetInfo[iJetColl].nJet]  = DeepCSVbN ;
    JetInfo[iJetColl].Jet_DeepCSVcN[JetInfo[iJetColl].nJet]  = DeepCSVcN ;
    JetInfo[iJetColl].Jet_DeepCSVlN[JetInfo[iJetColl].nJet]  = DeepCSVlN ;
    JetInfo[iJetColl].Jet_DeepCSVbbN[JetInfo[iJetColl].nJet] = DeepCSVbbN;
    JetInfo[iJetColl].Jet_DeepCSVccN[JetInfo[iJetColl].nJet] = DeepCSVccN;
    JetInfo[iJetColl].Jet_DeepCSVbP[JetInfo[iJetColl].nJet]  = DeepCSVbP ;
    JetInfo[iJetColl].Jet_DeepCSVcP[JetInfo[iJetColl].nJet]  = DeepCSVcP ;
    JetInfo[iJetColl].Jet_DeepCSVlP[JetInfo[iJetColl].nJet]  = DeepCSVlP ;
    JetInfo[iJetColl].Jet_DeepCSVbbP[JetInfo[iJetColl].nJet] = DeepCSVbbP;
    JetInfo[iJetColl].Jet_DeepCSVccP[JetInfo[iJetColl].nJet] = DeepCSVccP;
    JetInfo[iJetColl].Jet_ProbaN[JetInfo[iJetColl].nJet]   = ProbaN;
    JetInfo[iJetColl].Jet_ProbaP[JetInfo[iJetColl].nJet]   = ProbaP;
    JetInfo[iJetColl].Jet_Proba[JetInfo[iJetColl].nJet]    = Proba;
    JetInfo[iJetColl].Jet_BprobN[JetInfo[iJetColl].nJet]   = BprobN;
    JetInfo[iJetColl].Jet_BprobP[JetInfo[iJetColl].nJet]   = BprobP;
    JetInfo[iJetColl].Jet_Bprob[JetInfo[iJetColl].nJet]    = Bprob;
    JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet]     = SvtxN;
    JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet]      = Svtx;
    JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet]   = SvtxNHP;
    JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet]    = SvtxHP;
    JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] = CombinedSvtxN;
    JetInfo[iJetColl].Jet_CombSvxP[JetInfo[iJetColl].nJet] = CombinedSvtxP;
    JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet]  = CombinedSvtx;
    JetInfo[iJetColl].Jet_CombIVF[JetInfo[iJetColl].nJet]     = CombinedIVF;
    JetInfo[iJetColl].Jet_CombIVF_P[JetInfo[iJetColl].nJet]   = CombinedIVF_P;
    JetInfo[iJetColl].Jet_CombIVF_N[JetInfo[iJetColl].nJet]   = CombinedIVF_N;
    JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet]  = SoftMN;
    JetInfo[iJetColl].Jet_SoftMuP[JetInfo[iJetColl].nJet]  = SoftMP;
    JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet]   = SoftM;
    JetInfo[iJetColl].Jet_SoftElN[JetInfo[iJetColl].nJet]  = SoftEN;
    JetInfo[iJetColl].Jet_SoftElP[JetInfo[iJetColl].nJet]  = SoftEP;
    JetInfo[iJetColl].Jet_SoftEl[JetInfo[iJetColl].nJet]   = SoftE;
    JetInfo[iJetColl].Jet_DoubleSV[JetInfo[iJetColl].nJet] = DoubleSV;
    JetInfo[iJetColl].Jet_cMVA[JetInfo[iJetColl].nJet] = cMVA;
    JetInfo[iJetColl].Jet_cMVAv2[JetInfo[iJetColl].nJet] = cMVAv2;
    JetInfo[iJetColl].Jet_cMVAv2N[JetInfo[iJetColl].nJet] = cMVAv2Neg;
    JetInfo[iJetColl].Jet_cMVAv2P[JetInfo[iJetColl].nJet] = cMVAv2Pos;

    // TagInfo TaggingVariables
    if ( storeTagVariables )
    {
      reco::TaggingVariableList ipVars = ipTagInfo->taggingVariables();
      reco::TaggingVariableList svVars = svTagInfo->taggingVariables();
      int nTracks = ipTagInfo->selectedTracks().size();
      int nSVs = svTagInfo->nVertices();

      // per jet
      JetInfo[iJetColl].TagVar_jetNTracks[JetInfo[iJetColl].nJet]                  = nTracks;
      JetInfo[iJetColl].TagVar_jetNSecondaryVertices[JetInfo[iJetColl].nJet]       = nSVs;
      //-------------
      JetInfo[iJetColl].TagVar_chargedHadronEnergyFraction[JetInfo[iJetColl].nJet] = pjet->chargedHadronEnergyFraction();
      JetInfo[iJetColl].TagVar_neutralHadronEnergyFraction[JetInfo[iJetColl].nJet] = pjet->neutralHadronEnergyFraction();
      JetInfo[iJetColl].TagVar_photonEnergyFraction[JetInfo[iJetColl].nJet]        = pjet->photonEnergyFraction();
      JetInfo[iJetColl].TagVar_electronEnergyFraction[JetInfo[iJetColl].nJet]      = pjet->electronEnergyFraction();
      JetInfo[iJetColl].TagVar_muonEnergyFraction[JetInfo[iJetColl].nJet]          = pjet->muonEnergyFraction();
      JetInfo[iJetColl].TagVar_chargedHadronMultiplicity[JetInfo[iJetColl].nJet]   = pjet->chargedHadronMultiplicity();
      JetInfo[iJetColl].TagVar_neutralHadronMultiplicity[JetInfo[iJetColl].nJet]   = pjet->neutralHadronMultiplicity();
      JetInfo[iJetColl].TagVar_photonMultiplicity[JetInfo[iJetColl].nJet]          = pjet->photonMultiplicity();
      JetInfo[iJetColl].TagVar_electronMultiplicity[JetInfo[iJetColl].nJet]        = pjet->electronMultiplicity();
      JetInfo[iJetColl].TagVar_muonMultiplicity[JetInfo[iJetColl].nJet]            = pjet->muonMultiplicity();

      // per jet per track
      JetInfo[iJetColl].Jet_nFirstTrkTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkTagVar;

      std::vector<float> tagValList = ipVars.getList(reco::btau::trackMomentum,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackMomentum[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackEta,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackEta[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackPhi,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackPhi[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackPtRel,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackPtRel[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackPPar,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackPPar[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackEtaRel,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackEtaRel[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackDeltaR,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackDeltaR[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackPtRatio,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackPtRatio[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackPParRatio,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackPParRatio[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackSip2dVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackSip2dVal[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackSip2dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackSip2dSig[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackSip3dVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackSip3dVal[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackSip3dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackSip3dSig[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackDecayLenVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackDecayLenVal[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackDecayLenSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackDecayLenSig[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackJetDistVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackJetDistVal[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackJetDistSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackJetDistSig[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackChi2,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackChi2[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackNTotalHits,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackNTotalHits[JetInfo[iJetColl].nTrkTagVar] );
      tagValList = ipVars.getList(reco::btau::trackNPixelHits,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_trackNPixelHits[JetInfo[iJetColl].nTrkTagVar] );

      JetInfo[iJetColl].nTrkTagVar += nTracks;
      JetInfo[iJetColl].Jet_nLastTrkTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkTagVar;

      // per jet per secondary vertex
      JetInfo[iJetColl].Jet_nFirstSVTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nSVTagVar;

      for(int svIdx=0; svIdx < nSVs; ++svIdx)
      {
        JetInfo[iJetColl].TagVar_vertexMass[JetInfo[iJetColl].nSVTagVar + svIdx]    = svTagInfo->secondaryVertex(svIdx).p4().mass();
        JetInfo[iJetColl].TagVar_vertexNTracks[JetInfo[iJetColl].nSVTagVar + svIdx] = vtxTracks(svTagInfo->secondaryVertex(svIdx));
      }
      tagValList = svVars.getList(reco::btau::vertexJetDeltaR,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_vertexJetDeltaR[JetInfo[iJetColl].nSVTagVar] );
      tagValList = svVars.getList(reco::btau::flightDistance2dVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_flightDistance2dVal[JetInfo[iJetColl].nSVTagVar] );
      tagValList = svVars.getList(reco::btau::flightDistance2dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_flightDistance2dSig[JetInfo[iJetColl].nSVTagVar] );
      tagValList = svVars.getList(reco::btau::flightDistance3dVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_flightDistance3dVal[JetInfo[iJetColl].nSVTagVar] );
      tagValList = svVars.getList(reco::btau::flightDistance3dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVar_flightDistance3dSig[JetInfo[iJetColl].nSVTagVar] );

      JetInfo[iJetColl].nSVTagVar += nSVs;
      JetInfo[iJetColl].Jet_nLastSVTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nSVTagVar;
    }

    // CSV TaggingVariables
    if ( storeCSVTagVariables )
    {
      std::vector<const reco::BaseTagInfo*>  baseTagInfos;
      JetTagComputer::TagInfoHelper helper(baseTagInfos);
      baseTagInfos.push_back( ipTagInfo );
      baseTagInfos.push_back( svTagInfo );
      // TaggingVariables
      reco::TaggingVariableList vars = computer->taggingVariables(helper);

      // per jet
      JetInfo[iJetColl].TagVarCSV_trackJetPt[JetInfo[iJetColl].nJet]                  = ( vars.checkTag(reco::btau::trackJetPt) ? vars.get(reco::btau::trackJetPt) : -9999 );
      JetInfo[iJetColl].TagVarCSV_vertexCategory[JetInfo[iJetColl].nJet]              = ( vars.checkTag(reco::btau::vertexCategory) ? vars.get(reco::btau::vertexCategory) : -9999 );
      JetInfo[iJetColl].TagVarCSV_jetNSecondaryVertices[JetInfo[iJetColl].nJet]       = ( vars.checkTag(reco::btau::jetNSecondaryVertices) ? vars.get(reco::btau::jetNSecondaryVertices) : 0 );
      JetInfo[iJetColl].TagVarCSV_trackSumJetEtRatio[JetInfo[iJetColl].nJet]          = ( vars.checkTag(reco::btau::trackSumJetEtRatio) ? vars.get(reco::btau::trackSumJetEtRatio) : -9999 );
      JetInfo[iJetColl].TagVarCSV_trackSumJetDeltaR[JetInfo[iJetColl].nJet]           = ( vars.checkTag(reco::btau::trackSumJetDeltaR) ? vars.get(reco::btau::trackSumJetDeltaR) : -9999 );
      JetInfo[iJetColl].TagVarCSV_trackSip2dValAboveCharm[JetInfo[iJetColl].nJet]     = ( vars.checkTag(reco::btau::trackSip2dValAboveCharm) ? vars.get(reco::btau::trackSip2dValAboveCharm) : -9999 );
      JetInfo[iJetColl].TagVarCSV_trackSip2dSigAboveCharm[JetInfo[iJetColl].nJet]     = ( vars.checkTag(reco::btau::trackSip2dSigAboveCharm) ? vars.get(reco::btau::trackSip2dSigAboveCharm) : -9999 );
      JetInfo[iJetColl].TagVarCSV_trackSip3dValAboveCharm[JetInfo[iJetColl].nJet]     = ( vars.checkTag(reco::btau::trackSip3dValAboveCharm) ? vars.get(reco::btau::trackSip3dValAboveCharm) : -9999 );
      JetInfo[iJetColl].TagVarCSV_trackSip3dSigAboveCharm[JetInfo[iJetColl].nJet]     = ( vars.checkTag(reco::btau::trackSip3dSigAboveCharm) ? vars.get(reco::btau::trackSip3dSigAboveCharm) : -9999 );
      JetInfo[iJetColl].TagVarCSV_vertexMass[JetInfo[iJetColl].nJet]                  = ( vars.checkTag(reco::btau::vertexMass) ? vars.get(reco::btau::vertexMass) : -9999 );
      JetInfo[iJetColl].TagVarCSV_vertexNTracks[JetInfo[iJetColl].nJet]               = ( vars.checkTag(reco::btau::vertexNTracks) ? vars.get(reco::btau::vertexNTracks) : 0 );
      JetInfo[iJetColl].TagVarCSV_vertexEnergyRatio[JetInfo[iJetColl].nJet]           = ( vars.checkTag(reco::btau::vertexEnergyRatio) ? vars.get(reco::btau::vertexEnergyRatio) : -9999 );
      JetInfo[iJetColl].TagVarCSV_vertexJetDeltaR[JetInfo[iJetColl].nJet]             = ( vars.checkTag(reco::btau::vertexJetDeltaR) ? vars.get(reco::btau::vertexJetDeltaR) : -9999 );
      JetInfo[iJetColl].TagVarCSV_flightDistance2dVal[JetInfo[iJetColl].nJet]         = ( vars.checkTag(reco::btau::flightDistance2dVal) ? vars.get(reco::btau::flightDistance2dVal) : -9999 );
      JetInfo[iJetColl].TagVarCSV_flightDistance2dSig[JetInfo[iJetColl].nJet]         = ( vars.checkTag(reco::btau::flightDistance2dSig) ? vars.get(reco::btau::flightDistance2dSig) : -9999 );
      JetInfo[iJetColl].TagVarCSV_flightDistance3dVal[JetInfo[iJetColl].nJet]         = ( vars.checkTag(reco::btau::flightDistance3dVal) ? vars.get(reco::btau::flightDistance3dVal) : -9999 );
      JetInfo[iJetColl].TagVarCSV_flightDistance3dSig[JetInfo[iJetColl].nJet]         = ( vars.checkTag(reco::btau::flightDistance3dSig) ? vars.get(reco::btau::flightDistance3dSig) : -9999 );

      // per jet per track
      JetInfo[iJetColl].Jet_nFirstTrkTagVarCSV[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkTagVarCSV;
      std::vector<float> tagValList = vars.getList(reco::btau::trackSip2dSig,false);
      JetInfo[iJetColl].TagVarCSV_jetNTracks[JetInfo[iJetColl].nJet] = tagValList.size();

      tagValList = vars.getList(reco::btau::trackMomentum,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackMomentum[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackEta,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackEta[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackPhi,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackPhi[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackPtRel,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackPtRel[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackPPar,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackPPar[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackDeltaR,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackDeltaR[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackPtRatio,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackPtRatio[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackPParRatio,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackPParRatio[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackSip2dVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackSip2dVal[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackSip2dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackSip2dSig[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackSip3dVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackSip3dVal[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackSip3dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackSip3dSig[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackDecayLenVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackDecayLenVal[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackDecayLenSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackDecayLenSig[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackJetDistVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackJetDistVal[JetInfo[iJetColl].nTrkTagVarCSV] );
      tagValList = vars.getList(reco::btau::trackJetDistSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackJetDistSig[JetInfo[iJetColl].nTrkTagVarCSV] );

      JetInfo[iJetColl].nTrkTagVarCSV += JetInfo[iJetColl].TagVarCSV_jetNTracks[JetInfo[iJetColl].nJet];
      JetInfo[iJetColl].Jet_nLastTrkTagVarCSV[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkTagVarCSV;
      //---------------------------
      JetInfo[iJetColl].Jet_nFirstTrkEtaRelTagVarCSV[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkEtaRelTagVarCSV;
      tagValList = vars.getList(reco::btau::trackEtaRel,false);
      JetInfo[iJetColl].TagVarCSV_jetNTracksEtaRel[JetInfo[iJetColl].nJet] = tagValList.size();

      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].TagVarCSV_trackEtaRel[JetInfo[iJetColl].nTrkEtaRelTagVarCSV] );

      JetInfo[iJetColl].nTrkEtaRelTagVarCSV += JetInfo[iJetColl].TagVarCSV_jetNTracksEtaRel[JetInfo[iJetColl].nJet];
      JetInfo[iJetColl].Jet_nLastTrkEtaRelTagVarCSV[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkEtaRelTagVarCSV;
      
    }

    if ( storeCTagVariables )
    {
      JetInfo[iJetColl].CTag_Jet_CvsB[JetInfo[iJetColl].nJet] = CvsB;
      JetInfo[iJetColl].CTag_Jet_CvsBN[JetInfo[iJetColl].nJet] = CvsBNeg;
      JetInfo[iJetColl].CTag_Jet_CvsBP[JetInfo[iJetColl].nJet] = CvsBPos;
      JetInfo[iJetColl].CTag_Jet_CvsL[JetInfo[iJetColl].nJet] = CvsL;
      JetInfo[iJetColl].CTag_Jet_CvsLN[JetInfo[iJetColl].nJet] = CvsLNeg;
      JetInfo[iJetColl].CTag_Jet_CvsLP[JetInfo[iJetColl].nJet] = CvsLPos;
  
      std::vector<const reco::BaseTagInfo*>  slbaseTagInfos;
      slbaseTagInfos.push_back( ipTagInfoCTag );
      slbaseTagInfos.push_back( svTagInfoCTag );
      slbaseTagInfos.push_back( softPFMuTagInfoCTag );
      slbaseTagInfos.push_back( softPFElTagInfoCTag ); 
      JetTagComputer::TagInfoHelper slhelper(slbaseTagInfos);
      // TaggingVariables
      reco::TaggingVariableList slvars = slcomputer->taggingVariables(slhelper);

      // per jet
      JetInfo[iJetColl].CTag_vertexCategory[JetInfo[iJetColl].nJet]              = ( slvars.checkTag(reco::btau::vertexCategory) ? slvars.get(reco::btau::vertexCategory) : -9999 );
      JetInfo[iJetColl].CTag_jetNSecondaryVertices[JetInfo[iJetColl].nJet]       = ( slvars.checkTag(reco::btau::jetNSecondaryVertices) ? slvars.get(reco::btau::jetNSecondaryVertices) : 0 );
      JetInfo[iJetColl].CTag_trackSumJetEtRatio[JetInfo[iJetColl].nJet]          = ( slvars.checkTag(reco::btau::trackSumJetEtRatio) ? slvars.get(reco::btau::trackSumJetEtRatio) : -9999 );
      JetInfo[iJetColl].CTag_trackSumJetDeltaR[JetInfo[iJetColl].nJet]           = ( slvars.checkTag(reco::btau::trackSumJetDeltaR) ? slvars.get(reco::btau::trackSumJetDeltaR) : -9999 );
      JetInfo[iJetColl].CTag_trackSip2dSigAboveCharm[JetInfo[iJetColl].nJet]     = ( slvars.checkTag(reco::btau::trackSip2dSigAboveCharm) ? slvars.get(reco::btau::trackSip2dSigAboveCharm) : -9999 );
      JetInfo[iJetColl].CTag_trackSip3dSigAboveCharm[JetInfo[iJetColl].nJet]     = ( slvars.checkTag(reco::btau::trackSip3dSigAboveCharm) ? slvars.get(reco::btau::trackSip3dSigAboveCharm) : -9999 );
      JetInfo[iJetColl].CTag_vertexMass[JetInfo[iJetColl].nJet]                  = ( slvars.checkTag(reco::btau::vertexMass) ? slvars.get(reco::btau::vertexMass) : -9999 );
      JetInfo[iJetColl].CTag_vertexNTracks[JetInfo[iJetColl].nJet]               = ( slvars.checkTag(reco::btau::vertexNTracks) ? slvars.get(reco::btau::vertexNTracks) : 0 );
      JetInfo[iJetColl].CTag_vertexEnergyRatio[JetInfo[iJetColl].nJet]           = ( slvars.checkTag(reco::btau::vertexEnergyRatio) ? slvars.get(reco::btau::vertexEnergyRatio) : -9999 );
      JetInfo[iJetColl].CTag_vertexJetDeltaR[JetInfo[iJetColl].nJet]             = ( slvars.checkTag(reco::btau::vertexJetDeltaR) ? slvars.get(reco::btau::vertexJetDeltaR) : -9999 );
      JetInfo[iJetColl].CTag_flightDistance2dSig[JetInfo[iJetColl].nJet]         = ( slvars.checkTag(reco::btau::flightDistance2dSig) ? slvars.get(reco::btau::flightDistance2dSig) : -9999 );
      JetInfo[iJetColl].CTag_flightDistance3dSig[JetInfo[iJetColl].nJet]         = ( slvars.checkTag(reco::btau::flightDistance3dSig) ? slvars.get(reco::btau::flightDistance3dSig) : -9999 );
      JetInfo[iJetColl].CTag_massVertexEnergyFraction[JetInfo[iJetColl].nJet]    = ( slvars.checkTag(reco::btau::massVertexEnergyFraction) ? slvars.get(reco::btau::massVertexEnergyFraction) : -0.1);
      JetInfo[iJetColl].CTag_vertexBoostOverSqrtJetPt[JetInfo[iJetColl].nJet]    = ( slvars.checkTag(reco::btau::vertexBoostOverSqrtJetPt) ? slvars.get(reco::btau::vertexBoostOverSqrtJetPt) : -0.1);
      JetInfo[iJetColl].CTag_vertexLeptonCategory[JetInfo[iJetColl].nJet]        = ( slvars.checkTag(reco::btau::vertexLeptonCategory) ? slvars.get(reco::btau::vertexLeptonCategory) : -1);      
 
      // per jet per track
      JetInfo[iJetColl].Jet_nFirstTrkCTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkCTagVar;
      std::vector<float> tagValList = slvars.getList(reco::btau::trackSip2dSig,false);
      JetInfo[iJetColl].CTag_jetNTracks[JetInfo[iJetColl].nJet] = tagValList.size(); 

      tagValList = slvars.getList(reco::btau::trackPtRel,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackPtRel[JetInfo[iJetColl].nTrkCTagVar] );  
      tagValList = slvars.getList(reco::btau::trackPPar,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackPPar[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackDeltaR,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackDeltaR[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackPtRatio,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackPtRatio[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackPParRatio,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackPParRatio[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackSip2dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackSip2dSig[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackSip3dSig,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackSip3dSig[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackDecayLenVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackDecayLenVal[JetInfo[iJetColl].nTrkCTagVar] );
      tagValList = slvars.getList(reco::btau::trackJetDistVal,false);
      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackJetDistVal[JetInfo[iJetColl].nTrkCTagVar] );
     
      JetInfo[iJetColl].nTrkCTagVar += JetInfo[iJetColl].CTag_jetNTracks[JetInfo[iJetColl].nJet];
      JetInfo[iJetColl].Jet_nLastTrkCTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkCTagVar;
      //---------------------------
      JetInfo[iJetColl].Jet_nFirstTrkEtaRelCTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkEtaRelCTagVar;
      tagValList = slvars.getList(reco::btau::trackEtaRel,false);
      JetInfo[iJetColl].CTag_jetNTracksEtaRel[JetInfo[iJetColl].nJet] = tagValList.size();

      if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), &JetInfo[iJetColl].CTag_trackEtaRel[JetInfo[iJetColl].nTrkEtaRelCTagVar] );

      JetInfo[iJetColl].nTrkEtaRelCTagVar += JetInfo[iJetColl].CTag_jetNTracksEtaRel[JetInfo[iJetColl].nJet];
      JetInfo[iJetColl].Jet_nLastTrkEtaRelCTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkEtaRelCTagVar;
      
      //per jet per lepton
      JetInfo[iJetColl].Jet_nFirstLepCTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nLeptons;
      std::vector<float> sltagValList = slvars.getList(reco::btau::leptonPtRel,false);
      JetInfo[iJetColl].CTag_jetNLeptons[JetInfo[iJetColl].nJet] = sltagValList.size();

      sltagValList = slvars.getList(reco::btau::leptonPtRel,false);
      if(sltagValList.size()>0) std::copy( sltagValList.begin(), sltagValList.end(), &JetInfo[iJetColl].CTag_leptonPtRel[JetInfo[iJetColl].nLeptons] );
      sltagValList = slvars.getList(reco::btau::leptonSip3d,false);
      if(sltagValList.size()>0) std::copy( sltagValList.begin(), sltagValList.end(), &JetInfo[iJetColl].CTag_leptonSip3d[JetInfo[iJetColl].nLeptons] );
      sltagValList = slvars.getList(reco::btau::leptonDeltaR,false);
      if(sltagValList.size()>0) std::copy( sltagValList.begin(), sltagValList.end(), &JetInfo[iJetColl].CTag_leptonDeltaR[JetInfo[iJetColl].nLeptons] );
      sltagValList = slvars.getList(reco::btau::leptonRatioRel,false);
      if(sltagValList.size()>0) std::copy( sltagValList.begin(), sltagValList.end(), &JetInfo[iJetColl].CTag_leptonRatioRel[JetInfo[iJetColl].nLeptons] );
      sltagValList = slvars.getList(reco::btau::leptonEtaRel,false);
      if(sltagValList.size()>0) std::copy( sltagValList.begin(), sltagValList.end(), &JetInfo[iJetColl].CTag_leptonEtaRel[JetInfo[iJetColl].nLeptons] );
      sltagValList = slvars.getList(reco::btau::leptonRatio,false);
      if(sltagValList.size()>0) std::copy( sltagValList.begin(), sltagValList.end(), &JetInfo[iJetColl].CTag_leptonRatio[JetInfo[iJetColl].nLeptons] );     

      JetInfo[iJetColl].nLeptons += JetInfo[iJetColl].CTag_jetNLeptons[JetInfo[iJetColl].nJet];
      JetInfo[iJetColl].Jet_nLastLepCTagVar[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nLeptons;
    }


    std::vector<reco::btag::TrackIPData>  ipdata = ipTagInfo->impactParameterData();
    std::vector<std::size_t> indexes( sortedIndexes(ipdata) );

    TrackCategories::Flags flags1P;
    TrackCategories::Flags flags2P;
    TrackCategories::Flags flags3P;
    TrackCategories::Flags flags1N;
    TrackCategories::Flags flags2N;
    TrackCategories::Flags flags3N;
    int idSize = 0;

    if ( useTrackHistory_ && indexes.size() != 0 && !isData_ ) {

      idSize = indexes.size();
      flags1P = classifier_.evaluate( toTrackRef(selectedTracks[indexes[0]]) ).flags();
      if(idSize > 1) flags2P = classifier_.evaluate( toTrackRef(selectedTracks[indexes[1]]) ).flags();
      if(idSize > 2) flags3P = classifier_.evaluate( toTrackRef(selectedTracks[indexes[2]]) ).flags();
      flags1N = classifier_.evaluate( toTrackRef(selectedTracks[indexes[idSize-1]]) ).flags();
      if(idSize > 1) flags2N = classifier_.evaluate( toTrackRef(selectedTracks[indexes[idSize-2]]) ).flags();
      if(idSize > 2) flags3N = classifier_.evaluate( toTrackRef(selectedTracks[indexes[idSize-3]]) ).flags();
    }

    JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCHEBJetTags_.c_str());
    JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCNegHEBJetTags_.c_str());

    JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCHPBJetTags_.c_str());
    JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCNegHPBJetTags_.c_str());

    //*****************************************************************
    // get track histories for 1st, 2nd and 3rd track (TC)
    //*****************************************************************
    JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] = 0;

    // Track Histosry
    if ( useTrackHistory_ && indexes.size()!=0 && !isData_ ) {
      if ( flags1P[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 1));
      if ( flags1P[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 2));
      if ( flags1P[TrackCategories::SignalEvent] )           JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 3));
      if ( flags1P[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 4));
      if ( flags1P[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 5));
      if ( flags1P[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 6));
      if ( flags1P[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 7));
      if ( flags1P[TrackCategories::Fake] )	          JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 8));
      if ( flags1P[TrackCategories::SharedInnerHits] )	  JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 9));
      if ( idSize > 1 ) {
        if ( flags2P[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 1));
        if ( flags2P[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 2));
        if ( flags2P[TrackCategories::SignalEvent] )           JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 3));
        if ( flags2P[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 4));
        if ( flags2P[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 5));
        if ( flags2P[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 6));
        if ( flags2P[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 7));
        if ( flags2P[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 8));
        if ( flags2P[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 9));
      }
      if ( idSize > 2 ) {
        if ( flags3P[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 1));
        if ( flags3P[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 2));
        if ( flags3P[TrackCategories::SignalEvent] )           JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 3));
        if ( flags3P[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 4));
        if ( flags3P[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 5));
        if ( flags3P[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 6));
        if ( flags3P[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 7));
        if ( flags3P[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 8));
        if ( flags3P[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 9));
      }
      if ( flags1N[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 1));
      if ( flags1N[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 2));
      if ( flags1N[TrackCategories::SignalEvent] )           JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 3));
      if ( flags1N[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 4));
      if ( flags1N[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 5));
      if ( flags1N[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 6));
      if ( flags1N[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 7));
      if ( flags1N[TrackCategories::Fake] )	          JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 8));
      if ( flags1N[TrackCategories::SharedInnerHits] )	  JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 9));
      if ( idSize > 1 ) {
        if ( flags2N[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 1));
        if ( flags2N[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 2));
        if ( flags2N[TrackCategories::SignalEvent] )           JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 3));
        if ( flags2N[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 4));
        if ( flags2N[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 5));
        if ( flags2N[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 6));
        if ( flags2N[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 7));
        if ( flags2N[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 8));
        if ( flags2N[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 9));
      }
      if ( idSize > 2 ) {
        if ( flags3N[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 1));
        if ( flags3N[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 2));
        if ( flags3N[TrackCategories::SignalEvent] )           JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 3));
        if ( flags3N[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 4));
        if ( flags3N[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 5));
        if ( flags3N[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 6));
        if ( flags3N[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 7));
        if ( flags3N[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 8));
        if ( flags3N[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 9));
      }
    }

    //*****************************************************************
    // get track Histories of tracks in jets (for Jet Proba)
    //*****************************************************************
    JetInfo[iJetColl].Jet_histJet[JetInfo[iJetColl].nJet] = 0;

    if(useTrackHistory_ && !isData_){
      const Tracks jetProbTracks( ipTagInfo->selectedTracks() );

      cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0;
      can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0;

      for (unsigned int i=0; i<jetProbTracks.size(); ++i) {
        reco::btag::TrackIPData ip = (ipTagInfo->impactParameterData())[i];

        if ( ip.ip3d.significance() > 0 ) {
          TrackCategories::Flags theFlag = classifier_.evaluate( toTrackRef(jetProbTracks[i]) ).flags();
          if ( theFlag[TrackCategories::BWeakDecay] )	      cap0 = 1;
          if ( theFlag[TrackCategories::CWeakDecay] )	      cap1 = 1;
          if ( theFlag[TrackCategories::SignalEvent] )	      cap2 = 1;
          if ( theFlag[TrackCategories::ConversionsProcess] ) cap3 = 1;
          if ( theFlag[TrackCategories::KsDecay] )	      cap4 = 1;
          if ( theFlag[TrackCategories::LambdaDecay] )        cap5 = 1;
          if ( theFlag[TrackCategories::HadronicProcess] )    cap6 = 1;
          if ( theFlag[TrackCategories::Fake] ) 	      cap7 = 1;
          if ( theFlag[TrackCategories::SharedInnerHits] )    cap8 = 1;
        }
        else {
          TrackCategories::Flags theFlag = classifier_.evaluate( toTrackRef(jetProbTracks[i]) ).flags();
          if ( theFlag[TrackCategories::BWeakDecay] )	      can0 = 2;
          if ( theFlag[TrackCategories::CWeakDecay] )	      can1 = 2;
          if ( theFlag[TrackCategories::SignalEvent] )	      can2 = 2;
          if ( theFlag[TrackCategories::ConversionsProcess] ) can3 = 2;
          if ( theFlag[TrackCategories::KsDecay] )	      can4 = 2;
          if ( theFlag[TrackCategories::LambdaDecay] )        can5 = 2;
          if ( theFlag[TrackCategories::HadronicProcess] )    can6 = 2;
          if ( theFlag[TrackCategories::Fake] ) 	      can7 = 2;
          if ( theFlag[TrackCategories::SharedInnerHits] )    can8 = 2;
        }
      }
      JetInfo[iJetColl].Jet_histJet[JetInfo[iJetColl].nJet] =  cap0+can0 + (cap1+can1)*10 + (cap2+can2)*100
        + (cap3+can3)*1000     + (cap4+can4)*10000
        + (cap5+can5)*100000   + (cap6+can6)*1000000
        + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }

    //*****************************************************************
    // get track histories associated to sec. vertex (for simple SV)
    //*****************************************************************
    JetInfo[iJetColl].Jet_histSvx[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_nFirstSV[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nSV;
    JetInfo[iJetColl].Jet_SV_multi[JetInfo[iJetColl].nJet]  = svTagInfo->nVertices();

    // if secondary vertices present
    reco::TrackKinematics tau1Kinematics;
    reco::TrackKinematics tau2Kinematics;

    size_t tau1_vtx = 0, tau2_vtx = 0;
    JetInfo[iJetColl].Jet_tau1_nSecondaryVertices[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_tau2_nSecondaryVertices[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] = -1;
    JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet] = -1;
    JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet]  = -1;
    JetInfo[iJetColl].Jet_tau2_vertexDeltaR[JetInfo[iJetColl].nJet]  = -1;
    JetInfo[iJetColl].Jet_tau1_vertexNTracks[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_tau2_vertexNTracks[JetInfo[iJetColl].nJet] = 0;

    for (size_t vtx = 0; vtx < (size_t)svTagInfo->nVertices(); ++vtx)
    {
	JetInfo[iJetColl].SV_x[JetInfo[iJetColl].nSV]    = position(svTagInfo->secondaryVertex(vtx)).x();
	JetInfo[iJetColl].SV_y[JetInfo[iJetColl].nSV]    = position(svTagInfo->secondaryVertex(vtx)).y();
	JetInfo[iJetColl].SV_z[JetInfo[iJetColl].nSV]    = position(svTagInfo->secondaryVertex(vtx)).z();
	JetInfo[iJetColl].SV_ex[JetInfo[iJetColl].nSV]   = xError(svTagInfo->secondaryVertex(vtx));
	JetInfo[iJetColl].SV_ey[JetInfo[iJetColl].nSV]   = yError(svTagInfo->secondaryVertex(vtx));
	JetInfo[iJetColl].SV_ez[JetInfo[iJetColl].nSV]   = zError(svTagInfo->secondaryVertex(vtx));
	JetInfo[iJetColl].SV_chi2[JetInfo[iJetColl].nSV] = chi2(svTagInfo->secondaryVertex(vtx));
	JetInfo[iJetColl].SV_ndf[JetInfo[iJetColl].nSV]  = ndof(svTagInfo->secondaryVertex(vtx));

	JetInfo[iJetColl].SV_flight[JetInfo[iJetColl].nSV]      = svTagInfo->flightDistance(vtx).value();
	JetInfo[iJetColl].SV_flightErr[JetInfo[iJetColl].nSV]   = svTagInfo->flightDistance(vtx).error();
	JetInfo[iJetColl].SV_flight2D[JetInfo[iJetColl].nSV]    = svTagInfo->flightDistance(vtx, true).value();
	JetInfo[iJetColl].SV_flight2DErr[JetInfo[iJetColl].nSV] = svTagInfo->flightDistance(vtx, true).error();
	JetInfo[iJetColl].SV_nTrk[JetInfo[iJetColl].nSV]        = vtxTracks(svTagInfo->secondaryVertex(vtx));

	const Vertex &vertex = svTagInfo->secondaryVertex(vtx);

	JetInfo[iJetColl].SV_vtx_pt[JetInfo[iJetColl].nSV]  = vertex.p4().pt();
	JetInfo[iJetColl].SV_vtx_eta[JetInfo[iJetColl].nSV] = vertex.p4().eta();
	JetInfo[iJetColl].SV_vtx_phi[JetInfo[iJetColl].nSV] = vertex.p4().phi();
	JetInfo[iJetColl].SV_mass[JetInfo[iJetColl].nSV]    = vertex.p4().mass();

	Int_t totcharge=0;
	reco::TrackKinematics vertexKinematics;

	// get the vertex kinematics and charge
	vertexKinematicsAndChange(vertex, vertexKinematics, totcharge);
        if (currentAxes.size() > 1)
        {
                if (reco::deltaR2(svTagInfo->flightDirection(vtx),currentAxes[1]) < reco::deltaR2(svTagInfo->flightDirection(vtx),currentAxes[0]))
                {
                        tau2Kinematics  = tau2Kinematics + vertexKinematics;
                        JetInfo[iJetColl].Jet_tau2_vertexNTracks[JetInfo[iJetColl].nJet] += vtxTracks(svTagInfo->secondaryVertex(vtx));
                        if( JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet] < 0 ) 
                        {
                          JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet] = svTagInfo->flightDistance(vtx,true).significance();
                          JetInfo[iJetColl].Jet_tau2_vertexDeltaR[JetInfo[iJetColl].nJet] = reco::deltaR(svTagInfo->flightDirection(vtx),currentAxes[1]);
                          tau2_vtx = vtx;
                        }
                        JetInfo[iJetColl].Jet_tau2_nSecondaryVertices[JetInfo[iJetColl].nJet] += 1.;
                }
                else
                {
                        tau1Kinematics = tau1Kinematics + vertexKinematics;
                        JetInfo[iJetColl].Jet_tau1_vertexNTracks[JetInfo[iJetColl].nJet] += vtxTracks(svTagInfo->secondaryVertex(vtx));
                        if( JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] < 0 )
                        {
                          JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] = svTagInfo->flightDistance(vtx,true).significance();
                          JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet] = reco::deltaR(svTagInfo->flightDirection(vtx),currentAxes[0]);
                          tau1_vtx = vtx;
                        }
                        JetInfo[iJetColl].Jet_tau1_nSecondaryVertices[JetInfo[iJetColl].nJet] += 1.;
                }

        }
        else if (currentAxes.size() > 0)
        {
                tau1Kinematics = tau1Kinematics + vertexKinematics;
                JetInfo[iJetColl].Jet_tau1_vertexNTracks[JetInfo[iJetColl].nJet] += vtxTracks(svTagInfo->secondaryVertex(vtx));
                if( JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] < 0 )
                {
                  JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] = svTagInfo->flightDistance(vtx,true).significance();
                  JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet] = reco::deltaR(svTagInfo->flightDirection(vtx),currentAxes[0]);
                  tau1_vtx = vtx;
                }
                JetInfo[iJetColl].Jet_tau1_nSecondaryVertices[JetInfo[iJetColl].nJet] += 1.;
        }

	// total charge at the secondary vertex
	JetInfo[iJetColl].SV_totCharge[JetInfo[iJetColl].nSV]=totcharge;

	math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
	edm::RefToBase<reco::Jet> jet = ipTagInfo->jet();
	math::XYZVector jetDir = jet->momentum().Unit();
	GlobalVector flightDir = svTagInfo->flightDirection(vtx);

	JetInfo[iJetColl].SV_deltaR_jet[JetInfo[iJetColl].nSV]     = ( reco::deltaR(flightDir, jetDir) );
	JetInfo[iJetColl].SV_deltaR_sum_jet[JetInfo[iJetColl].nSV] = ( reco::deltaR(vertexSum, jetDir) );
	JetInfo[iJetColl].SV_deltaR_sum_dir[JetInfo[iJetColl].nSV] = ( reco::deltaR(vertexSum, flightDir) );

	Line::PositionType pos(GlobalPoint(position(vertex).x(),position(vertex).y(),position(vertex).z()));
	Line trackline(pos,flightDir);
	// get the Jet  line
	Line::PositionType pos2(GlobalPoint(pv->x(),pv->y(),pv->z()));
	Line::DirectionType dir2(GlobalVector(jetDir.x(),jetDir.y(),jetDir.z()));
	Line jetline(pos2,dir2);
	// now compute the distance between the two lines
	JetInfo[iJetColl].SV_vtxDistJetAxis[JetInfo[iJetColl].nSV] = (jetline.distance(trackline)).mag();

	JetInfo[iJetColl].SV_EnergyRatio[JetInfo[iJetColl].nSV]= vertexSum.E() / allSum.E();
	JetInfo[iJetColl].SV_dir_x[JetInfo[iJetColl].nSV]= flightDir.x();
	JetInfo[iJetColl].SV_dir_y[JetInfo[iJetColl].nSV]= flightDir.y();
	JetInfo[iJetColl].SV_dir_z[JetInfo[iJetColl].nSV]= flightDir.z();


	++JetInfo[iJetColl].nSV;

    } //// if secondary vertices present
    JetInfo[iJetColl].Jet_nLastSV[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nSV;

    JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet] = -1;
    JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet] = -1.;
    JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet] = -1;
    JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet] = -1.;
    JetInfo[iJetColl].Jet_tau1_vertexMass_corrected[JetInfo[iJetColl].nJet] = -1;
    JetInfo[iJetColl].Jet_tau2_vertexMass_corrected[JetInfo[iJetColl].nJet] = -1.;

    if ( JetInfo[iJetColl].Jet_tau1_nSecondaryVertices[JetInfo[iJetColl].nJet] > 0. )
    {
      math::XYZTLorentzVector tau1_vertexSum = tau1Kinematics.weightedVectorSum();
      JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet] = tau1_vertexSum.E() / allSum.E();
      if ( JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet] > 50. ) JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet] = 50.;

      double tau1_vertexMass = tau1_vertexSum.M();
      JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet] = tau1_vertexMass;

      GlobalVector dir = svTagInfo->flightDirection(tau1_vtx);
      double tau1_vertexPt2 = math::XYZVector(dir.x(), dir.y(), dir.z()).Cross(tau1_vertexSum).Mag2() / dir.mag2();
      tau1_vertexMass = std::sqrt(tau1_vertexMass * tau1_vertexMass + tau1_vertexPt2) + std::sqrt(tau1_vertexPt2);
      JetInfo[iJetColl].Jet_tau1_vertexMass_corrected[JetInfo[iJetColl].nJet] = tau1_vertexMass;
    }

    if ( JetInfo[iJetColl].Jet_tau2_nSecondaryVertices[JetInfo[iJetColl].nJet] > 0. )
    {
      math::XYZTLorentzVector tau2_vertexSum = tau2Kinematics.weightedVectorSum();
      JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet] = tau2_vertexSum.E() / allSum.E();
      if ( JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet] > 50. ) JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet] = 50.;

      double tau2_vertexMass = tau2_vertexSum.M();
      JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet] = tau2_vertexMass;

      GlobalVector dir = svTagInfo->flightDirection(tau2_vtx);
      double tau2_vertexPt2 = math::XYZVector(dir.x(), dir.y(), dir.z()).Cross(tau2_vertexSum).Mag2() / dir.mag2();
      tau2_vertexMass = std::sqrt(tau2_vertexMass * tau2_vertexMass + tau2_vertexPt2) + std::sqrt(tau2_vertexPt2);
      JetInfo[iJetColl].Jet_tau2_vertexMass_corrected[JetInfo[iJetColl].nJet] = tau2_vertexMass;
    }

    if ( runFatJets_ && iJetColl == 0 )
    {
      // when only one tau axis has SVs assigned, they are all assigned to the 1st tau axis
      // in the special case below need to swap values
      if( (JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet]<0 && JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet]>0) )
      {
        float temp = JetInfo[iJetColl].Jet_tau1_nSecondaryVertices[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_nSecondaryVertices[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_nSecondaryVertices[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_nSecondaryVertices[JetInfo[iJetColl].nJet] = temp;
        //-----
        temp = JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet] = temp;
        //-----
        temp = JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_vertexDeltaR[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_vertexDeltaR[JetInfo[iJetColl].nJet] = temp;
        //-----
        temp = JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet] = temp;
        //-----
        temp = JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet] = temp;
        //-----
        temp = JetInfo[iJetColl].Jet_tau1_vertexMass_corrected[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_vertexMass_corrected[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_vertexMass_corrected[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_vertexMass_corrected[JetInfo[iJetColl].nJet] = temp;
        //-----
        temp = JetInfo[iJetColl].Jet_tau1_vertexNTracks[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau1_vertexNTracks[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].Jet_tau2_vertexNTracks[JetInfo[iJetColl].nJet];
        JetInfo[iJetColl].Jet_tau2_vertexNTracks[JetInfo[iJetColl].nJet] = temp;
      }

      // get the TaggingVariables
      const reco::TaggingVariableList vars = bdsvTagInfo->taggingVariables();

      //--------------------------
      JetInfo[iJetColl].Jet_z_ratio[JetInfo[iJetColl].nJet] = vars.get(reco::btau::z_ratio);

      JetInfo[iJetColl].Jet_trackSip3dSig_3[JetInfo[iJetColl].nJet]=vars.get(reco::btau::trackSip3dSig_3);
      JetInfo[iJetColl].Jet_trackSip3dSig_2[JetInfo[iJetColl].nJet]=vars.get(reco::btau::trackSip3dSig_2);
      JetInfo[iJetColl].Jet_trackSip3dSig_1[JetInfo[iJetColl].nJet]=vars.get(reco::btau::trackSip3dSig_1);
      JetInfo[iJetColl].Jet_trackSip3dSig_0[JetInfo[iJetColl].nJet]=vars.get(reco::btau::trackSip3dSig_0);

      JetInfo[iJetColl].Jet_tau2_trackSip3dSig_0[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_trackSip3dSig_0);
      JetInfo[iJetColl].Jet_tau1_trackSip3dSig_0[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_trackSip3dSig_0);
      JetInfo[iJetColl].Jet_tau2_trackSip3dSig_1[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_trackSip3dSig_1);
      JetInfo[iJetColl].Jet_tau1_trackSip3dSig_1[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_trackSip3dSig_1);

      JetInfo[iJetColl].Jet_trackSip2dSigAboveCharm_0[JetInfo[iJetColl].nJet] = vars.get(reco::btau::trackSip2dSigAboveCharm);
      JetInfo[iJetColl].Jet_trackSip2dSigAboveBottom_0[JetInfo[iJetColl].nJet] = vars.get(reco::btau::trackSip2dSigAboveBottom_0);
      JetInfo[iJetColl].Jet_trackSip2dSigAboveBottom_1[JetInfo[iJetColl].nJet] = vars.get(reco::btau::trackSip2dSigAboveBottom_1);    

      JetInfo[iJetColl].Jet_tau2_trackEtaRel_0[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_trackEtaRel_0);
      JetInfo[iJetColl].Jet_tau2_trackEtaRel_1[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_trackEtaRel_1);
      JetInfo[iJetColl].Jet_tau2_trackEtaRel_2[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_trackEtaRel_2);
      JetInfo[iJetColl].Jet_tau1_trackEtaRel_0[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_trackEtaRel_0);
      JetInfo[iJetColl].Jet_tau1_trackEtaRel_1[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_trackEtaRel_1);
      JetInfo[iJetColl].Jet_tau1_trackEtaRel_2[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_trackEtaRel_2);

      JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_vertexMass);
      JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_vertexEnergyRatio);
      JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_vertexDeltaR);
      JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau1_flightDistance2dSig);
      JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_vertexMass);
      JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_vertexEnergyRatio);
      JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet] = vars.get(reco::btau::tau2_flightDistance2dSig);

      JetInfo[iJetColl].Jet_nTracks_fat[JetInfo[iJetColl].nJet] = vars.get(reco::btau::jetNTracks);
      JetInfo[iJetColl].Jet_nSV_fat[JetInfo[iJetColl].nJet] = vars.get(reco::btau::jetNSecondaryVertices);
      //--------------------------

      std::map<std::string,float> variables;
      variables["z_ratio"] = JetInfo[iJetColl].Jet_z_ratio[JetInfo[iJetColl].nJet];
      variables["trackSipdSig_3"] = JetInfo[iJetColl].Jet_trackSip3dSig_3[JetInfo[iJetColl].nJet]; 
      variables["trackSipdSig_2"] = JetInfo[iJetColl].Jet_trackSip3dSig_2[JetInfo[iJetColl].nJet]; 
      variables["trackSipdSig_1"] = JetInfo[iJetColl].Jet_trackSip3dSig_1[JetInfo[iJetColl].nJet]; 
      variables["trackSipdSig_0"] = JetInfo[iJetColl].Jet_trackSip3dSig_0[JetInfo[iJetColl].nJet];
      variables["trackSipdSig_1_0"] = JetInfo[iJetColl].Jet_tau2_trackSip3dSig_0[JetInfo[iJetColl].nJet];
      variables["trackSipdSig_0_0"] = JetInfo[iJetColl].Jet_tau1_trackSip3dSig_0[JetInfo[iJetColl].nJet];
      variables["trackSipdSig_1_1"] = JetInfo[iJetColl].Jet_tau2_trackSip3dSig_1[JetInfo[iJetColl].nJet];
      variables["trackSipdSig_0_1"] = JetInfo[iJetColl].Jet_tau1_trackSip3dSig_1[JetInfo[iJetColl].nJet]; 
      variables["trackSip2dSigAboveCharm_0"] = JetInfo[iJetColl].Jet_trackSip2dSigAboveCharm_0[JetInfo[iJetColl].nJet];
      variables["trackSip2dSigAboveBottom_0"] = JetInfo[iJetColl].Jet_trackSip2dSigAboveBottom_0[JetInfo[iJetColl].nJet];
      variables["trackSip2dSigAboveBottom_1"] = JetInfo[iJetColl].Jet_trackSip2dSigAboveBottom_1[JetInfo[iJetColl].nJet];
      variables["tau1_trackEtaRel_0"] = JetInfo[iJetColl].Jet_tau2_trackEtaRel_0[JetInfo[iJetColl].nJet];
      variables["tau1_trackEtaRel_1"] = JetInfo[iJetColl].Jet_tau2_trackEtaRel_1[JetInfo[iJetColl].nJet];
      variables["tau1_trackEtaRel_2"] = JetInfo[iJetColl].Jet_tau2_trackEtaRel_2[JetInfo[iJetColl].nJet];
      variables["tau0_trackEtaRel_0"] = JetInfo[iJetColl].Jet_tau1_trackEtaRel_0[JetInfo[iJetColl].nJet];
      variables["tau0_trackEtaRel_1"] = JetInfo[iJetColl].Jet_tau1_trackEtaRel_1[JetInfo[iJetColl].nJet];
      variables["tau0_trackEtaRel_2"] = JetInfo[iJetColl].Jet_tau1_trackEtaRel_2[JetInfo[iJetColl].nJet];
      variables["tau_vertexMass_0"] = JetInfo[iJetColl].Jet_tau1_vertexMass[JetInfo[iJetColl].nJet];
      variables["tau_vertexEnergyRatio_0"] = JetInfo[iJetColl].Jet_tau1_vertexEnergyRatio[JetInfo[iJetColl].nJet];
      variables["tau_vertexDeltaR_0"] = JetInfo[iJetColl].Jet_tau1_vertexDeltaR[JetInfo[iJetColl].nJet];
      variables["tau_flightDistance2dSig_0"] = JetInfo[iJetColl].Jet_tau1_flightDistance2dSig[JetInfo[iJetColl].nJet];
      variables["tau_vertexMass_1"] = JetInfo[iJetColl].Jet_tau2_vertexMass[JetInfo[iJetColl].nJet];
      variables["tau_vertexEnergyRatio_1"] = JetInfo[iJetColl].Jet_tau2_vertexEnergyRatio[JetInfo[iJetColl].nJet];
      variables["tau_flightDistance2dSig_1"] = JetInfo[iJetColl].Jet_tau2_flightDistance2dSig[JetInfo[iJetColl].nJet];
      variables["jetNTracks"] = JetInfo[iJetColl].Jet_nTracks_fat[JetInfo[iJetColl].nJet];
      variables["nSV"] = JetInfo[iJetColl].Jet_nSV_fat[JetInfo[iJetColl].nJet];

      JetInfo[iJetColl].Jet_BDTG_SV[JetInfo[iJetColl].nJet] = evaluator_SV_->evaluate(variables);

//       std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
//       std::cout << "Jet pt, eta, phi: " << pjet->pt() << ", " << pjet->eta() << ", " << pjet->phi() << std::endl;
//       for(std::map<std::string,float>::const_iterator it = variables.begin(); it != variables.end(); ++it)
//       {
//         std::cout << it->first << ": " << it->second << std::endl;
//       }
//       std::cout << "Discriminator: " << JetInfo[iJetColl].Jet_BDTG_SV[JetInfo[iJetColl].nJet] << std::endl;
    }


    cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0;
    can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0;

    //*****************************************************************
    // for Mistag studies
    //*****************************************************************
    const Tracks svxPostracks( svTagInfo->vertexTracks(0) );

    if ( useTrackHistory_ && !isData_ ) {
      for (unsigned int i=0; i<svxPostracks.size(); ++i) {
        TrackCategories::Flags theFlag = classifier_.evaluate( toTrackRef(svxPostracks[i]) ).flags();
        if ( theFlag[TrackCategories::BWeakDecay] )         cap0 = 1;
        if ( theFlag[TrackCategories::CWeakDecay] )         cap1 = 1;
        if ( theFlag[TrackCategories::SignalEvent] )           cap2 = 1;
        if ( theFlag[TrackCategories::ConversionsProcess] ) cap3 = 1;
        if ( theFlag[TrackCategories::KsDecay] )            cap4 = 1;
        if ( theFlag[TrackCategories::LambdaDecay] )        cap5 = 1;
        if ( theFlag[TrackCategories::HadronicProcess] )    cap6 = 1;
        if ( theFlag[TrackCategories::Fake] )	            cap7 = 1;
        if ( theFlag[TrackCategories::SharedInnerHits] )    cap8 = 1;
      }
    }

    const Tracks svxNegtracks( svNegTagInfo->vertexTracks(0) );

    if ( useTrackHistory_ && !isData_ ) {
      for (unsigned int i=0; i<svxNegtracks.size(); ++i) {
        TrackCategories::Flags theFlag = classifier_.evaluate( toTrackRef(svxNegtracks[i]) ).flags();
        if ( theFlag[TrackCategories::BWeakDecay] )         can0 = 2;
        if ( theFlag[TrackCategories::CWeakDecay] )         can1 = 2;
        if ( theFlag[TrackCategories::SignalEvent] )           can2 = 2;
        if ( theFlag[TrackCategories::ConversionsProcess] ) can3 = 2;
        if ( theFlag[TrackCategories::KsDecay] )            can4 = 2;
        if ( theFlag[TrackCategories::LambdaDecay] )        can5 = 2;
        if ( theFlag[TrackCategories::HadronicProcess] )    can6 = 2;
        if ( theFlag[TrackCategories::Fake] )	            can7 = 2;
        if ( theFlag[TrackCategories::SharedInnerHits] )    can8 = 2;
      }
      JetInfo[iJetColl].Jet_histSvx[JetInfo[iJetColl].nJet] =  cap0+can0 + (cap1+can1)*10 + (cap2+can2)*100
        + (cap3+can3)*1000     + (cap4+can4)*10000
        + (cap5+can5)*100000   + (cap6+can6)*1000000
        + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }

    //*********************************
    // Jet analysis

    Histos[iJetColl]->hData_All_NTracks->Fill( ntagtracks ) ;
    Histos[iJetColl]->hData_All_JetPt->Fill( ptjet );
    Histos[iJetColl]->hData_All_JetEta->Fill( TMath::Abs( etajet ) );

    //*****************************************************************
    //define positive and negative tags
    //*****************************************************************
    float varpos = -1000.;
    float varneg = -1000.;

    if ( selTagger_ == 0 ) {       // jet proba
      if ( ProbaP > 0 ) varpos = 20.*ProbaP;
      if ( ProbaN > 0 ) varneg = 20.*ProbaN;
    }
    else if ( selTagger_ == 2 ) {     // TC High Eff.
      if ( JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet] > 0. ) varpos = JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet];
      if ( JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet] > 0. ) varneg = JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 3 ) {     // TC High Pure.
      if ( JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet] > 0. ) varpos = JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet];
      if ( JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet] > 0. ) varneg = JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 4 ) {    // SSVHE
      if  (JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet] > 1)   varpos =  10.*JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet] - 10.;
      if  (JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] < -1) varneg = -10.*JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] - 10.;
    }
    else if ( selTagger_ == 5 ) {    // SV combined
      if  (JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet] > 0)  varpos = 50.*JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet];
      if  (JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] > 0) varneg = 50.*JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 6 ) {    // soft muon ptrel
      if  ( JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet]  > 0) varpos =  5*JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet];
      if  ( JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet] < 0) varneg = -5*JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 7 ) {    // SSVHP
      if  (JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet] > 1)   varpos =  10.*JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet] - 10.;
      if  (JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] < -1) varneg = -10.*JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] - 10.;
    }

    Histos[iJetColl]->hData_NTracks->Fill( ntagtracks ) ;
    Histos[iJetColl]->hData_JetPt->Fill( ptjet );
    Histos[iJetColl]->hData_JetEta->Fill( TMath::Abs( etajet ) );

    //*********************************
    // Tagging

    if ( varneg > 0 ) Histos[iJetColl]->hData_Tagger->Fill(-varneg );
    if ( varpos > 0 ) Histos[iJetColl]->hData_Tagger->Fill( varpos );

    if ( JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHE->Fill( JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHE->Fill(-JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHP->Fill( JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHP->Fill(-JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_ProbaP[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->  hData_Tagger_JP->Fill( 10*JetInfo[iJetColl].Jet_ProbaP[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_ProbaN[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->  hData_Tagger_JP->Fill(-10*JetInfo[iJetColl].Jet_ProbaN[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet]  >  1 ) Histos[iJetColl]->   hData_Tagger_SSVHE->Fill( 5*JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] < -1 ) Histos[iJetColl]->   hData_Tagger_SSVHE->Fill( 5*JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet]  >  1 ) Histos[iJetColl]-> hData_Tagger_SSVHP->Fill( 5*JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] < -1 ) Histos[iJetColl]-> hData_Tagger_SSVHP->Fill( 5*JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet] > 0  ) Histos[iJetColl]->hData_Tagger_CSV->Fill( 25*JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->hData_Tagger_CSV->Fill(-25*JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet]  > 0 ) Histos[iJetColl]-> hData_Tagger_MU->Fill( 2.5*JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet] < 0 ) Histos[iJetColl]-> hData_Tagger_MU->Fill( 2.5*JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet] );


    if ( !isData_ )
    {
      Histos[iJetColl]->hAllFlav_Flavour->Fill( flavour );
      if ( varneg > 0 ) Histos[iJetColl]->hAllFlav_Tagger->Fill(-varneg );
      if ( varpos > 0 ) Histos[iJetColl]->hAllFlav_Tagger->Fill( varpos );
      if ( flavour == 1 || flavour == 21 ) { // light q+gluon
        if ( varneg > 0 ) Histos[iJetColl]->hLightFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hLightFlav_Tagger->Fill( varpos );
      } // light q+gluon
      if ( flavour == 21 ) { // gluon jets
        if ( varneg > 0 ) Histos[iJetColl]->hGluonFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hGluonFlav_Tagger->Fill( varpos );
      } // gluon jets
      else if ( flavour == 1 ) { // uds jets
        if ( varneg > 0 ) Histos[iJetColl]->hUDSFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hUDSFlav_Tagger->Fill( varpos );
      } // uds jets
      else if ( flavour == 4 ) { // c jets
        if ( varneg > 0 ) Histos[iJetColl]->hCFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hCFlav_Tagger->Fill( varpos );
      } // c jets
      else if ( flavour == 5 ) { // b jets
        if ( varneg > 0 ) Histos[iJetColl]->hBFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hBFlav_Tagger->Fill( varpos );
      } // b jets
    }



    // Track History
    if ( useTrackHistory_ && indexes.size()!=0 && !isData_ ) {

      int cat1P = 0, cat2P = 0, cat3P = 0, catP = 0;
      int cat1N = 0, cat2N = 0, cat3N = 0, catN = 0;
      flags1P[TrackCategories::ConversionsProcess] ;
      if ( flags1P[TrackCategories::BWeakDecay] )              cat1P = 1;
      else if ( flags1P[TrackCategories::CWeakDecay] )         cat1P = 2;
      else if ( flags1P[TrackCategories::SignalEvent] )           cat1P = 3;
      else if ( flags1P[TrackCategories::ConversionsProcess] ) cat1P = 4;
      else if ( flags1P[TrackCategories::KsDecay] )            cat1P = 5;
      else if ( flags1P[TrackCategories::LambdaDecay] )        cat1P = 6;
      else if ( flags1P[TrackCategories::HadronicProcess] )    cat1P = 7;
      else if ( flags1P[TrackCategories::Fake] )               cat1P = 8;
      else if ( flags1P[TrackCategories::SharedInnerHits] )    cat1P = 9;
      if(idSize > 1){
        if ( flags2P[TrackCategories::BWeakDecay] )              cat2P = 1;
        else if ( flags2P[TrackCategories::CWeakDecay] )         cat2P = 2;
        else if ( flags2P[TrackCategories::SignalEvent] )           cat2P = 3;
        else if ( flags2P[TrackCategories::ConversionsProcess] ) cat2P = 4;
        else if ( flags2P[TrackCategories::KsDecay] )            cat2P = 5;
        else if ( flags2P[TrackCategories::LambdaDecay] )        cat2P = 6;
        else if ( flags2P[TrackCategories::HadronicProcess] )    cat2P = 7;
        else if ( flags2P[TrackCategories::Fake] )               cat2P = 8;
        else if ( flags2P[TrackCategories::SharedInnerHits] )    cat2P = 9;
      }
      if(idSize > 2){
        if ( flags3P[TrackCategories::BWeakDecay] )              cat3P = 1;
        else if ( flags3P[TrackCategories::CWeakDecay] )         cat3P = 2;
        else if ( flags3P[TrackCategories::SignalEvent] )           cat3P = 3;
        else if ( flags3P[TrackCategories::ConversionsProcess] ) cat3P = 4;
        else if ( flags3P[TrackCategories::KsDecay] )            cat3P = 5;
        else if ( flags3P[TrackCategories::LambdaDecay] )        cat3P = 6;
        else if ( flags3P[TrackCategories::HadronicProcess] )    cat3P = 7;
        else if ( flags3P[TrackCategories::Fake] )               cat3P = 8;
        else if ( flags3P[TrackCategories::SharedInnerHits] )    cat3P = 9;
      }
      if ( flags1N[TrackCategories::BWeakDecay] )              cat1N = 1;
      else if ( flags1N[TrackCategories::CWeakDecay] )         cat1N = 2;
      else if ( flags1N[TrackCategories::SignalEvent] )           cat1N = 3;
      else if ( flags1N[TrackCategories::ConversionsProcess] ) cat1N = 4;
      else if ( flags1N[TrackCategories::KsDecay] )            cat1N = 5;
      else if ( flags1N[TrackCategories::LambdaDecay] )        cat1N = 6;
      else if ( flags1N[TrackCategories::HadronicProcess] )    cat1N = 7;
      else if ( flags1N[TrackCategories::Fake] )               cat1N = 8;
      else if ( flags1N[TrackCategories::SharedInnerHits] )    cat1N = 9;
      if(idSize > 1){
        if ( flags2N[TrackCategories::BWeakDecay] )              cat2N = 1;
        else if ( flags2N[TrackCategories::CWeakDecay] )         cat2N = 2;
        else if ( flags2N[TrackCategories::SignalEvent] )           cat2N = 3;
        else if ( flags2N[TrackCategories::ConversionsProcess] ) cat2N = 4;
        else if ( flags2N[TrackCategories::KsDecay] )            cat2N = 5;
        else if ( flags2N[TrackCategories::LambdaDecay] )        cat2N = 6;
        else if ( flags2N[TrackCategories::HadronicProcess] )    cat2N = 7;
        else if ( flags2N[TrackCategories::Fake] )               cat2N = 8;
        else if ( flags2N[TrackCategories::SharedInnerHits] )    cat2N = 9;
      }
      if(idSize > 2){
        if ( flags3N[TrackCategories::BWeakDecay] )              cat3N = 1;
        else if ( flags3N[TrackCategories::CWeakDecay] )         cat3N = 2;
        else if ( flags3N[TrackCategories::SignalEvent] )           cat3N = 3;
        else if ( flags3N[TrackCategories::ConversionsProcess] ) cat3N = 4;
        else if ( flags3N[TrackCategories::KsDecay] )            cat3N = 5;
        else if ( flags3N[TrackCategories::LambdaDecay] )        cat3N = 6;
        else if ( flags3N[TrackCategories::HadronicProcess] )    cat3N = 7;
        else if ( flags3N[TrackCategories::Fake] )               cat3N = 8;
        else if ( flags3N[TrackCategories::SharedInnerHits] )    cat3N = 9;
      }

      if ( selTagger_ == 0 || selTagger_ == 3 ) {
        catP = cat3P;
        catN = cat3N;
      }
      else if ( selTagger_ == 2 ) {
        catP = cat2P;
        catN = cat2N;
      }
      else if ( selTagger_ == 1 ) {
        catP = cat1P;
        catN = cat1N;
      }

      if ( varneg > 0 ) {
        if      ( catN == 1 ) Histos[iJetColl]->hAllFlav_Tagger_Bwd->Fill(-varneg );
        else if ( catN == 2 ) Histos[iJetColl]->hAllFlav_Tagger_Cwd->Fill(-varneg );
        else if ( catN == 3 ) Histos[iJetColl]->hAllFlav_Tagger_Tau->Fill(-varneg );
        else if ( catN == 4 ) Histos[iJetColl]->hAllFlav_Tagger_Gam->Fill(-varneg );
        else if ( catN == 5 ) Histos[iJetColl]->hAllFlav_Tagger_K0s->Fill(-varneg );
        else if ( catN == 6 ) Histos[iJetColl]->hAllFlav_Tagger_Lam->Fill(-varneg );
        else if ( catN == 7 ) Histos[iJetColl]->hAllFlav_Tagger_Int->Fill(-varneg );
        else if ( catN == 8 ) Histos[iJetColl]->hAllFlav_Tagger_Fak->Fill(-varneg );
        else if ( catN == 9 ) Histos[iJetColl]->hAllFlav_Tagger_Bad->Fill(-varneg );
        else                  Histos[iJetColl]->hAllFlav_Tagger_Oth->Fill(-varneg );
      }
      if ( varpos > 0 ) {
        if      ( catP == 1 ) Histos[iJetColl]->hAllFlav_Tagger_Bwd->Fill( varpos );
        else if ( catP == 2 ) Histos[iJetColl]->hAllFlav_Tagger_Cwd->Fill( varpos );
        else if ( catP == 3 ) Histos[iJetColl]->hAllFlav_Tagger_Tau->Fill( varpos );
        else if ( catP == 4 ) Histos[iJetColl]->hAllFlav_Tagger_Gam->Fill( varpos );
        else if ( catP == 5 ) Histos[iJetColl]->hAllFlav_Tagger_K0s->Fill( varpos );
        else if ( catP == 6 ) Histos[iJetColl]->hAllFlav_Tagger_Lam->Fill( varpos );
        else if ( catP == 7 ) Histos[iJetColl]->hAllFlav_Tagger_Int->Fill( varpos );
        else if ( catP == 8 ) Histos[iJetColl]->hAllFlav_Tagger_Fak->Fill( varpos );
        else if ( catP == 9 ) Histos[iJetColl]->hAllFlav_Tagger_Bad->Fill( varpos );
        else                  Histos[iJetColl]->hAllFlav_Tagger_Oth->Fill( varpos );
      }
    }

    ++numjet;
    ++JetInfo[iJetColl].nJet;
  } // end loop on jet

  Histos[iJetColl]->hData_All_NJets->Fill( JetInfo[iJetColl].nJet );
  Histos[iJetColl]->hData_NJets->Fill( numjet );

  return;
} // BTagAnalyzerT:: processJets

template<typename IPTI,typename VTX>
float BTagAnalyzerT<IPTI,VTX>::calculPtRel(const reco::Track& theMuon, const pat::Jet& theJet )
{
  double pmu = TMath::Sqrt( theMuon.px()*theMuon.px() + theMuon.py()*theMuon.py()  + theMuon.pz()*theMuon.pz() );

  double jetpx = 0 ;
  double jetpy = 0 ;
  double jetpz = 0 ;

  if( theJet.isCaloJet() ){
    jetpx = theJet.px() + theMuon.px();
    jetpy = theJet.py() + theMuon.py();
    jetpz = theJet.pz() + theMuon.pz();
  }else{
    jetpx = theJet.px();
    jetpy = theJet.py();
    jetpz = theJet.pz();
  }

  double jetp = TMath::Sqrt(jetpx*jetpx + jetpy*jetpy + jetpz*jetpz);

  double ptrel  = ( jetpx * theMuon.px()  + jetpy * theMuon.py() + jetpz * theMuon.pz() ) / jetp;
  ptrel = TMath::Sqrt( pmu * pmu  - ptrel * ptrel );

  return ptrel;
}

template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::setTracksPVBase( const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight )
{
  iPV = -1;
  PVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::VertexCollection::const_iterator IV;
  typedef reco::Vertex::trackRef_iterator IT;

  for(IV iv=pvHandle->begin(); iv!=pvHandle->end(); ++iv)
  {
    const reco::Vertex & vtx = *iv;
    // loop over tracks in vertices
    for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
    {
      const reco::TrackBaseRef & baseRef = *it;
      // one of the tracks in the vertex is the same as the track considered in the function
      if( baseRef == trackBaseRef )
      {
        float w = vtx.trackWeight(baseRef);
        // select the vertex for which the track has the highest weight
        if( w > PVweight )
        {
          PVweight = w;
          iPV = ( iv - pvHandle->begin() );
          break;
        }
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::beginJob() {
  //cat 0
  cat0.etaMax = 2.5;
  cat0.etaMin = 0;
  cat0.nHitsMax= 50;
  cat0.nHitsMin= 1;
  cat0.nPixelHitsMax = 1;
  cat0.nPixelHitsMin = 1;
  cat0.pMax= 5000;
  cat0.pMin= 0;
  cat0.chiMin= 0;
  cat0.chiMax= 5;
  cat0.withFirstPixel = 0;

  //cat 1
  cat1.etaMax        = 2.5;
  cat1.etaMin        = 0;
  cat1.nHitsMax      = 50;
  cat1.nHitsMin      = 1;
  cat1.nPixelHitsMax = 8;
  cat1.nPixelHitsMin = 2;
  cat1.pMax          = 5000;
  cat1.pMin          = 0;
  cat1.chiMin        = 2.5;
  cat1.chiMax        = 5;
  cat1.withFirstPixel = 0;

  //cat 2
  cat2.etaMax        = 0.8;
  cat2.etaMin        = 0;
  cat2.nHitsMax      = 50;
  cat2.nHitsMin      = 1;
  cat2.nPixelHitsMax = 8;
  cat2.nPixelHitsMin = 3;
  cat2.pMax          = 8;
  cat2.pMin          = 0;
  cat2.chiMin        = 0;
  cat2.chiMax        = 2.5;
  cat2.withFirstPixel = 0;

  //cat 3
  cat3.etaMax        = 1.6;
  cat3.etaMin        = 0.8;
  cat3.nHitsMax      = 50;
  cat3.nHitsMin      = 1;
  cat3.nPixelHitsMax = 8;
  cat3.nPixelHitsMin = 3;
  cat3.pMax          = 8;
  cat3.pMin          = 0;
  cat3.chiMin        = 0;
  cat3.chiMax        = 2.5;
  cat3.withFirstPixel = 0;

  //cat 4
  cat4.etaMax        = 2.5;
  cat4.etaMin        = 1.6;
  cat4.nHitsMax      = 50;
  cat4.nHitsMin      = 1;
  cat4.nPixelHitsMax = 8;
  cat4.nPixelHitsMin = 3;
  cat4.pMax          = 8;
  cat4.pMin          = 0;
  cat4.chiMin        = 0;
  cat4.chiMax        = 2.5;
  cat4.withFirstPixel = 0;

  //cat 5
  cat5.etaMax        = 2.5;
  cat5.etaMin        = 0;
  cat5.nHitsMax      = 50;
  cat5.nHitsMin      = 1;
  cat5.nPixelHitsMax = 2;
  cat5.nPixelHitsMin = 2;
  cat5.pMax          = 8;
  cat5.pMin          = 0;
  cat5.chiMin        = 0;
  cat5.chiMax        = 2.5;
  cat5.withFirstPixel = 0;

  //cat 6
  cat6.etaMax        = 0.8;
  cat6.etaMin        = 0;
  cat6.nHitsMax      = 50;
  cat6.nHitsMin      = 1;
  cat6.nPixelHitsMax = 8;
  cat6.nPixelHitsMin = 3;
  cat6.pMax          = 5000;
  cat6.pMin          = 8;
  cat6.chiMin        = 0;
  cat6.chiMax        = 2.5;
  cat6.withFirstPixel = 0;

  //cat 7
  cat7.etaMax        = 1.6;
  cat7.etaMin        = 0.8;
  cat7.nHitsMax      = 50;
  cat7.nHitsMin      = 1;
  cat7.nPixelHitsMax = 8;
  cat7.nPixelHitsMin = 3;
  cat7.pMax          = 5000;
  cat7.pMin          = 8;
  cat7.chiMin        = 0;
  cat7.chiMax        = 2.5;
  cat7.withFirstPixel = 0;

  //cat 8
  cat8.etaMax        = 2.5;
  cat8.etaMin        = 1.6;
  cat8.nHitsMax      = 50;
  cat8.nHitsMin      = 1;
  cat8.nPixelHitsMax = 8;
  cat8.nPixelHitsMin = 3;
  cat8.pMax          = 5000;
  cat8.pMin          = 8;
  cat8.chiMin        = 0;
  cat8.chiMax        = 2.5;
  cat8.withFirstPixel = 0;

  //cat 9
  cat9.etaMax        = 2.5;
  cat9.etaMin        = 0;
  cat9.nHitsMax      = 50;
  cat9.nHitsMin      = 1;
  cat9.nPixelHitsMax = 2;
  cat9.nPixelHitsMin = 2;
  cat9.pMax          = 5000;
  cat9.pMin          = 8;
  cat9.chiMin        = 0;
  cat9.chiMax        = 2.5;
  cat9.withFirstPixel = 0;
}


// ------------ method called once each job just after ending the event loop  ------------
template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::endJob() {
}

template<typename IPTI,typename VTX>
std::vector< float > BTagAnalyzerT<IPTI,VTX>::getTrackProbabilies(std::vector< float > v, const int ipType){

  std::vector< float > vectTrackProba;

  for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); ++q) {
    //positives and negatives tracks
    double p3d = 0;
    if(ipType == 0){
      if ( *q>=0){p3d= (*q)/2.;}else{p3d=1.+(*q)/2.;}
      //if(-log(p3d)> 5) p3d=exp(-5.0);
      vectTrackProba.push_back(p3d);
    }

    //positives tracks only
    if(ipType == 1 && *q >=0 ){
      vectTrackProba.push_back(*q);
    }

    //negatives tracks only
    if(ipType == 2 && *q <0){
      vectTrackProba.push_back(-(*q));
    }
  }

  return vectTrackProba;
}

template<typename IPTI,typename VTX>
double BTagAnalyzerT<IPTI,VTX>::calculProbability(std::vector< float > v) {

  int ngoodtracks=v.size();
  double SumJet=0.;
  double m_minTrackProb = 0.005;
  for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); ++q) {
    SumJet += (*q>m_minTrackProb)?log(*q):log(m_minTrackProb);
  }

  double ProbJet;
  double Loginvlog = 0;

  if ( SumJet < 0. ) {
    if ( ngoodtracks >= 2 ) {
      Loginvlog = log(-SumJet);
    }
    double Prob = 1.;
    double lfact = 1.;
    for (int l=1; l!=ngoodtracks; ++l) {
      lfact *= l;
      Prob += exp( l*Loginvlog-log(1.*lfact) );
    }
    double LogProb = log(Prob);
    ProbJet = std::min(exp( std::max(LogProb+SumJet,-30.) ), 1.);
  } else {
    ProbJet = 1.;
  }
  return ProbJet;
}

template<typename IPTI,typename VTX>
bool BTagAnalyzerT<IPTI,VTX>::findCat(const reco::Track* track, CategoryFinder& d) {
  //int numcat=-1;

  double   p = track->p();
  double eta = track->eta();
  double chi = track->normalizedChi2();
  int   nhit = track->numberOfValidHits();
  int   npix = track->hitPattern().numberOfValidPixelHits();

  bool result = ( p > d.pMin  &&  p  < d.pMax		  &&
      fabs(eta) > d.etaMin    &&  fabs(eta) < d.etaMax    &&
      nhit >= d.nHitsMin      &&  nhit <= d.nHitsMax	  &&
      npix >= d.nPixelHitsMin &&  npix <= d.nPixelHitsMax &&
      chi >= d.chiMin         &&  chi <= d.chiMax );

  return result;
}

template<typename IPTI,typename VTX>
const edm::Ptr<reco::Muon> BTagAnalyzerT<IPTI,VTX>::matchMuon(const edm::Ptr<reco::Candidate>& theMuon, const edm::View<reco::Muon>& muons ){

  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(theMuon.get());

  if(pcand) // MiniAOD case
  {
    for(edm::View<reco::Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon )
    {
       const pat::Muon * patmuon = dynamic_cast<const pat::Muon *>(&(*muon));

       if(patmuon)
       {
         if(patmuon->originalObjectRef()==theMuon)
           return muons.ptrAt(muon - muons.begin());
       }
    }
    return edm::Ptr<reco::Muon>();
  }
  else
  {
    const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(theMuon.get());

    return edm::refToPtr( pfcand->muonRef() );
  }
}

template<typename IPTI,typename VTX>
std::vector<simPrimaryVertex> BTagAnalyzerT<IPTI,VTX>::getSimPVs(const edm::Handle<edm::HepMCProduct>& evtMC){

  std::vector<simPrimaryVertex> simpv;
  const HepMC::GenEvent* evt=evtMC->GetEvent();
  if (evt) {

    for ( HepMC::GenEvent::vertex_const_iterator vitr= evt->vertices_begin();
        vitr != evt->vertices_end(); ++vitr ) { // loop for vertex ...
      HepMC::FourVector pos = (*vitr)->position();
      //HepLorentzVector pos = (*vitr)->position();

      bool hasMotherVertex=false;

      for ( HepMC::GenVertex::particle_iterator
          mother  = (*vitr)->particles_begin(HepMC::parents);
          mother != (*vitr)->particles_end(HepMC::parents);
          ++mother ) {
        HepMC::GenVertex * mv=(*mother)->production_vertex();
        if (mv) {
          hasMotherVertex=true;
          break; //if verbose_, print all particles of gen vertices
        }
      }

      if (hasMotherVertex) {continue;}

      // could be a new vertex, check  all primaries found so far to avoid multiple entries
      const double mm=0.1;
      simPrimaryVertex sv(pos.x()*mm,pos.y()*mm,pos.z()*mm);
      simPrimaryVertex *vp=NULL;  // will become non-NULL if a vertex is found and then point to it
      for (std::vector<simPrimaryVertex>::iterator v0=simpv.begin();
          v0!=simpv.end(); ++v0) {
        if ( (fabs(sv.x-v0->x)<1e-5) && (fabs(sv.y-v0->y)<1e-5)
            && (fabs(sv.z-v0->z)<1e-5) ) {
          vp=&(*v0);
          break;
        }
      }

      if (!vp) {	  // this is a new vertex
        simpv.push_back(sv);
        vp=&simpv.back();
      }
      vp->genVertex.push_back((*vitr)->barcode());
    }
  }
  return simpv;
}

template<typename IPTI,typename VTX>
int BTagAnalyzerT<IPTI,VTX>::isFromGSP(const reco::Candidate* c)
{
  int isFromGSP = 0;

  if( c->numberOfMothers() == 1 ) {
    const reco::Candidate* dau = c;
    const reco::Candidate* mom = c->mother();
    while( dau->numberOfMothers() == 1 && !( isHardProcess(mom->status()) && (abs(mom->pdgId())==4 || abs(mom->pdgId())==5) ) ) {
      if( abs(mom->pdgId())==21 )
      {
        isFromGSP = 1;
        break;
      }
      dau = mom;
      mom = dau->mother();
    }
  }

  return isFromGSP;
}

template<typename IPTI,typename VTX>
bool BTagAnalyzerT<IPTI,VTX>::isHardProcess(const int status)
{
  // if Pythia8
  if( hadronizerType_ & (1 << 1) )
  {
    if( status>=21 && status<=29 )
      return true;
  }
  else // assuming Pythia6
  {
    if( status==3 )
      return true;
  }

  return false;
}

// -------------------------------------------------------------------------
// NameCompatible
// -------------------------------------------------------------------------
template<typename IPTI,typename VTX>
bool BTagAnalyzerT<IPTI,VTX>::NameCompatible(const std::string& pattern, const std::string& name)
{
  const boost::regex regexp(edm::glob2reg(pattern));

  return boost::regex_match(name,regexp);
}

// ------------ method that matches groomed and original jets based on minimum dR ------------
template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::matchGroomedJets(const edm::Handle<PatJetCollection>& jets,
 const edm::Handle<PatJetCollection>& groomedJets,
 std::vector<int>& matchedIndices)
{
   std::vector<bool> jetLocks(jets->size(),false);
   std::vector<int>  jetIndices;

   for(size_t gj=0; gj<groomedJets->size(); ++gj)
   {
     double matchedDR = 1e9;
     int matchedIdx = -1;

     for(size_t j=0; j<jets->size(); ++j)
     {
       if( jetLocks.at(j) ) continue; // skip jets that have already been matched

       double tempDR = reco::deltaR( jets->at(j).rapidity(), jets->at(j).phi(), groomedJets->at(gj).rapidity(), groomedJets->at(gj).phi() );
       if( tempDR < matchedDR )
       {
         matchedDR = tempDR;
         matchedIdx = j;
       }
     }

     if( matchedIdx>=0 ) jetLocks.at(matchedIdx) = true;
     jetIndices.push_back(matchedIdx);
   }

   if( std::find( jetIndices.begin(), jetIndices.end(), -1 ) != jetIndices.end() )
     edm::LogError("JetMatchingFailed") << "Matching groomed to original jets failed. Please check that the two jet collections belong to each other.";

   for(size_t j=0; j<jets->size(); ++j)
   {
     std::vector<int>::iterator matchedIndex = std::find( jetIndices.begin(), jetIndices.end(), j );

     matchedIndices.push_back( matchedIndex != jetIndices.end() ? std::distance(jetIndices.begin(),matchedIndex) : -1 );
   }
}

// -------------- etaRelToTauAxis ----------------
template<>
void BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::etaRelToTauAxis(const Vertex & vertex, const fastjet::PseudoJet & tauAxis, std::vector<float> & tau_trackEtaRel)
{
  math::XYZVector direction(tauAxis.px(), tauAxis.py(), tauAxis.pz());
  bool hasRefittedTracks = vertex.hasRefittedTracks();
  for(reco::Vertex::trackRef_iterator track = vertex.tracks_begin();
                                      track != vertex.tracks_end(); ++track)
  {
    double w = vertex.trackWeight(*track);
    if (w < 0.5)
      continue;
    if (hasRefittedTracks)
    {
      const reco::Track actualTrack = vertex.refittedTrack(*track);
      tau_trackEtaRel.push_back(std::abs(reco::btau::etaRel(direction.Unit(), actualTrack.momentum())));
    }
    else
      tau_trackEtaRel.push_back(std::abs(reco::btau::etaRel(direction.Unit(), (*track)->momentum())));
  }
}

template<>
void BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::etaRelToTauAxis(const Vertex & vertex, const fastjet::PseudoJet & tauAxis, std::vector<float> & tau_trackEtaRel)
{
  math::XYZVector direction(tauAxis.px(), tauAxis.py(), tauAxis.pz());
  const std::vector<reco::CandidatePtr> & tracks = vertex.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track)
    tau_trackEtaRel.push_back(std::abs(reco::btau::etaRel(direction.Unit(), (*track)->momentum())));
}



// -------------- recalcNsubjettiness ----------------
template<typename IPTI,typename VTX>
void BTagAnalyzerT<IPTI,VTX>::recalcNsubjettiness(const pat::Jet & jet, float & tau1, float & tau2, std::vector<fastjet::PseudoJet> & currentAxes) const
{
  std::vector<fastjet::PseudoJet> fjParticles;

  // loop over jet constituents and push them in the vector of FastJet constituents
  for(const reco::CandidatePtr & daughter : jet.daughterPtrVector())
  {
    if ( daughter.isNonnull() && daughter.isAvailable() )
    {
      const reco::Jet * subjet = dynamic_cast<const reco::Jet *>(daughter.get());
      // if the daughter is actually a subjet
      if( subjet && daughter->numberOfDaughters() > 1 )
      {
        // loop over subjet constituents and push them in the vector of FastJet constituents
        for(size_t i=0; i<daughter->numberOfDaughters(); ++i)
        {
          const reco::Candidate * constit = daughter->daughter(i);

          if ( constit )
            fjParticles.push_back( fastjet::PseudoJet( constit->px(), constit->py(), constit->pz(), constit->energy() ) );
          else
            edm::LogWarning("MissingJetConstituent") << "Jet constituent required for N-subjettiness computation is missing!";
        }
      }
      else
        fjParticles.push_back( fastjet::PseudoJet( daughter->px(), daughter->py(), daughter->pz(), daughter->energy() ) );
    }
    else
      edm::LogWarning("MissingJetConstituent") << "Jet constituent required for N-subjettiness computation is missing!";
  }

  // N-subjettiness calculator
  fastjet::contrib::Njettiness njettiness(fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta_,R0_));

  // calculate N-subjettiness
  tau1 = njettiness.getTau(1, fjParticles);
  tau2 = njettiness.getTau(2, fjParticles);
  currentAxes = njettiness.currentAxes();
}

// -------------- template specializations --------------------
// -------------- fillHelpers ----------------
template<>
void BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::fillHelpers(const edm::Event& iEvent)
{
  // do nothing
}

template<>
void BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::fillHelpers(const edm::Event& iEvent)
{
  std::vector<reco::JetTagInfo> jetTagInfos = m_helper.makeBaseVector(iEvent);
  if ( runSubJets_ ) {
    for ( size_t i = 0; i < SubJetLabels_.size(); ++i )
      std::vector<reco::JetTagInfo> subJetTagInfos = m_subjetHelper[i].makeBaseVector(iEvent);
  }
}

// -------------- toIPTagInfo ----------------
template<>
const BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::IPTagInfo *
BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::toIPTagInfo(const pat::Jet & jet, const std::string & tagInfos)
{
  return jet.tagInfoTrackIP(tagInfos.c_str());
}

template<>
const BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::IPTagInfo *
BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::toIPTagInfo(const pat::Jet & jet, const std::string & tagInfos)
{
  return jet.tagInfoCandIP(tagInfos.c_str());
}

// -------------- toSVTagInfo ----------------
template<>
const BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::SVTagInfo *
BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::toSVTagInfo(const pat::Jet & jet, const std::string & tagInfos)
{
  return jet.tagInfoSecondaryVertex(tagInfos.c_str());
}

template<>
const BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::SVTagInfo *
BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::toSVTagInfo(const pat::Jet & jet, const std::string & tagInfos)
{
  return jet.tagInfoCandSecondaryVertex(tagInfos.c_str());
}

// -------------- toAllTracks ----------------
template<>
const BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::Tracks
BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::toAllTracks(const pat::Jet & jet, const std::string & tagInfos, const reco::JetTagInfo & jetTagInfo, const int & iJetColl)
{
  return toIPTagInfo(jet,tagInfos)->tracks();
}

template<>
const BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::Tracks
BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::toAllTracks(const pat::Jet & jet, const std::string & tagInfos, const reco::JetTagInfo & jetTagInfo, const int & iJetColl)
{
  if( iJetColl > 0)
    return m_subjetHelper[iJetColl-1].tracks(jetTagInfo);
  else
    return m_helper.tracks(jetTagInfo);
}

// -------------- setTracksPV ----------------
template<>
void BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight)
{
  setTracksPVBase(trackRef, pvHandle, iPV, PVweight);
}

template<>
void BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight)
{
  iPV = -1;
  PVweight = 0.;

  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(trackRef.get());

  if(pcand) // MiniAOD case
  {
    if( pcand->fromPV() == pat::PackedCandidate::PVUsedInFit )
    {
      iPV = 0;
      PVweight = 1.;
    }
  }
  else
  {
    const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(trackRef.get());

    setTracksPVBase(pfcand->trackRef(), pvHandle, iPV, PVweight);
  }
}

// -------------- setTracksSV ----------------
template<>
void BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::setTracksSV(const TrackRef & trackRef, const SVTagInfo * svTagInfo, int & isFromSV, int & iSV, float & SVweight)
{
  isFromSV = 0;
  iSV = -1;
  SVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::Vertex::trackRef_iterator IT;

  size_t nSV = svTagInfo->nVertices();
  for(size_t iv=0; iv<nSV; ++iv)
  {
    const reco::Vertex & vtx = svTagInfo->secondaryVertex(iv);
    // loop over tracks in vertices
    for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
    {
      const reco::TrackBaseRef & baseRef = *it;
      // one of the tracks in the vertex is the same as the track considered in the function
      if( baseRef == trackBaseRef )
      {
        float w = vtx.trackWeight(baseRef);
        // select the vertex for which the track has the highest weight
        if( w > SVweight )
        {
          SVweight = w;
          isFromSV = 1;
          iSV = iv;
          break;
        }
      }
    }
  }
}

template<>
void BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::setTracksSV(const TrackRef & trackRef, const SVTagInfo * svTagInfo, int & isFromSV, int & iSV, float & SVweight)
{
  isFromSV = 0;
  iSV = -1;
  SVweight = 0.;

  typedef std::vector<reco::CandidatePtr>::const_iterator IT;

  size_t nSV = svTagInfo->nVertices();
  for(size_t iv=0; iv<nSV; ++iv)
  {
    const Vertex & vtx = svTagInfo->secondaryVertex(iv);
    const std::vector<reco::CandidatePtr> & tracks = vtx.daughterPtrVector();

    // one of the tracks in the vertex is the same as the track considered in the function
    if( std::find(tracks.begin(),tracks.end(),trackRef) != tracks.end() )
    {
      SVweight = 1.;
      isFromSV = 1;
      iSV = iv;
    }

    // select the first vertex for which the track is used in the fit
    // (reco::VertexCompositePtrCandidate does not store track weights so can't select the vertex for which the track has the highest weight)
    if(iSV>=0)
      break;
  }
}

// -------------- vertexKinematicsAndChange ----------------
template<>
void BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex>::vertexKinematicsAndChange(const Vertex & vertex, reco::TrackKinematics & vertexKinematics, Int_t & charge)
{
  Bool_t hasRefittedTracks = vertex.hasRefittedTracks();

  for(reco::Vertex::trackRef_iterator track = vertex.tracks_begin();
      track != vertex.tracks_end(); ++track) {
    Double_t w = vertex.trackWeight(*track);
    if (w < 0.5)
      continue;
    if (hasRefittedTracks) {
      reco::Track actualTrack = vertex.refittedTrack(*track);
      vertexKinematics.add(actualTrack, w);
      charge+=actualTrack.charge();
    }
    else {
      const reco::Track& mytrack = **track;
      vertexKinematics.add(mytrack, w);
      charge+=mytrack.charge();
    }
  }
}

template<>
void BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate>::vertexKinematicsAndChange(const Vertex & vertex, reco::TrackKinematics & vertexKinematics, Int_t & charge)
{
  const std::vector<reco::CandidatePtr> & tracks = vertex.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track) {
    const reco::Track& mytrack = *(*track)->bestTrack();
    vertexKinematics.add(mytrack, 1.0);
    charge+=mytrack.charge();
  }
}

// define specific instances of the templated BTagAnalyzer
typedef BTagAnalyzerT<reco::TrackIPTagInfo,reco::Vertex> BTagAnalyzerLegacy;
typedef BTagAnalyzerT<reco::CandIPTagInfo,reco::VertexCompositePtrCandidate> BTagAnalyzer;
//define plugins
DEFINE_FWK_MODULE(BTagAnalyzerLegacy);
DEFINE_FWK_MODULE(BTagAnalyzer);

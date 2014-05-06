#ifndef BTagAnalyzer_h
#define BTagAnalyzer_h

// -*- C++ -*-
//
// Package:    BTagAnalyzer
// Class:      BTagAnalyzer
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
// FastJet
#include "fastjet/PseudoJet.hh"
// N-subjettiness
#include "RecoBTag/PerformanceMeasurements/interface/Njettiness.hh"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "RecoBTag/PerformanceMeasurements/interface/CategoryFinder.h"

// reco track and vertex
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

#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackHistory/interface/TrackCategories.h"
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"

#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// reconstruct IP
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "FWCore/Utilities/interface/RegexMatch.h"
#include <boost/regex.hpp>

#include "RecoBTag/PerformanceMeasurements/interface/JetInfoBranches.h"
#include "RecoBTag/PerformanceMeasurements/interface/EventInfoBranches.h"
#include "RecoBTag/PerformanceMeasurements/interface/BookHistograms.h"

//
// constants, enums and typedefs
//
typedef std::vector<pat::Jet> PatJetCollection;
typedef std::map<const pat::Jet*,const pat::Jet*> JetToJetMap;

//
// class declaration
//

using namespace std;
using namespace reco;

const UInt_t MAX_JETCOLLECTIONS=2;

class BTagAnalyzer : public edm::EDAnalyzer
{
  public:
    explicit BTagAnalyzer(const edm::ParameterSet&);
    ~BTagAnalyzer();

    // auxiliary class holding simulated primary vertices
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

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    bool findCat(const reco::Track* ,CategoryFinder& );
    std::vector< float > getTrackProbabilies(std::vector< float > , const int );
    double calculProbability(std::vector< float > );

    float calculPtRel(const reco::Track& theMuon, const pat::Jet& theJet);

    int matchMuon(const edm::RefToBase<reco::Track>& theMuon, const edm::View<reco::Muon>& muons);

    void setTracksPV( const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight );

    void setTracksSV( const reco::TrackRef & trackRef, const reco::SecondaryVertexTagInfo *, int & isFromSV, int & iSV, float & SVweight );

    int getMuonTk(double pt, const int);
    bool NameCompatible(const std::string& pattern, const std::string& name);

    std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::HepMCProduct>& evtMC);

    void processTrig(const edm::Handle<edm::TriggerResults>&, const std::vector<std::string>&) ;

    void processJets(const edm::Handle<PatJetCollection>&, const edm::Handle<PatJetCollection>&,
        const edm::Event&, const edm::EventSetup&, const JetToJetMap&, const int) ;

    int isFromGSP(const reco::Candidate* c);

    bool isHardProcess(const int status);

    // ----------member data ---------------------------
    std::string outputFile_;
    //std::vector< std::string > moduleLabel_;

    bool runSubJets_ ;
    bool allowJetSkipping_ ;
    bool storeTagVariables_;
    bool storeCSVTagVariables_;

    edm::InputTag src_;  // Generator/handronizer module label
    edm::InputTag muonCollectionName_;
    edm::InputTag patMuonCollectionName_;
    edm::InputTag triggerTable_;
    edm::InputTag SVComputer_;

    edm::InputTag JetCollectionTag_;
    edm::InputTag FatJetCollectionTag_;
    edm::InputTag PrunedFatJetCollectionTag_;

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

    std::string combinedSVRetrainedBJetTags_;
    std::string combinedSVRetrainedNegBJetTags_;
    std::string combinedSVRetrainedPosBJetTags_;

    std::string combinedCSVJPBJetTags_;
    std::string combinedCSVJPNegBJetTags_;
    std::string combinedCSVJPPosBJetTags_;

    std::string combinedCSVSLBJetTags_;
    std::string combinedCSVSLNegBJetTags_;
    std::string combinedCSVSLPosBJetTags_;

    std::string combinedCSVJPSLBJetTags_;
    std::string combinedCSVJPSLNegBJetTags_;
    std::string combinedCSVJPSLPosBJetTags_;

    std::string simpleIVFSVHighPurBJetTags_;
    std::string simpleIVFSVHighEffBJetTags_;
    std::string doubleIVFSVHighEffBJetTags_;
    std::string combinedIVFSVBJetTags_;
    std::string combinedIVFSVPosBJetTags_;

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

    std::string softPFMuonTagInfos_;
    std::string softPFElectronTagInfos_;

    edm::InputTag primaryVertexColl_;
    edm::InputTag tracksColl_;

    bool useTrackHistory_;
    TFile*  rootFile_;
    double minJetPt_;
    double maxJetEta_;

    int selTagger_;
    bool isData_;
    bool useSelectedTracks_;
    bool produceJetTrackTree_;
    bool produceAllTrackTree_;
    bool producePtRelTemplate_;

    bool use_ttbar_filter_;
    edm::InputTag channel_;

    // trigger list
    std::vector<std::string> triggerPathNames_;

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
    JetInfoBranches JetInfo[MAX_JETCOLLECTIONS] ;
    BookHistograms* Histos[MAX_JETCOLLECTIONS] ;

    std::vector<std::string> ttbarTriggerPathNames_;

    const  reco::Vertex  *pv;

    bool PFJet80 ;
    std::vector<std::string> PFJet80TriggerPathNames_;

    TLorentzVector thelepton1;
    TLorentzVector thelepton2;

    const GenericMVAJetTagComputer *computer ;

    Njettiness nsubjettinessCalculator;

    edm::View<reco::Muon> muons ;

    edm::ESHandle<TransientTrackBuilder> trackBuilder ;
    edm::Handle<reco::VertexCollection> primaryVertex ;

    int cap0, cap1, cap2, cap3, cap4, cap5, cap6, cap7, cap8;
    int can0, can1, can2, can3, can4, can5, can6, can7, can8;

    // Generator/hadronizer type (information stored bitwise)
    unsigned int hadronizerType_;
};

#endif

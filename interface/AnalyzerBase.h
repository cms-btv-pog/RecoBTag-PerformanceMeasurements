#ifndef ANALYZERBASE_H
#define ANALYZERBASE_H

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

//#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
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
#include "DataFormats/BTauReco/interface/DeepFlavourFeatures.h"
#include "DataFormats/BTauReco/interface/DeepDoubleXFeatures.h"
#include "DataFormats/BTauReco/interface/DeepBoostedJetFeatures.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

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

// Math clusters to TrackingParticles
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoBTag/SecondaryVertex/interface/V0Filter.h"
// #include "RecoBTag/ImpactParameter/plugins/IPProducer.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "FWCore/Utilities/interface/RegexMatch.h"
#include <boost/regex.hpp>

#include "RecoBTag/PerformanceMeasurements/interface/JetInfoBranches.h"
#include "RecoBTag/PerformanceMeasurements/interface/EventInfoBranches.h"
#include "RecoBTag/PerformanceMeasurements/interface/BookHistograms.h"
#include "RecoBTag/PerformanceMeasurements/interface/VariableParser.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Common/interface/Provenance.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "fastjet/contrib/Njettiness.hh"

// #include "RecoBTag/SecondaryVertex/interface/CombinedSVSoftLeptonComputer.h"
#include "DataFormats/BTauReco/interface/DeepFlavourTagInfo.h"
#include "DataFormats/BTauReco/interface/DeepDoubleXTagInfo.h"
#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"



namespace analyzerBase{

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

    const reco::TrackBaseRef toTrackRef(const reco::TrackRef & trk);
    const reco::TrackBaseRef toTrackRef(const edm::Ptr<reco::Candidate> & cnd);

    const math::XYZPoint & position(const reco::Vertex & sv);
    const math::XYZPoint & position(const reco::VertexCompositePtrCandidate & sv);
    const double xError(const reco::Vertex & sv);
    const double xError(const reco::VertexCompositePtrCandidate & sv);
    const double yError(const reco::Vertex & sv);
    const double yError(const reco::VertexCompositePtrCandidate & sv);
    const double zError(const reco::Vertex & sv);
    const double zError(const reco::VertexCompositePtrCandidate & sv);
    const double chi2(const reco::Vertex & sv);
    const double chi2(const reco::VertexCompositePtrCandidate & sv);
    const double ndof(const reco::Vertex & sv);
    const double ndof(const reco::VertexCompositePtrCandidate & sv);
    const unsigned int vtxTracks(const reco::Vertex & sv);
    const unsigned int vtxTracks(const reco::VertexCompositePtrCandidate & sv);

    // This is needed to get a TrackingParticle --> Cluster match (instead of Cluster-->TP) (only needed in processJets)
    using P = std::pair<OmniClusterRef, TrackingParticleRef>;
    bool compare(const P& i, const P& j);


    enum JetFlavor {UNDEFINED, G, UD, S, C, GCC, CC, B, GBB, BB, LeptonicB, LeptonicB_C, TAU};

    JetFlavor jet_flavour(const pat::Jet& jet,
        const std::vector<reco::GenParticle>& gToBB,
        const std::vector<reco::GenParticle>& gToCC,
        const std::vector<reco::GenParticle>& neutrinosLepB,
        const std::vector<reco::GenParticle>& neutrinosLepB_C,
        const std::vector<reco::GenParticle>& alltaus,
        bool usePhysForLightAndUndefined=false);

}
#endif

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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
//#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

//#include "RecoBTag/MCTools/interface/JetFlavour.h"
//#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CategoryFinder.h"

// reco track and vertex
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
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
#include "TNtuple.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
//#include "RecoBTag/MCTools/interface/JetFlavour.h"
//#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

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

//residual jet corrections
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//reconstruct IP
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

struct ltstr
{
    bool operator()(const edm::RefToBase<reco::Jet> s1, edm::RefToBase<reco::Jet> s2) const
    {
        if (s1.id() != s2.id()) return s1.id()<s2.id();
        return s1.key()< s2.key();
    }
};


//
// class declaration
//
//using BTagMCTools::JetFlavour;
using namespace edm;
using namespace reco;
//using namespace BTagMCTools;

const UInt_t nMaxTrk_=100000;
const UInt_t nMaxJets_=10000;


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

    float calculPtRel();

    bool findCat(const reco::Track* ,CategoryFinder& );
    std::vector< float > getTrackProbabilies(std::vector< float > , int );
    double calculProbability(std::vector< float > );

    reco::JetFlavour getMatchedParton(const reco::Jet &jet);
    float calculPtRel(reco::Track theMuon, reco::Jet theJet, double JES , string jetcoll );
    
    int matchMuon(const edm::RefToBase<reco::Track>& theMuon, edm::View<reco::Muon>& muons);

    void setTracksPV( const reco::Vertex *pv, bool isPV );
    
    double getGenJetPt(reco::Jet theJet, GenJetCollection &theJets);
    int getMuonTk(double pt);
    bool NameCompatible(const std::string& pattern, const std::string& name);
    
    // ----------member data ---------------------------
    std::string outputFile_;
    //std::vector< std::string > moduleLabel_;

    std::string flavourMatchOptionf;
    edm::InputTag flavourSourcef;
    edm::InputTag muonCollectionName_;
    edm::InputTag triggerTable_;
    edm::InputTag SVComputer_;
    
    //JetFlavourIdentifier jetFlavourIdentifier_;

    std::string CaloJetCollectionTags_;
    std::string jetCorrector_;

    std::string jetPModuleName_;
    std::string jetPPosModuleName_;
    std::string jetPNegModuleName_;

    std::string jetBModuleName_;
    std::string jetBNegModuleName_;
    std::string jetBPosModuleName_;

    std::string trackCHEModuleName_;
    std::string trackCNegHEModuleName_;

    std::string trackCHPModuleName_;
    std::string trackCNegHPModuleName_;

    std::string combinedSvtxModuleName_;
    std::string combinedSvtxPosModuleName_;
    std::string combinedSvtxNegModuleName_;

    std::string svtxModuleNameHighPur_;
    std::string svtxNegModuleNameHighPur_;
    std::string svtxModuleNameHighEff_;
    std::string svtxNegModuleNameHighEff_;

    std::string softMuonModuleName_;
    std::string softMuonNegModuleName_;
    std::string softMuonTagInfoName_;

    std::string primaryVertexColl_;

    bool useTrackHistory_;
    TFile*  rootFile_;
    double minJetPt_;
    double maxJetEta_;


    int selTagger_;
    double tagCut_;
    double vetoPos_;
    int ntrackMin_;
    bool isData_;
    bool use_selected_tracks_;
    bool produceJetProbaTree_;
    bool producePtRelTemplate_;

    // trigger list
    std::vector<std::string> triggernames_;
    bool TriggerInfo_;
    std::string genJetCollection_;
    
    std::map<edm::RefToBase<reco::Jet>, unsigned int, ltstr> flavoursMapf;
    edm::Handle<reco::JetFlavourMatchingCollection> theJetPartonMapf;
    int TaggedJet(reco::Jet , edm::Handle<reco::JetTagCollection >  );

    TNtuple* nTuplesJets;
    TrackClassifier classifier_;

    ///////////////
    // Some Histograms
    
    edm::Service<TFileService> fs;

    TH1F* hData_All_NJets       ;
    TH1F* hData_All_NTracks     ;
    TH1F* hData_All_JetPt       ;
    TH1F* hData_All_JetEta      ;
    TH1F* hData_NJets           ;
    TH1F* hData_NTracks         ;
    TH1F* hData_JetPt           ;
    TH1F* hData_JetEta          ;
    TH1F* hData_Tagger          ;
    TH1F* hData_Tagger_TCHE     ;
    TH1F* hData_Tagger_TCHP     ;
    TH1F* hData_Tagger_JP	;
    TH1F* hData_Tagger_SSVHE    ;
    TH1F* hData_Tagger_SSVHP    ;
    TH1F* hData_Tagger_CSV      ;
    TH1F* hData_Tagger_MU	;

    TH1F* hAllFlav_Flavour         ;
    TH1F* hAllFlav_Tagger          ;
    TH1F* hAllFlav_Tagger_Gam      ;
    TH1F* hAllFlav_Tagger_K0s      ;
    TH1F* hAllFlav_Tagger_Lam      ;
    TH1F* hAllFlav_Tagger_Bwd      ;
    TH1F* hAllFlav_Tagger_Cwd      ;
    TH1F* hAllFlav_Tagger_Tau      ;
    TH1F* hAllFlav_Tagger_Int      ;
    TH1F* hAllFlav_Tagger_Fak      ;
    TH1F* hAllFlav_Tagger_Bad      ;
    TH1F* hAllFlav_Tagger_Oth      ;

    TH1F* hLightFlav_Tagger        ;
    TH1F* hGluonFlav_Tagger        ;
    TH1F* hUDSFlav_Tagger          ;
    TH1F* hCFlav_Tagger            ;
    TH1F* hBFlav_Tagger            ;

    TH1F*  IPSign_cat0  ;
    TH1F*  IPSign_cat1  ;
    TH1F*  IPSign_cat2  ;
    TH1F*  IPSign_cat3  ;
    TH1F*  IPSign_cat4  ;
    TH1F*  IPSign_cat5  ;
    TH1F*  IPSign_cat6  ;
    TH1F*  IPSign_cat7  ;
    TH1F*  IPSign_cat8  ;
    TH1F*  IPSign_cat9  ;
    TH1F*  TrackProbaNeg ;
    TH1F*  TrackProbaNeg_Cat0 ;
    TH1F*  TrackProbaNeg_Cat1 ;
    TH1F*  TrackProbaNeg_Cat2 ;
    TH1F*  TrackProbaNeg_Cat3 ;
    TH1F*  TrackProbaNeg_Cat4 ;
    TH1F*  TrackProbaNeg_Cat5 ;
    TH1F*  TrackProbaNeg_Cat6 ;
    TH1F*  TrackProbaNeg_Cat7 ;
    TH1F*  TrackProbaNeg_Cat8 ;
    TH1F*  TrackProbaNeg_Cat9 ;

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
    
    int   nTrack;
    float Track_dxy[nMaxTrk_];
    float Track_dz[nMaxTrk_];
    float Track_zIP[nMaxTrk_];
    float Track_length[nMaxTrk_];
    float Track_dist[nMaxTrk_];
    float Track_IP2D[nMaxTrk_];
    float Track_IP2Dsig[nMaxTrk_];
    float Track_IP2Derr[nMaxTrk_];
    float Track_IP[nMaxTrk_];
    float Track_IPsig[nMaxTrk_];
    float Track_IPerr[nMaxTrk_];
    float Track_Proba[nMaxTrk_];
    float Track_p[nMaxTrk_];
    float Track_pt[nMaxTrk_];
    float Track_eta[nMaxTrk_];
    float Track_phi[nMaxTrk_];
    float Track_chi2[nMaxTrk_];
    int   Track_charge[nMaxTrk_];
    int   Track_history[nMaxTrk_];
    int   Track_nHitStrip[nMaxTrk_];
    int   Track_nHitPixel[nMaxTrk_];
    int   Track_nHitAll[nMaxTrk_];
    int   Track_nHitTIB[nMaxTrk_];
    int   Track_nHitTID[nMaxTrk_];
    int   Track_nHitTOB[nMaxTrk_];
    int   Track_nHitTEC[nMaxTrk_];
    int   Track_nHitPXB[nMaxTrk_];
    int   Track_nHitPXF[nMaxTrk_];
    int   Track_isHitL1[nMaxTrk_];
    int   Track_PV[nMaxTrk_];
    int   Track_SV[nMaxTrk_];
    int   Track_isfromSV[nMaxTrk_];
    float Track_PVweight[nMaxTrk_];
    float Track_SVweight[nMaxTrk_];
    int   Track_category[nMaxTrk_];
    
    
    
    

    int nJet;
    float Jet_pt[nMaxJets_];
    float Jet_genpt[nMaxJets_];
    float Jet_residual[nMaxJets_];
    float Jet_jes[nMaxJets_];
    float Jet_eta[nMaxJets_];
    float Jet_phi[nMaxJets_];
    float Jet_Ip1N[nMaxJets_];
    float Jet_Ip1P[nMaxJets_];
    float Jet_Ip2N[nMaxJets_];
    float Jet_Ip2P[nMaxJets_];
    float Jet_Ip3N[nMaxJets_];
    float Jet_Ip3P[nMaxJets_];
    float Jet_Ip4N[nMaxJets_];
    float Jet_Ip4P[nMaxJets_];
    float Jet_Mass4N[nMaxJets_];
    float Jet_Mass4P[nMaxJets_];
    float Jet_ProbaN[nMaxJets_];
    float Jet_ProbaP[nMaxJets_];
    float Jet_Proba[nMaxJets_];
    float Jet_BprobN[nMaxJets_];
    float Jet_Bprob[nMaxJets_];
    float Jet_BprobP[nMaxJets_];
    float Jet_SvxN[nMaxJets_];
    float Jet_Svx[nMaxJets_];
    int   Jet_SvxNTracks[nMaxJets_];
    int   Jet_SvxTracks[nMaxJets_];
    float Jet_SvxNHP[nMaxJets_];
    float Jet_SvxHP[nMaxJets_];
    float Jet_SvxMass[nMaxJets_];
    float Jet_CombSvxN[nMaxJets_];
    float Jet_CombSvxP[nMaxJets_];
    float Jet_CombSvx[nMaxJets_];
    float Jet_SoftMuN[nMaxJets_];
    float Jet_SoftMu[nMaxJets_];
    int   Jet_hist1[nMaxJets_];
    int   Jet_hist2[nMaxJets_];
    int   Jet_hist3[nMaxJets_];
    int   Jet_histJet[nMaxJets_];
    int   Jet_histSvx[nMaxJets_];
    int   Jet_ntracks[nMaxJets_];
    int   Jet_flavour[nMaxJets_];
    int   Jet_nFirstTrack[nMaxJets_];
    int   Jet_nLastTrack[nMaxJets_]; 
    int   Jet_nFirstSV[nMaxJets_];
    int   Jet_nLastSV[nMaxJets_];
    int   Jet_nFirstTrkInc[nMaxJets_];
    int   Jet_nLastTrkInc[nMaxJets_];
    int   Jet_SV_multi[nMaxJets_]; 
    
    int   nMuon;
    int   Muon_IdxJet[10000];
    int   Muon_nMuHit[10000];
    int   Muon_nTkHit[10000];
    int   Muon_nPixHit[10000];
    int   Muon_nOutHit[10000];
    int   Muon_isGlobal[10000];
    int   Muon_nMatched[10000];
    float Muon_chi2[10000];
    float Muon_chi2Tk[10000];
    float Muon_pt[10000];
    float Muon_eta[10000];
    float Muon_phi[10000];
    float Muon_ptrel[10000];
    float Muon_vz[10000];
    int   Muon_hist[10000];
    int   Muon_TrackIdx[10000];
    float Muon_IPsig[10000];
    float Muon_IP[10000];
    float Muon_Proba[10000];
    float Muon_IP2D[10000];
    float Muon_IP2Dsig[10000];
    float Muon_deltaR[10000];
    float Muon_ratio[10000]; 
    
    int   nTrkInc;
    float TrkInc_pt[nMaxTrk_];
    float TrkInc_ptrel[nMaxTrk_];
    float TrkInc_IPsig[nMaxTrk_];
    float TrkInc_IP[nMaxTrk_];
    
    int nPV;
    float PV_x[10000];
    float PV_y[10000];
    float PV_z[10000];
    float PV_ex[10000];
    float PV_ey[10000];
    float PV_ez[10000];
    float PV_chi2[10000];
    float PV_ndf[10000];
    int   PV_isgood[10000];
    int   PV_isfake[10000];
    
    int nSV;
    float SV_x[10000];
    float SV_y[10000];
    float SV_z[10000];
    float SV_ex[10000];
    float SV_ey[10000];
    float SV_ez[10000];
    float SV_chi2[10000];
    float SV_ndf[10000];
    float SV_flight[10000];
    float SV_flightErr[10000];
    float SV_deltaR_jet[10000]; 
    float SV_deltaR_sum_jet[10000];
    float SV_deltaR_sum_dir[10000];
    float SV_energy_ratio[10000];
    float SV_aboveC[10000];	 
    float SV_vtx_pt[10000];
    float SV_flight2D[10000];
    float SV_flight2DErr[10000];
    float SV_totCharge[10000]; 
    float SV_vtxDistJetAxis[10000]; 
    int   SV_nTrk[10000]; 
    int   SV_nTrk_firstVxt[10000];
    float SV_mass[10000];	 
    float SV_vtx_eta[10000];
    float SV_vtx_phi[10000];   
    
    int nPUtrue;                // the true number of pileup interactions that have been added to the event
    int nPU;                    // the number of pileup interactions that have been added to the event
    int   PU_bunch[10000];      // 0 if on time pileup, -1 or +1 if out-of-time
    float PU_z[10000];          // the true primary vertex position along the z axis for each added interaction
    float PU_sumpT_low[10000];  // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
    float PU_sumpT_high[10000]; // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
    int   PU_ntrks_low[10000];  // the number of tracks originating from each interaction, where track pT > low_cu
    int   PU_ntrks_high[10000]; // the number of tracks originating from each interaction, where track pT > high_cut 
    float mcweight;
    math::XYZVector jetVertex[1000]; 
    
    int nCFromGSplit;
    float cFromGSplit_pT[10000];
    float cFromGSplit_eta[10000];
    float cFromGSplit_phi[10000];
    int   cFromGSplit_pdgID[10000];
    
    int nBFromGSplit;
    float bFromGSplit_pT[10000];
    float bFromGSplit_eta[10000];
    float bFromGSplit_phi[10000];
    int   bFromGSplit_pdgID[10000];
    
    int BitTrigger;
    int Run;
    int Evt;
    int LumiBlock;
    float PVz;
    float PVzSim;
    float pthat;
    FactorizedJetCorrector *resJEC_PF ;
    FactorizedJetCorrector *resJEC_JPT;
    FactorizedJetCorrector *resJEC_Calo;  

    std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::HepMCProduct> evtMC);
 
};

#endif

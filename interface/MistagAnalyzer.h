#ifndef MistagAnalyzer_h
#define MistagAnalyzer_h

// -*- C++ -*-
//
// Package:    MistagAnalyzer
// Class:      MistagAnalyzer
//
/**\class MistagAnalyzer MistagAnalyzer.cc RecoBTag/PerformanceMeasurements/plugins/MistagAnalyzer.cc

Description: <one line class summary>

Implementation:
    <Notes on implementation>
*/
//
// Original Author:  Andrea Jeremy
//         Created:  Tue Jul 15 16:55:19 CEST 2008
// $Id: MistagAnalyzer.h,v 1.12 2010/04/01 11:46:51 jandrea Exp $
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

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

//#include "RecoBTag/MCTools/interface/JetFlavour.h"
//#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"


// reco track and vertex
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
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


// trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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

class MistagAnalyzer : public edm::EDAnalyzer
{
public:
    explicit MistagAnalyzer(const edm::ParameterSet&);
    ~MistagAnalyzer();


private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    float calculPtRel();

    reco::JetFlavour getMatchedParton(const reco::Jet &jet);
    float calculPtRel(reco::Track theMuon, reco::Jet theJet, double JES );
//$$
    void setTracksPV( const reco::Vertex *pv, bool isPV );
//$$

    // ----------member data ---------------------------
    std::string outputFile_;
    //std::vector< std::string > moduleLabel_;

    std::string flavourMatchOptionf;
    edm::InputTag flavourSourcef;
    //JetFlavourIdentifier jetFlavourIdentifier_;

//$$    std::string primaryVertexColl_;

    std::string CaloJetCollectionTags_;
    std::string jetCorrector_;
    std::string jetPModuleName_;
    std::string jetPPosModuleName_;
    std::string jetPNegModuleName_;

    std::string trackCHEModuleName_;
    std::string trackCNegHEModuleName_;

    std::string trackCHPModuleName_;
    std::string trackCNegHPModuleName_;

    std::string combinedSvtxModuleName_;
    std::string combinedSvtxNegModuleName_;

    std::string svtxModuleName_;
    std::string svtxNegModuleName_;

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
    bool produceJetProbaTree_;
    // trigger list
    std::vector<std::string> triggernames_;
    bool          TriggerInfo_;

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
    TH1F* hData_Tagger_SSV      ;
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

    TTree *smalltree;
    
    int   nTrack;
    float Track_dxy[10000];
    float Track_dz[10000];
    float Track_zIP[10000];
//$$
    float Track_length[10000];
    float Track_dist[10000];
//$$
    float Track_IP2D[10000];
//$$    float Track_IP2Dsig[10000];
//$$
    float Track_IP[10000];
//$$
    float Track_IPsig[10000];
    float Track_Proba[10000];
    float Track_pt[10000];
    float Track_eta[10000];
    float Track_phi[10000];
    float Track_chi2[10000];
    float Track_charge[10000];
//$$
    float Track_PVweight[10000];
    float Track_SVweight[10000];

    int   Track_history[10000];
    int   Track_category[10000];
    int   Track_PV[10000];
    int   Track_SV[10000];
//$$
    int   Track_nHitStrip[10000];
    int   Track_nHitPixel[10000];
    int   Track_nHitAll[10000];
    int   Track_nHitTIB[10000];
    int   Track_nHitTID[10000];
    int   Track_nHitTOB[10000];
    int   Track_nHitTEC[10000];
    int   Track_nHitPXB[10000];
    int   Track_nHitPXF[10000];
    int   Track_isHitL1[10000];

    float Track_SharedMuSimHits[10000];
    float Track_MatchSimTrackID[10000];
  
    int nJet;
    float Jet_pt[10000];
    float Jet_jes[10000];
    float Jet_eta[10000];
    float Jet_phi[10000];
    float Jet_Ip1N[10000];
    float Jet_Ip1P[10000];
    float Jet_Ip2N[10000];
    float Jet_Ip2P[10000];
    float Jet_Ip3N[10000];
    float Jet_Ip3P[10000];
    float Jet_ProbaN[10000];
    float Jet_ProbaP[10000];
    float Jet_Proba[10000];
    float Jet_SvxN[10000];
    float Jet_Svx[10000];
    float Jet_CombSvxN[10000];
    float Jet_CombSvx[10000];
    float Jet_SoftMuN[10000];
    float Jet_SoftMu[10000];
    float Jet_hist1[10000];
    float Jet_hist2[10000];
    float Jet_hist3[10000];
    float Jet_histJet[10000];
    float Jet_histSvx[10000];
    float Jet_histMuon[10000];
    float Jet_mu_nHit[10000];
    float Jet_mu_chi2[10000];
    float Jet_mu_pt[10000];
    float Jet_mu_ptrel[10000];
    int   Jet_ntracks[10000];
    int   Jet_flavour[10000];
    int   Jet_nFirstTrack[10000];
    int   Jet_nLastTrack[10000]; 
    int   Jet_nFirstSV[10000];
    int   Jet_nLastSV[10000];
    
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
    
//$$    int nSelJets;
    int BitTrigger;
    int Run;
    int Evt;
    int LumiBlock;
//$$
};

#endif

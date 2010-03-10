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
// $Id: MistagAnalyzer.h,v 1.9 2009/10/03 20:00:33 yumiceva Exp $
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


//for triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"



struct ltstr
{
    bool operator()(const edm::RefToBase<reco::Jet> s1, edm::RefToBase<reco::Jet> s2) const
    {
        if (s1.id() != s2.id()) return s1.id()<s2.id();
        return s1.key()< s2.key();
    }
};











//
// class decleration
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
    virtual void beginJob(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    float calculPtRel();

    reco::JetFlavour getMatchedParton(const reco::CaloJet &jet);
    float calculPtRel(reco::Track theMuon, reco::Jet theJet, double JES );
    // ----------member data ---------------------------
    std::string outputFile_;
    //std::vector< std::string > moduleLabel_;


    std::string flavourMatchOptionf;
    edm::InputTag flavourSourcef;
    //JetFlavourIdentifier jetFlavourIdentifier_;

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


    //trigger list
    std::vector<std::string> triggernames_;
    bool          TriggerInfo_;

    std::map<edm::RefToBase<reco::Jet>, unsigned int, ltstr> flavoursMapf;
    edm::Handle<reco::JetFlavourMatchingCollection> theJetPartonMapf;
    int TaggedJet(reco::CaloJet , edm::Handle<reco::JetTagCollection >  );

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
    TH1F* hData_Tagger_JP	      ;
    TH1F* hData_Tagger_SSV      ;
    TH1F* hData_Tagger_CSV      ;
    TH1F* hData_Tagger_MU	      ;

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

    TH1F* hLightFlav_Tagger         ;
    TH1F* hGluonFlav_Tagger          ;
    TH1F* hUDSFlav_Tagger         ;
    TH1F* hCFlav_Tagger             ;
    TH1F* hBFlav_Tagger             ;

    TTree *smalltree;
    
    
    
    
    int   nTrackCand;
    float TrackCand_phi[10000];
    float TrackCand_eta[10000];
    float TrackCand_e[10000];
    float TrackCand_p[10000];
    float TrackCand_charge[10000];
    float TrackCand_dz[10000];
    float TrackCand_d0[10000];
    float TrackCand_pt[10000];
    float TrackCand_chi2[10000];
    float TrackCand_Normchi2[10000];
    
    
      
    float TrackCand_IPsignificance[10000];
    float TrackCand_TransverseIPsignificance[10000];
    float TrackCand_TransverseIP[10000];
    float TrackCand_Proba[10000];
    float TrackCand_DecayLength[10000];
    float TrackCand_DistJetAxis[10000];
    float TrackCand_zIP[10000];
    int   TrackCand_isHitL1[10000];
    int   TrackCand_nHitPixel[10000];
    int   TrackCand_nHitTracker[10000];
    int   TrackCand_nHitTOB[10000];
    int   TrackCand_nHitTIB[10000];
    int   TrackCand_nHitTEC[10000];
    int   TrackCand_nHitTID[10000];
    int   TrackCand_nHitPixelEC[10000];
    int   TrackCand_nHitPixelBL[10000];
    float TrackCand_SharedMuSimHits[10000];
    float TrackCand_MatchSimTrackID[10000];
  
    
    
    
    int nJetCand;
    float Ntagtracks[10000];
    float JetCand_pt[10000];
    float JetCand_jes[10000];
    float JetCand_eta[10000];
    float JetCand_phi[10000];
    int   JetCand_multiplicity[10000];
    int   JetCand_flavor[10000];
    int   JetCand_nFirstTrack[10000];
    int   JetCand_nLastTrack[10000]; 
    float JetCand_Ip1N[10000];
    float JetCand_Ip1P[10000];
    float JetCand_Ip2N[10000];
    float JetCand_Ip2P[10000];
    float JetCand_Ip3N[10000];
    float JetCand_Ip3P[10000];
    float JetCand_ProbaN[10000];
    float JetCand_ProbaP[10000];
    float JetCand_Proba[10000];
    float JetCand_SvtxN[10000];
    float JetCand_Svtx[10000];
    float JetCand_CombinedSvtxN[10000];
    float JetCand_CombinedSvtx[10000];
    float JetCand_SoftMN[10000];
    float JetCand_SoftM[10000];
    float JetCand_Category1[10000];
    float JetCand_Category2[10000];
    float JetCand_Category3[10000];
    float JetCand_CategoryJet[10000];
    float JetCand_CategorySVx[10000];
    float JetCand_CategoryMuon[10000];
    float JetCand_mu_NHits_tracker[10000];
    float JetCand_mu_Chi2[10000];
    float JetCand_mu_pT[10000];
    float JetCand_mu_d0[10000];
    float JetCand_ptRel[10000];
    int   JetCand_nFirstSV[10000];
    int   JetCand_nLastSV[10000];
    
    
    
  
    int nPrimaryV;
    float PrimaryV_x[10000];
    float PrimaryV_y[10000];
    float PrimaryV_z[10000];
    float PrimaryV_ex[10000];
    float PrimaryV_ey[10000];
    float PrimaryV_ez[10000];
    float PrimaryV_chi2[10000];
    float PrimaryV_ndf[10000];
    int   PrimaryV_isgood[10000];
    int   PrimaryV_isfake[10000];
    
    
    
    
     
    int nSecondaryV;
    float SecondaryV_x[10000];
    float SecondaryV_y[10000];
    float SecondaryV_z[10000];
    float SecondaryV_ex[10000];
    float SecondaryV_ey[10000];
    float SecondaryV_ez[10000];
    float SecondaryV_chi2[10000];
    float SecondaryV_ndf[10000];
    float SecondaryV_flightDistance[10000];
    float SecondaryV_flightDistanceError[10000];
    
    
    
    
    
    
    
    
    
    int nSelJets;
    int BitTrigger;
    
    
};

#endif

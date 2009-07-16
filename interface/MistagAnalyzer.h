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
// $Id: MistagAnalyzer.h,v 1.4 2009/07/16 00:47:42 jindal Exp $
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

struct ltstr{
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

class MistagAnalyzer : public edm::EDAnalyzer {
 public:
  explicit MistagAnalyzer(const edm::ParameterSet&);
  ~MistagAnalyzer();
  
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  float calculPtRel();
  
  reco::JetFlavour getMatchedParton(const reco::CaloJet &jet);
  // ----------member data ---------------------------
  std::string outputFile_;
  //std::vector< std::string > moduleLabel_;
  
  
  std::string flavourMatchOptionf;
  edm::InputTag flavourSourcef;
  //  JetFlavourIdentifier jetFlavourIdentifier_;
  
  std::string CaloJetCollectionTags_;
  
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
  
  
  
  std::map<edm::RefToBase<reco::Jet>, unsigned int, ltstr> flavoursMapf;
  edm::Handle<reco::JetFlavourMatchingCollection> theJetPartonMapf;
  int TaggedJet(reco::CaloJet , edm::Handle<reco::JetTagCollection >  );
  
  TNtuple* nTuplesJets;
  TrackClassifier classifier_;
  
  //**********************************
  // Data
  //**********************************
  TH1F* hData_All_NJets       ;
  TH1F* hData_All_NTracks     ;
  TH1F* hData_All_JetPt       ;
  TH1F* hData_All_JetEta      ;
  TH1F* hData_NJets           ;
  TH1F* hData_NTracks         ;
  TH1F* hData_JetPt           ;
  TH1F* hData_JetEta          ;
  TH1F* hData_NegTag_NTracks  ;
  TH1F* hData_NegTag_JetPt    ;
  TH1F* hData_NegTag_JetEta   ;
  TH1F* hData_PosTag_NTracks  ;
  TH1F* hData_PosTag_JetPt    ;
  TH1F* hData_PosTag_JetEta   ;
  TH1F* hData_Tagger          ;
  TH1F* hData_Veto_NTracks        ;
  TH1F* hData_Veto_JetPt          ;
  TH1F* hData_Veto_JetEta         ;
  TH1F* hData_Veto_NegTag_NTracks ;
  TH1F* hData_Veto_NegTag_JetPt   ;
  TH1F* hData_Veto_NegTag_JetEta  ;
  TH1F* hData_Veto_Tagger      ;
  
  //**********************************
  // All flavours in Monte Carlo
  //**********************************
  TH1F* hAllFlav_All_NJets       ;
  TH1F* hAllFlav_All_NTracks     ;
  TH1F* hAllFlav_All_JetPt       ;
  TH1F* hAllFlav_All_JetEta      ;
  TH1F* hAllFlav_All_Flavour     ;
  TH1F* hAllFlav_NJets           ;
  TH1F* hAllFlav_NTracks         ;
  TH1F* hAllFlav_JetPt           ;
  TH1F* hAllFlav_JetEta          ;
  TH1F* hAllFlav_Flavour         ;
  TH1F* hAllFlav_NegTag_NTracks  ;
  TH1F* hAllFlav_NegTag_JetPt    ;
  TH1F* hAllFlav_NegTag_JetEta   ;
  TH1F* hAllFlav_PosTag          ;
  TH1F* hAllFlav_PosTag_NTracks  ;
  TH1F* hAllFlav_PosTag_JetPt    ;
  TH1F* hAllFlav_PosTag_JetEta   ;
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
  
  
  TH1F* hAllFlav_Veto_NTracks        ;
  TH1F* hAllFlav_Veto_JetPt          ;
  TH1F* hAllFlav_Veto_JetEta         ;
  TH1F* hAllFlav_Veto_Flavour        ;
  TH1F* hAllFlav_Veto_NegTag_NTracks ;
  TH1F* hAllFlav_Veto_NegTag_JetPt   ;
  TH1F* hAllFlav_Veto_NegTag_JetEta  ;
  TH1F* hAllFlav_Veto_Tagger         ;
  
  TH1F* hAllFlav_K0s_NTracks        ;
  TH1F* hAllFlav_K0s_JetPt        ;
  TH1F* hAllFlav_K0s_JetEta        ;
  TH1F* hAllFlav_K0s_PosTag_NTracks ;
  TH1F* hAllFlav_K0s_PosTag_JetPt   ;
  TH1F* hAllFlav_K0s_PosTag_JetEta  ;
  TH1F* hAllFlav_K0s_Veto_NegTag_NTracks ;
  TH1F* hAllFlav_K0s_Veto_NegTag_JetPt     ;
  TH1F* hAllFlav_K0s_Veto_NegTag_JetEta  ;
  TH1F* hAllFlav_Lam_NTracks        ;
  TH1F* hAllFlav_Lam_JetPt        ;
  TH1F* hAllFlav_Lam_JetEta        ;
  TH1F* hAllFlav_Lam_PosTag_NTracks ;
  TH1F* hAllFlav_Lam_PosTag_JetPt   ;
  TH1F* hAllFlav_Lam_PosTag_JetEta  ;
  TH1F* hAllFlav_Lam_Veto_NegTag_NTracks ;
  TH1F* hAllFlav_Lam_Veto_NegTag_JetPt     ;
  TH1F* hAllFlav_Lam_Veto_NegTag_JetEta  ;
  TH1F* hAllFlav_Gam_NTracks        ;
  TH1F* hAllFlav_Gam_JetPt        ;
  TH1F* hAllFlav_Gam_JetEta        ;
  TH1F* hAllFlav_Gam_PosTag_NTracks ;
  TH1F* hAllFlav_Gam_PosTag_JetPt   ;
  TH1F* hAllFlav_Gam_PosTag_JetEta  ;
  TH1F* hAllFlav_Gam_Veto_NegTag_NTracks ;
  TH1F* hAllFlav_Gam_Veto_NegTag_JetPt     ;
  TH1F* hAllFlav_Gam_Veto_NegTag_JetEta  ;
  TH1F* hAllFlav_Fak_NTracks        ;
  TH1F* hAllFlav_Fak_JetPt        ;
  TH1F* hAllFlav_Fak_JetEta        ;
  TH1F* hAllFlav_Fak_PosTag_NTracks ;
  TH1F* hAllFlav_Fak_PosTag_JetPt   ;
  TH1F* hAllFlav_Fak_PosTag_JetEta  ;
  TH1F* hAllFlav_Fak_Veto_NegTag_NTracks ;
  TH1F* hAllFlav_Fak_Veto_NegTag_JetPt     ;
  TH1F* hAllFlav_Fak_Veto_NegTag_JetEta  ;
  
  //**********************************
  // udsg-jets
  //**********************************
  TH1F* hLightFlav_All_NTracks    ;
  TH1F* hLightFlav_All_JetPt      ;
  TH1F* hLightFlav_All_JetEta     ;
  TH1F* hLightFlav_NTracks        ;
  TH1F* hLightFlav_JetPt          ;
  TH1F* hLightFlav_JetEta         ;
  TH1F* hLightFlav_NegTag_NTracks ;
  TH1F* hLightFlav_NegTag_JetPt   ;
  TH1F* hLightFlav_NegTag_JetEta  ;
  TH1F* hLightFlav_PosTag         ;
  TH1F* hLightFlav_PosTag_NTracks ;
  TH1F* hLightFlav_PosTag_JetPt   ;
  TH1F* hLightFlav_PosTag_JetEta  ;
  TH1F* hLightFlav_Tagger         ;
  TH1F* hLightFlav_Tagger_Gam      ;
  TH1F* hLightFlav_Tagger_K0s      ;
  TH1F* hLightFlav_Tagger_Lam      ;
  TH1F* hLightFlav_Tagger_Bwd      ;
  TH1F* hLightFlav_Tagger_Cwd      ;
  TH1F* hLightFlav_Tagger_Tau      ;
  TH1F* hLightFlav_Tagger_Int      ;
  TH1F* hLightFlav_Tagger_Fak      ;
  TH1F* hLightFlav_Tagger_Bad      ;
  TH1F* hLightFlav_Tagger_Oth      ;
  TH1F* hLightFlav_Veto_NTracks        ;
  TH1F* hLightFlav_Veto_JetPt          ;
  TH1F* hLightFlav_Veto_JetEta         ;
  TH1F* hLightFlav_Veto_NegTag_NTracks ;
  TH1F* hLightFlav_Veto_NegTag_JetPt   ;
  TH1F* hLightFlav_Veto_NegTag_JetEta  ;
  TH1F* hLightFlav_Veto_Tagger           ;
  
  TH1F* hLightFlav_K0s_NTracks          ;
  TH1F* hLightFlav_K0s_JetPt          ;
  TH1F* hLightFlav_K0s_JetEta          ;
  TH1F* hLightFlav_K0s_PosTag_NTracks ;
  TH1F* hLightFlav_K0s_PosTag_JetPt   ;
  TH1F* hLightFlav_K0s_PosTag_JetEta  ;
  TH1F* hLightFlav_K0s_Veto_Tagger    ;
  TH1F* hLightFlav_Lam_NTracks          ;
  TH1F* hLightFlav_Lam_JetPt          ;
  TH1F* hLightFlav_Lam_JetEta          ;
  TH1F* hLightFlav_Lam_PosTag_NTracks ;
  TH1F* hLightFlav_Lam_PosTag_JetPt   ;
  TH1F* hLightFlav_Lam_PosTag_JetEta  ;
  TH1F* hLightFlav_Lam_Veto_Tagger    ;
  TH1F* hLightFlav_Gam_NTracks          ;
  TH1F* hLightFlav_Gam_JetPt          ;
  TH1F* hLightFlav_Gam_JetEta          ;
  TH1F* hLightFlav_Gam_PosTag_NTracks ;
  TH1F* hLightFlav_Gam_PosTag_JetPt   ;
  TH1F* hLightFlav_Gam_PosTag_JetEta  ;
  TH1F* hLightFlav_Gam_Veto_Tagger    ;
  TH1F* hLightFlav_Fak_NTracks          ;
  TH1F* hLightFlav_Fak_JetPt          ;
  TH1F* hLightFlav_Fak_JetEta          ;
  TH1F* hLightFlav_Fak_PosTag_NTracks ;
  TH1F* hLightFlav_Fak_PosTag_JetPt   ;
  TH1F* hLightFlav_Fak_PosTag_JetEta  ;
  TH1F* hLightFlav_Fak_Veto_Tagger    ;
  
  //**********************************
  // gluon-jets
  //**********************************
  TH1F* hGluonFlav_All_NTracks     ;
  TH1F* hGluonFlav_All_JetPt       ;
  TH1F* hGluonFlav_All_JetEta      ;
  TH1F* hGluonFlav_NTracks         ;
  TH1F* hGluonFlav_JetPt           ;
  TH1F* hGluonFlav_JetEta          ;
  TH1F* hGluonFlav_PosTag          ;
  TH1F* hGluonFlav_PosTag_NTracks  ;
  TH1F* hGluonFlav_PosTag_JetPt    ;
  TH1F* hGluonFlav_PosTag_JetEta   ;
  TH1F* hGluonFlav_Tagger          ;
  TH1F* hGluonFlav_Veto_NTracks        ;
  TH1F* hGluonFlav_Veto_JetPt          ;
  TH1F* hGluonFlav_Veto_JetEta         ;
  TH1F* hGluonFlav_Veto_NegTag_NTracks ;
  TH1F* hGluonFlav_Veto_NegTag_JetPt   ;
  TH1F* hGluonFlav_Veto_NegTag_JetEta  ;
  TH1F* hGluonFlav_Veto_Tagger           ;
  
  //**********************************
  // uds-jets
  //**********************************
  TH1F* hUDSFlav_All_NTracks    ;
  TH1F* hUDSFlav_All_JetPt      ;
  TH1F* hUDSFlav_All_JetEta     ;
  TH1F* hUDSFlav_NTracks        ;
  TH1F* hUDSFlav_JetPt          ;
  TH1F* hUDSFlav_JetEta         ;
  TH1F* hUDSFlav_PosTag         ;
  TH1F* hUDSFlav_PosTag_NTracks ;
  TH1F* hUDSFlav_PosTag_JetPt   ;
  TH1F* hUDSFlav_PosTag_JetEta  ;
  TH1F* hUDSFlav_Tagger         ;
  TH1F* hUDSFlav_Veto_NTracks        ;
  TH1F* hUDSFlav_Veto_JetPt          ;
  TH1F* hUDSFlav_Veto_JetEta         ;
  TH1F* hUDSFlav_Veto_NegTag_NTracks ;
  TH1F* hUDSFlav_Veto_NegTag_JetPt   ;
  TH1F* hUDSFlav_Veto_NegTag_JetEta  ;
  TH1F* hUDSFlav_Veto_Tagger         ;
  
  //**********************************
  // c-jets
  //**********************************
  TH1F* hCFlav_All_NTracks        ;
  TH1F* hCFlav_All_JetPt          ;
  TH1F* hCFlav_All_JetEta         ;
  TH1F* hCFlav_NTracks            ;
  TH1F* hCFlav_JetPt              ;
  TH1F* hCFlav_JetEta             ;
  TH1F* hCFlav_PosTag             ;
  TH1F* hCFlav_PosTag_NTracks     ;
  TH1F* hCFlav_PosTag_JetPt       ;
  TH1F* hCFlav_PosTag_JetEta      ;
  TH1F* hCFlav_Tagger             ;
  TH1F* hCFlav_Veto_NTracks        ;
  TH1F* hCFlav_Veto_JetPt          ;
  TH1F* hCFlav_Veto_JetEta         ;
  TH1F* hCFlav_Veto_NegTag_NTracks ;
  TH1F* hCFlav_Veto_NegTag_JetPt   ;
  TH1F* hCFlav_Veto_NegTag_JetEta  ;
  TH1F* hCFlav_Veto_Tagger       ;
  
  //**********************************
  // b-jets
  //**********************************
  TH1F* hBFlav_All_NTracks        ;
  TH1F* hBFlav_All_JetPt          ;
  TH1F* hBFlav_All_JetEta         ;
  TH1F* hBFlav_NTracks            ;
  TH1F* hBFlav_JetPt              ;
  TH1F* hBFlav_JetEta             ;
  TH1F* hBFlav_PosTag             ;
  TH1F* hBFlav_PosTag_NTracks     ;
  TH1F* hBFlav_PosTag_JetPt       ;
  TH1F* hBFlav_PosTag_JetEta      ;
  TH1F* hBFlav_Tagger             ;
  TH1F* hBFlav_Veto_NTracks        ;
  TH1F* hBFlav_Veto_JetPt          ;
  TH1F* hBFlav_Veto_JetEta         ;
  TH1F* hBFlav_Veto_NegTag_NTracks ;
  TH1F* hBFlav_Veto_NegTag_JetPt   ;
  TH1F* hBFlav_Veto_NegTag_JetEta  ;
  TH1F* hBFlav_Veto_Tagger       ;
  
  
  
  
};

#endif

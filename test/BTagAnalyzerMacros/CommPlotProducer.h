//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 13 15:05:58 2012 by ROOT version 5.32/00
// from TTree ttree/ttree
// found on file: JetTreeTP_Pt-80to120.root
//////////////////////////////////////////////////////////

#ifndef CommPlotProducer_h
#define CommPlotProducer_h
#define ntrack_max 10000

#include "TH1D.h"
#include "TH2D.h"
#include <iomanip>
#include <TROOT.h>
#include <TChain.h>
#include <string.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <map>
#include "../../../../PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
 
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class CommPlotProducer {
public :
   bool isData;
   bool use_selected_tracks;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   bool            produceJetProbaTree;
   bool            produceNewAlgoTree;
   double           n15_20, n20_30, n30_50,n50_80,n80_120,n120_170,n170_300,n300_470,n470_600,n600_800, n800_1000, n1000_inf;
   double           n15_30, n120_150, n150_inf;
   bool            use15_20,use20_30, use30_50,use50_80,use80_120,use120_170,use170_300,use300_470,use470_600,use600_800, use800_1000, use1000_inf;
   bool            use15_30, use120_150, use150_inf; 
   std::vector<TH1D*>   HistoBtag;  
   std::vector<TH2D*>   HistoBtag2D;  
    
//   std::map<TString, int>   HistoBtag_map; 
//   std::map<TString, int>   HistoBtag2D_map; 
   std::map<std::string, int>   HistoBtag_map;
   std::map<std::string, int>   HistoBtag2D_map;

   
   double x_section[12]; 
   double nmc_evt_vect[12];
   float WeightXS;
   float sum_xs;

   TString gentype;
   bool qcdtype;
   int sqrtstev;
     
   int numb_histo;
   int numb_histo2D;   
   
   reweight::LumiReWeighting LumiWeights;
   bool puweight;
   
   // Declaration of leaf types

   Int_t           nBitTrigger;
   Int_t           BitTrigger[3];   //[nBitTrigger]
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Float_t         pthat;
   Float_t         mcweight;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         nPUtrue;
   Int_t           nPU;

   Int_t           ncQuarks;
   Float_t         cQuark_pT[50];   //[ncQuarks]
   Float_t         cQuark_eta[50];   //[ncQuarks]
   Float_t         cQuark_phi[50];   //[ncQuarks]
   Int_t           cQuark_pdgID[50];   //[ncQuarks]
   Int_t           cQuark_status[50];   //[ncQuarks]
   Int_t           cQuark_fromGSP[50];   //[ncQuarks]
   Int_t           nbQuarks;
   Float_t         bQuark_pT[50];   //[nbQuarks]
   Float_t         bQuark_eta[50];   //[nbQuarks]
   Float_t         bQuark_phi[50];   //[nbQuarks]
   Int_t           bQuark_pdgID[50];   //[nbQuarks]
   Int_t           bQuark_status[50];   //[nbQuarks]
   Int_t           bQuark_fromGSP[50];   //[nbQuarks]
   Int_t           nBHadrons;
   Float_t         BHadron_pT[50];   //[nBHadrons]
   Float_t         BHadron_eta[50];   //[nBHadrons]
   Float_t         BHadron_phi[50];   //[nBHadrons]
   Float_t         BHadron_mass[50];   //[nBHadrons]
   Int_t           BHadron_pdgID[50];   //[nBHadrons]
   //Int_t           BHadron_status[50];   //[nBHadrons]
   Int_t           BHadron_mother[50];   //[nBHadrons]
   Int_t           BHadron_hasBdaughter[50];   //[nBHadrons]
   Int_t           nDHadrons;
   Int_t           nDaughters;
   Float_t         DHadron_pT[50];   //[nDHadrons]
   Float_t         DHadron_eta[50];   //[nDHadrons]
   Float_t         DHadron_phi[50];   //[nDHadrons]
   Float_t         DHadron_mass[50];   //[nDHadrons]
   //Float_t         DHadron_vx[50];   //[nDHadrons]
   //Float_t         DHadron_vy[50];   //[nDHadrons]
   //Float_t         DHadron_vz[50];   //[nDHadrons]
   //Float_t         DHadron_daughterVx[50];   //[nDHadrons]
   //Float_t         DHadron_daughterVy[50];   //[nDHadrons]
   //Float_t         DHadron_daughterVz[50];   //[nDHadrons]
   Int_t           DHadron_pdgID[50];   //[nDHadrons]
   Int_t           DHadron_nDaughters[50];   //[nDHadrons]
   Int_t           DHadron_DaughtersPdgID[100];   //[nDaughters]
   Int_t           DHadron_nChargedDaughters[50];   //[nDHadrons]
   Int_t           nGenlep;
   Float_t         Genlep_pT[50];   //[nGenlep]
   Float_t         Genlep_eta[50];   //[nGenlep]
   Float_t         Genlep_phi[50];   //[nGenlep]
   Int_t           Genlep_pdgID[50];   //[nGenlep]
   Int_t           Genlep_status[50];   //[nGenlep]
   Int_t           Genlep_mother[50];   //[nGenlep]
   Int_t           nGenquark;
   Float_t         Genquark_pT[50];   //[nGenquark]
   Float_t         Genquark_eta[50];   //[nGenquark]
   Float_t         Genquark_phi[50];   //[nGenquark]
   Int_t           Genquark_pdgID[50];   //[nGenquark]
   Int_t           Genquark_mother[50];   //[nGenquark]
   Int_t           nJet;
   Float_t         Jet_pt[1000];   //[nJet]
   Float_t         Jet_genpt[1000];   //[nJet]
   Float_t         Jet_residual[1000];   //[nJet]
   Float_t         Jet_jes[1000];   //[nJet]
   Float_t         Jet_eta[1000];   //[nJet]
   Float_t         Jet_phi[1000];   //[nJet]
   Float_t         Jet_mass[1000];   //[nJet]
   Int_t           Jet_ntracks[1000];   //[nJet]
   Int_t           Jet_flavour[1000];   //[nJet]
   Float_t         Jet_Ip2N[1000];   //[nJet]
   Float_t         Jet_Ip2P[1000];   //[nJet]
   Float_t         Jet_Ip3N[1000];   //[nJet]
   Float_t         Jet_Ip3P[1000];   //[nJet]
   Float_t         Jet_ProbaN[1000];   //[nJet]
   Float_t         Jet_ProbaP[1000];   //[nJet]
   Float_t         Jet_Proba[1000];   //[nJet]
   Float_t         Jet_BprobN[1000];   //[nJet]
   Float_t         Jet_BprobP[1000];   //[nJet]
   Float_t         Jet_Bprob[1000];   //[nJet]
   Float_t         Jet_SvxN[1000];   //[nJet]
   Float_t         Jet_Svx[1000];   //[nJet]
   Float_t         Jet_SvxNHP[1000];   //[nJet]
   Float_t         Jet_SvxHP[1000];   //[nJet]
   Float_t         Jet_SvxMass[1000];   //[nJet]
   Float_t         Jet_CombSvxN[1000];   //[nJet]
   Float_t         Jet_CombSvxP[1000];   //[nJet]
   Float_t         Jet_CombSvx[1000];   //[nJet]
   Float_t         Jet_RetCombSvxN[1000];   //[nJet]
   Float_t         Jet_RetCombSvxP[1000];   //[nJet]
   Float_t         Jet_RetCombSvx[1000];   //[nJet]
   //Float_t         Jet_CombCSVJP_N[1000];   //[nJet]
   //Float_t         Jet_CombCSVJP_P[1000];   //[nJet]
   //Float_t         Jet_CombCSVJP[1000];   //[nJet]
   Float_t         Jet_CombCSVSL_N[1000];   //[nJet]
   Float_t         Jet_CombCSVSL_P[1000];   //[nJet]
   Float_t         Jet_CombCSVSL[1000];   //[nJet]
   //Float_t         Jet_CombCSVJPSL_N[1000];   //[nJet]
   //Float_t         Jet_CombCSVJPSL_P[1000];   //[nJet]
   //Float_t         Jet_CombCSVJPSL[1000];   //[nJet]
   //Float_t         Jet_SimpIVF_HP[1000];   //[nJet]
   //Float_t         Jet_SimpIVF_HE[1000];   //[nJet]
   //Float_t         Jet_DoubIVF_HE[1000];   //[nJet]
   Float_t         Jet_CombIVF[1000];   //[nJet]
   //Float_t         Jet_CombIVF_P[1000];   //[nJet]
   Float_t         Jet_SoftMuN[1000];   //[nJet]
   Float_t         Jet_SoftMuP[1000];   //[nJet]
   Float_t         Jet_SoftMu[1000];   //[nJet]
   Float_t         Jet_SoftElN[1000];   //[nJet]
   Float_t         Jet_SoftElP[1000];   //[nJet]
   Float_t         Jet_SoftEl[1000];   //[nJet]
   Int_t           Jet_hist1[1000];   //[nJet]
   Int_t           Jet_hist2[1000];   //[nJet]
   Int_t           Jet_hist3[1000];   //[nJet]
   Int_t           Jet_histJet[1000];   //[nJet]
   Int_t           Jet_histSvx[1000];   //[nJet]
   Int_t           Jet_nFirstTrack[1000];   //[nJet]
   Int_t           Jet_nLastTrack[1000];   //[nJet]
   Int_t           Jet_nFirstSV[1000];   //[nJet]
   Int_t           Jet_nLastSV[1000];   //[nJet]
   Int_t           Jet_SV_multi[1000];   //[nJet]
   Int_t           Jet_nFirstTrkInc[1000];   //[nJet]
   Int_t           Jet_nLastTrkInc[1000];   //[nJet]
   Int_t           Jet_VtxCat[1000];   //[nJet]
   Int_t           Jet_looseID[1000];   //[nJet]
   Int_t           Jet_tightID[1000];   //[nJet]
   Int_t           nTrkInc;
   Float_t         TrkInc_pt[50];   //[nTrkInc]
   Float_t         TrkInc_eta[50];   //[nTrkInc]
   Float_t         TrkInc_phi[50];   //[nTrkInc]
   Float_t         TrkInc_ptrel[50];   //[nTrkInc]
   Float_t         TrkInc_IPsig[50];   //[nTrkInc]
   Float_t         TrkInc_IP[50];   //[nTrkInc]
   Int_t           nMuon;
   Int_t           Muon_IdxJet[50];   //[nMuon]
   Int_t           Muon_nMuHit[50];   //[nMuon]
   Int_t           Muon_nTkHit[50];   //[nMuon]
   Int_t           Muon_nPixHit[50];   //[nMuon]
   Int_t           Muon_nOutHit[50];   //[nMuon]
   Int_t           Muon_isGlobal[50];   //[nMuon]
   Int_t           Muon_nMatched[50];   //[nMuon]
   Float_t         Muon_chi2[50];   //[nMuon]
   Float_t         Muon_chi2Tk[50];   //[nMuon]
   Float_t         Muon_pt[50];   //[nMuon]
   Float_t         Muon_eta[50];   //[nMuon]
   Float_t         Muon_phi[50];   //[nMuon]
   Float_t         Muon_ptrel[50];   //[nMuon]
   Float_t         Muon_vz[50];   //[nMuon]
   Int_t           Muon_hist[50];   //[nMuon]
   Int_t           Muon_TrackIdx[50];   //[nMuon]
   Float_t         Muon_IPsig[50];   //[nMuon]
   Float_t         Muon_IP[50];   //[nMuon]
   Float_t         Muon_IP2Dsig[50];   //[nMuon]
   Float_t         Muon_IP2D[50];   //[nMuon]
   Float_t         Muon_Proba[50];   //[nMuon]
   Float_t         Muon_deltaR[50];   //[nMuon]
   Float_t         Muon_ratio[50];   //[nMuon]
   Float_t         Muon_ratioRel[50];   //[nMuon]
   Int_t           nPFElectron;
   Int_t           PFElectron_IdxJet[50];   //[nPFElectron]
   Float_t         PFElectron_pt[50];   //[nPFElectron]
   Float_t         PFElectron_eta[50];   //[nPFElectron]
   Float_t         PFElectron_phi[50];   //[nPFElectron]
   Float_t         PFElectron_ptrel[50];   //[nPFElectron]
   Float_t         PFElectron_deltaR[50];   //[nPFElectron]
   Float_t         PFElectron_ratio[50];   //[nPFElectron]
   Float_t         PFElectron_ratioRel[50];   //[nPFElectron]
   Float_t         PFElectron_IP[50];   //[nPFElectron]
   Float_t         PFElectron_IP2D[50];   //[nPFElectron]
   Int_t           nPFMuon;
   Int_t           PFMuon_IdxJet[50];   //[nPFMuon]
   Float_t         PFMuon_pt[50];   //[nPFMuon]
   Float_t         PFMuon_eta[50];   //[nPFMuon]
   Float_t         PFMuon_phi[50];   //[nPFMuon]
   Float_t         PFMuon_ptrel[50];   //[nPFMuon]
   Float_t         PFMuon_deltaR[50];   //[nPFMuon]
   Float_t         PFMuon_ratio[50];   //[nPFMuon]
   Float_t         PFMuon_ratioRel[50];   //[nPFMuon]
   Float_t         PFMuon_IP[50];   //[nPFMuon]
   Float_t         PFMuon_IP2D[50];   //[nPFMuon]
   Int_t           PFMuon_GoodQuality[50];   //[nPFMuon]
   Int_t           nSV;
   Float_t         SV_x[50];   //[nSV]
   Float_t         SV_y[50];   //[nSV]
   Float_t         SV_z[50];   //[nSV]
   Float_t         SV_ex[50];   //[nSV]
   Float_t         SV_ey[50];   //[nSV]
   Float_t         SV_ez[50];   //[nSV]
   Float_t         SV_chi2[50];   //[nSV]
   Float_t         SV_ndf[50];   //[nSV]
   Float_t         SV_flight[50];   //[nSV]
   Float_t         SV_flightErr[50];   //[nSV]
   Float_t         SV_deltaR_jet[50];   //[nSV]
   Float_t         SV_deltaR_sum_jet[50];   //[nSV]
   Float_t         SV_deltaR_sum_dir[50];   //[nSV]
   Float_t         SV_energy_ratio[50];   //[nSV]
   Float_t         SV_aboveC[50];   //[nSV]
   Float_t         SV_vtx_pt[50];   //[nSV]
   Float_t         SV_flight2D[50];   //[nSV]
   Float_t         SV_flight2DErr[50];   //[nSV]
   Float_t         SV_totCharge[50];   //[nSV]
   Float_t         SV_vtxDistJetAxis[50];   //[nSV]
   Int_t           SV_nTrk[50];   //[nSV]
   Int_t           SV_nTrk_firstVxt[50];   //[nSV]
   Float_t         SV_mass[50];   //[nSV]
   Float_t         SV_vtx_eta[50];   //[nSV]
   Float_t         SV_vtx_phi[50];   //[nSV]

   
   
   // List of branches for probatree
   
   
   
   
   
   Int_t nTrack; 
   Float_t Track_dxy[ntrack_max];            //[nTrack]
   Float_t Track_dz[ntrack_max];             //[nTrack]
   Float_t Track_zIP[ntrack_max];            //[nTrack]
   Float_t Track_length[ntrack_max];           //[nTrack]
   Float_t Track_dist[ntrack_max];           //[nTrack]
   Float_t Track_IP2D[ntrack_max];           //[nTrack]
   Float_t Track_IP2Dsig[ntrack_max];           //[nTrack]
   Float_t Track_IP2Derr[ntrack_max];        //[nTrack]
   Float_t Track_IP[ntrack_max];           //[nTrack]
   Float_t Track_IPsig[ntrack_max];           //[nTrack]
   Float_t Track_IPerr[ntrack_max];         //[nTrack]
   Float_t Track_Proba[ntrack_max];           //[nTrack]
   Float_t Track_p[ntrack_max];              //[nTrack]  
   Float_t Track_pt[ntrack_max];           //[nTrack]
   Float_t Track_eta[ntrack_max];           //[nTrack]
   Float_t Track_phi[ntrack_max];           //[nTrack]
   Float_t Track_chi2[ntrack_max];           //[nTrack]
   Int_t Track_charge[ntrack_max];           //[nTrack]
   Int_t Track_history[ntrack_max];           //[nTrack]
   Int_t Track_nHitStrip[ntrack_max];           //[nTrack]
   Int_t Track_nHitPixel[ntrack_max];           //[nTrack]
   Int_t Track_nHitAll[ntrack_max];           //[nTrack]
   Int_t Track_nHitTIB[ntrack_max];           //[nTrack]
   Int_t Track_nHitTID[ntrack_max];           //[nTrack]
   Int_t Track_nHitTOB[ntrack_max];           //[nTrack]
   Int_t Track_nHitTEC[ntrack_max];           //[nTrack]
   Int_t Track_nHitPXB[ntrack_max];           //[nTrack]
   Int_t Track_nHitPXF[ntrack_max];           //[nTrack]
   Int_t Track_isHitL1[ntrack_max];           //[nTrack]
   Int_t Track_PV[ntrack_max];               //[nTrack] 
   Int_t Track_SV[ntrack_max];            //[nTrack]     
   Float_t Track_PVweight[ntrack_max];           //[nTrack]
   Float_t Track_SVweight[ntrack_max];           //[nTrack]
   Int_t Track_category[ntrack_max];           //[nTrack]
   Int_t Track_isfromSV[ntrack_max];            //[nTrack]     

   Float_t PV_x[1000];          //[nPV]
   Float_t PV_y[1000];          //[nPV]
   Float_t PV_z[1000];          //[nPV]
   Float_t PV_ex[1000];          //[nPV]
   Float_t PV_ey[1000];          //[nPV]
   Float_t PV_ez[1000];          //[nPV]
   Float_t PV_chi2[1000];          //[nPV]
   Float_t PV_ndf[1000];          //[nPV]
   Int_t PV_isgood[1000];          //[nPV]
   Int_t PV_isfake[1000];          //[nPV]
   

   // List of branches
   TBranch        *b_nBitTrigger;   //!
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!  
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_ncQuarks;   //!
   TBranch        *b_cQuark_pT;   //!
   TBranch        *b_cQuark_eta;   //!
   TBranch        *b_cQuark_phi;   //!
   TBranch        *b_cQuark_pdgID;   //!
   TBranch        *b_cQuark_status;   //!
   TBranch        *b_cQuark_fromGSP;   //!
   TBranch        *b_nbQuarks;   //!
   TBranch        *b_bQuark_pT;   //!
   TBranch        *b_bQuark_eta;   //!
   TBranch        *b_bQuark_phi;   //!
   TBranch        *b_bQuark_pdgID;   //!
   TBranch        *b_bQuark_status;   //!
   TBranch        *b_bQuark_fromGSP;   //!
   TBranch        *b_nBHadrons;   //!
   TBranch        *b_BHadron_pT;   //!
   TBranch        *b_BHadron_eta;   //!
   TBranch        *b_BHadron_phi;   //!
   TBranch        *b_BHadron_mass;   //!
   TBranch        *b_BHadron_pdgID;   //!
   //TBranch        *b_BHadron_status;   //!
   TBranch        *b_BHadron_mother;   //!
   TBranch        *b_BHadron_hasBdaughter;   //!
   TBranch        *b_nDHadrons;   //!
   TBranch        *b_nDaughters;   //!
   TBranch        *b_DHadron_pT;   //!
   TBranch        *b_DHadron_eta;   //!
   TBranch        *b_DHadron_phi;   //!
   TBranch        *b_DHadron_mass;   //!
   //TBranch        *b_DHadron_vx;   //!
   //TBranch        *b_DHadron_vy;   //!
   //TBranch        *b_DHadron_vz;   //!
   //TBranch        *b_DHadron_daughterVx;   //!
   //TBranch        *b_DHadron_daughterVy;   //!
   //TBranch        *b_DHadron_daughterVz;   //!
   TBranch        *b_DHadron_pdgID;   //!
   TBranch        *b_DHadron_nDaughters;   //!
   TBranch        *b_DHadron_DaughtersPdgID;   //!
   TBranch        *b_DHadron_nChargedDaughters;   //!
   TBranch        *b_nGenlep;   //!  
   TBranch        *b_Genlep_pT;   //!
   TBranch        *b_Genlep_eta;   //!
   TBranch        *b_Genlep_phi;   //!
   TBranch        *b_Genlep_pdgID;   //!
   TBranch        *b_Genlep_status;   //!
   TBranch        *b_Genlep_mother;   //!
   TBranch        *b_nGenquark;   //!
   TBranch        *b_Genquark_pT;   //!
   TBranch        *b_Genquark_eta;   //!
   TBranch        *b_Genquark_phi;   //!
   TBranch        *b_Genquark_pdgID;   //!
   TBranch        *b_Genquark_mother;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_genpt;   //!
   TBranch        *b_Jet_residual;   //!
   TBranch        *b_Jet_jes;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_ntracks;   //!
   TBranch        *b_Jet_flavour;   //!
   TBranch        *b_Jet_Ip2N;   //!
   TBranch        *b_Jet_Ip2P;   //!
   TBranch        *b_Jet_Ip3N;   //! 
   TBranch        *b_Jet_Ip3P;   //!  
   TBranch        *b_Jet_ProbaN;   //!
   TBranch        *b_Jet_ProbaP;   //!
   TBranch        *b_Jet_Proba;   //!
   TBranch        *b_Jet_BprobN;   //!
   TBranch        *b_Jet_BprobP;   //!
   TBranch        *b_Jet_Bprob;   //!  
   TBranch        *b_Jet_SvxN;   //!
   TBranch        *b_Jet_Svx;   //!
   TBranch        *b_Jet_SvxNHP;   //!
   TBranch        *b_Jet_SvxHP;   //!
   TBranch        *b_Jet_SvxMass;   //!
   TBranch        *b_Jet_CombSvxN;   //!
   TBranch        *b_Jet_CombSvxP;   //!
   TBranch        *b_Jet_CombSvx;   //!
   TBranch        *b_Jet_RetCombSvxN;   //!
   TBranch        *b_Jet_RetCombSvxP;   //!
   TBranch        *b_Jet_RetCombSvx;   //!
   //TBranch        *b_Jet_CombCSVJP_N;   //!
   //TBranch        *b_Jet_CombCSVJP_P;   //!
   //TBranch        *b_Jet_CombCSVJP;   //!
   TBranch        *b_Jet_CombCSVSL_N;   //!
   TBranch        *b_Jet_CombCSVSL_P;   //!
   TBranch        *b_Jet_CombCSVSL;   //!
   //TBranch        *b_Jet_CombCSVJPSL_N;   //!
   //TBranch        *b_Jet_CombCSVJPSL_P;   //!
   //TBranch        *b_Jet_CombCSVJPSL;   //!
   //TBranch        *b_Jet_SimpIVF_HP;   //!
   //TBranch        *b_Jet_SimpIVF_HE;   //!
   //TBranch        *b_Jet_DoubIVF_HE;   //!
   TBranch        *b_Jet_CombIVF;   //!
   //TBranch        *b_Jet_CombIVF_P;   //!
   TBranch        *b_Jet_SoftMuN;   //!
   TBranch        *b_Jet_SoftMuP;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_SoftElN;   //!
   TBranch        *b_Jet_SoftElP;   //!
   TBranch        *b_Jet_SoftEl;   //!
   TBranch        *b_Jet_hist1;   //!  
   TBranch        *b_Jet_hist2;   //!   
   TBranch        *b_Jet_hist3;   //!   
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //! 
   TBranch        *b_Jet_nFirstTrack;   //!
   TBranch        *b_Jet_nLastTrack;   //!
   TBranch        *b_Jet_nFirstSV;   //!  
   TBranch        *b_Jet_nLastSV;   //!   
   TBranch        *b_Jet_SV_multi;   //!
   TBranch        *b_Jet_nFirstTrkInc;   //!
   TBranch        *b_Jet_nLastTrkInc;   //!
   TBranch        *b_Jet_VtxCat;   //! 
   TBranch        *b_Jet_looseID;   //!
   TBranch        *b_Jet_tightID;   //!
   TBranch        *b_nTrkInc;   //! 
   TBranch        *b_TrkInc_pt;   //!
   TBranch        *b_TrkInc_eta;   //!
   TBranch        *b_TrkInc_phi;   //!
   TBranch        *b_TrkInc_ptrel;   //!
   TBranch        *b_TrkInc_IPsig;   //!
   TBranch        *b_TrkInc_IP;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_IdxJet;   //!
   TBranch        *b_Muon_nMuHit;   //!
   TBranch        *b_Muon_nTkHit;   //!
   TBranch        *b_Muon_nPixHit;   //!
   TBranch        *b_Muon_nOutHit;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_nMatched;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_chi2Tk;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_ptrel;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_hist;   //!
   TBranch        *b_Muon_TrackIdx;   //!
   TBranch        *b_Muon_IPsig;   //!
   TBranch        *b_Muon_IP;   //!
   TBranch        *b_Muon_IP2Dsig;   //!
   TBranch        *b_Muon_IP2D;   //!
   TBranch        *b_Muon_Proba;   //!
   TBranch        *b_Muon_deltaR;   //!
   TBranch        *b_Muon_ratio;   //!
   TBranch        *b_Muon_ratioRel;   //!
   TBranch        *b_nPFElectron;   //!
   TBranch        *b_PFElectron_IdxJet;   //!
   TBranch        *b_PFElectron_pt;   //!
   TBranch        *b_PFElectron_eta;   //!
   TBranch        *b_PFElectron_phi;   //!
   TBranch        *b_PFElectron_ptrel;   //!
   TBranch        *b_PFElectron_deltaR;   //!
   TBranch        *b_PFElectron_ratio;   //!
   TBranch        *b_PFElectron_ratioRel;   //!
   TBranch        *b_PFElectron_IP;   //!
   TBranch        *b_PFElectron_IP2D;   //!
   TBranch        *b_nPFMuon;   //!
   TBranch        *b_PFMuon_IdxJet;   //!
   TBranch        *b_PFMuon_pt;   //!
   TBranch        *b_PFMuon_eta;   //!
   TBranch        *b_PFMuon_phi;   //!
   TBranch        *b_PFMuon_ptrel;   //!
   TBranch        *b_PFMuon_deltaR;   //!
   TBranch        *b_PFMuon_ratio;   //!
   TBranch        *b_PFMuon_ratioRel;   //!
   TBranch        *b_PFMuon_IP;   //!
   TBranch        *b_PFMuon_IP2D;   //!
   TBranch        *b_PFMuon_GoodQuality;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_SV_ex;   //!
   TBranch        *b_SV_ey;   //!
   TBranch        *b_SV_ez;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_ndf;   //!
   TBranch        *b_SV_flight;   //!
   TBranch        *b_SV_flightErr;   //!
   TBranch        *b_SV_deltaR_jet;   //!
   TBranch        *b_SV_deltaR_sum_jet;   //!
   TBranch        *b_SV_deltaR_sum_dir;   //!
   TBranch        *b_SV_energy_ratio;   //!
   TBranch        *b_SV_aboveC;   //! 
   TBranch        *b_SV_vtx_pt;   //!
   TBranch        *b_SV_flight2D;   //! 
   TBranch        *b_SV_flight2DErr;   //!
   TBranch        *b_SV_totCharge;   //!
   TBranch        *b_SV_vtxDistJetAxis;   //!
   TBranch        *b_SV_nTrk;   //!
   TBranch        *b_SV_nTrk_firstVxt;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_vtx_eta;   //!
   TBranch        *b_SV_vtx_phi;   //!

   
  

//   //--------------------------------------
//   // track information 
//   //--------------------------------------
  TBranch *b_nTrack;
  TBranch *b_Trackdxy;
  TBranch *b_Track_dz;   //!
  TBranch *b_Track_zIP;   //!
  TBranch *b_Tracklength;
  TBranch *b_Trackdist;
  TBranch *b_TrackIP2D;
  TBranch *b_TrackIP2Dsig;
  TBranch *b_TrackIP;
  TBranch *b_TrackIPsig;
  TBranch *b_TrackIP2Derr;
  TBranch *b_TrackIPerr;  
  TBranch *b_TrackProba;
  TBranch *b_Trackp;
  TBranch *b_Trackpt;
  TBranch *b_Tracketa;
  TBranch *b_Trackphi;
  TBranch *b_Trackchi2;
  TBranch *b_Trackcharge;
  TBranch *b_Trackhistory;
  TBranch *b_TracknHitStrip;
  TBranch *b_TracknHitPixel;
  TBranch *b_TracknHitAll;
  TBranch *b_TracknHitTIB;
  TBranch *b_TracknHitTID;
  TBranch *b_TracknHitTOB;
  TBranch *b_TracknHitTEC;
  TBranch *b_TracknHitPXB;
  TBranch *b_TracknHitPXF;
  TBranch *b_TrackisHitL1;
  TBranch *b_TrackPV;
  TBranch *b_TrackSV;
  TBranch *b_TrackPVweight;
  TBranch *b_TrackSVweight;
  TBranch *b_Trackcategory;
  TBranch *b_TrackisfromSV;

  //--------------------------------------
  // primary vertex information 
  //--------------------------------------
  TBranch *b_PV_x;
  TBranch *b_PV_y;
  TBranch *b_PV_z;
  TBranch *b_PV_ex;
  TBranch *b_PV_ey;
  TBranch *b_PV_ez;
  TBranch *b_PV_chi2;
  TBranch *b_PV_ndf;
  TBranch *b_PV_isgood;
  TBranch *b_PV_isfake;
  
   CommPlotProducer(TChain *supertree=0, bool infotree1=true, bool infotree2=false, int sqrts=13);
   virtual ~CommPlotProducer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(int trigger, float PtMin_Cut, float PtMax_Cut, TString outputname);
   virtual void     Loop(TString trignam, int trigger, float PtMin_Cut, float PtMax_Cut, TString outputname);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     AddHisto(TString name, TString title,  int nbins, float min, float max);
   virtual void     AddHisto2D(TString name,TString title,int nbins,float min,float max,int nbins2,float min2,float max2);   
   virtual void     FillHisto_int(int flavour, bool isGS, int number, int value, double weight);
   virtual void     FillHisto_float(int flavour, bool isGS, int number, float value, double weight);
   virtual void     FillHisto_floatFromMap(TString name, int flavour, bool isGS, float value, double weight);
   virtual void     FillHisto_intFromMap(TString name, int flavour, bool isGS, int value, double weight);
   virtual void     FillHisto2D_int_floatFromMap(TString name, int flavour, bool isGS, int value, float value2, double weight);
   virtual void     FillHisto2D_float_floatFromMap(TString name, int flavour, bool isGS, float value, float value2, double weight);
   virtual void     SetPU(vector<float> PUvector, TString PUdataFile);
   virtual void     SetPU2012_S10(TString PUdataFile);
   virtual void     SetPU2012_S7(TString PUdataFile);   
   virtual void     SetPV();
   virtual void     SetXS();
   virtual void     Counter();
   virtual double    GetEvtWeight() ; 
   virtual void     SetXS(TString generator, bool MuEnriched, int TeV) ;
   virtual void     SetInfo(TString generator, bool qcd, int TeV) ;
   virtual void     Fill_nevent(double n15,double n20, double n30,double n50,double n80,double n120,double n170,double n300,double n470,double n600, double n800, double n1000); 
   virtual void     SetSumXS();
   
   bool             passMuonSelection(int muidx, int ijet);
   bool             passTrigger(TString trigger, int pttrig);
};
#endif

#ifdef CommPlotProducer_cxx
CommPlotProducer::CommPlotProducer(TChain *superTree, bool infotree1, bool infotree2, int sqrts)
{  

   numb_histo = 0;
   numb_histo2D = 0;
   use15_20   =false;
   use20_30   =false;
   use30_50   =false;
   use50_80   =false;
   use80_120  =false;
   use120_170 =false;
   use170_300 =false;
   use300_470 =false;
   use470_600 =false;
   use600_800 =false;
   use800_1000=false;
   use1000_inf=false;
   use15_30   =false;
   use120_150 =false;
   use150_inf =false;
   
   n15_20   =0;
   n20_30   =0;
   n30_50   =0;
   n50_80   =0;
   n80_120  =0;
   n120_170 =0;
   n170_300 =0;
   n300_470 =0;
   n470_600 =0;
   n600_800 =0;
   n800_1000=0;
   n1000_inf=0;
   n15_30   =0;
   n120_150 =0;
   n150_inf =0;
   sqrtstev=sqrts;
   produceJetProbaTree=infotree1;
   produceNewAlgoTree=infotree2;

   if (produceJetProbaTree) use_selected_tracks=false;
   else use_selected_tracks=true;

   puweight=false;
   
   if (superTree==0) {
      TChain *newchain = new TChain("btagana/ttree");
      newchain->Add("/opt/sbg/data/data1/cms/cbeluffi/BTag_2013_01_10/Commissioning/CMSSW_5_3_2_patch4/src/RecoBTag/PerformanceMeasurements/test/NTuples_80_120/*.root");
      superTree=newchain;
   }


   Init(superTree);
}

CommPlotProducer::~CommPlotProducer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CommPlotProducer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CommPlotProducer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CommPlotProducer::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("nBitTrigger", &nBitTrigger, &b_nBitTrigger);
   fChain->SetBranchAddress("BitTrigger", BitTrigger, &b_BitTrigger);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("ncQuarks", &ncQuarks, &b_ncQuarks);
   fChain->SetBranchAddress("cQuark_pT", cQuark_pT, &b_cQuark_pT);
   fChain->SetBranchAddress("cQuark_eta", cQuark_eta, &b_cQuark_eta);
   fChain->SetBranchAddress("cQuark_phi", cQuark_phi, &b_cQuark_phi);
   fChain->SetBranchAddress("cQuark_pdgID", cQuark_pdgID, &b_cQuark_pdgID);
   fChain->SetBranchAddress("cQuark_status", cQuark_status, &b_cQuark_status);
   fChain->SetBranchAddress("cQuark_fromGSP", cQuark_fromGSP, &b_cQuark_fromGSP);
   fChain->SetBranchAddress("nbQuarks", &nbQuarks, &b_nbQuarks);
   fChain->SetBranchAddress("nbQuarks", &nbQuarks, &b_nbQuarks);
   fChain->SetBranchAddress("bQuark_pT", bQuark_pT, &b_bQuark_pT);
   fChain->SetBranchAddress("bQuark_eta", bQuark_eta, &b_bQuark_eta);
   fChain->SetBranchAddress("bQuark_phi", bQuark_phi, &b_bQuark_phi);
   fChain->SetBranchAddress("bQuark_pdgID", bQuark_pdgID, &b_bQuark_pdgID);
   fChain->SetBranchAddress("bQuark_status", bQuark_status, &b_bQuark_status);
   fChain->SetBranchAddress("bQuark_fromGSP", bQuark_fromGSP, &b_bQuark_fromGSP);
   fChain->SetBranchAddress("nBHadrons", &nBHadrons, &b_nBHadrons);
   fChain->SetBranchAddress("BHadron_pT", BHadron_pT, &b_BHadron_pT);
   fChain->SetBranchAddress("BHadron_eta", BHadron_eta, &b_BHadron_eta);
   fChain->SetBranchAddress("BHadron_phi", BHadron_phi, &b_BHadron_phi);
   fChain->SetBranchAddress("BHadron_mass", BHadron_mass, &b_BHadron_mass);
   fChain->SetBranchAddress("BHadron_pdgID", BHadron_pdgID, &b_BHadron_pdgID);
   //fChain->SetBranchAddress("BHadron_status", BHadron_status, &b_BHadron_status);
   fChain->SetBranchAddress("BHadron_mother", BHadron_mother, &b_BHadron_mother);
   fChain->SetBranchAddress("BHadron_hasBdaughter", BHadron_hasBdaughter, &b_BHadron_hasBdaughter);
   fChain->SetBranchAddress("nDHadrons", &nDHadrons, &b_nDHadrons);
   fChain->SetBranchAddress("nDaughters", &nDaughters, &b_nDaughters);
   fChain->SetBranchAddress("DHadron_pT", DHadron_pT, &b_DHadron_pT);
   fChain->SetBranchAddress("DHadron_eta", DHadron_eta, &b_DHadron_eta);
   fChain->SetBranchAddress("DHadron_phi", DHadron_phi, &b_DHadron_phi);
   fChain->SetBranchAddress("DHadron_mass", DHadron_mass, &b_DHadron_mass);
   //fChain->SetBranchAddress("DHadron_vx", DHadron_vx, &b_DHadron_vx);
   //fChain->SetBranchAddress("DHadron_vy", DHadron_vy, &b_DHadron_vy);
   //fChain->SetBranchAddress("DHadron_vz", DHadron_vz, &b_DHadron_vz);
   //fChain->SetBranchAddress("DHadron_daughterVx", DHadron_daughterVx, &b_DHadron_daughterVx);
   //fChain->SetBranchAddress("DHadron_daughterVy", DHadron_daughterVy, &b_DHadron_daughterVy);
   //fChain->SetBranchAddress("DHadron_daughterVz", DHadron_daughterVz, &b_DHadron_daughterVz);
   fChain->SetBranchAddress("DHadron_pdgID", DHadron_pdgID, &b_DHadron_pdgID);
   fChain->SetBranchAddress("DHadron_nDaughters", DHadron_nDaughters, &b_DHadron_nDaughters);
   fChain->SetBranchAddress("DHadron_DaughtersPdgID", DHadron_DaughtersPdgID, &b_DHadron_DaughtersPdgID);
   fChain->SetBranchAddress("DHadron_nChargedDaughters", DHadron_nChargedDaughters, &b_DHadron_nChargedDaughters);
   fChain->SetBranchAddress("nGenlep", &nGenlep, &b_nGenlep);
   fChain->SetBranchAddress("Genlep_pT", Genlep_pT, &b_Genlep_pT);
   fChain->SetBranchAddress("Genlep_eta", Genlep_eta, &b_Genlep_eta);
   fChain->SetBranchAddress("Genlep_phi", Genlep_phi, &b_Genlep_phi);
   fChain->SetBranchAddress("Genlep_pdgID", Genlep_pdgID, &b_Genlep_pdgID);
   fChain->SetBranchAddress("Genlep_status", Genlep_status, &b_Genlep_status);
   fChain->SetBranchAddress("Genlep_mother", Genlep_mother, &b_Genlep_mother);
   fChain->SetBranchAddress("nGenquark", &nGenquark, &b_nGenquark);
   fChain->SetBranchAddress("Genquark_pT", Genquark_pT, &b_Genquark_pT);
   fChain->SetBranchAddress("Genquark_eta", Genquark_eta, &b_Genquark_eta);
   fChain->SetBranchAddress("Genquark_phi", Genquark_phi, &b_Genquark_phi);
   fChain->SetBranchAddress("Genquark_pdgID", Genquark_pdgID, &b_Genquark_pdgID);
   fChain->SetBranchAddress("Genquark_mother", Genquark_mother, &b_Genquark_mother);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_genpt", Jet_genpt, &b_Jet_genpt);
   fChain->SetBranchAddress("Jet_residual", Jet_residual, &b_Jet_residual);
   fChain->SetBranchAddress("Jet_jes", Jet_jes, &b_Jet_jes);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_ntracks", Jet_ntracks, &b_Jet_ntracks);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   fChain->SetBranchAddress("Jet_Ip2N", Jet_Ip2N, &b_Jet_Ip2N);
   fChain->SetBranchAddress("Jet_Ip2P", Jet_Ip2P, &b_Jet_Ip2P);
   fChain->SetBranchAddress("Jet_Ip3N", Jet_Ip3N, &b_Jet_Ip3N);
   fChain->SetBranchAddress("Jet_Ip3P", Jet_Ip3P, &b_Jet_Ip3P);
   fChain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN, &b_Jet_ProbaN);
   fChain->SetBranchAddress("Jet_ProbaP", Jet_ProbaP, &b_Jet_ProbaP);
   fChain->SetBranchAddress("Jet_Proba", Jet_Proba, &b_Jet_Proba);
   fChain->SetBranchAddress("Jet_BprobN", Jet_BprobN, &b_Jet_BprobN);
   fChain->SetBranchAddress("Jet_BprobP", Jet_BprobP, &b_Jet_BprobP);
   fChain->SetBranchAddress("Jet_Bprob", Jet_Bprob, &b_Jet_Bprob);
   fChain->SetBranchAddress("Jet_SvxN", Jet_SvxN, &b_Jet_SvxN); 
   fChain->SetBranchAddress("Jet_Svx", Jet_Svx, &b_Jet_Svx);
   fChain->SetBranchAddress("Jet_SvxNHP", Jet_SvxNHP, &b_Jet_SvxNHP);
   fChain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP, &b_Jet_SvxHP);
//   fChain->SetBranchAddress("Jet_SvxMass", Jet_SvxMass, &b_Jet_SvxMass);
   fChain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN, &b_Jet_CombSvxN);
   fChain->SetBranchAddress("Jet_CombSvxP", Jet_CombSvxP, &b_Jet_CombSvxP);
   fChain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx, &b_Jet_CombSvx);
   fChain->SetBranchAddress("Jet_RetCombSvxN", Jet_RetCombSvxN, &b_Jet_RetCombSvxN);
   fChain->SetBranchAddress("Jet_RetCombSvxP", Jet_RetCombSvxP, &b_Jet_RetCombSvxP);
   fChain->SetBranchAddress("Jet_RetCombSvx", Jet_RetCombSvx, &b_Jet_RetCombSvx);
   //fChain->SetBranchAddress("Jet_CombCSVJP_N", Jet_CombCSVJP_N, &b_Jet_CombCSVJP_N);
   //fChain->SetBranchAddress("Jet_CombCSVJP_P", Jet_CombCSVJP_P, &b_Jet_CombCSVJP_P);
   //fChain->SetBranchAddress("Jet_CombCSVJP", Jet_CombCSVJP, &b_Jet_CombCSVJP);
   fChain->SetBranchAddress("Jet_CombCSVSL_N", Jet_CombCSVSL_N, &b_Jet_CombCSVSL_N);
   fChain->SetBranchAddress("Jet_CombCSVSL_P", Jet_CombCSVSL_P, &b_Jet_CombCSVSL_P);
   fChain->SetBranchAddress("Jet_CombCSVSL", Jet_CombCSVSL, &b_Jet_CombCSVSL);
   //fChain->SetBranchAddress("Jet_CombCSVJPSL_N", Jet_CombCSVJPSL_N, &b_Jet_CombCSVJPSL_N);
   //fChain->SetBranchAddress("Jet_CombCSVJPSL_P", Jet_CombCSVJPSL_P, &b_Jet_CombCSVJPSL_P);
   //fChain->SetBranchAddress("Jet_CombCSVJPSL", Jet_CombCSVJPSL, &b_Jet_CombCSVJPSL);
   //fChain->SetBranchAddress("Jet_SimpIVF_HP", Jet_SimpIVF_HP, &b_Jet_SimpIVF_HP);
   //fChain->SetBranchAddress("Jet_SimpIVF_HE", Jet_SimpIVF_HE, &b_Jet_SimpIVF_HE);
   //fChain->SetBranchAddress("Jet_DoubIVF_HE", Jet_DoubIVF_HE, &b_Jet_DoubIVF_HE);
   if (sqrtstev==13) fChain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF, &b_Jet_CombIVF);
   //fChain->SetBranchAddress("Jet_CombIVF_P", Jet_CombIVF_P, &b_Jet_CombIVF_P);
   fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
   fChain->SetBranchAddress("Jet_SoftMuP", Jet_SoftMuP, &b_Jet_SoftMuP);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_SoftElN", Jet_SoftElN, &b_Jet_SoftElN);
   fChain->SetBranchAddress("Jet_SoftElP", Jet_SoftElP, &b_Jet_SoftElP);
   fChain->SetBranchAddress("Jet_SoftEl", Jet_SoftEl, &b_Jet_SoftEl);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   fChain->SetBranchAddress("Jet_nFirstTrack", Jet_nFirstTrack, &b_Jet_nFirstTrack);
   fChain->SetBranchAddress("Jet_nLastTrack", Jet_nLastTrack, &b_Jet_nLastTrack);
   fChain->SetBranchAddress("Jet_nFirstSV", Jet_nFirstSV, &b_Jet_nFirstSV);
   fChain->SetBranchAddress("Jet_nLastSV", Jet_nLastSV, &b_Jet_nLastSV);
   fChain->SetBranchAddress("Jet_SV_multi", Jet_SV_multi, &b_Jet_SV_multi);
   fChain->SetBranchAddress("Jet_nFirstTrkInc", Jet_nFirstTrkInc, &b_Jet_nFirstTrkInc);
   fChain->SetBranchAddress("Jet_nLastTrkInc", Jet_nLastTrkInc, &b_Jet_nLastTrkInc);
//   fChain->SetBranchAddress("Jet_VtxCat", Jet_VtxCat, &b_Jet_VtxCat);
   fChain->SetBranchAddress("Jet_looseID", Jet_looseID, &b_Jet_looseID);
   fChain->SetBranchAddress("Jet_tightID", Jet_tightID, &b_Jet_tightID);
   fChain->SetBranchAddress("nTrkInc", &nTrkInc, &b_nTrkInc);
   fChain->SetBranchAddress("TrkInc_pt", &TrkInc_pt, &b_TrkInc_pt);
   fChain->SetBranchAddress("TrkInc_eta", &TrkInc_eta, &b_TrkInc_eta);
   fChain->SetBranchAddress("TrkInc_phi", &TrkInc_phi, &b_TrkInc_phi);
   fChain->SetBranchAddress("TrkInc_ptrel", &TrkInc_ptrel, &b_TrkInc_ptrel);
   fChain->SetBranchAddress("TrkInc_IPsig", &TrkInc_IPsig, &b_TrkInc_IPsig);
   fChain->SetBranchAddress("TrkInc_IP", &TrkInc_IP, &b_TrkInc_IP);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_IdxJet", Muon_IdxJet, &b_Muon_IdxJet);
   fChain->SetBranchAddress("Muon_nMuHit", Muon_nMuHit, &b_Muon_nMuHit);
   fChain->SetBranchAddress("Muon_nTkHit", Muon_nTkHit, &b_Muon_nTkHit);
   fChain->SetBranchAddress("Muon_nPixHit", Muon_nPixHit, &b_Muon_nPixHit);
   fChain->SetBranchAddress("Muon_nOutHit", Muon_nOutHit, &b_Muon_nOutHit);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_nMatched", Muon_nMatched, &b_Muon_nMatched);
   fChain->SetBranchAddress("Muon_chi2", Muon_chi2, &b_Muon_chi2);
   fChain->SetBranchAddress("Muon_chi2Tk", Muon_chi2Tk, &b_Muon_chi2Tk);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_ptrel", Muon_ptrel, &b_Muon_ptrel);
   fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_hist", Muon_hist, &b_Muon_hist);
   fChain->SetBranchAddress("Muon_TrackIdx", Muon_TrackIdx, &b_Muon_TrackIdx);
   fChain->SetBranchAddress("Muon_IPsig", Muon_IPsig, &b_Muon_IPsig);
   fChain->SetBranchAddress("Muon_IP", Muon_IP, &b_Muon_IP);
   fChain->SetBranchAddress("Muon_IP2Dsig", Muon_IP2Dsig, &b_Muon_IP2Dsig);
   fChain->SetBranchAddress("Muon_IP2D", Muon_IP2D, &b_Muon_IP2D);
   fChain->SetBranchAddress("Muon_Proba", Muon_Proba, &b_Muon_Proba);
   fChain->SetBranchAddress("Muon_deltaR", Muon_deltaR, &b_Muon_deltaR);
   fChain->SetBranchAddress("Muon_ratio", Muon_ratio, &b_Muon_ratio);
   fChain->SetBranchAddress("Muon_ratioRel", Muon_ratioRel, &b_Muon_ratioRel);
   fChain->SetBranchAddress("nPFElectron", &nPFElectron, &b_nPFElectron);
   fChain->SetBranchAddress("PFElectron_IdxJet", PFElectron_IdxJet, &b_PFElectron_IdxJet);
   fChain->SetBranchAddress("PFElectron_pt", PFElectron_pt, &b_PFElectron_pt);
   fChain->SetBranchAddress("PFElectron_eta", PFElectron_eta, &b_PFElectron_eta);
   fChain->SetBranchAddress("PFElectron_phi", PFElectron_phi, &b_PFElectron_phi);
   fChain->SetBranchAddress("PFElectron_ptrel", PFElectron_ptrel, &b_PFElectron_ptrel);
   fChain->SetBranchAddress("PFElectron_deltaR", PFElectron_deltaR, &b_PFElectron_deltaR);
   fChain->SetBranchAddress("PFElectron_ratio", PFElectron_ratio, &b_PFElectron_ratio);
   fChain->SetBranchAddress("PFElectron_ratioRel", PFElectron_ratioRel, &b_PFElectron_ratioRel);
   fChain->SetBranchAddress("PFElectron_IP", PFElectron_IP, &b_PFElectron_IP);
   fChain->SetBranchAddress("PFElectron_IP2D", PFElectron_IP2D, &b_PFElectron_IP2D);
   fChain->SetBranchAddress("nPFMuon", &nPFMuon, &b_nPFMuon);
   fChain->SetBranchAddress("PFMuon_IdxJet", PFMuon_IdxJet, &b_PFMuon_IdxJet);
   fChain->SetBranchAddress("PFMuon_pt", PFMuon_pt, &b_PFMuon_pt);
   fChain->SetBranchAddress("PFMuon_eta", PFMuon_eta, &b_PFMuon_eta);
   fChain->SetBranchAddress("PFMuon_phi", PFMuon_phi, &b_PFMuon_phi);
   fChain->SetBranchAddress("PFMuon_ptrel", PFMuon_ptrel, &b_PFMuon_ptrel);
   fChain->SetBranchAddress("PFMuon_deltaR", PFMuon_deltaR, &b_PFMuon_deltaR);
   fChain->SetBranchAddress("PFMuon_ratio", PFMuon_ratio, &b_PFMuon_ratio);
   fChain->SetBranchAddress("PFMuon_ratioRel", PFMuon_ratioRel, &b_PFMuon_ratioRel);
   fChain->SetBranchAddress("PFMuon_IP", PFMuon_IP, &b_PFMuon_IP);
   fChain->SetBranchAddress("PFMuon_IP2D", PFMuon_IP2D, &b_PFMuon_IP2D);
   fChain->SetBranchAddress("PFMuon_GoodQuality", PFMuon_GoodQuality, &b_PFMuon_GoodQuality);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("SV_ex", SV_ex, &b_SV_ex);
   fChain->SetBranchAddress("SV_ey", SV_ey, &b_SV_ey);
   fChain->SetBranchAddress("SV_ez", SV_ez, &b_SV_ez);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_ndf", SV_ndf, &b_SV_ndf);
   fChain->SetBranchAddress("SV_flight", SV_flight, &b_SV_flight);
   fChain->SetBranchAddress("SV_flightErr", SV_flightErr, &b_SV_flightErr);
   fChain->SetBranchAddress("SV_deltaR_jet", SV_deltaR_jet, &b_SV_deltaR_jet);
   fChain->SetBranchAddress("SV_deltaR_sum_jet", SV_deltaR_sum_jet, &b_SV_deltaR_sum_jet);
   fChain->SetBranchAddress("SV_deltaR_sum_dir", SV_deltaR_sum_dir, &b_SV_deltaR_sum_dir);
//   fChain->SetBranchAddress("SV_energy_ratio", SV_energy_ratio, &b_SV_energy_ratio);
//   fChain->SetBranchAddress("SV_aboveC", SV_aboveC, &b_SV_aboveC);
   fChain->SetBranchAddress("SV_vtx_pt", SV_vtx_pt, &b_SV_vtx_pt);
   fChain->SetBranchAddress("SV_flight2D", SV_flight2D, &b_SV_flight2D);
   fChain->SetBranchAddress("SV_flight2DErr", SV_flight2DErr, &b_SV_flight2DErr);
   fChain->SetBranchAddress("SV_totCharge", SV_totCharge, &b_SV_totCharge);
   fChain->SetBranchAddress("SV_vtxDistJetAxis", SV_vtxDistJetAxis, &b_SV_vtxDistJetAxis);
   fChain->SetBranchAddress("SV_nTrk", SV_nTrk, &b_SV_nTrk);
//   fChain->SetBranchAddress("SV_nTrk_firstVxt", SV_nTrk_firstVxt, &b_SV_nTrk_firstVxt);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_vtx_eta", SV_vtx_eta, &b_SV_vtx_eta);
   fChain->SetBranchAddress("SV_vtx_phi", SV_vtx_phi, &b_SV_vtx_phi);

   
   
   
  if ( produceJetProbaTree ) {
//      
//     //--------------------------------------
//   // track information 
//   //--------------------------------------
  fChain->SetBranchAddress("nTrack",	      &nTrack, 	     &b_nTrack);
  fChain->SetBranchAddress("Track_dxy",       Track_dxy,     &b_Trackdxy);	
  fChain->SetBranchAddress("Track_dz",        Track_dz, &b_Track_dz);
  fChain->SetBranchAddress("Track_zIP",       Track_zIP, &b_Track_zIP);
  fChain->SetBranchAddress("Track_length",    Track_length,  &b_Tracklength);
  fChain->SetBranchAddress("Track_dist",      Track_dist,    &b_Trackdist);
  fChain->SetBranchAddress("Track_IP2D",      Track_IP2D,    &b_TrackIP2D);
  fChain->SetBranchAddress("Track_IP2Dsig",   Track_IP2Dsig, &b_TrackIP2Dsig);
  fChain->SetBranchAddress("Track_IP",        Track_IP,	     &b_TrackIP	);
  fChain->SetBranchAddress("Track_IP2Derr",   Track_IP2Derr, &b_TrackIP2Derr);
  fChain->SetBranchAddress("Track_IPerr",     Track_IPerr,   &b_TrackIPerr);  
  fChain->SetBranchAddress("Track_IPsig",     Track_IPsig,   &b_TrackIPsig);
  fChain->SetBranchAddress("Track_Proba",     Track_Proba,   &b_TrackProba);
  fChain->SetBranchAddress("Track_p",         Track_p,	     &b_Trackp );   
  fChain->SetBranchAddress("Track_pt",        Track_pt,	     &b_Trackpt);
  fChain->SetBranchAddress("Track_eta",       Track_eta,     &b_Tracketa);
  fChain->SetBranchAddress("Track_phi",       Track_phi,     &b_Trackphi);
  fChain->SetBranchAddress("Track_chi2",      Track_chi2,    &b_Trackchi2);
  fChain->SetBranchAddress("Track_charge",    Track_charge,  &b_Trackcharge);
  fChain->SetBranchAddress("Track_history",   Track_history, &b_Trackhistory);
  fChain->SetBranchAddress("Track_nHitStrip", Track_nHitStrip,  &b_TracknHitStrip);
  fChain->SetBranchAddress("Track_nHitPixel", Track_nHitPixel,  &b_TracknHitPixel);
  fChain->SetBranchAddress("Track_nHitAll",   Track_nHitAll,    &b_TracknHitAll);
  fChain->SetBranchAddress("Track_nHitTIB",   Track_nHitTIB,    &b_TracknHitTIB);
  fChain->SetBranchAddress("Track_nHitTID",   Track_nHitTID,    &b_TracknHitTID);
  fChain->SetBranchAddress("Track_nHitTOB",   Track_nHitTOB,    &b_TracknHitTOB);
  fChain->SetBranchAddress("Track_nHitTEC",   Track_nHitTEC,    &b_TracknHitTEC);
  fChain->SetBranchAddress("Track_nHitPXB",   Track_nHitPXB,    &b_TracknHitPXB);
  fChain->SetBranchAddress("Track_nHitPXF",   Track_nHitPXF,    &b_TracknHitPXF);
  fChain->SetBranchAddress("Track_isHitL1",   Track_isHitL1,    &b_TrackisHitL1);
  fChain->SetBranchAddress("Track_PV",        Track_PV,         &b_TrackPV );    
  fChain->SetBranchAddress("Track_SV",        Track_SV,         &b_TrackSV);     
  fChain->SetBranchAddress("Track_PVweight",  Track_PVweight,   &b_TrackPVweight );
  fChain->SetBranchAddress("Track_SVweight",  Track_SVweight,   &b_TrackSVweight);
  fChain->SetBranchAddress("Track_category",  Track_category,   &b_Trackcategory );
  fChain->SetBranchAddress("Track_isfromSV",  Track_isfromSV,   &b_TrackisfromSV);     

  //--------------------------------------
  // primary vertex information 
  //--------------------------------------
  fChain->SetBranchAddress("PV_x"	   ,PV_x  ,  &b_PV_x );
  fChain->SetBranchAddress("PV_y"	   ,PV_y  ,  &b_PV_y );	
  fChain->SetBranchAddress("PV_z"	   ,PV_z  ,  &b_PV_z );	
  fChain->SetBranchAddress("PV_ex"	   ,PV_ex ,  &b_PV_ex );	
  fChain->SetBranchAddress("PV_ey"	   ,PV_ey ,  &b_PV_ey );	
  fChain->SetBranchAddress("PV_ez"	   ,PV_ez ,  &b_PV_ez );	
  fChain->SetBranchAddress("PV_chi2"	   ,PV_chi2 ,&b_PV_chi2 );
  fChain->SetBranchAddress("PV_ndf"	   ,PV_ndf , &b_PV_ndf );	
  fChain->SetBranchAddress("PV_isgood"    ,PV_isgood,&b_PV_isgood	 );
  fChain->SetBranchAddress("PV_isfake"    ,PV_isfake,&b_PV_isfake );
  
  
  }
   Notify();
}

Bool_t CommPlotProducer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CommPlotProducer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   cout << " Show entry = " << entry << endl;
   fChain->Show(entry);
}
Int_t CommPlotProducer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) {
 double DeltaPhi = TMath::Abs(phi2 - phi1);
 if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
 return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}


#endif // #ifdef CommPlotProducer_cxx



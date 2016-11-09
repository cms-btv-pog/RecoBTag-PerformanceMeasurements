//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 13 15:05:58 2012 by ROOT version 5.32/00
// from TTree ttree/ttree
// found on file: JetTreeTP_Pt-80to120.root
//////////////////////////////////////////////////////////

#ifndef CommPlotProducer4ttbar_h
#define CommPlotProducer4ttbar_h
#define ntrack_max 2000
#define nMaxJets_ 1000
#define nMaxLeptons_ 1000
#define nMaxTrk_ 1000

#include "TH1D.h"
#include "TH2D.h"
#include <iomanip>
#include <TROOT.h>
#include <TChain.h>
#include <string.h>
#include <TFile.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <map>
#include "../TTbarSelector.h"
 
// Header file for the classes stored in the TTree if any.


class CommPlotProducer4ttbar 
{
public :


   bool isData;
   bool use_selected_tracks;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   bool            produceJetProbaTree;
   bool            produceNewAlgoTree;
   bool            produceCTagTree;

   std::vector<TH1D*>   HistoBtag;  
   std::vector<TH1D*>   HistoTTbar;  
   std::vector<TH2D*>   HistoBtag2D;  
    
   std::map<std::string, int>   HistoBtag_map;
   std::map<std::string, int>   HistoTTbar_map; 
   std::map<std::string, int>   HistoBtag2D_map;

   
   double x_section[10]; 
     
   int numb_histo;
   int numb_histo2D;   
   int numb_histo2;
   
   
   TTbarSelector thettbarselector_;
   
   // Declaration of leaf types

   Int_t           nBitTrigger;
   Int_t           BitTrigger[3];   //[nBitTrigger]
   Int_t           nJet;
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         pthat;
   Float_t         mcweight;
   Float_t         nPUtrue;
   Int_t           nPU;
   Int_t           ttbar_chan;
   Int_t           ttbar_trigWord;
   Int_t           ttbar_metfilterWord;
   Int_t           ttbar_allmepartons;
   Int_t           ttbar_matchmepartons;
   Int_t           ttbar_nl;
   Float_t         ttbar_lpt[2];   //[ttbar_nl]
   Float_t         ttbar_leta[2];   //[ttbar_nl]
   Float_t         ttbar_lphi[2];   //[ttbar_nl]
   Float_t         ttbar_lm[2];   //[ttbar_nl]
   Int_t           ttbar_lid[2];   //[ttbar_nl]
   Int_t           ttbar_lgid[2];   //[ttbar_nl]
   Int_t           ttbar_lch[2];   //[ttbar_nl]
   Float_t         ttbar_metpt;
   Float_t         ttbar_metphi;
   Float_t         ttbar_rho;
   Int_t           ttbar_nw;
   Float_t         ttbar_w[250];   //[ttbar_nw]   
   Float_t         lepton1_pT;
   Float_t         lepton1_eta;
   Float_t         lepton1_phi;
   Float_t         lepton2_pT;
   Float_t         lepton2_eta;
   Float_t         lepton2_phi;
   Float_t         met;
   Float_t         mll;
   Int_t           trig_ttbar;

   
   
   Float_t         Jet_pt[1000];   //[nJet]
   Float_t         Jet_genpt[1000];   //[nJet]
   Float_t         Jet_residual[1000];   //[nJet]
   Float_t         Jet_jes[1000];   //[nJet]
   Float_t         Jet_eta[1000];   //[nJet]
   Float_t         Jet_phi[1000];   //[nJet]
   Float_t         Jet_mass[1000];   //[nJet]
   Int_t           Jet_ntracks[1000];   //[nJet]
   Int_t           Jet_flavour[1000];   //[nJet]
   Float_t         Jet_Ip1N[1000];   //[nJet]
   Float_t         Jet_Ip1P[1000];   //[nJet] 
   Float_t         Jet_Ip2N[1000];   //[nJet]
   Float_t         Jet_Ip2P[1000];   //[nJet]
   Float_t         Jet_Ip3N[1000];   //[nJet]
   Float_t         Jet_Ip3P[1000];   //[nJet]
   Float_t         Jet_ProbaN[1000];   //[nJet]
   Float_t         Jet_ProbaP[1000];   //[nJet]
   Float_t         Jet_BprobN[1000];   //[nJet]
   Float_t         Jet_BprobP[1000];   //[nJet]
   Float_t         Jet_SvxN[1000];   //[nJet]
   Float_t         Jet_Svx[1000];   //[nJet]
   Float_t         Jet_SvxNHP[1000];   //[nJet]
   Float_t         Jet_SvxHP[1000];   //[nJet]
   Float_t         Jet_SvxMass[1000];   //[nJet]
   //Float_t         Jet_CombSvxN[1000];   //[nJet]
   //Float_t         Jet_CombSvxP[1000];   //[nJet]
   Float_t         Jet_CombSvx[1000];   //[nJet]
   //Float_t         Jet_SimpIVF_HP[1000];   //[nJet]
   //Float_t         Jet_SimpIVF_HE[1000];   //[nJet]
   //Float_t         Jet_DoubIVF_HE[1000];   //[nJet]
   Float_t         Jet_cMVAv2[1000];   //[nJet]
   Float_t         Jet_CombIVF[1000];   //[nJet]
   Float_t         Jet_CombIVF_P[1000];   //[nJet]
   //Float_t         Jet_RetCombSvxN[1000];   //[nJet]
   //Float_t         Jet_RetCombSvxP[1000];   //[nJet]
   //Float_t         Jet_RetCombSvx[1000];   //[nJet]
   //Float_t         Jet_CombCSVJP_N[1000];   //[nJet]
   //Float_t         Jet_CombCSVJP_P[1000];   //[nJet]
   //Float_t         Jet_CombCSVJP[1000];   //[nJet]
   //Float_t         Jet_CombCSVSL_N[1000];   //[nJet]
   //Float_t         Jet_CombCSVSL_P[1000];   //[nJet]
   //Float_t         Jet_CombCSVSL[1000];   //[nJet]
   //Float_t         Jet_CombCSVJPSL_N[1000];   //[nJet]
   //Float_t         Jet_CombCSVJPSL_P[1000];   //[nJet]
   //Float_t         Jet_CombCSVJPSL[1000];   //[nJet]
   Float_t         Jet_SoftMuN[1000];   //[nJet]
   Float_t         Jet_SoftMuP[1000];   //[nJet]
   Float_t         Jet_SoftMu[1000];   //[nJet]
   Float_t         Jet_SoftElN[1000];   //[nJet]
   Float_t         Jet_SoftElP[1000];   //[nJet]
   Float_t         Jet_SoftEl[1000];   //[nJet]

   Float_t         Jet_nFirstTrkEtaRelTagVarCSV[1000];   //[nJet]
   Float_t         Jet_nLastTrkEtaRelTagVarCSV[1000];   //[nJet]
   Float_t         Jet_nFirstTrkTagVarCSV[1000];   //[nJet]
   Float_t         Jet_nLastTrkTagVarCSV[1000];   //[nJet]
   Float_t         TagVarCSV_trackSip2dSigAboveCharm[1000];   //[nJet]
   Float_t         TagVarCSV_trackSumJetEtRatio[1000];   //[nJet]
   Float_t         TagVarCSV_trackSumJetDeltaR[1000];   //[nJet]
   Float_t         TagVarCSV_trackEtaRel[1000];   //[nJet]
   Float_t         TagVarCSV_trackSip3dSig[1000];   //[nJet]
   Float_t         TagVarCSV_vertexEnergyRatio[1000];   //[nJet]
   Float_t         TagVarCSV_vertexCategory[1000];   //[nJet]
   Float_t         TagVarCSV_vertexMass[1000];   //[nJet]
   Float_t         TagVarCSV_vertexNTracks[1000];   //[nJet]
   Float_t         TagVarCSV_flightDistance2dSig[1000];   //[nJet]
   Float_t         TagVarCSV_vertexJetDeltaR[1000];   //[nJet]

   Float_t         CTag_Jet_CvsB[1000];   //[nJet] 
   Float_t         CTag_Jet_CvsL[1000];   //[nJet]

   Int_t           Jet_hist1[1000];   //[nJet]
   Int_t           Jet_hist2[1000];   //[nJet]
   Int_t           Jet_hist3[1000];   //[nJet]
   Int_t           Jet_histJet[1000];   //[nJet]
   Int_t           Jet_histSvx[1000];   //[nJet]
   Int_t           Jet_nFirstTrack[1000];   //[nJet]
   Int_t           Jet_nLastTrack[1000];   //[nJet]
   Int_t           Jet_nFirstTrkInc[1000];   //[nJet]
   Int_t           Jet_nLastTrkInc[1000];   //[nJet]
   Int_t           Jet_nFirstSV[1000];   //[nJet]
   Int_t           Jet_nLastSV[1000];   //[nJet]
   Int_t           Jet_SV_multi[1000];   //[nJet]
   
   
   Int_t           nTrkInc;
   Float_t         TrkInc_pt[1000];   //[nTrkInc]
   Float_t         TrkInc_eta[1000];   //[nTrkInc]
   Float_t         TrkInc_phi[1000];   //[nTrkInc]
   Float_t         TrkInc_ptrel[1000];   //[nTrkInc]
   Float_t         TrkInc_IPsig[1000];   //[nTrkInc]
   Float_t         TrkInc_IP[1000];   //[nTrkInc]
   
   Int_t           nCFromGSplit;
   Float_t         cFromGSplit_pT[1000];   //[nCFromGSplit]
   Float_t         cFromGSplit_eta[1000];   //[nCFromGSplit]
   Float_t         cFromGSplit_phi[1000];   //[nCFromGSplit]
   Int_t           nBFromGSplit;
   Float_t         bFromGSplit_pT[1000];   //[nBFromGSplit]
   Float_t         bFromGSplit_eta[1000];   //[nBFromGSplit]
   Float_t         bFromGSplit_phi[1000];   //[nBFromGSplit]


   
   Int_t           nBHadrons;
   Float_t         BHadron_pT[1000];   //[nBHadrons]
   Float_t         BHadron_eta[1000];   //[nBHadrons]
   Float_t         BHadron_phi[1000];   //[nBHadrons]
   Float_t         BHadron_mass[1000];   //[nBHadrons]
   Int_t           BHadron_pdgID[1000];   //[nBHadrons]

   
   Int_t           nPFElectron;
   Int_t           PFElectron_IdxJet[100];   //[nPFElectron]
   Float_t         PFElectron_pt[100];   //[nPFElectron]
   Float_t         PFElectron_eta[100];   //[nPFElectron]
   Float_t         PFElectron_phi[100];   //[nPFElectron]
   Float_t         PFElectron_ptrel[100];   //[nPFElectron]
   Float_t         PFElectron_deltaR[100];   //[nPFElectron]
   Float_t         PFElectron_ratio[100];   //[nPFElectron]
   Float_t         PFElectron_ratioRel[100];   //[nPFElectron]
   Int_t           nPFMuon;
   Int_t           PFMuon_IdxJet[100];   //[nPFMuon]
   Float_t         PFMuon_pt[100];   //[nPFMuon]
   Float_t         PFMuon_eta[100];   //[nPFMuon]
   Float_t         PFMuon_phi[100];   //[nPFMuon]
   Float_t         PFMuon_ptrel[100];	//[nPFMuon]
   Float_t         PFMuon_deltaR[100];   //[nPFMuon]
   Float_t         PFMuon_ratio[100];	//[nPFMuon]
   Float_t         PFMuon_ratioRel[100];   //[nPFMuon]
   Float_t         PFMuon_IPsig[100];   //[nPFMuon]
   Int_t           PFMuon_GoodQuality[100];   //[nPFMuon]
   
   
   
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
   
   Int_t nSV;	    
   Float_t SV_x[1000];            //[nSV] 
   Float_t SV_y[1000];              //[nSV]   
   Float_t SV_z[1000];             //[nSV] 	    
   Float_t SV_ex[1000];             //[nSV] 	    
   Float_t SV_ey[1000];             //[nSV] 	    
   Float_t SV_ez[1000];             //[nSV] 	    
   Float_t SV_chi2 [1000];              //[nSV]    
   Float_t SV_ndf[1000];             //[nSV]       
   Float_t SV_flight[1000];              //[nSV]   
   Float_t SV_flightErr[1000];             //[nSV] 
   Float_t SV_deltaR_jet[1000];             //[nSV]
   Float_t SV_deltaR_sum_jet[1000];             //[nSV]
   Float_t SV_deltaR_sum_dir[1000];             //[nSV]
   Float_t SV_EnergyRatio[1000];             //[nSV]
   Float_t SV_vtx_pt[1000];             //[nSV]
   Float_t SV_flight2D[1000];             //[nSV]
   Float_t SV_flight2DErr[1000];             //[nSV]
   Float_t SV_totCharge[1000];             //[nSV]
   Float_t SV_vtxDistJetAxis[1000];             //[nSV]
   Int_t SV_nTrk[1000];              //[nSV]
   Float_t SV_mass[1000];             //[nSV]	 
   Float_t SV_vtx_eta[1000];             //[nSV]
   Float_t SV_vtx_phi[1000];              //[nSV]  
       
   //List of branches for CTag tree
   //per jet
   Int_t Jet_nFirstTrkCTagVar[nMaxJets_];
   Int_t Jet_nLastTrkCTagVar[nMaxJets_];
   Int_t Jet_nFirstTrkEtaRelCTagVar[nMaxJets_];
   Int_t Jet_nLastTrkEtaRelCTagVar[nMaxJets_];
   Int_t Jet_nFirstLepCTagVar[nMaxJets_];
   Int_t Jet_nLastLepCTagVar[nMaxJets_]; 
   Float_t CTag_Jet_CvsBN[nMaxJets_];
   Float_t CTag_Jet_CvsBP[nMaxJets_]; 
   Float_t CTag_Jet_CvsLN[nMaxJets_];
   Float_t CTag_Jet_CvsLP[nMaxJets_];
   Float_t CTag_jetNTracks[nMaxJets_];
   Float_t CTag_jetNTracksEtaRel[nMaxJets_];
   Float_t CTag_jetNLeptons[nMaxJets_];
   Float_t CTag_trackSumJetEtRatio[nMaxJets_];
   Float_t CTag_trackSumJetDeltaR[nMaxJets_];
   Float_t CTag_trackSip2dSigAboveCharm[nMaxJets_];
   Float_t CTag_trackSip3dSigAboveCharm[nMaxJets_];
   Float_t CTag_vertexCategory[nMaxJets_];
   Float_t CTag_jetNSecondaryVertices[nMaxJets_];
   Float_t CTag_vertexMass[nMaxJets_];
   Float_t CTag_vertexNTracks[nMaxJets_];
   Float_t CTag_vertexEnergyRatio[nMaxJets_];
   Float_t CTag_vertexJetDeltaR[nMaxJets_];
   //Float_t CTag_flightDistance2dVal[nMaxJets_];
   Float_t CTag_flightDistance2dSig[nMaxJets_];
   //Float_t CTag_flightDistance3dVal[nMaxJets_];
   Float_t CTag_flightDistance3dSig[nMaxJets_];
   //Float_t CTag_vertexFitProb[nMaxJets_];
   Float_t CTag_massVertexEnergyFraction[nMaxJets_];
   Float_t CTag_vertexBoostOverSqrtJetPt[nMaxJets_];
   Float_t CTag_vertexLeptonCategory[nMaxJets_];
   //per jet per track
   Int_t nTrkCTagVar;
   Int_t nTrkEtaRelCTagVar;
   Float_t CTag_trackPtRel[nMaxTrk_];                            
   Float_t CTag_trackPPar[nMaxTrk_];                            
   Float_t CTag_trackDeltaR[nMaxTrk_];                           
   Float_t CTag_trackPtRatio[nMaxTrk_];                          
   Float_t CTag_trackPParRatio[nMaxTrk_];                        
   //Float_t CTag_trackSip2dVal[nMaxTrk_];                         
   Float_t CTag_trackSip2dSig[nMaxTrk_];                         
   //Float_t CTag_trackSip3dVal[nMaxTrk_];                         
   Float_t CTag_trackSip3dSig[nMaxTrk_];                         
   Float_t CTag_trackDecayLenVal[nMaxTrk_];                      
   Float_t CTag_trackJetDistVal[nMaxTrk_];                       
   Float_t CTag_trackEtaRel[nMaxTrk_]; 
   //per jet per lepton
   Int_t   nLeptons;
   Float_t CTag_leptonPtRel[nMaxLeptons_];
   Float_t CTag_leptonSip3d[nMaxLeptons_];
   Float_t CTag_leptonDeltaR[nMaxLeptons_];
   Float_t CTag_leptonRatioRel[nMaxLeptons_];
   Float_t CTag_leptonEtaRel[nMaxLeptons_];
   Float_t CTag_leptonRatio[nMaxLeptons_];
   
   
   // List of branches
   TBranch        *b_nBitTrigger;   //!
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_ttbar_chan;   //!
   TBranch        *b_ttbar_trigWord;   //!
   TBranch        *b_ttbar_metfilterWord;   //!
   TBranch        *b_ttbar_allmepartons;   //!
   TBranch        *b_ttbar_matchmepartons;   //!
   TBranch        *b_ttbar_nl;   //!
   TBranch        *b_ttbar_lpt;   //!
   TBranch        *b_ttbar_leta;   //!
   TBranch        *b_ttbar_lphi;   //!
   TBranch        *b_ttbar_lm;   //!
   TBranch        *b_ttbar_lid;   //!
   TBranch        *b_ttbar_lgid;   //!
   TBranch        *b_ttbar_lch;   //!
   TBranch        *b_ttbar_metpt;   //!
   TBranch        *b_ttbar_metphi;   //!
   TBranch        *b_ttbar_rho;   //!
   TBranch        *b_ttbar_nw;   //!
   TBranch        *b_ttbar_w;   //!   
   //TBranch        *b_lepton1_pT;   //!
   //TBranch        *b_lepton1_eta;   //!
   //TBranch        *b_lepton1_phi;   //!
   //TBranch        *b_lepton2_pT;   //!
   //TBranch        *b_lepton2_eta;   //!
   //TBranch        *b_lepton2_phi;   //!
   //TBranch        *b_met;   //!
   //TBranch        *b_mll;   //!
   //TBranch        *b_trig_ttbar;   //!

   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_genpt;   //!
   TBranch        *b_Jet_residual;   //!
   TBranch        *b_Jet_jes;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_ntracks;   //!
   TBranch        *b_Jet_flavour;   //!
   TBranch        *b_Jet_Ip1N;   //!
   TBranch        *b_Jet_Ip1P;   //!
   TBranch        *b_Jet_Ip2N;   //!
   TBranch        *b_Jet_Ip2P;   //!
   TBranch        *b_Jet_Ip3N;   //!
   TBranch        *b_Jet_Ip3P;   //!
   TBranch        *b_Jet_ProbaN;   //!
   TBranch        *b_Jet_ProbaP;   //!
   TBranch        *b_Jet_BprobN;   //!
   TBranch        *b_Jet_BprobP;   //!
   TBranch        *b_Jet_SvxN;   //!
   TBranch        *b_Jet_Svx;   //!
   TBranch        *b_Jet_SvxNHP;   //!
   TBranch        *b_Jet_SvxHP;   //!
   TBranch        *b_Jet_SvxMass;   //!
   //TBranch        *b_Jet_CombSvxP;   //!
   TBranch        *b_Jet_CombSvx;   //!
   //TBranch        *b_Jet_SimpIVF_HP;   //!
   //TBranch        *b_Jet_SimpIVF_HE;   //!
   //TBranch        *b_Jet_DoubIVF_HE;   //!
   TBranch        *b_Jet_cMVAv2;   //!
   TBranch        *b_Jet_CombIVF;   //!
   TBranch        *b_Jet_CombIVF_P;   //!
   //TBranch        *b_Jet_RetCombSvxN;   //!
   //TBranch        *b_Jet_RetCombSvxP;   //!
   //TBranch        *b_Jet_RetCombSvx;   //!
   //TBranch        *b_Jet_CombCSVJP_N;   //!
   //TBranch        *b_Jet_CombCSVJP_P;   //!
   //TBranch        *b_Jet_CombCSVJP;   //!
   //TBranch        *b_Jet_CombCSVSL_N;   //!
   //TBranch        *b_Jet_CombCSVSL_P;   //!
   //TBranch        *b_Jet_CombCSVSL;   //!
   //TBranch        *b_Jet_CombCSVJPSL_N;   //!
   //TBranch        *b_Jet_CombCSVJPSL_P;   //!
   //TBranch        *b_Jet_CombCSVJPSL;   //!
   TBranch        *b_Jet_SoftMuN;   //!
   TBranch        *b_Jet_SoftMuP;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_SoftElN;   //!
   TBranch        *b_Jet_SoftElP;   //!
   TBranch        *b_Jet_SoftEl;   //!

   TBranch        *b_Jet_nFirstTrkEtaRelTagVarCSV;
   TBranch        *b_Jet_nLastTrkEtaRelTagVarCSV;         
   TBranch        *b_Jet_nFirstTrkTagVarCSV;              
   TBranch        *b_Jet_nLastTrkTagVarCSV;               
   TBranch        *b_TagVarCSV_trackSip2dSigAboveCharm;   
   TBranch        *b_TagVarCSV_trackSumJetEtRatio;        
   TBranch        *b_TagVarCSV_trackSumJetDeltaR;         
   TBranch        *b_TagVarCSV_trackEtaRel;               
   TBranch        *b_TagVarCSV_trackSip3dSig;             
   TBranch        *b_TagVarCSV_vertexEnergyRatio;         
   TBranch        *b_TagVarCSV_vertexCategory;            
   TBranch        *b_TagVarCSV_vertexMass;                
   TBranch        *b_TagVarCSV_vertexNTracks;             
   TBranch        *b_TagVarCSV_flightDistance2dSig;       
   TBranch        *b_TagVarCSV_vertexJetDeltaR;           

   TBranch        *b_CTag_Jet_CvsB;           
   TBranch        *b_CTag_Jet_CvsL;  

   TBranch        *b_Jet_hist1;   //!
   TBranch        *b_Jet_hist2;   //!
   TBranch        *b_Jet_hist3;   //!
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //!
   TBranch        *b_Jet_nFirstTrack;   //!
   TBranch        *b_Jet_nLastTrack;   //!
   TBranch        *b_Jet_nFirstTrkInc;   //!
   TBranch        *b_Jet_nLastTrkInc;   //!
   TBranch        *b_Jet_nFirstSV;
   TBranch        *b_Jet_nLastSV;
   TBranch        *b_Jet_SV_multi;
   TBranch        *b_nTrkInc;   //!
   TBranch        *b_TrkInc_pt;   //!
   TBranch        *b_TrkInc_eta;   //!
   TBranch        *b_TrkInc_phi;   //!
   TBranch        *b_TrkInc_ptrel;   //!
   TBranch        *b_TrkInc_IPsig;   //!
   TBranch        *b_TrkInc_IP;   //!
   TBranch        *b_nCFromGSplit;   //!
   TBranch        *b_cFromGSplit_pT;   //!
   TBranch        *b_cFromGSplit_eta;   //!
   TBranch        *b_cFromGSplit_phi;   //!
   TBranch        *b_nBFromGSplit;   //!
   TBranch        *b_bFromGSplit_pT;   //!
   TBranch        *b_bFromGSplit_eta;   //!
   TBranch        *b_bFromGSplit_phi;   //!
   TBranch        *b_nBHadrons;   //!
   TBranch        *b_BHadron_pT;   //!
   TBranch        *b_BHadron_eta;   //!
   TBranch        *b_BHadron_phi;   //!
   TBranch        *b_BHadron_mass;   //!
   TBranch        *b_BHadron_pdgID;   //!


   TBranch        *b_nPFElectron;   //!
   TBranch        *b_PFElectron_IdxJet;   //!
   TBranch        *b_PFElectron_pt;   //!
   TBranch        *b_PFElectron_eta;   //!
   TBranch        *b_PFElectron_phi;   //!
   TBranch        *b_PFElectron_ptrel;   //!
   TBranch        *b_PFElectron_deltaR;   //!
   TBranch        *b_PFElectron_ratio;   //!
   TBranch        *b_PFElectron_ratioRel;   //!
   TBranch        *b_nPFMuon;   //!
   TBranch        *b_PFMuon_IdxJet;   //!
   TBranch        *b_PFMuon_pt;   //!
   TBranch        *b_PFMuon_eta;   //!
   TBranch        *b_PFMuon_phi;   //!
   TBranch        *b_PFMuon_ptrel;   //!
   TBranch        *b_PFMuon_deltaR;   //!
   TBranch        *b_PFMuon_ratio;   //!
   TBranch        *b_PFMuon_ratioRel;   //!
   TBranch        *b_PFMuon_IPsig;   //!
   TBranch        *b_PFMuon_GoodQuality;   //!   
   
   //Lists of CTag brahch
   TBranch        *c_Jet_nFirstTrkCTagVar;
   TBranch        *c_Jet_nLastTrkCTagVar;
   TBranch        *c_Jet_nFirstTrkEtaRelCTagVar;
   TBranch        *c_Jet_nLastTrkEtaRelCTagVar;
   TBranch        *c_Jet_nFirstLepCTagVar;
   TBranch        *c_Jet_nLastLepCTagVar;
   TBranch        *c_CTag_Jet_CvsB;
   TBranch        *c_CTag_Jet_CvsBN;
   TBranch        *c_CTag_Jet_CvsBP;
   TBranch        *c_CTag_Jet_CvsL;
   TBranch        *c_CTag_Jet_CvsLN;
   TBranch        *c_CTag_Jet_CvsLP;
   TBranch        *c_CTag_jetNTracks;
   TBranch        *c_CTag_jetNTracksEtaRel;
   TBranch        *c_CTag_jetNLeptons;
   TBranch        *c_CTag_trackSumJetEtRatio;
   TBranch        *c_CTag_trackSumJetDeltaR;
   TBranch        *c_CTag_trackSip2dSigAboveCharm;
   TBranch        *c_CTag_trackSip3dSigAboveCharm;
   TBranch        *c_CTag_vertexCategory;
   TBranch        *c_CTag_jetNSecondaryVertices;
   TBranch        *c_CTag_vertexMass;
   TBranch        *c_CTag_vertexNTracks;
   TBranch        *c_CTag_vertexEnergyRatio;
   TBranch        *c_CTag_vertexJetDeltaR;
   //TBranch        *c_CTag_flightDistance2dVal;
   TBranch        *c_CTag_flightDistance2dSig;
   //TBranch        *c_CTag_flightDistance3dVal;
   TBranch        *c_CTag_flightDistance3dSig;
   //TBranch        *c_CTag_vertexFitProb;
   TBranch        *c_CTag_massVertexEnergyFraction;
   TBranch        *c_CTag_vertexBoostOverSqrtJetPt;
   TBranch        *c_CTag_vertexLeptonCategory;
   //track information
   TBranch        *c_nTrkCTagVar;
   TBranch        *c_nTrkEtaRelCTagVar;
   TBranch        *c_CTag_trackPtRel;
   TBranch        *c_CTag_trackPPar;
   TBranch        *c_CTag_trackDeltaR;
   TBranch        *c_CTag_trackPtRatio;
   TBranch        *c_CTag_trackPParRatio;
   //TBranch        *c_CTag_trackSip2dVal;
   TBranch        *c_CTag_trackSip2dSig;
   //TBranch        *c_CTag_trackSip3dVal;
   TBranch        *c_CTag_trackSip3dSig;
   TBranch        *c_CTag_trackDecayLenVal;
   TBranch        *c_CTag_trackJetDistVal;
   TBranch        *c_CTag_trackEtaRel;
   //lepton information
   TBranch        *c_nLeptons;
   TBranch        *c_CTag_leptonPtRel;
   TBranch        *c_CTag_leptonSip3d;
   TBranch        *c_CTag_leptonDeltaR;
   TBranch        *c_CTag_leptonRatioRel;
   TBranch        *c_CTag_leptonEtaRel;
   TBranch        *c_CTag_leptonRatio;

   //--------------------------------------
   // track information 
   //--------------------------------------
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
  
  //--------------------------------------
  // secondary vertex information 
  //--------------------------------------
  TBranch *b_nSV;
  TBranch *b_SV_x;
  TBranch *b_SV_y;
  TBranch *b_SV_z;
  TBranch *b_SV_ex;
  TBranch *b_SV_ey;
  TBranch *b_SV_ez;
  TBranch *b_SV_chi2;
  TBranch *b_SV_ndf;
  TBranch *b_SV_flight;
  TBranch *b_SV_flightErr;
  TBranch *b_SV_deltaR_jet;             
  TBranch *b_SV_deltaR_sum_jet;  
  TBranch *b_SV_deltaR_sum_dir;  
  TBranch *b_SV_EnergyRatio;    
  TBranch *b_SV_vtx_pt;          
  TBranch *b_SV_flight2D;        
  TBranch *b_SV_flight2DErr;     
  TBranch *b_SV_totCharge;       
  TBranch *b_SV_vtxDistJetAxis;  
  TBranch *b_SV_nTrk;  
  TBranch *b_SV_mass;  
  TBranch *b_SV_vtx_eta;
  TBranch *b_SV_vtx_phi;
  
   CommPlotProducer4ttbar(TChain *supertree=0, bool infotree1=true, bool infotree2=true, bool infotreeCTag=false);
   virtual ~CommPlotProducer4ttbar();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(int datatype, int trig_data, float PtMin_Cut, float PtMax_Cut, TString outputname, TH1F* wgtcounter, TString syst);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     AddHisto(TString name, TString title,  int nbins, float min, float max, TString syst);
   virtual void     AddHistottbar(TString name, TString title,  int nbins, float min, float max, TString syst);
   virtual void     AddHisto2D(TString name,TString title,int nbins,float min,float max,int nbins2,float min2,float max2, TString syst);   
   virtual void     FillHisto_int(int flavour, bool isGS, int number, int value, double weight);
   virtual void     FillHisto_float(int flavour, bool isGS, int number, float value, double weight);
   virtual void     FillHisto_floatFromMap(TString name, int flavour, bool isGS, float value, double weight);
   virtual void     FillHistottbar_intFromMap(TString name, int flavour, bool isGS, int value, double weight);
   virtual void     FillHistottbar_floatFromMap(TString name, int flavour, bool isGS, float value, double weight);
   virtual void     FillHisto_intFromMap(TString name, int flavour, bool isGS, int value, double weight);
   virtual void     FillHisto2D_int_floatFromMap(TString name, int flavour, bool isGS, int value, float value2, double weight);
   virtual void     FillHisto2D_float_floatFromMap(TString name, int flavour, bool isGS, float value, float value2, double weight);
   virtual void     SetNorm(float xnorm);

//private :
   TGraph *puWgtGr_,*puWgtDownGr_,*puWgtUpGr_;

};
#endif

#ifdef CommPlotProducer4ttbar_cxx
CommPlotProducer4ttbar::CommPlotProducer4ttbar(TChain *superTree, bool infotree1, bool infotree2, bool infotreeCTag)
{  

   numb_histo = 0;
   numb_histo2 = 0;
   numb_histo2D = 0;

   produceJetProbaTree=infotree1;
   produceNewAlgoTree=infotree2;
   produceCTagTree=infotreeCTag;

   if (produceJetProbaTree) use_selected_tracks=false;
   else use_selected_tracks=true;

   
   Init(superTree);
}

CommPlotProducer4ttbar::~CommPlotProducer4ttbar()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CommPlotProducer4ttbar::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CommPlotProducer4ttbar::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);

   if (centry < 0) return centry;

   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }

   return centry;
}

void CommPlotProducer4ttbar::Init(TChain *tree)
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

   fChain->SetBranchAddress("BitTrigger", &BitTrigger, &b_BitTrigger);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("ttbar_chan", &ttbar_chan, &b_ttbar_chan);
   fChain->SetBranchAddress("ttbar_trigWord", &ttbar_trigWord, &b_ttbar_trigWord);
   fChain->SetBranchAddress("ttbar_metfilterWord", &ttbar_metfilterWord, &b_ttbar_metfilterWord);
   fChain->SetBranchAddress("ttbar_allmepartons", &ttbar_allmepartons, &b_ttbar_allmepartons);
   fChain->SetBranchAddress("ttbar_matchmepartons", &ttbar_matchmepartons, &b_ttbar_matchmepartons);
   fChain->SetBranchAddress("ttbar_nl", &ttbar_nl, &b_ttbar_nl);
   fChain->SetBranchAddress("ttbar_lpt", ttbar_lpt, &b_ttbar_lpt);
   fChain->SetBranchAddress("ttbar_leta", ttbar_leta, &b_ttbar_leta);
   fChain->SetBranchAddress("ttbar_lphi", ttbar_lphi, &b_ttbar_lphi);
   fChain->SetBranchAddress("ttbar_lm", ttbar_lm, &b_ttbar_lm);
   fChain->SetBranchAddress("ttbar_lid", ttbar_lid, &b_ttbar_lid);
   fChain->SetBranchAddress("ttbar_lgid", ttbar_lgid, &b_ttbar_lgid);
   fChain->SetBranchAddress("ttbar_lch", ttbar_lch, &b_ttbar_lch);
   fChain->SetBranchAddress("ttbar_metpt", &ttbar_metpt, &b_ttbar_metpt);
   fChain->SetBranchAddress("ttbar_metphi", &ttbar_metphi, &b_ttbar_metphi);
   fChain->SetBranchAddress("ttbar_rho", &ttbar_rho, &b_ttbar_rho);
   fChain->SetBranchAddress("ttbar_nw", &ttbar_nw, &b_ttbar_nw);
   fChain->SetBranchAddress("ttbar_w", ttbar_w, &b_ttbar_w);

   //fChain->SetBranchAddress("lepton1_pT", &lepton1_pT, &b_lepton1_pT);
   //fChain->SetBranchAddress("lepton1_eta", &lepton1_eta, &b_lepton1_eta);
   //fChain->SetBranchAddress("lepton1_phi", &lepton1_phi, &b_lepton1_phi);
   //fChain->SetBranchAddress("lepton2_pT", &lepton2_pT, &b_lepton2_pT);
   //fChain->SetBranchAddress("lepton2_eta", &lepton2_eta, &b_lepton2_eta);
   //fChain->SetBranchAddress("lepton2_phi", &lepton2_phi, &b_lepton2_phi);
   //fChain->SetBranchAddress("met", &met, &b_met);
   //fChain->SetBranchAddress("mll", &mll, &b_mll);
   //fChain->SetBranchAddress("trig_ttbar", &trig_ttbar, &b_trig_ttbar);

   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_genpt", Jet_genpt, &b_Jet_genpt);
   fChain->SetBranchAddress("Jet_residual", Jet_residual, &b_Jet_residual);
   fChain->SetBranchAddress("Jet_jes", Jet_jes, &b_Jet_jes);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_ntracks", Jet_ntracks, &b_Jet_ntracks);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   //fChain->SetBranchAddress("Jet_Ip1N", Jet_Ip1N, &b_Jet_Ip1N);
   //fChain->SetBranchAddress("Jet_Ip1P", Jet_Ip1P, &b_Jet_Ip1P);
   fChain->SetBranchAddress("Jet_Ip2N", Jet_Ip2N, &b_Jet_Ip2N);
   fChain->SetBranchAddress("Jet_Ip2P", Jet_Ip2P, &b_Jet_Ip2P);
   fChain->SetBranchAddress("Jet_Ip3N", Jet_Ip3N, &b_Jet_Ip3N);
   fChain->SetBranchAddress("Jet_Ip3P", Jet_Ip3P, &b_Jet_Ip3P);
   fChain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN, &b_Jet_ProbaN);
   fChain->SetBranchAddress("Jet_ProbaP", Jet_ProbaP, &b_Jet_ProbaP);
   fChain->SetBranchAddress("Jet_BprobN", Jet_BprobN, &b_Jet_BprobN);
   fChain->SetBranchAddress("Jet_BprobP", Jet_BprobP, &b_Jet_BprobP);
   fChain->SetBranchAddress("Jet_SvxN", Jet_SvxN, &b_Jet_SvxN);
   fChain->SetBranchAddress("Jet_Svx", Jet_Svx, &b_Jet_Svx);
   fChain->SetBranchAddress("Jet_SvxNHP", Jet_SvxNHP, &b_Jet_SvxNHP);
   fChain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP, &b_Jet_SvxHP);
   //fChain->SetBranchAddress("Jet_SvxMass", Jet_SvxMass, &b_Jet_SvxMass);
   //fChain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN, &b_Jet_CombSvxN);
   //fChain->SetBranchAddress("Jet_CombSvxP", Jet_CombSvxP, &b_Jet_CombSvxP);
   fChain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx, &b_Jet_CombSvx);
   //fChain->SetBranchAddress("Jet_SimpIVF_HP", Jet_SimpIVF_HP, &b_Jet_SimpIVF_HP);
   //fChain->SetBranchAddress("Jet_SimpIVF_HE", Jet_SimpIVF_HE, &b_Jet_SimpIVF_HE);
   //fChain->SetBranchAddress("Jet_DoubIVF_HE", Jet_DoubIVF_HE, &b_Jet_DoubIVF_HE);
   fChain->SetBranchAddress("Jet_cMVAv2", Jet_cMVAv2, &b_Jet_cMVAv2);
   fChain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF, &b_Jet_CombIVF);
   fChain->SetBranchAddress("Jet_CombIVF_P", Jet_CombIVF_P, &b_Jet_CombIVF_P);
     fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
     fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   if(!produceCTagTree){
   fChain->SetBranchAddress("Jet_nFirstTrack", Jet_nFirstTrack, &b_Jet_nFirstTrack);
   fChain->SetBranchAddress("Jet_nLastTrack", Jet_nLastTrack, &b_Jet_nLastTrack);
   //fChain->SetBranchAddress("Jet_nFirstTrkInc", Jet_nFirstTrkInc, &b_Jet_nFirstTrkInc);
   //fChain->SetBranchAddress("Jet_nLastTrkInc", Jet_nLastTrkInc, &b_Jet_nLastTrkInc);
     fChain->SetBranchAddress("Jet_nFirstSV", Jet_nFirstSV, &b_Jet_nFirstSV);
     fChain->SetBranchAddress("Jet_nLastSV", Jet_nLastSV, &b_Jet_nLastSV);
   }
   fChain->SetBranchAddress("Jet_SV_multi",Jet_SV_multi , &b_Jet_SV_multi);
   //fChain->SetBranchAddress("nTrkInc", &nTrkInc, &b_nTrkInc);
   //fChain->SetBranchAddress("TrkInc_pt", TrkInc_pt, &b_TrkInc_pt);
   //fChain->SetBranchAddress("TrkInc_eta", TrkInc_eta, &b_TrkInc_eta);
   //fChain->SetBranchAddress("TrkInc_ptrel", TrkInc_ptrel, &b_TrkInc_ptrel);
   //fChain->SetBranchAddress("TrkInc_IPsig", TrkInc_IPsig, &b_TrkInc_IPsig);
   //fChain->SetBranchAddress("TrkInc_IP", TrkInc_IP, &b_TrkInc_IP);
   //fChain->SetBranchAddress("nCFromGSplit", &nCFromGSplit, &b_nCFromGSplit);
   //fChain->SetBranchAddress("cFromGSplit_pT", cFromGSplit_pT, &b_cFromGSplit_pT);
   //fChain->SetBranchAddress("cFromGSplit_eta", cFromGSplit_eta, &b_cFromGSplit_eta);
   //fChain->SetBranchAddress("cFromGSplit_phi", cFromGSplit_phi, &b_cFromGSplit_phi);
   //fChain->SetBranchAddress("nBFromGSplit", &nBFromGSplit, &b_nBFromGSplit);
   //fChain->SetBranchAddress("bFromGSplit_pT", bFromGSplit_pT, &b_bFromGSplit_pT);
   //fChain->SetBranchAddress("bFromGSplit_eta", bFromGSplit_eta, &b_bFromGSplit_eta);
   //fChain->SetBranchAddress("bFromGSplit_phi", bFromGSplit_phi, &b_bFromGSplit_phi);   
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);

    if(produceCTagTree){
   //--------------------------------------
   //   CTag information 
   //--------------------------------------
   fChain->SetBranchAddress("Jet_nFirstTrkCTagVar"        ,Jet_nFirstTrkCTagVar        ,&c_Jet_nFirstTrkCTagVar);
   fChain->SetBranchAddress("Jet_nLastTrkCTagVar"         ,Jet_nLastTrkCTagVar         ,&c_Jet_nLastTrkCTagVar);
   fChain->SetBranchAddress("Jet_nFirstTrkEtaRelCTagVar"        ,Jet_nFirstTrkEtaRelCTagVar        ,&c_Jet_nFirstTrkEtaRelCTagVar);
   fChain->SetBranchAddress("Jet_nLastTrkEtaRelCTagVar"        ,Jet_nLastTrkEtaRelCTagVar        ,&c_Jet_nLastTrkEtaRelCTagVar);
   fChain->SetBranchAddress("Jet_nFirstLepCTagVar"        ,Jet_nFirstLepCTagVar        ,&c_Jet_nFirstLepCTagVar);
   fChain->SetBranchAddress("Jet_nLastLepCTagVar"         ,Jet_nLastLepCTagVar         ,&c_Jet_nLastLepCTagVar);   
   fChain->SetBranchAddress("CTag_Jet_CvsB", CTag_Jet_CvsB, &c_CTag_Jet_CvsB);
   fChain->SetBranchAddress("CTag_Jet_CvsBN", CTag_Jet_CvsBN, &c_CTag_Jet_CvsBN);
   fChain->SetBranchAddress("CTag_Jet_CvsBP", CTag_Jet_CvsBP, &c_CTag_Jet_CvsBP);
   fChain->SetBranchAddress("CTag_Jet_CvsL", CTag_Jet_CvsL, &c_CTag_Jet_CvsL);
   fChain->SetBranchAddress("CTag_Jet_CvsLN", CTag_Jet_CvsLN, &c_CTag_Jet_CvsLN);   
   fChain->SetBranchAddress("CTag_Jet_CvsLP", CTag_Jet_CvsLP, &c_CTag_Jet_CvsLP);
   fChain->SetBranchAddress("CTag_jetNTracks"               ,CTag_jetNTracks               ,&c_CTag_jetNTracks);
   fChain->SetBranchAddress("CTag_jetNTracksEtaRel"         ,CTag_jetNTracksEtaRel         ,&c_CTag_jetNTracksEtaRel);
   fChain->SetBranchAddress("CTag_trackSumJetEtRatio"       ,CTag_trackSumJetEtRatio       ,&c_CTag_trackSumJetEtRatio);
   fChain->SetBranchAddress("CTag_trackSumJetDeltaR"        ,CTag_trackSumJetDeltaR        ,&c_CTag_trackSumJetDeltaR);
   fChain->SetBranchAddress("CTag_trackSip2dSigAboveCharm"  ,CTag_trackSip2dSigAboveCharm  ,&c_CTag_trackSip2dSigAboveCharm);
   fChain->SetBranchAddress("CTag_trackSip3dSigAboveCharm"  ,CTag_trackSip3dSigAboveCharm  ,&c_CTag_trackSip3dSigAboveCharm);
   fChain->SetBranchAddress("CTag_vertexCategory"           ,CTag_vertexCategory           ,&c_CTag_vertexCategory);    
   fChain->SetBranchAddress("CTag_jetNSecondaryVertices"    ,CTag_jetNSecondaryVertices    ,&c_CTag_jetNSecondaryVertices);
   fChain->SetBranchAddress("CTag_vertexMass"               ,CTag_vertexMass               ,&c_CTag_vertexMass);
   fChain->SetBranchAddress("CTag_vertexNTracks"            ,CTag_vertexNTracks            ,&c_CTag_vertexNTracks);
   fChain->SetBranchAddress("CTag_vertexEnergyRatio"        ,CTag_vertexEnergyRatio        ,&c_CTag_vertexEnergyRatio);
   fChain->SetBranchAddress("CTag_vertexJetDeltaR"          ,CTag_vertexJetDeltaR          ,&c_CTag_vertexJetDeltaR);
   //fChain->SetBranchAddress("CTag_flightDistance2dVal"      ,CTag_flightDistance2dVal      ,&c_CTag_flightDistance2dVal);
   fChain->SetBranchAddress("CTag_flightDistance2dSig"      ,CTag_flightDistance2dSig      ,&c_CTag_flightDistance2dSig);
   //fChain->SetBranchAddress("CTag_flightDistance3dVal"      ,CTag_flightDistance3dVal      ,&c_CTag_flightDistance3dVal);
   fChain->SetBranchAddress("CTag_flightDistance3dSig"      ,CTag_flightDistance3dSig      ,&c_CTag_flightDistance3dSig);  
   //fChain->SetBranchAddress("CTag_vertexFitProb"            ,CTag_vertexFitProb         ,&c_CTag_vertexFitProb);
   fChain->SetBranchAddress("CTag_massVertexEnergyFraction" ,CTag_massVertexEnergyFraction     ,&c_CTag_massVertexEnergyFraction);
   fChain->SetBranchAddress("CTag_vertexBoostOverSqrtJetPt" ,CTag_vertexBoostOverSqrtJetPt     ,&c_CTag_vertexBoostOverSqrtJetPt);
   fChain->SetBranchAddress("CTag_vertexLeptonCategory"     ,CTag_vertexLeptonCategory     ,&c_CTag_vertexLeptonCategory);
   fChain->SetBranchAddress("nTrkCTagVar"               ,&nTrkCTagVar              ,&c_nTrkCTagVar);
   fChain->SetBranchAddress("nTrkEtaRelCTagVar"         ,&nTrkEtaRelCTagVar        ,&c_nTrkEtaRelCTagVar);
   fChain->SetBranchAddress("CTag_trackPtRel"        ,CTag_trackPtRel        ,&c_CTag_trackPtRel);
   fChain->SetBranchAddress("CTag_trackPPar"         ,CTag_trackPPar         ,&c_CTag_trackPPar);
   fChain->SetBranchAddress("CTag_trackDeltaR"       ,CTag_trackDeltaR       ,&c_CTag_trackDeltaR);
   fChain->SetBranchAddress("CTag_trackPtRatio"      ,CTag_trackPtRatio      ,&c_CTag_trackPtRatio);
   fChain->SetBranchAddress("CTag_trackPParRatio"    ,CTag_trackPParRatio    ,&c_CTag_trackPParRatio);
   //fChain->SetBranchAddress("CTag_trackSip2dVal"     ,CTag_trackSip2dVal     ,&c_CTag_trackSip2dVal);
   fChain->SetBranchAddress("CTag_trackSip2dSig"     ,CTag_trackSip2dSig     ,&c_CTag_trackSip2dSig);
   //fChain->SetBranchAddress("CTag_trackSip3dVal"     ,CTag_trackSip3dVal     ,&c_CTag_trackSip3dVal);
   fChain->SetBranchAddress("CTag_trackSip3dSig"     ,CTag_trackSip3dSig     ,&c_CTag_trackSip3dSig);
   fChain->SetBranchAddress("CTag_trackDecayLenVal"  ,CTag_trackDecayLenVal  ,&c_CTag_trackDecayLenVal);
   fChain->SetBranchAddress("CTag_trackJetDistVal"   ,CTag_trackJetDistVal   ,&c_CTag_trackJetDistVal);
   fChain->SetBranchAddress("CTag_trackEtaRel"       ,CTag_trackEtaRel       ,&c_CTag_trackEtaRel);
   fChain->SetBranchAddress("CTag_jetNLeptons"               ,CTag_jetNLeptons               ,&c_CTag_jetNLeptons);  
   fChain->SetBranchAddress("nLeptons"               ,&nLeptons              ,&c_nLeptons);
   fChain->SetBranchAddress("CTag_leptonPtRel"         ,CTag_leptonPtRel         ,&c_CTag_leptonPtRel);
   fChain->SetBranchAddress("CTag_leptonSip3d"         ,CTag_leptonSip3d         ,&c_CTag_leptonSip3d);
   fChain->SetBranchAddress("CTag_leptonDeltaR"         ,CTag_leptonDeltaR         ,&c_CTag_leptonDeltaR);
   fChain->SetBranchAddress("CTag_leptonRatioRel"         ,CTag_leptonRatioRel         ,&c_CTag_leptonRatioRel);
   fChain->SetBranchAddress("CTag_leptonEtaRel"         ,CTag_leptonEtaRel         ,&c_CTag_leptonEtaRel);
   fChain->SetBranchAddress("CTag_leptonRatio"         ,CTag_leptonRatio         ,&c_CTag_leptonRatio);

   fChain->SetBranchAddress("Jet_nFirstTrkEtaRelTagVarCSV",      Jet_nFirstTrkEtaRelTagVarCSV,        &b_Jet_nFirstTrkEtaRelTagVarCSV);
   fChain->SetBranchAddress("Jet_nLastTrkEtaRelTagVarCSV",       Jet_nLastTrkEtaRelTagVarCSV,         &b_Jet_nLastTrkEtaRelTagVarCSV);
   fChain->SetBranchAddress("Jet_nFirstTrkTagVarCSV",            Jet_nFirstTrkTagVarCSV,              &b_Jet_nFirstTrkTagVarCSV);
   fChain->SetBranchAddress("Jet_nLastTrkTagVarCSV",             Jet_nLastTrkTagVarCSV,               &b_Jet_nLastTrkTagVarCSV);
   fChain->SetBranchAddress("TagVarCSV_trackSip2dSigAboveCharm", TagVarCSV_trackSip2dSigAboveCharm,   &b_TagVarCSV_trackSip2dSigAboveCharm);
   fChain->SetBranchAddress("TagVarCSV_trackSumJetEtRatio",      TagVarCSV_trackSumJetEtRatio,        &b_TagVarCSV_trackSumJetEtRatio);
   fChain->SetBranchAddress("TagVarCSV_trackSumJetDeltaR",       TagVarCSV_trackSumJetDeltaR,         &b_TagVarCSV_trackSumJetDeltaR);
   fChain->SetBranchAddress("TagVarCSV_trackEtaRel",             TagVarCSV_trackEtaRel,               &b_TagVarCSV_trackEtaRel);
   fChain->SetBranchAddress("TagVarCSV_trackSip3dSig",           TagVarCSV_trackSip3dSig,             &b_TagVarCSV_trackSip3dSig);
   fChain->SetBranchAddress("TagVarCSV_vertexEnergyRatio",       TagVarCSV_vertexEnergyRatio,         &b_TagVarCSV_vertexEnergyRatio);
   fChain->SetBranchAddress("TagVarCSV_vertexCategory",          TagVarCSV_vertexCategory,            &b_TagVarCSV_vertexCategory);
   fChain->SetBranchAddress("TagVarCSV_vertexMass",              TagVarCSV_vertexMass,                &b_TagVarCSV_vertexMass);
   fChain->SetBranchAddress("TagVarCSV_vertexNTracks",           TagVarCSV_vertexNTracks,             &b_TagVarCSV_vertexNTracks);
   fChain->SetBranchAddress("TagVarCSV_flightDistance2dSig",     TagVarCSV_flightDistance2dSig,       &b_TagVarCSV_flightDistance2dSig);
   fChain->SetBranchAddress("TagVarCSV_vertexJetDeltaR",         TagVarCSV_vertexJetDeltaR,           &b_TagVarCSV_vertexJetDeltaR);
   } 
   
   
  if ( produceJetProbaTree ) {
      
   //--------------------------------------
   // track information 
   //--------------------------------------
  fChain->SetBranchAddress("nTrack",	      &nTrack, 	     &b_nTrack);
  fChain->SetBranchAddress("Track_dxy",       Track_dxy,     &b_Trackdxy);	
  fChain->SetBranchAddress("Track_dz",        Track_dz, &b_Track_dz);
  //fChain->SetBranchAddress("Track_zIP",       Track_zIP, &b_Track_zIP);
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
  
  //--------------------------------------
  // secondary vertex information 
  //--------------------------------------
  fChain->SetBranchAddress("nSV"	   ,&nSV	  ,&b_nSV );
  fChain->SetBranchAddress("SV_x"	   ,SV_x  ,&b_SV_x );	
  fChain->SetBranchAddress("SV_y"	   ,SV_y  ,&b_SV_y );	
  fChain->SetBranchAddress("SV_z"	   ,SV_z  ,&b_SV_z );	
  fChain->SetBranchAddress("SV_ex"	   ,SV_ex ,&b_SV_ex );	
  fChain->SetBranchAddress("SV_ey"	   ,SV_ey ,&b_SV_ey );	
  fChain->SetBranchAddress("SV_ez"	   ,SV_ez ,&b_SV_ez );	
  fChain->SetBranchAddress("SV_chi2"	   ,SV_chi2 ,&b_SV_chi2	 );
  fChain->SetBranchAddress("SV_ndf"	   ,SV_ndf  ,&b_SV_ndf );	
  fChain->SetBranchAddress("SV_flight"    ,SV_flight,&b_SV_flight	 );
  fChain->SetBranchAddress("SV_flightErr" ,SV_flightErr,&b_SV_flightErr );
  fChain->SetBranchAddress("SV_deltaR_jet" ,SV_deltaR_jet,&b_SV_deltaR_jet);
  fChain->SetBranchAddress("SV_deltaR_sum_jet" ,SV_deltaR_sum_jet,&b_SV_deltaR_sum_jet);
  fChain->SetBranchAddress("SV_deltaR_sum_dir" ,SV_deltaR_sum_dir,&b_SV_deltaR_sum_dir);
  fChain->SetBranchAddress("SV_EnergyRatio" ,SV_EnergyRatio,&b_SV_EnergyRatio);
  fChain->SetBranchAddress("SV_vtx_pt" ,SV_vtx_pt,&b_SV_vtx_pt);
  fChain->SetBranchAddress("SV_flight2D" ,SV_flight2D,&b_SV_flight2D);
  fChain->SetBranchAddress("SV_flight2DErr" ,SV_flight2DErr,&b_SV_flight2DErr);
  fChain->SetBranchAddress("SV_totCharge" ,SV_totCharge,&b_SV_totCharge);
  fChain->SetBranchAddress("SV_vtxDistJetAxis" ,SV_vtxDistJetAxis,&b_SV_vtxDistJetAxis);
  fChain->SetBranchAddress("SV_nTrk" ,SV_nTrk,&b_SV_nTrk);
  fChain->SetBranchAddress("SV_mass" ,SV_mass,&b_SV_mass);
  fChain->SetBranchAddress("SV_vtx_eta" ,SV_vtx_eta,&b_SV_vtx_eta);
  fChain->SetBranchAddress("SV_vtx_phi" ,SV_vtx_phi,&b_SV_vtx_phi);

  
  }
  if (produceNewAlgoTree) {
   //fChain->SetBranchAddress("Jet_RetCombSvxN", Jet_RetCombSvxN, &b_Jet_RetCombSvxN);
   //fChain->SetBranchAddress("Jet_RetCombSvxP", Jet_RetCombSvxP, &b_Jet_RetCombSvxP);
   //fChain->SetBranchAddress("Jet_RetCombSvx", Jet_RetCombSvx, &b_Jet_RetCombSvx);
   //fChain->SetBranchAddress("Jet_CombCSVJP_N", Jet_CombCSVJP_N, &b_Jet_CombCSVJP_N);
   //fChain->SetBranchAddress("Jet_CombCSVJP_P", Jet_CombCSVJP_P, &b_Jet_CombCSVJP_P);
   //fChain->SetBranchAddress("Jet_CombCSVJP", Jet_CombCSVJP, &b_Jet_CombCSVJP);
   //fChain->SetBranchAddress("Jet_CombCSVSL_N", Jet_CombCSVSL_N, &b_Jet_CombCSVSL_N);
   //fChain->SetBranchAddress("Jet_CombCSVSL_P", Jet_CombCSVSL_P, &b_Jet_CombCSVSL_P);
   //fChain->SetBranchAddress("Jet_CombCSVSL", Jet_CombCSVSL, &b_Jet_CombCSVSL);
   //fChain->SetBranchAddress("Jet_CombCSVJPSL_N", Jet_CombCSVJPSL_N, &b_Jet_CombCSVJPSL_N);
   //fChain->SetBranchAddress("Jet_CombCSVJPSL_P", Jet_CombCSVJPSL_P, &b_Jet_CombCSVJPSL_P);
   //fChain->SetBranchAddress("Jet_CombCSVJPSL", Jet_CombCSVJPSL, &b_Jet_CombCSVJPSL);

   fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
   fChain->SetBranchAddress("Jet_SoftMuP", Jet_SoftMuP, &b_Jet_SoftMuP);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_SoftElN", Jet_SoftElN, &b_Jet_SoftElN);
   fChain->SetBranchAddress("Jet_SoftElP", Jet_SoftElP, &b_Jet_SoftElP);
   fChain->SetBranchAddress("Jet_SoftEl", Jet_SoftEl, &b_Jet_SoftEl);
   
   fChain->SetBranchAddress("Jet_nFirstTrkEtaRelTagVarCSV",      Jet_nFirstTrkEtaRelTagVarCSV,        &b_Jet_nFirstTrkEtaRelTagVarCSV);
   fChain->SetBranchAddress("Jet_nLastTrkEtaRelTagVarCSV",       Jet_nLastTrkEtaRelTagVarCSV,         &b_Jet_nLastTrkEtaRelTagVarCSV);         
   fChain->SetBranchAddress("Jet_nFirstTrkTagVarCSV",            Jet_nFirstTrkTagVarCSV,              &b_Jet_nFirstTrkTagVarCSV);              
   fChain->SetBranchAddress("Jet_nLastTrkTagVarCSV",             Jet_nLastTrkTagVarCSV,               &b_Jet_nLastTrkTagVarCSV);               
   fChain->SetBranchAddress("TagVarCSV_trackSip2dSigAboveCharm", TagVarCSV_trackSip2dSigAboveCharm,   &b_TagVarCSV_trackSip2dSigAboveCharm);   
   fChain->SetBranchAddress("TagVarCSV_trackSumJetEtRatio",      TagVarCSV_trackSumJetEtRatio,        &b_TagVarCSV_trackSumJetEtRatio);        
   fChain->SetBranchAddress("TagVarCSV_trackSumJetDeltaR",       TagVarCSV_trackSumJetDeltaR,         &b_TagVarCSV_trackSumJetDeltaR);         
   fChain->SetBranchAddress("TagVarCSV_trackEtaRel",             TagVarCSV_trackEtaRel,               &b_TagVarCSV_trackEtaRel);               
   fChain->SetBranchAddress("TagVarCSV_trackSip3dSig",           TagVarCSV_trackSip3dSig,             &b_TagVarCSV_trackSip3dSig);             
   fChain->SetBranchAddress("TagVarCSV_vertexEnergyRatio",       TagVarCSV_vertexEnergyRatio,         &b_TagVarCSV_vertexEnergyRatio);         
   fChain->SetBranchAddress("TagVarCSV_vertexCategory",          TagVarCSV_vertexCategory,            &b_TagVarCSV_vertexCategory);            
   fChain->SetBranchAddress("TagVarCSV_vertexMass",              TagVarCSV_vertexMass,                &b_TagVarCSV_vertexMass);                
   fChain->SetBranchAddress("TagVarCSV_vertexNTracks",           TagVarCSV_vertexNTracks,             &b_TagVarCSV_vertexNTracks);             
   fChain->SetBranchAddress("TagVarCSV_flightDistance2dSig",     TagVarCSV_flightDistance2dSig,       &b_TagVarCSV_flightDistance2dSig);       
   fChain->SetBranchAddress("TagVarCSV_vertexJetDeltaR",         TagVarCSV_vertexJetDeltaR,           &b_TagVarCSV_vertexJetDeltaR);           

   fChain->SetBranchAddress("CTag_Jet_CvsB",CTag_Jet_CvsB , &b_CTag_Jet_CvsB);
   fChain->SetBranchAddress("CTag_Jet_CvsL",CTag_Jet_CvsL , &b_CTag_Jet_CvsL);
   
   fChain->SetBranchAddress("nPFElectron", &nPFElectron, &b_nPFElectron);
   fChain->SetBranchAddress("PFElectron_IdxJet", PFElectron_IdxJet, &b_PFElectron_IdxJet);
   fChain->SetBranchAddress("PFElectron_pt", PFElectron_pt, &b_PFElectron_pt);
   fChain->SetBranchAddress("PFElectron_eta", PFElectron_eta, &b_PFElectron_eta);
   fChain->SetBranchAddress("PFElectron_phi", PFElectron_phi, &b_PFElectron_phi);
   fChain->SetBranchAddress("PFElectron_ptrel", PFElectron_ptrel, &b_PFElectron_ptrel);
   fChain->SetBranchAddress("PFElectron_deltaR", PFElectron_deltaR, &b_PFElectron_deltaR);
   fChain->SetBranchAddress("PFElectron_ratio", PFElectron_ratio, &b_PFElectron_ratio);
   fChain->SetBranchAddress("PFElectron_ratioRel", PFElectron_ratioRel, &b_PFElectron_ratioRel);
   fChain->SetBranchAddress("nPFMuon", &nPFMuon, &b_nPFMuon);
   fChain->SetBranchAddress("PFMuon_IdxJet", PFMuon_IdxJet, &b_PFMuon_IdxJet);
   fChain->SetBranchAddress("PFMuon_pt", PFMuon_pt, &b_PFMuon_pt);
   fChain->SetBranchAddress("PFMuon_eta", PFMuon_eta, &b_PFMuon_eta);
   fChain->SetBranchAddress("PFMuon_phi", PFMuon_phi, &b_PFMuon_phi);
   fChain->SetBranchAddress("PFMuon_ptrel", PFMuon_ptrel, &b_PFMuon_ptrel);
   fChain->SetBranchAddress("PFMuon_deltaR", PFMuon_deltaR, &b_PFMuon_deltaR);
   fChain->SetBranchAddress("PFMuon_ratio", PFMuon_ratio, &b_PFMuon_ratio);
   fChain->SetBranchAddress("PFMuon_ratioRel", PFMuon_ratioRel, &b_PFMuon_ratioRel);
   fChain->SetBranchAddress("PFMuon_IPsig", PFMuon_IPsig, &b_PFMuon_IPsig);
   fChain->SetBranchAddress("PFMuon_GoodQuality", PFMuon_GoodQuality, &b_PFMuon_GoodQuality);   
  }


   //pileup weights
   puWgtGr_ = 0, puWgtDownGr_ = 0, puWgtUpGr_ = 0;
   
   TString puWgtUrl("${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test_ttbarSelector_Moriond16/BTagAnalyzerMacros/pileupWgts.root");
   gSystem->ExpandPathName(puWgtUrl);
   TFile *fIn=TFile::Open(puWgtUrl);
   if(fIn)
   {
           puWgtGr_     = (TGraph *)fIn->Get("puwgts_nom");
           puWgtDownGr_ = (TGraph *)fIn->Get("puwgts_down");
           puWgtUpGr_   = (TGraph *)fIn->Get("puwgts_up");
           fIn->Close();
   }
   else
   {
           std::cout << "Unable to find pileupWgts.root, no PU reweighting will be applied" << std::endl;
   }
   

   Notify();
}

Bool_t CommPlotProducer4ttbar::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CommPlotProducer4ttbar::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#endif // #ifdef CommPlotProducer4ttbar_cxx

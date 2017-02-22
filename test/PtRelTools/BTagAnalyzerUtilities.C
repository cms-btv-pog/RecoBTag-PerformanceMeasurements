#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"

#include <fstream>
#include <iostream>

#include "BTagAnalyzerUtilities/Ntuples/RunIISpring16MiniAODv2-80X.h"
#include "BTagAnalyzerUtilities/Taggers/RunIISpring16MiniAODv2-80X.h"

int Run, Evt, nPU, nPV, nTrkInc, nBHadrons, nDHadrons, nGenlep, nGenquark, nBitTrigger;
int BitTrigger[100], Jet_flavour[50], Jet_nFirstTrkInc[50], Jet_nLastTrkInc[50], Muon_IdxJet[50]; 
float pthat, PVz, nPUtrue, mcweight;
int LumiBlock;
int nFatJet;
float FatJet_pt[50], FatJet_eta[50], FatJet_phi[50], FatJet_mass[50], FatJet_residual[50], FatJet_tau1[50], FatJet_tau2[50];
int FatJet_looseID[50], FatJet_nSubJets[50], FatJet_FirstSubJet[50], FatJet_LastSubJet[50];
int nJet, Jet_looseID[50], Jet_nbHadrons[50];
float Jet_pt[50], Jet_genpt[50], Jet_eta[50], Jet_phi[50], Jet_Ip1P[50], Jet_Ip2P[50], Jet_Ip3P[50], Jet_Ip3N[50], Jet_jes[50];
float Jet_Svx[50], Jet_SvxHP[50], Jet_Proba[50], Jet_ProbaP[50], Jet_ProbaN[50], Jet_Bprob[50], Jet_CombSvxP[50], Jet_CombSvxN[50], Jet_CombSvx[50];
float Jet_SimpIVF_HE[50], Jet_SimpIVF_HP[50], Jet_CombIVF[50], Jet_CombIVF_P[50], Jet_cMVAv2[50], Jet_residual[50];
float Jet_RetCombSvxP[50], Jet_RetCombSvxN[50], Jet_RetCombSvx[50], Jet_CombCSVJP_P[50], Jet_CombCSVJP_N[50], Jet_CombCSVJP[50];
float Jet_CombCSVJPSL_P[50], Jet_CombCSVJPSL_N[50], Jet_CombCSVJPSL[50], Jet_CombCSVSL_P[50], Jet_CombCSVSL_N[50], Jet_CombCSVSL[50];
float Jet_SoftMu[50], Jet_SoftMuP[50], Jet_SoftMuN[50], Jet_SoftEl[50], Jet_SoftElP[50], Jet_SoftElN[50];
int Jet_VtxCat[50], Jet_nseltracks[50];
int nMuon, nPFMuon;
float Muon_pt[50], Muon_eta[50], Muon_phi[50], Muon_vz[50];
float PFMuon_pt[50], PFMuon_eta[50], PFMuon_phi[50], PFMuon_ptrel[50], PFMuon_IP[50], PFMuon_IPsig[50]; 
int PFMuon_GoodQuality[50], PFMuon_IdxJet[50];
float Muon_ptrel[50], Muon_IP[50], Muon_IPsig[50], Muon_Proba[50], Muon_chi2[50], Muon_chi2Tk[50]; 
int Muon_nMuHit[50], Muon_nTkHit[50], Muon_nPixHit[50], Muon_nOutHit[50], Muon_isGlobal[50], Muon_nMatched[50];
float TrkInc_pt[200], TrkInc_eta[200], TrkInc_phi[200], TrkInc_ptrel[200], TrkInc_IP[200], TrkInc_IPsig[200];
float BHadron_pT[100], BHadron_eta[100], BHadron_phi[100], BHadron_mass[100]; int BHadron_pdgID[100], BHadron_hasBdaughter[100], BHadron_DHadron2[100], BHadron_DHadron1[100];
float DHadron_pT[100], DHadron_eta[100], DHadron_phi[100], DHadron_mass[100]; int DHadron_pdgID[100];
float Genlep_pT[20], Genlep_eta[20], Genlep_phi[20]; int Genlep_pdgID[20], Genlep_mother[20];
float Genquark_pT[20], Genquark_eta[20], Genquark_phi[20]; int Genquark_pdgID[20], Genquark_mother[20];
int ncQuarks, nbQuarks;
float cQuark_eta[50], cQuark_phi[50], bQuark_eta[50], bQuark_phi[50];
int cQuark_status[50], cQuark_pdgID[50], bQuark_status[50], bQuark_pdgID[50], cQuark_fromGSP[50], bQuark_fromGSP[50];
int nStdJet, StdJet_flavour[50];
float StdJet_pt[50], StdJet_eta[50], StdJet_phi[50], StdJet_Proba[50], StdJet_Bprob[50], StdJet_CombIVF[50], StdJet_cMVAv2[50];

TTree *GetChain(TFile *ThisTree, bool ReadLightTracks = false, TString TreeName = "btagana", bool ReadOnlyEventInfo = false) {

  TTree *tchain = (TTree*)ThisTree->Get(TreeName + "/ttree"); 

  if (TreeName=="btagana") {

    tchain->SetBranchAddress("Run", &Run);
    tchain->SetBranchAddress("Evt", &Evt);
    tchain->SetBranchAddress("nBitTrigger", &nBitTrigger);
    tchain->SetBranchAddress("BitTrigger", BitTrigger);
    tchain->SetBranchAddress("LumiBlock", &LumiBlock);
    tchain->SetBranchAddress("nPU", &nPU);
    tchain->SetBranchAddress("nPUtrue", &nPUtrue);
    tchain->SetBranchAddress("nPV", &nPV);
    tchain->SetBranchAddress("pthat", &pthat);
    tchain->SetBranchAddress("mcweight", &mcweight);
    tchain->SetBranchAddress("PVz", &PVz);

  }

  if (TreeName=="btaganaFatJets") {

    TString SubjetType = "SoftDropSubJetInfo";
    
    tchain->SetBranchAddress("FatJetInfo.nJet", &nFatJet);
    tchain->SetBranchAddress("FatJetInfo.Jet_pt", FatJet_pt);
    tchain->SetBranchAddress("FatJetInfo.Jet_eta", FatJet_eta);
    tchain->SetBranchAddress("FatJetInfo.Jet_phi", FatJet_phi);
    tchain->SetBranchAddress("FatJetInfo.Jet_mass", FatJet_mass);
    tchain->SetBranchAddress("FatJetInfo.Jet_looseID", FatJet_looseID);
    tchain->SetBranchAddress("FatJetInfo.Jet_residual", FatJet_residual);
    tchain->SetBranchAddress("FatJetInfo.Jet_tau2", FatJet_tau2);
    tchain->SetBranchAddress("FatJetInfo.Jet_tau1", FatJet_tau1);
    if (SubjetType=="PrunedSubJetInfo") {
      tchain->SetBranchAddress("FatJetInfo.Jet_nSubJets_Pruned", FatJet_nSubJets);
      tchain->SetBranchAddress("FatJetInfo.Jet_nFirstSJ_Pruned", FatJet_FirstSubJet);
      tchain->SetBranchAddress("FatJetInfo.Jet_nLastSJ_Pruned", FatJet_LastSubJet);
    } else if (SubjetType=="SoftDropSubJetInfo") {
      tchain->SetBranchAddress("FatJetInfo.Jet_nSubJets_SoftDrop", FatJet_nSubJets);
      tchain->SetBranchAddress("FatJetInfo.Jet_nFirstSJ_SoftDrop", FatJet_FirstSubJet);
      tchain->SetBranchAddress("FatJetInfo.Jet_nLastSJ_SoftDrop", FatJet_LastSubJet);
    }
    tchain->SetBranchAddress(SubjetType + ".nJet", &nJet);
    tchain->SetBranchAddress(SubjetType + ".Jet_pt", Jet_pt);
    tchain->SetBranchAddress(SubjetType + ".Jet_genpt", Jet_genpt);
    tchain->SetBranchAddress(SubjetType + ".Jet_eta", Jet_eta);
    tchain->SetBranchAddress(SubjetType + ".Jet_phi", Jet_phi);
    tchain->SetBranchAddress(SubjetType + ".Jet_flavour", Jet_flavour);
    //tchain->SetBranchAddress(SubjetType + ".Jet_Ip1P", Jet_Ip1P);
    tchain->SetBranchAddress(SubjetType + ".Jet_Ip2P", Jet_Ip2P);
    tchain->SetBranchAddress(SubjetType + ".Jet_Ip3P", Jet_Ip3P);
    tchain->SetBranchAddress(SubjetType + ".Jet_Ip3N", Jet_Ip3N);
    tchain->SetBranchAddress(SubjetType + ".Jet_Svx", Jet_Svx);
    tchain->SetBranchAddress(SubjetType + ".Jet_SvxHP", Jet_SvxHP);
    //tchain->SetBranchAddress(SubjetType + ".Jet_Proba",  Jet_Proba);
    tchain->SetBranchAddress(SubjetType + ".Jet_ProbaP", Jet_Proba);
    tchain->SetBranchAddress(SubjetType + ".Jet_ProbaN", Jet_ProbaN);
    tchain->SetBranchAddress(SubjetType + ".Jet_BprobP", Jet_Bprob);
    tchain->SetBranchAddress(SubjetType + ".Jet_CombSvxP", Jet_CombSvxP);
    tchain->SetBranchAddress(SubjetType + ".Jet_CombSvxN", Jet_CombSvxN);
    tchain->SetBranchAddress(SubjetType + ".Jet_CombSvx", Jet_CombSvx);
    tchain->SetBranchAddress(SubjetType + ".Jet_SoftMuP", Jet_SoftMuP);
    tchain->SetBranchAddress(SubjetType + ".Jet_SoftMuN", Jet_SoftMuN);
    tchain->SetBranchAddress(SubjetType + ".Jet_SoftMu", Jet_SoftMu);
    tchain->SetBranchAddress(SubjetType + ".Jet_SoftElP", Jet_SoftElP);
    tchain->SetBranchAddress(SubjetType + ".Jet_SoftElN", Jet_SoftElN);
    tchain->SetBranchAddress(SubjetType + ".Jet_SoftEl", Jet_SoftEl);
    tchain->SetBranchAddress(SubjetType + ".Jet_CombIVF", Jet_CombIVF);
    tchain->SetBranchAddress(SubjetType + ".Jet_CombIVF_P", Jet_CombIVF_P);
    tchain->SetBranchAddress(SubjetType + ".Jet_cMVAv2", Jet_cMVAv2);
    tchain->SetBranchAddress(SubjetType + ".Jet_nseltracks", Jet_nseltracks);
    if (ReadLightTracks) tchain->SetBranchAddress("Jet_nFirstTrkInc", Jet_nFirstTrkInc);
    if (ReadLightTracks) tchain->SetBranchAddress("Jet_nLastTrkInc", Jet_nLastTrkInc);
    tchain->SetBranchAddress(SubjetType + ".nPFMuon", &nPFMuon);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_GoodQuality", PFMuon_GoodQuality);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_pt", PFMuon_pt);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_eta", PFMuon_eta);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_phi", PFMuon_phi);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_ptrel", PFMuon_ptrel);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_IP", PFMuon_IP);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_IPsig", PFMuon_IPsig);
    tchain->SetBranchAddress(SubjetType + ".PFMuon_IdxJet", PFMuon_IdxJet);
    if (ReadLightTracks) {
      tchain->SetBranchAddress(SubjetType + ".nTrkInc", &nTrkInc);
      tchain->SetBranchAddress(SubjetType + ".TrkInc_pt", TrkInc_pt);
      tchain->SetBranchAddress(SubjetType + ".TrkInc_eta", TrkInc_eta);
      tchain->SetBranchAddress(SubjetType + ".TrkInc_phi", TrkInc_phi);
      tchain->SetBranchAddress(SubjetType + ".TrkInc_ptrel", TrkInc_ptrel);
      tchain->SetBranchAddress(SubjetType + ".TrkInc_IP", TrkInc_IP);
      tchain->SetBranchAddress(SubjetType + ".TrkInc_IPsig", TrkInc_IPsig);
    }
   
  } else if (TreeName=="btagana") {

    if (!ReadOnlyEventInfo) {

      /*tchain->SetBranchAddress("Run", &Run);
	tchain->SetBranchAddress("nBitTrigger", &nBitTrigger);
	tchain->SetBranchAddress("BitTrigger", BitTrigger);
	tchain->SetBranchAddress("nPU", &nPU);
	tchain->SetBranchAddress("nPUtrue", &nPUtrue);
	tchain->SetBranchAddress("nPV", &nPV);
	tchain->SetBranchAddress("pthat", &pthat);
	tchain->SetBranchAddress("PVz", &PVz);*/
      tchain->SetBranchAddress("nJet", &nJet);
      tchain->SetBranchAddress("Jet_pt", Jet_pt);
      tchain->SetBranchAddress("Jet_jes", Jet_jes);
      tchain->SetBranchAddress("Jet_genpt", Jet_genpt);
      tchain->SetBranchAddress("Jet_eta", Jet_eta);
      tchain->SetBranchAddress("Jet_phi", Jet_phi);
      tchain->SetBranchAddress("Jet_flavour", Jet_flavour);
      tchain->SetBranchAddress("Jet_looseID", Jet_looseID);
      tchain->SetBranchAddress("Jet_nbHadrons", Jet_nbHadrons);
      //tchain->SetBranchAddress("Jet_Ip1P", Jet_Ip1P);
      tchain->SetBranchAddress("Jet_Ip2P", Jet_Ip2P);
      tchain->SetBranchAddress("Jet_Ip3P", Jet_Ip3P);
      tchain->SetBranchAddress("Jet_Ip3N", Jet_Ip3N);
      tchain->SetBranchAddress("Jet_Svx", Jet_Svx);
      tchain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP);
      //tchain->SetBranchAddress("Jet_Proba",  Jet_Proba);
      tchain->SetBranchAddress("Jet_ProbaP", Jet_Proba);
      tchain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN);
      tchain->SetBranchAddress("Jet_BprobP", Jet_Bprob);
      tchain->SetBranchAddress("Jet_CombSvxP", Jet_CombSvxP);
      tchain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN);
      tchain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx);
      /*tchain->SetBranchAddress("Jet_RetCombSvxP", Jet_RetCombSvxP);
	tchain->SetBranchAddress("Jet_RetCombSvxN", Jet_RetCombSvxN);
	tchain->SetBranchAddress("Jet_RetCombSvx", Jet_RetCombSvx);
	tchain->SetBranchAddress("Jet_CombCSVJP_P", Jet_CombCSVJP_P);
	tchain->SetBranchAddress("Jet_CombCSVJP_N", Jet_CombCSVJP_N);
	tchain->SetBranchAddress("Jet_CombCSVJP", Jet_CombCSVJP);
	tchain->SetBranchAddress("Jet_CombCSVJPSL_P", Jet_CombCSVJPSL_P);
	tchain->SetBranchAddress("Jet_CombCSVJPSL_N", Jet_CombCSVJPSL_N);
	tchain->SetBranchAddress("Jet_CombCSVJPSL", Jet_CombCSVJPSL);
	tchain->SetBranchAddress("Jet_CombCSVSL_P", Jet_CombCSVSL_P);
	tchain->SetBranchAddress("Jet_CombCSVSL_N", Jet_CombCSVSL_N);
	tchain->SetBranchAddress("Jet_CombCSVSL", Jet_CombCSVSL);*/
      tchain->SetBranchAddress("Jet_SoftMuP", Jet_SoftMuP);
      tchain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN);
      tchain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu);
      tchain->SetBranchAddress("Jet_SoftElP", Jet_SoftElP);
      tchain->SetBranchAddress("Jet_SoftElN", Jet_SoftElN);
      tchain->SetBranchAddress("Jet_SoftEl", Jet_SoftEl);
      //tchain->SetBranchAddress("Jet_SimpIVF_HP", Jet_SimpIVF_HP);
      //tchain->SetBranchAddress("Jet_SimpIVF_HE", Jet_SimpIVF_HE);
      tchain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF);
      tchain->SetBranchAddress("Jet_CombIVF_P", Jet_CombIVF_P);
      tchain->SetBranchAddress("Jet_cMVAv2", Jet_cMVAv2);
      //tchain->SetBranchAddress("Jet_VtxCat", Jet_VtxCat);
      if (ReadLightTracks) tchain->SetBranchAddress("Jet_nFirstTrkInc", Jet_nFirstTrkInc);
      if (ReadLightTracks) tchain->SetBranchAddress("Jet_nLastTrkInc", Jet_nLastTrkInc);
      tchain->SetBranchAddress("nPFMuon", &nPFMuon);
      tchain->SetBranchAddress("PFMuon_GoodQuality", PFMuon_GoodQuality);
      tchain->SetBranchAddress("PFMuon_pt", PFMuon_pt);
      tchain->SetBranchAddress("PFMuon_eta", PFMuon_eta);
      tchain->SetBranchAddress("PFMuon_phi", PFMuon_phi);
      tchain->SetBranchAddress("PFMuon_ptrel", PFMuon_ptrel);
      tchain->SetBranchAddress("PFMuon_IP", PFMuon_IP);
      //tchain->SetBranchAddress("PFMuon_IPsig", PFMuon_IPsig);
      tchain->SetBranchAddress("PFMuon_IdxJet", PFMuon_IdxJet);
      /*tchain->SetBranchAddress("nMuon", &nMuon);
	tchain->SetBranchAddress("Muon_IdxJet", Muon_IdxJet);
	tchain->SetBranchAddress("Muon_pt", Muon_pt);
	tchain->SetBranchAddress("Muon_eta", Muon_eta);
	tchain->SetBranchAddress("Muon_phi", Muon_phi);
	tchain->SetBranchAddress("Muon_vz", Muon_vz);
	tchain->SetBranchAddress("Muon_ptrel", Muon_ptrel);
	tchain->SetBranchAddress("Muon_IP", Muon_IP);
	tchain->SetBranchAddress("Muon_IPsig", Muon_IPsig);
	tchain->SetBranchAddress("Muon_Proba", Muon_Proba);
	tchain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal);
	tchain->SetBranchAddress("Muon_nMatched", Muon_nMatched);
	tchain->SetBranchAddress("Muon_nTkHit", Muon_nTkHit);
	tchain->SetBranchAddress("Muon_nPixHit", Muon_nPixHit);
	tchain->SetBranchAddress("Muon_nMuHit", Muon_nMuHit);
	tchain->SetBranchAddress("Muon_nOutHit", Muon_nOutHit);
	tchain->SetBranchAddress("Muon_chi2", Muon_chi2);
	tchain->SetBranchAddress("Muon_chi2Tk", Muon_chi2Tk);
	tchain->SetBranchAddress("nBFromGSplit", &nBFromGSplit);
	tchain->SetBranchAddress("bFromGSplit_pT", bFromGSplit_pT);
	tchain->SetBranchAddress("bFromGSplit_eta", bFromGSplit_eta);
	tchain->SetBranchAddress("bFromGSplit_phi", bFromGSplit_phi);
	tchain->SetBranchAddress("nCFromGSplit", &nCFromGSplit);
	tchain->SetBranchAddress("cFromGSplit_pT", cFromGSplit_pT);
	tchain->SetBranchAddress("cFromGSplit_eta", cFromGSplit_eta);
	tchain->SetBranchAddress("cFromGSplit_phi", cFromGSplit_phi);*/
      if (ReadLightTracks) {
	tchain->SetBranchAddress("nTrkInc", &nTrkInc);
	tchain->SetBranchAddress("TrkInc_pt", TrkInc_pt);
	tchain->SetBranchAddress("TrkInc_eta", TrkInc_eta);
	tchain->SetBranchAddress("TrkInc_phi", TrkInc_phi);
	tchain->SetBranchAddress("TrkInc_ptrel", TrkInc_ptrel);
	tchain->SetBranchAddress("TrkInc_IP", TrkInc_IP);
	tchain->SetBranchAddress("TrkInc_IPsig", TrkInc_IPsig);
      }

    } else {

      tchain->SetBranchAddress("nJet", &nStdJet);
      tchain->SetBranchAddress("Jet_pt",StdJet_pt);
      tchain->SetBranchAddress("Jet_eta",StdJet_eta);
      tchain->SetBranchAddress("Jet_phi",StdJet_phi);
      tchain->SetBranchAddress("Jet_flavour",StdJet_flavour);
      tchain->SetBranchAddress("Jet_ProbaP",StdJet_Proba);
      tchain->SetBranchAddress("Jet_BprobP",StdJet_Bprob);
      tchain->SetBranchAddress("Jet_CombIVF", StdJet_CombIVF);
      tchain->SetBranchAddress("Jet_cMVAv2", StdJet_cMVAv2);

    }

  }

  if (TreeName=="btagana") {

    tchain->SetBranchAddress("nBHadrons", &nBHadrons);
    tchain->SetBranchAddress("BHadron_pT", BHadron_pT);
    tchain->SetBranchAddress("BHadron_eta", BHadron_eta);
    tchain->SetBranchAddress("BHadron_phi", BHadron_phi);
    tchain->SetBranchAddress("BHadron_mass", BHadron_mass);
    tchain->SetBranchAddress("BHadron_pdgID", BHadron_pdgID);
    tchain->SetBranchAddress("BHadron_hasBdaughter", BHadron_hasBdaughter);
    tchain->SetBranchAddress("BHadron_DHadron1", BHadron_DHadron1);
    tchain->SetBranchAddress("BHadron_DHadron2", BHadron_DHadron2);
    tchain->SetBranchAddress("nDHadrons", &nDHadrons);
    tchain->SetBranchAddress("DHadron_pT", DHadron_pT);
    tchain->SetBranchAddress("DHadron_eta", DHadron_eta);
    tchain->SetBranchAddress("DHadron_phi", DHadron_phi);
    tchain->SetBranchAddress("DHadron_mass", DHadron_mass);
    tchain->SetBranchAddress("DHadron_pdgID", DHadron_pdgID);
    tchain->SetBranchAddress("nGenlep", &nGenlep);
    tchain->SetBranchAddress("Genlep_pT", Genlep_pT);
    tchain->SetBranchAddress("Genlep_phi", Genlep_phi);
    tchain->SetBranchAddress("Genlep_eta", Genlep_eta);
    tchain->SetBranchAddress("Genlep_pdgID", Genlep_pdgID);
    tchain->SetBranchAddress("Genlep_mother", Genlep_mother);
    tchain->SetBranchAddress("nGenquark", &nGenquark);
    tchain->SetBranchAddress("Genquark_pT", Genquark_pT);
    tchain->SetBranchAddress("Genquark_phi", Genquark_phi);
    tchain->SetBranchAddress("Genquark_eta", Genquark_eta);
    tchain->SetBranchAddress("Genquark_pdgID", Genquark_pdgID);
    tchain->SetBranchAddress("Genquark_mother", Genquark_mother);
    tchain->SetBranchAddress("ncQuarks", &ncQuarks);
    tchain->SetBranchAddress("cQuark_eta", cQuark_eta);
    tchain->SetBranchAddress("cQuark_phi", cQuark_phi);
    tchain->SetBranchAddress("cQuark_status", cQuark_status);
    tchain->SetBranchAddress("cQuark_pdgID", cQuark_pdgID);
    tchain->SetBranchAddress("cQuark_fromGSP", cQuark_fromGSP);
    tchain->SetBranchAddress("nbQuarks", &nbQuarks);
    tchain->SetBranchAddress("bQuark_eta", bQuark_eta);
    tchain->SetBranchAddress("bQuark_phi", bQuark_phi);
    tchain->SetBranchAddress("bQuark_status", bQuark_status);
    tchain->SetBranchAddress("bQuark_pdgID", bQuark_pdgID);
    tchain->SetBranchAddress("bQuark_fromGSP", bQuark_fromGSP);

  }

  return tchain;
 
}

int GetNumberOfDataRanges(TString DataType) {

  int nDataRanges = nBTagMuRanges;
  if (DataType=="JetHT") nDataRanges = nJetRunRanges;
  if (DataType=="QCDMu") nDataRanges = nMonteCarloPtHatRanges;
  if (DataType=="QCD") nDataRanges = nMCInclusivePtHatRanges;
  //if (DataType=="TTP") nDataRanges = nTTbarPowhegRanges;
  //if (DataType=="TTM") nDataRanges = nTTbarMadGraphRanges;
  //if (DataType=="Hbb") nDataRanges = nGluGluHbbRanges;

  return nDataRanges;

}

TString GetDataRangeName(TString DataType, int DataRange) {
  
  TString DataRangeName = BTagMuRangeName[DataRange];
  if (DataType=="QCDMu") DataRangeName = MonteCarloPtHatRange[DataRange];
  else if (DataType=="JetHT") DataRangeName = JetRunRangeName[DataRange];
  else if (DataType=="QCD") DataRangeName = MCInclusivePtHatRange[DataRange];
  //else if (DataType=="TTP") DataRangeName = TTbarPowhegRange[DataRange];
  //else if (DataType=="TTM") DataRangeName = TTbarMadGraphRange[DataRange];
  //else if (DataType=="Hbb") DataRangeName = GluGluHbbRange[DataRange];

  return DataRangeName;

}

float GetPtHatWeight(TString DataType, int PtHatRange) {

  if (DataType=="QCDMu") return CrossSection[PtHatRange]/GeneratedEvents[PtHatRange];
  if (DataType=="QCD") return CrossSectionInclusive[PtHatRange]/GeneratedEventsInclusive[PtHatRange];
  return 1.;

}

int GetNumberOfTrees(TString DataType, int DataRange) {

  int nTrees = nBTagMuTrees[DataRange];
  if (DataType=="QCDMu") nTrees = nMonteCarloTrees[DataRange];
  else if (DataType=="JetHT") nTrees = nJetTrees[DataRange]; 
  else if (DataType=="QCD") nTrees = nMCInclusiveTrees[DataRange];
  //else if (DataType=="TTP") nTrees = nTTbarPowhegTrees[DataRange];
  //else if (DataType=="TTM") nTrees = nTTbarMadGraphTrees[DataRange];
  //else if (DataType=="Hbb") nTrees = nGluGluHbbTrees[DataRange];
  
  return nTrees;

}

void GetDataRangeInformation(TString DataType, int DataRange, TString *DataRangeName, float *PtHatWeight, TString *FileDirectoryName, int *NumberOfTrees, int *FirstTree) {

  *DataRangeName = GetDataRangeName(DataType, DataRange);

  *PtHatWeight = GetPtHatWeight(DataType, DataRange);

  TString DataRangeDirectoryName = *DataRangeName;
  if (DataRangeDirectoryName.Contains(":"))
    DataRangeDirectoryName.Remove(DataRangeDirectoryName.First(":"));

  TString ThisFileDirectoryName = "EOSPath" + DataRangeDirectoryName + "/JetTree_data" + TreeContentFlag;

  if (DataType=="BTagMu") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathBTagMu);
  if (DataType=="QCDMu") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathQCDMu);
  if (DataType=="JetHT") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathJetHT);
  if (DataType=="QCD") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathQCD);
  //if (DataType=="TTP") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathTTbarPowheg);
  //if (DataType=="TTM") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathTTbarMadGraph);
  //if (DataType=="Hbb") ThisFileDirectoryName.ReplaceAll("EOSPath", EOSPathGluGluHbb);

  if (ThisFileDirectoryName.Contains("QCD")) 
    ThisFileDirectoryName.ReplaceAll("/JetTree_data", "/JetTree_mc");

  if ((*DataRangeName).Contains("July09")) ThisFileDirectoryName.ReplaceAll("CMSSW_8_0_12", "CMSSW_8_0_8_patch1");

  *FileDirectoryName = ThisFileDirectoryName;
   
  *NumberOfTrees = GetNumberOfTrees(DataType, DataRange);

  *FirstTree = 1;

  if ((*DataRangeName).Contains(":")) {

    TString PartitionFlag = *DataRangeName;
    PartitionFlag.Replace(0, PartitionFlag.First(":")+1, "");

    TString SectionFlag = PartitionFlag;
    SectionFlag.Remove(SectionFlag.First("o"));
    int Section = SectionFlag.Atoi();

    TString nSectionsFlag = PartitionFlag;
    nSectionsFlag.Replace(0, nSectionsFlag.First("f")+1, "");
    int nSections = nSectionsFlag.Atoi();
    
    int NumberOfTreesPerSection = (*NumberOfTrees)/nSections; 
    *FirstTree = 1 + (Section - 1)*NumberOfTreesPerSection; 
    if (Section<nSections) *NumberOfTrees = Section*NumberOfTreesPerSection;
    
  }
  
}

TString GetFileName(TString FileDirectoryName, int FileNumber) {

  TString FileName = FileDirectoryName + "_"; FileName += FileNumber; FileName += ".root";

  return FileName;

}

double DeltaPhi(double phi1, double phi2) {

  double DeltaPhi = fabs(phi2 - phi1);
  if (DeltaPhi>3.141593) DeltaPhi -= 2.*3.141593;
  return fabs(DeltaPhi);
  
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) {

  double DeltaPhi = fabs(phi2 - phi1);
  if (DeltaPhi>3.141593) DeltaPhi -= 2.*3.141593;
  return sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
  
}

bool IsTaggedJet(int ijet, TString ThisTaggerName, TString JetType = "Default") {

  if (ThisTaggerName=="NONE") 
    return true; 

  if (ThisTaggerName=="NOTAG") 
    return true;

  float ThisDiscriminatorCut = -1;
  for (int tg = 0; tg<nBTagAnalyzerTaggers; tg++)
    if (ThisTaggerName==BTagAnalyzerTaggerName[tg])
      ThisDiscriminatorCut = BTagAnalyzerTaggerCut[tg];

  if (ThisDiscriminatorCut==-1)
    std::cout << "Wrong discriminator cut " << ThisTaggerName << std::endl;

  double ThisDiscriminatorValue = -1;

  if (JetType=="Default") {

    if (ThisTaggerName.Contains("JBP"))
      ThisDiscriminatorValue = Jet_Bprob[ijet];
    else if (ThisTaggerName.Contains("JP"))
      ThisDiscriminatorValue = Jet_Proba[ijet]; 
    else if (ThisTaggerName.Contains("CSVv2"))
      ThisDiscriminatorValue = Jet_CombIVF[ijet];  
    else if (ThisTaggerName.Contains("CMVAv2"))
      ThisDiscriminatorValue = Jet_cMVAv2[ijet]; 
    else if (ThisTaggerName.Contains("TCHP")) 
      ThisDiscriminatorValue = Jet_Ip3P[ijet];
    else if (ThisTaggerName.Contains("TCHE")) 
      ThisDiscriminatorValue = Jet_Ip2P[ijet];
    else if (ThisTaggerName.Contains("CSV"))
      ThisDiscriminatorValue = Jet_CombSvx[ijet];
    else if (ThisTaggerName.Contains("SSVHE"))
      ThisDiscriminatorValue = Jet_Svx[ijet];

  } else if (JetType=="Standard") {

    if (ThisTaggerName.Contains("JBP"))
      ThisDiscriminatorValue = StdJet_Bprob[ijet];
    else if (ThisTaggerName.Contains("JP"))
      ThisDiscriminatorValue = StdJet_Proba[ijet]; 
    else if (ThisTaggerName.Contains("CSVv2"))
      ThisDiscriminatorValue = StdJet_CombIVF[ijet];
    else if (ThisTaggerName.Contains("cMVAv2"))
      ThisDiscriminatorValue = StdJet_cMVAv2[ijet];

  } 

  if (ThisDiscriminatorValue==-1)
    std::cout << "Wrong discriminator value " << ThisTaggerName << std::endl;

  return ThisDiscriminatorValue>ThisDiscriminatorCut;
  
  std::cout << "Wrong tagger name " << ThisTaggerName << std::endl;

  return false;

}

bool IsTaggedJetOld(int ijet, TString ThisTaggerName) {

  if (ThisTaggerName=="NONE") 
    return true; 

  if (ThisTaggerName=="NOTAG") 
    return true;

  if (ThisTaggerName=="TCHPL")
    return (Jet_Ip3P[ijet]>1.19);

  if (ThisTaggerName=="TCHET") 
    return Jet_Ip2P[ijet]>10.2;

  if (ThisTaggerName=="TCHPM")
    return Jet_Ip3P[ijet]>1.93;

  if (ThisTaggerName=="TCHPT")
    return Jet_Ip3P[ijet]>3.41;

  if (ThisTaggerName=="CSVL")
    return Jet_CombSvx[ijet]>0.244;

  if (ThisTaggerName=="CSVM")
    return Jet_CombSvx[ijet]>0.679;

  if (ThisTaggerName=="CSVT")
    return Jet_CombSvx[ijet]>0.898;

  if (ThisTaggerName=="CSVS")
    return Jet_CombSvx[ijet]>0.930;

  if (ThisTaggerName=="CSVU")
    return Jet_CombSvx[ijet]>0.960;

  if (ThisTaggerName=="CSLL")
    return Jet_CombCSVJPSL[ijet]>0.531;

  if (ThisTaggerName=="CSLM")
    return Jet_CombCSVJPSL[ijet]>0.758;

  if (ThisTaggerName=="CSLT")
    return Jet_CombCSVJPSL[ijet]>0.849;

  if (ThisTaggerName=="CSLS")
    return Jet_CombCSVJPSL[ijet]>0.900;

  if (ThisTaggerName=="CSLU")
    return Jet_CombCSVJPSL[ijet]>0.950;

  if (ThisTaggerName=="RCVL")
    return Jet_RetCombSvx[ijet]>0.405;

  if (ThisTaggerName=="RCVM")
    return Jet_RetCombSvx[ijet]>0.783;

  if (ThisTaggerName=="RCVT")
    return Jet_RetCombSvx[ijet]>0.920;

  if (ThisTaggerName=="RCVS")
    return Jet_RetCombSvx[ijet]>0.945;

  if (ThisTaggerName=="RCVU")
    return Jet_RetCombSvx[ijet]>0.970;

  if (ThisTaggerName=="SSVHEL")
    return Jet_Svx[ijet]> 3.05;

  if (ThisTaggerName=="TCHEL")
    return Jet_Ip2P[ijet]> 1.7;
 
  if (ThisTaggerName=="JBPT")
    return Jet_Bprob[ijet]>3.74;
 
  if (ThisTaggerName=="JBPT1")
    return Jet_Bprob[ijet]>3.14;
 
  if (ThisTaggerName=="JBPM")
    return Jet_Bprob[ijet]>2.55;
 
  if (ThisTaggerName=="JBPM1")
    return Jet_Bprob[ijet]>2.43;
 
  if (ThisTaggerName=="JBPM2")
    return Jet_Bprob[ijet]>2.31;
 
  if (ThisTaggerName=="JBPM3")
    return Jet_Bprob[ijet]>2.18;
 
  if (ThisTaggerName=="JBPM4")
    return Jet_Bprob[ijet]>2.06;
 
  if (ThisTaggerName=="JBPM5")
    return Jet_Bprob[ijet]>1.94;
 
  if (ThisTaggerName=="JBPM6")
    return Jet_Bprob[ijet]>1.82;
 
  if (ThisTaggerName=="JBPM7")
    return Jet_Bprob[ijet]>1.70;
 
  if (ThisTaggerName=="JBPM8")
    return Jet_Bprob[ijet]>1.57;
 
  if (ThisTaggerName=="JBPM9")
    return Jet_Bprob[ijet]>1.45;
 
  if (ThisTaggerName=="JBPL")
    return Jet_Bprob[ijet]>1.33;

  if (ThisTaggerName=="JPL")
    return Jet_Proba[ijet]>0.245;
 
  if (ThisTaggerName=="JPM")
    return Jet_Proba[ijet]>0.515;
 
  if (ThisTaggerName=="JPT")
    return Jet_Proba[ijet]>0.760;
 
  if (ThisTaggerName=="CSVv2L")
    return Jet_CombIVF[ijet]>0.460;
 
  if (ThisTaggerName=="CSVv2M")
    return Jet_CombIVF[ijet]>0.800;
 
  if (ThisTaggerName=="CSVv2T")
    return Jet_CombIVF[ijet]>0.935;
 
  if (ThisTaggerName=="CMVAv2L")
    return Jet_CombIVF[ijet]>-0.715;
 
  if (ThisTaggerName=="CMVAv2M")
    return Jet_CombIVF[ijet]>0.185;
 
  if (ThisTaggerName=="CMVAv2T")
    return Jet_CombIVF[ijet]>0.875;
  
  std::cout << "Wrong tagger name " << ThisTaggerName << std::endl;

  return false;

}

int GetMuonJet(int *iMu) {

  int jMu = -1;

  int nSoftMuon  = 0;

  for (int nmu = 0; nmu<nMuon; nmu++) {

    if (Muon_pt[nmu]<=5.) continue; 
    if (fabs(Muon_eta[nmu])>=2.4) continue; 
    if (Muon_isGlobal[nmu]==0) continue; 
    if (Muon_nMuHit[nmu]==0) continue; 
    if (Muon_nMatched[nmu]<=1) continue; 
    if (Muon_nTkHit[nmu]<=10) continue; 
    if (Muon_nPixHit[nmu]<=1) continue;
    if (Muon_nOutHit[nmu]>=3) continue;
    if (Muon_chi2[nmu]>=10.) continue;
    if (Muon_chi2Tk[nmu]>=10.) continue;
    if (fabs(Muon_vz[nmu]-PVz)>=1.) continue;

    jMu = Muon_IdxJet[nmu];
    
    *iMu = nmu;
    
    nSoftMuon++;
      
  }

  if (nSoftMuon!=1) return -1;

  if (jMu<0) return -1;

  return jMu;

}

int GetPFMuonJet(int *yMu, bool IsUnique = true) {

  int kMu = -1;

  int nSoftMuon  = 0;

  float SoftMuonPt = -1.;
  
  for (int nmu = 0; nmu<nPFMuon; nmu++) {

    if (PFMuon_pt[nmu]<=5.) continue; 
    if (fabs(PFMuon_eta[nmu])>=2.4) continue; 
    if (PFMuon_GoodQuality[nmu]<2) continue; 
    
    /*int JetIdx = -1; float MinDR = 0.4;
    for (int ijet = 0; ijet<nJet; ijet++) 
      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], PFMuon_eta[nmu], PFMuon_phi[nmu])<MinDR) {

	MinDR = DeltaR(Jet_eta[ijet], Jet_phi[ijet], PFMuon_eta[nmu], PFMuon_phi[nmu]);
	JetIdx = ijet;

      }
    */

    int JetIdx = PFMuon_IdxJet[nmu];
    
    if (JetIdx>=0) {

      if (PFMuon_pt[nmu]>SoftMuonPt) {

	kMu = JetIdx;
	
	*yMu = nmu;

	SoftMuonPt = PFMuon_pt[nmu];
	
      }
      
      nSoftMuon++;
      
    }

  }

  if (IsUnique && nSoftMuon!=1) return -1;

  if (kMu<0) return -1;

  return kMu;

}

int IsMuonJet(int ijet, bool TightMuon = true) {

  int nSoftMuon  = 0;

  int mu_index = -1;

  for (int nmu = 0; nmu<nMuon; nmu++) 
    if (ijet==Muon_IdxJet[nmu]) {

      if (TightMuon) {
	if (Muon_pt[nmu]<=5.) continue; 
	if (fabs(Muon_eta[nmu])>=2.4) continue; 
	if (Muon_isGlobal[nmu]==0) continue; 
	if (Muon_nMuHit[nmu]==0) continue; 
	if (Muon_nMatched[nmu]<=1) continue; 
	if (Muon_nTkHit[nmu]<=10) continue; 
	if (Muon_nPixHit[nmu]<=1) continue;
	if (Muon_nOutHit[nmu]>=3) continue;
	if (Muon_chi2[nmu]>=10.) continue;
	if (Muon_chi2Tk[nmu]>=10.) continue;
	if (fabs(Muon_vz[nmu]-PVz)>=1.) continue;
      }

      //int ptBin = -1;
      //for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
      //if (Jet_pt[Muon_IdxJet[nmu]]>PtRelPtEdge[ptb]) ptBin = ptb;
      //if (Muon_pt[nmu]<=MuonPtCut[ptBin] && TightMuon) continue;  
          
      mu_index = nmu;

      nSoftMuon++;
      
    }
  
  //if (nSoftMuon!=1) return -1;
  
  return mu_index;

}

bool PassTriggerBit(TString ThisCode) {

  if (ThisCode=="NONE") return true;

  if (ThisCode=="") {
    
    std::cout << "Please select a trigger code" << std::endl;
    return false;

  }

  int triggerIdx = -1;
 
  if (ThisCode=="_DiJet20X") triggerIdx =   2;
  if (ThisCode=="_DiJet40X") triggerIdx =   3;
  if (ThisCode=="_Jet30_v")  triggerIdx =   1;
  if (ThisCode=="_Jet60_v")  triggerIdx =   4;
  if (ThisCode=="_Jet110_v") triggerIdx =   9;
  if (ThisCode=="_Jet300_v") triggerIdx =  18;
  if (ThisCode=="_PFJet40")  triggerIdx =   2;
  if (ThisCode=="_PFJet60")  triggerIdx =  87;
  if (ThisCode=="_PFJet80")  triggerIdx =   7;
  if (ThisCode=="_PFJet140") triggerIdx =  12;
  if (ThisCode=="_PFJet200") triggerIdx =  15;
  if (ThisCode=="_PFJet260") triggerIdx =  17;
  if (ThisCode=="_PFJet320") triggerIdx =  19;
  if (ThisCode=="_PFJet400") triggerIdx =  21;
  if (ThisCode=="_PFJet60")  triggerIdx =  87;
  if (ThisCode=="_PFJet450") triggerIdx =  88;
  if (ThisCode=="_PFJet500") triggerIdx =  89;
  if (ThisCode=="_DiJet20")  triggerIdx =  35;
  if (ThisCode=="_DiJet40")  triggerIdx =  41;
  if (ThisCode=="_DiJet70")  triggerIdx =  44;
  if (ThisCode=="_DiJet110") triggerIdx =  47;
  if (ThisCode=="_Jet300")   triggerIdx =  50; // 18 only for 7 TeV data

  if (triggerIdx==-1) {
    
    std::cout << "Please select an existing trigger code" << std::endl;
    return false;

  }

  int bitIdx = int(triggerIdx/32);
  if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) return true;
  
  return false;

}

bool PassTriggerBit(int ThisBit, int ThisCode) {

  if (ThisCode==-1) return true;

  int ThisMagnitude = 0;
  for (int mg = 0; mg<10; mg++) 
    if (ThisCode/pow(10, mg)>=1.) ThisMagnitude = mg;

  int ThisBitRounded = ThisBit/pow(10, ThisMagnitude+1);
  int ThisBitCut = ThisBit/pow(10, ThisMagnitude);
  int ThisBitCode = ThisBitCut - 10*ThisBitRounded;

  int ThisCodeCode = ThisCode/pow(10, ThisMagnitude);

  return (ThisBitCode & ThisCodeCode);

}

bool IsFromGluonSplitting(int ijet, int SplitFlavour = 5, bool isPythia8 = false) {

  int qStatus = 2;
  if (isPythia8) qStatus = 71;

  if (SplitFlavour==4) {
    
    for (int k = 0; k<ncQuarks-1; k++) 	
      if (cQuark_status[k]==qStatus)       
	if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], cQuark_eta[k], cQuark_phi[k])<=0.5) 
	  for (int l = k+1; l<ncQuarks; l++) 
	    if (cQuark_status[l]==qStatus) 
	      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], cQuark_eta[l], cQuark_phi[l])<=0.5) 		  
		if (cQuark_pdgID[k]*cQuark_pdgID[l]<0) return true;
 
  } else if (SplitFlavour==5 || SplitFlavour==21) {
 
    for (int k = 0; k<nbQuarks-1; k++) 
      if (bQuark_status[k]==qStatus) 
	if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], bQuark_eta[k], bQuark_phi[k])<=0.5)  
	  for (int l = k+1; l<nbQuarks; l++)
	    if (bQuark_status[l]==qStatus) 
	      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], bQuark_eta[l], bQuark_phi[l])<=0.5)
		if (bQuark_pdgID[k]*bQuark_pdgID[l]<0) return true;
 
  } else if (SplitFlavour==21) {

    for (int k = 0; k<nbQuarks; k++)
      if (bQuark_status[k]==qStatus && bQuark_fromGSP[k]==1) 	
	if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], bQuark_eta[k], bQuark_phi[k])<=0.5) return true;

  }

  return false;
  
}

bool IsFromGluonSplittingFromHadron(int ijet, int SplitFlavour = 5) {

  int nHadrons = 0;

  if (SplitFlavour==4) {

    for (int k = 0; k<nDHadrons; k++)
      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], DHadron_eta[k], DHadron_phi[k])<=0.4) {

	bool GoodDHadron = true;

	for (int l = 0; l<nBHadrons; l++)
	  if (BHadron_DHadron1[l]==k || BHadron_DHadron2[l]==k) 
	    GoodDHadron = false;

	if (GoodDHadron) nHadrons++;
	
      }
 
  } else if (SplitFlavour==5 || SplitFlavour==21) {

    for (int k = 0; k<nBHadrons; k++)
      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], BHadron_eta[k], BHadron_phi[k])<=0.4) 
	if (BHadron_hasBdaughter[k]==0)
	  nHadrons++;

  } 

  if (nHadrons>=2) return true;

  return false;
  
}

bool HasMuon(int ThisJet, bool TightCut) {

  if (TightCut && nMuon>0) return true;

  for (int nm = 0; nm<nMuon; nm++)
    if (Muon_IdxJet[nm]==ThisJet) return true;

  return false;

}

bool HasPFMuon(int ThisJet, bool TightCut) {

  if (TightCut && nPFMuon>0) return true;

  for (int nm = 0; nm<nPFMuon; nm++) 
    if (PFMuon_GoodQuality[nm]>=1)
      if (PFMuon_IdxJet[nm]==ThisJet) return true;/* {
    
      float MinDR = 0.4;
      for (int ijet = 0; ijet<nJet; ijet++) 
	if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], PFMuon_eta[nm], PFMuon_phi[nm])<MinDR) 
	  return true;
      
	  }*/
  
  return false;

}

bool IsPFMuonEvent(bool TightCut) {

  if (TightCut && nPFMuon>0) return true;

  for (int nm = 0; nm<nPFMuon; nm++) 
    if (PFMuon_GoodQuality[nm]>=1) return true;
  
  return false;

}

bool HasTaggedJet(int ThisJet, bool inclusive = false, TString TaggerVeto = "TCHEL") {

  for (int ijet = 0; ijet<nJet; ijet++)
    if (ijet!=ThisJet || inclusive)
      if (IsTaggedJet(ijet, TaggerVeto)) return true;

  return false;

}

bool PassPtHat(TString QCDSample, int ijet) {

  bool ThisPass = true;

  if (QCDSample=="Inclusive") {

    // Fall 11
//     if (pthat<  30. && Jet_pt[ijet]>  85.) ThisPass = false;
//     if (pthat<  50. && Jet_pt[ijet]> 100.) ThisPass = false;
//     if (pthat<  80. && Jet_pt[ijet]> 150.) ThisPass = false;
//     if (pthat< 120. && Jet_pt[ijet]> 250.) ThisPass = false;
//     if (pthat< 170. && Jet_pt[ijet]> 320.) ThisPass = false;
//     //if (pthat< 300. && Jet_pt[ijet]> 400.) ThisPass = false;
//     if (pthat< 470. && Jet_pt[ijet]> 620.) ThisPass = false;
//     if (pthat< 600. && Jet_pt[ijet]> 720.) ThisPass = false;
//     if (pthat< 800. && Jet_pt[ijet]> 920.) ThisPass = false;
//     if (pthat> 800. && Jet_pt[ijet]< 600.) ThisPass = false;

    // Summer 12
    if (pthat<  30. && Jet_pt[ijet]> 200.) ThisPass = false;
    //if (pthat<  50. && Jet_pt[ijet]> 140.) ThisPass = false;
    //if (pthat<  80. && Jet_pt[ijet]> 160.) ThisPass = false;
    if (pthat<  50. && Jet_pt[ijet]> 200.) ThisPass = false;
    if (pthat<  80. && Jet_pt[ijet]> 200.) ThisPass = false;
    if (pthat< 120. && Jet_pt[ijet]> 250.) ThisPass = false;
    if (pthat< 170. && Jet_pt[ijet]> 340.) ThisPass = false;
    if (pthat< 300. && Jet_pt[ijet]> 520.) ThisPass = false;
    //if (pthat< 470. && Jet_pt[ijet]> 620.) ThisPass = false;
    //if (pthat< 600. && Jet_pt[ijet]> 720.) ThisPass = false;
    //if (pthat< 800. && Jet_pt[ijet]> 920.) ThisPass = false;
    //if (pthat> 800. && Jet_pt[ijet]< 600.) ThisPass = false;
    
  } else if (QCDSample=="MuEnrichedPt5") {

    // Fall 11
//     if (pthat<  20. && Jet_pt[ijet]>  60.) ThisPass = false;
//     //if (pthat<  30. && Jet_pt[ijet]>  80.) ThisPass = false;
//     if (pthat<  30. && Jet_pt[ijet]>  85.) ThisPass = false;
//     if (pthat<  50. && Jet_pt[ijet]> 100.) ThisPass = false;
//     if (pthat<  80. && Jet_pt[ijet]> 150.) ThisPass = false;
//     //if (pthat< 120. && Jet_pt[ijet]> 210.) ThisPass = false;
//     if (pthat< 120. && Jet_pt[ijet]> 250.) ThisPass = false;
//     if (pthat< 170. && Jet_pt[ijet]> 320.) ThisPass = false;
//     //if (pthat< 300. && Jet_pt[ijet]> 400.) ThisPass = false;
//     if (pthat< 470. && Jet_pt[ijet]> 620.) ThisPass = false;
//     if (pthat< 600. && Jet_pt[ijet]> 720.) ThisPass = false;
//     if (pthat< 800. && Jet_pt[ijet]> 920.) ThisPass = false;
//     if (pthat> 800. && Jet_pt[ijet]< 600.) ThisPass = false;

    // Summer 12
    if (pthat<  20. && Jet_pt[ijet]>  60.) ThisPass = false;
    if (pthat<  30. && Jet_pt[ijet]>  85.) ThisPass = false;
    if (pthat<  50. && Jet_pt[ijet]> 120.) ThisPass = false;
    if (pthat<  80. && Jet_pt[ijet]> 160.) ThisPass = false;
    if (pthat< 120. && Jet_pt[ijet]> 220.) ThisPass = false;
    if (pthat< 170. && Jet_pt[ijet]> 320.) ThisPass = false;
    if (pthat< 300. && Jet_pt[ijet]> 440.) ThisPass = false;
    if (pthat< 470. && Jet_pt[ijet]> 620.) ThisPass = false;
    if (pthat< 600. && Jet_pt[ijet]> 720.) ThisPass = false;
    if (pthat< 800. && Jet_pt[ijet]> 920.) ThisPass = false;
    if (pthat> 800. && Jet_pt[ijet]< 600.) ThisPass = false;
    
  }

  return ThisPass;

}

int GetFatJet(int ijet,int *ajet) {
  
  *ajet = -1;

  for (int fatjet = 0; fatjet<nFatJet; fatjet++) 
    for (int subjet = FatJet_FirstSubJet[fatjet]; subjet<FatJet_FirstSubJet[fatjet]+FatJet_nSubJets[fatjet]; subjet++) 
      if (subjet==ijet) {
      
	if (FatJet_pt[fatjet]<425.)       return -1;
	if (fabs(FatJet_eta[fatjet])>2.4) return -2;
	if (FatJet_looseID[fatjet]==0)    return -3;
	
	double dRfatfatMin = 1000.;
	      
	for (int gatjet = 0; gatjet<nFatJet; gatjet++) 
	  if (gatjet!=fatjet) {
	    
	    double dRfatfat = DeltaR(FatJet_eta[fatjet], FatJet_phi[fatjet],
				     FatJet_eta[gatjet], FatJet_phi[gatjet]);
	    if (dRfatfat<dRfatfatMin ) dRfatfatMin = dRfatfat;
	    
	  }
				
	if (dRfatfatMin<0.8) return -4;

	float fatptjet = FatJet_pt[fatjet];
	float fatptmass = FatJet_mass[fatjet]/fatptjet;
	      
	if (FatJet_nSubJets[fatjet]!=2) return -5;

	int ijet1 = FatJet_FirstSubJet[fatjet]; 
	int ijet2 = ijet1 + 1;

	if (ijet==ijet1) *ajet = ijet2;
	else if (ijet==ijet2) *ajet = ijet1;
 
	if (ijet1==ijet2) return -6;
	if (Jet_pt[ijet1]==0. || Jet_pt[ijet2]==0.) return -7;
	//if (Jet_pt[ijet1]<20 || Jet_pt[ijet2]<20) return -8;
       
	double dRsubsub = DeltaR( Jet_eta[ijet1], Jet_phi[ijet1],
				  Jet_eta[ijet2], Jet_phi[ijet2] );
 
	// Infrared-safe cut
	//if (dRsubsub<fatptmass) return -9;// Still needed now???

	// additional delta R cut
	if ( dRsubsub > 0.825 ) return -10;
     
	//if ( dRsubsub < 0.5 ) continue; // study in DR bins

	if (FatJet_tau2[fatjet]/FatJet_tau1[fatjet]>=0.5) return -11;

	return fatjet;

      }

  return -999;

}

double TrackPhiSolution(double SinPhi, double alpha, double beta, double gamma, double *DS1, double *DS2) {

  double phi_trk_1 = asin(SinPhi);
  double phi_trk_2 = 3.14159265358979312 - phi_trk_1;

  double cos_phi_trk_1 = cos(phi_trk_1);
  double cos_phi_trk_2 = cos(phi_trk_2);
  
  double D1 = alpha - beta*cos_phi_trk_1 - gamma*SinPhi;
  double D2 = alpha - beta*cos_phi_trk_2 - gamma*SinPhi;

  *DS1 = D1;
  *DS2 = D2;

  //std::cout << "DS " << D1 << " " << D2 << std::endl;

  if (fabs(D1)<fabs(D2)) return phi_trk_1;
  else return phi_trk_2;

}

double GetTrackPhi(double pt_rel, double pt_trk, double eta_jet, double phi_jet, double eta_trk) {
  
  double sin_rel = pt_rel/(pt_trk*cosh(eta_trk));
  double cos_rel = sqrt(1- sin_rel*sin_rel);

  double theta_jet = 2.*atan(exp(-1.*eta_jet));
  double theta_trk = 2.*atan(exp(-1.*eta_trk));

  double alpha = cos_rel - cos(theta_jet)*cos(theta_trk);
  double beta  = sin(theta_jet)*cos(phi_jet)*sin(theta_trk);
  double gamma = sin(theta_jet)*sin(phi_jet)*sin(theta_trk);

  double A  = 2*alpha*gamma;
  double B  = 4*beta*beta*(gamma*gamma + beta*beta - alpha*alpha);
  double C  = 2*(gamma*gamma + beta*beta);

  if (fabs(asin(sin_rel))<fabs(theta_jet-theta_trk)) return phi_jet;

//   if (fabs(B)<0.0001) { 

//     B = 0.;

//     if (fabs(A/C)>1.)
//       if (fabs(A/C)-1.<0.0005) A = A/fabs(A)*fabs(C);

//   }

  double S1 = (A + sqrt(B))/C;
  double S2 = (A - sqrt(B))/C;

  double DS1, DS2;

  if (fabs(S1)<=1. && fabs(S2)>1.) {

    return TrackPhiSolution(S1, alpha, beta, gamma, &DS1, &DS2);

  } else if (fabs(S1)>1. && fabs(S2)<=1.) {
    
    return TrackPhiSolution(S2, alpha, beta, gamma, &DS1, &DS2);

  } else if (fabs(S1)<=1. && fabs(S2)<=1.) {

    double phi_trk_1 = TrackPhiSolution(S1, alpha, beta, gamma, &DS1, &DS2);
    double phi_trk_2 = TrackPhiSolution(S2, alpha, beta, gamma, &DS1, &DS2);
    
    double DeltaPhi_1 = fabs(phi_jet - phi_trk_1);
    if (DeltaPhi_1>3.141593) DeltaPhi_1 -= 2.*3.141593;
    
    double DeltaPhi_2 = fabs(phi_jet - phi_trk_2);
    if (DeltaPhi_2>3.141593) DeltaPhi_2 -= 2.*3.141593;

    if (fabs(DeltaPhi_1)<0.4 && fabs(DeltaPhi_2)>=0.4) {

      return phi_trk_1;

    } if (fabs(DeltaPhi_1)>=0.4 && fabs(DeltaPhi_2)<0.4) {

      return phi_trk_2;

    } if (fabs(DeltaPhi_1)<0.4 && fabs(DeltaPhi_2)<0.4) { 
    
      //std::cout << "Utilities::GetTrackPhi -> two solutions " << DeltaPhi_1 << " " << DeltaPhi_2 << std::endl;
      if (fabs(DeltaPhi_1)<fabs(DeltaPhi_2)) return phi_trk_1;
      else return phi_trk_2;

    } else {

      std::cout << "Utilities::GetTrackPhi -> no good solutions " << DeltaPhi_1 << " " << DeltaPhi_2 << std::endl;
      std::cout << "      " << A << " " << B << " " << C << " " << beta << " " << gamma << " " << alpha << " " << " " << phi_jet << " " << theta_jet-theta_trk << " (" << theta_jet << ")" << std::endl << "      " << pt_rel << " " << pt_trk << " " << sin_rel << " " << cos_rel << " " << DS1 << " " << DS2 << std::endl;
      return -999.;

    }

  } else {

    std::cout << "Utilities::GetTrackPhi -> no solutions " << S1 << " " << S2 << std::endl;
    std::cout << "      " << A << " " << B << " " << C << " " << beta << " " << gamma << " " << alpha << " " << " " << phi_jet << " " << theta_jet-theta_trk << " (" << theta_jet << ", " << eta_trk << ")" << std::endl << "      " << pt_rel << " " << pt_trk << " " << sin_rel << " " << cos_rel << " " << DS1 << " " << DS2 << std::endl;
    return -999.;
    
  }
  
}

double DeltaRfromPtRel(double pt_rel, double pt_trk, double eta_jet, double phi_jet, double eta_trk) {
	
  double TrkPhi = GetTrackPhi(pt_rel, pt_trk, eta_jet, phi_jet, eta_trk);
  
  if (TrkPhi<-100.) return 9999.;
  
  return DeltaR(eta_jet, phi_jet, eta_trk, TrkPhi);

}

int GetBHadron(int ijet, float DRCut = 0.5) {

  int iBHadron = -1;
  
  float B_Mass = 0.;
  
  for (int ib = 0; ib<nBHadrons; ib++) 
    if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], BHadron_eta[ib], BHadron_phi[ib])<DRCut) {

      if (BHadron_mass[ib]>B_Mass) {

	B_Mass = BHadron_mass[ib];
	iBHadron = ib;
      
      }

    }

  return iBHadron;

}

float GetEnergy(float Pt, float Mass, float Eta) {

  float Momentum = Pt*cosh(Eta);

  return sqrt(Mass*Mass + Momentum*Momentum);

}

// Till Spring16_25nsV6
/*const int nJEUEtaBins = 40, nJEUPtPoints = 44;
float JEUEtaEdge[nJEUEtaBins+1] = {-5.4, -5.0, -4.4, -4.0, -3.5, -3.0, -2.8, -2.6, -2.4, 
				   -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 
				   0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 
				   2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.4};
float JEUPtPoint[nJEUPtPoints] = {9.0, 11.0, 13.5, 16.5, 19.5, 22.5, 26.0, 30.0, 34.5, 40., 46.0, 52.5, 60.0, 
				  69.0, 79.0, 90.5, 105.5, 123.5, 143.0, 163.5, 185.0, 208.0, 232.5, 258.5, 286.0,
				  331.0, 396.0, 468.5, 549.5, 639.0, 738.0, 847.5, 968.5, 1102.0, 1249.5, 1412.0, 
				  1590.5, 1787.0, 1945.0, 2119.0, 2369.0, 2643.5, 2945.0, 3276.5};
float Jet_JEU[2][nJEUEtaBins][nJEUPtPoints];
int BhoFactor = 132;*/
//From Spring16_25nsV6
const int nJEUEtaBins = 40, nJEUPtPoints = 50;
float JEUEtaEdge[nJEUEtaBins+1] = {-5.4, -5.0, -4.4, -4.0, -3.5, -3.0, -2.8, -2.6, -2.4, 
				   -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 
				   0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 
				   2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.4};
float JEUPtPoint[nJEUPtPoints] = {9.0, 11.0, 13.5, 16.5, 19.5, 22.5, 26.0, 30.0, 34.5, 40., 46.0, 52.5, 60.0, 
				  69.0, 79.0, 90.5, 105.5, 123.5, 143.0, 163.5, 185.0, 208.0, 232.5, 258.5, 286.0,
				  331.0, 396.0, 468.5, 549.5, 639.0, 738.0, 847.5, 968.5, 1102.0, 1249.5, 1412.0, 
				  1590.5, 1787.0, 2003.0, 2241.0, 2503.0, 2790.5, 3107.0, 3455.0, 3837.0, 4257.0, 
				  4719.0, 5226.5, 5784.0, 6538.0};
float Jet_JEU[2][nJEUEtaBins][nJEUPtPoints];
int BhoFactor = 150;

float GetJEU(int DataType, float jeteta, float jetpt) {

  int ThisEtaBin = -1, ThisPtPoint = -1;
  
  for (int eb = 0; eb<nJEUEtaBins; eb++)
    if (jeteta>JEUEtaEdge[eb]) ThisEtaBin = eb;

  for (int pp = 0; pp<nJEUPtPoints; pp++) 
    if (jetpt>JEUPtPoint[pp]) ThisPtPoint = pp;

  float DY = Jet_JEU[DataType][ThisEtaBin][ThisPtPoint+1] - Jet_JEU[DataType][ThisEtaBin][ThisPtPoint];
  float DX = JEUPtPoint[ThisPtPoint+1] - JEUPtPoint[ThisPtPoint];
  float Offset = Jet_JEU[DataType][ThisEtaBin][ThisPtPoint];

  return DY/DX*(jetpt - JEUPtPoint[ThisPtPoint]) + Offset;
  
}

void ReadJetEnergyUncertainty(TString DataType, TString CorrectionVersion, TString JetType) {

  int dt = 0;
  if (DataType=="_MC") dt = 1;

  TString JEUFileName = "./JEU/" + CorrectionVersion + DataType + "/" + CorrectionVersion + DataType + "_Uncertainty" + JetType + ".txt";
  ifstream JEUFile; JEUFile.open(JEUFileName);
  
  if (!JEUFile)
    throw std::invalid_argument("JEU file not found!");
  
  TString DummyString;
  for (int dd = 0; dd<7; dd++) 
    JEUFile >> DummyString;
  
  for (int eb = 0; eb<nJEUEtaBins; eb++) {
    
    float etamin, etamax; int Bho;
    JEUFile >> etamin >> etamax >> Bho;
    
    if (etamin!=JEUEtaEdge[eb]) cout << "ReadJetEnergyUncertainty: warning -> etamin = " << etamin << endl;
    if (Bho!=BhoFactor) cout << "ReadJetEnergyUncertainty: warning -> Bho = " << Bho << endl;
    
    for (int pp = 0; pp<nJEUPtPoints; pp++) {
      
      float ptpoint, jeudw, jeuup;
      JEUFile >> ptpoint >> jeudw >> jeuup;
      
      if (ptpoint!=JEUPtPoint[pp]) cout << "ReadJetEnergyUncertainty: warning -> ptpoint = " << ptpoint << endl;
      if (jeudw!=jeuup) cout << "ReadJetEnergyUncertainty: warning -> jeudw = " << jeudw << ", jeuup = " << jeuup << endl;
      
      Jet_JEU[dt][eb][pp] = jeudw;
      
    }
    
  }
  
}

/*
TH1D *ReshapeHistogram(TH1D *THISHISTO, TString ThisVariable) {

  int nBins = THISHISTO->GetNbinsX();
  float xi = THISHISTO->GetBinLowEdge(1);
  float xf = THISHISTO->GetBinLowEdge(nBins) + THISHISTO->GetBinWidth(1);

  TH1D *INTHISTO = new TH1D("INTHISTO", "", nBins, xi, xf);

  for (int ib = 1; ib<=nBins+1; ib++) {

    float ThisIntegral = THISHISTO->Integral(ib, nBins+1);
    INTHISTO->SetBinContent(ib, ThisIntegral);

  }

  float ThisCut[4], ThisSF[4] = {0.98, 0.96, 0.92, 1.};
  if (ThisVariable=="CSV") { 
    ThisCut[0] = 0.244; ThisCut[1] = 0.679; ThisCut[2] = 0.898; ThisCut[3] = 1.001; 
  } else if (ThisVariable=="JP") {
    ThisCut[0] = 0.275; ThisCut[1] = 0.545; ThisCut[2] = 0.790; ThisCut[3] = 2.5001; 
  } else if (ThisVariable=="PtRel") {
    ThisCut[0] = 0.2; ThisCut[1] = 1.; ThisCut[2] = 2.; ThisCut[3] = 5.001; 
  }

  float LowSF = 1.; int FirstBin = 1;
  for (int icut = 0; icut<4; icut++) {

    float HighSF = ThisSF[icut];

    int LastBin = INTHISTO->FindBin(ThisCut[icut]);
    if (icut<3)
      if (fabs(ThisCut[icut]-INTHISTO->GetBinLowEdge(LastBin))<(INTHISTO->GetBinWidth(LastBin)/2.)) LastBin -= 1.;
    
    float PassSF = (HighSF - LowSF)/(LastBin-FirstBin+1);
    
    for (int ib = FirstBin; ib<=LastBin; ib++) {

      float ThisScale = LowSF + (ib-FirstBin)*PassSF;
      float NewContent = ThisScale*INTHISTO->GetBinContent(ib);
      INTHISTO->SetBinContent(ib, NewContent);

    }

    LowSF = HighSF;
    
    FirstBin = LastBin + 1;

  }

  TH1D *RESHISTO = new TH1D("RESHISTO", "", nBins, xi, xf);

  float LastContent = 0.;
  for (int ib = nBins+1; ib>=1; ib--) {
  
    float ThisContent = INTHISTO->GetBinContent(ib) - LastContent;
    RESHISTO->SetBinContent(ib, ThisContent);

    LastContent = INTHISTO->GetBinContent(ib);
    
  }

//   THISHISTO->Rebin(10);
//   RESHISTO->Rebin(10);

//   THISHISTO->Draw();
//   RESHISTO->SetLineColor(2);  
//   RESHISTO->Draw("same");

//   float restotal = RESHISTO->Integral(1, nBins+1); 
//   float thistotal = THISHISTO->Integral(1, nBins+1); 
//   cout << restotal << " " << thistotal << endl;
//   for (int icut = 0; icut<3; icut++) {

//     int tt = RESHISTO->FindBin(ThisCut[icut]);
//     cout << RESHISTO->Integral(tt, nBins+1)/restotal << " " << THISHISTO->Integral(tt, nBins+1)/thistotal << endl; 

//   }

  return RESHISTO;

}


*/

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr  7 10:47:28 2011 by ROOT version 5.27/06b
// from TTree ttree/ttree
// found on file: JetTree.root
//////////////////////////////////////////////////////////

#ifndef JetTree_h
#define JetTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class JetTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           BitTrigger;
   Int_t           nJet;
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         pthat;
   Int_t           nPU;
   Int_t           PU_bunch[50];   //[nPU]
   Float_t         PU_z[50];   //[nPU]
   Float_t         PU_sumpT_low[50];   //[nPU]
   Float_t         PU_sumpT_high[50];   //[nPU]
   Int_t           PU_ntrks_low[50];   //[nPU]
   Int_t           PU_ntrks_high[50];   //[nPU]
   Float_t         Jet_pt[50];   //[nJet]
   Float_t         Jet_et[50];   //[nJet]
   Float_t         Jet_residual[50];   //[nJet]
   Float_t         Jet_jes[50];   //[nJet]
   Float_t         Jet_eta[50];   //[nJet]
   Float_t         Jet_phi[50];   //[nJet]
   Int_t           Jet_ntracks[50];   //[nJet]
   Int_t           Jet_flavour[50];   //[nJet]
   Float_t         Jet_Ip1N[50];   //[nJet]
   Float_t         Jet_Ip1P[50];   //[nJet]
   Float_t         Jet_Ip2N[50];   //[nJet]
   Float_t         Jet_Ip2P[50];   //[nJet]
   Float_t         Jet_Ip3N[50];   //[nJet]
   Float_t         Jet_Ip3P[50];   //[nJet]
   Float_t         Jet_ProbaN[50];   //[nJet]
   Float_t         Jet_ProbaP[50];   //[nJet]
   Float_t         Jet_Proba[50];   //[nJet]
   Float_t         Jet_BprobN[50];   //[nJet]
   Float_t         Jet_Bprob[50];   //[nJet]
   Float_t         Jet_SvxN[50];   //[nJet]
   Float_t         Jet_Svx[50];   //[nJet]
   Int_t           Jet_SvxNTracks[50];   //[nJet]
   Int_t           Jet_SvxTracks[50];   //[nJet]
   Float_t         Jet_SvxNHP[50];   //[nJet]
   Float_t         Jet_SvxHP[50];   //[nJet]
   Float_t         Jet_SvxMass[50];   //[nJet]
   Float_t         Jet_CombSvxN[50];   //[nJet]
   Float_t         Jet_CombSvx[50];   //[nJet]
   Float_t         Jet_SoftMuN[50];   //[nJet]
   Float_t         Jet_SoftMu[50];   //[nJet]
   Int_t           Jet_hist1[50];   //[nJet]
   Int_t           Jet_hist2[50];   //[nJet]
   Int_t           Jet_hist3[50];   //[nJet]
   Int_t           Jet_histJet[50];   //[nJet]
   Int_t           Jet_histSvx[50];   //[nJet]
   Int_t           Jet_nFirstTrack[50];   //[nJet]
   Int_t           Jet_nLastTrack[50];   //[nJet]
   Int_t           Jet_nFirstSV[50];   //[nJet]
   Int_t           Jet_nLastSV[50];   //[nJet]
   Int_t           nMuon;
   Int_t           Muon_IdxJet[10];   //[nMuon]
   Int_t           Muon_nMuHit[10];   //[nMuon]
   Int_t           Muon_nTkHit[10];   //[nMuon]
   Int_t           Muon_nPixHit[10];   //[nMuon]
   Int_t           Muon_nOutHit[10];   //[nMuon]
   Int_t           Muon_isGlobal[10];   //[nMuon]
   Int_t           Muon_nMatched[10];   //[nMuon]
   Float_t         Muon_chi2[10];   //[nMuon]
   Float_t         Muon_chi2Tk[10];   //[nMuon]
   Float_t         Muon_pt[10];   //[nMuon]
   Float_t         Muon_eta[10];   //[nMuon]
   Float_t         Muon_ptrel[10];   //[nMuon]
   Float_t         Muon_vz[10];   //[nMuon]
   Int_t           Muon_hist[10];   //[nMuon]
   Int_t           nCFromGSplit;
   Float_t         cFromGSplit_pT[40];   //[nCFromGSplit]
   Float_t         cFromGSplit_eta[40];   //[nCFromGSplit]
   Float_t         cFromGSplit_phi[40];   //[nCFromGSplit]
   Int_t           nBFromGSplit;
   Float_t         bFromGSplit_pT[40];   //[nBFromGSplit]
   Float_t         bFromGSplit_eta[40];   //[nBFromGSplit]
   Float_t         bFromGSplit_phi[40];   //[nBFromGSplit]

   // List of branches
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_PU_bunch;   //!
   TBranch        *b_PU_z;   //!
   TBranch        *b_PU_sumpT_low;   //!
   TBranch        *b_PU_sumpT_high;   //!
   TBranch        *b_PU_ntrks_low;   //!
   TBranch        *b_PU_ntrks_high;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_et;   //!
   TBranch        *b_Jet_residual;   //!
   TBranch        *b_Jet_jes;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
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
   TBranch        *b_Jet_Proba;   //!
   TBranch        *b_Jet_BprobN;   //!
   TBranch        *b_Jet_Bprob;   //!
   TBranch        *b_Jet_SvxN;   //!
   TBranch        *b_Jet_Svx;   //!
   TBranch        *b_Jet_SvxNTracks;   //!
   TBranch        *b_Jet_SvxTracks;   //!
   TBranch        *b_Jet_SvxNHP;   //!
   TBranch        *b_Jet_SvxHP;   //!
   TBranch        *b_Jet_SvxMass;   //!
   TBranch        *b_Jet_CombSvxN;   //!
   TBranch        *b_Jet_CombSvx;   //!
   TBranch        *b_Jet_SoftMuN;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_hist1;   //!
   TBranch        *b_Jet_hist2;   //!
   TBranch        *b_Jet_hist3;   //!
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //!
   TBranch        *b_Jet_nFirstTrack;   //!
   TBranch        *b_Jet_nLastTrack;   //!
   TBranch        *b_Jet_nFirstSV;   //!
   TBranch        *b_Jet_nLastSV;   //!
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
   TBranch        *b_Muon_ptrel;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_hist;   //!
   TBranch        *b_nCFromGSplit;   //!
   TBranch        *b_cFromGSplit_pT;   //!
   TBranch        *b_cFromGSplit_eta;   //!
   TBranch        *b_cFromGSplit_phi;   //!
   TBranch        *b_nBFromGSplit;   //!
   TBranch        *b_bFromGSplit_pT;   //!
   TBranch        *b_bFromGSplit_eta;   //!
   TBranch        *b_bFromGSplit_phi;   //!

   JetTree(TTree *tree=0);
   virtual ~JetTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//$$
//$$   void   Loop();
   void   Loop(Int_t    aNN = 0, 
               Float_t  aTagCut = 0.01, 
               Float_t  aPtMin = 20., 
               Float_t  aPtMax = 1000., 
               Float_t  aEtaMin = 0., 
               Float_t  aEtaMax = 10., 
               Float_t  aFreeCut = 0., 
               Int_t    aIntCut = 0, 
               TString  afilename = "output.root");
//$$
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JetTree_cxx
JetTree::JetTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("JetTree.root");
      if (!f) {
         f = new TFile("JetTree.root");
         f->cd("JetTree.root:/mistag");
      }
      tree = (TTree*)gDirectory->Get("ttree");

   }
   Init(tree);
}

JetTree::~JetTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetTree::Init(TTree *tree)
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
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("PU_bunch", PU_bunch, &b_PU_bunch);
   fChain->SetBranchAddress("PU_z", PU_z, &b_PU_z);
   fChain->SetBranchAddress("PU_sumpT_low", PU_sumpT_low, &b_PU_sumpT_low);
   fChain->SetBranchAddress("PU_sumpT_high", PU_sumpT_high, &b_PU_sumpT_high);
   fChain->SetBranchAddress("PU_ntrks_low", PU_ntrks_low, &b_PU_ntrks_low);
   fChain->SetBranchAddress("PU_ntrks_high", PU_ntrks_high, &b_PU_ntrks_high);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_et", Jet_et, &b_Jet_et);
   fChain->SetBranchAddress("Jet_residual", Jet_residual, &b_Jet_residual);
   fChain->SetBranchAddress("Jet_jes", Jet_jes, &b_Jet_jes);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_ntracks", Jet_ntracks, &b_Jet_ntracks);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   fChain->SetBranchAddress("Jet_Ip1N", Jet_Ip1N, &b_Jet_Ip1N);
   fChain->SetBranchAddress("Jet_Ip1P", Jet_Ip1P, &b_Jet_Ip1P);
   fChain->SetBranchAddress("Jet_Ip2N", Jet_Ip2N, &b_Jet_Ip2N);
   fChain->SetBranchAddress("Jet_Ip2P", Jet_Ip2P, &b_Jet_Ip2P);
   fChain->SetBranchAddress("Jet_Ip3N", Jet_Ip3N, &b_Jet_Ip3N);
   fChain->SetBranchAddress("Jet_Ip3P", Jet_Ip3P, &b_Jet_Ip3P);
   fChain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN, &b_Jet_ProbaN);
   fChain->SetBranchAddress("Jet_ProbaP", Jet_ProbaP, &b_Jet_ProbaP);
   fChain->SetBranchAddress("Jet_Proba", Jet_Proba, &b_Jet_Proba);
   fChain->SetBranchAddress("Jet_BprobN", Jet_BprobN, &b_Jet_BprobN);
   fChain->SetBranchAddress("Jet_Bprob", Jet_Bprob, &b_Jet_Bprob);
   fChain->SetBranchAddress("Jet_SvxN", Jet_SvxN, &b_Jet_SvxN);
   fChain->SetBranchAddress("Jet_Svx", Jet_Svx, &b_Jet_Svx);
   fChain->SetBranchAddress("Jet_SvxNTracks", Jet_SvxNTracks, &b_Jet_SvxNTracks);
   fChain->SetBranchAddress("Jet_SvxTracks", Jet_SvxTracks, &b_Jet_SvxTracks);
   fChain->SetBranchAddress("Jet_SvxNHP", Jet_SvxNHP, &b_Jet_SvxNHP);
   fChain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP, &b_Jet_SvxHP);
   fChain->SetBranchAddress("Jet_SvxMass", Jet_SvxMass, &b_Jet_SvxMass);
   fChain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN, &b_Jet_CombSvxN);
   fChain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx, &b_Jet_CombSvx);
   fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   fChain->SetBranchAddress("Jet_nFirstTrack", Jet_nFirstTrack, &b_Jet_nFirstTrack);
   fChain->SetBranchAddress("Jet_nLastTrack", Jet_nLastTrack, &b_Jet_nLastTrack);
   fChain->SetBranchAddress("Jet_nFirstSV", Jet_nFirstSV, &b_Jet_nFirstSV);
   fChain->SetBranchAddress("Jet_nLastSV", Jet_nLastSV, &b_Jet_nLastSV);
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
   fChain->SetBranchAddress("Muon_ptrel", Muon_ptrel, &b_Muon_ptrel);
   fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_hist", Muon_hist, &b_Muon_hist);
   fChain->SetBranchAddress("nCFromGSplit", &nCFromGSplit, &b_nCFromGSplit);
   fChain->SetBranchAddress("cFromGSplit_pT", cFromGSplit_pT, &b_cFromGSplit_pT);
   fChain->SetBranchAddress("cFromGSplit_eta", cFromGSplit_eta, &b_cFromGSplit_eta);
   fChain->SetBranchAddress("cFromGSplit_phi", cFromGSplit_phi, &b_cFromGSplit_phi);
   fChain->SetBranchAddress("nBFromGSplit", &nBFromGSplit, &b_nBFromGSplit);
   fChain->SetBranchAddress("bFromGSplit_pT", bFromGSplit_pT, &b_bFromGSplit_pT);
   fChain->SetBranchAddress("bFromGSplit_eta", bFromGSplit_eta, &b_bFromGSplit_eta);
   fChain->SetBranchAddress("bFromGSplit_phi", bFromGSplit_phi, &b_bFromGSplit_phi);
   Notify();
}

Bool_t JetTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetTree_cxx

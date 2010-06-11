//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 10 18:18:33 2010 by ROOT version 5.20/00
// from TTree ttree/ttree
// found on file: JetTree_All.root
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
   Float_t         Jet_pt[50];   //[nJet]
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
   Float_t         Jet_SvxN[50];   //[nJet]
   Float_t         Jet_Svx[50];   //[nJet]
   Int_t           Jet_SvxNTracks[50];   //[nJet]
   Int_t           Jet_SvxTracks[50];   //[nJet]
   Float_t         Jet_SvxNHP[50];   //[nJet]
   Float_t         Jet_SvxHP[50];   //[nJet]
   Float_t         Jet_CombSvxN[50];   //[nJet]
   Float_t         Jet_CombSvx[50];   //[nJet]
   Float_t         Jet_SoftMuN[50];   //[nJet]
   Float_t         Jet_SoftMu[50];   //[nJet]
   Int_t           Jet_hist1[50];   //[nJet]
   Int_t           Jet_hist2[50];   //[nJet]
   Int_t           Jet_hist3[50];   //[nJet]
   Int_t           Jet_histJet[50];   //[nJet]
   Int_t           Jet_histSvx[50];   //[nJet]
   Int_t           Jet_histMuon[50];   //[nJet]
   Int_t           Jet_mu_nHit[50];   //[nJet]
   Float_t         Jet_mu_chi2[50];   //[nJet]
   Float_t         Jet_mu_pt[50];   //[nJet]
   Float_t         Jet_mu_ptrel[50];   //[nJet]
   Int_t           Jet_nFirstTrack[50];   //[nJet]
   Int_t           Jet_nLastTrack[50];   //[nJet]
   Int_t           Jet_nFirstSV[50];   //[nJet]
   Int_t           Jet_nLastSV[50];   //[nJet]

   // List of branches
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Jet_pt;   //!
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
   TBranch        *b_Jet_SvxN;   //!
   TBranch        *b_Jet_Svx;   //!
   TBranch        *b_Jet_SvxNTracks;   //!
   TBranch        *b_Jet_SvxTracks;   //!
   TBranch        *b_Jet_SvxNHP;   //!
   TBranch        *b_Jet_SvxHP;   //!
   TBranch        *b_Jet_CombSvxN;   //!
   TBranch        *b_Jet_CombSvx;   //!
   TBranch        *b_Jet_SoftMuN;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_hist1;   //!
   TBranch        *b_Jet_hist2;   //!
   TBranch        *b_Jet_hist3;   //!
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //!
   TBranch        *b_Jet_histMuon;   //!
   TBranch        *b_Jet_mu_nHit;   //!
   TBranch        *b_Jet_mu_chi2;   //!
   TBranch        *b_Jet_mu_pt;   //!
   TBranch        *b_Jet_mu_ptrel;   //!
   TBranch        *b_Jet_nFirstTrack;   //!
   TBranch        *b_Jet_nLastTrack;   //!
   TBranch        *b_Jet_nFirstSV;   //!
   TBranch        *b_Jet_nLastSV;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("JetTree_All.root");
      if (!f) {
         f = new TFile("JetTree_All.root");
         f->cd("JetTree_All.root:/mistag");
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
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
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
   fChain->SetBranchAddress("Jet_SvxN", Jet_SvxN, &b_Jet_SvxN);
   fChain->SetBranchAddress("Jet_Svx", Jet_Svx, &b_Jet_Svx);
   fChain->SetBranchAddress("Jet_SvxNTracks", Jet_SvxNTracks, &b_Jet_SvxNTracks);
   fChain->SetBranchAddress("Jet_SvxTracks", Jet_SvxTracks, &b_Jet_SvxTracks);
   fChain->SetBranchAddress("Jet_SvxNHP", Jet_SvxNHP, &b_Jet_SvxNHP);
   fChain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP, &b_Jet_SvxHP);
   fChain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN, &b_Jet_CombSvxN);
   fChain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx, &b_Jet_CombSvx);
   fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   fChain->SetBranchAddress("Jet_histMuon", Jet_histMuon, &b_Jet_histMuon);
   fChain->SetBranchAddress("Jet_mu_nHit", Jet_mu_nHit, &b_Jet_mu_nHit);
   fChain->SetBranchAddress("Jet_mu_chi2", Jet_mu_chi2, &b_Jet_mu_chi2);
   fChain->SetBranchAddress("Jet_mu_pt", Jet_mu_pt, &b_Jet_mu_pt);
   fChain->SetBranchAddress("Jet_mu_ptrel", Jet_mu_ptrel, &b_Jet_mu_ptrel);
   fChain->SetBranchAddress("Jet_nFirstTrack", Jet_nFirstTrack, &b_Jet_nFirstTrack);
   fChain->SetBranchAddress("Jet_nLastTrack", Jet_nLastTrack, &b_Jet_nLastTrack);
   fChain->SetBranchAddress("Jet_nFirstSV", Jet_nFirstSV, &b_Jet_nFirstSV);
   fChain->SetBranchAddress("Jet_nLastSV", Jet_nLastSV, &b_Jet_nLastSV);
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

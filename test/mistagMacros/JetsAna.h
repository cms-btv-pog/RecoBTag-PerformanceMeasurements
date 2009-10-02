//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr  6 17:30:18 2009 by ROOT version 5.20/00
// from TTree Jets/All Jets
// found on file: results_80_120.root
//////////////////////////////////////////////////////////

#ifndef Jets_h
#define Jets_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Jets {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         Event_Njets;
   Float_t         Event_Trig;
   Float_t         Jet_Number;
   Float_t         Jet_Flavour;
   Float_t         Jet_Ntracks;
   Float_t         Jet_Pt;
   Float_t         Jet_Eta;
   Float_t         Jet_Phi;
   Float_t         TagN_TC1;
   Float_t         TagP_TC1;
   Float_t         TagN_TC2;
   Float_t         TagP_TC2;
   Float_t         TagN_TC3;
   Float_t         TagP_TC3;
   Float_t         TagN_JP;
   Float_t         TagP_JP;
   Float_t         Tag_JP;
   Float_t         TagN_SSV;
   Float_t         TagP_SSV;
   Float_t         TagN_CSV;
   Float_t         TagP_CSV;
   Float_t         TagN_MU;
   Float_t         TagP_MU;
   Float_t         Cat_TC1;
   Float_t         Cat_TC2;
   Float_t         Cat_TC3;
   Float_t         Cat_JP;
   Float_t         Cat_SSV;
   Float_t         Cat_MU;
   Float_t         Muon_hits;
   Float_t         Muon_chi2;
   Float_t         Muon_d0;
   Float_t         Muon_Pt;
   Float_t         Muon_PtRel;

   // List of branches
   TBranch        *b_Event_Njets;   //!
   TBranch        *b_Event_Trig;   //!
   TBranch        *b_Jet_Number;   //!
   TBranch        *b_Jet_Flavour;   //!
   TBranch        *b_Jet_Ntracks;   //!
   TBranch        *b_Jet_Pt;   //!
   TBranch        *b_Jet_Eta;   //!
   TBranch        *b_Jet_Phi;   //!
   TBranch        *b_TagN_TC1;   //!
   TBranch        *b_TagP_TC1;   //!
   TBranch        *b_TagN_TC2;   //!
   TBranch        *b_TagP_TC2;   //!
   TBranch        *b_TagN_TC3;   //!
   TBranch        *b_TagP_TC3;   //!
   TBranch        *b_TagN_JP;   //!
   TBranch        *b_TagP_JP;   //!
   TBranch        *b_Tag_JP;   //!
   TBranch        *b_TagN_SSV;   //!
   TBranch        *b_TagP_SSV;   //!
   TBranch        *b_TagN_CSV;   //!
   TBranch        *b_TagP_CSV;   //!
   TBranch        *b_TagN_MU;   //!
   TBranch        *b_TagP_MU;   //!
   TBranch        *b_Cat_TC1;   //!
   TBranch        *b_Cat_TC2;   //!
   TBranch        *b_Cat_TC3;   //!
   TBranch        *b_Cat_JP;   //!
   TBranch        *b_Cat_SSV;   //!
   TBranch        *b_Cat_MU;   //!
   TBranch        *b_Muon_hits;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_d0;   //!
   TBranch        *b_Muon_Pt;   //!
   TBranch        *b_Muon_PtRel;   //!

   Jets(TTree *tree=0);
   virtual ~Jets();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//$$
//$$   void   Loop();
   void   Loop(Int_t    aNN = 0, 
               Double_t aTagCut = 0.01, 
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

#ifdef Jets_cxx
Jets::Jets(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results_80_120.root");
      if (!f) {
         f = new TFile("results_80_120.root");
      }
      tree = (TTree*)gDirectory->Get("Jets");

   }
   Init(tree);
}

Jets::~Jets()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Jets::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Jets::LoadTree(Long64_t entry)
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

void Jets::Init(TTree *tree)
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

   fChain->SetBranchAddress("Event_Njets", &Event_Njets, &b_Event_Njets);
   fChain->SetBranchAddress("Event_Trig", &Event_Trig, &b_Event_Trig);
   fChain->SetBranchAddress("Jet_Number", &Jet_Number, &b_Jet_Number);
   fChain->SetBranchAddress("Jet_Flavour", &Jet_Flavour, &b_Jet_Flavour);
   fChain->SetBranchAddress("Jet_Ntracks", &Jet_Ntracks, &b_Jet_Ntracks);
   fChain->SetBranchAddress("Jet_Pt", &Jet_Pt, &b_Jet_Pt);
   fChain->SetBranchAddress("Jet_Eta", &Jet_Eta, &b_Jet_Eta);
   fChain->SetBranchAddress("Jet_Phi", &Jet_Phi, &b_Jet_Phi);
   fChain->SetBranchAddress("TagN_TC1", &TagN_TC1, &b_TagN_TC1);
   fChain->SetBranchAddress("TagP_TC1", &TagP_TC1, &b_TagP_TC1);
   fChain->SetBranchAddress("TagN_TC2", &TagN_TC2, &b_TagN_TC2);
   fChain->SetBranchAddress("TagP_TC2", &TagP_TC2, &b_TagP_TC2);
   fChain->SetBranchAddress("TagN_TC3", &TagN_TC3, &b_TagN_TC3);
   fChain->SetBranchAddress("TagP_TC3", &TagP_TC3, &b_TagP_TC3);
   fChain->SetBranchAddress("TagN_JP", &TagN_JP, &b_TagN_JP);
   fChain->SetBranchAddress("TagP_JP", &TagP_JP, &b_TagP_JP);
   fChain->SetBranchAddress("Tag_JP", &Tag_JP, &b_Tag_JP);
   fChain->SetBranchAddress("TagN_SSV", &TagN_SSV, &b_TagN_SSV);
   fChain->SetBranchAddress("TagP_SSV", &TagP_SSV, &b_TagP_SSV);
   fChain->SetBranchAddress("TagN_CSV", &TagN_CSV, &b_TagN_CSV);
   fChain->SetBranchAddress("TagP_CSV", &TagP_CSV, &b_TagP_CSV);
   fChain->SetBranchAddress("TagN_MU", &TagN_MU, &b_TagN_MU);
   fChain->SetBranchAddress("TagP_MU", &TagP_MU, &b_TagP_MU);
   fChain->SetBranchAddress("Cat_TC1", &Cat_TC1, &b_Cat_TC1);
   fChain->SetBranchAddress("Cat_TC2", &Cat_TC2, &b_Cat_TC2);
   fChain->SetBranchAddress("Cat_TC3", &Cat_TC3, &b_Cat_TC3);
   fChain->SetBranchAddress("Cat_JP", &Cat_JP, &b_Cat_JP);
   fChain->SetBranchAddress("Cat_SSV", &Cat_SSV, &b_Cat_SSV);
   fChain->SetBranchAddress("Cat_MU", &Cat_MU, &b_Cat_MU);
   fChain->SetBranchAddress("Muon_hits", &Muon_hits, &b_Muon_hits);
   fChain->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
   fChain->SetBranchAddress("Muon_d0", &Muon_d0, &b_Muon_d0);
   fChain->SetBranchAddress("Muon_Pt", &Muon_Pt, &b_Muon_Pt);
   fChain->SetBranchAddress("Muon_PtRel", &Muon_PtRel, &b_Muon_PtRel);
   Notify();
}

Bool_t Jets::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Jets::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Jets::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Jets_cxx

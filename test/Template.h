#ifndef Template_h
#define Template_h
/** \class Template
 *  
 * Analyze ROOT files produced by analyzer and create plots
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: Template.h,v 1.1 2007/09/24 18:26:48 yumiceva Exp $
 *
 */

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>

#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"


class Template {

  public:
	
	Template(TString filename);
	virtual ~Template();
	BTagEvent *fS8evt;
	
	TChain            *fChain;
	Int_t             fCurrent;
	TFile            *ffile;
	TFile            *foutfile;
	
	virtual Int_t  GetEntry(Long64_t entry);
	virtual void     Init();
	virtual void     Loop();
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Show(Long64_t entry = -1);
	virtual Int_t    Cut(Long64_t entry);
	virtual Bool_t   Notify();
	void Add(TString filename);
	
};

#endif

#ifdef Template_cxx
Template::Template(TString filename)
{
			
	fChain = new TChain("summary");
	fChain->Add(filename);
	
	fS8evt = new BTagEvent();
	
	//cout << " file loaded and objects created" << endl;
	Init();
	
}

Template::~Template()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Template::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t Template::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Template::Add(TString filename)
{
	fChain->Add(filename);
}

void Template::Init()
{

	//if (!tree) return;

	//fChain = tree;
	fCurrent = -1;
	
	fChain->SetBranchAddress("s8.",&fS8evt);
	
}

Bool_t Template::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Template::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Template::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Template_cxx

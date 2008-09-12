/** \class BasePlotter
 *  
 * Wrapper class to give access to the BTagEvent content
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 * \author Victor E. Bazterra, UIC
 *
 * \version $Id: BasePlotter.h,v 1.1 2008/05/17 11:06:19 bazterra Exp $
 *
 */

#ifndef BasePlotter_h
#define BasePlotter_h

#include "TChain.h"
#include "TFile.h"
#include "TObject.h"
#include "TString.h"

#ifdef NOSCRAMV
#include "BTagEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"
#endif

class BasePlotter : public TObject {

  public:

   BasePlotter(TString filename="")
   {
     fChain_ = new TChain("summary");
     if ( filename != "" ) fChain_->Add(filename);
     event_ = new BTagEvent();
   }
    
   ~BasePlotter()
   {
     if (!fChain_) return;
     delete fChain_->GetCurrentFile();
   }

   void Add(TString filename)
   {
     fChain_->Add(filename);
   }
    
   void Loop(Long64_t max = 0);

   virtual void Write() {}  

   ClassDef(BasePlotter,1);

  protected:
  
   virtual void Book() = 0;  

   virtual void Fill(BTagEvent *) = 0;
    
  private:

    BTagEvent * event_;

    TChain * fChain_;

    Int_t fCurrent_;

    Int_t GetEntry(Long64_t entry)
    {
      if (!fChain_) return 0;
      return fChain_->GetEntry(entry);
    }

    // Set the environment to read one entry
	Long64_t LoadTree(Long64_t entry);
};


#endif 

#ifndef BTagHistograms_h
#define BTagHistograms_h

/**_________________________________________________________________
   class:   BTagHistograms.h

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagHistograms.h,v 1.2 2008/06/11 08:46:58 tboccali Exp $

________________________________________________________________**/


#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"

#include <map>

class BTagHistograms
{

public:

    BTagHistograms();
    ~BTagHistograms();

    void Init(TString type, TString suffix1="", TString suffix2="");
    void Fill1d(TString name, Double_t x , Double_t weight =1);
    void Fill2d(TString name, Double_t x, Double_t y, Double_t weight=1);
    void Save();
    void SaveToFile(TString filename);
    void Fit(TString name, Double_t mean);
    void DeleteHisto()
    {
        for (std::map<TString,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih)
        {
            TH1 *htemp = ih->second;
            delete htemp;
        }
        for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
        {
            TH2 *htemp = ih->second;
            delete htemp;
        }
    };

private:

    std::map<TString, TCanvas*> cv_map;
    std::map<TString, TH1*> h1;
    std::map<TString, TH2*> h2;
    TFile            *ffile;
    TFile            *foutfile;

};

#endif

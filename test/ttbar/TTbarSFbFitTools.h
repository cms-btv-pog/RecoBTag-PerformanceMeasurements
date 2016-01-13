#ifndef _TTbarSFbFitTools_h_
#define _TTbarSFbFitTools_h_

#include <vector>

#include "TH1F.h"
#include "TObjArray.h"
#include "TString.h"

struct TTbarFracFitterResult_t
{
  float eff,effUnc,sf,sfUnc,effExp,effExpUnc;
  int minuitStatus;
};

class TTbarFracFitter
{
 public:
  TTbarFracFitter();
  TTbarFracFitterResult_t fit(TObjArray &fracTempl,TH1F *data,Int_t idxOfInterest=0,TString saveResultIn="");
  TTbarFracFitterResult_t fit(TObjArray &passTemplates, TH1F *passDataH,
			      TObjArray &failTemplates, TH1F *failDataH,
			      Int_t idxOfInterest=0,TString saveResultIn="",
			      Float_t lumi=2.444);
  ~TTbarFracFitter();
};

#endif

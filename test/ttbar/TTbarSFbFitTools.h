#ifndef _TTbarSFbFitTools_h_
#define _TTbarSFbFitTools_h_

#include <vector>

#include "TH1F.h"
#include "TObjArray.h"
#include "TString.h"

struct TTbarFracFitterResult_t
{
  float nExp,nExpUnc,nObs,nObsUnc,sf,sfUnc;
  int minuitStatus;
};

class TTbarFracFitter
{
 public:
  TTbarFracFitter();
  TTbarFracFitterResult_t fit(TObjArray &fracTempl,TH1F *data,Int_t idxOfInterest=0,TString saveResultIn="");
  ~TTbarFracFitter();
};

#endif

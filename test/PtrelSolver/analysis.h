#ifndef _analysis_h_
#define _analysis_h_

#include "TCut.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TChain.h"
#include "TString.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TFile.h"
#include "TGraphErrors.h"

Double_t    pdf(Double_t *xx, Double_t *par);
Double_t    pdf1(Double_t *xx, Double_t *par);
Double_t    combined_pdf(Double_t *xx, Double_t *par);



void        formatHist1(TH1 *hh, const char *xtitle, const char *ytitle);
const char *saveAsEPS(TCanvas *cc, TH1 *hist, const char *opt = "", TLegend *leg = 0);
double      weightSum(TH1F *temp, TProfile *w);
void        Divide(TH1F *result, TH1 *h1, TH1 *h2);
TH1F       *Divide(TH1 *h1, TH1 *h2);

#endif

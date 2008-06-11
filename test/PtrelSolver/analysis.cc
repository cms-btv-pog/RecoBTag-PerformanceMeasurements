#include "analysis.h"
#include <math.h>

/**************************************************************************
 *
 * templates: d * [x^a * e( bx^2) +c]
 * d is the normalization factor. 
 *
 **************************************************************************/
Double_t pdf(Double_t *xx, Double_t *par) {


  Double_t x = xx[0];

  return par[3] * (pow(x, par[0]) * exp(par[1] * x * x) + par[2]);
}


// this is for c flavor
Double_t pdf1(Double_t *xx, Double_t *par) {


  Double_t x = xx[0];

  return par[4] * (pow(x, par[0]) * exp(par[1] * pow(x, par[3])) + par[2]);
}



// make combined functions, assuming that each has parameter size 4
Double_t combined_pdf(Double_t *xx, Double_t *par) {

  Double_t f= 0;

  f += pdf1(xx, &par[0]); 
  f += pdf1(xx, &par[5]); 

  return f;
}
/**************************************************************************/




void formatHist1(TH1 *hh, const char *xtitle, const char *ytitle) {

  if (!hh) return;
  hh->UseCurrentStyle();
  hh->GetXaxis()->SetTitle(xtitle);
  hh->GetXaxis()->SetNdivisions(505);

  hh->GetYaxis()->SetTitle(ytitle);
  hh->GetYaxis()->SetTitleOffset(1.4);
  hh->GetYaxis()->SetNdivisions(510);

  if (hh->GetZaxis()) {
    hh->GetZaxis()->SetNdivisions(510);
  }
}


const char *saveAsEPS(TCanvas *cc, TH1 *hist, const char *opt, TLegend *leg) {

  if (!cc || !hist) return 0;
  cc->cd();

  TString *epsname = new  TString(hist->GetName());
  epsname->Append(".eps");
  hist->Draw(opt);

  if (leg) leg->Draw();
  cc->SaveAs(epsname->Data());
  return epsname->Data();
}


double weightSum(TH1F *temp, TProfile *w) {
  if (!temp || !w) return 0;
  double value = 0;
  int nbins = temp->GetNbinsX(); 
 for (int jj = 1; jj <= nbins; jj++) {

    value += temp->GetBinContent(jj) * w->GetBinContent(jj);
  }

  return value;
}


void Divide(TH1F *result, TH1 *h1, TH1 *h2) {
  if (!h1 || !h2|| !result) return;

  int numbins = h1->GetNbinsX();
  for (int ii = 1; ii <= numbins; ii++) {

    float val = 1.0, err = 0;
    float n1 = h1->GetBinContent(ii);
    float n2 = h2->GetBinContent(ii);

    if (h2->GetBinContent(ii)) {


      val = n1/n2;
      err = sqrt( pow(h1->GetBinError(ii), 2) /pow(n2, 2) + pow(n1, 2)/pow(n2, 4) * pow(h2->GetBinError(ii), 2)  );

    }
    result->SetBinContent(ii, val);
    result->SetBinError(ii, err);      
  }
}




TH1F *Divide(TH1 *h1, TH1 *h2) {

  if (!h1 || !h2) return 0;

  TH1F * result = new TH1F( *((TH1F *)h1) );
  TString name(h1->GetName()); name.Append("_ratio");
  result->SetName(name.Data());


  int numbins = h1->GetNbinsX();
  for (int ii = 1; ii <= numbins; ii++) {

    float val = 1.0, err = 0;
    float n1 = h1->GetBinContent(ii);
    float n2 = h2->GetBinContent(ii);

    if (h2->GetBinContent(ii)) {


      val = n1/n2;
      err = sqrt( pow(h1->GetBinError(ii), 2) /pow(n2, 2) + pow(n1, 2)/pow(n2, 4) * pow(h2->GetBinError(ii), 2)  );

    }
    result->SetBinContent(ii, val);
    result->SetBinError(ii, err);      
  }

  return result;
}

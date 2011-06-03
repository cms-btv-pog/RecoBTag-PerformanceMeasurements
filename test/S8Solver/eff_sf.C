#include <vector>
#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include "TLatex.h"
#include "TFile.h"

void eff_sf(const char *filename="s8.root", const char * OP= "", const char *datalegend="Data") {

  //gROOT->SetStyle("Plain");
  //gStyle->SetFillColor(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetPadBorderSize(0);  
  //gStyle->SetOptStat(1111111);

  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  Double_t xp,yp;
  Double_t xpe,ype;
  
  
  TFile* file = TFile::Open(filename);
  
  if (  ! file->IsOpen() ) {
    cerr << "Could not open the input file "  << endl;
    
  }
  
  
  
  gDirectory->cd("mcefficiency");
  
  
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
  TPad *d1, *d2;
  
  c1->Divide(1,2,0,0);
  d1 = (TPad*)c1->GetPad(1);
  d1->SetPad(0.01,0.35,0.95,0.95);
  d2 = (TPad*)c1->GetPad(2);
  d2->SetPad(0.01,0.02,0.95,0.35);
	
	
  d1->cd();
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.19);
  gPad->SetGrid();

  g1 = (TGraphErrors*) eff_tag_b;

  Double_t effmc[20];
  Double_t effdata[20];
  Double_t effmc_xerr[20];
  Double_t effdata_xerr[20];
  Double_t effmc_yerr[20];
  Double_t effdata_yerr[20];
  
  Double_t pt[20];
  
  Double_t sf[20] ;
  
  Double_t sf_y_error[20] ;
  
  Double_t sf_x_error[20];
	
	
	
	
  for(int i = 0;  g1->GetN() > i; ++i) {
    g1->GetPoint(i,xp,yp);
    xpe=g1->GetErrorX(i);
    ype=g1->GetErrorY(i);
    pt[i]=xp;
    effmc[i]=yp;
    //  cout << i << " " << yp << endl;
    effmc_xerr[i]=xpe;
    effmc_yerr[i]=ype;
  }
  
  g1->SetTitle("Efficiency");
  g1->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  TString tmplabel("b-tag Efficiency #epsilon_{b}");
  tmplabel = TString(OP) + TString(" ")+ tmplabel;
  g1->GetYaxis()->SetTitle(tmplabel);
  g1->GetYaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->SetLabelSize(0.03);
  g1->GetYaxis()->SetTitleOffset(0.9);
  //  g1->GetYaxis()->CenterTitle(true);
  g1->SetMarkerStyle(22);

  int n = g1->GetN() ;
  double low = 0.22; //(effmc[0]- 0.1);
  double high = 0.88; //(1.1* effmc[n-1]);
  std::cout << "low/high " << low <<"/" << high <<std::endl;
  g1->GetYaxis()->SetRangeUser(low,high);
  gPad->SetTickx();
  g1->Draw("AP");
  




  gDirectory->cd("../s8efficiency");
  
  g2 = (TGraphErrors*) eff_tag_b;
  
  
  for(int i = 0;  g2->GetN() > i; ++i) {
    g2->GetPoint(i,xp,yp);
    xpe=g2->GetErrorX(i);
    ype=g2->GetErrorY(i);
    pt[i]=xp;
    effdata[i]=yp;
    // cout << i << " " << yp << endl;
    effdata_xerr[i]=xpe;
    effdata_yerr[i]=ype;
  }
  
  
  int npoints=g2->GetN();
  g2->SetMarkerStyle(23);
  g2->Draw("P");
  
  std::ostringstream title;
  title << "#splitline{CMS Preliminary 2011}{ at #sqrt{s} = 7 TeV}";
  //title << "16 pb^{-1}";
  //title << " at #sqrt{s} = 7 TeV}";
  TLatex *label = new TLatex(3.570061, 23.08044, title.str().c_str());
  label->SetNDC();
  label->SetTextAlign(13);
  label->SetX(0.5);
  label->SetY(0.970);
  double yl = low - 0.05 ;
  double yh = yl + 0.15;
  TLegend* leg = new TLegend(0.4,yl,0.7,yh); //use inspector
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g1, "Simulation", "P");
  leg->AddEntry(g2, datalegend, "P");
  leg->Draw();
  label->Draw();
  
  
  d2->cd();
  gPad->SetLeftMargin(0.19);
  gPad->SetTopMargin(0);
  gPad->SetRightMargin(0.04);
  gPad->SetFrameBorderSize(0);
  gPad->SetBottomMargin(0.3);
  gPad->SetTickx();
  gPad->SetRightMargin(0.04); 
  gPad->SetGrid();
  
  
  for(int i = 0;  npoints > i; ++i) {
    sf[i]=effdata[i]/effmc[i];
    sf_x_error[i]=effdata_xerr[i];
    sf_y_error[i]=((effdata_yerr[i]*effdata_yerr[i])/(effdata[i]*effdata[i]))+
      ((effmc_yerr[i]*effmc_yerr[i])/(effmc[i]*effmc[i]));
    sf_y_error[i]=sqrt(sf_y_error[i])*sf[i];
    cout << effdata[i] << " " << effmc[i] << " " << sf[i] << 
      " " << sf_y_error[i] << " " << sf_y_error[i]/sf[i] << endl;
    
  }
  
  
  TGraphErrors* e0 = new TGraphErrors(npoints, pt, sf, sf_x_error, sf_y_error);
  //e0->SetMarkerColor(7);
  e0->SetMarkerStyle(20);
  e0->SetMarkerSize(1.3);
  e0->SetTitle("Scale Factors");
  e0->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  e0->GetXaxis()->SetTitleOffset(1.2);

  //e0->GetYaxis()->SetRangeUser(0, 1);
  //e0->GetXaxis()->CenterTitle(kTRUE);
  e0->GetYaxis()->SetTitle("Data/MC Scale Factor SF_{b}");
  e0->GetYaxis()->SetTitleSize(0.06);
  e0->GetYaxis()->SetTitleOffset(0.9);
  e0->GetYaxis()->CenterTitle(true);
  e0->GetXaxis()->SetLabelSize(0.07); 
  e0->GetYaxis()->SetLabelSize(0.06); 
  e0->GetXaxis()->SetTitleSize(0.07);
  
  //    e0->GetXaxis()->SetLabelOffset(0.04);
  
  e0->GetYaxis()->SetRangeUser(0.48,1.52);

  //e0->GetYaxis()->CenterTitle(kTRUE);
  e0->Draw("AP");
  
  std::ostringstream title;   
  title << "#splitline{CMS Preliminary 2011}{at #sqrt{s} = 7 TeV}";
  
  TLatex *label = new TLatex(3.570061, 23.08044, title.str().c_str());
  label->SetNDC();
  label->SetTextAlign(13);
  label->SetX(0.4);
  label->SetY(0.890);
  //label->Draw();
  c1->Update();
  
  c1->Print("eff_tag_b.pdf","pdf");
  

}

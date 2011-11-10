#include <string>
#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TROOT.h"
#include <iomanip>


void eff_sf(const char *filename="s8.root", const char * OP= "", const char *datalegend="Data", const char *ptrelfilename="") {

  //gROOT->SetStyle("Plain");
  //gStyle->SetFillColor(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetPadBorderSize(0);  
  //gStyle->SetOptStat(0);

  // use the TDR CMS style
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptFit(0);

  // define some colors
  int s8color = 2; // red
  int ptrelcolor = 4; // blue
  int s8marker_data = 20;
  int s8marker_mc = 24;
  int ptrelmarker_data = 21;
  int ptrelmarker_mc = 25;

  Double_t xp,yp;
  Double_t xpe,ype;
  
  bool HasPtrel = false;

  TFile *file, *ptrelfile;
  
  file = new TFile(filename);
  
  if ( TString(ptrelfilename) != "" ) 
    {
      ptrelfile = new TFile(ptrelfilename);
      HasPtrel = true;
    }

  if ( file->IsZombie() ) 
    {    
      cout << "Could not open input file " << filename  << " exit now. " << endl;
      return;
    }
  
  
  // get plots.
  file->cd();
  file->cd("mcefficiency");
  // get MC Efficiencies
  TGraphErrors *g1 = (TGraphErrors*) gDirectory->Get("eff_tag_b");

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
  //gPad->SetGrid();

  g1->SetTitle("Efficiency");
  g1->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  TString tmplabel("b-tag Efficiency #epsilon_{b}");
  tmplabel = TString(OP) + TString(" ")+ tmplabel;
  g1->GetYaxis()->SetTitle(tmplabel);
  g1->GetYaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->SetLabelSize(0.05);
  g1->GetYaxis()->SetTitleOffset(1.2);
  //  g1->GetYaxis()->CenterTitle(true);
  g1->SetMarkerStyle(s8marker_mc);

  int n = g1->GetN() ;
  double low = (effmc[0]- 0.1);
  double high = (1.1* effmc[n-1]);
  //std::cout << "low/high " << low <<"/" << high <<std::endl;
  g1->GetYaxis()->SetRangeUser(low,high);
  gPad->SetTickx();
  g1->SetMarkerColor(s8color);
  g1->SetLineColor(s8color);
  g1->Draw("AP");
  
  // get S8 efficiencies
  gDirectory->cd("../s8efficiency");
  
  TGraphErrors *g2 = (TGraphErrors*) gDirectory->Get("eff_tag_b");
   
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
  g2->SetMarkerStyle(s8marker_data);
  g2->SetMarkerColor(s8color);
  g2->SetLineColor(s8color);

  g2->Draw("P");
  
  std::string title = "CMS prelim. at #sqrt{s} = 7 TeV"; 
  TLatex *label = new TLatex(3.570061, 23.08044, title.c_str());
//  TLatex *label = new TLatex(3.570061, 23.08044, title.c_str());
  label->SetNDC();
  label->SetTextAlign(13);
  label->SetX(0.25);
  label->SetY(0.98);
 
  // get ptrel results if available
  TGraphErrors *ptrel_data_g;
  TGraphErrors *ptrel_mc_g;
  TGraphErrors *ptrel_sf_g;

  if (HasPtrel) 
    {
      ptrelfile->cd();
      //TH1F *htmp = (TH1F*) gDirectory->Get("Data_"+TString(OP)+"__KinWeights_PVWeighting_DW");
      ptrel_data_g = new TGraphErrors( (TH1*) gDirectory->Get("Data_"+TString(OP)+"__KinWeights_PVWeighting_80GeVLarge") );
      ptrel_mc_g   = new TGraphErrors( (TH1*) gDirectory->Get("MC_"+TString(OP)+"__KinWeights_PVWeighting_80GeVLarge") );
      ptrel_sf_g   = new TGraphErrors( (TH1*) gDirectory->Get("SF_"+TString(OP)+"__KinWeights_PVWeighting_80GeVLarge") );

      ptrel_data_g->SetLineColor(ptrelcolor);
      ptrel_data_g->SetMarkerColor(ptrelcolor);
      ptrel_mc_g->SetLineColor(ptrelcolor);
      ptrel_mc_g->SetMarkerColor(ptrelcolor);
      
      ptrel_data_g->SetMarkerStyle(ptrelmarker_data);
      ptrel_mc_g->SetMarkerStyle(ptrelmarker_mc);
      ptrel_data_g->SetMarkerSize(1.5);
      ptrel_mc_g->SetMarkerSize(1.5);

      ptrel_data_g->Draw("P");
      ptrel_mc_g->Draw("P");
      
    }

  double yl = low ;
  double yh = yl + 0.2;
  TLegend* leg = new TLegend(0.62,yl,0.92,yh); //use inspector
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  if (HasPtrel)
    {
      leg->AddEntry(g1, "Simulation System8", "P");
      leg->AddEntry(g2, TString(datalegend)+" System8", "P");
      leg->AddEntry(ptrel_mc_g, "Simulation PtRel", "P");
      leg->AddEntry(ptrel_data_g, TString(datalegend)+" Ptrel", "P");
    }
  else
    {
      leg->AddEntry(g1, "Simulation", "P");
      leg->AddEntry(g2, datalegend, "P");
    }

  leg->Draw();
  label->Draw();
  TPaveText p6(0.62,yh,0.92,yh+0.05,"brNDC");
  p6.SetBorderSize(0);
  p6.SetFillColor(kWhite);
  p6.SetTextSize(0.04);
  std::ostringstream os;
  os << " for 1.2 < #eta < 2.4 ";
  TString buffer(os.str()); 
  p6.AddText(buffer);
  p6.Draw();
  c1->Update();
 
  
  d2->cd();
  
  gPad->SetLeftMargin(0.19);
  gPad->SetTopMargin(0);
  gPad->SetRightMargin(0.04);
  gPad->SetFrameBorderSize(0);
  gPad->SetBottomMargin(0.3);
  gPad->SetTickx();
  gPad->SetRightMargin(0.04); 
  //gPad->SetGrid();
  
  
  for(int i = 0;  npoints > i; ++i) {
    sf[i]=effdata[i]/effmc[i];
    sf_x_error[i]=effdata_xerr[i];
    sf_y_error[i]=((effdata_yerr[i]*effdata_yerr[i])/(effdata[i]*effdata[i]))+
      ((effmc_yerr[i]*effmc_yerr[i])/(effmc[i]*effmc[i]));
    sf_y_error[i]=sqrt(sf_y_error[i])*sf[i];
    //cout << effdata[i] << " " << effmc[i] << " " << sf[i] << 
    //  " " << sf_y_error[i] << " " << sf_y_error[i]/sf[i] << endl;
    
  }
  
  
  TGraphErrors* e0 = new TGraphErrors(npoints, pt, sf, sf_x_error, sf_y_error);
  e0->SetMarkerColor(s8color);
  e0->SetLineColor(s8color);
  e0->SetMarkerStyle(s8marker_data);
  e0->SetMarkerSize(1.3);
  e0->SetTitle("Scale Factors");
  e0->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  e0->GetXaxis()->SetTitleOffset(1.2);
  e0->GetXaxis()->SetTitleSize(0.09);
  e0->GetXaxis()->SetLabelSize(0.09);
  //e0->GetYaxis()->SetRangeUser(0, 1);
  //e0->GetXaxis()->CenterTitle(kTRUE);
  e0->GetYaxis()->SetTitle("Data/Sim. SF_{b}");
  e0->GetYaxis()->SetTitleSize(0.09);
  e0->GetYaxis()->SetTitleOffset(0.7);
  e0->GetYaxis()->CenterTitle(true);
  e0->GetYaxis()->SetLabelSize(0.09); 
  e0->GetYaxis()->SetRangeUser(0.48,1.52);

  TMultiGraph *the_sf = new TMultiGraph();
  the_sf->Add(e0,"p");
  //e0->Draw("AP");
  
  // plot ptrel SF if available
  if (HasPtrel)
    {
      ptrel_sf_g->SetMarkerSize(1.3);
      ptrel_sf_g->SetMarkerStyle(ptrelmarker_data);
      ptrel_sf_g->SetLineColor(ptrelcolor);
      ptrel_sf_g->SetMarkerColor(ptrelcolor);

      //ptrel_sf_g->Draw("P");
      the_sf->Add(ptrel_sf_g,"p");

      // get last bin content
      Double_t last_bin_ptrel,last_bin_s8;
      Double_t tmpx;
      //ptrel_sf_g->Print("all");
      ptrel_sf_g->GetPoint(149,tmpx,last_bin_ptrel);
      e0->GetPoint(e0->GetN() - 1,tmpx,last_bin_s8);
      cout << last_bin_ptrel << endl;
      cout << last_bin_s8 << endl;
      cout << "Averaged SF for last bin: SF = " << (last_bin_ptrel+last_bin_s8)/2 << endl;

    }

  the_sf->Draw("a");
  the_sf->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  the_sf->GetXaxis()->SetTitleOffset(1.2);
  the_sf->GetXaxis()->SetTitleSize(0.09);
  the_sf->GetXaxis()->SetLabelSize(0.09);
  the_sf->GetYaxis()->SetTitle("Data/Sim. SF_{b}");
  the_sf->GetYaxis()->SetTitleSize(0.09);
  the_sf->GetYaxis()->SetTitleOffset(0.7);
  the_sf->GetYaxis()->CenterTitle(true);
  the_sf->GetYaxis()->SetLabelSize(0.09);
  the_sf->GetYaxis()->SetRangeUser(0.48,1.52);

  if (HasPtrel)
    {
      TLegend* leg2 = new TLegend(0.62,0.76,0.92,0.76+(yh-yl)*0.5); //use inspector
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);
      leg2->AddEntry(e0, "System8", "P");
      leg2->AddEntry(ptrel_sf_g, "PtRel", "P");
      leg2->Draw();
    }
 // Fit individually ptrel and s8 scale factors
  /*
  TF1 *s8_fit = new TF1("s8_fit","pol1",20,120);
  s8_fit->SetLineColor(s8color);
  TF1 *ptrel_fit =new TF1("ptrel_fit","pol1",20,120);
  ptrel_fit->SetLineColor(ptrelcolor);

  e0->Fit("s8_fit","RE");
  cout << "System8 fit chi2 = "<< s8_fit->GetChisquare() << endl;
  if (HasPtrel) 
    {
      ptrel_sf_g->Fit("ptrel_fit","RE");
      cout << "PtRel fit chi2 = "<< ptrel_fit->GetChisquare() << endl;
    }
  */
  // Fit simultaneously ptrel and s8
  TF1 *sf_fit = new TF1("sf_fit","pol0",20,240);
  sf_fit->SetLineColor(1);
  sf_fit->SetLineWidth(2);
  the_sf->Fit("sf_fit","RE");
  cout << "fit chi2 = " << sf_fit->GetChisquare() << endl;
  cout << "fit p0 = " << sf_fit->GetParameter(0) << endl;
  cout << "fit err = " << sqrt(2)*sf_fit->GetParError(0) << endl;
  TPaveText p5(0.25,0.86,0.45,0.96,"brNDC");
  p5.SetBorderSize(0);
  p5.SetFillColor(kWhite);
  p5.SetTextSize(0.06);
  std::ostringstream os;
  os << "<SF> = " << setprecision(2) <<sf_fit->GetParameter(0) << " #pm "<< setprecision(0) << sqrt(2)*sf_fit->GetParError(0) <<" (stat.)";
  TString buffer(os.str()); 
  p5.AddText(buffer);
  p5.Draw();
  c1->Update();

  // Print file
  TString pdfname = "s8_eff_sf_"+TString(OP)+".pdf";
  c1->Print(pdfname,"pdf");
  

}

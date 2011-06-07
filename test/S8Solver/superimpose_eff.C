#include <vector>
#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include "TLatex.h"
#include "TFile.h"

void superimpose_eff(const char *filenominal="s8.root",
		     const char *filelow = "s8.root",
		     const char *filemed = "s8.root",
		     const char *filehigh = "s8.root",
		     const char * OP= "", bool etaplot= false) {

  //gROOT->SetStyle("Plain");
  //gStyle->SetFillColor(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetPadBorderSize(0);  
  //gStyle->SetOptStat(1111111);

  // use the TDR CMS style
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  Int_t colors[4] = {1,2,6,4};

  Int_t markers[4] = { 22,24,25,26};

  TString labels[4] = { "all","1 #leq NPV #leq 3","4 #leq NPV #leq 6","7 #leq NPV "};
  
  if (etaplot)
    {
      labels[0] = "|#eta| #leq 2.4";
      labels[1] = "|#eta| #leq 2.1";
      labels[2] = "2.1 #leq |#eta| #leq 2.4";
    }

  cout << "Trying to open input files\n" << endl;

  TFile* tfilenominal = new TFile(filenominal);
  TFile* tfilelow     = new TFile(filelow);
  TFile* tfilemed     = new TFile(filemed);
  TFile* tfilehigh    = new TFile(filehigh);

  int nfiles = 0;

  if ( tfilenominal->IsZombie() || tfilelow->IsZombie() ) {
    cout << "Could not open the input file nominal or low"  << endl;
    return;
  }

  
  tfilenominal->cd(); 
  gDirectory->cd("mcefficiency");
    
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
  
  TGraphErrors *gnominal = (TGraphErrors*) eff_tag_b;

  TGraphErrors *glow;
  TGraphErrors *gmed;
  TGraphErrors *ghigh;

  gnominal->SetTitle("Efficiency");
  gnominal->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  TString tmplabel("b-tag Efficiency #epsilon_{b}");
  tmplabel = TString(OP) + TString(" ")+ tmplabel;
  gnominal->GetYaxis()->SetTitle(tmplabel);
  gnominal->GetYaxis()->SetTitleSize(0.05);
  gnominal->GetYaxis()->SetLabelSize(0.03);
  gnominal->GetYaxis()->SetTitleOffset(0.9);
  //  gnominal->GetYaxis()->CenterTitle(true);
  gnominal->SetMarkerStyle(22);

  int n = gnominal->GetN() ;
  double low = 0.18; //(effmc[0]- 0.1);
  double high = 0.98; //(1.1* effmc[n-1]);
  //std::cout << "low/high " << low <<"/" << high <<std::endl;
  gnominal->GetYaxis()->SetRangeUser(low,high);
  gnominal->SetMarkerStyle( markers[0] );
  gnominal->SetMarkerColor(colors[0]);
  gnominal->SetLineColor(colors[0]);

  gnominal->Draw("AP");
  

  if ( ! tfilelow->IsZombie() )
    {
      nfiles++;
      tfilelow->cd();
      gDirectory->cd("mcefficiency");
      glow = (TGraphErrors*) eff_tag_b;
      glow->SetMarkerStyle(markers[1]);
      glow->SetMarkerColor(colors[1]);
      glow->SetLineColor(colors[1]);
      glow->Draw("P");
      
    }

  if ( ! tfilemed->IsZombie() )
    {
      nfiles++;
      tfilemed->cd();
      gDirectory->cd("mcefficiency");
      gmed = (TGraphErrors*) eff_tag_b;
      gmed->SetMarkerStyle(markers[2]);
      gmed->SetMarkerColor(colors[2]);
      gmed->SetLineColor(colors[2]);

      gmed->Draw("P");
    }

  if ( ! tfilehigh->IsZombie() )
    {
      nfiles++;
      tfilehigh->cd();
      gDirectory->cd("mcefficiency");
      ghigh = (TGraphErrors*) eff_tag_b;
      ghigh->SetMarkerStyle(markers[3]);
      ghigh->SetMarkerColor(colors[3]);
      ghigh->SetLineColor(colors[3]);

      ghigh->Draw("P");
    }

  std::ostringstream title;
  title << "CMS Preliminary 2011 at #sqrt{s} = 7 TeV";

  TLatex *label = new TLatex(3.570061, 23.08044, title.str().c_str());
  label->SetNDC();
  label->SetTextAlign(13);
  label->SetX(0.2);
  label->SetY(0.990);
  double yl = low - 0.0 ;
  double yh = yl + 0.15;
  TLegend* leg = new TLegend(0.44,yl,0.74,yh); //use inspector                                                                                                                                                            
  leg->SetFillColor(0);
  leg->SetBorderSize(0);

  leg->AddEntry(gnominal,labels[0], "P");
  leg->AddEntry(glow,labels[1], "P");
  if (nfiles>=2)
    leg->AddEntry(gmed,labels[2], "P");
  if (nfiles>=3)
    leg->AddEntry(ghigh,labels[3], "P");

  leg->Draw();
  label->Draw();
  c1->SetGrid();

  c1->Update();
  
  TString suffix = "NPVregions";
  if (etaplot) suffix = "etaregions";

  TString pdfname = "superimpose_effb_"+suffix+"_"+TString(OP)+".pdf";
  //if TString(datalegend) !="Data" ) pdfname = "s8_eff_sf_"+TString(OP)+"_closure.pdf";
  c1->Print(pdfname,"pdf");
  
  tfilenominal->Close();
  tfilelow->Close();
  tfilemed->Close();
  tfilehigh->Close();
 
}

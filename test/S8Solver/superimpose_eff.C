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

void superimpose_eff(const char *filenominal="s8.root",
		     const char *filelow = "s8.root",
		     const char *filemed = "s8.root",
		     const char *filehigh = "s8.root",
		     const char * OP= "", bool etaplot= true) {

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

  if (etaplot)
    {
      labels[0] = "|#eta| #leq 2.4";
      labels[1] = "|#eta| #leq 1.2";
      labels[2] = "1.2 #leq |#eta| #leq 2.4";
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

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);

  
  tfilenominal->cd(); 
  gDirectory->cd("s8efficiency");
      
  TGraphErrors *gnominal = (TGraphErrors*) gDirectory->Get("eff_tag_b");

  gDirectory->cd("../mcefficiency");
  
  TGraphErrors *gnominalmc = (TGraphErrors*) gDirectory->Get("eff_tag_b");
  TGraphErrors *gnominalsf;

  TGraphErrors *glow, *glowmc, *glowsf;
  TGraphErrors *gmed, *gmedmc, *gmedsf;
  TGraphErrors *ghigh, *ghighmc, *ghighsf;

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

  //int n = gnominal->GetN() ;
  double low = 0.08; //(effmc[0]- 0.1);
  double high = 0.58; //(1.1* effmc[n-1]);
  //std::cout << "low/high " << low <<"/" << high <<std::endl;
  gnominal->GetYaxis()->SetRangeUser(low,high);
  gnominal->SetMarkerStyle( markers[0] );
  gnominal->SetMarkerColor(colors[0]);
  gnominal->SetLineColor(colors[0]);

  gnominal->Draw("AP");

  for(int i = 0;  gnominalmc->GetN() > i; ++i) {
    Double_t xp,yp,xpe,ype;
    gnominalmc->GetPoint(i,xp,yp);
    xpe=gnominalmc->GetErrorX(i);
    ype=gnominalmc->GetErrorY(i);
    pt[i]=xp;
    effmc[i]=yp;
    effmc_xerr[i]=xpe;
    effmc_yerr[i]=ype;
  }
  for(int i = 0;  gnominal->GetN() > i; ++i) {
    Double_t xp,yp,xpe,ype;
    gnominal->GetPoint(i,xp,yp);
    xpe=gnominal->GetErrorX(i);
    ype=gnominal->GetErrorY(i);
    pt[i]=xp;
    effdata[i]=yp;
    effdata_xerr[i]=xpe;
    effdata_yerr[i]=ype;
  }
  for(int i = 0;  gnominal->GetN() > i; ++i) {
    sf[i]=effdata[i]/effmc[i];
    sf_x_error[i]=effdata_xerr[i];
    sf_y_error[i]=((effdata_yerr[i]*effdata_yerr[i])/(effdata[i]*effdata[i]))+
      ((effmc_yerr[i]*effmc_yerr[i])/(effmc[i]*effmc[i]));
    sf_y_error[i]=sqrt(sf_y_error[i])*sf[i];
    //cout << effdata[i] << " " << effmc[i] << " " << sf[i] <<
    //  " " << sf_y_error[i] << " " << sf_y_error[i]/sf[i] << endl;

  }
  
  gnominalsf = new TGraphErrors(gnominal->GetN(), pt, sf, sf_x_error, sf_y_error);


  if ( ! tfilelow->IsZombie() )
    {
      nfiles++;
      tfilelow->cd();
      gDirectory->cd("s8efficiency");

      glow = (TGraphErrors*) gDirectory->Get("eff_tag_b");

      gDirectory->cd("../mcefficiency");
      glowmc = (TGraphErrors*) gDirectory->Get("eff_tag_b");
      
      glow->SetMarkerStyle(markers[1]);
      glow->SetMarkerColor(colors[1]);
      glow->SetLineColor(colors[1]);
      glow->Draw("P");

      for(int i = 0;  glowmc->GetN() > i; ++i) {
	Double_t xp,yp,xpe,ype;
	glowmc->GetPoint(i,xp,yp);
	xpe=glowmc->GetErrorX(i);
	ype=glowmc->GetErrorY(i);
	pt[i]=xp;
	effmc[i]=yp;
	effmc_xerr[i]=xpe;
	effmc_yerr[i]=ype;
      }
      for(int i = 0;  glow->GetN() > i; ++i) {
	Double_t xp,yp,xpe,ype;
	glow->GetPoint(i,xp,yp);
	xpe=glow->GetErrorX(i);
	ype=glow->GetErrorY(i);
	pt[i]=xp;
	effdata[i]=yp;
	effdata_xerr[i]=xpe;
	effdata_yerr[i]=ype;
      }
      for(int i = 0;  glow->GetN() > i; ++i) {
	sf[i]=effdata[i]/effmc[i];
	sf_x_error[i]=effdata_xerr[i];
	sf_y_error[i]=((effdata_yerr[i]*effdata_yerr[i])/(effdata[i]*effdata[i]))+
	  ((effmc_yerr[i]*effmc_yerr[i])/(effmc[i]*effmc[i]));
	sf_y_error[i]=sqrt(sf_y_error[i])*sf[i];
      }

      glowsf = new TGraphErrors(glow->GetN(), pt, sf, sf_x_error, sf_y_error);
      glowsf->SetMarkerStyle(markers[1]);
      glowsf->SetMarkerColor(colors[1]);
      glowsf->SetLineColor(colors[1]);

    }

  if ( ! tfilemed->IsZombie() )
    {
      nfiles++;
      tfilemed->cd();
      gDirectory->cd("s8efficiency");

      gmed = (TGraphErrors*) gDirectory->Get("eff_tag_b");

      gDirectory->cd("../mcefficiency");
      gmedmc = (TGraphErrors*) gDirectory->Get("eff_tag_b");

      gmed->SetMarkerStyle(markers[2]);
      gmed->SetMarkerColor(colors[2]);
      gmed->SetLineColor(colors[2]);

      gmed->Draw("P");

      for(int i = 0;  gmedmc->GetN() > i; ++i) {
	Double_t xp,yp,xpe,ype;
	gmedmc->GetPoint(i,xp,yp);
	xpe=gmedmc->GetErrorX(i);
	ype=gmedmc->GetErrorY(i);
	pt[i]=xp;
	effmc[i]=yp;
	effmc_xerr[i]=xpe;
	effmc_yerr[i]=ype;
      }
      for(int i = 0;  gmed->GetN() > i; ++i) {
	Double_t xp,yp,xpe,ype;
	gmed->GetPoint(i,xp,yp);
	xpe=gmed->GetErrorX(i);
	ype=gmed->GetErrorY(i);
	pt[i]=xp;
	effdata[i]=yp;
	effdata_xerr[i]=xpe;
	effdata_yerr[i]=ype;
      }
      for(int i = 0;  gmed->GetN() > i; ++i) {
	sf[i]=effdata[i]/effmc[i];
	sf_x_error[i]=effdata_xerr[i];
	sf_y_error[i]=((effdata_yerr[i]*effdata_yerr[i])/(effdata[i]*effdata[i]))+
	  ((effmc_yerr[i]*effmc_yerr[i])/(effmc[i]*effmc[i]));
	sf_y_error[i]=sqrt(sf_y_error[i])*sf[i];
      }

      gmedsf = new TGraphErrors(gmed->GetN(), pt, sf, sf_x_error, sf_y_error);
      gmedsf->SetMarkerStyle(markers[2]);
      gmedsf->SetMarkerColor(colors[2]);
      gmedsf->SetLineColor(colors[2]);

    }

  if ( ! tfilehigh->IsZombie() )
    {
      nfiles++;
      tfilehigh->cd();
      gDirectory->cd("s8efficiency");

      ghigh = (TGraphErrors*) gDirectory->Get("eff_tag_b");

      gDirectory->cd("../mcefficiency");
      ghighmc = (TGraphErrors*) gDirectory->Get("eff_tag_b");

      ghigh->SetMarkerStyle(markers[3]);
      ghigh->SetMarkerColor(colors[3]);
      ghigh->SetLineColor(colors[3]);

      ghigh->Draw("P");

      for(int i = 0;  ghighmc->GetN() > i; ++i) {
	Double_t xp,yp,xpe,ype;
	ghighmc->GetPoint(i,xp,yp);
	xpe=ghighmc->GetErrorX(i);
	ype=ghighmc->GetErrorY(i);
	pt[i]=xp;
	effmc[i]=yp;
	effmc_xerr[i]=xpe;
	effmc_yerr[i]=ype;
      }
      for(int i = 0;  ghigh->GetN() > i; ++i) {
	Double_t xp,yp,xpe,ype;
	ghigh->GetPoint(i,xp,yp);
	xpe=ghigh->GetErrorX(i);
	ype=ghigh->GetErrorY(i);
	pt[i]=xp;
	effdata[i]=yp;
	effdata_xerr[i]=xpe;
	effdata_yerr[i]=ype;
      }
      for(int i = 0;  ghigh->GetN() > i; ++i) {
	sf[i]=effdata[i]/effmc[i];
	sf_x_error[i]=effdata_xerr[i];
	sf_y_error[i]=((effdata_yerr[i]*effdata_yerr[i])/(effdata[i]*effdata[i]))+
	  ((effmc_yerr[i]*effmc_yerr[i])/(effmc[i]*effmc[i]));
	sf_y_error[i]=sqrt(sf_y_error[i])*sf[i];
      }

      ghighsf = new TGraphErrors(ghigh->GetN(), pt, sf, sf_x_error, sf_y_error);
      ghighsf->SetMarkerStyle(markers[3]);
      ghighsf->SetMarkerColor(colors[3]);
      ghighsf->SetLineColor(colors[3]);

    }

  //std::ostringstream title;
  string title ="CMS Preliminary 2011 at #sqrt{s} = 7 TeV";

  TLatex *label = new TLatex(3.570061, 23.08044, title.c_str());
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
  

  TCanvas* c2 = new TCanvas("c2", "c2", 600, 800);

  gnominalsf->GetYaxis()->SetRangeUser(0.5,1.5);

  gnominalsf->GetYaxis()->SetTitle("Data/MC Scale Factor SF_{b}");
  gnominalsf->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  gnominalsf->Draw("AP");

  if (glowsf) glowsf->Draw("P");
  if (gmedsf) gmedsf->Draw("P");
  if (ghighsf) ghighsf->Draw("P");

  c2->SetGrid();
  label->Draw();
  TLegend* leg2 = new TLegend(0.75,.75,0.9,.91); 
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);

  leg2->AddEntry(gnominal,labels[0], "P");
  leg2->AddEntry(glow,labels[1], "P");
  if (nfiles>=2)
    leg2->AddEntry(gmed,labels[2], "P");
  if (nfiles>=3)
    leg2->AddEntry(ghigh,labels[3], "P");

  leg2->Draw();

  pdfname = "superimpose_sfb_"+suffix+"_"+TString(OP)+".pdf";
  c2->Print(pdfname,"pdf");
  
  tfilenominal->Close();
  tfilelow->Close();
  tfilemed->Close();
  tfilehigh->Close();
 
}

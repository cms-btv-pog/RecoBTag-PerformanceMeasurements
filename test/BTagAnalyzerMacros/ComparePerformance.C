#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "THStack.h"
#include "TFile.h"
#include "TROOT.h"
#include "TColor.h"
#include "TF1.h"
#include "TMath.h"
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TString.h>

using namespace std;

TString adddir = "root_genPT30to80/";
TString filename_1=adddir+"QCD_8P8NoPU.root";
TString labelfile1="QCD 8TeV Pyth8 No PU";
TString filename_2=adddir+"QCD_8P8wPU.root";
TString labelfile2="QCD 8TeV Pyth8 w/ PU";

bool boolfile3=false;
TString filename_3=adddir+"QCD_8P6NoPU.root";
TString labelfile3="QCD 8TeV Pyth6 No PU";

bool boolfile4=true;
TString filename_4=adddir+"QCD_8P6wPU.root";
TString labelfile4="QCD 8TeV Pyth6 w/ PU";

bool boolfile5=false;
TString filename_5=adddir+"QCD_13P8NoPU.root";
TString labelfile5="QCD 13TeV Pyth8 No PU";

bool boolfile6=true;
TString filename_6=adddir+"QCD_13P8wPU.root";
TString labelfile6="QCD 13TeV Pyth8 w/ PU";

bool boolfile7=true;
TString filename_7=adddir+"phys14.root";
TString labelfile7="QCD 13TeV Pyth8 72X";
bool boolfile8=true;
TString filename_8=adddir+"8tev_P8.root";
TString labelfile8="QCD 8TeV Pyth8 53X";
bool boolfile9=true;
TString filename_9=adddir+"8tev_P6.root";
TString labelfile9="QCD 8TeV Pyth6 53X";

bool merge_g=true;

TString dir4plots="perf_plots";

//TString title = "CMS Simulation, HLT_BTagMu_DiJet40_Mu5, Jet pT > 60 GeV/c";
TString title = "test QCD50-80, Jet gen pT : 50-200 GeV/c";
TString titlePhys = "Performance for PHYS14 QCD samples";
TString titlePhys1 = "Performance for PHYS14 QCD samples (Jet gen pT 50-200 GeV/c)";
TString titlePhys2 = "Performance for PHYS14 QCD samples (Jet gen pT 500-1000 GeV/c)";
//TString title = "test QCD50-80, Jet reco pT : 30-120 GeV/c";


TString format=".gif"; // .png or .pdf or .gif

void DrawPlot(TString name, bool log, TString titlename=titlePhys);
void DrawEff(TString name, bool log, TString titlename=titlePhys);
void DrawEffvsPt(TString name, bool log, TString titlename=titlePhys, int icase1=1);
void DrawAllPlot(bool log, int case1=0);
void DrawTest(bool log, int case1=0);
void ComputeEffAndErr(TH1D* disc_b, TH1D* disc_c, TH1D* disc_l, Double_t * eff_b, Double_t * err_b, Double_t * eff_c, Double_t * err_c, Double_t * eff_l, Double_t * err_l);
void DefinePerfTGraph( TString  filename_x, TString name, TGraphErrors* & BvsUDSG_1, TString ytitle_for_BvsUDSG_1, TGraphErrors* &  BvsC_1, TString ytitle_for_BvsC_1, int iMstyle, float xMsize, int iColor);
void DefineEffHist(TString filename_x, TString name, TH1D * & Eff1_b, TH1D * & Eff1_c, TH1D * & Eff1_l, int iColor);

//--------

void ComparePerformance(){

    TString action = "mkdir "+dir4plots+"/";
    system(action);

    DrawPlot("TCHE"	      ,1);
    DrawPlot("TCHP"	      ,1);
    DrawPlot("JP"	      ,1);
    DrawPlot("JBP"	      ,1);
    DrawPlot("SSV"	      ,1);
    DrawPlot("SSVHP"      ,1);
    DrawPlot("CSV"	      ,1);
    DrawPlot("CSVIVF"     ,1);

}

void DrawPlot(TString name, bool log, TString titlename){

    cout << "Name of histo:" << name << '\n';

    TString mistag;
    if (!merge_g) mistag = "uds-jet misid. probability";
    else mistag = "udsg-jet misid. probability";

    TGraphErrors* BvsUDSG_1 = 0;
    TGraphErrors* BvsC_1 = 0;
    DefinePerfTGraph(filename_1, name, BvsUDSG_1, mistag, BvsC_1, "c-jet misid. probability", 20, 0.75, kBlue); 

    TGraphErrors* BvsUDSG_2 =0;
    TGraphErrors* BvsC_2 =0;
    DefinePerfTGraph(filename_2, name, BvsUDSG_2, mistag, BvsC_2, "c-jet misid. probability", 24, 0.75, kOrange); 

    TGraphErrors* BvsUDSG_3 =0;
    TGraphErrors* BvsC_3 =0;
    if (boolfile3) DefinePerfTGraph(filename_3, name, BvsUDSG_3, mistag, BvsC_3, "c-jet misid. probability", 21, 0.75, kGray+2); 

    TGraphErrors* BvsUDSG_4 =0;
    TGraphErrors* BvsC_4 =0;
    if (boolfile4) DefinePerfTGraph(filename_4, name, BvsUDSG_4, mistag, BvsC_4, "c-jet misid. probability", 25, 0.75, kViolet); 

    TGraphErrors* BvsUDSG_5 =0;
    TGraphErrors* BvsC_5 =0;
    if (boolfile5) DefinePerfTGraph(filename_5, name, BvsUDSG_5, mistag, BvsC_5, "c-jet misid. probability", 22, 0.75, kCyan-6); 

    TGraphErrors* BvsUDSG_6 =0;
    TGraphErrors* BvsC_6 =0;
    if (boolfile6) DefinePerfTGraph(filename_6, name, BvsUDSG_6, mistag, BvsC_6, "c-jet misid. probability", 26, 0.75, 46); 

    TGraphErrors* BvsUDSG_7 =0;
    TGraphErrors* BvsC_7 =0;
    if (boolfile7) DefinePerfTGraph(filename_7, name, BvsUDSG_7, mistag, BvsC_7, "c-jet misid. probability", 33, 0.75, 2); 
    TGraphErrors* BvsUDSG_8 =0;
    TGraphErrors* BvsC_8 =0;
    if (boolfile8) DefinePerfTGraph(filename_8, name, BvsUDSG_8, mistag, BvsC_8, "c-jet misid. probability", 33, 0.75, 9); 
    TGraphErrors* BvsUDSG_9 =0;
    TGraphErrors* BvsC_9 =0;
    if (boolfile9) DefinePerfTGraph(filename_9, name, BvsUDSG_9, mistag, BvsC_9, "c-jet misid. probability", 33, 0.75, kGreen+2); 

    if (boolfile7) {
        BvsUDSG_7->SetLineWidth(2);
        BvsC_7->SetLineWidth(2);
    }
    if (boolfile8) {
        BvsUDSG_8->SetLineWidth(2);
        BvsC_8->SetLineWidth(2);
    }
    if (boolfile9) {
        BvsUDSG_9->SetLineWidth(2);
        BvsC_9->SetLineWidth(2);
    }

    // CREATE CANVAS
    gStyle->SetOptTitle(0);
  
    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->cd();

    TPad *c1_1 = new TPad("canvas_1", "canvas_1",0,0,1.0,1.0);
    c1_1->Draw();
    c1_1->cd(); 
    c1_1->SetLogy(log);

    TLegend* qw = 0;
    qw = new TLegend(0.15,0.7,0.45,0.87);
    qw->SetFillColor(kWhite);
    qw->AddEntry(BvsUDSG_1, labelfile1, "p");
    qw->AddEntry(BvsUDSG_2, labelfile2, "p");
    if (boolfile3) qw->AddEntry(BvsUDSG_3, labelfile3, "p");
    if (boolfile4) qw->AddEntry(BvsUDSG_4, labelfile4, "p");
    if (boolfile5) qw->AddEntry(BvsUDSG_5, labelfile5, "p");
    if (boolfile6) qw->AddEntry(BvsUDSG_6, labelfile6, "p");
    if (boolfile7) qw->AddEntry(BvsUDSG_7, labelfile7, "p");
    if (boolfile8) qw->AddEntry(BvsUDSG_8, labelfile8, "p");
    if (boolfile9) qw->AddEntry(BvsUDSG_9, labelfile9, "p");

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.045);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);
    

    BvsUDSG_1->Draw("ALP");
    BvsUDSG_2->Draw("LPsame");
    if (boolfile3) BvsUDSG_3->Draw("LPsame");
    if (boolfile4) BvsUDSG_4->Draw("LPsame");
    if (boolfile5) BvsUDSG_5->Draw("LPsame");
    if (boolfile6) BvsUDSG_6->Draw("LPsame");
    if (boolfile7) BvsUDSG_7->Draw("LPsame");
    if (boolfile8) BvsUDSG_8->Draw("LPsame");
    if (boolfile9) BvsUDSG_9->Draw("LPsame");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw->Draw();

    TString nameuuu="bvsuds_";
    if (merge_g) nameuuu="bvsudsg_";
    TString name_plot=nameuuu+name+"_Linear"+format;
    if(log) name_plot=nameuuu+name+"_Log"+format;
    c1->SaveAs(dir4plots+"/"+name_plot);


    TCanvas *c2 = new TCanvas("c2", "c2",10,32,782,552);
    c2->cd();
    TPad *c2_1 = new TPad("canvas_2", "canvas_2",0,0,1.0,1.0);
    c2_1->Draw();
    c2_1->cd();
    c2_1->SetLogy(log);

    BvsC_1->Draw("ALP");
    BvsC_2->Draw("LPsame");
    if (boolfile3) BvsC_3->Draw("LPsame");
    if (boolfile4) BvsC_4->Draw("LPsame");
    if (boolfile5) BvsC_5->Draw("LPsame");
    if (boolfile6) BvsC_6->Draw("LPsame");
    if (boolfile7) BvsC_7->Draw("LPsame");
    if (boolfile8) BvsC_8->Draw("LPsame");
    if (boolfile9) BvsC_9->Draw("LPsame");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw->Draw();

    name_plot="bvsc_"+name+"_Linear"+format;
    if(log) name_plot="bvsc_"+name+"_Log"+format;
    c2->SaveAs(dir4plots+"/"+name_plot);

  
}

void DrawEff(TString name, bool log, TString titlename){

    cout << "Name of histo:" << name << '\n';

    TH1D * Eff1_b = 0;
    TH1D * Eff1_c = 0;
    TH1D * Eff1_l = 0;
    DefineEffHist(filename_1, name, Eff1_b, Eff1_c, Eff1_l, kBlue);

    TH1D * Eff2_b = 0;
    TH1D * Eff2_c = 0;
    TH1D * Eff2_l = 0;
    DefineEffHist(filename_2, name, Eff2_b, Eff2_c, Eff2_l, kRed);

    TH1D * Eff3_b = 0;
    TH1D * Eff3_c = 0;
    TH1D * Eff3_l = 0;
    if (boolfile3) DefineEffHist(filename_3, name, Eff3_b, Eff3_c, Eff3_l, kGreen+2);

    TH1D * Eff4_b = 0;
    TH1D * Eff4_c = 0;
    TH1D * Eff4_l = 0;
    if (boolfile4) DefineEffHist(filename_4, name, Eff4_b, Eff4_c, Eff4_l, kOrange);


    TH1D * Eff5_b = 0;
    TH1D * Eff5_c = 0;
    TH1D * Eff5_l = 0;
    if (boolfile5) DefineEffHist(filename_5, name, Eff5_b, Eff5_c, Eff5_l, 13);

    TH1D * Eff6_b = 0;
    TH1D * Eff6_c = 0;
    TH1D * Eff6_l = 0;
    if (boolfile6) DefineEffHist(filename_6, name, Eff6_b, Eff6_c, Eff6_l, 28);

    TH1D * Eff7_b = 0;
    TH1D * Eff7_c = 0;
    TH1D * Eff7_l = 0;
    if (boolfile7) DefineEffHist(filename_7, name, Eff7_b, Eff7_c, Eff7_l, 2);
    TH1D * Eff8_b = 0;
    TH1D * Eff8_c = 0;
    TH1D * Eff8_l = 0;
    if (boolfile8) DefineEffHist(filename_8, name, Eff8_b, Eff8_c, Eff8_l, 9);
    TH1D * Eff9_b = 0;
    TH1D * Eff9_c = 0;
    TH1D * Eff9_l = 0;
    if (boolfile9) DefineEffHist(filename_9, name, Eff9_b, Eff9_c, Eff9_l, kGreen+2);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // CREATE CANVAS
  
    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->cd();

    TPad *c1_1 = new TPad("canvas_1", "canvas_1",0,0,1.0,1.0);
    c1_1->Draw();
    c1_1->cd(); 
    c1_1->SetLogy(log);

    TLegend* qw = 0;
    qw = new TLegend(0.15,0.15,0.45,0.32);
    qw->SetFillColor(kWhite);
    qw->AddEntry(Eff1_b, labelfile1, "l");
    qw->AddEntry(Eff2_b, labelfile2, "l");
    if (boolfile3) qw->AddEntry(Eff3_b, labelfile3, "l");
    if (boolfile4) qw->AddEntry(Eff4_b, labelfile4, "l");
    if (boolfile5) qw->AddEntry(Eff5_b, labelfile5, "l");
    if (boolfile6) qw->AddEntry(Eff6_b, labelfile6, "l");
    if (boolfile7) qw->AddEntry(Eff7_b, labelfile7, "l");
    if (boolfile8) qw->AddEntry(Eff8_b, labelfile8, "l");
    if (boolfile9) qw->AddEntry(Eff9_b, labelfile9, "l");

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.045);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);
    
    Eff1_b->Draw();
    Eff2_b->Draw("same");
    if (boolfile3) Eff3_b->Draw("same");
    if (boolfile4) Eff4_b->Draw("same");
    if (boolfile5) Eff5_b->Draw("same");
    if (boolfile6) Eff6_b->Draw("same");
    if (boolfile7) Eff7_b->Draw("same");
    if (boolfile8) Eff8_b->Draw("same");
    if (boolfile9) Eff9_b->Draw("same");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw->Draw();

    TString name_plot="effb_"+name+"_Linear"+format;
    if(log) name_plot="effb_"+name+"_Log"+format;
    c1->SaveAs(dir4plots+"/"+name_plot);

    TCanvas *c2 = new TCanvas("c2", "c2",10,32,782,552);
    c2->cd();
    TPad *c2_1 = new TPad("canvas_2", "canvas_2",0,0,1.0,1.0);
    c2_1->Draw();
    c2_1->cd();
    c2_1->SetLogy(log);

    Eff1_c->Draw();
    Eff2_c->Draw("same");
    if (boolfile3) Eff3_c->Draw("same");
    if (boolfile4) Eff4_c->Draw("same");
    if (boolfile5) Eff5_c->Draw("same");
    if (boolfile6) Eff6_c->Draw("same");
    if (boolfile7) Eff7_c->Draw("same");
    if (boolfile8) Eff8_c->Draw("same");
    if (boolfile9) Eff9_c->Draw("same");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw->Draw();
    name_plot="effc_"+name+"_Linear"+format;
    if(log) name_plot="effc_"+name+"_Log"+format;
    c2->SaveAs(dir4plots+"/"+name_plot);
  
    TCanvas *c3 = new TCanvas("c3", "c3",10,32,782,552);
    c3->cd();
    TPad *c3_1 = new TPad("canvas_3", "canvas_3",0,0,1.0,1.0);
    c3_1->Draw();
    c3_1->cd();
    c3_1->SetLogy(log);

    Eff1_l->Draw();
    Eff2_l->Draw("same");
    if (boolfile3) Eff3_l->Draw("same");
    if (boolfile4) Eff4_l->Draw("same");
    if (boolfile5) Eff5_l->Draw("same");
    if (boolfile6) Eff6_l->Draw("same");
    if (boolfile7) Eff7_l->Draw("same");
    if (boolfile8) Eff8_l->Draw("same");
    if (boolfile9) Eff9_l->Draw("same");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw->Draw();
    TString nameuuu="effl_";
    if (merge_g) nameuuu="effl_wg_";
    name_plot=nameuuu+name+"_Linear"+format;
    if(log) name_plot=nameuuu+name+"_Log"+format;
    c3->SaveAs(dir4plots+"/"+name_plot);
}

void DrawEffvsPt(TString name, bool log, TString titlename, int icase1){

    cout << "Name of histo: jet" << icase1<< "_pt_" << name << '\n';
    TString scase1="";
    if (icase1==1) scase1="1";
    else if (icase1==2) scase1="2";
    cout << "jet"+ scase1 +"_ptgen_b" << " and " << "jet"+ scase1 +"_pt_" +name+"_b" << endl;

    // File1
    TH1D* hist_b_r1;
    TH1D* hist_c_r1;
    TH1D* hist_l_r1; 
    TH1D* hist_b_1;
    TH1D* hist_c_1;
    TH1D* hist_l_1; 

    TFile *myFile_1     = new TFile(filename_1);
    myFile_1->cd();
    
    hist_b_r1         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_b");
    hist_c_r1         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_c");
    hist_l_r1         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_l");
    hist_b_1         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_b");
    hist_c_1         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_c");
    hist_l_1         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_l");


    int nbin_max_1= hist_b_1->GetNbinsX();
    float minx=hist_b_1->GetXaxis()->GetXmin();
    float maxx=hist_b_1->GetXaxis()->GetXmax();

    TH1D * Eff1_b= new TH1D ("Eff1_b",name,nbin_max_1,minx,maxx);
    TH1D * Eff1_c = new TH1D ("Eff1_c",name,nbin_max_1,minx,maxx);
    TH1D * Eff1_l = new TH1D ("Eff1_l",name,nbin_max_1,minx,maxx);

    // equivalent numbers (= number without weights with the same statistical error than weighted events)

    for (int ii=1; ii<nbin_max_1+1; ii++) {
          double val_b=0., val_c=0., val_l=0.;
          double error_b=0., error_c=0., error_l=0.;
          double neq_b=0;
          double neq_c=0;
          double neq_l=0;

          if (hist_b_r1->GetBinError(ii)>0)  neq_b= (hist_b_r1->GetBinContent(ii)*hist_b_r1->GetBinContent(ii)) / (hist_b_r1->GetBinError(ii)*hist_b_r1->GetBinError(ii));
          if (hist_c_r1->GetBinError(ii)>0)  neq_c= (hist_c_r1->GetBinContent(ii)*hist_c_r1->GetBinContent(ii)) / (hist_c_r1->GetBinError(ii)*hist_c_r1->GetBinError(ii));
          if (hist_l_r1->GetBinError(ii)>0)  neq_l= (hist_l_r1->GetBinContent(ii)*hist_l_r1->GetBinContent(ii)) / (hist_l_r1->GetBinError(ii)*hist_l_r1->GetBinError(ii));

          if (hist_b_r1->GetBinContent(ii)>0) val_b = (hist_b_1->GetBinContent(ii))/hist_b_r1->GetBinContent(ii);
          if (neq_b > 0.) error_b = sqrt(neq_b *val_b*(1.-val_b)) / neq_b;
          if (hist_c_r1->GetBinContent(ii)>0) val_c = (hist_c_1->GetBinContent(ii))/hist_c_r1->GetBinContent(ii);
          if (neq_c > 0.) error_c = sqrt(neq_c *val_c*(1.-val_c)) / neq_c;
          if (hist_l_r1->GetBinContent(ii)>0) val_l = (hist_l_1->GetBinContent(ii))/hist_l_r1->GetBinContent(ii);
          if (neq_l > 0.) error_l = sqrt(neq_l *val_l*(1.-val_l)) / neq_l;

          Eff1_b->SetBinContent(ii, val_b);
          Eff1_b->SetBinError(ii,   error_b);
          Eff1_c->SetBinContent(ii, val_c);
          Eff1_c->SetBinError(ii,   error_c);
          Eff1_l->SetBinContent(ii, val_l);
          Eff1_l->SetBinError(ii,   error_l);
    }

    Eff1_b->SetLineColor(kBlue);
    Eff1_c->SetLineColor(kBlue);
    Eff1_l->SetLineColor(kBlue);
    Eff1_b->SetMarkerColor(kBlue);
    Eff1_c->SetMarkerColor(kBlue);
    Eff1_l->SetMarkerColor(kBlue);

    Eff1_b->GetYaxis()->SetTitle("b-jet efficiency");  
    Eff1_l->GetYaxis()->SetTitle("udsg-jet misid. probability");
    Eff1_c->GetYaxis()->SetTitle("c-jet misid. probability");

    Eff1_b->SetMaximum(1.);        
    Eff1_c->SetMaximum(1.);        
    Eff1_l->SetMaximum(1.);        
    
    // File2
    TH1D* hist_b_r2;
    TH1D* hist_c_r2;
    TH1D* hist_l_r2; 
    TH1D* hist_b_2;
    TH1D* hist_c_2;
    TH1D* hist_l_2; 

    TFile *myFile_2     = new TFile(filename_2);
    myFile_2->cd();
    
    hist_b_r2         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_b");
    hist_c_r2         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_c");
    hist_l_r2         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_l");
    hist_b_2         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_b");
    hist_c_2         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_c");
    hist_l_2         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_l");


    int nbin_max_2= hist_b_2->GetNbinsX();

    TH1D * Eff2_b= new TH1D ("Eff2_b",name,nbin_max_2,minx,maxx);
    TH1D * Eff2_c = new TH1D ("Eff2_c",name,nbin_max_2,minx,maxx);
    TH1D * Eff2_l = new TH1D ("Eff2_l",name,nbin_max_2,minx,maxx);

    // equivalent numbers (= number without weights with the same statistical error than weighted events)

    for (int ii=1; ii<nbin_max_2+1; ii++) {
          double val_b=0., val_c=0., val_l=0.;
          double error_b=0., error_c=0., error_l=0.;
          double neq_b=0;
          double neq_c=0;
          double neq_l=0;

          if (hist_b_r2->GetBinError(ii)>0)  neq_b= (hist_b_r2->GetBinContent(ii)*hist_b_r2->GetBinContent(ii)) / (hist_b_r2->GetBinError(ii)*hist_b_r2->GetBinError(ii));
          if (hist_c_r2->GetBinError(ii)>0)  neq_c= (hist_c_r2->GetBinContent(ii)*hist_c_r2->GetBinContent(ii)) / (hist_c_r2->GetBinError(ii)*hist_c_r2->GetBinError(ii));
          if (hist_l_r2->GetBinError(ii)>0)  neq_l= (hist_l_r2->GetBinContent(ii)*hist_l_r2->GetBinContent(ii)) / (hist_l_r2->GetBinError(ii)*hist_l_r2->GetBinError(ii));

          if (hist_b_r2->GetBinContent(ii)>0) val_b = (hist_b_2->GetBinContent(ii))/hist_b_r2->GetBinContent(ii);
          if (neq_b > 0.) error_b = sqrt(neq_b *val_b*(1.-val_b)) / neq_b;
          if (hist_c_r2->GetBinContent(ii)>0) val_c = (hist_c_2->GetBinContent(ii))/hist_c_r2->GetBinContent(ii);
          if (neq_c > 0.) error_c = sqrt(neq_c *val_c*(1.-val_c)) / neq_c;
          if (hist_l_r2->GetBinContent(ii)>0) val_l = (hist_l_2->GetBinContent(ii))/hist_l_r2->GetBinContent(ii);
          if (neq_l > 0.) error_l = sqrt(neq_l *val_l*(1.-val_l)) / neq_l;

          Eff2_b->SetBinContent(ii, val_b);
          Eff2_b->SetBinError(ii,   error_b);
          Eff2_c->SetBinContent(ii, val_c);
          Eff2_c->SetBinError(ii,   error_c);
          Eff2_l->SetBinContent(ii, val_l);
          Eff2_l->SetBinError(ii,   error_l);
    }

    Eff2_b->SetLineColor(kRed);
    Eff2_c->SetLineColor(kRed);
    Eff2_l->SetLineColor(kRed);
    Eff2_b->SetMarkerColor(kRed);
    Eff2_c->SetMarkerColor(kRed);
    Eff2_l->SetMarkerColor(kRed);

    Eff2_b->GetYaxis()->SetTitle("b-jet efficiency");  
    Eff2_l->GetYaxis()->SetTitle("udsg-jet misid. probability");
    Eff2_c->GetYaxis()->SetTitle("c-jet misid. probability");

    Eff2_b->SetMaximum(1.);        
    Eff2_c->SetMaximum(1.);        
    Eff2_l->SetMaximum(1.);        

    // File3
    TH1D* hist_b_r3;
    TH1D* hist_c_r3;
    TH1D* hist_l_r3; 
    TH1D* hist_b_3;
    TH1D* hist_c_3;
    TH1D* hist_l_3; 

    if (boolfile3) {
     TFile *myFile_3     = new TFile(filename_3);
     myFile_3->cd();
     
     hist_b_r3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_b");
     hist_c_r3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_c");
     hist_l_r3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_l");
     hist_b_3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_b");
     hist_c_3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_c");
     hist_l_3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_l");
    }
    else {
     myFile_2->cd();
    
     hist_b_r3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_b");
     hist_c_r3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_c");
     hist_l_r3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_ptgen_l");
     hist_b_3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_b");
     hist_c_3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_c");
     hist_l_3         = (TH1D*)gROOT->FindObject("jet"+ scase1 +"_pt_" +name+"_l");
    }



    int nbin_max_3= hist_b_3->GetNbinsX();

    TH1D * Eff3_b= new TH1D ("Eff3_b",name,nbin_max_3,minx,maxx);
    TH1D * Eff3_c = new TH1D ("Eff3_c",name,nbin_max_3,minx,maxx);
    TH1D * Eff3_l = new TH1D ("Eff3_l",name,nbin_max_3,minx,maxx);

    // equivalent numbers (= number without weights with the same statistical error than weighted events)

    for (int ii=1; ii<nbin_max_3+1; ii++) {
          double val_b=0., val_c=0., val_l=0.;
          double error_b=0., error_c=0., error_l=0.;
          double neq_b=0;
          double neq_c=0;
          double neq_l=0;

          if (hist_b_r3->GetBinError(ii)>0)  neq_b= (hist_b_r3->GetBinContent(ii)*hist_b_r3->GetBinContent(ii)) / (hist_b_r3->GetBinError(ii)*hist_b_r3->GetBinError(ii));
          if (hist_c_r3->GetBinError(ii)>0)  neq_c= (hist_c_r3->GetBinContent(ii)*hist_c_r3->GetBinContent(ii)) / (hist_c_r3->GetBinError(ii)*hist_c_r3->GetBinError(ii));
          if (hist_l_r3->GetBinError(ii)>0)  neq_l= (hist_l_r3->GetBinContent(ii)*hist_l_r3->GetBinContent(ii)) / (hist_l_r3->GetBinError(ii)*hist_l_r3->GetBinError(ii));

          if (hist_b_r3->GetBinContent(ii)>0) val_b = (hist_b_3->GetBinContent(ii))/hist_b_r3->GetBinContent(ii);
          if (neq_b > 0.) error_b = sqrt(neq_b *val_b*(1.-val_b)) / neq_b;
          if (hist_c_r3->GetBinContent(ii)>0) val_c = (hist_c_3->GetBinContent(ii))/hist_c_r3->GetBinContent(ii);
          if (neq_c > 0.) error_c = sqrt(neq_c *val_c*(1.-val_c)) / neq_c;
          if (hist_l_r3->GetBinContent(ii)>0) val_l = (hist_l_3->GetBinContent(ii))/hist_l_r3->GetBinContent(ii);
          if (neq_l > 0.) error_l = sqrt(neq_l *val_l*(1.-val_l)) / neq_l;

          Eff3_b->SetBinContent(ii, val_b);
          Eff3_b->SetBinError(ii,   error_b);
          Eff3_c->SetBinContent(ii, val_c);
          Eff3_c->SetBinError(ii,   error_c);
          Eff3_l->SetBinContent(ii, val_l);
          Eff3_l->SetBinError(ii,   error_l);
    }


    Eff3_b->SetLineColor(kGreen+2);
    Eff3_c->SetLineColor(kGreen+2);
    Eff3_l->SetLineColor(kGreen+2);
    Eff3_b->SetMarkerColor(kGreen+2);
    Eff3_c->SetMarkerColor(kGreen+2);
    Eff3_l->SetMarkerColor(kGreen+2);

    Eff3_b->GetYaxis()->SetTitle("b-jet efficiency");  
    Eff3_l->GetYaxis()->SetTitle("udsg-jet misid. probability");
    Eff3_c->GetYaxis()->SetTitle("c-jet misid. probability");

    Eff3_b->SetMaximum(1.);        
    Eff3_c->SetMaximum(1.);        
    Eff3_l->SetMaximum(1.);        

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // CREATE CANVAS
  
    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->cd();

    // FIRST MC & DATA
    TPad *c1_1 = new TPad("canvas_1", "canvas_1",0,0,1.0,1.0);
    c1_1->Draw();
    c1_1->cd(); 
    c1_1->SetLogy(log);

    TLegend* qw = 0;
    qw = new TLegend(0.15,0.15,0.45,0.32);
    qw->SetFillColor(kWhite);
    qw->AddEntry(Eff1_b, "13TeV 20PU25ns Phys14", "l");
    qw->AddEntry(Eff2_b, "13TeV 20PU25ns CSA14", "l");
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")) qw->AddEntry(Eff3_b, "8TeV nPUtrue17-23 53X", "l");
    TLegend* qw2 = 0;
    qw2 = new TLegend(0.65,0.68,0.95,0.95);
    qw2->SetFillColor(kWhite);
    qw2->AddEntry(Eff1_b, "13TeV 20PU25ns Phys14", "l");
    qw2->AddEntry(Eff2_b, "13TeV 20PU25ns CSA14", "l");
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")) qw2->AddEntry(Eff3_b, "8TeV nPUtrue17-23 53X", "l");

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.045);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);
    
    Eff1_b->Draw();
    Eff2_b->Draw("same");
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")) Eff3_b->Draw("same");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw->Draw();
    TString name_plot="effvspt_b_"+name+"_Linear_"+scase1+format;
    if(log) name_plot="effvspt_b_"+name+"_Log_"+scase1+format;
    c1->SaveAs(dir4plots+"/"+name_plot);

    TCanvas *c2 = new TCanvas("c2", "c2",10,32,782,552);
    c2->cd();
    TPad *c2_1 = new TPad("canvas_2", "canvas_2",0,0,1.0,1.0);
    c2_1->Draw();
    c2_1->cd();
    c2_1->SetLogy(log);

    if (Eff1_c->GetMinimum()> Eff2_c->GetMinimum()) Eff1_c->SetMinimum(Eff2_c->GetMinimum()*0.5) ;
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")  and Eff1_c->GetMinimum()> Eff3_c->GetMinimum()) Eff1_c->SetMinimum(Eff3_c->GetMinimum()*0.5) ;
    Eff1_c->Draw();
    Eff2_c->Draw("same");
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")) Eff3_c->Draw("same");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw2->Draw();
    name_plot="effvspt_c_"+name+"_Linear_"+scase1+format;
    if(log) name_plot="effvspt_c_"+name+"_Log_"+scase1+format;
    c2->SaveAs(dir4plots+"/"+name_plot);
  
    TCanvas *c3 = new TCanvas("c3", "c3",10,32,782,552);
    c3->cd();
    TPad *c3_1 = new TPad("canvas_3", "canvas_3",0,0,1.0,1.0);
    c3_1->Draw();
    c3_1->cd();
    c3_1->SetLogy(log);

    if (Eff1_l->GetMinimum()> Eff2_l->GetMinimum()) Eff1_l->SetMinimum(Eff2_l->GetMinimum()*0.5) ;
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")  and Eff1_l->GetMinimum()> Eff3_l->GetMinimum()) Eff1_l->SetMinimum(Eff3_l->GetMinimum()*0.5) ;
    Eff1_l->Draw();
    Eff2_l->Draw("same");
    if (boolfile3 and (name != "CSVIVF" && name != "CSVIVF_1" && name != "CSVIVF_2" && name != "csvivfl" && name != "csvivfm")) Eff3_l->Draw("same");
    latex->DrawLatex(0.1, 0.94, titlename);
    qw2->Draw();
    name_plot="effvspt_l_"+name+"_Linear_"+scase1+format;
    if(log) name_plot="effvspt_l_"+name+"_Log_"+scase1+format;
    c3->SaveAs(dir4plots+"/"+name_plot);
}



void DrawAllPlot(bool log, int case1){


    TFile *myFile     = new TFile(filename_1);
    myFile->cd();

    TGraph* EffBvsUDSG[100];
    TGraph* EffBvsC[100];

    // CREATE CANVAS

    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->cd();

    // FIRST MC & DATA
    TPad *c1_1 = new TPad("canvas_1", "canvas_1",0,0,1.0,1.0);
    c1_1->Draw();
    c1_1->cd();
    c1_1->SetLogy(log);

    TCanvas *c2 = new TCanvas("c2", "c2",10,32,782,552);
    c2->cd();

    // FIRST MC & DATA
    TPad *c2_1 = new TPad("canvas_2", "canvas_2",0,0,1.0,1.0);
    c2_1->Draw();
    c2_1->cd();
    c2_1->SetLogy(log);

    gStyle->SetOptTitle(0);

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.045);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);

    TLegend* qw = 0;
    qw = new TLegend(0.15,0.6,0.25,0.87);
    qw->SetFillColor(kWhite);

    const char *disc_array[] = {"TCHE","TCHP","JP","JBP","SSV","SSVHP","CSV","CSVIVF"};
    const char *disc_array_1[] = {"TCHE_1","TCHP_1","JP_1","JBP_1","SSV_1","SSVHP_1","CSV_1","CSVIVF_1"};
    const char *disc_array_2[] = {"TCHE_2","TCHP_2","JP_2","JBP_2","SSV_2","SSVHP_2","CSV_2","CSVIVF_2"};
    int markerstyle[] = {20,22,33,34,23,29,31,21};
    //int markercolor[] = {46,2,6,4,7,3,1,8};
    int markercolor[] = {kOrange,kOrange+2,kBlue,kBlue-2,8,kGreen+2,1,2};

    std::vector<std::string> disc_names1(disc_array, std::end(disc_array));
    std::vector<std::string> disc_names2(disc_array_1, std::end(disc_array_1));
    std::vector<std::string> disc_names3(disc_array_2, std::end(disc_array_2));
   
    std::vector<std::string> disc_names;
    if (case1==0) disc_names=disc_names1;
    else if (case1==1) disc_names=disc_names2;
    else if (case1==2) disc_names=disc_names3;

    int count = 0;

    cout << "Looping over disc. names\n";

    for(vector<string>::const_iterator iter = disc_names.begin(); iter != disc_names.end(); ++iter) {

//        if (*iter == "CSVIVF" && boolfile3) continue;
//        if (*iter == "CSVIVF_1" && boolfile3) continue;
//        if (*iter == "CSVIVF_2" && boolfile3) continue;

        TH1D* disc_b;
        TH1D* disc_c;
        TH1D* disc_l;

        disc_b         = (TH1D*)gROOT->FindObject((*iter+"_b").c_str());
        disc_c         = (TH1D*)gROOT->FindObject((*iter+"_c").c_str());
        disc_l         = (TH1D*)gROOT->FindObject((*iter+"_l").c_str());

        int nbin_max= disc_b->GetNbinsX();

        Double_t eff_b[nbin_max];
        Double_t eff_l[nbin_max];
        Double_t eff_c[nbin_max];

        Double_t err_b[nbin_max];
        Double_t err_l[nbin_max];
        Double_t err_c[nbin_max];

        ComputeEffAndErr(disc_b, disc_c, disc_l, eff_b, err_b, eff_c, err_c, eff_l, err_l);

        TGraphErrors* BvsUDSG = new TGraphErrors(nbin_max,eff_b,eff_l,err_b,err_l);
        TGraphErrors* BvsC    = new TGraphErrors(nbin_max,eff_b,eff_c,err_b,err_c);

        BvsUDSG->GetXaxis()->SetLimits(0.,1.);
        BvsUDSG->GetHistogram()->SetMaximum(1.);
        BvsUDSG->GetHistogram()->SetMinimum(0.0001);

        BvsC->GetXaxis()->SetLimits(0.,1.);
        BvsC->GetHistogram()->SetMaximum(1.);
        BvsC->GetHistogram()->SetMinimum(0.001);

        // SET COSMETICS

        BvsUDSG->SetMarkerStyle(markerstyle[count]);
        BvsUDSG->SetMarkerColor(markercolor[count]);
        BvsUDSG->SetLineColor(markercolor[count]);
        BvsUDSG->SetMarkerSize(0.75);
        BvsUDSG->GetXaxis()->SetTitle("b-jet efficiency");
        BvsUDSG->GetYaxis()->SetTitle("udsg-jet misid. probability");
        cout << "Cloning " << (*iter).c_str() << '\n';
        EffBvsUDSG[count] = (TGraphErrors*)BvsUDSG->Clone();

        BvsC->SetMarkerStyle(markerstyle[count]);
        BvsC->SetMarkerSize(0.75);
        BvsC->SetMarkerColor(markercolor[count]);
        BvsC->SetLineColor(markercolor[count]);
        BvsC->GetXaxis()->SetTitle("b-jet efficiency");
        BvsC->GetYaxis()->SetTitle("c-jet misid. probability");

        EffBvsC[count] = (TGraphErrors*)BvsC->Clone();

        qw->AddEntry(BvsUDSG, (*iter).c_str() ,"p");
        ++count;
    }

    c1->cd();
    c1_1->cd();
    cout << "Drawing the graphs\n";
    for (int tgs = 0; tgs < count; ++tgs){
        if (tgs == 0) EffBvsUDSG[tgs]->Draw("ALP");
        else EffBvsUDSG[tgs]->Draw("LPsame");
    }
    TString scase1="";
    if (case1==1) scase1="_1";
    else if (case1==2) scase1="_2";
   
    if (case1==0) latex->DrawLatex(0.1, 0.94, titlePhys);
    else if (case1==1) latex->DrawLatex(0.1, 0.94, titlePhys1);
    else if (case1==2) latex->DrawLatex(0.1, 0.94, titlePhys2);
    qw->Draw();
    TString name_plot="bvsudsg_all_Linear"+scase1+format;
    if(log) name_plot="bvsudsg_all_Log"+scase1+format;
    c1->SaveAs(dir4plots+"/"+name_plot);

    c2->cd();
    c2_1->cd();
    for (int tgs = 0; tgs < count; ++tgs){
        if (tgs == 0) EffBvsC[tgs]->Draw("ALP");
        else EffBvsC[tgs]->Draw("LPsame");
    }
    if (case1==0) latex->DrawLatex(0.1, 0.94, titlePhys);
    else if (case1==1) latex->DrawLatex(0.1, 0.94, titlePhys1);
    else if (case1==2) latex->DrawLatex(0.1, 0.94, titlePhys2);
    qw->Draw();
    name_plot="bvsc_all_Linear"+scase1+format;
    if(log) name_plot="bvsc_all_Log"+scase1+format;
    c2->SaveAs(dir4plots+"/"+name_plot);

}


void ComputeEffAndErr(TH1D* disc_b, TH1D* disc_c, TH1D* disc_l, Double_t * eff_b, Double_t * err_b, Double_t * eff_c, Double_t * err_c, Double_t * eff_l, Double_t * err_l)
{

        int nbin_max= disc_b->GetNbinsX();
        double tot_b=disc_b->Integral(0,nbin_max+1);
        double tot_c=disc_c->Integral(0,nbin_max+1);
        double tot_udsg=disc_l->Integral(0,nbin_max+1);

        // equivalent numbers (= number without weights with the same statistical error than weighted events)
        double neq_b=0;
        double neq_c=0;
        double neq_l=0;

        for (int ii=0; ii<nbin_max+2; ii++) {
          if (disc_b->GetBinError(ii)>0)  neq_b+=(disc_b->GetBinContent(ii)*disc_b->GetBinContent(ii)) / (disc_b->GetBinError(ii)*disc_b->GetBinError(ii));
          if (disc_c->GetBinError(ii)>0)  neq_c+=(disc_c->GetBinContent(ii)*disc_c->GetBinContent(ii)) / (disc_c->GetBinError(ii)*disc_c->GetBinError(ii));
          if (disc_l->GetBinError(ii)>0)  neq_l+=(disc_l->GetBinContent(ii)*disc_l->GetBinContent(ii)) / (disc_l->GetBinError(ii)*disc_l->GetBinError(ii));
        }

        for (int ii=0; ii<nbin_max; ii++) {

            double val_b=0., val_c=0., val_udsg=0.;
            double error_b=0., error_c=0., error_udsg=0.;


            if (tot_b > 0.) val_b = disc_b->Integral(ii+1,nbin_max+1) / tot_b;
            if (neq_b > 0.) error_b = sqrt(neq_b *val_b*(1.-val_b)) / neq_b;

            if (tot_c > 0.) val_c = disc_c->Integral(ii+1,nbin_max+1) / tot_c;
            if (neq_c > 0.) error_c = sqrt(neq_c *val_c*(1.-val_c)) / neq_c;

            if (tot_udsg > 0.) val_udsg = disc_l->Integral(ii+1,nbin_max+1) / tot_udsg;
            if (neq_l > 0.) error_udsg = sqrt(neq_l *val_udsg*(1.-val_udsg)) / neq_l;

            eff_b[ii] = val_b;
            err_b[ii] = error_b;

            eff_c[ii] = val_c;
            err_c[ii] = error_c;

            eff_l[ii] = val_udsg;
            err_l[ii] = error_udsg;

        }

}
void DrawTest(bool log, int case1){


    TFile *myFile     = new TFile(filename_1);
    myFile->cd();

    TGraph* EffBvsUDSG[100];
    TGraph* EffBvsC[100];

    // CREATE CANVAS

    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->cd();

    // FIRST MC & DATA
    TPad *c1_1 = new TPad("canvas_1", "canvas_1",0,0,1.0,1.0);
    c1_1->Draw();
    c1_1->cd();
    c1_1->SetLogy(log);

    TCanvas *c2 = new TCanvas("c2", "c2",10,32,782,552);
    c2->cd();

    // FIRST MC & DATA
    TPad *c2_1 = new TPad("canvas_2", "canvas_2",0,0,1.0,1.0);
    c2_1->Draw();
    c2_1->cd();
    c2_1->SetLogy(log);

    gStyle->SetOptTitle(0);

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.045);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);

    TLegend* qw = 0;
    qw = new TLegend(0.15,0.6,0.25,0.87);
    qw->SetFillColor(kWhite);

    const char *disc_array_1[] = {"CSV_1","CSVIVF_1"};
    int markerstyle[] = {31,21};
    int markercolor[] = {1,2};

    std::vector<std::string> disc_names(disc_array_1, std::end(disc_array_1));
   
    int count = 0;
    int count2 = 0;

    cout << "Looping over disc. names\n";

    for(vector<string>::const_iterator iter = disc_names.begin(); iter != disc_names.end(); ++iter) {

//        if (*iter == "CSVIVF" && boolfile3) continue;
//        if (*iter == "CSVIVF_1" && boolfile3) continue;
//        if (*iter == "CSVIVF_2" && boolfile3) continue;

        TH1D* disc_b;
        TH1D* disc_c;
        TH1D* disc_l;

        disc_b         = (TH1D*)gROOT->FindObject((*iter+"_b").c_str());
        disc_c         = (TH1D*)gROOT->FindObject((*iter+"_c").c_str());
        disc_l         = (TH1D*)gROOT->FindObject((*iter+"_l").c_str());

        int nbin_max= disc_b->GetNbinsX();

        Double_t eff_b[nbin_max];
        Double_t eff_l[nbin_max];
        Double_t eff_c[nbin_max];

        Double_t err_b[nbin_max];
        Double_t err_l[nbin_max];
        Double_t err_c[nbin_max];

        ComputeEffAndErr(disc_b, disc_c, disc_l, eff_b, err_b, eff_c, err_c, eff_l, err_l);

        TGraphErrors* BvsUDSG = new TGraphErrors(nbin_max,eff_b,eff_l,err_b,err_l);
        TGraphErrors* BvsC    = new TGraphErrors(nbin_max,eff_b,eff_c,err_b,err_c);

        BvsUDSG->GetXaxis()->SetLimits(0.,1.);
        BvsUDSG->GetHistogram()->SetMaximum(1.);
        BvsUDSG->GetHistogram()->SetMinimum(0.0001);

        BvsC->GetXaxis()->SetLimits(0.,1.);
        BvsC->GetHistogram()->SetMaximum(1.);
        BvsC->GetHistogram()->SetMinimum(0.001);

        // SET COSMETICS

        BvsUDSG->SetMarkerStyle(markerstyle[count]);
        BvsUDSG->SetMarkerColor(markercolor[count]);
        BvsUDSG->SetLineColor(markercolor[count]);
        BvsUDSG->SetMarkerSize(0.75);
        BvsUDSG->GetXaxis()->SetTitle("b-jet efficiency");
        BvsUDSG->GetYaxis()->SetTitle("udsg-jet misid. probability");
        cout << "Cloning " << (*iter).c_str() << '\n';
        EffBvsUDSG[count] = (TGraphErrors*)BvsUDSG->Clone();

        BvsC->SetMarkerStyle(markerstyle[count]);
        BvsC->SetMarkerSize(0.75);
        BvsC->SetMarkerColor(markercolor[count]);
        BvsC->SetLineColor(markercolor[count]);
        BvsC->GetXaxis()->SetTitle("b-jet efficiency");
        BvsC->GetYaxis()->SetTitle("c-jet misid. probability");

        EffBvsC[count] = (TGraphErrors*)BvsC->Clone();

        qw->AddEntry(BvsUDSG, (*iter).c_str() ,"p");
        ++count;
    }
    TFile *myFile2     = new TFile(filename_3);
    myFile2->cd();
    const char *disc_array_2[] = {"CSV_1"};
    int markerstyle2[] = {20};
    int markercolor2[] = {8};

    std::vector<std::string> disc_names2(disc_array_2, std::end(disc_array_2));
   
    cout << "Looping over disc. names\n";

    for(vector<string>::const_iterator iter = disc_names2.begin(); iter != disc_names2.end(); ++iter) {


        TH1D* disc_b;
        TH1D* disc_c;
        TH1D* disc_l;

        disc_b         = (TH1D*)gROOT->FindObject((*iter+"_b").c_str());
        disc_c         = (TH1D*)gROOT->FindObject((*iter+"_c").c_str());
        disc_l         = (TH1D*)gROOT->FindObject((*iter+"_l").c_str());

        int nbin_max= disc_b->GetNbinsX();

        Double_t eff_b[nbin_max];
        Double_t eff_l[nbin_max];
        Double_t eff_c[nbin_max];

        Double_t err_b[nbin_max];
        Double_t err_l[nbin_max];
        Double_t err_c[nbin_max];

        ComputeEffAndErr(disc_b, disc_c, disc_l, eff_b, err_b, eff_c, err_c, eff_l, err_l);

        TGraphErrors* BvsUDSG = new TGraphErrors(nbin_max,eff_b,eff_l,err_b,err_l);
        TGraphErrors* BvsC    = new TGraphErrors(nbin_max,eff_b,eff_c,err_b,err_c);

        BvsUDSG->GetXaxis()->SetLimits(0.,1.);
        BvsUDSG->GetHistogram()->SetMaximum(1.);
        BvsUDSG->GetHistogram()->SetMinimum(0.0001);

        BvsC->GetXaxis()->SetLimits(0.,1.);
        BvsC->GetHistogram()->SetMaximum(1.);
        BvsC->GetHistogram()->SetMinimum(0.001);

        // SET COSMETICS

        BvsUDSG->SetMarkerStyle(markerstyle2[count2]);
        BvsUDSG->SetMarkerColor(markercolor2[count2]);
        BvsUDSG->SetLineColor(markercolor2[count2]);
        BvsUDSG->SetMarkerSize(0.75);
        BvsUDSG->GetXaxis()->SetTitle("b-jet efficiency");
        BvsUDSG->GetYaxis()->SetTitle("udsg-jet misid. probability");
        cout << "Cloning " << (*iter).c_str() << '\n';
        EffBvsUDSG[count] = (TGraphErrors*)BvsUDSG->Clone();

        BvsC->SetMarkerStyle(markerstyle2[count2]);
        BvsC->SetMarkerSize(0.75);
        BvsC->SetMarkerColor(markercolor2[count2]);
        BvsC->SetLineColor(markercolor2[count2]);
        BvsC->GetXaxis()->SetTitle("b-jet efficiency");
        BvsC->GetYaxis()->SetTitle("c-jet misid. probability");

        EffBvsC[count] = (TGraphErrors*)BvsC->Clone();

        TString name = *iter + " 8 TeV";
        qw->AddEntry(BvsUDSG, name ,"p");
        ++count;
        ++count2;
    }

    c1->cd();
    c1_1->cd();
    cout << "Drawing the graphs\n";
    for (int tgs = 0; tgs < count; ++tgs){
        if (tgs == 0) EffBvsUDSG[tgs]->Draw("ALP");
        else EffBvsUDSG[tgs]->Draw("LPsame");
    }
   
    if (case1==0) latex->DrawLatex(0.1, 0.94, titlePhys);
    else if (case1==1) latex->DrawLatex(0.1, 0.94, titlePhys1);
    else if (case1==2) latex->DrawLatex(0.1, 0.94, titlePhys2);
    qw->Draw();
    TString name_plot="bvsudsg_test_Linear"+format;
    if(log) name_plot="bvsudsg_test_Log"+format;
    c1->SaveAs(dir4plots+"/"+name_plot);

    c2->cd();
    c2_1->cd();
    for (int tgs = 0; tgs < count; ++tgs){
        if (tgs == 0) EffBvsC[tgs]->Draw("ALP");
        else EffBvsC[tgs]->Draw("LPsame");
    }
    if (case1==0) latex->DrawLatex(0.1, 0.94, titlePhys);
    else if (case1==1) latex->DrawLatex(0.1, 0.94, titlePhys1);
    else if (case1==2) latex->DrawLatex(0.1, 0.94, titlePhys2);
    qw->Draw();
    name_plot="bvsc_test_Linear"+format;
    if(log) name_plot="bvsc_test_Log"+format;
    c2->SaveAs(dir4plots+"/"+name_plot);

}

void DefinePerfTGraph( TString  filename_x, TString name, TGraphErrors* & BvsUDSG_1, TString ytitle_for_BvsUDSG_1, TGraphErrors* &  BvsC_1, TString ytitle_for_BvsC_1, int iMstyle, float xMsize, int iColor) 
{
    
    // File1
    TH1D* hist_b_1;
    TH1D* hist_c_1;
    TH1D* hist_l_1; 
    TH1D* hist_g_1; 

    TFile *myFile_1     = new TFile(filename_x);
    myFile_1->cd();
    
    hist_b_1         = (TH1D*)gROOT->FindObject(name+"_b");
    hist_c_1         = (TH1D*)gROOT->FindObject(name+"_c");
    hist_l_1         = (TH1D*)gROOT->FindObject(name+"_l");
    if (merge_g) {
     hist_g_1         = (TH1D*)gROOT->FindObject(name+"_g");
     hist_l_1->Add(hist_g_1);
    }


    int nbin_max_1= hist_b_1->GetNbinsX();


    Double_t eff_b_1[nbin_max_1];
    Double_t eff_l_1[nbin_max_1];
    Double_t eff_c_1[nbin_max_1];
    Double_t err_b_1[nbin_max_1];
    Double_t err_l_1[nbin_max_1];
    Double_t err_c_1[nbin_max_1];
    ComputeEffAndErr(hist_b_1, hist_c_1, hist_l_1, eff_b_1, err_b_1, eff_c_1, err_c_1, eff_l_1, err_l_1);

    BvsUDSG_1 = new TGraphErrors(nbin_max_1,eff_b_1,eff_l_1,err_b_1,err_l_1);
    BvsC_1    = new TGraphErrors(nbin_max_1,eff_b_1,eff_c_1,err_b_1,err_c_1);

    // SET COSMETICS
    BvsUDSG_1->SetMarkerStyle(iMstyle);
    BvsUDSG_1->SetMarkerSize(xMsize);
    BvsUDSG_1->SetMarkerColor(iColor);
    BvsUDSG_1->SetLineColor(iColor);
    BvsUDSG_1->GetXaxis()->SetTitle("b-jet efficiency");  
    BvsUDSG_1->GetYaxis()->SetTitle(ytitle_for_BvsUDSG_1);

    BvsC_1->SetMarkerStyle(iMstyle);
    BvsC_1->SetMarkerSize(xMsize);
    BvsC_1->SetMarkerColor(iColor);
    BvsC_1->SetLineColor(iColor);
    BvsC_1->GetXaxis()->SetTitle("b-jet efficiency");  
    BvsC_1->GetYaxis()->SetTitle(ytitle_for_BvsC_1);

    BvsUDSG_1->GetXaxis()->SetLimits(0.,1.);
    BvsUDSG_1->GetHistogram()->SetMaximum(1.);        
    BvsUDSG_1->GetHistogram()->SetMinimum(0.0001);
    
    BvsC_1->GetXaxis()->SetLimits(0.,1.);
    BvsC_1->GetHistogram()->SetMaximum(1.);        
    BvsC_1->GetHistogram()->SetMinimum(0.001);

}

void DefineEffHist(TString filename_x, TString name, TH1D * & Eff1_b, TH1D * & Eff1_c, TH1D * & Eff1_l, int iColor)
{
    // File1
    TH1D* hist_b_1;
    TH1D* hist_c_1;
    TH1D* hist_l_1; 
    TH1D* hist_g_1; 

    TFile *myFile_1     = new TFile(filename_x);
    myFile_1->cd();
    
    hist_b_1         = (TH1D*)gROOT->FindObject(name+"_b");
    hist_c_1         = (TH1D*)gROOT->FindObject(name+"_c");
    hist_l_1         = (TH1D*)gROOT->FindObject(name+"_l");
    if (merge_g) {
     hist_g_1         = (TH1D*)gROOT->FindObject(name+"_g");
     hist_l_1->Add(hist_g_1);
    }


    int nbin_max_1= hist_b_1->GetNbinsX();
    float minx=hist_b_1->GetXaxis()->GetXmin();
    float maxx=hist_b_1->GetXaxis()->GetXmax();

    Eff1_b= new TH1D ("Eff1_b",name,nbin_max_1,minx,maxx);
    Eff1_c = new TH1D ("Eff1_c",name,nbin_max_1,minx,maxx);
    Eff1_l = new TH1D ("Eff1_l",name,nbin_max_1,minx,maxx);

    Double_t eff_b_1[nbin_max_1];
    Double_t eff_l_1[nbin_max_1];
    Double_t eff_c_1[nbin_max_1];
    Double_t err_b_1[nbin_max_1];
    Double_t err_l_1[nbin_max_1];
    Double_t err_c_1[nbin_max_1];
    ComputeEffAndErr(hist_b_1, hist_c_1, hist_l_1, eff_b_1, err_b_1, eff_c_1, err_c_1, eff_l_1, err_l_1);

    for (int ii=0; ii<nbin_max_1; ii++) {
        Eff1_b->SetBinContent(ii+1, eff_b_1[ii]);
        Eff1_b->SetBinError(ii+1,   err_b_1[ii]);
        Eff1_c->SetBinContent(ii+1, eff_c_1[ii]);
        Eff1_c->SetBinError(ii+1,   err_c_1[ii]);
        Eff1_l->SetBinContent(ii+1, eff_l_1[ii]);
        Eff1_l->SetBinError(ii+1,   err_l_1[ii]);
    }

    Eff1_b->SetLineColor(iColor);
    Eff1_c->SetLineColor(iColor);
    Eff1_l->SetLineColor(iColor);
    Eff1_b->SetMarkerColor(iColor);
    Eff1_c->SetMarkerColor(iColor);
    Eff1_l->SetMarkerColor(iColor);

    Eff1_b->GetYaxis()->SetTitle("b-jet efficiency");  
    if (!merge_g) Eff1_l->GetYaxis()->SetTitle("uds-jet misid. probability");
    else Eff1_l->GetYaxis()->SetTitle("udsg-jet misid. probability");
    Eff1_c->GetYaxis()->SetTitle("c-jet misid. probability");

    Eff1_b->SetMaximum(1.);        
    Eff1_c->SetMaximum(1.);        
    Eff1_l->SetMaximum(1.);        
}

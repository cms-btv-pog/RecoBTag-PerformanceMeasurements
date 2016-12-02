#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "THStack.h"
#include "TFile.h"
#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TLegendEntry.h>
#include <TColor.h>

#include "tdrstyle.C"

using namespace std;

TString Jettypetitle="AK4 jets (p_{T} > 20 GeV)";
TString filename="ttbar/output_all.root";
TString output="Commissioning_plots_JetPt20/";
TString CMStitle= "CMS";
TString Preliminarytitle= "Preliminary";
TString Selectiontitle= "e#mu channel, #geq 2 jets";
TString Lumititle= "#sqrt{s} = 13 TeV, 25ns";
//TString title= "CMS 2015 preliminary, #sqrt{s} = 13 TeV,  2.44 fb^{-1}";
TString format1=".pdf"; // .png or .pdf or .gif
TString format2=".png"; // .png or .pdf or .gif
bool bOverflow=true;
bool b_ordering=true;
bool web = true;

// Configuration for CTag commissioning //
bool c_ordering=false;

void Draw(TString name, TString histotitle, bool log, int move_legend=0);
void DrawTTbar(TString name, TString histotitle, bool log, int move_legend=0);
void DrawTagRate(TString name, TString histotitle, bool log);
void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log);
void OverFlowBinFix(TH1D* );

//--------

void DrawCommPlot(bool Draw_track_plots, bool Draw_Nminus1_plots, bool Draw_sv_plots, bool Draw_muons_plots, bool Draw_discriminator_plots , 
bool Draw_newdiscriminator_plots, bool Draw_tagRate_plots, bool Draw_2D_plots)
{

  TString action = "mkdir ttbar/Commissioning_plots";
  system(action);
  
  int move=1;

  DrawTTbar("nPV", "number of PV",0);
  DrawTTbar("met", "MET (GeV)",0);
  DrawTTbar("mll", "M_{ll} (GeV)",0);
  DrawTTbar("njet","number of jets",0);
  DrawTTbar("nbtag","number of btag jets",0);
  DrawTTbar("pt_e","Leading electron P_{T} (GeV)",0);
  DrawTTbar("pt_mu","Leading muon P_{T} (GeV)",0);
  DrawTTbar("pt_jet","Leading jet P_{T} (GeV)",0);

}

//--------------------------

void Draw(TString name, TString histotitle, bool log, int move_legend)
{

 if(c_ordering){
    //Muon channel//
   filename="output_allMuEFG.root";
   output="Commissioning_plots_singleMuEFG_JetPt25_103116/";
   Selectiontitle= "single muon channel, #geq 4 jets";
   //Eletron Channel//
   //filename="output_allElecEFG.root";
   //output="Commissioning_plots_singleElecEFG_JetPt25_103116/";
   //Selectiontitle= "single electron channel, #geq 4 jets";

   TString MakeDir = "mkdir ttbar/"+output;
   TString action = MakeDir;
   system(action);
   Selectiontitle= "single muon channel, #geq 4 jets";
   bOverflow=false;
   b_ordering=false; 
 }

 TH1D* hist_b;
 TH1D* hist_c;
 TH1D* hist_pu;
 //TH1D* hist_gsplit;
 TH1D* hist_l;
 TH1D* hist_data;
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_b         = (TH1D*)gROOT->FindObject(name+"_b");
 hist_c         = (TH1D*)gROOT->FindObject(name+"_c");
 hist_pu        = (TH1D*)gROOT->FindObject(name+"_pu");
 //hist_gsplit    = (TH1D*)gROOT->FindObject(name+"_bfromg");
 hist_l         = (TH1D*)gROOT->FindObject(name+"_l");
 hist_data      = (TH1D*)gROOT->FindObject(name+"_data");
 

 if (bOverflow && name!="SSV" && name!="SSVHP") 
 {
   OverFlowBinFix(hist_b);
   OverFlowBinFix(hist_c);
   OverFlowBinFix(hist_pu);
   //OverFlowBinFix(hist_gsplit);
   OverFlowBinFix(hist_l);
   OverFlowBinFix(hist_data);
 }

 TH1D* histo_tot = (TH1D*) hist_b->Clone();
 histo_tot->Sumw2();
 histo_tot ->Add(hist_c);
 histo_tot ->Add(hist_l);  
 histo_tot ->Add(hist_pu);
 //histo_tot ->Add(hist_gsplit);

 
 //float scale_f = 1.;
 float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_pu->Integral() + hist_l->Integral());
 //float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_gsplit->Integral() + hist_l->Integral());

 hist_b       ->Scale(scale_f);
 hist_c       ->Scale(scale_f);
 hist_pu      ->Scale(scale_f);
 //hist_gsplit  ->Scale(scale_f);
 hist_l       ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);
  
 TH1D* histo_ratio;
 histo_ratio = (TH1D*) hist_data->Clone();
 histo_ratio->SetName("histo_ratio");
 histo_ratio->SetTitle("");
  
 histo_ratio->Divide(histo_tot);
  
 hist_data  ->SetLineWidth(2);
 hist_data  ->SetMarkerStyle(20);  
 hist_data  ->SetMarkerSize(0.75); 

 hist_c     ->SetFillColor(8);
 hist_b     ->SetFillColor(2);
 hist_pu    ->SetFillColor(kYellow);
 //hist_gsplit->SetFillColor(7);
 hist_l     ->SetFillColor(4);
  
 histo_tot  ->SetLineColor(2);
  

 THStack *stack = new THStack("stack","stack");
  
 if (b_ordering)
 {
   stack      ->Add(hist_b);
   stack      ->Add(hist_c);
   stack      ->Add(hist_pu);
   stack      ->Add(hist_l);
   //stack      ->Add(hist_gsplit);
 }
 else 
 {
   stack      ->Add(hist_pu);
   //stack      ->Add(hist_gsplit);
   stack      ->Add(hist_l);
   stack      ->Add(hist_c);
   stack      ->Add(hist_b);
 }

 setTDRStyle();
 gStyle->SetErrorX(0.);

 TCanvas *c1 = new TCanvas(name, name,10,32,782,552);
 c1->SetFillColor(10);
 c1->SetBorderMode(0);
 c1->SetBorderSize(2);
 c1->SetFrameFillColor(0);
 c1->SetFrameBorderMode(0);

 c1->cd();

 TPad* canvas_1 = new TPad("canvas"+name, "canvas"+name,0,0.23,1.0,0.98);
 canvas_1 ->Draw();
 canvas_1->SetFillColor(0);
 canvas_1->SetBorderMode(0);
 canvas_1->SetBorderSize(2);
 canvas_1->SetTopMargin(0.065);
 canvas_1->SetFrameBorderMode(0);
 canvas_1 ->cd();
 canvas_1->SetLogy(log);
 
 //if (hist_data->GetMaximum() > stack->GetMaximum() ) stack->SetMaximum( hist_data->GetMaximum()*1.1) ;

 if (!log) stack->SetMaximum( stack->GetMaximum()*2);
 else stack->SetMaximum( stack->GetMaximum()*5000);

 if (name=="trk_multi_sel" && log) stack->SetMaximum( stack->GetMaximum()*50000);
 
 if (log) stack->SetMinimum(0.1);
 else  stack->SetMinimum(0.);


 float xmov=1;
 if (move_legend==2)  
 {
    if (log)  xmov=5;
    else xmov=1.1;
    stack->SetMaximum( xmov*stack->GetMaximum() ) ;
 }

 if (log && name == "CSVv2") stack->SetMaximum(2*10e7);

 stack    ->Draw("hist");  
  
 stack    ->GetHistogram()->SetTitleSize(0.08,"Y");
 stack    ->GetHistogram()->SetTitleOffset(0.81,"Y"); 

 stack    ->GetHistogram()->GetXaxis()->SetLabelSize(0);
 stack    ->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
 stack    ->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
 stack    ->GetHistogram()->GetXaxis()->SetTitle(name);

 if(name == "CSV")        stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.02");
 else if(name == "CSVv2") stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.02");
 else if(name == "cMVAv2") stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.04");
 else if(name == "track_IPs") stack    ->GetHistogram()->GetYaxis()->SetTitle("Tracks/0.7");
 else if(name == "sv_flight3DSig") stack    ->GetHistogram()->GetYaxis()->SetTitle("SVs/1.6");
 else if(name == "tagvarCSV_vertexmass_cat0") stack    ->GetHistogram()->GetYaxis()->SetTitle("SVs/0.16 GeV");
 else if(name == "JP") stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.05");
 else if(name == "JBP") stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.16");
 else if(name == "SoftMu") stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.02");
 else if(name == "SoftEl") stack    ->GetHistogram()->GetYaxis()->SetTitle("Jets/0.02");
 else                     stack    ->GetHistogram()->GetYaxis()->SetTitle("entries");


 hist_data->Draw("same E1X0");

 TLegend* leg;
 if (move_legend==1) 
 {
   leg =  new TLegend(0.35,0.15,0.70,0.42, NULL, "brNDC");
 }
 else if (move_legend==3) 
 {
   leg =  new TLegend(0.35,0.63,0.70,0.90, NULL, "brNDC");
 }
 else leg =  new TLegend(0.65,0.60,0.95,0.90, NULL, "brNDC");

  
  //Legend
 leg->AddEntry(hist_data,     "Data"     ,  "ep");
 leg->AddEntry(hist_b,        "b"        ,   "f");
 leg->AddEntry(hist_c,        "c"        ,   "f");
 leg->AddEntry(hist_l,        "udsg"     ,   "f");
 leg->AddEntry(hist_pu,       "pile-up"  ,   "f");
 //leg->AddEntry(hist_gsplit,   "b from gluon splitting"     ,"f");

 leg->SetBorderSize(1);
 leg->SetTextFont(62);
 leg->SetLineColor(1);
 leg->SetLineStyle(1);
 leg->SetLineWidth(1);
 leg->SetFillStyle(1001);
 leg->SetBorderSize(0);
 
 leg->SetFillColor(0);
 leg->Draw();

  
 TLatex *latex = new TLatex(0.20, 0.89, CMStitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(62); //22
 latex->SetTextSize(0.065);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.82, Preliminarytitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(52); //22
 latex->SetTextSize(0.05681);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.98, 0.95125, Lumititle);
 latex->SetNDC();
 latex->SetTextAlign(31);
 latex->SetTextFont(42); //22
 latex->SetTextSize(0.04875);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.74, Selectiontitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(42); //22
 latex->SetTextSize(0.055);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.68, Jettypetitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(42); //22
 latex->SetTextSize(0.055);
 latex->SetLineWidth(2);
 latex->Draw();

 canvas_1->Modified();

 c1->cd();  
  
 TPad* canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.32);
 canvas_2->Draw();
 canvas_2->SetFillColor(0);
 canvas_2->SetBorderMode(0);
 canvas_2->SetBorderSize(2);
 canvas_2->SetGridy();
 canvas_2->SetBottomMargin(0.31);
 canvas_2->SetFrameBorderMode(0);

 canvas_2->cd();
 gPad->SetBottomMargin(0.375);
 gPad->SetGridy();
  
 histo_ratio->SetMarkerStyle(20);
 histo_ratio->SetMarkerSize(0.75);
 histo_ratio->SetLineWidth(2);
 
 histo_ratio->GetYaxis()->SetTitle("Data/MC");
 histo_ratio->GetYaxis()->SetTitleFont(62);
 histo_ratio->GetYaxis()->SetLabelFont(62);
 histo_ratio->GetXaxis()->SetTitle(histotitle);
 histo_ratio->GetXaxis()->SetTitleFont(62);
 histo_ratio->GetXaxis()->SetLabelFont(62);
 histo_ratio->GetYaxis()->SetNdivisions( 505 );

 histo_ratio->SetTitleOffset(0.9, "X");
 histo_ratio->SetTitleOffset(0.31, "Y");
  
 double labelsizex=0.12;
 double labelsizey=0.12; 
 double titlesizex=0.15;
 double titlesizey=0.14;
  
   
 histo_ratio->GetXaxis()->SetLabelSize( labelsizex );
 histo_ratio->GetXaxis()->SetTitleSize( titlesizex );
 histo_ratio->GetYaxis()->SetLabelSize( labelsizey );
 histo_ratio->GetYaxis()->SetTitleSize( titlesizey );

 histo_ratio->SetMinimum(0.4);
 histo_ratio->SetMaximum(1.6);
 histo_ratio->Draw("E1X0");

 c1->cd();  
  
 TString name_plot1=name+"_Linear"+format1; 
 if(log) name_plot1=name+"_Log"+format1;
 TString name_plot2=name+"_Linear"+format2; 
 if(log) name_plot2=name+"_Log"+format2;
 c1->SaveAs("ttbar/"+output+name_plot1);
 c1->SaveAs("ttbar/"+output+name_plot2);
}

//--------------------------

void DrawTTbar(TString name, TString histotitle, bool log, int move_legend)
{

 
 if(name.SubString("afterJetSel") == "afterJetSel") Selectiontitle= "e#mu channel, = 2 jets";
 else                                               Selectiontitle= "e#mu channel, #geq 2 jets";

 Lumititle= "2.5 fb^{-1} (13 TeV, 25ns)";
 Jettypetitle="AK4 jets (p_{T} > 30 GeV)";

 TH1D* hist_ttbar;
 TH1D* hist_dy;
 TH1D* hist_st;
 TH1D* hist_ww;
 TH1D* hist_wz;
 TH1D* hist_zz;
 TH1D* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_ttbar     = (TH1D*)gROOT->FindObject(name+"_ttbar");
 hist_dy        = (TH1D*)gROOT->FindObject(name+"_dy");
 hist_st        = (TH1D*)gROOT->FindObject(name+"_st");
 hist_ww        = (TH1D*)gROOT->FindObject(name+"_ww");
 hist_wz        = (TH1D*)gROOT->FindObject(name+"_wz");
 hist_zz        = (TH1D*)gROOT->FindObject(name+"_zz");
 hist_data      = (TH1D*)gROOT->FindObject(name+"_data");

 if (bOverflow) {
  OverFlowBinFix(hist_ttbar);
  OverFlowBinFix(hist_dy);
  OverFlowBinFix(hist_st);
  OverFlowBinFix(hist_ww);
  OverFlowBinFix(hist_wz);
  OverFlowBinFix(hist_zz);
  OverFlowBinFix(hist_data);
 }



 TH1D* histo_tot = (TH1D*) hist_ttbar->Clone();
 histo_tot->Sumw2();
 histo_tot ->Add(hist_dy);
 histo_tot ->Add(hist_st);
 histo_tot ->Add(hist_ww);
 histo_tot ->Add(hist_wz);
 histo_tot ->Add(hist_zz);

 

 float scale_f = 1.;
 //float scale_f = (hist_data->Integral())/(histo_tot->Integral());

 hist_ttbar   ->Scale(scale_f);
 hist_dy      ->Scale(scale_f);
 hist_st      ->Scale(scale_f);
 hist_ww      ->Scale(scale_f);
 hist_wz      ->Scale(scale_f);
 hist_zz      ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);

  
 TH1D* histo_ratio;
 histo_ratio = (TH1D*) hist_data->Clone();
 histo_ratio->SetName("histo_ratio");
 histo_ratio->SetTitle("");
  
 histo_ratio->Divide(histo_tot);
  
 hist_data  ->SetLineWidth(2);
 hist_data  ->SetMarkerStyle(20);  
 hist_data  ->SetMarkerSize(0.75); 

 hist_dy     ->SetFillColor(kAzure-2);
 hist_ttbar  ->SetFillColor(kRed+1);
 hist_st     ->SetFillColor(kMagenta);
 hist_ww     ->SetFillColor(kGreen+4);
 hist_wz     ->SetFillColor(kOrange+2);
 hist_zz     ->SetFillColor(kBlue+2);
  
 THStack *stack = new THStack("stack","stack");
  
 stack      ->Add(hist_ttbar);
 stack      ->Add(hist_dy);
 stack      ->Add(hist_st);
 stack      ->Add(hist_ww);
 stack      ->Add(hist_wz);
 stack      ->Add(hist_zz);

 setTDRStyle();
 gStyle->SetErrorX(0.);

 TCanvas *c1 = new TCanvas(name, name,10,32,782,552);
 c1->SetFillColor(10);
 c1->SetBorderMode(0);
 c1->SetBorderSize(2);
 c1->SetFrameFillColor(0);
 c1->SetFrameBorderMode(0);

 c1->cd();

 TPad* canvas_1 = new TPad("canvas"+name, "canvas"+name,0,0.23,1.0,0.98);
 canvas_1 ->Draw();
 canvas_1->SetFillColor(0);
 canvas_1->SetBorderMode(0);
 canvas_1->SetBorderSize(2);
 canvas_1->SetTopMargin(0.065);
 canvas_1->SetFrameBorderMode(0);
 canvas_1 ->cd();
 canvas_1->SetLogy(log);
 
 //if (hist_data->GetMaximum() > stack->GetMaximum() ) stack->SetMaximum( hist_data->GetMaximum()*1.1) ;

 if (!log) stack->SetMaximum( stack->GetMaximum()*2);
 else stack->SetMaximum( stack->GetMaximum()*5000);

 if (log) stack->SetMinimum(0.1);
 else  stack->SetMinimum(0.);

 stack    ->Draw("hist");  
  
 stack    ->GetHistogram()->SetTitleSize(0.08,"Y");
 stack    ->GetHistogram()->SetTitleOffset(0.81,"Y"); 

 stack    ->GetHistogram()->GetXaxis()->SetLabelSize(0);
 stack    ->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
 stack    ->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
 stack    ->GetHistogram()->GetXaxis()->SetTitle(name);
 stack    ->GetHistogram()->GetYaxis()->SetTitle("Events");


 hist_data->Draw("same E1X0");

 TLegend* leg;
 leg =  new TLegend(0.70,0.50,0.95,0.90, NULL, "brNDC");

  
  //Legend
 leg->AddEntry(hist_data,     "Data",                       "p");
 leg->AddEntry(hist_ttbar,    "t#bar{t}"           ,         "f");
 leg->AddEntry(hist_dy,       "DY"     ,"f");
 leg->AddEntry(hist_st,       "tW"           ,         "f");
 leg->AddEntry(hist_ww,       "WW"           ,         "f");
 leg->AddEntry(hist_wz,       "WZ"           ,         "f");
 leg->AddEntry(hist_zz,       "ZZ"           ,         "f");

 leg->SetBorderSize(1);
 leg->SetTextFont(62);
 leg->SetLineColor(1);
 leg->SetLineStyle(1);
 leg->SetLineWidth(1);
 leg->SetFillStyle(1001);
 leg->SetBorderSize(0);
 
 leg->SetFillColor(0);
 leg->Draw();

  
 TLatex *latex = new TLatex(0.20, 0.89, CMStitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(62); //22
 latex->SetTextSize(0.065);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.82, Preliminarytitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(52); //22
 latex->SetTextSize(0.05681);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.98, 0.95125, Lumititle);
 latex->SetNDC();
 latex->SetTextAlign(31);
 latex->SetTextFont(42); //22
 latex->SetTextSize(0.04875);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.74, Selectiontitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(42); //22
 latex->SetTextSize(0.055);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.68, Jettypetitle);
 latex->SetNDC();
 latex->SetTextAlign(13);
 latex->SetTextFont(42); //22
 latex->SetTextSize(0.055);
 latex->SetLineWidth(2);
 latex->Draw();

 if(name.SubString("SFapplied") == "SFapplied") 
 {
        latex = new TLatex(0.20, 0.62, "closure test");
        latex->SetNDC();
        latex->SetTextAlign(13);
        latex->SetTextFont(42); //22
        latex->SetTextSize(0.055);
        latex->SetLineWidth(2);
        latex->Draw();
 } 

 canvas_1->Modified();

 c1->cd();  

 TPad* canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.32);
 canvas_2->Draw();
 canvas_2->SetFillColor(0);
 canvas_2->SetBorderMode(0);
 canvas_2->SetBorderSize(2);
 canvas_2->SetGridy();
 canvas_2->SetBottomMargin(0.31);
 canvas_2->SetFrameBorderMode(0);

 canvas_2->cd();
 gPad->SetBottomMargin(0.375);
 gPad->SetGridy();
  
 histo_ratio->SetMarkerStyle(20);
 histo_ratio->SetMarkerSize(0.75);
 histo_ratio->SetLineWidth(2);
 
 histo_ratio->GetYaxis()->SetTitle("Data/MC");
 histo_ratio->GetYaxis()->SetTitleFont(62);
 histo_ratio->GetYaxis()->SetLabelFont(62);
 histo_ratio->GetXaxis()->SetTitle(histotitle);
 histo_ratio->GetXaxis()->SetTitleFont(62);
 histo_ratio->GetXaxis()->SetLabelFont(62);
 histo_ratio->GetYaxis()->SetNdivisions( 505 );

 histo_ratio->SetTitleOffset(0.9, "X");
 histo_ratio->SetTitleOffset(0.31, "Y");
  
 double labelsizex=0.12;
 double labelsizey=0.12; 
 double titlesizex=0.15;
 double titlesizey=0.14;
  
   
 histo_ratio->GetXaxis()->SetLabelSize( labelsizex );
 histo_ratio->GetXaxis()->SetTitleSize( titlesizex );
 histo_ratio->GetYaxis()->SetLabelSize( labelsizey );
 histo_ratio->GetYaxis()->SetTitleSize( titlesizey );

 histo_ratio->SetMinimum(0.4);
 histo_ratio->SetMaximum(1.6);
 histo_ratio->Draw("E1X0");

  
 c1->cd();  
  
 TString name_plot1="ttbar_"+name+"_Linear"+format1; 
 if(log) name_plot1="ttbar_"+name+"_Log"+format1;
 TString name_plot2="ttbar_"+name+"_Linear"+format2; 
 if(log) name_plot2="ttbar_"+name+"_Log"+format2;
 c1->SaveAs("ttbar/Commissioning_plots/"+name_plot1);
 c1->SaveAs("ttbar/Commissioning_plots/"+name_plot2);


}



//--------------------------


void DrawTagRate(TString name, TString histotitle, bool log)
{


 TH1D* hist_b;
 TH1D* hist_c;
 TH1D* hist_gsplit;
 TH1D* hist_l;
 TH1D* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_b         = (TH1D*)gROOT->FindObject(name+"_b");
 hist_c         = (TH1D*)gROOT->FindObject(name+"_c");
 hist_gsplit    = (TH1D*)gROOT->FindObject(name+"_bfromg");
 hist_l         = (TH1D*)gROOT->FindObject(name+"_l");
 hist_data      = (TH1D*)gROOT->FindObject(name+"_data");
 
/*
 if (bOverflow) {
  OverFlowBinFix(hist_b);
  OverFlowBinFix(hist_c);
  OverFlowBinFix(hist_gsplit);
  OverFlowBinFix(hist_l);
  OverFlowBinFix(hist_data);
 }
*/


 TH1D* histo_tot = (TH1D*) hist_b->Clone();
 histo_tot->Sumw2();
 histo_tot ->Add(hist_c);
 histo_tot ->Add(hist_gsplit);
 histo_tot ->Add(hist_l);  
  
 //float scale_f = 1.;
 float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_gsplit->Integral() + hist_l->Integral());

 hist_b       ->Scale(scale_f);
 hist_c       ->Scale(scale_f);
 hist_gsplit  ->Scale(scale_f);
 hist_l       ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);
 
 int   nbinx=hist_data->GetNbinsX();
 float minx=hist_data->GetXaxis()->GetXmin();
 float maxx=hist_data->GetXaxis()->GetXmax();
 
 TH1D * TagRate_Data = new TH1D ("TagRate_Data",histotitle,nbinx,minx,maxx);
 TH1D * TagRate_MC   = new TH1D ("TagRate_MC",histotitle,nbinx,minx,maxx);
 TH1D * TagRate_MC_b = new TH1D ("TagRate_MC_b",histotitle,nbinx,minx,maxx);
 TH1D * TagRate_MC_c = new TH1D ("TagRate_MC_c",histotitle,nbinx,minx,maxx);
 TH1D * TagRate_MC_udsg = new TH1D ("TagRate_MC_udsg",histotitle,nbinx,minx,maxx);
 TH1D * TagRate_MC_gspl = new TH1D ("TagRate_MC_gspl",histotitle,nbinx,minx,maxx);

 int nbin_max= hist_data->GetNbinsX();

 for (int ii=0; ii<nbin_max; ii++) {
    float totdata=hist_data->Integral(0,nbin_max+1);
    float totmc=histo_tot->Integral(0,nbin_max+1);

    float val = hist_data->Integral(ii+1,nbin_max+1) / totdata;
    float err = sqrt(totdata *val*(1-val))/ totdata;
      
    float valMC = histo_tot->Integral(ii+1,nbin_max+1)/ totmc;
    float errMC = sqrt(totmc *valMC*(1-valMC))/ totmc;

    TagRate_Data->SetBinContent(ii+1, val);
    TagRate_Data->SetBinError(ii+1, err);
    TagRate_MC->SetBinContent(ii+1,   histo_tot->Integral(ii+1,nbin_max+1) / totmc);
    TagRate_MC->SetBinError(ii+1, errMC  );
    TagRate_MC_b->SetBinContent(ii+1, hist_b->Integral(ii+1,nbin_max+1) / totmc);
    TagRate_MC_c->SetBinContent(ii+1, hist_c->Integral(ii+1,nbin_max+1) / totmc);
    TagRate_MC_udsg->SetBinContent(ii+1, hist_l->Integral(ii+1,nbin_max+1) / totmc);
    TagRate_MC_gspl->SetBinContent(ii+1, hist_gsplit->Integral(ii+1,nbin_max+1) / totmc);

  }

   double titleoffsety=0.2;
   double titlesizex=0.17;
   double titlesizey=0.08;
   double labelsizex=0.14;
   double labelsizey=0.12;

   TagRate_Data  ->GetYaxis()->SetLabelSize(labelsizey);
   TagRate_Data  ->GetYaxis()->SetTitleSize(titlesizey);
   TagRate_Data  ->GetYaxis()->SetTitleOffset(titleoffsety);
  
   TagRate_MC     ->GetYaxis()->SetLabelSize(labelsizey);
   TagRate_MC     ->GetYaxis()->SetTitleSize(titlesizey);
   TagRate_MC     ->GetYaxis()->SetTitleOffset(titleoffsety);
   TagRate_MC_b     ->GetYaxis()->SetLabelSize(labelsizey);
   TagRate_MC_b     ->GetYaxis()->SetTitleSize(titlesizey);
   TagRate_MC_b     ->GetYaxis()->SetTitleOffset(titleoffsety);


  // MAKE DATA/MC RATIO
  
  TH1D* histo_ratio;
  histo_ratio = (TH1D*) TagRate_Data->Clone();
  histo_ratio->SetName("histo0_ratio");
  histo_ratio->SetTitle("");
  histo_ratio->Divide(TagRate_MC);


  // SET COLORS
  TagRate_MC->SetLineColor(2);
  TagRate_MC_b->SetFillColor(2);
  TagRate_MC_c->SetFillColor(8);
  TagRate_MC_gspl->SetFillColor(7);
  TagRate_MC_udsg->SetFillColor(4);

  // DO STACK
  THStack* hs= new THStack();

  if (b_ordering)
  {
    hs->Add(TagRate_MC_b);
    hs->Add(TagRate_MC_gspl);  
    hs->Add(TagRate_MC_c);
    hs->Add(TagRate_MC_udsg);
  }
  else 
  {
    hs->Add(TagRate_MC_udsg);
    hs->Add(TagRate_MC_c);
    hs->Add(TagRate_MC_gspl);  
    hs->Add(TagRate_MC_b);
  }
  
  // SET COSMETICS
  TagRate_Data->SetMarkerStyle(20);
  TagRate_Data->SetMarkerSize(0.75);
  TagRate_Data->GetXaxis()->SetTitle();  


  gStyle->SetOptTitle(0);

  // CREATE CANVAS
  
  TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
  c1->cd();

  // FIRST MC & DATA
  TPad *c1_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
  c1_1->Draw();
  c1_1->cd(); 

  if (TagRate_Data->GetMaximum() > hs->GetMaximum() ) 
  {
    hs->SetMaximum(TagRate_Data->GetMaximum()*1.1 );
  }
  hs->Draw("hist");

  hs->GetHistogram()->SetTitleSize(0.08,"Y");
  hs->GetHistogram()->SetTitleOffset(0.55,"Y"); 
  hs->GetHistogram()->GetYaxis()->SetTitle("Tag Rate");
  hs->GetHistogram()->GetXaxis()->SetTitle(histotitle);

  TagRate_Data->Draw("e same");  
    // ADD LEGEND
  TLegend* qw = 0;
  qw = new TLegend(0.6,0.73,0.95,1.);
  qw->AddEntry(TagRate_Data,        "e#mu ttbar data"                     ,"p");
  qw->AddEntry(TagRate_MC_b,        "b quark"                  ,"f");
  qw->AddEntry(TagRate_MC_gspl,     "b from gluon splitting"   ,"f");
  qw->AddEntry(TagRate_MC_c,        "c quark"                  ,"f");
  qw->AddEntry(TagRate_MC_udsg,     "uds quark or gluon"     ,"f");
 
  qw->SetFillColor(0);
  qw->Draw();

  
 TLatex *latex = new TLatex(0.20, 0.89, CMStitle);
 latex->SetNDC();
 latex->SetTextSize(0.065);
 latex->SetTextFont(62); //22
 latex->SetTextAlign(13);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.20, 0.82, Preliminarytitle);
 latex->SetNDC();
 latex->SetTextSize(0.05681);
 latex->SetTextFont(52); //22
 latex->SetTextAlign(13);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.98, 0.95125, Lumititle);
 latex->SetNDC();
 latex->SetTextSize(0.04875);
 latex->SetTextFont(42); //22
 latex->SetTextAlign(31);
 latex->SetLineWidth(2);
 latex->Draw();

 latex = new TLatex(0.55, 0.98, Selectiontitle);
 latex->SetNDC();
 latex->SetTextSize(0.08);
 latex->SetTextFont(42); //22
 latex->SetTextAlign(13);
 latex->SetLineWidth(2);
 latex->SetX(0.63);
 latex->SetY(0.90);
 latex->Draw();

  c1->cd();  
  
  TPad* canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.32);
  canvas_2->Draw();
  canvas_2->cd();

  gPad->SetBottomMargin(0.375);
  gPad->SetGridy();
  
 histo_ratio->SetMarkerStyle(20);
 histo_ratio->SetMarkerSize(0.75);
 histo_ratio->SetLineWidth(2);
 
 histo_ratio->GetYaxis()->SetTitle("Data/MC");
 histo_ratio->GetXaxis()->SetTitle(histotitle);
 histo_ratio->GetYaxis()->SetNdivisions( 505 );

 histo_ratio->GetXaxis()->SetLabelSize( labelsizex);
 histo_ratio->GetXaxis()->SetTitleSize( titlesizex );

 histo_ratio->SetMinimum(0.4);
 histo_ratio->SetMaximum(1.6);
 histo_ratio->Draw("E1X0");


  c1->cd();
  
  TString name_plot1="tag_"+name+"_Linear"+format1;
  if(log) name_plot1="tag_"+name+"_Log"+format1;
  TString name_plot2="tag_"+name+"_Linear"+format2;
  if(log) name_plot2="tag_"+name+"_Log"+format2;
  c1->SaveAs("ttbar/Commissioning_plots/"+name_plot1);
  c1->SaveAs("ttbar/Commissioning_plots/"+name_plot2);
}


//--------------------------


void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log)
{

 TH2D* hist_b;
 TH2D* hist_c;
 TH2D* hist_gsplit;
 TH2D* hist_l;
 TH2D* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_b         = (TH2D*)gROOT->FindObject(name+"_b");
 hist_c         = (TH2D*)gROOT->FindObject(name+"_c");
 hist_gsplit    = (TH2D*)gROOT->FindObject(name+"_bfromg");
 hist_l         = (TH2D*)gROOT->FindObject(name+"_l");
 hist_data      = (TH2D*)gROOT->FindObject(name+"_data");
 

 TH2D* histo_tot = (TH2D*) hist_b->Clone();
 histo_tot ->Add(hist_c);
 histo_tot ->Add(hist_gsplit);
 histo_tot ->Add(hist_l); 
 
 //float scale_f = 1.;
 float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_gsplit->Integral() + hist_l->Integral());

 hist_b       ->Scale(scale_f);
 hist_c       ->Scale(scale_f);
 hist_gsplit  ->Scale(scale_f);
 hist_l       ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);

  
 TProfile * pro_mc = histo_tot->ProfileX(name+"_tot");
 TProfile * pro_mc_b = hist_b->ProfileX();
 TProfile * pro_mc_c = hist_c->ProfileX();
 TProfile * pro_mc_udsg = hist_l->ProfileX();
 TProfile * pro_mc_gspl = hist_gsplit->ProfileX();
 TProfile * pro_data = hist_data->ProfileX();
 
   // SET COLORS
  pro_mc->SetLineColor(1);
  pro_mc_b->SetLineColor(2);
  pro_mc_c->SetLineColor(8);
  pro_mc_udsg->SetLineColor(4);
  pro_mc_gspl->SetLineColor(7);

  pro_data->SetMarkerStyle(20);
  pro_data->SetMarkerSize(0.75);
  
  pro_mc_gspl->GetXaxis()->SetTitle(titleX);  
  pro_mc_gspl->GetYaxis()->SetTitle(titleY);  
  
  
  // SET COSMETICS
  pro_data->SetMarkerStyle(20);
  pro_data->SetMarkerSize(0.75);
  //pro_mc_gspl->GetXaxis()->SetTitle();  
  //pro_mc_gspl->GetYaxis()->SetTitle();  


  gStyle->SetOptTitle(0);


  Double_t titleoffsetx=0.8;
  Double_t titleoffsety=0.8;
  Double_t titlesizex=0.05;
  Double_t titlesizey=0.05;
  Double_t labelsizex=0.035;
  Double_t labelsizey=0.035;

  pro_data->GetYaxis()->SetLabelSize(labelsizey);
  pro_data->GetYaxis()->SetTitleSize(titlesizey);
  pro_data->GetYaxis()->SetTitleOffset(titleoffsety);
  pro_mc->GetYaxis()->SetLabelSize(labelsizey);
  pro_mc->GetYaxis()->SetTitleSize(titlesizey);
  pro_mc->GetYaxis()->SetTitleOffset(titleoffsety);
  pro_mc_b->GetYaxis()->SetLabelSize(labelsizey);
  pro_mc_b->GetYaxis()->SetTitleSize(titlesizey);
  pro_mc_b->GetYaxis()->SetTitleOffset(titleoffsety);
  pro_mc_c->GetYaxis()->SetLabelSize(labelsizey);
  pro_mc_c->GetYaxis()->SetTitleSize(titlesizey);
  pro_mc_c->GetYaxis()->SetTitleOffset(titleoffsety);

  pro_mc_gspl->GetYaxis()->SetLabelSize(labelsizey);
  pro_mc_gspl->GetYaxis()->SetTitleSize(titlesizey);
  pro_mc_gspl->GetYaxis()->SetTitleOffset(titleoffsety);

  pro_mc_udsg->GetYaxis()->SetLabelSize(labelsizey);
  pro_mc_udsg->GetYaxis()->SetTitleSize(titlesizey);
  pro_mc_udsg->GetYaxis()->SetTitleOffset(titleoffsety);



  pro_data->GetXaxis()->SetLabelSize(labelsizex);
  pro_data->GetXaxis()->SetTitleSize(titlesizex);
  pro_data->GetXaxis()->SetTitleOffset(titleoffsetx);
  pro_mc->GetXaxis()->SetLabelSize(labelsizex);
  pro_mc->GetXaxis()->SetTitleSize(titlesizex);
  pro_mc->GetXaxis()->SetTitleOffset(titleoffsetx);
  pro_mc_b->GetXaxis()->SetLabelSize(labelsizex);
  pro_mc_b->GetXaxis()->SetTitleSize(titlesizex);
  pro_mc_b->GetXaxis()->SetTitleOffset(titleoffsetx);
  pro_mc_c->GetXaxis()->SetLabelSize(labelsizex);
  pro_mc_c->GetXaxis()->SetTitleSize(titlesizex);
  pro_mc_c->GetXaxis()->SetTitleOffset(titleoffsetx);

  pro_mc_gspl->GetXaxis()->SetLabelSize(labelsizex);
  pro_mc_gspl->GetXaxis()->SetTitleSize(titlesizex);
  pro_mc_gspl->GetXaxis()->SetTitleOffset(titleoffsetx);

  pro_mc_udsg->GetXaxis()->SetLabelSize(labelsizex);
  pro_mc_udsg->GetXaxis()->SetTitleSize(titlesizex);
  pro_mc_udsg->GetXaxis()->SetTitleOffset(titleoffsetx);

  TCanvas *canvas = new TCanvas("c1", "c1",10,32,782,552);
  canvas->cd();

  float maxhist= pro_mc_gspl->GetMaximum();
  if (pro_mc_b->GetMaximum() > maxhist) maxhist = pro_mc_b->GetMaximum()*1.1;
  if (pro_mc_c->GetMaximum() > maxhist) maxhist = pro_mc_c->GetMaximum()*1.1;
  if (pro_mc_udsg->GetMaximum() > maxhist) maxhist = pro_mc_udsg->GetMaximum()*1.1;
  if (pro_mc->GetMaximum() > maxhist) maxhist = pro_mc->GetMaximum()*1.1;
  if (pro_data->GetMaximum() > maxhist) maxhist = pro_data->GetMaximum()*1.1;
  if (maxhist> pro_mc_gspl->GetMaximum()) pro_mc_gspl->SetMaximum(maxhist);

  pro_mc_gspl->Draw("hist");
  pro_mc_b->Draw("hist,same");
  pro_mc_c->Draw("hist,same");
  pro_mc_udsg->Draw("hist,same");
  pro_mc->Draw("hist,same");
  pro_data->Draw("e,same");

  TLegend* qw = 0;
  qw =  new TLegend(0.6,0.73,0.95,1.);
  qw->AddEntry(pro_data,        "e#mu ttbar data"                   ,"p");
  qw->AddEntry(pro_mc,          "total "                 ,"l");
  qw->AddEntry(pro_mc_b,        "b quark"                ,"l");
  qw->AddEntry(pro_mc_gspl,     "b from gluon splitting" ,"l");
  qw->AddEntry(pro_mc_c,        "c quark"                ,"l");
  qw->AddEntry(pro_mc_udsg,     "uds quark or gluon"     ,"l");

  qw->SetFillColor(0);
  qw->Draw();



  TString name_plot1=name+"_Linear"+format1;
  if(log) name_plot1=name+"_Log"+format1;
  TString name_plot2=name+"_Linear"+format2;
  if(log) name_plot2=name+"_Log"+format2;
  canvas->SaveAs("ttbar/Commissioning_plots/"+name_plot1);
  canvas->SaveAs("ttbar/Commissioning_plots/"+name_plot2);

}

//------------

void OverFlowBinFix(TH1D* histo)
{

  Float_t val, errval;

  val=histo->GetBinContent(1)+histo->GetBinContent(0);
  errval=0;


  errval=pow(histo->GetBinError(0),2)+pow(histo->GetBinError(1),2);
  errval=sqrt(errval);
  histo->SetBinContent(1,val);
  histo->SetBinError(1,errval);
  histo->SetBinContent(0,0);
  histo->SetBinError(0,0);

  Int_t highbin=histo->GetNbinsX();

  val=histo->GetBinContent(highbin)+histo->GetBinContent(highbin+1);
  errval=0;


  errval=pow(histo->GetBinError(highbin),2)+pow(histo->GetBinError(highbin+1),2);
  errval=sqrt(errval);
  histo->SetBinContent(highbin,val);
  histo->SetBinError(highbin,errval);
  histo->SetBinContent(highbin+1,0);
  histo->SetBinError(highbin+1,0);
}


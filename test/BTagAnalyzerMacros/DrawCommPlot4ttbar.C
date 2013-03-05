#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
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

using namespace std;

TString filename="output_test.root";
TString title= "CMS 2012 preliminary, #sqrt{s} = 8 TeV,  12.2 fb^{-1}";
TString format=".gif"; // .png or .pdf
bool bOverflow=true;
bool b_ordering=false;


void Draw(TString name, TString histotitle, bool log, int move_legend=0);
void DrawTTbar(TString name, TString histotitle, bool log, int move_legend=0);
void DrawTagRate(TString name, TString histotitle, bool log);
void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log);
void OverFlowBinFix(TH1F* );

//--------

void DrawCommPlot(bool Draw_track_plots, bool Draw_Nminus1_plots, bool Draw_sv_plots, bool Draw_muons_plots, bool Draw_discriminator_plots , bool Draw_tagRate_plots, bool Draw_2D_plots){

  TString action = "mkdir Commissioning_plots";
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



//Draw("jet_multi"    ,"number of jets",1);	    
Draw("jet_pt_all"   ,"pT of all jets",1);	    
Draw("jet_pt_sv"    ,"pT of jets containing a SV",1);
Draw("jet_eta"      ,"eta of all jets",	1, move);    
Draw("jet_phi"      ,"phi of all jets",1, move);

if (Draw_track_plots){	      
  Draw("track_multi"  ,      "number of tracks in the jets",1, move+1);		   
  Draw("trk_multi_sel"  ,    "number of selected tracks in the jets",1);	   
  Draw("track_chi2"   ,      "normalized #chi^{2} of the tracks"	,1);	   
  Draw("track_nHit" ,      "number of hits",1);		   
  Draw("track_HPix"   ,      "number of hits in the Pixel",1);		   
  Draw("track_IPs"    ,      "3D IP significance of all tracks",1);	   
  Draw("track_IPs1tr" ,      "3D IP significance of the first track",1, move+1);	   
  Draw("track_IPs2tr" ,      "3D IP significance of the second track",1);     
  Draw("track_IP"     ,      "3D IP of all tracks",1, move+1);			   
  Draw("track_IP1tr"  ,      "3D IP of the first track"	,1, move+1);	   
  Draw("track_IP2tr"  ,      "3D IP of the second track",1, move+1);		   
  Draw("track_IP2Ds"	,    "2D IP significance of all tracks"	,1);   
  Draw("track_IP2Ds1tr" ,    "2D IP significance of the first track",1, move+1);	   
  Draw("track_IP2Ds2tr" ,    "2D IP significance of the second track" ,1);    
  Draw("track_IP2D"    ,     "2D IP of all tracks",1);			   
  Draw("track_IP2D1tr" ,     "2D IP of the first track"	,1, move+1);	   
  Draw("track_IP2D2tr" ,     "2D IP of the second track",1);		   
  Draw("track_IP2Derr1tr" ,  "2D IP error of the first track",1);  	   
  Draw("track_IPerr1tr"   ,  "3D IP error of the first track" ,1); 	   
  Draw("track_IP2Derr2tr" ,  "2D IP error of the second track" ,1);	   
  Draw("track_IPerr2tr"   ,  "3D IP error of the second track" ,1);	   
  Draw("track_IP2Derr" ,     "2D IP error of all tracks",1);				   
  Draw("track_IPerr"   ,     "3D IP error of all tracks",1);				   
  Draw("track_len"     ,     "Track decay length",1);				   
  Draw("track_dist"    ,     "Track distance to the jet axis"	,1);	   
  Draw("track_dz"     ,      "Track transverse IP",1);				   
  Draw("track_isfromSV",     "Track is from SV",1);			   
  Draw("track_pt"     ,      "pT of all the tracks",1);			   
	
  Draw("track_IPs3tr" ,    "3D IP significance of the third track",1);	
  Draw("track_IP3tr"  ,      "3D IP of the third track"	,1);	
  Draw("track_IPerr3tr"   ,  "3D IP error of the third track" ,1); 	
  Draw("track_IP2Ds3tr" ,    "2D IP significance of the second track" ,1);
  Draw("track_IP2D3tr" ,     "2D IP of the third track"	,1);	
  Draw("track_IP2Derr3tr" ,  "2D IP error of the third track",1);  	
  }  
if (Draw_Nminus1_plots){  
  Draw("track_chi2_cut"    ,"#chi^{2} normalized  without this cut",1);	    
  Draw("track_nHit_cut"  ,"number of hits without this cut",1);
  Draw("track_HPix_cut"    ,"number of hits in the Pixel without this cut",1);  
  Draw("track_len_cut"     ,"decay length without this cut",1);		    
  Draw("track_dist_cut"    ,"distance to the jet axis without this cut" ,1);  
  Draw("track_dz_cut"     ,"transverse IP without this cut",1);		    
  Draw("track_pt_cut"	   ,"pT without this cut",1);
  }
if (Draw_sv_plots){

  Draw("sv_mass","SV mass",1);
  Draw("sv_multi","nr. of SV",1);
  Draw("sv_chi2norm","SV norm. #chi^{2}",1, move+1);
  Draw("sv_deltaR_sumJet","#Delta R between the jet and the SV",1);
  Draw("sv_deltaR_sumDir","#Delta R between the jet direction and the SV",1);
  Draw("sv_en_ratio","SV energy ratio",1, move+1);	
  Draw("sv_aboveC","IP2D of the first track above the charm threshold",1, move+1);	
  Draw("sv_pt","SV p_{T}",1); 	
  Draw("sv_eta","SV #eta",1, move);	
  Draw("sv_phi","SV #phi",1, move);	
  Draw("sv_flightSig2D","SV 2D flight distance significance",1, move+1);
  Draw("sv_flight2D","SV 2D flight distance",1);	
  Draw("sv_flight3D","SV 3D flight distance",1);	
  Draw("sv_flight3DSig","SV 3D flight distance significance",1, move+1);
  Draw("sv_multi_0","nr. of SV including bin 0",1); 
  //Draw("sv_tot_charge ","SV charge",1);
  Draw("svnTrk","nr. of tracks from a SV",1);	
  Draw("svnTrk_firstVxt","nr. of tracks from the first SV",1);
  Draw("sv_flight3Derr","SV 3D flight distance error",1, move+2);
  Draw("sv_flight2Derr","SV 2D flight distance error",1, move+1);
  Draw("sv_mass_3trk","SV mass if more than 3 tracks attached to an SV",1);
  }
if (Draw_muons_plots){  
  Draw("muon_multi"   ,      "number of muons", 1);       
  Draw("muon_multi_sel"   ,  "number of selected muons",1);
  Draw("mu_ptrel"     ,      "p_{T} rel. of the muon",1);    
  Draw("mu_chi2"      ,      "norm. #chi^{2} of the muon", 1);
  Draw("muon_Pt",	     "Muon p_{T}",1);	       
  Draw("muon_eta",	     "Muon #eta",1, move);	       
  Draw("muon_phi",	     "Muon #phi",1, move);	       
  Draw("muon_Ip3d",	     "Muon 3D IP",1, move+1);	       
  Draw("muon_Ip2d",	     "Muon 2D IP",1, move+1);	       
  Draw("muon_Sip3d",	     "Muon 3D IP significance",1, move+1);
  Draw("muon_Sip2d",	     "Muon 2D IP significance",1, move+1);
  Draw("muon_DeltaR",	     "Muon1 #Delta R",1);
  
  } 
if (Draw_discriminator_plots){  
  Draw("TCHE_extended1"       ,"TCHE (extended1)",1, move+1);    
  Draw("TCHP_extended1"       ,"TCHP (extended1)",1);    
  Draw("TCHE_extended2"       ,"TCHE (extended2)", 1);   
  Draw("TCHP_extended2"       ,"TCHP (extended2)",1);    
  Draw("discri_ssche0",      "SSVHE Discriminator", 1);
  Draw("discri_sschp0",      "SSVHP Discriminator", 1);

  Draw("TCHE"	      ,"TCHE Discriminator", 1); 		    
  Draw("TCHP"	      ,"TCHP Discriminator",1);  		    
  Draw("JP"	      ,"JP Discriminator",1);			    
  Draw("JBP"	      ,"JBP Discriminator",1);			    
  Draw("SSV"	      ,"SSVHE Discriminator",1);			    
  Draw("SSVHP"        ,"SSVHP Discriminator",1); 		    
  Draw("CSV"	      ,"CSV Discriminator",1, move+2);
  
  }
  			    
if (Draw_tagRate_plots){ 
  DrawTagRate("TCHE_extended1","TCHE (extended1)", 0);
  DrawTagRate("TCHP_extended1"," TCHP (extended1)", 0);
  DrawTagRate("TCHE_extended2","TCHE (extended2)", 0);
  DrawTagRate("TCHP_extended2","TCHP (extended2)", 0);
  DrawTagRate("discri_ssche0","SSVHE (extended)", 0);
  DrawTagRate("discri_sschp0","SSVHP (extended)", 0);

  DrawTagRate("TCHE"	      ,"TCHE Discriminator", 0);
  DrawTagRate("TCHP"	      ,"TCHP Discriminator", 0);
  DrawTagRate("JP"	      ,"JP Discriminator", 0);
  DrawTagRate("JBP"	      ,"JBP Discriminator", 0);
  DrawTagRate("SSV"	      ,"SSVHE Discriminator", 0);
  DrawTagRate("SSVHP"         ,"SSVHP Discriminator", 0);
  DrawTagRate("CSV"	      ,"CSV Discriminator", 0);

  
  }
  
if (Draw_2D_plots){
  Draw2DPlot("seltrack_vs_jetpt", "nr. of selected tracks as a function of the jet p_{T}", "jet p_{T}","nr. of selected tracks",0);
  Draw2DPlot("sv_mass_vs_flightDist3D", " SV mass as a function of the SV 3D flight distance ","SV 3D flight distance","SV mass",  0);
  Draw2DPlot("avg_sv_mass_vs_jetpt","Avg SV mass as a function of the jet p_{T}","jet p_{T}","Avg SV mass", 0);	       
  Draw2DPlot("sv_deltar_jet_vs_jetpt","#Delta R between the SV and the jet as a function of the jet p_{T}","jet p_{T}","#Delta R between the SV and the jet", 0); 
  Draw2DPlot("sv_deltar_sum_jet_vs_jetpt","#Delta R between the SV and the jet sum as a function of the jet p_{T}","jet p_{T}","#Delta R between the SV and the jet sum", 0);     
  Draw2DPlot("sv_deltar_sum_dir_vs_jetpt","#Delta R between the SV and the jet direction as a function of the jet p_{T}", "jet p_{T}","#Delta R between the SV and the jet direction", 0); 
  Draw2DPlot("muon_ptrel_vs_jetpt","Muon_p{T}^{rel} as a function of the jet p_{T}","jet p_{T}","Muon_p{T}^{rel}",	 0);       
  Draw2DPlot("muon_DeltaR_vs_jetpt","Muon #Delta R as a function of the jet p_{T}","jet p_{T}","Muon #Delta R",  0);	         
  }

}

//--------------------------

void Draw(TString name, TString histotitle, bool log, int move_legend)

{

 
 TH1F* hist_b;
 TH1F* hist_c;
 TH1F* hist_gsplit;
 TH1F* hist_l;
 TH1F* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_b         = (TH1F*)gROOT->FindObject(name+"_b");
 hist_c         = (TH1F*)gROOT->FindObject(name+"_c");
 hist_gsplit    = (TH1F*)gROOT->FindObject(name+"_bfromg");
 hist_l         = (TH1F*)gROOT->FindObject(name+"_l");
 hist_data      = (TH1F*)gROOT->FindObject(name+"_data");
 

 if (bOverflow) {
  OverFlowBinFix(hist_b);
  OverFlowBinFix(hist_c);
  OverFlowBinFix(hist_gsplit);
  OverFlowBinFix(hist_l);
  OverFlowBinFix(hist_data);
 }



 TH1F* histo_tot = (TH1F*) hist_b->Clone();
 histo_tot->Sumw2();
 histo_tot ->Add(hist_c);
 histo_tot ->Add(hist_gsplit);
 histo_tot ->Add(hist_l);  

 

 float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_gsplit->Integral() + hist_l->Integral());

 hist_b       ->Scale(scale_f);
 hist_c       ->Scale(scale_f);
 hist_gsplit  ->Scale(scale_f);
 hist_l       ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);
  
 double titleoffsety=0.2;
 double titlesizex=0.17;
 double titlesizey=0.2;
 double labelsizex=0.14;
 double labelsizey=0.12; 
  
   
 hist_data  ->GetYaxis()->SetLabelSize(labelsizey);
 hist_data  ->GetYaxis()->SetTitleSize(titlesizey);
 hist_data  ->GetYaxis()->SetTitleOffset(titleoffsety);
  
 hist_b     ->GetYaxis()->SetLabelSize(labelsizey);
 hist_b     ->GetYaxis()->SetTitleSize(titlesizey);
 hist_b     ->GetYaxis()->SetTitleOffset(titleoffsety);
  
 TH1F* histo_ratio;
 histo_ratio = (TH1F*) hist_data->Clone();
 histo_ratio->SetName("histo_ratio");
 histo_ratio->SetTitle("");
  
 histo_ratio->Divide(histo_tot);
  
 hist_data  ->SetLineWidth(2);
 hist_data  ->SetMarkerStyle(20);  
 hist_data  ->SetMarkerSize(0.75); 

 hist_c     ->SetFillColor(8);
 hist_b     ->SetFillColor(2);
 hist_gsplit->SetFillColor(7);
 hist_l     ->SetFillColor(4);
  
 histo_tot  ->SetLineColor(2);
  

 THStack *stack = new THStack("stack","stack");
  
 if (b_ordering){
  stack      ->Add(hist_b);
  stack      ->Add(hist_gsplit);
  stack      ->Add(hist_c);
  stack      ->Add(hist_l);
 }
 else {
  stack      ->Add(hist_l);
  stack      ->Add(hist_c);
  stack      ->Add(hist_gsplit);
  stack      ->Add(hist_b);
 }


 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);  
 //gStyle->SetLogy(log);
  
 TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
 c1->SetFillColor(10);
 c1->  cd();   
  

 TPad* canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
 canvas_1 ->Draw();
 canvas_1 ->cd();
 
 canvas_1->SetLogy(log);
 
 if (hist_data->GetMaximum() > stack->GetMaximum() ) stack->SetMaximum( hist_data->GetMaximum()*1.1) ;

 if (move_legend==2) {
    if (log)  stack->SetMaximum( 5.*stack->GetMaximum() ) ;
    else stack->SetMaximum( 1.1*stack->GetMaximum() ) ;
 }

 stack    ->Draw("hist");  
  
 stack    ->GetHistogram()->GetXaxis()->SetTitle(name);
 stack    ->GetHistogram()->GetYaxis()->SetTitle("entries");

 stack    ->GetHistogram()->SetTitleSize(0.08,"Y");
 stack    ->GetHistogram()->SetTitleOffset(0.65,"Y"); 

 hist_data->Draw("same e");

// TLegend* qw =  new TLegend(0.54,0.63,0.88,0.9);
 TLegend* qw;
 if (move_legend==1) {
   qw =  new TLegend(0.35,0.15,0.70,0.42);
 }
 else if (move_legend==3) {
  qw =  new TLegend(0.35,0.63,0.70,0.90);
 }
 else qw =  new TLegend(0.6,0.73,0.95,1.);
  
  //Legend
 qw->AddEntry(hist_data,     "e#mu ttbar data",                       "p");
 qw->AddEntry(hist_b,        "b quark"           ,         "f");
 qw->AddEntry(hist_gsplit,   "b from gluon splitting"     ,"f");
 qw->AddEntry(hist_c,        "c quark"           ,         "f");
 qw->AddEntry(hist_l,        "uds quark or gluon"     ,    "f");

 
 qw->SetFillColor(0);
 qw->Draw();
  
  
 TLatex *latex = new TLatex();
 latex->SetNDC();
 latex->SetTextSize(0.055);
 latex->SetTextFont(42); //22

 latex->SetTextAlign(13);
 latex->DrawLatex(0.08, 0.96, title);


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
  
 TString name_plot=name+"_Linear"+format; 
 if(log) name_plot=name+"_Log"+format;
 c1->SaveAs("Commissioning_plots/"+name_plot);


}

//--------------------------

void DrawTTbar(TString name, TString histotitle, bool log, int move_legend)

{

 
 TH1F* hist_ttbar;
 TH1F* hist_dy;
 TH1F* hist_st;
 TH1F* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_ttbar     = (TH1F*)gROOT->FindObject(name+"_ttbar");
 hist_dy        = (TH1F*)gROOT->FindObject(name+"_dy");
 hist_st        = (TH1F*)gROOT->FindObject(name+"_st");
 hist_data      = (TH1F*)gROOT->FindObject(name+"_data");
 

 if (bOverflow) {
  OverFlowBinFix(hist_ttbar);
  OverFlowBinFix(hist_dy);
  OverFlowBinFix(hist_st);
  OverFlowBinFix(hist_data);
 }



 TH1F* histo_tot = (TH1F*) hist_ttbar->Clone();
 histo_tot->Sumw2();
 histo_tot ->Add(hist_dy);
 histo_tot ->Add(hist_st);

 

/*
 float scale_f = (hist_data->Integral())/(hist_tot->Integral());

 hist_ttbar   ->Scale(scale_f);
 hist_dy      ->Scale(scale_f);
 hist_st      ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);
*/
  
 double titleoffsety=0.2;
 double titlesizex=0.17;
 double titlesizey=0.2;
 double labelsizex=0.14;
 double labelsizey=0.12; 
  
   
 hist_data  ->GetYaxis()->SetLabelSize(labelsizey);
 hist_data  ->GetYaxis()->SetTitleSize(titlesizey);
 hist_data  ->GetYaxis()->SetTitleOffset(titleoffsety);
  
 hist_ttbar     ->GetYaxis()->SetLabelSize(labelsizey);
 hist_ttbar     ->GetYaxis()->SetTitleSize(titlesizey);
 hist_ttbar     ->GetYaxis()->SetTitleOffset(titleoffsety);
  
 TH1F* histo_ratio;
 histo_ratio = (TH1F*) hist_data->Clone();
 histo_ratio->SetName("histo_ratio");
 histo_ratio->SetTitle("");
  
 histo_ratio->Divide(histo_tot);
  
 hist_data  ->SetLineWidth(2);
 hist_data  ->SetMarkerStyle(20);  
 hist_data  ->SetMarkerSize(0.75); 

 hist_dy     ->SetFillColor(kAzure-2);
 hist_ttbar  ->SetFillColor(kRed+1);
 hist_st     ->SetFillColor(kMagenta);
  
  

 THStack *stack = new THStack("stack","stack");
  
 stack      ->Add(hist_ttbar);
 stack      ->Add(hist_st);
 stack      ->Add(hist_dy);


 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);  
 //gStyle->SetLogy(log);
  
 TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
 c1->SetFillColor(10);
 c1->  cd();   
  

 TPad* canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
 canvas_1 ->Draw();
 canvas_1 ->cd();
 
 canvas_1->SetLogy(log);
 
 if (hist_data->GetMaximum() > stack->GetMaximum() ) stack->SetMaximum( hist_data->GetMaximum()*1.1) ;

 if (move_legend==2) {
    if (log)  stack->SetMaximum( 5.*stack->GetMaximum() ) ;
    else stack->SetMaximum( 1.1*stack->GetMaximum() ) ;
 }

 stack    ->Draw("hist");  
  
 stack    ->GetHistogram()->GetXaxis()->SetTitle(name);
 stack    ->GetHistogram()->GetYaxis()->SetTitle("entries");

 stack    ->GetHistogram()->SetTitleSize(0.08,"Y");
 stack    ->GetHistogram()->SetTitleOffset(0.65,"Y"); 

 hist_data->Draw("same e");

// TLegend* qw =  new TLegend(0.54,0.63,0.88,0.9);
 TLegend* qw;
 if (move_legend==1) {
   qw =  new TLegend(0.35,0.15,0.70,0.42);
 }
 else if (move_legend==3) {
  qw =  new TLegend(0.35,0.63,0.70,0.90);
 }
 else qw =  new TLegend(0.6,0.73,0.95,1.);
  
  //Legend
 qw->AddEntry(hist_data,     "e#mu data",                       "p");
 qw->AddEntry(hist_ttbar,    "t#bar{t} MC"           ,         "f");
 qw->AddEntry(hist_dy,       "DY"     ,"f");
 qw->AddEntry(hist_st,       "tW"           ,         "f");

 
 qw->SetFillColor(0);
 qw->Draw();
  
  
 TLatex *latex = new TLatex();
 latex->SetNDC();
 latex->SetTextSize(0.055);
 latex->SetTextFont(42); //22

 latex->SetTextAlign(13);
 latex->DrawLatex(0.08, 0.96, title);


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
  
 TString name_plot="ttbar_"+name+"_Linear"+format; 
 if(log) name_plot="ttbar_"+name+"_Log"+format;
 c1->SaveAs("Commissioning_plots/"+name_plot);


}



//--------------------------


void DrawTagRate(TString name, TString histotitle, bool log){


 TH1F* hist_b;
 TH1F* hist_c;
 TH1F* hist_gsplit;
 TH1F* hist_l;
 TH1F* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_b         = (TH1F*)gROOT->FindObject(name+"_b");
 hist_c         = (TH1F*)gROOT->FindObject(name+"_c");
 hist_gsplit    = (TH1F*)gROOT->FindObject(name+"_bfromg");
 hist_l         = (TH1F*)gROOT->FindObject(name+"_l");
 hist_data      = (TH1F*)gROOT->FindObject(name+"_data");
 
 if (bOverflow) {
  OverFlowBinFix(hist_b);
  OverFlowBinFix(hist_c);
  OverFlowBinFix(hist_gsplit);
  OverFlowBinFix(hist_l);
  OverFlowBinFix(hist_data);
 }


 TH1F* histo_tot = (TH1F*) hist_b->Clone();
 histo_tot->Sumw2();
 histo_tot ->Add(hist_c);
 histo_tot ->Add(hist_gsplit);
 histo_tot ->Add(hist_l);  
  
 float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_gsplit->Integral() + hist_l->Integral());

 hist_b       ->Scale(scale_f);
 hist_c       ->Scale(scale_f);
 hist_gsplit  ->Scale(scale_f);
 hist_l       ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);
 
 int   nbinx=hist_data->GetNbinsX();
 float minx=hist_data->GetXaxis()->GetXmin();
 float maxx=hist_data->GetXaxis()->GetXmax();
 
 TH1F * TagRate_Data = new TH1F ("TagRate_Data",histotitle,nbinx,minx,maxx);
 TH1F * TagRate_MC   = new TH1F ("TagRate_MC",histotitle,nbinx,minx,maxx);
 TH1F * TagRate_MC_b = new TH1F ("TagRate_MC_b",histotitle,nbinx,minx,maxx);
 TH1F * TagRate_MC_c = new TH1F ("TagRate_MC_c",histotitle,nbinx,minx,maxx);
 TH1F * TagRate_MC_udsg = new TH1F ("TagRate_MC_udsg",histotitle,nbinx,minx,maxx);
 TH1F * TagRate_MC_gspl = new TH1F ("TagRate_MC_gspl",histotitle,nbinx,minx,maxx);

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
   double titlesizey=0.2;
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
  
  TH1F* histo_ratio;
  histo_ratio = (TH1F*) TagRate_Data->Clone();
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

  if (b_ordering){
   hs->Add(TagRate_MC_b);
   hs->Add(TagRate_MC_gspl);  
   hs->Add(TagRate_MC_c);
   hs->Add(TagRate_MC_udsg);
  }
  else {
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

  if (TagRate_Data->GetMaximum() > hs->GetMaximum() ) {
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

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.055);
  latex->SetTextFont(42); //22

  latex->SetTextAlign(13);
 latex->DrawLatex(0.08, 0.96, title);
  

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
  
  TString name_plot="tag_"+name+"_Linear"+format;
  if(log) name_plot="tag_"+name+"_Log"+format;
  c1->SaveAs("Commissioning_plots/"+name_plot);
}


//--------------------------


void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log){



 TH2F* hist_b;
 TH2F* hist_c;
 TH2F* hist_gsplit;
 TH2F* hist_l;
 TH2F* hist_data;
 
 
 TFile *myFile     = new TFile(filename);
 
 myFile->cd();
 hist_b         = (TH2F*)gROOT->FindObject(name+"_b");
 hist_c         = (TH2F*)gROOT->FindObject(name+"_c");
 hist_gsplit    = (TH2F*)gROOT->FindObject(name+"_bfromg");
 hist_l         = (TH2F*)gROOT->FindObject(name+"_l");
 hist_data      = (TH2F*)gROOT->FindObject(name+"_data");
 

 TH2F* histo_tot = (TH2F*) hist_b->Clone();
 histo_tot ->Add(hist_c);
 histo_tot ->Add(hist_gsplit);
 histo_tot ->Add(hist_l); 
 
 float scale_f = (hist_data->Integral())/(hist_b->Integral() + hist_c ->Integral()+ hist_gsplit->Integral() + hist_l->Integral());

 hist_b       ->Scale(scale_f);
 hist_c       ->Scale(scale_f);
 hist_gsplit  ->Scale(scale_f);
 hist_l       ->Scale(scale_f);
 histo_tot    ->Scale(scale_f);

  
 TProfile * pro_mc = histo_tot->ProfileX();
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

  TString name_plot=name+"_Linear"+format;
  if(log) name_plot=name+"_Log"+format;
  canvas->SaveAs("Commissioning_plots/"+name_plot);

}

//------------

void OverFlowBinFix(TH1F* histo){

  Float_t val, errval;

  val=histo->GetBinContent(1)+histo->GetBinContent(0);
  errval=0;
  if(histo->GetBinContent(1)!=0.)
    errval+=pow(histo->GetBinError(1)/histo->GetBinContent(1),2);
  if(histo->GetBinContent(0)!=0.)
    errval+=pow(histo->GetBinError(0)/histo->GetBinContent(0),2);
  errval=sqrt(errval)*val;
  histo->SetBinContent(1,val);
  histo->SetBinError(1,errval);
  histo->SetBinContent(0,0);
  histo->SetBinError(0,0);

  Int_t highbin=histo->GetNbinsX();

  val=histo->GetBinContent(highbin)+histo->GetBinContent(highbin+1);
  errval=0;
  if(histo->GetBinContent(highbin)!=0.)
    errval+=pow(histo->GetBinError(highbin)/histo->GetBinContent(highbin),2);
  if(histo->GetBinContent(highbin+1)!=0.)
    errval+=pow(histo->GetBinError(highbin+1)/histo->GetBinContent(highbin+1),2);
  errval=sqrt(errval)*val;
  histo->SetBinContent(highbin,val);
  histo->SetBinError(highbin,errval);
  histo->SetBinContent(highbin+1,0);
  histo->SetBinError(highbin+1,0);
}



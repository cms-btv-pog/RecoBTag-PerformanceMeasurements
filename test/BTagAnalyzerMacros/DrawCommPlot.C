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
#include <TString.h>

using namespace std;

//TString filename="qcd_pythia/jet_40_60_500_bugfix/output_all.root";
//TString filename="qcd_pythia/jet_40_60_500/output_newalgo1.root";
TString filename="qcd_pythia/jet_40_60_500_bugfix/output_newpf1.root";
//TString filename="qcd_herwig/jet40/output_herwig_jet40.root";
TString datacaption = "HLT_PFJet40, jet pT>60 ";

bool btag = false;
//TString filename="qcd_pythia/btag_40_65_300_bugfix/output_btag_all.root";
//TString filename="qcd_pythia/btag_40_65_300/output_btag_newalgo1.root";
//TString filename="qcd_pythia/btag_40_binning2/output_btag_all.root";
//TString filename="qcd_pythia/btag_40_binning2/output_btag_mcinclu.root";
//TString datacaption = "HLT_BTag_Mu5_Dijet40, jet pT>65 ";

//TString filename="qcd_pythia/btag_20_45_300/output_btag_all.root";

//TString filename="qcd_herwig/btag20/output_btag20_herwig.root";
//TString filename="qcd_pythia/btag_20_45_300/output_btag_newalgo1.root";
bool test_bin= false;
//TString filename="qcd_pythia/btag_40_65_300/output_btag_oldprocess.root";

//TString filename="test_proof/save_proof_btag40_65_300_w_trigMC.root";
//TString filename="qcd_pythia/btag_40_65_300/output_btag_newalgo.root";
//TString filename="qcd_pythia/btag_40_65_300/output_btag_oldprocess.root";
//TString filename="test_proof/save_proof_trig40_60_500.root";
//TString filename="test/abcd_new/output_all.root";
//TString filename="test/abcd_80_80_470/output_all.root";
//TString filename="test/compaNt/final_out.root";
//TString filename="test/camille_abcd/test_cami_80_100_500.root";
//TString dir4plots="test/Commissioning_plots_pthat";
//TString dir4plots="test_proof/Commissioning_btag40_65";
TString dir4plots="Commissioning_plots";
TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  19.7 fb^{-1}";
//TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  old Process";
//TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  14. fb^{-1}";
//TString title= "OLD #sqrt{s} = 8 TeV";
//TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  A:0.89 fb^{-1}";
//TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  B:4.43 fb^{-1}";
//TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  C:7.15 fb^{-1}";
//TString title= "CMS 2013 preliminary, #sqrt{s} = 8 TeV,  D:7.27 fb^{-1}";
//TString datacaption = "HLT_BTag_Mu5_Dijet20, jet pT>45 ";
//TString datacaption = "HLT_PFJet80, jet pT>100 ";
//TString datacaption = "HLT_PFJet80, jet pT>80 ";
TString format=".pdf"; // .png or .pdf or .gif
bool bOverflow=true;
bool web = false;



void Draw(TString name, TString histotitle, bool log);
void DrawMC(TString name, TString histotitle, bool log);
void DrawTagRate(TString name, TString histotitle, bool log);
void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log);
void OverFlowBinFix(TH1D* );
void anplot();
void anplot2();

//--------

void DrawCommPlot(bool Draw_track_plots, bool Draw_Nminus1_plots, bool Draw_sv_plots, bool Draw_muons_plots, bool Draw_discriminator_plots , 
bool Draw_newdiscriminator_plots, bool Draw_tagRate_plots, bool Draw_2D_plots, bool Draw_pflepton){

  TString action = "mkdir "+dir4plots+"/";
  system(action);
  
  
//Draw("jet_multi"    ,"number of jets",1);	    
Draw("jet_pt_all"   ,"pT of all jets",1);	    
Draw("jet_eta"      ,"eta of all jets",	0);    
Draw("jet_phi"      ,"phi of all jets",0);
DrawMC("nPV"      ,"# of PV",0);

if (Draw_track_plots){	      
  Draw("track_multi"  ,      "number of tracks in the jets",0);		   
  Draw("trk_multi_sel"  ,    "number of selected tracks in the jets",0);	   
  Draw("track_chi2"   ,      "normalized #chi^{2} of the tracks"	,1);	   
  Draw("track_nHit" ,      "number of hits",1);		   
  Draw("track_HPix"   ,      "number of hits in the Pixel",1);		   
  Draw("track_len"     ,     "Track decay length",1);				   
  Draw("track_dist"    ,     "Track distance to the jet axis"	,1);	   
  Draw("track_dz"     ,      "Track IP_dz",1);				   
  Draw("track_pt"     ,      "pT of all the tracks",1);			   
  Draw("track_pt15"     ,      "pT of all the tracks",1);			   
  //Draw("track_isfromSV",     "Track is from SV",1);			   

  Draw("track_IPs"    ,      "3D IP significance of all tracks",1);	   
  Draw("track_IPs1tr" ,      "3D IP significance of the first track",1);	   
  Draw("track_IPs2tr" ,      "3D IP significance of the second track",1);     
  Draw("track_IPs3tr" ,      "3D IP significance of the third track",1);	
  Draw("track_IP"     ,      "3D IP of all tracks",1);			   
  Draw("track_IP1tr"  ,      "3D IP of the first track"	,1);	   
  Draw("track_IP2tr"  ,      "3D IP of the second track",1);		   
  Draw("track_IP3tr"  ,      "3D IP of the third track"	,1);	
  Draw("track_IP2Ds"	,    "2D IP significance of all tracks"	,1);   
  Draw("track_IP2Ds1tr" ,    "2D IP significance of the first track",1);	   
  Draw("track_IP2Ds2tr" ,    "2D IP significance of the second track" ,1);    
  Draw("track_IP2Ds3tr" ,    "2D IP significance of the third track" ,1);
  Draw("track_IP2D"    ,     "2D IP of all tracks",1);			   
  Draw("track_IP2D1tr" ,     "2D IP of the first track"	,1);	   
  Draw("track_IP2D2tr" ,     "2D IP of the second track",1);		   
  Draw("track_IP2D3tr" ,     "2D IP of the third track"	,1);	
  Draw("track_IP2Derr" ,     "2D IP error of all tracks",1);				   
  Draw("track_IP2Derr1tr" ,  "2D IP error of the first track",1);  	   
  Draw("track_IP2Derr2tr" ,  "2D IP error of the second track" ,1);	   
  Draw("track_IP2Derr3tr" ,  "2D IP error of the third track",1);  	
  Draw("track_IPerr"   ,     "3D IP error of all tracks",1);				   
  Draw("track_IPerr1tr"   ,  "3D IP error of the first track" ,1); 	   
  Draw("track_IPerr2tr"   ,  "3D IP error of the second track" ,1);	   
  Draw("track_IPerr3tr"   ,  "3D IP error of the third track" ,1); 	
	
  }  
if (Draw_Nminus1_plots){  
  Draw("track_chi2_cut"    ,"Normalized #chi^{2} @N-1 step",1);	    
  Draw("track_nHit_cut"  ," Number of hits @N-1 step",1);
  Draw("track_HPix_cut"    ,"Number of hits in the Pixel @N-1 step",1);  
  Draw("track_len_cut"     ,"Decay length @N-1 step",1);		    
  Draw("track_dist_cut"    ,"Distance to the jet axis @N-1 step" ,1);  
  Draw("track_dz_cut"     , "Transverse IP @N-1 step",1);		    
  Draw("track_pt_cut"	   ,"Track pT @N-1 step",1);
  Draw("track_pt15_cut"	   ,"Track pT @N-1 step",1);
  }
if (Draw_sv_plots){

  Draw("jet_pt_sv"    ,"pT of jets containing a SV",1);
  Draw("sv_multi_0","nr. of SV including bin 0",1); 
  Draw("sv_multi","nr. of SV",1);
  Draw("sv_mass","SV mass",0);
  Draw("sv_mass_3trk","SV mass if #tracks@SV >=3",0);
  Draw("sv_chi2norm","SV norm. #chi^{2}",1);
  Draw("sv_deltaR_jet","Delta R between the jet and the SV direction.",0);
  Draw("sv_deltaR_sumJet","#Delta R between the jet and the SV {p}",0);
  Draw("sv_deltaR_sumDir","#Delta R between the SV direction and the SV {p}",0);
  Draw("sv_en_ratio","SV energy ratio",0);	
  Draw("sv_aboveC","IP2D of the first track above the charm threshold",1);	
  Draw("sv_pt","SV p_{T}",1); 	
  Draw("sv_eta","SV #eta",0);	
  Draw("sv_phi","SV #phi",0);	
  Draw("sv_flight3D","SV 3D flight distance",1);	
  Draw("sv_flight2D","SV 2D flight distance",1);	
  Draw("sv_flight2D_3trk","SV 2D flight dist. if #tracks@SV >=3",1);	
  Draw("sv_flight2DSig_3trk","SV 2D flight dist. sig. if #tracks@SV >=3",1);	
  Draw("sv_flight3DSig","SV 3D flight distance significance",1);
  Draw("sv_flightSig2D","SV 2D flight distance significance",1);
  Draw("sv_flight3Derr","SV 3D flight distance error",1);
  Draw("sv_flight2Derr","SV 2D flight distance error",1);
  Draw("svnTrk","nr. of tracks from a SV",1);	
  Draw("svnTrk_firstVxt","nr. of tracks from the first SV",1);
  }
if (Draw_muons_plots){  
  Draw("muon_multi"   ,      "number of muons", 1);       
  Draw("muon_multi_sel"   ,  "number of selected muons",1);
  Draw("mu_ptrel"     ,      "p_{T} rel. of the muon",0);    
  Draw("mu_chi2"      ,      "norm. #chi^{2} of the muon", 1);
  Draw("muon_Pt",	     "Muon p_{T}",1);	       
  Draw("muon_eta",	     "Muon #eta",0);	       
  Draw("muon_phi",	     "Muon #phi",0);	       
  Draw("muon_Ip3d",	     "Muon 3D IP",1);	       
  Draw("muon_Ip2d",	     "Muon 2D IP",1);	       
  Draw("muon_Sip3d",	     "Muon 3D IP significance",1);
  Draw("muon_Sip2d",	     "Muon 2D IP significance",1);
  Draw("muon_DeltaR",	     "Muon1 #Delta R",0);
  
  } 
if (Draw_discriminator_plots){  
  Draw("TCHE_extended1"       ,"TCHE (extended)",1);    
  Draw("TCHP_extended1"       ,"TCHP (extended)",1);    
  Draw("discri_ssche0",      "SSVHE Discriminator", 1);
  Draw("discri_sschp0",      "SSVHP Discriminator", 1);

  Draw("TCHE"	      ,"TCHE Discriminator", 1); 		    
  Draw("TCHP"	      ,"TCHP Discriminator",1);  		    
  Draw("JP"	      ,"JP Discriminator",1);			    
  Draw("JBP"	      ,"JBP Discriminator",1);			    
  Draw("SSV"	      ,"SSVHE Discriminator",1);			    
  Draw("SSVHP"        ,"SSVHP Discriminator",1); 		    
  Draw("CSV"	      ,"CSV Discriminator",1);
  
  }
if (Draw_newdiscriminator_plots) {
  Draw("RetCombSvx"         ,"CSV_v1 Discriminator",1);
  Draw("RetCombSvx"         ,"CSV_v1 Discriminator",0);
  Draw("CombCSVJP"          ,"CSVJP Discriminator",1);
  Draw("CombCSVSL"          ,"CSVSL_v1 Discriminator",1);
  Draw("CombCSVJPSL"        ,"CSVJPSL Discriminator",1);
  Draw("SoftMu"        ,"Soft Muon Discriminator",1);
  Draw("SoftMu"        ,"Soft Muon Discriminator",0);
  Draw("SoftEl"        ,"Soft Electron Discriminator",1);
}
  			    
if (Draw_tagRate_plots){ 
  DrawTagRate("TCHE_extended1","TCHE (extended)", 1);
  DrawTagRate("TCHP_extended1"," TCHP (extended)", 1);
  DrawTagRate("discri_ssche0","SSVHE (extended)", 1);
  DrawTagRate("discri_sschp0","SSVHP (extended)", 1);

  DrawTagRate("TCHE"	      ,"TCHE Discriminator", 1);
  DrawTagRate("TCHP"	      ,"TCHP Discriminator", 1);
  DrawTagRate("JP"	      ,"JP Discriminator", 1);
  DrawTagRate("JBP"	      ,"JBP Discriminator", 1);
  DrawTagRate("SSV"	      ,"SSVHE Discriminator", 1);
  DrawTagRate("SSVHP"         ,"SSVHP Discriminator", 1);
  DrawTagRate("CSV"	      ,"CSV Discriminator", 1);

  
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

if (Draw_pflepton) {

   Draw("pfelectron_multi","number of pfelectron",1);
   Draw("pfelectron_pt","pfelectron p_{T}",1);
   Draw("pfelectron_eta","pfelectron #eta",0);
   Draw("pfelectron_phi","pfelectron #phi", 0);
   Draw("pfelectron_Sip","3D IP significance of pfelectron",1);
   Draw("pfelectron_ptrel","pT rel. of the pfelectron",0);
   Draw("pfelectron_ratio","ratio of pfelectron", 1);
   Draw("pfelectron_ratiorel",   "ratioRel of pfelectron",1);
   Draw("pfelectron_deltar", "#DeltaR(pfelectron,jet)",0);

   Draw("pfmuon_multi","number of pfmuon",1);
   Draw("pfmuon_goodquality","quality of pfmuons",0);
   Draw("pfmuon_pt","pfmuon p_{T}",1);
   Draw("pfmuon_eta","pfmuon #eta",0);
   Draw("pfmuon_phi","pfmuon #phi", 0);
   Draw("pfmuon_Sip","3D IP significance of pfmuon",1);
   Draw("pfmuon_ptrel","pT rel. of the pfmuon",0);
   Draw("pfmuon_ratio","ratio of pfmuon", 1);
   Draw("pfmuon_ratiorel",   "ratioRel of pfmuon",1);
   Draw("pfmuon_deltar", "#DeltaR(pfmuon,jet)",0);

}

if (test_bin) {
  for (int i=1;i<10;i++) {

    char bbb[10];
    sprintf(bbb,"bin%d",i);
    TString bini = bbb;
    cout <<" i " << i << "bin" << i << " "  << bini << endl;
    Draw(bini+"_SSVHP",              "SSVHP ("+bini+")",                      1);
    Draw(bini+"_JP",                 "JP ("+bini+")",                         1);
    Draw(bini+"_CSV",                "CSV ("+bini+")",                        1);
    Draw(bini+"_SVmass",             "SV_mass ("+bini+")",                    0);
    Draw(bini+"_mu_ptrel",           "pT rel. of the muon ("+bini+")",        0);
    Draw(bini+"_muon_DeltaR",        "DeltaR(jet,mu) ("+bini+")",             0);
    Draw(bini+"_sv_flight3DSig",     "3D flight dist. sig ("+bini+")",        1);
    Draw(bini+"_sv_pt",              "Vtx p_{T} ("+bini+")",                  0);
    Draw(bini+"_sv_deltaR_jet",      "DeltaR(jet,SV) ("+bini+")",             0);
   }


}
}
void anplot() {

  TString action = "mkdir "+dir4plots+"/";
  system(action);
  
  
Draw("jet_pt_all"   ,"pT of all jets",1);	    

  Draw("track_IPs"    ,      "3D IP significance of all tracks",1);	   
  Draw("track_IP"     ,      "3D IP of all tracks",1);			   
  Draw("track_IPerr"   ,     "3D IP error of all tracks",1);				   
  Draw("track_chi2_cut"    ,"Normalized #chi^{2} @N-1 step",1);	    
  Draw("track_nHit_cut"  ," Number of hits @N-1 step",1);
  Draw("track_HPix_cut"    ,"Number of hits in the Pixel @N-1 step",1);  
  Draw("track_len_cut"     ,"Decay length @N-1 step",1);		    
  Draw("track_dist_cut"    ,"Distance to the jet axis @N-1 step" ,1);  
  Draw("track_dz_cut"     , "IP_dz @N-1 step",1);		    
  Draw("track_pt_cut"	   ,"Track pT @N-1 step",1);
  Draw("track_pt15_cut"	   ,"Track pT @N-1 step",1);
  Draw("track_IP2D_cut"    ,"IP2D @N-1 step",1);

/*
  Draw("sv_multi_0","nr. of SV including bin 0",1); 
  Draw("sv_mass","SV mass",0);
  Draw("sv_deltaR_jet","Delta R between the jet and the SV direction.",0);
  Draw("sv_pt","SV p_{T}",1); 	

  Draw("mu_ptrel"     ,      "p_{T} rel. of the muon",0);    
  Draw("muon_Pt",	     "Muon p_{T}",1);	       
  Draw("muon_Sip3d",	     "Muon 3D IP significance",1);
  Draw("muon_DeltaR",	     "Muon1 #Delta R",0);
  
  Draw("discri_ssche0",      "SSVHE Discriminator", 1);
  Draw("discri_sschp0",      "SSVHP Discriminator", 1);

  Draw("TCHE"	      ,"TCHE Discriminator", 1); 		    
  Draw("TCHP"	      ,"TCHP Discriminator",1);  		    
  Draw("JP"	      ,"JP Discriminator",1);			    
  Draw("JBP"	      ,"JBP Discriminator",1);			    
  Draw("CSV"	      ,"CSV Discriminator",1);
*/
}
void anplot2() {

  TString action = "mkdir "+dir4plots+"/";
  system(action);
  
  
  Draw("CSV"	      ,"CSV Discriminator",1);
  Draw("CombCSVSL"          ,"CSVSL_v1 Discriminator",1);
/*
  Draw("RetCombSvx"         ,"CSV_v1 Discriminator",1);
  Draw("RetCombSvx"         ,"CSV_v1 Discriminator",0);
  Draw("CombCSVJP"          ,"CSVJP Discriminator",1);
  Draw("CombCSVJPSL"        ,"CSVJPSL Discriminator",1);
  Draw("SoftMu"        ,"Soft Muon Discriminator",1);
  Draw("SoftEl"        ,"Soft Electron Discriminator",1);
*/
  			    

}

//--------------------------

void Draw(TString name, TString histotitle, bool log)

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
 
 bool ovbin = true;
  for (int i=1;i<10;i++) {
    char ccc[10];
    sprintf(ccc,"bin%d",i);
    TString biccc = ccc;
    if (name==biccc+"_SVmass" ||
        name==biccc+"_sv_flight3DSig" ||
        name==biccc+"_sv_pt" ||
        name==biccc+"_sv_deltaR_jet" ) ovbin= false;
  }

 if (bOverflow && name!="SSV" && name!="SSVHP" && name!="sv_mass" 
     && name!="sv_mass_3trk" && ovbin) {
  OverFlowBinFix(hist_b);
  OverFlowBinFix(hist_c);
  OverFlowBinFix(hist_gsplit);
  OverFlowBinFix(hist_l);
  OverFlowBinFix(hist_data);
 }



 TH1D* histo_tot = (TH1D*) hist_b->Clone();
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
 hist_gsplit->SetFillColor(7);
 hist_l     ->SetFillColor(4);
  
 histo_tot  ->SetLineColor(2);
  

 THStack *stack = new THStack("stack","stack");
  
 stack      ->Add(hist_b);
 stack      ->Add(hist_gsplit);
 stack      ->Add(hist_c);
 stack      ->Add(hist_l);

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
 if ((name=="CombCSVJP" || name=="CombCSVJPSL") && btag) stack->SetMaximum( hist_data->GetMaximum()*3.);
 if (name=="jet_phi" || name=="sv_phi" || name=="muon_phi") {
  if (log) stack->SetMinimum(0.01);
  else  stack->SetMinimum(0.); 
 }
 if (test_bin) {
 float minmin=hist_b->GetMinimum();
 if (minmin<hist_l->GetMinimum()) minmin=hist_l->GetMinimum();
 if (minmin<hist_c->GetMinimum()) minmin=hist_c->GetMinimum();
 if (minmin<hist_gsplit->GetMinimum()) minmin=hist_gsplit->GetMinimum();
 if (minmin>0) stack->SetMinimum(minmin/2);
 else {
  if (log) stack->SetMinimum(0.01);
  else  stack->SetMinimum(0.); 
 }
 }


 if (name=="CSV" || name=="CombCSVSL") {
// stack->SetMaximum(2*stack->GetMaximum());
 if (stack->GetMinimum()> hist_b->GetMinimum()) {
 stack->SetMinimum(hist_b->GetMinimum());
 }}

 if (name=="pfelectron_eta" || name=="pfmuon_eta") stack->SetMaximum( hist_data->GetMaximum()*1.5);

 stack    ->Draw("hist");  
  
 stack    ->GetHistogram()->GetXaxis()->SetTitle(name);
 stack    ->GetHistogram()->GetYaxis()->SetTitle("entries");

 stack    ->GetHistogram()->SetTitleSize(0.08,"Y");
 stack    ->GetHistogram()->SetTitleOffset(0.65,"Y"); 

 hist_data->Draw("same e");


 int move_legend=0;
 if (name=="jet_phi" || name=="sv_phi" || name=="muon_phi" || name=="pfelectron_phi" || name=="pfmuon_phi") move_legend=1;
 if (btag) {
  if (name=="CombCSVJPSL" || name=="CombCSVSL" || (name=="SoftMu" && !log) ) move_legend=2;
  if (name=="CSV" || name=="CombCSVJP" || (name=="RetCombSvx" && !log) ) move_legend=3;
  if ((name=="RetCombSvx" || name=="SoftMu" ) && log) move_legend=1;
 }
 if (log && name=="sv_en_ratio" ) move_legend=1;
// TLegend* qw =  new TLegend(0.54,0.63,0.88,0.9);
 TLegend* qw;
 if (move_legend==1) {
   qw =  new TLegend(0.35,0.25,0.70,0.52);
 }
 else if (move_legend==2) {
     qw =  new TLegend(0.10,0.63,0.45,0.90);
 }
 else if (move_legend==3) {
     qw =  new TLegend(0.35,0.63,0.70,0.90);
 }
 else qw =  new TLegend(0.6,0.73,0.95,1.);

  
  //Legend
 qw->AddEntry(hist_data,     datacaption,                       "p");
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
 c1->SaveAs(dir4plots+"/"+name_plot);

 if (log && web) {  // save also _Linear for web
  canvas_1 ->cd();
  canvas_1->SetLogy(false);
  c1->cd();
  c1->SaveAs(dir4plots+"/"+name+"_Linear"+format);
 }




}


//--------------------------


void DrawTagRate(TString name, TString histotitle, bool log){


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
 

 TH1D* histo_tot = (TH1D*) hist_b->Clone();
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

  hs->Add(TagRate_MC_b);
  hs->Add(TagRate_MC_gspl);  
  hs->Add(TagRate_MC_c);
  hs->Add(TagRate_MC_udsg);
  
  
  
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
  c1_1->SetLogy(log);

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
  qw->AddEntry(TagRate_Data,        datacaption                     ,"p");
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
  c1->SaveAs(dir4plots+"/"+name_plot);

  if (log && web) {  // save also _Linear for web
   c1_1 ->cd();
   c1_1->SetLogy(false);
   c1->cd();
   c1->SaveAs(dir4plots+"/tag_"+name+"_Linear"+format);
  }
}


//--------------------------


void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log){



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

  float minhist= pro_mc_gspl->GetMinimum();
  if (pro_mc_b->GetMinimum() < minhist) minhist = pro_mc_b->GetMinimum()*0.9;
  if (pro_mc_c->GetMinimum() < minhist) minhist = pro_mc_c->GetMinimum()*0.9;
  if (pro_mc_udsg->GetMinimum() < minhist) minhist = pro_mc_udsg->GetMinimum()*0.9;
  if (pro_mc->GetMinimum() < minhist) minhist = pro_mc->GetMinimum()*0.9;
  if (pro_data->GetMinimum() < minhist) minhist = pro_data->GetMinimum()*0.9;

  if (maxhist> pro_mc_gspl->GetMaximum()) pro_mc_gspl->SetMaximum(maxhist);
  if (pro_mc_gspl->GetMinimum() >minhist) pro_mc_gspl->SetMinimum(minhist);

  pro_mc_gspl->Draw("hist");
  pro_mc_b->Draw("hist,same");
  pro_mc_c->Draw("hist,same");
  pro_mc_udsg->Draw("hist,same");
  pro_mc->Draw("hist,same");
  pro_data->Draw("e,same");

  TLegend* qw = 0;
  qw =  new TLegend(0.6,0.73,0.95,1.);
  qw->AddEntry(pro_data,        datacaption                   ,"p");
  qw->AddEntry(pro_mc,          "total "                 ,"l");
  qw->AddEntry(pro_mc_b,        "b quark"                ,"l");
  qw->AddEntry(pro_mc_gspl,     "b from gluon splitting" ,"l");
  qw->AddEntry(pro_mc_c,        "c quark"                ,"l");
  qw->AddEntry(pro_mc_udsg,     "uds quark or gluon"     ,"l");

  qw->SetFillColor(0);
  qw->Draw();


  TString name_plot=name+"_Linear"+format;
  if(log) name_plot=name+"_Log"+format;
  canvas->SaveAs(dir4plots+"/"+name_plot);

}

//------------

void OverFlowBinFix(TH1D* histo){

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


void DrawMC(TString name, TString histotitle, bool log)
 
{
 TFile *myFile     = new TFile(filename);
 TH1D* hist_mc;     
 TH1D* hist_data;  
 if (name=="nPV") {
//   hist_mc       = (TH1D*)gROOT->FindObject(name+"_mc_unw");
   hist_mc       = (TH1D*)gROOT->FindObject(name+"_mc");
  }
 else {
   hist_mc       = (TH1D*)gROOT->FindObject(name+"_mc");
 }
 hist_data     = (TH1D*)gROOT->FindObject(name+"_data");
 float scale_f = (hist_data->Integral())/(hist_mc->Integral());
 hist_mc       ->Scale(scale_f);

 double titleoffsety=0.2;
 double titlesizex=0.17;
 double titlesizey=0.2;
 double labelsizex=0.14;
  
   
 hist_data  ->GetYaxis()->SetTitleSize(titlesizey);
 hist_data  ->GetYaxis()->SetTitleOffset(titleoffsety);
  
 hist_mc     ->GetYaxis()->SetTitleSize(titlesizey);
 hist_mc     ->GetYaxis()->SetTitleOffset(titleoffsety);
  
 TH1D* histo_ratio;
 histo_ratio = (TH1D*) hist_data->Clone();
 histo_ratio->Sumw2();
 histo_ratio->SetName("histo_ratio");
 histo_ratio->SetTitle("");
  
 histo_ratio->Divide(hist_mc);
  
 hist_data  ->SetLineWidth(2);
 hist_data  ->SetMarkerStyle(20);  
 hist_data  ->SetMarkerSize(0.75); 

 hist_mc     ->SetFillColor(2);

 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);  
  
 TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
 c1->SetFillColor(10);
 c1->  cd();   
  

 TPad* canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
 canvas_1 ->Draw();
 canvas_1 ->cd();
 
 canvas_1->SetLogy(log);
 
 hist_mc    ->Draw("hist");  
  
 hist_mc    ->GetXaxis()->SetTitle(name);
 hist_mc    ->GetYaxis()->SetTitle("entries");

 hist_mc    ->SetTitleSize(0.08,"Y");
 hist_mc    ->SetTitleOffset(0.65,"Y"); 

 hist_data->Draw("same e");


 TLegend* qw =  new TLegend(0.6,0.73,0.95,1.);
  //Legend
 qw->AddEntry(hist_data,     datacaption,                       "p");
 qw->AddEntry(hist_mc,        "MC "           ,         "f");
 
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

/*
 if (name=="nPV") {
   cout << " nPV ratio Data/MC " << endl;
   for (int i=1; i<61; i++) {
     cout <<i-1 << "  " << histo_ratio->GetBinContent(i) <<endl;
   }
 }
*/

 c1->cd();  
  
 TString name_plot=name+"_Linear"+format; 
 if(log) name_plot=name+"_Log"+format;
 c1->SaveAs(dir4plots+"/"+name_plot);

 

}
//--------------------------


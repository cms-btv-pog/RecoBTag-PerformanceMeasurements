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

TString adddir = "./";
TString filename_1=adddir+"QCD_8P8NoPU.root";
TString labelfile1="QCD 8TeV Pyth8 No PU";
TString filename_2=adddir+"QCD_8P8wPU.root";
TString labelfile2="QCD 8TeV Pyth8 w/ PU";

bool boolfile3=false;
TString filename_3=adddir+"QCD_8P6NoPU.root";
TString labelfile3="QCD 8TeV Pyth6 No PU";

bool boolfile4=false;
TString filename_4=adddir+"QCD_8P6wPU.root";
TString labelfile4="QCD 8TeV Pyth6 w/ PU";


bool boolfile5=false;
TString filename_5=adddir+"QCD_13P8NoPU.root";
TString labelfile5="QCD 13TeV Pyth8 No PU";

bool boolfile6=false;
TString filename_6=adddir+"QCD_13P8wPU.root";
TString labelfile6="QCD 13TeV Pyth8 w/ PU";



bool displayRatio = false;


TString dir4plots="plots";

//TString datacaption = "HLT_BTagMu_DiJet40_Mu5, jet pT>60";
//TString datacaption = "HLT_PFJet40, jet pT>60";
TString datacaption = "data";
TString filename = "qcd120_13tev_pu20bx25.root";
bool btag = false;


TString title = "CMS Simulation";


TString format=".gif"; // .png or .pdf or .gif
bool bOverflow=true;
bool web = true;



void Draw(TString name, TString histotitle, bool log);
void DrawMC(TString name, TString histotitle, bool log);
void DrawMCCompare(TString name, TString histotitle, bool log);
void DrawMCgen(TString name, TString histotitle, bool log);
void DrawTagRate(TString name, TString histotitle, bool log);
void Draw2DPlot(TString name, TString histotitle, TString titleX, TString titleY, bool log);
void OverFlowBinFix(TH1D* );
void DefineHstack(TString filename_x, TString name, TH1D*& hist_fill, THStack * & hs_fill, TH1D* & hist_fill_b, TH1D* & hist_fill_c, TH1D* & hist_fill_l, TH1D* & hist_fill_g, TH1D* & hist_fill_bfromg, TH1D* & hist_fill_cfromg, TH1D* & hist_fill_puu );
void DefineHline(TString filename_x, TString name, TH1D* & hist_line, int iMcolor, float xMsize, int iMstyle, int iLwidth);
void DefineHratio(TString nameRatio, TH1D* & histo_ratio1, TH1D* hnum, TH1D* hdenum, int iMcolor, float xMsize, int iMstyle, int iLwidth );
void DefineHgen(TString filename_x, TString name, TH1D* & hist_fill, int iColor);


void anplot();
void anplot2();

//--------

void DrawCommPlot(bool Draw_track_plots, bool Draw_Nminus1_plots, bool Draw_sv_plots, bool Draw_muons_plots, bool Draw_discriminator_plots , 
                  bool Draw_newdiscriminator_plots, bool Draw_tagRate_plots, bool Draw_2D_plots, bool Draw_pflepton){

    TString action = "mkdir "+dir4plots+"/";
    system(action);
  
    //cout << "Drawing pt, eta, phi plots\n";
    //Draw("jet_multi"    ,"number of jets",1);	    
    DrawMCCompare("jet_pt_all"   ,"pT of all jets",1);	    
    DrawMCCompare("jet_eta"      ,"eta of all jets",	0);    
    DrawMCCompare("jet_phi"      ,"phi of all jets",0);
    DrawMCCompare("track_multi"  ,"number of tracks in the jets", 0);
    //DrawMC("nPV"      ,"# of PV",0);

    //cout << "Drawing tack plots\n";
    if (Draw_track_plots){	      
        DrawMCCompare("track_multi"  ,      "number of tracks in the jets",0);		   
        DrawMCCompare("trk_multi_sel"  ,    "number of selected tracks in the jets",0);	   
        DrawMCCompare("track_chi2"   ,      "normalized #chi^{2} of the tracks"	,1);	   
        DrawMCCompare("track_nHit" ,      "number of hits",1);		   
        DrawMCCompare("track_HPix"   ,      "number of hits in the Pixel",1);		   
        DrawMCCompare("track_len"     ,     "Track decay length",1);				   
        DrawMCCompare("track_dist"    ,     "Track distance to the jet axis"	,1);	   
        DrawMCCompare("track_dz"     ,      "Track IP_dz",1);				   
        DrawMCCompare("track_pt"     ,      "pT of all the tracks",1);			   
        DrawMCCompare("track_pt15"     ,      "pT of all the tracks",1);			   
        //Draw("track_isfromSV",     "Track is from SV",1);			   

        DrawMCCompare("track_IPs"    ,      "3D IP significance of all tracks",1);	   
        DrawMCCompare("track_IPs1tr" ,      "3D IP significance of the first track",1);	   
        DrawMCCompare("track_IPs2tr" ,      "3D IP significance of the second track",1);     
        DrawMCCompare("track_IPs3tr" ,      "3D IP significance of the third track",1);	
        DrawMCCompare("track_IP"     ,      "3D IP of all tracks",1);			   
        DrawMCCompare("track_IP1tr"  ,      "3D IP of the first track"	,1);	   
        DrawMCCompare("track_IP2tr"  ,      "3D IP of the second track",1);		   
        DrawMCCompare("track_IP3tr"  ,      "3D IP of the third track"	,1);	
        DrawMCCompare("track_IP2Ds"	,    "2D IP significance of all tracks"	,1);   
        DrawMCCompare("track_IP2Ds1tr" ,    "2D IP significance of the first track",1);	   
        DrawMCCompare("track_IP2Ds2tr" ,    "2D IP significance of the second track" ,1);    
        DrawMCCompare("track_IP2Ds3tr" ,    "2D IP significance of the third track" ,1);
        DrawMCCompare("track_IP2D"    ,     "2D IP of all tracks",1);			   
        DrawMCCompare("track_IP2D1tr" ,     "2D IP of the first track"	,1);	   
        DrawMCCompare("track_IP2D2tr" ,     "2D IP of the second track",1);		   
        DrawMCCompare("track_IP2D3tr" ,     "2D IP of the third track"	,1);	
        DrawMCCompare("track_IP2Derr" ,     "2D IP error of all tracks",1);				   
        DrawMCCompare("track_IP2Derr1tr" ,  "2D IP error of the first track",1);  	   
        DrawMCCompare("track_IP2Derr2tr" ,  "2D IP error of the second track" ,1);	   
        DrawMCCompare("track_IP2Derr3tr" ,  "2D IP error of the third track",1);  	
        DrawMCCompare("track_IPerr"   ,     "3D IP error of all tracks",1);				   
        DrawMCCompare("track_IPerr1tr"   ,  "3D IP error of the first track" ,1); 	   
        DrawMCCompare("track_IPerr2tr"   ,  "3D IP error of the second track" ,1);	   
        DrawMCCompare("track_IPerr3tr"   ,  "3D IP error of the third track" ,1); 	
	
    }  
    if (Draw_Nminus1_plots){  
        DrawMCCompare("track_chi2_cut"    ,"Normalized #chi^{2} @N-1 step",1);	    
        DrawMCCompare("track_nHit_cut"  ," Number of hits @N-1 step",1);
        DrawMCCompare("track_HPix_cut"    ,"Number of hits in the Pixel @N-1 step",1);  
        DrawMCCompare("track_len_cut"     ,"Decay length @N-1 step",1);		    
        DrawMCCompare("track_dist_cut"    ,"Distance to the jet axis @N-1 step" ,1);  
        DrawMCCompare("track_dz_cut"     , "Transverse IP @N-1 step",1);		    
        DrawMCCompare("track_pt_cut"	   ,"Track pT @N-1 step",1);
        DrawMCCompare("track_pt15_cut"	   ,"Track pT @N-1 step",1);
    }
    if (Draw_sv_plots){

        DrawMCCompare("jet_pt_sv"    ,"pT of jets containing a SV",1);
        DrawMCCompare("sv_multi_0","nr. of SV including bin 0",1); 
        DrawMCCompare("sv_multi","nr. of SV",1);
        DrawMCCompare("sv_mass","SV mass",0);
        DrawMCCompare("sv_mass_3trk","SV mass if #tracks@SV >=3",0);
        DrawMCCompare("sv_chi2norm","SV norm. #chi^{2}",1);
        DrawMCCompare("sv_deltaR_jet","Delta R between the jet and the SV direction.",0);
        DrawMCCompare("sv_deltaR_sumJet","#Delta R between the jet and the SV {p}",0);
        DrawMCCompare("sv_deltaR_sumDir","#Delta R between the SV direction and the SV {p}",0);
        DrawMCCompare("sv_en_ratio","SV energy ratio",0);	
        DrawMCCompare("sv_aboveC","IP2D of the first track above the charm threshold",1);	
        DrawMCCompare("sv_pt","SV p_{T}",1); 	
        DrawMCCompare("sv_eta","SV #eta",0);	
        DrawMCCompare("sv_phi","SV #phi",0);	
        DrawMCCompare("sv_flight3D","SV 3D flight distance",1);	
        DrawMCCompare("sv_flight2D","SV 2D flight distance",1);	
        DrawMCCompare("sv_flight2D_3trk","SV 2D flight dist. if #tracks@SV >=3",1);	
        DrawMCCompare("sv_flight2DSig_3trk","SV 2D flight dist. sig. if #tracks@SV >=3",1);	
        DrawMCCompare("sv_flight3DSig","SV 3D flight distance significance",1);
        DrawMCCompare("sv_flightSig2D","SV 2D flight distance significance",1);
        DrawMCCompare("sv_flight3Derr","SV 3D flight distance error",1);
        DrawMCCompare("sv_flight2Derr","SV 2D flight distance error",1);
        DrawMCCompare("svnTrk","nr. of tracks from a SV",1);	
        DrawMCCompare("svnTrk_firstVxt","nr. of tracks from the first SV",1);
    }
    if (Draw_muons_plots){  
        DrawMCCompare("muon_multi"   ,      "number of muons", 1);       
        DrawMCCompare("muon_multi_sel"   ,  "number of selected muons",1);
        DrawMCCompare("mu_ptrel"     ,      "p_{T} rel. of the muon",0);    
        DrawMCCompare("mu_chi2"      ,      "norm. #chi^{2} of the muon", 1);
        DrawMCCompare("muon_Pt",	     "Muon p_{T}",1);	       
        DrawMCCompare("muon_eta",	     "Muon #eta",0);	       
        DrawMCCompare("muon_phi",	     "Muon #phi",0);	       
        DrawMCCompare("muon_Ip3d",	     "Muon 3D IP",1);	       
        DrawMCCompare("muon_Ip2d",	     "Muon 2D IP",1);	       
        DrawMCCompare("muon_Sip3d",	     "Muon 3D IP significance",1);
        DrawMCCompare("muon_Sip2d",	     "Muon 2D IP significance",1);
        DrawMCCompare("muon_DeltaR",	     "Muon1 #Delta R",0);
  
    } 
    if (Draw_discriminator_plots){  
        DrawMCCompare("TCHE_extended1"       ,"TCHE (extended)",1);    
        DrawMCCompare("TCHP_extended1"       ,"TCHP (extended)",1);    
        DrawMCCompare("discri_ssche0",      "SSVHE Discriminator", 1);
        DrawMCCompare("discri_sschp0",      "SSVHP Discriminator", 1);

        DrawMCCompare("TCHE"	      ,"TCHE Discriminator", 1); 		    
        DrawMCCompare("TCHP"	      ,"TCHP Discriminator",1);  		    
        DrawMCCompare("JP"	      ,"JP Discriminator",1);			    
        DrawMCCompare("JBP"	      ,"JBP Discriminator",1);			    
        DrawMCCompare("SSV"	      ,"SSVHE Discriminator",0);			    
        DrawMCCompare("SSVHP"        ,"SSVHP Discriminator",0); 		    
        DrawMCCompare("CSV"	      ,"CSV Discriminator",1);
  
    }
    if (Draw_newdiscriminator_plots) {
        DrawMCCompare("RetCombSvx"         ,"CSV_v1 Discriminator",1);
        DrawMCCompare("RetCombSvx"         ,"CSV_v1 Discriminator",0);
        DrawMCCompare("CombCSVJP"          ,"CSVJP Discriminator",1);
        DrawMCCompare("CombCSVSL"          ,"CSVSL_v1 Discriminator",1);
        DrawMCCompare("CombCSVJPSL"        ,"CSVJPSL Discriminator",1);
        DrawMCCompare("SoftMu"        ,"Soft Muon Discriminator",1);
        DrawMCCompare("SoftMu"        ,"Soft Muon Discriminator",0);
        DrawMCCompare("SoftEl"        ,"Soft Electron Discriminator",1);
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
        DrawTagRate("SSV"	      ,"SSVHE Discriminator", 0);
        DrawTagRate("SSVHP"         ,"SSVHP Discriminator", 0);
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

        DrawMCCompare("pfelectron_multi","number of pfelectron",1);
        DrawMCCompare("pfelectron_pt","pfelectron p_{T}",1);
        DrawMCCompare("pfelectron_eta","pfelectron #eta",0);
        DrawMCCompare("pfelectron_phi","pfelectron #phi", 0);
        DrawMCCompare("pfelectron_Sip","3D IP significance of pfelectron",1);
        DrawMCCompare("pfelectron_ptrel","pT rel. of the pfelectron",0);
        DrawMCCompare("pfelectron_ratio","ratio of pfelectron", 1);
        DrawMCCompare("pfelectron_ratiorel",   "ratioRel of pfelectron",1);
        DrawMCCompare("pfelectron_deltar", "#DeltaR(pfelectron,jet)",0);

        DrawMCCompare("pfmuon_multi","number of pfmuon",1);
        DrawMCCompare("pfmuon_goodquality","quality of pfmuons",0);
        DrawMCCompare("pfmuon_pt","pfmuon p_{T}",1);
        DrawMCCompare("pfmuon_eta","pfmuon #eta",0);
        DrawMCCompare("pfmuon_phi","pfmuon #phi", 0);
        DrawMCCompare("pfmuon_Sip","3D IP significance of pfmuon",1);
        DrawMCCompare("pfmuon_ptrel","pT rel. of the pfmuon",0);
        DrawMCCompare("pfmuon_ratio","ratio of pfmuon", 1);
        DrawMCCompare("pfmuon_ratiorel",   "ratioRel of pfmuon",1);
        DrawMCCompare("pfmuon_deltar", "#DeltaR(pfmuon,jet)",0);

    }

}
void anplot() {

    TString action = "mkdir "+dir4plots+"/";
    system(action);
  
  
    DrawMCCompare("jet_pt_all"   ,"pT of all jets",1);	    

    DrawMCCompare("track_IPs"    ,      "3D IP significance of all tracks",1);	   
    DrawMCCompare("track_IP"     ,      "3D IP of all tracks",1);			   
    DrawMCCompare("track_IPerr"   ,     "3D IP error of all tracks",1);				   
    DrawMCCompare("track_chi2_cut"    ,"Normalized #chi^{2} @N-1 step",1);	    
    DrawMCCompare("track_nHit_cut"  ," Number of hits @N-1 step",1);
    DrawMCCompare("track_HPix_cut"    ,"Number of hits in the Pixel @N-1 step",1);  
    DrawMCCompare("track_len_cut"     ,"Decay length @N-1 step",1);		    
    DrawMCCompare("track_dist_cut"    ,"Distance to the jet axis @N-1 step" ,1);  
    DrawMCCompare("track_dz_cut"     , "IP_dz @N-1 step",1);		    
    DrawMCCompare("track_pt_cut"	   ,"Track pT @N-1 step",1);
    DrawMCCompare("track_pt15_cut"	   ,"Track pT @N-1 step",1);
    DrawMCCompare("track_IP2D_cut"    ,"IP2D @N-1 step",1);

    /*
      DrawMCCompare("sv_multi_0","nr. of SV including bin 0",1); 
      DrawMCCompare("sv_mass","SV mass",0);
      DrawMCCompare("sv_deltaR_jet","Delta R between the jet and the SV direction.",0);
      DrawMCCompare("sv_pt","SV p_{T}",1); 	

      DrawMCCompare("mu_ptrel"     ,      "p_{T} rel. of the muon",0);    
      DrawMCCompare("muon_Pt",	     "Muon p_{T}",1);	       
      DrawMCCompare("muon_Sip3d",	     "Muon 3D IP significance",1);
      DrawMCCompare("muon_DeltaR",	     "Muon1 #Delta R",0);
  
      DrawMCCompare("discri_ssche0",      "SSVHE Discriminator", 1);
      DrawMCCompare("discri_sschp0",      "SSVHP Discriminator", 1);

      DrawMCCompare("TCHE"	      ,"TCHE Discriminator", 1); 		    
      DrawMCCompare("TCHP"	      ,"TCHP Discriminator",1);  		    
      DrawMCCompare("JP"	      ,"JP Discriminator",1);			    
      DrawMCCompare("JBP"	      ,"JBP Discriminator",1);			    
      DrawMCCompare("CSV"	      ,"CSV Discriminator",1);
    */
}
void anplot2() {

    TString action = "mkdir "+dir4plots+"/";
    system(action);
  
  
    DrawMCCompare("CSV"	      ,"CSV Discriminator",1);
    DrawMCCompare("CombCSVSL"          ,"CSVSL_v1 Discriminator",1);
    /*
      DrawMCCompare("RetCombSvx"         ,"CSV_v1 Discriminator",1);
      DrawMCCompare("RetCombSvx"         ,"CSV_v1 Discriminator",0);
      DrawMCCompare("CombCSVJP"          ,"CSVJP Discriminator",1);
      DrawMCCompare("CombCSVJPSL"        ,"CSVJPSL Discriminator",1);
      DrawMCCompare("SoftMu"        ,"Soft Muon Discriminator",1);
      DrawMCCompare("SoftEl"        ,"Soft Electron Discriminator",1);
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



    cout << histotitle << endl;
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
/*
    errval=0;
    if(histo->GetBinContent(1)!=0.)
        errval+=pow(histo->GetBinError(1)/histo->GetBinContent(1),2);
    if(histo->GetBinContent(0)!=0.)
        errval+=pow(histo->GetBinError(0)/histo->GetBinContent(0),2);
    errval=sqrt(errval)*val;
*/
    errval=sqrt( histo->GetBinError(1)*histo->GetBinError(1) + histo->GetBinError(0)*histo->GetBinError(0));
    histo->SetBinContent(1,val);
    histo->SetBinError(1,errval);
    histo->SetBinContent(0,0);
    histo->SetBinError(0,0);

    Int_t highbin=histo->GetNbinsX();

    val=histo->GetBinContent(highbin)+histo->GetBinContent(highbin+1);
/*  
    errval=0;
    if(histo->GetBinContent(highbin)!=0.)
        errval+=pow(histo->GetBinError(highbin)/histo->GetBinContent(highbin),2);
    if(histo->GetBinContent(highbin+1)!=0.)
        errval+=pow(histo->GetBinError(highbin+1)/histo->GetBinContent(highbin+1),2);
    errval=sqrt(errval)*val;
*/
    errval=sqrt( histo->GetBinError(highbin)*histo->GetBinError(highbin) + histo->GetBinError(highbin+1)*histo->GetBinError(highbin+1));
    histo->SetBinContent(highbin,val);
    histo->SetBinError(highbin,errval);
    histo->SetBinContent(highbin+1,0);
    histo->SetBinError(highbin+1,0);
}


void DrawMC(TString name, TString histotitle, bool log)
 
{
    TFile *myFile     = new TFile(filename);
    myFile->cd();

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


void DrawMCCompare(TString name, TString histotitle, bool log)
 
{


    TH1D* hist_fill;
    THStack *hs_fill;
    TH1D* hist_fill_b;
    TH1D* hist_fill_c;
    TH1D* hist_fill_l;
    TH1D* hist_fill_g;
    TH1D* hist_fill_bfromg;
    TH1D* hist_fill_cfromg;
    TH1D* hist_fill_puu;
    DefineHstack(filename_1, name, hist_fill, hs_fill, hist_fill_b, hist_fill_c, hist_fill_l, hist_fill_g, hist_fill_bfromg, hist_fill_cfromg, hist_fill_puu);




    TH1D* hist_line;
    DefineHline(filename_2, name, hist_line,  kOrange, 0.75, 20, 2) ;


    TH1D* hist_line2;
    if (boolfile3) DefineHline(filename_3, name, hist_line2, kGray+2, 1,1,1 );


    TH1D* hist_line4;
    if (boolfile4) DefineHline(filename_4, name, hist_line4, kViolet, 1,1,1);


    TH1D* hist_line5;
    if (boolfile5) DefineHline(filename_5, name, hist_line5, kCyan-6, 1,1,1);

    TH1D* hist_line6;
    if (boolfile6) DefineHline(filename_6, name, hist_line6, 46, 1,1,1);

    TH1D* histo_ratio1;
    DefineHratio("histo_ratio1", histo_ratio1, hist_line, hist_fill,   kOrange, 0.75, 20, 2) ;



    TH1D* histo_ratio2;
    if (boolfile3) DefineHratio("histo_ratio2", histo_ratio2, hist_line2, hist_fill,  kGray+2, 1,7,1 );
    TH1D* histo_ratio4;
    if (boolfile4) DefineHratio("histo_ratio4", histo_ratio4, hist_line4, hist_fill, kViolet, 1,7,1);
    TH1D* histo_ratio5;
    if (boolfile5)DefineHratio("histo_ratio5", histo_ratio5, hist_line5, hist_fill, kCyan-6, 1,7,1);
    TH1D* histo_ratio6;
    if (boolfile6)DefineHratio("histo_ratio6", histo_ratio6, hist_line6, hist_fill, 46, 1,7,1);



    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);  
  
    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->SetFillColor(10);
    c1->  cd();   
  

    TPad* canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
    canvas_1 ->Draw();
    canvas_1 ->cd();
 
    canvas_1->SetLogy(log);
 

    if (hist_line->GetMaximum() > hs_fill->GetMaximum()) hs_fill->SetMaximum(hist_line->GetMaximum()*1.1);

    hs_fill->Draw("hist");
    hs_fill->GetXaxis()->SetTitle(name);
    hs_fill->GetYaxis()->SetTitle("entries");
    hs_fill->GetYaxis()->SetLabelSize(0.05);
    hs_fill->GetYaxis()->SetTitleSize(0.08);
    hs_fill->GetYaxis()->SetTitleOffset(0.65); 

  

    hist_line->Draw("hist error same");
    if (boolfile3) hist_line2->Draw("hist same");
    if (boolfile4) hist_line4->Draw("hist same");
    if (boolfile5) hist_line5->Draw("hist same");
    if (boolfile6) hist_line6->Draw("hist same");


    TLegend* qw =  new TLegend(0.7,0.65,0.975,1.);
    //Legend
    qw->AddEntry(hist_fill_b, "bottom", "f");
    qw->AddEntry(hist_fill_bfromg, "bfromg", "f");
    qw->AddEntry(hist_fill_c, "charm", "f");
    qw->AddEntry(hist_fill_cfromg, "cfromg", "f");
    qw->AddEntry(hist_fill_l, "light", "f");
//    qw->AddEntry(hist_fill_puu, "PU", "f");

    qw->AddEntry(hist_line, labelfile2,"f");
    if (boolfile3) qw->AddEntry(hist_line2, labelfile3,"f");
    if (boolfile4) qw->AddEntry(hist_line4, labelfile4,"f");
    if (boolfile5) qw->AddEntry(hist_line5, labelfile5,"f");
    if (boolfile6) qw->AddEntry(hist_line6, labelfile6,"f");
 
    qw->SetFillColor(0);
    qw->Draw();
  
  
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.055);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);
    latex->DrawLatex(0.1, 0.96, title);

    TLatex *mclabels = new TLatex();
    mclabels->SetNDC();
    mclabels->SetTextSize(0.04);
    mclabels->SetTextFont(42);
    mclabels->SetTextAlign(13);

    c1->cd();  
  
    TPad* canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.32);
    canvas_2->Draw();
    canvas_2->cd();
    gPad->SetBottomMargin(0.375);
    gPad->SetGridy();

    histo_ratio1->GetXaxis()->SetTitle(histotitle);
    histo_ratio1->GetYaxis()->SetTitle("ratio");
    histo_ratio1->SetMinimum(0.5);
    histo_ratio1->SetMaximum(2.0);

    histo_ratio1->Draw("E1X0");
    
    if(boolfile3) histo_ratio2->Draw("E1X0 same");
    if(boolfile4) histo_ratio4->Draw("E1X0 same");
    if(boolfile5) histo_ratio5->Draw("E1X0 same");
    if(boolfile6) histo_ratio6->Draw("E1X0 same");

    TLegend* hrl =  new TLegend(0.7,0.7,0.975,1.);
    //Legend
    hrl->AddEntry(histo_ratio1, labelfile2+"over ref", "pl");
    if (boolfile3) hrl->AddEntry(histo_ratio2, labelfile3+"over ref", "pl");
    if (boolfile4) hrl->AddEntry(histo_ratio4, labelfile4+"over ref", "pl");
    if (boolfile5) hrl->AddEntry(histo_ratio5, labelfile5+"over ref", "pl");
    if (boolfile6) hrl->AddEntry(histo_ratio6, labelfile6+"over ref", "pl");
    hrl->SetFillColor(0);
    hrl->Draw();


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

    if (web && log) {
      canvas_1->cd();
      canvas_1->SetLogy(0);
      c1->cd();
      c1->SaveAs(dir4plots+"/"+name+"_Linear"+format);
    }


    if (displayRatio) {
      TCanvas *c2 = new TCanvas("c2", "c2",10,32,782,552);
      c2->SetFillColor(10);
      int howmuch=1;
      if (boolfile3) howmuch++;
      if (boolfile4) howmuch++;
      if (boolfile5) howmuch++;
      if (boolfile6) howmuch++;
      c2->Divide(1,howmuch);
      c2->  cd(1);   
          gPad->SetGridy();
      histo_ratio1->Draw("E1X0");
      howmuch=1;
      if (boolfile3) { 
          howmuch++;
          c2->cd(howmuch);
          gPad->SetGridy();
          histo_ratio2->SetMinimum(0.5);
	  histo_ratio2->SetMaximum(2.0);
          histo_ratio2->Draw("E1X0");
      }
      if (boolfile4) { 
          howmuch++;
          c2->cd(howmuch);
          gPad->SetGridy();
          histo_ratio4->SetMinimum(0.5);
	  histo_ratio4->SetMaximum(2.0);
          histo_ratio4->Draw("E1X0");
      }
      if (boolfile5) { 
          howmuch++;
          c2->cd(howmuch);
          gPad->SetGridy();
          histo_ratio5->SetMinimum(0.5);
	  histo_ratio5->SetMaximum(2.0);
          histo_ratio5->Draw("E1X0");
      }
      if (boolfile6) { 
          howmuch++;
          c2->cd(howmuch);
          gPad->SetGridy();
	  histo_ratio6->SetMinimum(0.5);
	  histo_ratio6->SetMaximum(2.0);
          histo_ratio6->Draw("E1X0");
      }
      c2->cd();
      c2->SaveAs(dir4plots+"/ratio_"+name+format);
    }

 

}
//--------------------------

void DrawMCgen(TString name, TString histotitle, bool log)
{
    //cout << "Getting 2025 file\n";
    
    TH1D* hist_fill;
    DefineHgen(filename_1, name, hist_fill,  2);


    //cout << "Getting 4050 file\n";
    TH1D* hist_line;
    DefineHgen(filename_2, name, hist_line, kOrange);

    TH1D* hist_line_2;
    if (boolfile3) DefineHgen(filename_3, name, hist_line_2, kGray+2);

    TH1D* hist_line_4;
    if (boolfile4) DefineHgen(filename_4, name, hist_line_4, kViolet);

    TH1D* hist_line_5;
    if (boolfile5) DefineHgen(filename_5, name, hist_line_5, kCyan-6);

    TH1D* hist_line_6;
    if (boolfile6) DefineHgen(filename_6, name, hist_line_6, 46);

    TH1D* histo_ratio1;
    DefineHratio("histo_ratio1", histo_ratio1, hist_line, hist_fill,   kOrange, 0.75, 20, 1) ;

    TH1D* histo_ratio2;
    if (boolfile3) DefineHratio("histo_ratio2", histo_ratio2, hist_line_2, hist_fill,   kGray+2, 1, 7, 1) ;

    TH1D* histo_ratio3;
    if (boolfile5) DefineHratio("histo_ratio3", histo_ratio3, hist_fill, hist_line_5,   1, 1, 7, 1) ;

    TH1D* histo_ratio4;
    if (boolfile4) DefineHratio("histo_ratio4", histo_ratio4, hist_line_4, hist_fill,   kViolet, 1, 7, 1) ;

    TH1D* histo_ratio5;
    if (boolfile5) DefineHratio("histo_ratio5", histo_ratio5, hist_line_5, hist_fill,   kCyan-6, 1, 7, 1) ;

    TH1D* histo_ratio6;
    if (boolfile6) DefineHratio("histo_ratio6", histo_ratio6, hist_line_6, hist_fill,   46, 1, 7, 1) ;

/*
    cout << "histo_ratio1" << endl;
    for (int i=1; i<histo_ratio2->GetNbinsX(); i++) {
    cout << i << "  " << histo_ratio1->GetBinContent(i) << endl;
    cout << "ratio 13/8 pthat" << endl;
    for (int i=1; i<histo_ratio3->GetNbinsX()+1; i++) {
    cout << i << "  " << histo_ratio3->GetBinContent(i) << endl;
    }
    */


    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);  
  
    TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
    c1->SetFillColor(10);
    c1->  cd();   
  

    TPad* canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
    canvas_1 ->Draw();
    canvas_1 ->cd();
 
    canvas_1->SetLogy(log);

    float max = hist_fill->GetMaximum();
    if (hist_line->GetMaximum() > max) max= hist_line->GetMaximum();
    if (boolfile3 && hist_line_2->GetMaximum() > max) max= hist_line_2->GetMaximum();
    hist_fill->SetMaximum(max*1.1);
 
    hist_fill->Draw("hist");
    hist_line->Draw("hist, same");
    if (boolfile3) hist_line_2->Draw("hist, same");
    if (boolfile4) hist_line_4->Draw("hist, same");
    if (boolfile5) hist_line_5->Draw("hist, same");
    if (boolfile6) hist_line_6->Draw("hist, same");
  
    TLegend* qw =  new TLegend(0.7,0.65,0.975,1.);
    //Legend
    qw->AddEntry(hist_fill, labelfile1, "f");
    qw->AddEntry(hist_line, labelfile2, "f");
    if (boolfile3) qw->AddEntry(hist_line_2, labelfile3, "f");
    if (boolfile4) qw->AddEntry(hist_line_4, labelfile4, "f");
    if (boolfile5) qw->AddEntry(hist_line_5, labelfile5, "f");
    if (boolfile6) qw->AddEntry(hist_line_6, labelfile6, "f");
 
    qw->SetFillColor(0);
    qw->Draw();
  
  
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.055);
    latex->SetTextFont(42); //22
    latex->SetTextAlign(13);
    latex->DrawLatex(0.1, 0.96, title);


    c1->cd();  
  
    TPad* canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.32);
    canvas_2->Draw();
    canvas_2->cd();
    gPad->SetBottomMargin(0.375);
    gPad->SetGridy();

    histo_ratio1->GetXaxis()->SetTitle(histotitle);
    histo_ratio1->SetMinimum(0.5);
    histo_ratio1->SetMaximum(2.0);

    histo_ratio1->Draw("E1X0");
    
    if (boolfile3) histo_ratio2->Draw("E1X0 same");
    if (boolfile4) histo_ratio4->Draw("E1X0 same");
    if (boolfile5) histo_ratio5->Draw("E1X0 same");
    if (boolfile6) histo_ratio6->Draw("E1X0 same");

    c1->cd();  
  
    TString name_plot=name+"_Linear"+format; 
    if(log) name_plot=name+"_Log"+format;
    c1->SaveAs(dir4plots+"/"+name_plot);

 

}
//--------------------------

void DefineHstack(TString filename_x, TString name, TH1D*& hist_fill, THStack * & hs_fill, TH1D* & hist_fill_b, TH1D* & hist_fill_c, TH1D* & hist_fill_l, TH1D* & hist_fill_g, TH1D* & hist_fill_bfromg, TH1D* & hist_fill_cfromg, TH1D* & hist_fill_puu )
{

    TFile *fileFill     = new TFile(filename_x);


    fileFill->cd();

    hist_fill_b       = (TH1D*)gROOT->FindObject(name+"_b");
    hist_fill_c       = (TH1D*)gROOT->FindObject(name+"_c");
    hist_fill_l       = (TH1D*)gROOT->FindObject(name+"_l");
    hist_fill_g       = (TH1D*)gROOT->FindObject(name+"_g");
    hist_fill_bfromg  = (TH1D*)gROOT->FindObject(name+"_bfromg");
    hist_fill_cfromg  = (TH1D*)gROOT->FindObject(name+"_cfromg");
    hist_fill_puu     = (TH1D*)gROOT->FindObject(name+"_puu");

    if (bOverflow) {
        OverFlowBinFix(hist_fill_b);
        OverFlowBinFix(hist_fill_c);
        OverFlowBinFix(hist_fill_l);
        OverFlowBinFix(hist_fill_g);
        OverFlowBinFix(hist_fill_bfromg);
        OverFlowBinFix(hist_fill_cfromg);
        OverFlowBinFix(hist_fill_puu);
    }

    hist_fill = (TH1D*) hist_fill_b->Clone();
    hist_fill->Add(hist_fill_bfromg);
    hist_fill->Add(hist_fill_c);
    hist_fill->Add(hist_fill_cfromg);
    hist_fill->Add(hist_fill_l);
    hist_fill->Add(hist_fill_g);
    hist_fill->Add(hist_fill_puu);

    float int_fill = hist_fill->Integral();
    if (int_fill > 0.){
        hist_fill->Scale(1./int_fill);
        hist_fill_b->Scale(1./int_fill);
        hist_fill_c->Scale(1./int_fill);
        hist_fill_l->Scale(1./int_fill);
        hist_fill_g->Scale(1./int_fill);
        hist_fill_puu->Scale(1./int_fill);
        hist_fill_bfromg->Scale(1./int_fill);
        hist_fill_cfromg->Scale(1./int_fill);
    }
/*
    cout << filename_x << "     " << hist_fill_b->Integral() 
    << "     " << hist_fill_c->Integral()  
    << "     " << hist_fill_l->Integral()
    << "     " << hist_fill_g->Integral()
    << "     " << hist_fill_bfromg->Integral()
    << "     " << hist_fill_cfromg->Integral() <<endl;
*/

    double titleoffsety=0.2;
    double titlesizey=0.2;
    hist_fill->GetYaxis()->SetTitleSize(titlesizey);
    hist_fill->GetYaxis()->SetTitleOffset(titleoffsety);

    hist_fill_b->SetFillColor(kRed);
    hist_fill_bfromg->SetFillColor(kCyan);
    hist_fill_c->SetFillColor(8);
    hist_fill_cfromg->SetFillColor(kMagenta);
    hist_fill_l->SetFillColor(kBlue);
    hist_fill_g->SetFillColor(kBlue);
    hist_fill_puu->SetFillColor(kYellow);
    hist_fill_b->SetLineColor(kRed);
    hist_fill_bfromg->SetLineColor(kCyan);
    hist_fill_c->SetLineColor(8);
    hist_fill_cfromg->SetLineColor(kMagenta);
    hist_fill_l->SetLineColor(kBlue);
    hist_fill_g->SetLineColor(kBlue);
    hist_fill_puu->SetLineColor(kYellow);

    hist_fill->SetLineWidth(2);
    hist_fill->SetLineColor(kBlack);
    hist_fill->SetMarkerStyle(20);  
    hist_fill->SetMarkerSize(0.75);
    hist_fill->SetMarkerColor(kBlack);

    hs_fill = new THStack("hs_fill","Stacked histograms AOD");
    hs_fill->Add(hist_fill_b);
    hs_fill->Add(hist_fill_bfromg);
    hs_fill->Add(hist_fill_c);
    hs_fill->Add(hist_fill_cfromg);
    hs_fill->Add(hist_fill_l);
    hs_fill->Add(hist_fill_g);
    hs_fill->Add(hist_fill_puu);


}

void DefineHline(TString filename_x, TString name, TH1D* & hist_line, int icolor, float xMsize, int iMstyle, int iLwidth )
{

    TFile *fileLine     = new TFile(filename_x);

    TH1D* hist_line_b;
    TH1D* hist_line_c;
    TH1D* hist_line_l;
    TH1D* hist_line_g;
    TH1D* hist_line_bfromg;
    TH1D* hist_line_cfromg;
    TH1D* hist_line_puu;

    fileLine->cd();

    hist_line_b       = (TH1D*)gROOT->FindObject(name+"_b");
    hist_line_c       = (TH1D*)gROOT->FindObject(name+"_c");
    hist_line_l       = (TH1D*)gROOT->FindObject(name+"_l");
    hist_line_g       = (TH1D*)gROOT->FindObject(name+"_g");
    hist_line_bfromg  = (TH1D*)gROOT->FindObject(name+"_bfromg");
    hist_line_cfromg  = (TH1D*)gROOT->FindObject(name+"_cfromg");
    hist_line_puu  = (TH1D*)gROOT->FindObject(name+"_puu");

    if (bOverflow) {
        OverFlowBinFix(hist_line_b);
        OverFlowBinFix(hist_line_c);
        OverFlowBinFix(hist_line_l);
        OverFlowBinFix(hist_line_g);
        OverFlowBinFix(hist_line_bfromg);
        OverFlowBinFix(hist_line_cfromg);
        OverFlowBinFix(hist_line_puu);
    }
    
    hist_line = (TH1D*) hist_line_b->Clone();
    hist_line->Add(hist_line_bfromg);
    hist_line->Add(hist_line_c);
    hist_line->Add(hist_line_cfromg);
    hist_line->Add(hist_line_l);
    hist_line->Add(hist_line_g);
    hist_line->Add(hist_line_puu);

    float int_line = hist_line->Integral();
    if (int_line > 0.){
        hist_line->Scale(1./int_line);
        hist_line_b->Scale(1./int_line);
        hist_line_c->Scale(1./int_line);
        hist_line_l->Scale(1./int_line);
        hist_line_g->Scale(1./int_line);
        hist_line_puu->Scale(1./int_line);
        hist_line_bfromg->Scale(1./int_line);
        hist_line_cfromg->Scale(1./int_line);
    }

/*
    cout << filename_x << "     " << hist_line_b->Integral() 
    << "     " << hist_line_c->Integral()  
    << "     " << hist_line_l->Integral()
    << "     " << hist_line_g->Integral()
    << "     " << hist_line_bfromg->Integral()
    << "     " << hist_line_cfromg->Integral() <<endl;
*/

    double titleoffsety=0.2;
    double titlesizey=0.2;
    hist_line->GetYaxis()->SetTitleSize(titlesizey);
    hist_line->GetYaxis()->SetTitleOffset(titleoffsety);

    hist_line->SetLineWidth(iLwidth);
    hist_line->SetLineColor(icolor);
    hist_line->SetMarkerStyle(iMstyle);  
    hist_line->SetMarkerSize(xMsize);
    hist_line->SetMarkerColor(icolor);

/*
    hist_line_b->SetLineColor(kRed+2);
    hist_line_b->SetMarkerColor(kRed+2);
    hist_line_b->SetLineWidth(2);
    hist_line_b->SetFillStyle(0);
    hist_line_bfromg->SetLineColor(kCyan+2);
    hist_line_bfromg->SetMarkerColor(kCyan+2);
    hist_line_bfromg->SetLineWidth(2);
    hist_line_bfromg->SetFillStyle(0);
    hist_line_c->SetLineColor(kGreen+2);
    hist_line_c->SetMarkerColor(kGreen+2);
    hist_line_c->SetLineWidth(2);
    hist_line_c->SetFillStyle(0);
    hist_line_cfromg->SetLineColor(kMagenta+3);
    hist_line_cfromg->SetMarkerColor(kMagenta+3);
    hist_line_cfromg->SetLineWidth(2);
    hist_line_cfromg->SetFillStyle(0);
    hist_line_l->SetLineColor(kBlack);
    hist_line_l->SetMarkerColor(kBlack);
    hist_line_l->SetLineWidth(2);
    hist_line_l->SetFillStyle(0);
    hist_line_g->SetLineColor(kBlack);
    hist_line_g->SetMarkerColor(kBlack);
    hist_line_g->SetLineWidth(2);
    hist_line_g->SetFillStyle(0);
    hist_line_puu->SetLineColor(kYellow);
    hist_line_puu->SetMarkerColor(kYellow);
    hist_line_puu->SetLineWidth(2);
    hist_line_puu->SetFillStyle(0);
*/

}

void DefineHratio(TString nameRatio, TH1D* & histo_ratio1, TH1D* hnum, TH1D* hdenum, int icolor, float xMsize, int iMstyle, int iLwidth )
{
    histo_ratio1 = (TH1D*) hnum->Clone();
    histo_ratio1->SetName(nameRatio);
    histo_ratio1->SetTitle(""); 
    histo_ratio1->Divide(hdenum);

    double titlesizex=0.17;
    double labelsizex=0.14;
    histo_ratio1->GetYaxis()->SetNdivisions( 505 );
    histo_ratio1->GetXaxis()->SetLabelSize( labelsizex);
    histo_ratio1->GetXaxis()->SetTitleSize( titlesizex );
    histo_ratio1->GetYaxis()->SetLabelSize( 0.1);

    histo_ratio1->SetLineWidth(iLwidth);
    histo_ratio1->SetLineColor(icolor);
    histo_ratio1->SetMarkerStyle(iMstyle);  
    histo_ratio1->SetMarkerSize(xMsize);
    histo_ratio1->SetMarkerColor(icolor);
}

void DefineHgen(TString filename_x, TString name, TH1D* & hist_fill, int iColor)
{

    TFile *fileFill     = new TFile(filename_x);
    fileFill->cd();
    hist_fill       = (TH1D*)gROOT->FindObject(name);
    if (bOverflow) {
        OverFlowBinFix(hist_fill);
    }

    float int_fill = hist_fill->Integral();
    hist_fill->Scale(1./int_fill);

    double titleoffsety=0.65;
    double titlesizey=0.08;
     
    hist_fill->GetYaxis()->SetTitleSize(titlesizey);
    hist_fill->GetYaxis()->SetTitleOffset(titleoffsety);
  
    hist_fill->SetLineColor(iColor);
}

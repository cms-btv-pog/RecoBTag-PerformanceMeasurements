#define CommPlotProducer4ttbar_cxx
#include "CommPlotProducer4ttbar.h"

#include <TH2.h>
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <iomanip>
#include "TLorentzVector.h"
#include <string.h>
#include "THStack.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>

#include "../common.h"

using namespace std;


//------------------------------------------------------------------------------//
 //--------------------   Set cross sections  ---------------------------------//
  //--------------------------------------------------------------------------//
void CommPlotProducer4ttbar::SetNorm(float xnorm) 
{
   x_section[0]=xnorm;
}


void CommPlotProducer4ttbar::Loop(int datatype, int trig_data, float PtMin_Cut, float PtMax_Cut, TString output_name, TH1F* wgtcounter, TString syst)
{ 
 
  if (syst != "") cout << "Running over " << syst << "..." << endl; 
  else            cout << "Running over nominal sample" << endl; 

 
  //---------------Configuration-----------------------------------------// 
  float PtMin  = PtMin_Cut;  
  float PtMax  = PtMax_Cut;  
  float EtaCut = 2.4; 
  double pi=acos(-1);
  
  // Cross check variables--------------------------------------------------//
  int   Nevent = 0; 
  float N_event_mc_before_sel   =0;
  float N_event_data_before_sel =0;
  float N_event_mc_after_sel    =0;
  float N_event_data_after_sel  =0;
  
  
  int njet_c     =0;
  int njet_b     =0;   
  int njet_bfromg=0;    
  int njet_l     =0;  
  int njet_pu    =0;  
  int njet_mc    =0; 
  int njet_data  =0; 
  int nSVbins    = 50; 
  int nTrackbins = 100; 

  
  bool passNhit;
  bool passPix    ; 
  bool passIPz    ;
  bool passPt     ; 
  bool passnormchi2; 
  bool passtrkdist ; 
  bool passtrklen  ;
  bool passTrackIP2D;
   
  TString tmpsyst = "";
  if(syst != "") tmpsyst = "_"+syst;

  //---------------------------------------------------------------------//
  TFile *myfile=new TFile(output_name+tmpsyst+".root",      "recreate");
  
  // --------------------------------------Histograms declaration------------------------------------------------//   
  TH1D* nPU_mc                  = new TH1D("nPU_mc",                "nPU_mc",                50,-0.5,49.5 );
  TH1D* nPU_data                = new TH1D("nPU_data",              "nPU_data",              50,-0.5,49.5 );
  TH1D* nPV_mc                  = new TH1D("nPV_mc",                "nPV_mc",                50,-0.5,49.5 );
  TH1D* pt_hat                  = new TH1D("pt_hat",                "pt_hat",                80,   0,800  );
  TH1D* jet_pt_mc               = new TH1D("jet_pt_mc",  	    "jet_pt_mc", 	     80,   0,PtMax);
  
  // --------------------------------------Histograms declaration -----------------------------------------//
  if(!produceCTagTree){ 
  AddHistottbar("nPV",            "number of PV",               50,-0.5,49.5, syst);
  AddHistottbar("nPV_unweighted", "unweighted number of PV",    50,-0.5,49.5, syst);
  AddHistottbar("met",            "MET",                        30,  0.,300., syst);
  AddHistottbar("mll",            "M_{ll}",                     60,  0.,300., syst);
  AddHistottbar("njet",           "number of jets",	        10,-0.5, 9.5, syst);
  AddHistottbar("njet_pt30",      "number of jets pt30",        10,-0.5, 9.5, syst);
  AddHistottbar("pt_e",           "P_{T}^{e}",                  50,  0.,200., syst);
  AddHistottbar("pt_mu",          "P_{T}^{#mu}",                50,  0.,200., syst);
  AddHistottbar("pt_jet",         "P_{T}^{leading jet}",        50,  0.,400., syst);

  // HIP check (as function of run range for Run2016B)
  AddHistottbar("nEvt_run",          "number of evt VS run",                  20,  0,  20, syst);
  AddHistottbar("nEvt_run_CSVv2L",   "number of evt VS run(b-jet csvl)",      20,  0,  20, syst);
  AddHistottbar("nEvt_run_CSVv2M",   "number of evt VS run(b-jet csvm)",      20,  0,  20, syst);
  AddHistottbar("nEvt_run_CSVv2T",   "number of evt VS run(b-jet csvt)",      20,  0,  20, syst);

  AddHistottbar("nbtag_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_CSVv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );


   
  AddHistottbar("nbtag_all_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_all_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_all_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_all_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

    
  AddHistottbar("nbtag_2b_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2b_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2b_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2b_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

     
  AddHistottbar("nbtag_1b1c_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1b1c_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1b1c_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1b1c_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

     
  AddHistottbar("nbtag_1b1l_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1b1l_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1b1l_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1b1l_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

     
  AddHistottbar("nbtag_2c_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2c_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2c_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2c_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

     
  AddHistottbar("nbtag_1c1l_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1c1l_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1c1l_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );
                           
  AddHistottbar("nbtag_1c1l_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

      
  AddHistottbar("nbtag_2l_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2l_Inf60_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_Inf60_afterJetSel_cMVAv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_Inf60_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2l_60-120_afterJetSel_cMVAv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_60-120_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_60-120_afterJetSel_cMVAv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );
                         
  AddHistottbar("nbtag_2l_120-320_afterJetSel_cMVAv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_120-320_afterJetSel_cMVAv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_120-320_afterJetSel_cMVAv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

 

  
  AddHistottbar("nbtag_all_afterJetSel_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_afterJetSel_CSVv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_afterJetSel_CSVv2M_SFapplied","number of btag jets (medium WP)",  	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_afterJetSel_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_all_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_all_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_all_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_all_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

 
  AddHistottbar("nbtag_2b_afterJetSel_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_afterJetSel_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2b_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2b_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2b_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2b_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  
  AddHistottbar("nbtag_1b1c_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1b1c_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1b1c_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",    6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1b1c_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",    6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1c_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",     6,-0.5,5.5 , syst   );

  
  AddHistottbar("nbtag_1b1l_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1b1l_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1b1l_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",    6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1b1l_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",    6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1b1l_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",     6,-0.5,5.5 , syst   );


  AddHistottbar("nbtag_2c_afterJetSel_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_afterJetSel_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2c_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2c_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2c_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2c_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

    
  AddHistottbar("nbtag_1c1l_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1c1l_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1c1l_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",    6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_1c1l_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",    6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_1c1l_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",     6,-0.5,5.5 , syst   );


  AddHistottbar("nbtag_2l_afterJetSel_CSVv2T","number of btag jets (tight WP)",		     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_afterJetSel_CSVv2L","number of btag jets (loose WP)",		     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2l_Inf60_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_Inf60_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_Inf60_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2l_60-120_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_60-120_afterJetSel_CSVv2M","number of btag jets (medium WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_60-120_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   );

  AddHistottbar("nbtag_2l_120-320_afterJetSel_CSVv2T","number of btag jets (tight WP)",	     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_120-320_afterJetSel_CSVv2M","number of btag jets (medium WP)",     6,-0.5,5.5 , syst   );
  AddHistottbar("nbtag_2l_120-320_afterJetSel_CSVv2L","number of btag jets (loose WP)",	     6,-0.5,5.5 , syst   ); 

  
  AddHisto("sv_deltaR_jet",      "sv_deltaR_jet",                                       50,0.,0.5  , syst  );
  AddHisto("sv_deltaR_sumJet",   "SVvtxSumJetDeltaR",                                   50,0.,0.5 , syst   );
  AddHisto("sv_deltaR_sumDir",   "SVvtxSumVtxDirDeltaR",                                50,0.,0.5 , syst   );
  AddHisto("sv_en_ratio",        "Fractional energy",                                   50,0.,1.  , syst   );  
  AddHisto("sv_pt",              "Vtx p_{T}",                                           50,0.,100. , syst  );
  AddHisto("sv_eta",             "Vtx #eta",                                            50, -2.5, 2.5, syst);
  AddHisto("sv_phi",             "Vtx #phi",                                            40, -1.*pi,pi, syst);
  AddHisto("sv_flightSig2D",     "Flight distance significance 2D",                     50,0.,80.   , syst );
  AddHisto("sv_flight2D",        "Flight distance 2D",                                  50,0.,2.5   , syst );
  AddHisto("sv_flight3D",        "Flight distance 3D",                                  50,0.,15.   , syst );  
  AddHisto("sv_flight3DSig" ,    "flight distance significance 3D",	                50,0.,80.   , syst );
  AddHisto("sv_multi_0"	  ,      "number of secondary vertex",                          6,-0.5,5.5  , syst );
  AddHisto("sv_multi"	  ,      "number of secondary vertex",                          6,-0.5,5.5  , syst );
  AddHisto("sv_mass"	  ,      "invariant mass of the secondary vertex",              50,0.,8.    , syst );
  AddHisto("sv_chi2norm"  ,      "normalized chi2 of the secondary vertex",             50,0.,10.   , syst );
  AddHisto("sv_tot_charge",      "Total charge",                                                        21,-10.5,10.5, syst );
  AddHisto("svnTrk",	         "Track multiplicity : SVnVertexTracks (centered)",                     13, -0.5,12.5, syst );
  AddHisto("sv_flight3Derr",     "Flight distance error 3D",                                            50,   0.,0.2 , syst );
  AddHisto("sv_flight2Derr",     "Flight distance error 2D",                                            50,   0.,0.05, syst );
  AddHisto("sv_mass_3trk"	,"invariant mass of the secondary vertex with at least 3 SV tracks",    50,   0.,8.  , syst );
  
  AddHisto("track_multi"  ,      "number of tracks in the jets",                40,-0.5,39.5 , syst );
  AddHisto("trk_multi_sel"  ,    "number of selected tracks in the jets",       40,-0.5,39.5 , syst );
  AddHisto("track_chi2"   ,      "normalized chi2 of the tracks",               100,0.,30.   , syst );
  AddHisto("track_nHit" ,        "number of hits ",               35,-0.5, 34.5 , syst);
  AddHisto("track_HPix"   ,      "number of hits in the Pixel",                 10,-0.5, 9.5 , syst );
  
  AddHisto("track_IPs"    ,      "3D IP significance of all tracks",	        100,-35.,35. , syst );
  AddHisto("track_IPs1tr" ,      "3D IP significance of the first track",       100,-35.,35. , syst );
  AddHisto("track_IPs2tr" ,      "3D IP significance of the second track",      100,-35.,35. , syst );
  AddHisto("track_IP"     ,      "3D IP of all tracks",	                        100,-0.1,0.1 , syst );
  AddHisto("track_IP1tr"  ,      "3D IP of the first track",                    100,-0.1,0.1 , syst );
  AddHisto("track_IP2tr"  ,      "3D IP of the second track",                   100,-0.1,0.1 , syst ); 
  AddHisto("track_IP2Ds"    ,    "2D IP significance of all tracks",	        100,-35.,35. , syst );
  AddHisto("track_IP2Ds1tr" ,    "2D IP significance of the first track",       100,-35.,35. , syst );
  AddHisto("track_IP2Ds2tr" ,    "2D IP significance of the second track",      100,-35.,35. , syst );
  AddHisto("track_IP2D"    ,     "2D IP of all tracks",	                        100,-0.1,0.1 , syst );
  AddHisto("track_IP2D1tr" ,     "2D IP of the first track",                    100,-0.1,0.1 , syst );
  AddHisto("track_IP2D2tr" ,     "2D IP of the second track",                   100,-0.1,0.1 , syst );
  AddHisto("track_IP2Derr1tr" ,  "2D IP error of the first track",              100,0,0.1    , syst );    
  AddHisto("track_IPerr1tr"   ,  "3D IP error of the first track",              100,0,0.1    , syst );  
  AddHisto("track_IP2Derr2tr" ,  "2D IP error of the second track",             100,0,0.1    , syst );    
  AddHisto("track_IPerr2tr"   ,  "3D IP error of the second track",             100,0,0.1    , syst ); 
  AddHisto("track_IP2Derr" ,     "2D IP error",                                 100,0,0.1    , syst );    
  AddHisto("track_IPerr"   ,     "3D IP error",                                 100,0,0.1    , syst ); 
  AddHisto("track_IPs3tr" ,      "3D IP significance of the third track",       100,-35.,35. , syst );
  AddHisto("track_IP3tr"  ,      "3D IP of the third track",                    100,-0.1,0.1 , syst );
  AddHisto("track_IPerr3tr"   ,  "3D IP error of the third track",              100,0,0.1    , syst ); 
  AddHisto("track_IP2Ds3tr" ,    "2D IP significance of the second track",      100,-35.,35. , syst );
  AddHisto("track_IP2D3tr" ,     "2D IP of the third track",                    100,-0.1,0.1 , syst );
  AddHisto("track_IP2Derr3tr" ,  "2D IP error of the third track",              100,0,0.1    , syst );   
   
  AddHisto("track_len"     ,     "decay length",		                100,0,25.    , syst );
  AddHisto("track_dist"    ,     "distance to the jet axis",                    100,0.,0.3   , syst );
  AddHisto("track_dz"     ,      "transverse IP",                               100,-20.,20. , syst );  
  AddHisto("track_isfromSV",     "Track is from SV",                            2,-0.5, 1.5  , syst );  
  AddHisto("track_pt"	  ,      "pT of all the tracks",	                80,0.,200.   , syst );
  AddHisto("track_chi2_cut"     ,"normalized chi2 ",  	                        100,0.,30.   , syst );
  AddHisto("track_nHit_cut"     ,"number of hits  ",                            35,-0.5, 34.5, syst );
  AddHisto("track_HPix_cut"     ,"number of hits in the Pixel ",                10,-0.5, 9.5 , syst );
  AddHisto("track_len_cut"      ,"decay length ",		                100,0,25.    , syst );
  AddHisto("track_dist_cut"     ,"distance to the jet axis ",                   100,0.,0.3   , syst );
  AddHisto("track_dz_cut"       ,"transverse IP ",		                100,-20., 20., syst );  
  AddHisto("track_pt_cut"       ,"pT ",	                                        80,0.,200.   , syst );
  AddHisto("track_IP2D_cut"     ,"IP2D ",	                                100,-0.1,0.1 , syst );
   
  AddHisto("TCHE_extended1"	  ,"TCHE_extended1",				70, -30.,30. , syst );
  AddHisto("TCHP_extended1"	  ,"TCHP_extended1",				70, -30.,30. , syst );
  AddHisto("TCHE_extended2"	  ,"TCHE_extended2",				100,-30.,30. , syst );
  AddHisto("TCHP_extended2"	  ,"TCHP_extended2",				100,-30.,30. , syst );
  AddHisto("discri_ssche0",	   "SSVHE Discriminator",                       80, -1., 7.  , syst ); 
  AddHisto("discri_sschp0",	   "SSVHP Discriminator",                       80, -1., 7.  , syst ); 
   
  AddHisto("TCHE"	  ,"TCHE",				     50,0.,30. , syst);
  AddHisto("TCHP"	  ,"TCHP",				     50,0.,30. , syst);  
  AddHisto("JP" 	  ,"JP",				     30,0.,1.5 , syst);
  AddHisto("JBP"	  ,"JBP",				     25,0.,4.  , syst);
  AddHisto("SSV"	  ,"SSVHE",				     70,0.,7.  , syst);
  AddHisto("SSVHP"	  ,"SSVHP",				     70,0.,7.  , syst);
  AddHisto("CSV"	  ,"CSV",				     50,0.,1.  , syst);
  AddHisto("CSVv2"	  ,"CSVv2",				     50,0.,1.  , syst);
  AddHisto("CSVv2_pu"	  ,"CSVv2_pu",				     50,0.,1.  , syst);
  AddHisto("cMVAv2"	  ,"cMVAv2",				     50,-1.,1. , syst);
  AddHisto("CvsB"     ,"CvsB",                                       50,-1.,1. , syst);
  AddHisto("CvsL"     ,"CvsL",                                       50,-1.,1. , syst);

  AddHisto("SoftMu"       ,"SoftMu",                                 50,0.,1.  , syst);
  AddHisto("SoftEl"       ,"SoftEl",                                 50,0.,1.  , syst);
  
  AddHisto("pfmuon_multi",      "number of pfmuons",	        7,    -0.5, 6.5, syst );
  AddHisto("pfmuon_goodquality","quality of pfmuons",           3,    -0.5, 2.5, syst );
  AddHisto("pfmuon_pt",		"pfmuon p_{T}",  	        50,      0, 100, syst );
  AddHisto("pfmuon_eta",   	"pfmuon #eta",  	        50,   -2.5, 2.5, syst );  
  AddHisto("pfmuon_phi",        "pfmuon #phi",                  40, -1.*pi,  pi, syst );
  AddHisto("pfmuon_Sip",	"3D IP significance of pfmuon", 50,    -35,  35, syst );
  AddHisto("pfmuon_ptrel",      "pT rel. of the muon",	        50,      0,   5, syst );
  AddHisto("pfmuon_ratio",      "ratio of pfmuon",              50,      0,   2, syst );  
  AddHisto("pfmuon_ratiorel",   "ratioRel of pfmuon",           50,      0,0.05, syst );  
  AddHisto("pfmuon_deltar",	"#DeltaR(pfmuon,jet)",          50,      0, 0.5, syst );
  
  AddHisto("pfelectron_multi",  "number of pfelectron",	        7,    -0.5, 6.5, syst );
  AddHisto("pfelectron_pt",	"pfelectron p_{T}",  	        50,      0, 100, syst );
  AddHisto("pfelectron_eta",   	"pfelectron #eta",  	        50,   -2.5, 2.5, syst );  
  AddHisto("pfelectron_phi",    "pfelectron #phi",              40, -1.*pi,  pi, syst );
  AddHisto("pfelectron_ptrel",  "pT rel. of the pfelectron",	50,      0,   5, syst );
  AddHisto("pfelectron_ratio",  "ratio of pfelectron",          50,      0,   2, syst );  
  AddHisto("pfelectron_ratiorel","ratioRel of pfelectron",      50,      0,0.05, syst );  
  AddHisto("pfelectron_deltar",	"#DeltaR(pfelectron,jet)",      50,      0, 0.5, syst );  

  
  AddHisto2D("seltrack_vs_jetpt", "sel track multiplicity vs jet pt",         30,60,1000, 100,-0.5,99.5,syst );
  AddHisto2D("sv_mass_vs_flightDist3D", " SVMass vs SV 3D flight distance ",  100,0, 10,100,0,6,        syst );			
  AddHisto2D("avg_sv_mass_vs_jetpt","Avg SVMass vs jet pt",                   30,60,1000, 100,0,6,      syst );
  AddHisto2D("sv_deltar_jet_vs_jetpt","SVJetDeltaR vs jet pt",                25,60,300, 50,0.,0.5,     syst );  
  AddHisto2D("sv_deltar_sum_jet_vs_jetpt","SVvtxSumJetDeltaR vs jet pt",      25,60,300, 50,0.,0.5,     syst );
  AddHisto2D("sv_deltar_sum_dir_vs_jetpt","SVvtxSumVtxDirDeltaR vs jet pt",   25,60,300, 50,0.,0.5,     syst ); 

  AddHisto("tagvarCSV_vertexCategory",          "vertex category",                      3, -0.5, 2.5, syst );
  AddHisto("tagvarCSV_Sig2dAboveCharm",         "IP significance 2D charm",       nSVbins, -35.,35. , syst );
  AddHisto("tagvarCSV_trackEtaRel",             "Track etaRel",                        40,   0.,8.  , syst );
  AddHisto("tagvarCSV_trackSumJetEtRatio",      "Track  SumJet ET ratio",              40,   0.,1.5 , syst );
  AddHisto("tagvarCSV_trackSumJetDeltaR",       "Track  SumJet Delta r",               40,   0.,0.5 , syst );

  AddHisto("tagvarCSV_vertexmass_cat0",         "SV mass",                        nSVbins,   0.,8.  , syst );
  AddHisto("tagvarCSV_vertexmass3trk_cat0",     "SV mass (at least 3 SV tracks)", nSVbins,   0.,8.  , syst );
  AddHisto("tagvarCSV_vertexNTracks_cat0",      "# SV tracks",                         13, -0.5,12.5, syst );
  AddHisto("tagvarCSV_energyratio",             "Fractional energy",              nSVbins,   0.,1.  , syst );
  AddHisto("tagvarCSV_trackSip3dSig",           "3D IP significance",          nTrackbins, -35.,35. , syst );
  AddHisto("tagvarCSV_2DsigFlightDist_cat0",    "Flight distance significance 2D",nSVbins,   0.,80. , syst );
  AddHisto("tagvarCSV_vertexJetDeltaR_cat0",    "DeltaR(SV,jet) ",                nSVbins,   0.,0.4 , syst );
  } //end !produceCTagTree
 
  AddHisto("jet_multi"    ,"number of jets",                 20,        0,      20,   syst);
  AddHisto("jet_pt_all"   ,"pT of all jets",                 PtMax/10,  0,      PtMax,syst);
  AddHisto("genjet_pt_all"        ,"genpT of all jets",         50     ,  -0.5,    49.5,syst);
  AddHisto("jet_pt_sv"    ,"pT of jets containing a SV",     PtMax/10,  0,      PtMax,syst);
  AddHisto("jet_eta"      ,"eta of all jets",                50,        -2.5,   2.5,  syst);
  AddHisto("jet_phi"      ,"phi of all jets",                40,        -1.*pi, pi,   syst);

  //CTag Comm//
  if(produceCTagTree){
  AddHisto("CTag_tagvarCSV_vertexCategory",          "vertex category CSV",                      3, -0.5, 2.5, syst );
  AddHisto("CTag_tagvarCSV_Sig2dAboveCharm",         "IP significance 2D charm CSV",       nSVbins, -35.,35. , syst );
  AddHisto("CTag_tagvarCSV_trackEtaRel",             "Track etaRel CSV",                        40,   0.,8.  , syst );
  AddHisto("CTag_tagvarCSV_trackSumJetEtRatio",      "Track  SumJet ET ratio CSV",              40,   0.,1.5 , syst );
  AddHisto("CTag_tagvarCSV_trackSumJetDeltaR",       "Track  SumJet Delta r CSV",               40,   0.,0.5 , syst );

  AddHisto("CTag_tagvarCSV_vertexmass_cat0",         "SV mass CSV",                        nSVbins,   0.,8.  , syst );
  AddHisto("CTag_tagvarCSV_vertexmass3trk_cat0",     "SV mass (at least 3 SV tracks) CSV", nSVbins,   0.,8.  , syst );
  AddHisto("CTag_tagvarCSV_vertexNTracks_cat0",      "# SV tracks CSV",                         13, -0.5,12.5, syst );
  AddHisto("CTag_tagvarCSV_energyratio",             "Fractional energy CSV",              nSVbins,   0.,1.  , syst );
  AddHisto("CTag_tagvarCSV_trackSip3dSig",           "3D IP significance CSV",          nTrackbins, -35.,35. , syst );
  AddHisto("CTag_tagvarCSV_2DsigFlightDist_cat0",    "Flight distance significance 2D CSV",nSVbins,   0.,80. , syst );
  AddHisto("CTag_tagvarCSV_vertexJetDeltaR_cat0",    "DeltaR(SV,jet) CSV",                nSVbins,   0.,0.4 , syst );

  AddHisto("JP"           ,"JP",                                     30,0.,1.5 , syst);
  AddHisto("CSVv2"        ,"CSVv2",                                  50,0.,1.  , syst);
  AddHisto("CSVv2_pu"     ,"CSVv2_pu",                               50,0.,1.  , syst);
  AddHisto("cMVAv2"       ,"cMVAv2",                                 50,-1.,1. , syst);

  AddHisto("CvsB"         ,"CvsB",                                   50,-1.,1. , syst);
  AddHisto("CvsBN"        ,"CvsBN",                                  50,-1.,1. , syst);
  AddHisto("CvsBP"        ,"CvsBP",                                  50,-1.,1. , syst);
  AddHisto("CvsL"         ,"CvsL",                                   50,-1.,1. , syst);
  AddHisto("CvsLN"        ,"CvsLN",                                  50,-1.,1. , syst);
  AddHisto("CvsLP"        ,"CvsLP",                                  50,-1.,1. , syst);
  AddHisto("CTag_jetNTracks"	        ,"CTag_jetNTracks",	 	40,-0.5,39.5 , syst );
  AddHisto("CTag_jetNTracksEtaRel"      ,"CTag_jetNTracksEtaRel",	40,-0.5,39.5 , syst );
  AddHisto("CTag_jetNLeptons"           ,"CTag_jetNLeptons",		7,-0.5, 6.5, syst );
  AddHisto("CTag_trackSumJetEtRatio"    ,"CTag_trackSumJetEtRatio",     40,   0.,1.5 , syst );
  AddHisto("CTag_trackSumJetDeltaR"     ,"CTag_trackSumJetDeltaR",	40,   0.,0.3 , syst );
  AddHisto("CTag_trackSip2dSigAboveCharm","CTag_trackSip2dSigAboveCharm",50, -35,35, syst );
  AddHisto("CTag_trackSip3dSigAboveCharm","CTag_trackSip3dSigAboveCharm",50, -35,35, syst );
  AddHisto("CTag_vertexCategory"	,"CTag_vertexCategory",		3, -0.5, 2.5, syst );
  AddHisto("CTag_jetNSecondaryVertices" ,"CTag_jetNSecondaryVertices",	6,-0.5,5.5  , syst );
  AddHisto("CTag_vertexMass"		,"CTag_vertexMass",		nSVbins,   0.,8.  , syst );
  AddHisto("CTag_vertexNTracks"		,"CTag_vertexNTracks",		13, -0.5,12.5, syst );
  AddHisto("CTag_vertexEnergyRatio"	,"CTag_vertexEnergyRatio",	nSVbins,   0.,1.  , syst );
  AddHisto("CTag_vertexJetDeltaR"	,"CTag_vertexJetDeltaR",	nSVbins,   0.,0.4 , syst );
  AddHisto("CTag_flightDistance2dSig"   ,"CTag_flightDistance2dSig",	nSVbins,   0.,80. , syst );
  AddHisto("CTag_flightDistance3dSig"   ,"CTag_flightDistance3dSig",    nSVbins,   0.,80. , syst );
  AddHisto("CTag_massVertexEnergyFraction","CTag_massVertexEnergyFraction", 50, 0, 1, syst );  
  AddHisto("CTag_vertexBoostOverSqrtJetPt","CTag_vertexBoostOverSqrtJetPt", 50, 0, 1, syst );
  AddHisto("CTag_vertexLeptonCategory"  ,"CTag_vertexLeptonCategory",   7,-0.5, 6.5, syst );
  AddHisto("CTag_trackPtRel"  ,"CTag_trackPtRel", 40, 0 ,10  , syst );
  AddHisto("CTag_trackPPar"  ,"CTag_trackPPar", 40, 0, 200  , syst );
  AddHisto("CTag_trackDeltaR"  ,"CTag_trackDeltaR", 40, 0, 0.4  , syst );
  AddHisto("CTag_trackPtRatio"  ,"CTag_trackPtRatio", 30, 0, 0.3  , syst );
  AddHisto("CTag_trackPParRatio"  ,"CTag_trackPParRatio", 40, 0.95, 1  , syst );
  AddHisto("CTag_trackSip2dSig"  ,"CTag_trackSip2dSig", 40, -30, 30  , syst );
  AddHisto("CTag_trackSip3dSig"  ,"CTag_trackSip3dSig", 40, -30, 30  , syst );
  AddHisto("CTag_trackDecayLenVal"  ,"CTag_trackDecayLenVal", 40, 0, 5  , syst );
  AddHisto("CTag_trackJetDistVal"  ,"CTag_trackJetDistVal", 40, -0.07, 0  , syst );
  AddHisto("CTag_trackEtaRel"  ,"CTag_trackEtaRel", 40, 0, 8  , syst );
  AddHisto("CTag_leptonPtRel"  ,"CTag_leptonPtRel", 40, 0 ,6  , syst );
  AddHisto("CTag_leptonSip3d"  ,"CTag_leptonSip3d", 50, -1, 1  , syst );
  AddHisto("CTag_leptonDeltaR"  ,"CTag_leptonDeltaR", 40 , 0, 0.4  , syst );
  AddHisto("CTag_leptonRatioRel"  ,"CTag_leptonRatioRel", 40, 0 , 0.02  , syst );
  AddHisto("CTag_leptonEtaRel"  ,"CTag_leptonEtaRel", 40, 0 , 0.1  , syst );
  AddHisto("CTag_leptonRatio"  ,"CTag_leptonRatio", 40 , 0 ,1  , syst );
  //With CTag_vertexCategory==0 //
  AddHisto("CTag_jetNSecondaryVertices_Vcat0" ,"CTag_jetNSecondaryVertices_Vcat0",  6,-0.5,5.5  , syst );
  AddHisto("CTag_vertexMass_Vcat0"            ,"CTag_vertexMass_Vcat0",             nSVbins,   0.,8.  , syst );
  AddHisto("CTag_vertexNTracks_Vcat0"         ,"CTag_vertexNTracks_Vcat0",          13, -0.5,12.5, syst );
  AddHisto("CTag_vertexEnergyRatio_Vcat0"     ,"CTag_vertexEnergyRatio_Vcat0",      nSVbins,   0.,1.  , syst );
  AddHisto("CTag_vertexJetDeltaR_Vcat0"       ,"CTag_vertexJetDeltaR_Vcat0",        nSVbins,   0.,0.4 , syst );
  AddHisto("CTag_flightDistance2dSig_Vcat0"   ,"CTag_flightDistance2dSig_Vcat0",    nSVbins,   0.,80. , syst );
  AddHisto("CTag_flightDistance3dSig_Vcat0"   ,"CTag_flightDistance3dSig_Vcat0",    nSVbins,   0.,80. , syst );
  AddHisto("CTag_massVertexEnergyFraction_Vcat0","CTag_massVertexEnergyFraction_Vcat0", 50, 0, 1, syst );
  AddHisto("CTag_vertexBoostOverSqrtJetPt_Vcat0","CTag_vertexBoostOverSqrtJetPt_Vcat0", 50, 0, 1, syst );
  }

 
  Nevent = 0;
  if (fChain == 0) return;
 
  double sumWeightWoPUreweighting       = 0;
  double sumWeightWithPUreweighting     = 0;
  double sumWeightWithPUreweightingUp   = 0;
  double sumWeightWithPUreweightingDown = 0;
 
  TRandom3 *theRand_ = new TRandom3(12345);  // used to apply btag efficiency (btag SF closure test)
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  //------------------------------------------------------------------------------------------------------------------//  
  //----------------------------------------EVENT LOOP ---------------------------------------------------------------// 
  //------------------------------------------------------------------------------------------------------------------//  
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (ientry % 10 == 0) printProgressBar(ientry, nentries);
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //if (ientry > 100000) break;

    //-----------------------------------
    //is data or MC ?
    //-----------------------------------
    if ( datatype == 0) isData=true;
    else                isData=false;
    
    if (isData) N_event_data_before_sel++;
    else        N_event_mc_before_sel++;

    //-----------------------------------
    //initialaze the weight at 1
    //-----------------------------------
    double ww=1; 
         
    vector<TLorentzVector>    theLeptColl;
    vector<Int_t>             theLeptIds;

    for( unsigned short int ilep = 0; ilep < ttbar_nl; ilep++)
    {
      TLorentzVector    theTmpLept;
      theTmpLept.SetPtEtaPhiM(ttbar_lpt[ilep], ttbar_leta[ilep], ttbar_lphi[ilep], ttbar_lm[ilep]);
      theLeptColl.push_back(theTmpLept);
      theLeptIds.push_back(ttbar_lid[ilep]*ttbar_lch[ilep]);
    }
    
    if(ttbar_nl >= 2)
    {
      thettbarselector_.lept1_.SetPtEtaPhiM(ttbar_lpt[0], ttbar_leta[0], ttbar_lphi[0], ttbar_lm[0]);
      thettbarselector_.lept2_.SetPtEtaPhiM(ttbar_lpt[1], ttbar_leta[1], ttbar_lphi[1], ttbar_lm[1]);
    }
 
    vector< pair< TLorentzVector, Float_t> > theJetColl;
    for (unsigned short int ijet = 0; ijet < nJet; ijet++)
    {
      TLorentzVector tmpJet_;
      tmpJet_.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], Jet_mass[ijet]);
      theJetColl.push_back(make_pair( tmpJet_, Jet_genpt[ijet]) );
    }  

    thettbarselector_.met_ = ttbar_metpt;

    //--------------------------------------------//  
    //-------------pile-up reweighting------------//  
    //--------------------------------------------//  
    Float_t puWgtLo(1.0), puWgtNom(1.0), puWgtHi(1.0);                           
    if(!isData)
    {
      if(puWgtGr_)     puWgtNom = puWgtGr_->Eval(nPUtrue);
      if(puWgtDownGr_) puWgtLo  = puWgtDownGr_->Eval(nPUtrue);
      if(puWgtUpGr_)   puWgtHi  = puWgtUpGr_->Eval(nPUtrue);
    }


    // to compute the event weight but no the event selection (for PU reweighting)
    bool computeEvtWgtOnly = true;
    bool computeWeights =  thettbarselector_.passTTbarSelection( isData, theLeptColl, theLeptIds, theJetColl, ttbar_trigWord, ttbar_w, ttbar_nw, wgtcounter, syst, computeEvtWgtOnly);

    // if something goes wrong in the event weight calculation
    if(!computeWeights)  continue;

    // compute the evt weight
    if(!isData) ww = thettbarselector_.evWgt*puWgtNom*x_section[0];

    if(     syst == "PU__minus") ww *= puWgtLo/puWgtNom;
    else if(syst == "PU__plus" ) ww *= puWgtHi/puWgtNom;

    //if(output_name == "output_dy1") cout << "evtwgt = " << ww << endl;

    // To rescale PU reweighting with ratio of sum of weights (at the end)
    if(!isData)
    {
        sumWeightWoPUreweighting        += thettbarselector_.evWgt*x_section[0];
        sumWeightWithPUreweighting      += thettbarselector_.evWgt*puWgtNom*x_section[0];
        sumWeightWithPUreweightingUp    += thettbarselector_.evWgt*puWgtLo*x_section[0];
        sumWeightWithPUreweightingDown  += thettbarselector_.evWgt*puWgtHi*x_section[0];
    }

    computeEvtWgtOnly = false;

    // do the event selection
    bool passsel = false;
    if(produceCTagTree){
    passsel =  thettbarselector_.passSemiLepTTbarSelection( isData, theLeptColl, theLeptIds, theJetColl, ttbar_trigWord, ttbar_w, ttbar_nw, wgtcounter, syst, computeEvtWgtOnly);
    } else{
    passsel = thettbarselector_.passTTbarSelection( isData, theLeptColl, theLeptIds, theJetColl, ttbar_trigWord, ttbar_w, ttbar_nw, wgtcounter, syst, computeEvtWgtOnly);
    }
    if(!passsel)  continue;

    bool fillCommissioningHistograms = false;
    // fill commissioning histograms only if running over nominal samples
    if(syst == "" && thettbarselector_.theSelJetColl.size() >= 2) fillCommissioningHistograms = true;


    //-----------------------------------
    //counter of events
    //-----------------------------------
    Nevent++;
    
    
    //-----------------------------------
    //Fill control plot
    //-----------------------------------
    if(fillCommissioningHistograms)
    {
        if(isData)
        {
                N_event_data_after_sel++;
                nPU_data       ->Fill(nPV);
        }
        else
        {
                N_event_mc_after_sel++;
                nPU_mc         ->Fill(nPUtrue,ww);
                pt_hat         ->Fill(pthat,ww);
                nPV_mc         ->Fill(nPV,ww);
        }
        if(!produceCTagTree){
        FillHistottbar_intFromMap("nPV",                datatype, 0 ,nPV                                ,ww);
        FillHistottbar_intFromMap("nPV_unweighted",     datatype, 0 ,nPV                                 ,1);
        FillHistottbar_floatFromMap("met",              datatype, 0 ,thettbarselector_.met_             ,ww);
        FillHistottbar_floatFromMap("mll",              datatype, 0 ,thettbarselector_.mll_             ,ww);
        FillHistottbar_floatFromMap("pt_e",             datatype, 0 ,thettbarselector_.lept1_.Pt()      ,ww);
        FillHistottbar_floatFromMap("pt_mu",            datatype, 0 ,thettbarselector_.lept2_.Pt()      ,ww);
        }
    }


    int nbjet_ttbar_CSVv2_TWP =0;
    int nbjet_ttbar_CSVv2_MWP =0;
    int nbjet_ttbar_CSVv2_LWP =0;
    int nbjet_ttbar_cMVAv2_TWP =0;
    int nbjet_ttbar_cMVAv2_MWP =0;
    int nbjet_ttbar_cMVAv2_LWP =0;
    float ptjet_ttbar=0;

    int nbjet_afterJetSel_ttbar_CSVv2_TWP =0;
    int nbjet_afterJetSel_ttbar_CSVv2_MWP =0;
    int nbjet_afterJetSel_ttbar_CSVv2_LWP =0;
    int nbjet_afterJetSel_ttbar_cMVAv2_TWP =0;
    int nbjet_afterJetSel_ttbar_cMVAv2_MWP =0;
    int nbjet_afterJetSel_ttbar_cMVAv2_LWP =0;

    int n_ttbar_bgenjet =0;
    int n_ttbar_cgenjet =0;
    int n_ttbar_lgenjet =0;
    int n_ttbar_applySF =0;

    bool isTTbarSelForSF    = false;
    bool isLowerThan60      = false;
    bool isBetween60and120  = false;
    bool isBetween120and320 = false;



    //-----------------------------------
    //Loop on jets 
    //-----------------------------------
    int nJets_pt30=0;
    for (unsigned int ijet = 0; ijet < thettbarselector_.theSelJetColl.size(); ijet++) 
    {

      int newJetIndex = thettbarselector_.theSelJetColl[ijet].second;

      float ptjet   = Jet_pt[newJetIndex];
      float etajet  = Jet_eta[newJetIndex];
      float phijet  = Jet_phi[newJetIndex];
      float ntrkjet  = Jet_ntracks[newJetIndex];  
      int   flav     = Jet_flavour[newJetIndex];

      if( ptjet >= 30 ) nJets_pt30++;
      
      //fill info for ttbar SF
      if(thettbarselector_.theSelJetColl.size() == 2)
      {
        if( thettbarselector_.theSelJetColl[0].first.Pt() >= 30 && thettbarselector_.theSelJetColl[1].first.Pt() >= 30 )
        {
            isTTbarSelForSF = true;

            // pT splitting
            if( thettbarselector_.theSelJetColl[0].first.Pt() < 60 && thettbarselector_.theSelJetColl[1].first.Pt() < 60)        isLowerThan60=true;
            else if( thettbarselector_.theSelJetColl[0].first.Pt() < 120 && thettbarselector_.theSelJetColl[1].first.Pt() < 120) isBetween60and120=true;
            else if( thettbarselector_.theSelJetColl[0].first.Pt() < 320 && thettbarselector_.theSelJetColl[1].first.Pt() < 320) isBetween120and320=true;

            if( !isData)
            {
                 if (fabs(flav) == 5)
                 {
                        n_ttbar_bgenjet++;
                 }
                 else if (fabs(flav) == 4)
                 {
                        n_ttbar_cgenjet++;
                 }
                 else if (fabs(flav) < 4 || fabs(flav) == 21) 
                 {
                        n_ttbar_lgenjet++;
                 }
            }

            if(!isData && fabs(flav) == 5)
            {
                //if( theRand_->Uniform(1.) <= 0.477369)  n_ttbar_applySF++;  //CSVv2T
                if( theRand_->Uniform(1.) <= 0.687265)  n_ttbar_applySF++; //CSVv2M
                //if( theRand_->Uniform(1.) <= 0.684196)  n_ttbar_applySF++; //CSVv2M
            }
            //else if (Jet_CombIVF[newJetIndex] > 0.935) n_ttbar_applySF++;  //CSVv2T
            else if (Jet_CombIVF[newJetIndex] > 0.800) n_ttbar_applySF++;  //CSVv2M

            if (Jet_CombIVF[newJetIndex] > 0.935) nbjet_afterJetSel_ttbar_CSVv2_TWP++; // Tight WP for 76X
            if (Jet_CombIVF[newJetIndex] > 0.800) nbjet_afterJetSel_ttbar_CSVv2_MWP++; // Medium WP for 76X
            if (Jet_CombIVF[newJetIndex] > 0.460) nbjet_afterJetSel_ttbar_CSVv2_LWP++; // Loose WP for 76X
            if (Jet_cMVAv2[newJetIndex] > 0.875) nbjet_afterJetSel_ttbar_cMVAv2_TWP++; // Tight WP for 76X
            if (Jet_cMVAv2[newJetIndex] > 0.185) nbjet_afterJetSel_ttbar_cMVAv2_MWP++; // Medium WP for 76X
            if (Jet_cMVAv2[newJetIndex] > -0.715) nbjet_afterJetSel_ttbar_cMVAv2_LWP++; // Loose WP for 76X

        }
      }

      if(!produceCTagTree){
      if (nTrack>1000) cout << " data nTrack " << nTrack << endl;
      }

      if (ptjet>ptjet_ttbar) ptjet_ttbar=ptjet;

      int n_sv             = 0.;
      float mass_sv        =-1.;
      float chi2norm_sv    =-1.;
      float flightSig_sv   =-1.;    
      float flight2DSig_sv =-1.;    
      float sv_dR_jet      =-1.;
      float sv_dR_dir_sum  =-1.; 
      float sv_dR_jet_sum  =-1.;
      float sv_en_rat      =-1.; 
      float sv_pt	   =-1.;      
      float sveta         =-1000.; 
      float svphi         =-1000.; 
      float sv_flight3D    =-1.;
      float sv_flight3Derr =-1.;
      float sv_flight2D    =-1.;
      float sv_flight2Derr =-1.;
      int   sv_totchar     =-1.;
      float sv_nTrk        = 0.;
      
      
      float tche     = Jet_Ip2P[newJetIndex];
      float tchp     = Jet_Ip3P[newJetIndex];
      float jetproba = Jet_ProbaP[newJetIndex];
      float jetbproba= Jet_BprobP[newJetIndex];
      float ssvhe    = Jet_Svx[newJetIndex] ;
      float ssvhp    = Jet_SvxHP[newJetIndex];
      float csv      = Jet_CombSvx[newJetIndex];
      float csv_v2   = Jet_CombIVF[newJetIndex];
      float cmva_v2  = Jet_cMVAv2[newJetIndex];
      float cvsB     = CTag_Jet_CvsB[newJetIndex]; 
      float cvsL     = CTag_Jet_CvsL[newJetIndex];

      float CvsB     = CTag_Jet_CvsB[newJetIndex];
      float CvsBN    = CTag_Jet_CvsBN[newJetIndex];
      float CvsBP    = CTag_Jet_CvsBP[newJetIndex];
      float CvsL     = CTag_Jet_CvsL[newJetIndex];
      float CvsLN    = CTag_Jet_CvsLN[newJetIndex];
      float CvsLP    = CTag_Jet_CvsLP[newJetIndex]; 

      bool isPU=false;
   
      if( fillCommissioningHistograms)
      {
        if(!isData && Jet_genpt[newJetIndex] <= 8) isPU = true;

        //fill jet multiplicity
        if (!isData) 
        {
                if      (fabs(flav)==5)                  njet_b++;
                else if (fabs(flav)==4)                  njet_c++;
	        else if (fabs(flav)<4 || fabs(flav)==21) njet_l++;
                njet_mc++;
                jet_pt_mc   ->Fill(ptjet,ww);
        }
        else njet_data++;
      
        if (csv_v2>0.935)   nbjet_ttbar_CSVv2_TWP++; // Tight WP for 76X
        if (csv_v2>0.800)   nbjet_ttbar_CSVv2_MWP++; // Medium WP for 76X
        if (csv_v2>0.460)   nbjet_ttbar_CSVv2_LWP++; // Loose WP for 76X
        if (cmva_v2>0.875)  nbjet_ttbar_cMVAv2_TWP++; // Tight WP for 76X
        if (cmva_v2>0.185)  nbjet_ttbar_cMVAv2_MWP++; // Medium WP for 76X
        if (cmva_v2>-0.715) nbjet_ttbar_cMVAv2_LWP++; // Loose WP for 76X

        //if(output_name == "output_dy1") cout << "csvv2 = " << csv_v2 << "ptjet= " << ptjet << endl;
      
        FillHisto_floatFromMap("jet_multi",                  flav, isPU ,thettbarselector_.theSelJetColl.size() , ww);
        FillHisto_floatFromMap("jet_pt_all",                 flav, isPU ,ptjet                                  , ww);
        FillHisto_floatFromMap("genjet_pt_all",              flav, isPU ,Jet_genpt[newJetIndex]                 , ww);
        
        if(isPU) FillHisto_floatFromMap("CSVv2_pu",          flav, isPU ,csv_v2                                 , ww);
        if(!produceCTagTree){
        if (nSV > 0)FillHisto_floatFromMap("jet_pt_sv",      flav, isPU ,ptjet                                  , ww);
        }
        FillHisto_floatFromMap("jet_eta",     flav, isPU ,etajet   , ww);
        FillHisto_floatFromMap("jet_phi",     flav, isPU ,phijet   , ww);
        FillHisto_intFromMap(  "track_multi", flav, isPU ,ntrkjet  , ww);

      }
      
      int n1_ip=-1;
      int n2_ip=-1;
      int n3_ip=-1;
      float sig1_ip=-9999;
      float sig2_ip=-9999;
      float sig3_ip=-9999;
      float sig12D_ip=-9999;
      float sig22D_ip=-9999;
      float sig32D_ip=-9999;   
         
      int ntracksel  =0;   

         
      if ( fillCommissioningHistograms && produceJetProbaTree ) 
      {
	
	for (int itrk=Jet_nFirstTrack[newJetIndex]; itrk<Jet_nLastTrack[newJetIndex] ; itrk++)
        {

	  //-------------------------//
	  //-----Track selection-----//
          //-------------------------//  
	  passNhit=false;
	  passPix= false;
	  passIPz=false;
	  passPt=false;
	  passnormchi2=false;
	  passtrkdist=false;
	  passtrklen=false;
	  passTrackIP2D=false;
	    
	  //if (Track_nHitAll[itrk]>=8)           passNhit=true;
	  //if (Track_nHitPixel[itrk]>=2)         passPix= true;
	  if (Track_nHitAll[itrk]>=0)           passNhit=true;     // HIP mitigation
	  if (Track_nHitPixel[itrk]>=1)         passPix= true;     // HIP mitigation
	  if (fabs(Track_dz[itrk])<17)          passIPz=true;
	  if (Track_pt[itrk]>1)                 passPt=true;
	  if (Track_chi2[itrk]<5)               passnormchi2=true;
	  if (fabs(Track_dist[itrk])<0.07)      passtrkdist=true;
	  if (Track_length[itrk]<5)             passtrklen=true;
	  if (fabs(Track_IP2D[itrk])<0.2)       passTrackIP2D=true;
	  
	  if (!use_selected_tracks)
          {
	    
	    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passTrackIP2D)
            {
	      FillHisto_floatFromMap("track_len_cut",   flav, isPU ,Track_length[itrk]          , ww);
	    }
	    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrklen && passTrackIP2D)
            {
	      FillHisto_floatFromMap("track_dist_cut",  flav, isPU ,fabs(Track_dist[itrk])      , ww);
	    }	    
	    if (passNhit && passPix && passIPz && passPt && passtrkdist && passtrklen && passTrackIP2D)
            {
	      FillHisto_floatFromMap("track_chi2_cut",  flav, isPU ,Track_chi2[itrk]	        , ww);
	    }	    
	    if (passNhit && passPix && passIPz && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D)
            {
	      FillHisto_floatFromMap("track_pt_cut",    flav, isPU ,Track_pt[itrk]              , ww);
	    }	    
	    if (passNhit && passPix && passPt && passnormchi2 && passtrkdist && passtrklen){
            
	      FillHisto_floatFromMap("track_dz_cut",    flav, isPU ,Track_dz[itrk]              , ww);
	    }
	    if (passNhit && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D)
            {
	      FillHisto_intFromMap(  "track_HPix_cut",  flav, isPU ,Track_nHitPixel[itrk]       , ww);
	    }	
	    if (passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D)
            {
	      FillHisto_intFromMap(  "track_nHit_cut",  flav, isPU ,Track_nHitAll[itrk]         , ww);	  
	    }
	    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen )
            {
	      FillHisto_intFromMap(  "track_IP2D_cut",  flav, isPU ,Track_IP2D[itrk]            , ww);	  
	    }	    
	  }
	  
          if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D)
          {
	    ntracksel++;
	    FillHisto_floatFromMap("track_chi2",    flav, isPU ,Track_chi2[itrk]	,ww);
	    FillHisto_intFromMap(  "track_nHit",    flav, isPU ,Track_nHitAll[itrk]     ,ww);
	    FillHisto_intFromMap(  "track_HPix",    flav, isPU ,Track_nHitPixel[itrk]   ,ww);
	    FillHisto_floatFromMap("track_IPs",     flav, isPU ,Track_IPsig[itrk]       ,ww);
	    FillHisto_floatFromMap("track_IP",      flav, isPU ,Track_IP[itrk]          ,ww);
	    FillHisto_floatFromMap("track_IP2Ds",   flav, isPU ,Track_IP2Dsig[itrk]     ,ww);
	    FillHisto_floatFromMap("track_IP2D",    flav, isPU ,Track_IP2D[itrk]        ,ww);
	    FillHisto_floatFromMap("track_IP2Derr", flav, isPU ,Track_IP2Derr[itrk]     ,ww);	  
	    FillHisto_floatFromMap("track_IPerr",   flav, isPU ,Track_IPerr[itrk]       ,ww);	  	  
	    FillHisto_floatFromMap("track_dz",      flav, isPU ,Track_dz[itrk]          ,ww);	  
	    FillHisto_intFromMap(  "track_isfromSV",flav, isPU ,Track_isfromSV[itrk]    ,ww);	  
	    FillHisto_floatFromMap("track_len",     flav, isPU ,Track_length[itrk]      ,ww);
	    FillHisto_floatFromMap("track_dist",    flav, isPU ,fabs(Track_dist[itrk])  ,ww);
	    FillHisto_floatFromMap("track_pt",      flav, isPU ,Track_pt[itrk]          ,ww);	  
	    
	  
            //Tracks sorted by IP
	    Float_t sig  =Track_IP[itrk]/Track_IPerr[itrk];
	    Float_t sig2D=Track_IP2D[itrk]/Track_IP2Derr[itrk];
            if (sig>sig1_ip) 
            {
              sig3_ip=sig2_ip;
	      sig2_ip=sig1_ip;
	      sig1_ip=sig;
              sig32D_ip=sig22D_ip;
	      sig22D_ip=sig12D_ip;
	      sig12D_ip=sig2D;	      
	      n3_ip=n2_ip;
	      n2_ip=n1_ip;
	      n1_ip=itrk;
	    }
	    else if (sig>sig2_ip)
            {
	      sig3_ip=sig2_ip;
	      sig2_ip=sig;	      
	      sig32D_ip=sig22D_ip;
	      sig22D_ip=sig2D;
	      n3_ip=n2_ip;
	      n2_ip=itrk;
	    }
	    else if (sig>sig3_ip) 
            {
	      sig3_ip=sig;
	      sig32D_ip=sig2D;
	      n3_ip=itrk;
            }	      

	  }//end selected tracks

	}//end tracks loop
	

	if (n1_ip>-1)  
        {
	  FillHisto_floatFromMap("track_IPs1tr",    flav, isPU ,sig1_ip               , ww);
	  FillHisto_floatFromMap("track_IP1tr",     flav, isPU ,Track_IP[n1_ip]       , ww);
	  FillHisto_floatFromMap("track_IPerr1tr",  flav, isPU ,Track_IPerr[n1_ip]    , ww);
	  FillHisto_floatFromMap("track_IP2Ds1tr",  flav, isPU ,sig12D_ip             , ww);
	  FillHisto_floatFromMap("track_IP2D1tr",   flav, isPU ,Track_IP2D[n1_ip]     , ww);
	  FillHisto_floatFromMap("track_IP2Derr1tr",flav, isPU ,Track_IP2Derr[n1_ip]  , ww);		  
	}

	if (n2_ip>-1) 
        {
	  FillHisto_floatFromMap("track_IPs2tr",    flav, isPU ,sig2_ip               , ww);
	  FillHisto_floatFromMap("track_IP2tr",     flav, isPU ,Track_IP[n2_ip]       , ww);
	  FillHisto_floatFromMap("track_IPerr2tr",  flav, isPU ,Track_IPerr[n2_ip]    , ww);	
	  FillHisto_floatFromMap("track_IP2Ds2tr",  flav, isPU ,sig22D_ip             , ww);
	  FillHisto_floatFromMap("track_IP2D2tr",   flav, isPU ,Track_IP2D[n2_ip]     , ww);
	  FillHisto_floatFromMap("track_IP2Derr2tr",flav, isPU ,Track_IP2Derr[n2_ip]  , ww);
	}
                
	if (n3_ip>-1) 
        {
	  FillHisto_floatFromMap("track_IPs3tr",    flav, isPU ,sig3_ip               , ww);
	  FillHisto_floatFromMap("track_IP3tr",     flav, isPU ,Track_IP[n3_ip]       , ww);
	  FillHisto_floatFromMap("track_IPerr3tr",  flav, isPU ,Track_IPerr[n3_ip]    , ww);	
	  FillHisto_floatFromMap("track_IP2Ds3tr",  flav, isPU ,sig32D_ip             , ww);
	  FillHisto_floatFromMap("track_IP2D3tr",   flav, isPU ,Track_IP2D[n3_ip]     , ww);
	  FillHisto_floatFromMap("track_IP2Derr3tr",flav, isPU ,Track_IP2Derr[n3_ip]  , ww);
	}

	FillHisto_intFromMap(        "trk_multi_sel",     flav, isPU ,ntracksel	         , ww);  
	FillHisto2D_int_floatFromMap("seltrack_vs_jetpt", flav, isPU ,ptjet ,  ntracksel , ww);


        FillHisto_floatFromMap("tagvarCSV_vertexCategory",      flav, isPU, TagVarCSV_vertexCategory[newJetIndex],              ww);
        FillHisto_floatFromMap("tagvarCSV_Sig2dAboveCharm",     flav, isPU, TagVarCSV_trackSip2dSigAboveCharm[newJetIndex],     ww);
        FillHisto_floatFromMap("tagvarCSV_trackSumJetEtRatio",  flav, isPU, TagVarCSV_trackSumJetEtRatio[newJetIndex],          ww);
        FillHisto_floatFromMap("tagvarCSV_trackSumJetDeltaR",   flav, isPU, TagVarCSV_trackSumJetDeltaR[newJetIndex],           ww);

        for (int inrel=Jet_nFirstTrkEtaRelTagVarCSV[newJetIndex]; inrel<Jet_nLastTrkEtaRelTagVarCSV[newJetIndex]; inrel++) 
        {
           FillHisto_floatFromMap("tagvarCSV_trackEtaRel",      flav, isPU, TagVarCSV_trackEtaRel[inrel],                       ww);
        }
        
        FillHisto_floatFromMap("tagvarCSV_energyratio",         flav, isPU, TagVarCSV_vertexEnergyRatio[newJetIndex],           ww);

        for (int inrel=Jet_nFirstTrkTagVarCSV[newJetIndex]; inrel<Jet_nLastTrkTagVarCSV[newJetIndex]; inrel++) 
        {
           FillHisto_floatFromMap("tagvarCSV_trackSip3dSig",    flav, isPU, TagVarCSV_trackSip3dSig[inrel],                     ww);
        }
        if (TagVarCSV_vertexCategory[newJetIndex]==0) 
        {
          FillHisto_floatFromMap("tagvarCSV_vertexmass_cat0",   flav, isPU, TagVarCSV_vertexMass[newJetIndex],                  ww);
          if (TagVarCSV_vertexNTracks[newJetIndex]>=3) 
          {
                FillHisto_floatFromMap("tagvarCSV_vertexmass3trk_cat0",     flav, isPU, TagVarCSV_vertexMass[newJetIndex],      ww);
          }
          FillHisto_floatFromMap("tagvarCSV_vertexNTracks_cat0",   flav, isPU, TagVarCSV_vertexNTracks[newJetIndex],            ww);
          FillHisto_floatFromMap("tagvarCSV_2DsigFlightDist_cat0", flav, isPU, TagVarCSV_flightDistance2dSig[newJetIndex],      ww);
          FillHisto_floatFromMap("tagvarCSV_vertexJetDeltaR_cat0", flav, isPU, TagVarCSV_vertexJetDeltaR[newJetIndex],          ww);
        }

	
        //---------------------------------
	//fill information related to SV
	//---------------------------------
	n_sv           = Jet_SV_multi[newJetIndex];	  
	FillHisto_intFromMap(  "sv_multi_0",      flav, isPU ,n_sv      ,ww);

	if (n_sv>0)
        {  
	  chi2norm_sv    = SV_chi2[Jet_nFirstSV[newJetIndex]]/SV_ndf[Jet_nFirstSV[newJetIndex]];
	  flightSig_sv   = SV_flight[Jet_nFirstSV[newJetIndex]]/SV_flightErr[Jet_nFirstSV[newJetIndex]];
	  flight2DSig_sv = SV_flight2D[Jet_nFirstSV[newJetIndex]]/SV_flight2DErr[Jet_nFirstSV[newJetIndex]];
          mass_sv        = SV_mass[Jet_nFirstSV[newJetIndex]];
	  sv_dR_jet      = SV_deltaR_jet[Jet_nFirstSV[newJetIndex]];
	  sv_dR_dir_sum  = SV_deltaR_sum_dir[Jet_nFirstSV[newJetIndex]];
	  sv_dR_jet_sum  = SV_deltaR_sum_jet[Jet_nFirstSV[newJetIndex]];
	  sv_en_rat      = SV_EnergyRatio[Jet_nFirstSV[newJetIndex]];
	  sv_pt          = SV_vtx_pt[Jet_nFirstSV[newJetIndex]];
	  sveta          = SV_vtx_eta[Jet_nFirstSV[newJetIndex]];
	  svphi          = SV_vtx_phi[Jet_nFirstSV[newJetIndex]];
	  sv_flight3D    = SV_flight[Jet_nFirstSV[newJetIndex]] ;  
          sv_flight3Derr = SV_flightErr[Jet_nFirstSV[newJetIndex]];
	  sv_flight2D    = SV_flight2D[Jet_nFirstSV[newJetIndex]];    
	  sv_flight2Derr = SV_flight2DErr[Jet_nFirstSV[newJetIndex]];    
	  sv_totchar     = SV_totCharge[Jet_nFirstSV[newJetIndex]] ;
	  sv_nTrk        = SV_nTrk[Jet_nFirstSV[newJetIndex]] ;  

          //-------------------------//
	  //-----SV histograms-------//
          //-------------------------//   
	  FillHisto_intFromMap(  "sv_multi",        flav, isPU ,n_sv               , ww);
	  FillHisto_floatFromMap("sv_chi2norm",     flav, isPU ,chi2norm_sv        , ww);
	  FillHisto_floatFromMap("sv_mass",         flav, isPU ,mass_sv            , ww);
	  FillHisto_floatFromMap("sv_deltaR_jet",   flav, isPU ,sv_dR_jet          , ww);
	  FillHisto_floatFromMap("sv_deltaR_sumJet",flav, isPU ,sv_dR_dir_sum      , ww);
	  FillHisto_floatFromMap("sv_deltaR_sumDir",flav, isPU ,sv_dR_jet_sum      , ww);
	  FillHisto_floatFromMap("sv_en_ratio",     flav, isPU ,sv_en_rat          , ww);
	  FillHisto_floatFromMap("sv_pt",           flav, isPU ,sv_pt              , ww);
	  FillHisto_floatFromMap("sv_flight2D",     flav, isPU ,sv_flight2D        , ww);
	  FillHisto_floatFromMap("sv_flight2Derr",  flav, isPU ,sv_flight2Derr     , ww);
	  FillHisto_floatFromMap("sv_flightSig2D",  flav, isPU ,flight2DSig_sv     , ww);
	  FillHisto_intFromMap("sv_tot_charge",     flav, isPU ,sv_totchar         , ww);
	  FillHisto_intFromMap(  "svnTrk",          flav, isPU ,sv_nTrk            , ww);
	  FillHisto_floatFromMap("sv_eta",          flav, isPU ,sveta              , ww);
	  FillHisto_floatFromMap("sv_phi",          flav, isPU ,svphi              , ww);	
	  FillHisto_floatFromMap("sv_flight3D",     flav, isPU ,sv_flight3D        , ww);
	  FillHisto_floatFromMap("sv_flight3Derr",  flav, isPU ,sv_flight3Derr     , ww);
	  FillHisto_floatFromMap("sv_flight3DSig",  flav, isPU ,flightSig_sv       , ww);
	
	  if (sv_nTrk >2)FillHisto_floatFromMap("sv_mass_3trk", flav, isPU ,mass_sv, ww);
	
	  FillHisto2D_float_floatFromMap("sv_mass_vs_flightDist3D"     ,flav,isPU ,sv_flight3D,mass_sv  , ww);	
	  FillHisto2D_float_floatFromMap("avg_sv_mass_vs_jetpt"        ,flav,isPU ,ptjet,mass_sv        , ww);
	  FillHisto2D_float_floatFromMap("sv_deltar_jet_vs_jetpt"      ,flav,isPU ,ptjet,sv_dR_jet      , ww);
	  FillHisto2D_float_floatFromMap("sv_deltar_sum_jet_vs_jetpt"  ,flav,isPU ,ptjet,sv_dR_dir_sum  , ww);
	  FillHisto2D_float_floatFromMap("sv_deltar_sum_dir_vs_jetpt"  ,flav,isPU ,ptjet,sv_dR_jet_sum  , ww);	    
	    
	}//end n_sv > 0 condition

      }//end produce jetProbaTree and fillCommissioningHistograms
	    

      if( fillCommissioningHistograms)
      {

        //Taggers
        FillHisto_floatFromMap("TCHE",  flav, isPU, tche	, ww);
        FillHisto_floatFromMap("TCHP",  flav, isPU, tchp	, ww);
        FillHisto_floatFromMap("JP",    flav, isPU, jetproba    , ww);
        FillHisto_floatFromMap("JBP",   flav, isPU, jetbproba   , ww);
        FillHisto_floatFromMap("SSV",   flav, isPU, ssvhe	, ww);
        FillHisto_floatFromMap("SSVHP", flav, isPU, ssvhp	, ww);
        FillHisto_floatFromMap("CSV",   flav, isPU, csv	        , ww);
        FillHisto_floatFromMap("CSVv2", flav, isPU, csv_v2      , ww);
        FillHisto_floatFromMap("cMVAv2",flav, isPU, cmva_v2     , ww);
        FillHisto_floatFromMap("CvsB",  flav, isPU, cvsB        , ww);
        FillHisto_floatFromMap("CvsL",  flav, isPU, cvsL        , ww);
      }

      if (fillCommissioningHistograms && produceNewAlgoTree) 
      {
        float softmu            = Jet_SoftMu[newJetIndex]  ;
        float solfel            = Jet_SoftEl[newJetIndex];

        FillHisto_floatFromMap("SoftMu",      flav, isPU, softmu        , ww);
        FillHisto_floatFromMap("SoftEl",      flav, isPU, solfel        , ww);

      }

      if( fillCommissioningHistograms)
      {
        //Taggers
        FillHisto_floatFromMap("TCHE_extended1",  flav, isPU, tche  , ww);
        FillHisto_floatFromMap("TCHP_extended1",  flav, isPU, tchp  , ww);
        FillHisto_floatFromMap("TCHE_extended2",  flav, isPU, tche  , ww);
        FillHisto_floatFromMap("TCHP_extended2",  flav, isPU, tchp  , ww);
        FillHisto_floatFromMap("discri_ssche0",   flav, isPU, ssvhe , ww);
        FillHisto_floatFromMap("discri_sschp0",   flav, isPU, ssvhp , ww);
      
      } // end of fillCommissioningHistograms     


      if (fillCommissioningHistograms && produceNewAlgoTree) 
      {
        // PFMuon
        int npfmu=0;
        int indpfmu=-1;

        TLorentzVector thejet;
        thejet.SetPtEtaPhiM(ptjet, etajet, phijet, 0);

        float minpf=0;
        for (int im=0; im<nPFMuon; im++) 
        {
          TLorentzVector thepfmu;
          thepfmu.SetPtEtaPhiM(PFMuon_pt[im],PFMuon_eta[im],PFMuon_phi[im],0);
          float drpfj=thepfmu.DeltaR(thejet);
          float diffdr = drpfj-PFMuon_deltaR[im];
          if (diffdr<0.) diffdr*=-1.;
          if (drpfj< 0.5 && diffdr<0.05 && PFMuon_GoodQuality[im]>0) 
          {
               if (PFMuon_pt[im]> minpf)
               {
                 indpfmu=im;
                 minpf=PFMuon_pt[im];
               }
               npfmu++;
          }
        }


        FillHisto_intFromMap(  "pfmuon_multi",  flav, isPU , npfmu   ,ww);
        
        if (indpfmu>-1)
        {
          FillHisto_intFromMap(  "pfmuon_goodquality", flav, isPU, PFMuon_GoodQuality[indpfmu]       , ww);
          FillHisto_floatFromMap("pfmuon_pt",          flav, isPU, PFMuon_pt[indpfmu]                , ww);
          FillHisto_floatFromMap("pfmuon_eta",         flav, isPU, PFMuon_eta[indpfmu]               , ww);
          FillHisto_floatFromMap("pfmuon_phi",         flav, isPU, PFMuon_phi[indpfmu]               , ww);
          FillHisto_floatFromMap("pfmuon_Sip",         flav, isPU, PFMuon_IPsig[indpfmu]             , ww);
          FillHisto_floatFromMap("pfmuon_ptrel",       flav, isPU, PFMuon_ptrel[indpfmu]             , ww);
          FillHisto_floatFromMap("pfmuon_ratio",       flav, isPU, PFMuon_ratio[indpfmu]             , ww);
          FillHisto_floatFromMap("pfmuon_ratiorel",    flav, isPU, PFMuon_ratioRel[indpfmu]          , ww);
          FillHisto_floatFromMap("pfmuon_deltar",      flav, isPU, PFMuon_deltaR[indpfmu]            , ww);
        }

        //PFElectron
        int npfel=0;
        int indpfel=-1;

        minpf=0;
        for (int im=0; im<nPFElectron; im++) 
        {
          TLorentzVector thepfel;
          thepfel.SetPtEtaPhiM(PFElectron_pt[im],PFElectron_eta[im],PFElectron_phi[im],0);
          float drpfj=thepfel.DeltaR(thejet);
          float diffdr = drpfj-PFElectron_deltaR[im];
          if (diffdr<0.) diffdr*=-1.;
          if (drpfj< 0.5 && diffdr<0.05 && PFElectron_pt[im]>2.)
          {
               if (PFElectron_pt[im]> minpf)
               {
                 indpfel=im;
                 minpf=PFElectron_pt[im];
               }
               npfel++;
          }
        }


        FillHisto_intFromMap(  "pfelectron_multi",  flav, isPU , npfel   ,ww);

        if (indpfel>-1)
        {
          FillHisto_floatFromMap("pfelectron_pt",        flav, isPU, PFElectron_pt[indpfel]           , ww);
          FillHisto_floatFromMap("pfelectron_eta",       flav, isPU, PFElectron_eta[indpfel]          , ww);
          FillHisto_floatFromMap("pfelectron_phi",       flav, isPU, PFElectron_phi[indpfel]          , ww);
          FillHisto_floatFromMap("pfelectron_ptrel",     flav, isPU, PFElectron_ptrel[indpfel]        , ww);
          FillHisto_floatFromMap("pfelectron_ratio",     flav, isPU, PFElectron_ratio[indpfel]        , ww);
          FillHisto_floatFromMap("pfelectron_ratiorel",  flav, isPU, PFElectron_ratioRel[indpfel]     , ww);
          FillHisto_floatFromMap("pfelectron_deltar",    flav, isPU, PFElectron_deltaR[indpfel]       , ww);
        } 


      }// end produceNewAlgoTree and fillCommissioningHistograms

      if ( fillCommissioningHistograms && produceCTagTree ){
        FillHisto_floatFromMap("CTag_tagvarCSV_vertexCategory",      flav, isPU, TagVarCSV_vertexCategory[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_tagvarCSV_Sig2dAboveCharm",     flav, isPU, TagVarCSV_trackSip2dSigAboveCharm[newJetIndex],     ww);
        FillHisto_floatFromMap("CTag_tagvarCSV_trackSumJetEtRatio",  flav, isPU, TagVarCSV_trackSumJetEtRatio[newJetIndex],          ww);
        FillHisto_floatFromMap("CTag_tagvarCSV_trackSumJetDeltaR",   flav, isPU, TagVarCSV_trackSumJetDeltaR[newJetIndex],           ww);

        for (int inrel=Jet_nFirstTrkEtaRelTagVarCSV[newJetIndex]; inrel<Jet_nLastTrkEtaRelTagVarCSV[newJetIndex]; inrel++)
        {
           FillHisto_floatFromMap("CTag_tagvarCSV_trackEtaRel",      flav, isPU, TagVarCSV_trackEtaRel[inrel],                       ww);
        }

        FillHisto_floatFromMap("CTag_tagvarCSV_energyratio",         flav, isPU, TagVarCSV_vertexEnergyRatio[newJetIndex],           ww);

        for (int inrel=Jet_nFirstTrkTagVarCSV[newJetIndex]; inrel<Jet_nLastTrkTagVarCSV[newJetIndex]; inrel++)
        {
           FillHisto_floatFromMap("CTag_tagvarCSV_trackSip3dSig",    flav, isPU, TagVarCSV_trackSip3dSig[inrel],                     ww);
        }
        if (TagVarCSV_vertexCategory[newJetIndex]==0)
        {
          FillHisto_floatFromMap("CTag_tagvarCSV_vertexmass_cat0",   flav, isPU, TagVarCSV_vertexMass[newJetIndex],                  ww);
          if (TagVarCSV_vertexNTracks[newJetIndex]>=3)
          {
                FillHisto_floatFromMap("CTag_tagvarCSV_vertexmass3trk_cat0",     flav, isPU, TagVarCSV_vertexMass[newJetIndex],      ww);
          }
          FillHisto_floatFromMap("CTag_tagvarCSV_vertexNTracks_cat0",   flav, isPU, TagVarCSV_vertexNTracks[newJetIndex],            ww);
          FillHisto_floatFromMap("CTag_tagvarCSV_2DsigFlightDist_cat0", flav, isPU, TagVarCSV_flightDistance2dSig[newJetIndex],      ww);
          FillHisto_floatFromMap("CTag_tagvarCSV_vertexJetDeltaR_cat0", flav, isPU, TagVarCSV_vertexJetDeltaR[newJetIndex],          ww);
        }
        
        FillHisto_floatFromMap("CTag_jetNTracks",               flav, isPU ,CTag_jetNTracks[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_jetNTracksEtaRel",         flav, isPU ,CTag_jetNTracksEtaRel[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_jetNLeptons",              flav, isPU ,CTag_jetNLeptons[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_trackSumJetEtRatio",       flav, isPU ,CTag_trackSumJetEtRatio[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_trackSumJetDeltaR",        flav, isPU ,CTag_trackSumJetDeltaR[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_trackSip2dSigAboveCharm",  flav, isPU ,CTag_trackSip2dSigAboveCharm[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_trackSip3dSigAboveCharm",  flav, isPU ,CTag_trackSip3dSigAboveCharm[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexCategory",           flav, isPU ,CTag_vertexCategory[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_jetNSecondaryVertices",    flav, isPU ,CTag_jetNSecondaryVertices[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexMass",               flav, isPU ,CTag_vertexMass[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexNTracks",            flav, isPU ,CTag_vertexNTracks[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexEnergyRatio",        flav, isPU ,CTag_vertexEnergyRatio[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexJetDeltaR",          flav, isPU ,CTag_vertexJetDeltaR[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_flightDistance2dSig",      flav, isPU ,CTag_flightDistance2dSig[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_flightDistance3dSig",      flav, isPU ,CTag_flightDistance3dSig[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_massVertexEnergyFraction", flav, isPU ,CTag_massVertexEnergyFraction[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexBoostOverSqrtJetPt", flav, isPU ,CTag_vertexBoostOverSqrtJetPt[newJetIndex],              ww);
        FillHisto_floatFromMap("CTag_vertexLeptonCategory",     flav, isPU ,CTag_vertexLeptonCategory[newJetIndex],              ww);

        if (CTag_vertexCategory[newJetIndex]==0){
          FillHisto_floatFromMap("CTag_jetNSecondaryVertices_Vcat0",    flav, isPU ,CTag_jetNSecondaryVertices[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_vertexMass_Vcat0",               flav, isPU ,CTag_vertexMass[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_vertexNTracks_Vcat0",            flav, isPU ,CTag_vertexNTracks[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_vertexEnergyRatio_Vcat0",        flav, isPU ,CTag_vertexEnergyRatio[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_vertexJetDeltaR_Vcat0",          flav, isPU ,CTag_vertexJetDeltaR[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_flightDistance2dSig_Vcat0",      flav, isPU ,CTag_flightDistance2dSig[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_flightDistance3dSig_Vcat0",      flav, isPU ,CTag_flightDistance3dSig[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_massVertexEnergyFraction_Vcat0", flav, isPU ,CTag_massVertexEnergyFraction[newJetIndex],              ww);
          FillHisto_floatFromMap("CTag_vertexBoostOverSqrtJetPt_Vcat0", flav, isPU ,CTag_vertexBoostOverSqrtJetPt[newJetIndex],              ww);
        }
         
        for (int inrel=Jet_nFirstTrkCTagVar[newJetIndex]; inrel<Jet_nLastTrkCTagVar[newJetIndex]; inrel++)
        {
           FillHisto_floatFromMap("CTag_trackPtRel",    flav, isPU, CTag_trackPtRel[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackPPar",    flav, isPU, CTag_trackPPar[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackDeltaR",    flav, isPU, CTag_trackDeltaR[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackPtRatio",    flav, isPU, CTag_trackPtRatio[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackPParRatio",    flav, isPU, CTag_trackPParRatio[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackSip2dSig",    flav, isPU, CTag_trackSip2dSig[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackSip3dSig",    flav, isPU, CTag_trackSip3dSig[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackDecayLenVal",    flav, isPU, CTag_trackDecayLenVal[inrel],                     ww);
           FillHisto_floatFromMap("CTag_trackJetDistVal",    flav, isPU, CTag_trackJetDistVal[inrel],                     ww);
        }
        for (int inrel=Jet_nFirstTrkEtaRelCTagVar[newJetIndex]; inrel<Jet_nLastTrkEtaRelCTagVar[newJetIndex]; inrel++)
        {
           FillHisto_floatFromMap("CTag_trackEtaRel",      flav, isPU, CTag_trackEtaRel[inrel],                       ww);
        }
        for (int inrel=Jet_nFirstLepCTagVar[newJetIndex]; inrel<Jet_nLastLepCTagVar[newJetIndex]; inrel++)
        {
           FillHisto_floatFromMap("CTag_leptonPtRel",      flav, isPU, CTag_leptonPtRel[inrel],                       ww);
           FillHisto_floatFromMap("CTag_leptonSip3d",      flav, isPU, CTag_leptonSip3d[inrel],                       ww);
           FillHisto_floatFromMap("CTag_leptonDeltaR",      flav, isPU, CTag_leptonDeltaR[inrel],                       ww);
           FillHisto_floatFromMap("CTag_leptonRatioRel",      flav, isPU, CTag_leptonRatioRel[inrel],                       ww);
           FillHisto_floatFromMap("CTag_leptonEtaRel",      flav, isPU, CTag_leptonEtaRel[inrel],                       ww);
           FillHisto_floatFromMap("CTag_leptonRatio",      flav, isPU, CTag_leptonRatio[inrel],                       ww);
        }

        FillHisto_floatFromMap("JP",    flav, isPU, jetproba    , ww);
        FillHisto_floatFromMap("CSVv2", flav, isPU, csv_v2      , ww);
        FillHisto_floatFromMap("cMVAv2",flav, isPU, cmva_v2     , ww);

        FillHisto_floatFromMap("CvsB",  flav, isPU, CvsB        ,ww);
        FillHisto_floatFromMap("CvsBN", flav, isPU, CvsBN       ,ww);
        FillHisto_floatFromMap("CvsBP", flav, isPU, CvsBP       ,ww);
        FillHisto_floatFromMap("CvsL",  flav, isPU, CvsL        ,ww);
        FillHisto_floatFromMap("CvsLN", flav, isPU, CvsLN       ,ww);
        FillHisto_floatFromMap("CvsLP", flav, isPU, CvsLP       ,ww);   

   }// End of fillCommissioningHistograms && produceCTagTree      

    } // End Loop on Jets

    if(!produceCTagTree){
    if( nJets_pt30 >= 2 ) FillHistottbar_intFromMap("njet_pt30", datatype, 0 , nJets_pt30  , ww);
 
    if(fillCommissioningHistograms)
    {
        FillHistottbar_intFromMap("njet",               datatype, 0 ,thettbarselector_.theSelJetColl.size()  , ww);
        FillHistottbar_intFromMap("nbtag_CSVv2T",       datatype, 0 ,nbjet_ttbar_CSVv2_TWP                   , ww);
        FillHistottbar_intFromMap("nbtag_CSVv2M",       datatype, 0 ,nbjet_ttbar_CSVv2_MWP                   , ww);
        FillHistottbar_intFromMap("nbtag_CSVv2L",       datatype, 0 ,nbjet_ttbar_CSVv2_LWP                   , ww);
        FillHistottbar_intFromMap("nbtag_cMVAv2T",      datatype, 0 ,nbjet_ttbar_cMVAv2_TWP                  , ww);
        FillHistottbar_intFromMap("nbtag_cMVAv2M",      datatype, 0 ,nbjet_ttbar_cMVAv2_MWP                  , ww);
        FillHistottbar_intFromMap("nbtag_cMVAv2L",      datatype, 0 ,nbjet_ttbar_cMVAv2_LWP                  , ww);
        FillHistottbar_floatFromMap("pt_jet",           datatype, 0 ,ptjet_ttbar                             , ww);

        // HIP check (as function of run range for Run2016B)
        if( isData && thettbarselector_.theSelJetColl.size() == 2 )//ALPHA
        {
            int cat=-1;
            if( Run < 273450 )      cat=0;
            else if( Run < 273730 ) cat=1;
            else if( Run < 274240 ) cat=2;
            else if( Run < 274284 ) cat=3;
            else if( Run < 274335 ) cat=4;
            else if( Run < 274382 ) cat=5;
            else if( Run < 274421 ) cat=6;
            else if( Run < 274440 ) cat=7;
            else if( Run < 274968 ) cat=8;
            else if( Run < 274970 ) cat=9;
            else if( Run < 275000 ) cat=10;
            else if( Run < 275068 ) cat=11;
            else if( Run < 275124 ) cat=12;
            else if( Run < 275292 ) cat=13;
            else if( Run < 275311 ) cat=14;
            else if( Run < 275345 ) cat=15;
            else if( Run < 275376 ) cat=16;
            else if( Run < 275657 ) cat=17;
            else                    cat=18; // fill events for later run ranges

            FillHistottbar_intFromMap("nEvt_run", datatype, 0, cat, 1);
            if ( nbjet_ttbar_CSVv2_LWP >= 2 ) FillHistottbar_intFromMap("nEvt_run_CSVv2L", datatype, 0, cat, 1); 
            if ( nbjet_ttbar_CSVv2_MWP >= 2 ) FillHistottbar_intFromMap("nEvt_run_CSVv2M", datatype, 0, cat, 1);
            if ( nbjet_ttbar_CSVv2_TWP >= 2 ) FillHistottbar_intFromMap("nEvt_run_CSVv2T", datatype, 0, cat, 1);
        }
    }

    TString tmp_ptbin = "";
    
    if(isTTbarSelForSF == true)
    {

        if(isLowerThan60)           tmp_ptbin = "_Inf60";        
        else if(isBetween60and120)  tmp_ptbin = "_60-120";        
        else if(isBetween120and320) tmp_ptbin = "_120-320";        

        FillHistottbar_intFromMap("nbtag_all_afterJetSel_CSVv2T",               datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP  ,ww);
        FillHistottbar_intFromMap("nbtag_all_afterJetSel_CSVv2M_SFapplied",     datatype, 0 ,n_ttbar_applySF                    ,ww);
        FillHistottbar_intFromMap("nbtag_all_afterJetSel_CSVv2M",               datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP  ,ww);
        FillHistottbar_intFromMap("nbtag_all_afterJetSel_CSVv2L",               datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP  ,ww);

        FillHistottbar_intFromMap("nbtag_all_afterJetSel_cMVAv2T",              datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP ,ww);
        FillHistottbar_intFromMap("nbtag_all_afterJetSel_cMVAv2M",              datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP ,ww);
        FillHistottbar_intFromMap("nbtag_all_afterJetSel_cMVAv2L",              datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP ,ww);

        if(tmp_ptbin != "")
        {

                FillHistottbar_intFromMap("nbtag_all"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_all"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_all"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_all"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_all"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_all"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

        }


        if( n_ttbar_bgenjet == 2)
        {

                FillHistottbar_intFromMap("nbtag_2b_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2b_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2b_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_2b_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2b_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2b_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                if(tmp_ptbin != "")
                {

                        FillHistottbar_intFromMap("nbtag_2b"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2b"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2b"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                        FillHistottbar_intFromMap("nbtag_2b"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2b"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2b"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                }
        }
        else if( n_ttbar_bgenjet == 1 && n_ttbar_cgenjet == 1 )
        {

                FillHistottbar_intFromMap("nbtag_1b1c_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1c_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1c_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_1b1c_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1c_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1c_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                if(tmp_ptbin != "")
                {

                        FillHistottbar_intFromMap("nbtag_1b1c"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1c"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1c"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                        FillHistottbar_intFromMap("nbtag_1b1c"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1c"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1c"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                }
        }
        else if( n_ttbar_bgenjet == 1 && n_ttbar_lgenjet == 1 )
        {

                FillHistottbar_intFromMap("nbtag_1b1l_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1l_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1l_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_1b1l_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1l_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1b1l_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                if(tmp_ptbin != "")
                {

                        FillHistottbar_intFromMap("nbtag_1b1l"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1l"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1l"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                        FillHistottbar_intFromMap("nbtag_1b1l"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1l"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1b1l"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                }
        }
        else if( n_ttbar_cgenjet == 2)
        {

                FillHistottbar_intFromMap("nbtag_2c_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2c_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2c_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_2c_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2c_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2c_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                if(tmp_ptbin != "")
                {

                        FillHistottbar_intFromMap("nbtag_2c"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2c"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2c"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                        FillHistottbar_intFromMap("nbtag_2c"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2c"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2c"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                }
        }
        else if( n_ttbar_cgenjet == 1 && n_ttbar_lgenjet == 1 )
        {

                FillHistottbar_intFromMap("nbtag_1c1l_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1c1l_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1c1l_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_1c1l_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1c1l_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_1c1l_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                if(tmp_ptbin != "")
                {

                        FillHistottbar_intFromMap("nbtag_1c1l"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1c1l"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1c1l"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                        FillHistottbar_intFromMap("nbtag_1c1l"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1c1l"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_1c1l"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                }
        }
        else if( n_ttbar_lgenjet == 2 )
        {

                FillHistottbar_intFromMap("nbtag_2l_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2l_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2l_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                FillHistottbar_intFromMap("nbtag_2l_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2l_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                FillHistottbar_intFromMap("nbtag_2l_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                if(tmp_ptbin != "")
                {

                        FillHistottbar_intFromMap("nbtag_2l"+tmp_ptbin+"_afterJetSel_CSVv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2l"+tmp_ptbin+"_afterJetSel_CSVv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2l"+tmp_ptbin+"_afterJetSel_CSVv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_CSVv2_LWP     ,ww);

                        FillHistottbar_intFromMap("nbtag_2l"+tmp_ptbin+"_afterJetSel_cMVAv2T",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_TWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2l"+tmp_ptbin+"_afterJetSel_cMVAv2M",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_MWP     ,ww);
                        FillHistottbar_intFromMap("nbtag_2l"+tmp_ptbin+"_afterJetSel_cMVAv2L",     datatype, 0 ,nbjet_afterJetSel_ttbar_cMVAv2_LWP     ,ww);

                }
        }
        else if( !isData ) cout << "There is an issue with the gen jets: see line: " << __LINE__ << " and n_ttbar[b,c,l]genjet = " << n_ttbar_bgenjet << " | " << n_ttbar_cgenjet << " | " << n_ttbar_lgenjet << endl;

    } // end of isTTbarSelForSF condition
   
    }//end of !produceCTagTree

    //cout << "---------------------End event---------------------------" << endl;
    
  }//----------------------------------------End loop events ----------------------------------------------------------------------------------------------//  
  
  cout << endl;
  cout << "----------------Configuration-------------------" <<endl;
  cout << " Run over "           << N_event_data_before_sel << " data events and " << N_event_mc_before_sel << " mc events."    << endl;
  cout << " After selection -> " << N_event_data_after_sel  << " data events and " << N_event_mc_after_sel  << "mc events left" << endl;   
  cout << endl;

  myfile->cd();

  // scale PU reweighting with ratio of sum of weights
  Double_t puScale = sumWeightWoPUreweighting/sumWeightWithPUreweighting;
  if(!isData) cout << "sumWeightWoPUreweighting = " << sumWeightWoPUreweighting << " | sumWeightWithPUreweighting = " << sumWeightWithPUreweighting << "  | scale (PUreweighting) = " << puScale << endl;

  for (unsigned int i=0; i<HistoBtag.size(); i++) 
  {
    if(!isData && sumWeightWithPUreweighting != 0) HistoBtag[i]->Scale(puScale);
    HistoBtag[i]->Write();
  }
  for (unsigned int i=0; i<HistoTTbar.size(); i++) 
  {
    if(!isData && sumWeightWithPUreweighting != 0) HistoTTbar[i]->Scale(puScale);
    HistoTTbar[i]->Write();
  }
  for (unsigned int i=0; i<HistoBtag2D.size(); i++) 
  {
    if(!isData && sumWeightWithPUreweighting != 0) HistoBtag2D[i]->Scale(puScale);
    HistoBtag2D[i]->Write();
  }  
  myfile->Close();

}

void CommPlotProducer4ttbar::AddHisto(TString name, TString title, int nbins, float min, float max, TString syst)
{

  TString syst_tmp = "";
  if(syst != "") syst_tmp = "_"+syst;
 
  TH1D* h_b      = new TH1D(name+"_b"+syst_tmp,         title+"_b"+syst_tmp,            nbins,min,max);
  TH1D* h_pu     = new TH1D(name+"_pu"+syst_tmp,        title+"_pu"+syst_tmp,           nbins,min,max);
  //TH1D* h_bfromg = new TH1D(name+"_bfromg"+syst_tmp,    title+"_bfromg"+syst_tmp,       nbins,min,max);  
  TH1D* h_c      = new TH1D(name+"_c"+syst_tmp,         title+"_c"+syst_tmp,            nbins,min,max);  
  TH1D* h_l      = new TH1D(name+"_l"+syst_tmp,         title+"_l"+syst_tmp,            nbins,min,max);
  TH1D* h_data   = new TH1D(name+"_data"+syst_tmp,      title+"_data"+syst_tmp,         nbins,min,max);
  
  h_b        ->Sumw2();
  h_pu       ->Sumw2(); 
  //h_bfromg   ->Sumw2();  
  h_c        ->Sumw2();  
  h_l        ->Sumw2(); 
  h_data     ->Sumw2();
  
  HistoBtag.push_back(h_b);
  HistoBtag.push_back(h_pu);  
  //HistoBtag.push_back(h_bfromg);  
  HistoBtag.push_back(h_c);  
  HistoBtag.push_back(h_l);  
  HistoBtag.push_back(h_data);  
  HistoBtag_map[name.Data()] = numb_histo;
  
  numb_histo++;
  
}

void CommPlotProducer4ttbar::AddHistottbar(TString name, TString title, int nbins, float min, float max, TString syst)
{
 
  TString syst_tmp = "";
  if(syst != "") syst_tmp = "_"+syst;
 
  TH1D* h_ttbar   = new TH1D(name+"_ttbar"+syst_tmp,    title+"_ttbar"+syst_tmp,        nbins,min,max);
  TH1D* h_dy      = new TH1D(name+"_dy"+syst_tmp,       title+"_dy"+syst_tmp,           nbins,min,max);  
  TH1D* h_st      = new TH1D(name+"_st"+syst_tmp,       title+"_st"+syst_tmp,           nbins,min,max);  
  TH1D* h_ww      = new TH1D(name+"_ww"+syst_tmp,       title+"_ww"+syst_tmp,           nbins,min,max);  
  TH1D* h_wz      = new TH1D(name+"_wz"+syst_tmp,       title+"_wz"+syst_tmp,           nbins,min,max);  
  TH1D* h_zz      = new TH1D(name+"_zz"+syst_tmp,       title+"_zz"+syst_tmp,           nbins,min,max);  
  TH1D* h_data    = new TH1D(name+"_data"+syst_tmp,     title+"_data"+syst_tmp,         nbins,min,max);
  
  
  h_ttbar     ->Sumw2();
  h_dy        ->Sumw2();  
  h_st        ->Sumw2(); 
  h_ww        ->Sumw2();  
  h_wz        ->Sumw2();  
  h_zz        ->Sumw2();  
  h_data      ->Sumw2();
  
  HistoTTbar.push_back(h_ttbar);
  HistoTTbar.push_back(h_dy);  
  HistoTTbar.push_back(h_st);  
  HistoTTbar.push_back(h_ww);  
  HistoTTbar.push_back(h_wz);  
  HistoTTbar.push_back(h_zz);  
  HistoTTbar.push_back(h_data);  
  HistoTTbar_map[name.Data()] = numb_histo2;
  
  numb_histo2++;
  
}

void CommPlotProducer4ttbar::FillHisto_float(int flavour, bool isPU, int number, float value, double weight)
{
  
  if (!isData)
  {
    if (isPU)                                       HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==5)                      HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
    
  }  
  else                                              HistoBtag[number*5 +4]->Fill(value);
  
  
}
void CommPlotProducer4ttbar::FillHisto_int(int flavour, bool isPU, int number, int value, double weight)
{
  
  if (!isData)
  {
    if (isPU)                                       HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==5)                      HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
  }  
  else                                              HistoBtag[number*5 +4]->Fill(value);
  
}


void CommPlotProducer4ttbar::FillHisto_floatFromMap(TString name, int flavour, bool isPU, float value, double weight)
{
  
  int number = HistoBtag_map[name.Data()] ;
  if (!isData)
  {
    if (isPU)                                       HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==5)                      HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
    
  }  
  else                                              HistoBtag[number*5 +4]->Fill(value);
   
}


void CommPlotProducer4ttbar::FillHisto_intFromMap(TString name, int flavour, bool isPU, int value, double weight)
{
  
  int number = HistoBtag_map[name.Data()] ;
  if (!isData)
  {
    if (isPU)                                       HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==5)                      HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
  }  
  else                                              HistoBtag[number*5 +4]->Fill(value);
  
}

void CommPlotProducer4ttbar::FillHistottbar_floatFromMap(TString name, int flavour, bool isPU, float value, double weight)
{
  
  int number = HistoTTbar_map[name.Data()] ;
  if (!isData){
    if (flavour==1)                  HistoTTbar[number*7 +0]->Fill(value,weight);
    else if (flavour==2)             HistoTTbar[number*7 +1]->Fill(value,weight);
    else if (flavour==3)             HistoTTbar[number*7 +2]->Fill(value,weight); 
    else if (flavour==4)             HistoTTbar[number*7 +3]->Fill(value,weight); 
    else if (flavour==5)             HistoTTbar[number*7 +4]->Fill(value,weight); 
    else if (flavour==6)             HistoTTbar[number*7 +5]->Fill(value,weight); 
  }  
  else                               HistoTTbar[number*7 +6]->Fill(value);
  
   
}


void CommPlotProducer4ttbar::FillHistottbar_intFromMap(TString name, int flavour, bool isPU, int value, double weight)
{
  
  int number = HistoTTbar_map[name.Data()] ;
  if (!isData)
  {
    if (flavour==1)                  HistoTTbar[number*7 +0]->Fill(value,weight);
    else if (flavour==2)             HistoTTbar[number*7 +1]->Fill(value,weight);
    else if (flavour==3)             HistoTTbar[number*7 +2]->Fill(value,weight); 
    else if (flavour==4)             HistoTTbar[number*7 +3]->Fill(value,weight); 
    else if (flavour==5)             HistoTTbar[number*7 +4]->Fill(value,weight); 
    else if (flavour==6)             HistoTTbar[number*7 +5]->Fill(value,weight); 
  }  
  else                               HistoTTbar[number*7 +6]->Fill(value);
  
}


//-----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------2D PLOTS-----------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
void CommPlotProducer4ttbar::AddHisto2D(TString name, TString title, int nbins, float min, float max, int nbins2, float
min2, float max2, TString syst)
{
 
  TString syst_tmp = "";
  if(syst != "") syst_tmp = "_"+syst;
 
  TH2D* h_b      = new TH2D(name+"_b"+syst_tmp,         title+"_b"+syst_tmp,            nbins,min,max,nbins2,min2,max2);
  TH2D* h_pu     = new TH2D(name+"_pu"+syst_tmp,        title+"_pu"+syst_tmp,           nbins,min,max,nbins2,min2,max2);
  //TH2D* h_bfromg = new TH2D(name+"_bfromg"+syst_tmp,    title+"_bfromg"+syst_tmp,       nbins,min,max,nbins2,min2,max2);  
  TH2D* h_c      = new TH2D(name+"_c"+syst_tmp,         title+"_c"+syst_tmp,            nbins,min,max,nbins2,min2,max2);  
  TH2D* h_l      = new TH2D(name+"_l"+syst_tmp,         title+"_l"+syst_tmp,            nbins,min,max,nbins2,min2,max2);
  TH2D* h_data   = new TH2D(name+"_data"+syst_tmp,      title+"_data"+syst_tmp,         nbins,min,max,nbins2,min2,max2);
  
  
  h_b        ->Sumw2();
  h_pu       ->Sumw2(); 
  //h_bfromg   ->Sumw2();  
  h_c        ->Sumw2();  
  h_l        ->Sumw2(); 
  h_data     ->Sumw2();
  
  HistoBtag2D.push_back(h_b);
  HistoBtag2D.push_back(h_pu);  
  //HistoBtag2D.push_back(h_bfromg);  
  HistoBtag2D.push_back(h_c);  
  HistoBtag2D.push_back(h_l);  
  HistoBtag2D.push_back(h_data);  
  HistoBtag2D_map[name.Data()] = numb_histo2D;
  numb_histo2D++;
  
}


void CommPlotProducer4ttbar::FillHisto2D_int_floatFromMap(TString name, int flavour, bool isPU, int value, float value2, double weight)
{
   
  int number = HistoBtag2D_map[name.Data()] ;
  if (!isData)
  {
    if (isPU)                                       HistoBtag2D[number*5 +1]->Fill(value,value2,weight);
    else if (fabs(flavour)==5)                      HistoBtag2D[number*5 +0]->Fill(value,value2,weight);
    else if (fabs(flavour)==4)                      HistoBtag2D[number*5 +2]->Fill(value,value2,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*5 +3]->Fill(value,value2,weight);
    
  }  
  else                                              HistoBtag2D[number*5 +4]->Fill(value,value2);
  
}

void CommPlotProducer4ttbar::FillHisto2D_float_floatFromMap(TString name, int flavour, bool isPU, float value, float value2, double weight)
{
  
  int number = HistoBtag2D_map[name.Data()] ;
  if (!isData){
    if (isPU)                                       HistoBtag2D[number*5 +1]->Fill(value,value2,weight);
    else if (fabs(flavour)==5)                      HistoBtag2D[number*5 +0]->Fill(value,value2,weight);
    else if (fabs(flavour)==4)                      HistoBtag2D[number*5 +2]->Fill(value,value2,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*5 +3]->Fill(value,value2,weight);
    
  }  
  else                                              HistoBtag2D[number*5 +4]->Fill(value,value2);
  
}


#define CommPlotProducer_cxx
#include "CommPlotProducer.h"

#include <TH2.h>
#include <TStyle.h>
#include "TH1F.h"
#include "TH2F.h"
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
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>

using namespace std;


//---------------------------------------------------------------------------------------//
//-------------Fill number of events for each sample--------------------------------------//
//----------------------------------------------------------------------------------------//

void CommPlotProducer::Fill_nevent(double n15,double n20,double n30,double n50,double n80,double n120,double n170,double n300,double
n470,double n600){

 
  nmc_evt_vect[0]=n15;
  nmc_evt_vect[1]=n20;
  nmc_evt_vect[2]=n30;  
  nmc_evt_vect[3]=n50;  
  nmc_evt_vect[4]=n80;  
  nmc_evt_vect[5]=n120;
  nmc_evt_vect[6]=n170;  
  nmc_evt_vect[7]=n300;  
  nmc_evt_vect[8]=n470;   
  nmc_evt_vect[9]=n600;
  
}

//---------------------------------------------------------------------------------------//
//-------------Use it if you don't know how many events per sample you have------------//
//----------------------------------------------------------------------------------------//

void CommPlotProducer::Counter(){    

  int   Nevent = 0;
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    
     
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    Nevent++;
    
    if(Nevent%50000 ==0 && Nevent!=0) cout << " number of processed events is " << Nevent << endl;    
    if      ( pthat >=  15. && pthat <  30. ){
      use15_30=true;
      n15_30++;
    } 
    else if ( pthat >=  30. && pthat <  50. ){
      use30_50=true;    
      n30_50++;
    } 
    else if ( pthat >=  50. && pthat <  80. ){
      use50_80=true;    
      n50_80++;
    }
    else if ( pthat >=  80. && pthat < 120. ){
      use80_120=true;    
      n80_120++;
    }
    else if ( pthat >= 120. && pthat < 170. ){
      use120_170=true;    
      n120_170++;
    }
    else if ( pthat >= 170. && pthat < 300. ){
      use170_300=true;    
      n170_300++;
    }
    else if ( pthat >= 300. && pthat < 470. ){
      use300_470=true;    
      n300_470++;
    }
    else if ( pthat >= 470. && pthat < 600. ){
      use470_600=true;    
      n470_600++;
    }
    else if ( pthat >= 600. && pthat <= 800. ){
      use600_800=true;    
      n600_800++;
    }

  }  
  
  if  (use15_30)  cout << "Run with QCD15_30   sample with " << n15_30  << " events" <<endl; 
  if  (use30_50)  cout << "Run with QCD30_50   sample with " << n30_50  << " events" <<endl;
  if  (use50_80)  cout << "Run with QCD50_80   sample with " << n50_80  << " events" <<endl;
  if  (use80_120) cout << "Run with QCD80_120  sample with " << n80_120 << " events" <<endl;
  if  (use120_170)cout << "Run with QCD120_170 sample with " << n120_170<< " events" <<endl;
  if  (use170_300)cout << "Run with QCD170_300 sample with " << n170_300<< " events" <<endl;
  if  (use300_470)cout << "Run with QCD300_470 sample with " << n300_470<< " events" <<endl;
  if  (use470_600)cout << "Run with QCD470_600 sample with " << n470_600<< " events" <<endl;
  if  (use600_800)cout << "Run with QCD600_800 sample with " << n600_800<< " events" <<endl;
}

//------------------------------------------------------------------------------//
//-------------Set PU--------------------------------------------------------//
//-----------------------------------------------------------------------------//
float CommPlotProducer::SetPU(vector<float> PUvector, TString PUdataFile){


  vector<float> mc_vect;
  vector<float> data_vect;
  
  
  TFile *filepuest = new TFile(PUdataFile,"READ");
  TH1F* npu_dat= (TH1F*) filepuest->Get("pileup");
    
  for( int i=0; i<60; ++i) {
    mc_vect.push_back(PUvector[i]);
  }
  for (int i=0; i<PUvector.size(); i++){
    data_vect.push_back(npu_dat->GetBinContent(i+1));
  }
  
  LumiWeights = reweight::LumiReWeighting(mc_vect,data_vect);
  
  WeightPU = LumiWeights.weight( nPUtrue );   
  
  return WeightPU;
}

float CommPlotProducer::SetPU2012_S7(TString PUdataFile){
 
  // Distribution used for Summer2012 MC.
  Double_t Summer2012_S7[60] = {2.344E-05,2.344E-05,2.344E-05,2.344E-05,4.687E-04,4.687E-04,7.032E-04,9.414E-04,1.234E-03,
  1.603E-03,2.464E-03,3.250E-03,5.021E-03,6.644E-03,8.502E-03,1.121E-02,1.518E-02,2.033E-02,2.608E-02,3.171E-02,
  3.667E-02,4.060E-02,4.338E-02,4.520E-02,4.641E-02,4.735E-02,4.816E-02,4.881E-02,4.917E-02,4.909E-02,4.842E-02,
  4.707E-02,4.501E-02,4.228E-02,3.896E-02,3.521E-02,3.118E-02,2.702E-02,2.287E-02,1.885E-02,1.508E-02,1.166E-02,
  8.673E-03,6.190E-03,4.222E-03,2.746E-03,1.698E-03,9.971E-04,5.549E-04,2.924E-04,1.457E-04,6.864E-05,3.054E-05,
  1.282E-05,5.081E-06,1.898E-06,6.688E-07,2.221E-07,6.947E-08,2.047E-08 };  

  vector<float> mc_vect;
  vector<float> data_vect;
  
  
  TFile *filepuest = new TFile(PUdataFile,"READ");
  TH1F* npu_dat= (TH1F*) filepuest->Get("pileup");
    
  for( int i=0; i<60; ++i) {
    mc_vect.push_back(Summer2012_S7[i]);
  }
  for (int i=0; i<60; i++){
    data_vect.push_back(npu_dat->GetBinContent(i+1));
  }
  
  LumiWeights = reweight::LumiReWeighting(mc_vect,data_vect);
  
  WeightPU = LumiWeights.weight( nPUtrue );   
  
  return WeightPU;


}

float CommPlotProducer::SetPU2012_S10(TString PUdataFile){
 
  // Distribution used for Summer2012 MC.
  Double_t Summer2012_S10[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 
  6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 
  5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 
  3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 
  8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 
  7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 
  2.322E-05, 1.570E-05, 5.005E-06};  

  vector<float> mc_vect;
  vector<float> data_vect;
  
  
  TFile *filepuest = new TFile(PUdataFile,"READ");
  TH1F* npu_dat= (TH1F*) filepuest->Get("pileup");
    
  for( int i=0; i<60; ++i) {
    mc_vect.push_back(Summer2012_S10[i]);
  }
  for (int i=0; i<60; i++){
    data_vect.push_back(npu_dat->GetBinContent(i+1));
  }
  
  LumiWeights = reweight::LumiReWeighting(mc_vect,data_vect);
  
  WeightPU = LumiWeights.weight( nPUtrue );   
  
  return WeightPU;


}
//------------------------------------------------------------------------------//
//-------------Set cross sections--------------------------------------------------------//
//-----------------------------------------------------------------------------//
int CommPlotProducer::SetXS(){

  double pythia_xs[9]={816000000.0,66285328.0, 8148778.0,1033680.0,156293.3,34138.15,1759.549,113.8791,27.01};

  for (int i=0; i<9; i++){
    x_section[i]=pythia_xs[i];
  }
  choice=0; 
  return choice;

}

int CommPlotProducer::SetXS(TString generator, bool MuEnriched, int TeV){

double pythia_xs8MU[8]={806298.0, 176187.6, 40448,7463.94,2299.752,151.8048, 11.79648, 2.690196};
double pythia_xs7MU[7]={1471168.0, 1224034.0,578463.0,144421.74,29048.7,4440.2215,2837.6712 };
double pythia_xs7  [10]={815900000.0,53120000.0,6359000.0,784300.0,115100.0,24260.0,1168.0,70.22,15.55,1.844 };
double herwig_xs8  [6]={5.32E7,6550000.0,836000.0,126600.0,27900.0,1579.0};

choice=0;
 
  if (generator=="pythia"){
    if (MuEnriched){
      if (TeV==8){
        choice=1;
        for (int i=0; i<8; i++){
          x_section[i]=pythia_xs8MU[i];
        }        
      }
      else {
        choice=2;
        for (int i=0; i<7; i++){
          x_section[i]=pythia_xs7MU[i];
        }       
      }
    }
    else{
      if (TeV==7){
        choice=3;
        for (int i=0; i<10; i++){
          x_section[i]=pythia_xs7[i];
        }         
      }
    }
  }
  if (generator=="herwig"){  
    choice=4;
    for (int i=0; i<6; i++){
      x_section[i]=herwig_xs8[i];
    }   
  }
  

  
  return choice;

}

void CommPlotProducer::SetSumXS(){
  sum_xs=0.0;
  
  if (nmc_evt_vect[0]>0)  sum_xs+=x_section[0];
  if (nmc_evt_vect[1]>0)  sum_xs+=x_section[1];
  if (nmc_evt_vect[2]>0)  sum_xs+=x_section[2];
  if (nmc_evt_vect[3]>0)  sum_xs+=x_section[3];
  if (nmc_evt_vect[4]>0) sum_xs+=x_section[4];
  if (nmc_evt_vect[5]>0) sum_xs+=x_section[5];
  if (nmc_evt_vect[6]>0) sum_xs+=x_section[6];
  if (nmc_evt_vect[7]>0) sum_xs+=x_section[7];
  if (nmc_evt_vect[8]>0) sum_xs+=x_section[8]; 
  if (nmc_evt_vect[9]>0) sum_xs+=x_section[9];

}


//-------------------------------------------------------------------------------------------//
//-------------Compute the event weight given the pthat-------------------------------------------//
//------------------------------------------------------------------------------------------//

float CommPlotProducer::GetEvtWeight(){

  WeightXS=0.0;
  float  nevt  =0.0;
  float  xs    =0.0;
  
  if (choice ==0){// Use pythia 8 TeV
  
    if ( pthat >=  15. && pthat <  30. ){
      nevt  =nmc_evt_vect[0];
      xs=x_section[0];
    } 
     
    if ( pthat >=  30. && pthat <  50. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    }
      
    if ( pthat >=  50. && pthat <  80. ){
      nevt  =nmc_evt_vect[2];
      xs=x_section[2];
    }     
    if ( pthat >=  80. && pthat <  120. ){
      nevt  =nmc_evt_vect[3];
      xs=x_section[3];
    } 
     
    if ( pthat >=  120. && pthat <  170. ){
      nevt  =nmc_evt_vect[4];
      xs=x_section[4];
    }
      
    if ( pthat >=  170. && pthat <  300. ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    }    
  
    if ( pthat >=  300. && pthat <  470. ){
      nevt  =nmc_evt_vect[6];
      xs=x_section[6];
    } 
     
    if ( pthat >=  470. && pthat <  600. ){
      nevt  =nmc_evt_vect[7];
      xs=x_section[7];
    }
      
    if ( pthat >=  600. && pthat <  800. ){
      nevt  =nmc_evt_vect[8];
      xs=x_section[8];
    }
  
  }
  if (choice ==1){// Use pythia 8 TeV MuEnriched

     
    if ( pthat >=  30. && pthat <  50. ){
      nevt  =nmc_evt_vect[0];
      xs=x_section[0];
    }
      
    if ( pthat >=  50. && pthat <  80. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    }     
    if ( pthat >=  80. && pthat <  120. ){
      nevt  =nmc_evt_vect[2];
      xs=x_section[2];
    } 
     
    if ( pthat >=  120. && pthat <  170. ){
      nevt  =nmc_evt_vect[3];
      xs=x_section[3];
    }
      
    if ( pthat >=  170. && pthat <  300. ){
      nevt  =nmc_evt_vect[4];
      xs=x_section[4];
    }    
  
    if ( pthat >=  300. && pthat <  470. ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    } 
     
    if ( pthat >=  470. && pthat <  600. ){
      nevt  =nmc_evt_vect[6];
      xs=x_section[6];
    }
      
    if ( pthat >=  600. && pthat <  800. ){
      nevt  =nmc_evt_vect[7];
      xs=x_section[7];
    }
  
  }  
  
  if (choice ==3){// Use pythia 7 TeV MuEnriched

    if ( pthat >=  15. && pthat <  30. ){
      nevt  =nmc_evt_vect[0];
      xs=x_section[0];
    } 
     
    if ( pthat >=  30. && pthat <  50. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    }
      
    if ( pthat >=  50. && pthat <  80. ){
      nevt  =nmc_evt_vect[2];
      xs=x_section[2];
    }     
    if ( pthat >=  80. && pthat <  120. ){
      nevt  =nmc_evt_vect[3];
      xs=x_section[3];
    } 
     
    if ( pthat >=  120. && pthat <  170. ){
      nevt  =nmc_evt_vect[4];
      xs=x_section[4];
    }
      
    if ( pthat >=  170. && pthat <  300. ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    }    
  
    if ( pthat >=  300. && pthat <  470. ){
      nevt  =nmc_evt_vect[6];
      xs=x_section[6];
    } 
     
    if ( pthat >=  470. && pthat <  600. ){
      nevt  =nmc_evt_vect[7];
      xs=x_section[7];
    }
      
    if ( pthat >=  600. && pthat <  800. ){
      nevt  =nmc_evt_vect[8];
      xs=x_section[8];
    }
    if ( pthat >=  800. && pthat <  1000. ){
      nevt  =nmc_evt_vect[9];
      xs=x_section[9];
    }  
  }
  
  if (choice ==2){// Use pythia 7 TeV MuEnriched

     
    if ( pthat >=  15. && pthat <  20. ){
      nevt  =nmc_evt_vect[0];
      xs=x_section[0];
    }
      
    if ( pthat >=  20. && pthat <  30. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    }     
    if ( pthat >=  30. && pthat < 50. ){
      nevt  =nmc_evt_vect[2];
      xs=x_section[2];
    } 
     
    if ( pthat >=  50. && pthat <  80. ){
      nevt  =nmc_evt_vect[3];
      xs=x_section[3];
    }
      
    if ( pthat >=  80. && pthat <  120. ){
      nevt  =nmc_evt_vect[4];
      xs=x_section[4];
    }    
  
    if ( pthat >=  120. && pthat <  150. ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    } 
    if ( pthat >=  150. ){
      nevt  =nmc_evt_vect[6];
      xs=x_section[6];
    }   
  }     
  
  if (choice ==4){// Use Herwig 8 TeV

     
    if ( pthat >=  30 && pthat <  50. ){
      nevt  =nmc_evt_vect[0];
      xs=x_section[0];
    }
      
    if ( pthat >=  50. && pthat <  80. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    }     
    if ( pthat >=  80. && pthat < 120. ){
      nevt  =nmc_evt_vect[2];
      xs=x_section[2];
    } 
     
    if ( pthat >=  120. && pthat <  170. ){
      nevt  =nmc_evt_vect[3];
      xs=x_section[3];
    }
      
    if ( pthat >=  170. && pthat <  300. ){
      nevt  =nmc_evt_vect[4];
      xs=x_section[4];
    }    
    if ( pthat >=  300 ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    }
   
  }  
  
  WeightXS =xs/(sum_xs*nevt);  
  
  if ( pthat <1. ) WeightXS=1. ;

  return  WeightXS;
   
  
}


void CommPlotProducer::Loop(int trigger, float PtMin_Cut, float PtMax_Cut, TString output_name)
  
  
{ 
  
  //---------------Configuration-----------------------------------------// 
  //produceJetProbaTree=false;
  int   IntCut = trigger;
  float PtMin  = PtMin_Cut;  
  float PtMax  = PtMax_Cut;  
  float EtaCut = 2.4; 
  
  int Year = 2012;  
  
  
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
  int njet_mc    =0; 
  int njet_data  =0; 

  
  bool passNhit;
  bool passPix    ; 
  bool passIPz    ;
  bool passPt     ; 
  bool passnormchi2; 
  bool passtrkdist ; 
  bool passtrklen  ;
  bool passTrackIP2D;
  
  //---------------------------------------------------------------------//
  TFile *myfile=new TFile(output_name+".root",      "recreate");
  
  // --------------------------------------Histograms declaration------------------------------------------------//   
  TH1F* nPU_mc                  = new TH1F("nPU_mc",                "nPU_mc",                50,-0.5,49.5);
  TH1F* nPU_data                = new TH1F("nPU_data",              "nPU_data",              50,-0.5,49.5);
  TH1F* nPV_mc                  = new TH1F("nPV_mc",                "nPV_mc",                50,-0.5,49.5);
  TH1F* pt_hat                  = new TH1F("pt_hat",                "pt_hat",                80,0,800);
  TH1F* jet_pt_mc               = new TH1F("jet_pt_mc",  	    "jet_pt_mc", 	     80,0,PtMax);
  
  // --------------------------------------Histograms declaration -----------------------------------------//
  
  AddHisto("nPU"          ,"number of PU event",	     60,-0.5,59.5);
 
  AddHisto("jet_multi"    ,"number of jets",		     20,0,20    );
  AddHisto("jet_pt_all"	  ,"pT of all jets",		     PtMax/10,0,PtMax );
  AddHisto("jet_pt_sv"	  ,"pT of jets containing a SV",     PtMax/10,0,PtMax);
  AddHisto("jet_eta"	  ,"eta of all jets",		     50,-2.5,2.5);
  AddHisto("jet_phi"	  ,"phi of all jets",		     20,-5,5    );
  
  AddHisto("muon_multi"   ,      "number of muons",	   7,-0.5,6.5    );
  AddHisto("muon_multi_sel"   ,  "number of selected muons",7,-0.5,6.5   );
  AddHisto("mu_ptrel"	  ,      "pT rel. of the muon",	   50,0,5        );
  AddHisto("mu_chi2"	  ,      "norm. chi2 of the muon", 50,0,10       );  
  AddHisto("muon_Pt",		 "Muon p_{T}",  	   200, 0, 100	 );
  AddHisto("muon_eta",		 "Muon #eta",  	           50, -2.5, 2.5 );  
  AddHisto("muon_phi",		 "Muon #phi",  	           80, -4, 4	 );
  AddHisto("muon_Ip3d", 	 "Muon 3D IP",  	   50, -0.1, 0.1 );
  AddHisto("muon_Ip2d", 	 "Muon 2D IP",  	   50, -0.1, 0.1 );
  AddHisto("muon_Sip3d",	 "Muon 3D IP significance",50, -35, 35   );
  AddHisto("muon_Sip2d",	 "Muon 2D IP significance",50, -35, 35   );
  AddHisto("muon_DeltaR",	 "Muon1 deltaR",50,0,0.5); //90
  
  AddHisto("sv_deltaR_sumJet",   "SVvtxSumJetDeltaR",                                   50,0.,0.5    );
  AddHisto("sv_deltaR_sumDir",   "SVvtxSumVtxDirDeltaR",                                50,0.,0.5    );
  AddHisto("sv_en_ratio",        "Fractional energy",                                   50,0.,1.     );  
  AddHisto("sv_aboveC",          "IP significance 2D charm",                            50,-35.,35.  );
  AddHisto("sv_pt",              "Vtx p_{T}",                                           50,0.,100.   );
  AddHisto("sv_eta",             "Vtx #eta",                                            50, -2.5, 2.5);
  AddHisto("sv_phi",             "Vtx #phi",                                            80, -4, 4    );
  AddHisto("sv_flightSig2D",     "Flight distance significance 2D",                     50,0.,80.    );
  AddHisto("sv_flight2D",        "Flight distance 2D",                                  50,0.,2.5    );
  AddHisto("sv_flight3D",        "Flight distance 3D",                                  50,0.,15.    );  
  AddHisto("sv_flight3DSig" ,    "flight distance significance 3D",	                50,0.,80.    );
  AddHisto("sv_multi_0"	  ,      "number of secondary vertex",                          6,-0.5,5.5   );
  AddHisto("sv_multi"	  ,      "number of secondary vertex",                          6,-0.5,5.5   );
  AddHisto("sv_mass"	  ,      "invariant mass of the secondary vertex",              50,0.,8.     );
  AddHisto("sv_chi2norm"  ,      "normalized chi2 of the secondary vertex",             50,0.,10.    );
  AddHisto("sv_tot_charge",      "Total charge",                                        21,-10.5,10.5);
  AddHisto("svnTrk",	         "Track multiplicity : SVnVertexTracks (centered)",     13,-0.5,12.5 );
  AddHisto("svnTrk_firstVxt",    "Track multiplicity : SVnFirstVertexTracks (centered)", 11,-0.5,10.5 ); 
  AddHisto("sv_flight3Derr",     "Flight distance error 3D",                            50,0.,0.2);
  AddHisto("sv_flight2Derr",     "Flight distance error 2D",                            50,0.,0.05);
  AddHisto("sv_mass_3trk"	,"invariant mass of the secondary vertex with at least 3 SV tracks",  50,0.,8.     );
  
  AddHisto("track_multi"  ,      "number of tracks in the jets",                40,-0.5,39.5  );
  AddHisto("trk_multi_sel"  ,    "number of selected tracks in the jets",       40,-0.5,39.5  );
  AddHisto("track_chi2"   ,      "normalized chi2 of the tracks",               100,0.,30.    );
  AddHisto("track_nHit" ,        "number of hits ",               35,-0.5, 34.5 );
  AddHisto("track_HPix"   ,      "number of hits in the Pixel",                 10,-0.5, 9.5  );
  
  AddHisto("track_IPs"    ,      "3D IP significance of all tracks",	        100,-35.,35.  );
  AddHisto("track_IPs1tr" ,      "3D IP significance of the first track",       100,-35.,35.  );
  AddHisto("track_IPs2tr" ,      "3D IP significance of the second track",      100,-35.,35.  );
  AddHisto("track_IP"     ,      "3D IP of all tracks",	                        100,-0.1,0.1);
  AddHisto("track_IP1tr"  ,      "3D IP of the first track",                    100,-0.1,0.1);
  AddHisto("track_IP2tr"  ,      "3D IP of the second track",                   100,-0.1,0.1); 
  AddHisto("track_IP2Ds"    ,    "2D IP significance of all tracks",	        100,-35.,35.  );
  AddHisto("track_IP2Ds1tr" ,    "2D IP significance of the first track",       100,-35.,35.  );
  AddHisto("track_IP2Ds2tr" ,    "2D IP significance of the second track",      100,-35.,35.  );
  AddHisto("track_IP2D"    ,     "2D IP of all tracks",	                        100,-0.1,0.1);
  AddHisto("track_IP2D1tr" ,     "2D IP of the first track",                    100,-0.1,0.1);
  AddHisto("track_IP2D2tr" ,     "2D IP of the second track",                   100,-0.1,0.1);
  AddHisto("track_IP2Derr1tr" ,  "2D IP error of the first track",              100,0,0.1);    
  AddHisto("track_IPerr1tr"   ,  "3D IP error of the first track",              100,0,0.1);  
  AddHisto("track_IP2Derr2tr" ,  "2D IP error of the second track",             100,0,0.1);    
  AddHisto("track_IPerr2tr"   ,  "3D IP error of the second track",             100,0,0.1); 
  AddHisto("track_IP2Derr" ,     "2D IP error",                                 100,0,0.1);    
  AddHisto("track_IPerr"   ,     "3D IP error",                                 100,0,0.1); 
  AddHisto("track_IPs3tr" ,      "3D IP significance of the third track",       100,-35.,35.  );
  AddHisto("track_IP3tr"  ,      "3D IP of the third track",                    100,-0.1,0.1);
  AddHisto("track_IPerr3tr"   ,  "3D IP error of the third track",              100,0,0.1); 
  AddHisto("track_IP2Ds3tr" ,    "2D IP significance of the second track",      100,-35.,35.  );
  AddHisto("track_IP2D3tr" ,     "2D IP of the third track",                    100,-0.1,0.1);
  AddHisto("track_IP2Derr3tr" ,  "2D IP error of the third track",              100,0,0.1);   
   
  AddHisto("track_len"     ,     "decay length",		                100,0,25.     );
  AddHisto("track_dist"    ,     "distance to the jet axis",                    100,0.,0.3    );
  AddHisto("track_dz"     ,     "transverse IP",                               100,-20,20  );  
  AddHisto("track_isfromSV",     "Track is from SV",                            2,-0.5, 1.5   );  
  AddHisto("track_pt"	  ,      "pT of all the tracks",	                80,0.,200.    );
  AddHisto("track_chi2_cut"     ,"normalized chi2 ",  	                        100,0.,30.    );
  AddHisto("track_nHit_cut"   ,"number of hits  ",               35,-0.5, 34.5 );
  AddHisto("track_HPix_cut"     ,"number of hits in the Pixel ",                 10,-0.5, 9.5  );
  AddHisto("track_len_cut"      ,"decay length ",		                100,0,25.     );
  AddHisto("track_dist_cut"     ,"distance to the jet axis ",                    100,0.,0.3    );
  AddHisto("track_dz_cut"      ,"transverse IP ",		                10,-0.5, 9.5  );  
  AddHisto("track_pt_cut"       ,"pT ",	                                        80,0.,200.);
  AddHisto("track_IP2D_cut"     ,"IP2D ",	                                100,-0.1,0.1);
   
  AddHisto("TCHE_extended1"	  ,"TCHE_extended1",				     70, -30.,30. );
  AddHisto("TCHP_extended1"	  ,"TCHP_extended1",				     70, -30.,30. );
  AddHisto("TCHE_extended2"	  ,"TCHE_extended2",				     100,-30.,30. );
  AddHisto("TCHP_extended2"	  ,"TCHP_extended2",				     100,-30.,30. );
  AddHisto("discri_ssche0",	 "SSVHE Discriminator",    80, -1., 7.   ); 
  AddHisto("discri_sschp0",	 "SSVHP Discriminator",    80, -1., 7.   ); 
   
  AddHisto("TCHE"	  ,"TCHE",				     50,0.,30. );
  AddHisto("TCHP"	  ,"TCHP",				     50,0.,30. );  
  AddHisto("JP" 	  ,"JP",				     50,0.,2.5 );
  AddHisto("JBP"	  ,"JBP",				     50,0.,8.  );
  AddHisto("SSV"	  ,"SSV",				     50,0.,7.  );
  AddHisto("SSVHP"	  ,"SSVHP",				     50,0.,7.  );
  AddHisto("CSV"	  ,"CSV",				     50,0.,1.  );
  
  AddHisto2D("seltrack_vs_jetpt", "sel track multiplicity vs jet pt",         30,60,1000, 100,-0.5,99.5);
  AddHisto2D("sv_mass_vs_flightDist3D", " SVMass vs SV 3D flight distance ",  100,0, 10,100,0,6);			
  AddHisto2D("avg_sv_mass_vs_jetpt","Avg SVMass vs jet pt",                   30,60,1000, 100,0,6);
  AddHisto2D("sv_deltar_jet_vs_jetpt","SVJetDeltaR vs jet pt",                25,60,300, 50,0.,0.5);  
  AddHisto2D("sv_deltar_sum_jet_vs_jetpt","SVvtxSumJetDeltaR vs jet pt",      25,60,300, 50,0.,0.5);
  AddHisto2D("sv_deltar_sum_dir_vs_jetpt","SVvtxSumVtxDirDeltaR vs jet pt",   25,60,300, 50,0.,0.5); 
  AddHisto2D("muon_ptrel_vs_jetpt","Muon_p{T}^{rel} vs jet pt",               30,60,1000,50,0,5);  
  AddHisto2D("muon_DeltaR_vs_jetpt","Muon1 DeltaR vs jet pt",                 30,60,1000,50,0,0.5);

  
  Nevent = 0;
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  //------------------------------------------------------------------------------------------------------------------//  
  //----------------------------------------EVENT LOOP ---------------------------------------------------------------// 
  //------------------------------------------------------------------------------------------------------------------//  
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    
    //-----------------------------------
    //is data or MC ?
    //-----------------------------------
    if ( pthat <1.)    isData=true;
    else              isData=false;
    
    if (isData) N_event_data_before_sel++;
    else        N_event_mc_before_sel++;
    //-----------------------------------
    //initialaze the weight at 1
    //-----------------------------------
    float ww=1; 
    
    //------------------------ cross section checks---------------------------------//
    if      ( pthat >=  15. && pthat <  30. ){
      n15_30++;
      use15_30=true;
    } 
    else if ( pthat >=  30. && pthat <  50. ){
      n30_50++;
      use30_50=true;
    } 
    else if ( pthat >=  50. && pthat <  80. ){
      n50_80++;
      use50_80=true;
    }
    else if ( pthat >=  80. && pthat < 120. ){
      n80_120++;
      use80_120=true;
    }
    else if ( pthat >= 120. && pthat < 170. ){
      n120_170++;
      use120_170=true;
    }
    else if ( pthat >= 170. && pthat < 300. ){
      n170_300++;
      use170_300=true;
    }
    else if ( pthat >= 300. && pthat < 470. ){
      n300_470++;
      use300_470=true;
    }
    else if ( pthat >= 470. && pthat < 600. ){
      n470_600++;
      use470_600=true;
    }
    else if ( pthat >= 600. && pthat <= 800. ){
      n600_800++;
      use600_800=true;
    }
    //-----------------------------------//
    //-------xs reweighting--------------//
    //-----------------------------------//

    ww=GetEvtWeight();
    //--------------------------------------------//  
    //-------------pile-up reweighting------------//  
    //--------------------------------------------//  
    
    if (!isData){
      if (nPUtrue >=59) nPUtrue=59;     
      ww*=WeightPU;
      
    }    
    
    if(!isData && WeightPU == 0  ) cout << "found null PU weights " <<  nPUtrue  << endl;
    if(!isData && WeightPU == 1  ) cout << "found null PU weights " <<  nPUtrue << endl;
    if(!isData && WeightPU > 10 ) cout << "diverging PU weights " << WeightPU << endl;
    
    
    //-----------------------------------
    //counter of events
    //-----------------------------------
    Nevent++;
    if(Nevent%50000 ==0 && Nevent!=0) cout << " number of processed events is " << Nevent << endl;
    
    //at least 1 jet in the event
    if (nJet<=0) continue;
    
    //-----------------------------------
    //Apply trigger selection
    //-----------------------------------
    // Triggers
    int trig100 = BitTrigger%100;
    int trig1000 = BitTrigger%1000;
    //int trig10000 = BitTrigger%10000;
    
    int trig1 = BitTrigger%10;
    int trig2 = (trig100 - trig1) / 10;
    int trig3 = (trig1000 - trig100) / 100;
    
    bool Jet30  = false, Jet40  = false, Jet60  = false, Jet80  = false;
    bool Jet110 = false, Jet140 = false, Jet150 = false, Jet190 = false;
    bool Jet200 = false, Jet240 = false, Jet260 = false;
    bool Jet300 = false, Jet320 = false;
    
    if ( trig2==1 || trig2==3 || trig2==5 || trig2==7 )  Jet40  = true;
    if ( trig2 >= 4 )				         Jet80  = true;
    if ( trig3==2 || trig3==3 || trig3==6 || trig3==7 )  Jet140 = true;
    if ( trig3 >= 4 )				         Jet200 = true;
    if ( trig1==1 || trig1==3 || trig1==5 || trig1==7 )  Jet260 = true;
    if ( trig1==2 || trig1==3 || trig1==6 || trig1==7 )  Jet320 = true;
    
    if ( IntCut ==  30 && !Jet30 )  continue;
    if ( IntCut ==  40 && !Jet40 )  continue;
    if ( IntCut ==  60 && !Jet60 )  continue;
    if ( IntCut ==  80 && !Jet80 )  continue;
    if ( IntCut == 110 && !Jet110 ) continue;
    if ( IntCut == 140 && !Jet140 ) continue;
    if ( IntCut == 150 && !Jet150 ) continue;
    if ( IntCut == 190 && !Jet190 ) continue;
    if ( IntCut == 200 && !Jet200 ) continue;
    if ( IntCut == 240 && !Jet240 ) continue;
    if ( IntCut == 260 && !Jet260 ) continue;
    if ( IntCut == 300 && !Jet300 ) continue;
    if ( IntCut == 320 && !Jet320 ) continue;
    
    if ( Year == 2012 && IntCut == 1 
	 && !Jet40 && !Jet80 && !Jet140
	 && !Jet200 && !Jet260 && !Jet320 ) continue;
    
    //-----------------------------------
    //Determine if there is at least 
    //on jet which pass the trigger
    //in the event => away from the TO
    //-----------------------------------
    
    bool JetPtCut = false;
    for (int ijet=0; ijet<nJet ; ijet++) {
      
      float ptjet = Jet_pt[ijet];
      //JetPtCut = true;
      //if (  ptjet > 100. && fabs(Jet_eta[ijet]) < EtaCut ) JetPtCut = true;
      //if(JetPtCut) cout << "ptjet 1 " <<   ptjet << "  fabs(Jet_eta[ijet]) " << fabs(Jet_eta[ijet]) << endl;
      
      if (      IntCut ==  40 && ptjet >  50. && fabs(Jet_eta[ijet])  < EtaCut ) JetPtCut = true;
      else if ( IntCut ==  80 && ptjet > 100. && fabs(Jet_eta[ijet])  < EtaCut ) JetPtCut = true;
      else if ( IntCut == 140 && ptjet > 160. && fabs(Jet_eta[ijet])  < EtaCut ) JetPtCut = true;
      else if ( IntCut == 200 && ptjet > 220. && fabs(Jet_eta[ijet])  < EtaCut ) JetPtCut = true;
      else if ( IntCut == 260 && ptjet > 300. && fabs(Jet_eta[ijet])  < EtaCut ) JetPtCut = true;
      else if ( IntCut == 320 && ptjet > 360. && fabs(Jet_eta[ijet])  < EtaCut ) JetPtCut = true;
      
    }
    if (!JetPtCut) continue;
    
    //-----------------------------------
    //Fill control plot
    //-----------------------------------
    
    if(isData){
      N_event_data_after_sel++;
      nPU_data       ->Fill(nPV);
    }
    else{
      //N_event_mc+=ww;      
      N_event_mc_after_sel++;
      nPU_mc         ->Fill(nPUtrue,ww);
      pt_hat         ->Fill(pthat,ww);
      nPV_mc         ->Fill(nPV,ww);
    }

      
    //-----------------------------------
    //Loop on jets 
    //-----------------------------------
    for (int ijet = 0; ijet < nJet; ijet++) {
      
      float ptjet    = Jet_pt[ijet];
      float etajet   = Jet_eta[ijet];
      float phijet   = Jet_phi[ijet];      
      float ntrkjet  = Jet_ntracks[ijet];  
      int   flav     = Jet_flavour[ijet];
      
          
      if (   ptjet  < PtMin  || ptjet  > PtMax   ) continue;
      if (   fabs(etajet) > EtaCut               ) continue;

      float mass_sv        =0.;
      int n_sv             =0.;
      float chi2norm_sv    =0.;
      float flightSig_sv   =0.;    
      float flight2DSig_sv =0.;    
      float sv_dR_jet      =0.;
      float sv_dR_dir_sum  =0.; 
      float sv_dR_jet_sum  =0.;
      float sv_en_rat      =0.; 
      float sv_abovC       =0.;  
      float sv_pt	   =0.;      
      float sveta         =0.; 
      float svphi         =0.; 
      float sv_flight3D    =0.;
      float sv_flight3Derr =0.;
      float sv_flight2D    =0.;
      float sv_flight2Derr =0.;
      int   sv_totchar     =0.;
      float sv_nTrk        =0.;
      float sv_1st_nTrk    =0.;
      
      int   idxFirstMuon = -1;
      
      float tche     = Jet_Ip2P[ijet];
      float tchp     = Jet_Ip3P[ijet];
      float jetproba = Jet_ProbaP[ijet];
      float jetbproba= Jet_BprobP[ijet];
      float ssvhe    = Jet_Svx[ijet] ;
      float ssvhp    = Jet_SvxHP[ijet];
      float csv      = Jet_CombSvx[ijet];
      
      //cout << "tche " << tche << endl;
      //cout << "tchp " << tchp << endl;
      
      bool isGluonSplit=false;
      
      
      //fill jet multiplicity
      if (!isData) {
        if (fabs(flav)==4)                  njet_c++;
	if (fabs(flav)<4 || fabs(flav)==21) njet_l++;
        njet_mc++;
        jet_pt_mc   ->Fill(ptjet,ww);
      }
      else njet_data++;
      
      //---------------------------------
      //test if b from gluon splitting
      //---------------------------------
      if (fabs(flav)==5){
	njet_b++;
	//----------------------//
	// is gluon splitting ?
	//----------------------//
	
	int nSplit=0;
	for (int k = 0; k < nBFromGSplit; k++) {
	  float dEta = Jet_eta[ijet] - bFromGSplit_eta[k];
	  float dPhi = Jet_phi[ijet] - bFromGSplit_phi[k];
	  if ( dPhi > 3.14159 ) dPhi = 2.*3.14159 - dPhi; 
	  float deltaR = TMath::Sqrt( dEta*dEta + dPhi*dPhi );
	  if ( deltaR < 0.5 ) nSplit++;
	}
	//----------------------//
	// is not gluon splitting  
	//----------------------//
	if  (nSplit >= 2 ) {
	  isGluonSplit=true;
	  njet_bfromg++;
	}
	else               isGluonSplit=false;
      }    
      
      FillHisto_floatFromMap("jet_multi",                  flav, isGluonSplit ,nJet     ,ww);
      FillHisto_floatFromMap("jet_pt_all",                 flav, isGluonSplit ,ptjet    ,ww);
      if (nSV > 0)FillHisto_floatFromMap("jet_pt_sv",      flav, isGluonSplit ,ptjet    ,ww);
      
      FillHisto_floatFromMap("jet_eta",     flav, isGluonSplit ,etajet   ,ww);
      FillHisto_floatFromMap("jet_phi",     flav, isGluonSplit ,phijet   ,ww);
      FillHisto_intFromMap(  "track_multi", flav, isGluonSplit ,ntrkjet  ,ww);
      
      //       cout << " ptjet:" << ptjet<<endl;
      //       cout << " etajet:" << etajet<<endl;
      //       cout << " phijet:" << phijet<<endl;
      //cout << " track_multi:" << ntrkjet<<"and last-first track: " << Jet_nLastTrack[ijet]-Jet_nFirstTrack[ijet] << endl;
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
         
      if ( produceJetProbaTree ) {
	
	for (int itrk=Jet_nFirstTrack[ijet]; itrk<Jet_nLastTrack[ijet] ; itrk++){

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
	    
	  if (Track_nHitAll[itrk]>=8)   passNhit=true;
	  if (Track_nHitPixel[itrk]>=2)   passPix= true;
	  if (fabs(Track_dz[itrk])<17)   passIPz=true;
	  if (Track_pt[itrk]>1)           passPt=true;
	  if (Track_chi2[itrk]<5)         passnormchi2=true;
	  if (fabs(Track_dist[itrk])<0.07)passtrkdist=true;
	  if (Track_length[itrk]<5)       passtrklen=true;
	  if (fabs(Track_IP2D[itrk])<0.2)   passTrackIP2D=true;
	  
	  if (!use_selected_tracks){
	    
	    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passTrackIP2D){
	      FillHisto_floatFromMap("track_len_cut",          flav, isGluonSplit ,Track_length[itrk] , ww);
	    }
	    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrklen && passTrackIP2D){
	      FillHisto_floatFromMap("track_dist_cut",         flav, isGluonSplit ,fabs(Track_dist[itrk])   , ww);
	    }	    
	    if (passNhit && passPix && passIPz && passPt && passtrkdist && passtrklen && passTrackIP2D){
	      FillHisto_floatFromMap("track_chi2_cut",         flav, isGluonSplit ,Track_chi2[itrk]	   ,ww);
	    }	    
	    if (passNhit && passPix && passIPz && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
	      FillHisto_floatFromMap("track_pt_cut",           flav, isGluonSplit ,Track_pt[itrk]     , ww);
	    }	    
	    if (passNhit && passPix && passPt && passnormchi2 && passtrkdist && passtrklen){
	      FillHisto_floatFromMap("track_dz_cut",          flav, isGluonSplit ,Track_dz[itrk]      ,ww);
	    }
	    if (passNhit && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
	      FillHisto_intFromMap(  "track_HPix_cut",         flav, isGluonSplit ,Track_nHitPixel[itrk],ww);
	    }	
	    if (passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
	      FillHisto_intFromMap(  "track_nHit_cut",       flav, isGluonSplit ,Track_nHitAll[itrk],ww);	  
	    }
	    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen ){
	      FillHisto_intFromMap(  "track_IP2D_cut",         flav, isGluonSplit ,Track_IP2D[itrk],ww);	  
	    }	    
	  }
	  if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
	    ntracksel++;
	    
	    FillHisto_floatFromMap("track_chi2",    flav, isGluonSplit ,Track_chi2[itrk]	   ,ww);
	    FillHisto_intFromMap(  "track_nHit",  flav, isGluonSplit ,Track_nHitAll[itrk],ww);
	    FillHisto_intFromMap(  "track_HPix",    flav, isGluonSplit ,Track_nHitPixel[itrk],ww);
	    FillHisto_floatFromMap("track_IPs",     flav, isGluonSplit ,Track_IPsig[itrk]    ,ww);
	    FillHisto_floatFromMap("track_IP",      flav, isGluonSplit ,Track_IP[itrk]       ,ww);
	    FillHisto_floatFromMap("track_IP2Ds",   flav, isGluonSplit ,Track_IP2Dsig[itrk]  ,ww);
	    FillHisto_floatFromMap("track_IP2D",    flav, isGluonSplit, Track_IP2D[itrk]     ,ww);
	    FillHisto_floatFromMap("track_IP2Derr", flav, isGluonSplit, Track_IP2Derr[itrk]  ,ww);	  
	    FillHisto_floatFromMap("track_IPerr",   flav, isGluonSplit, Track_IPerr[itrk]    ,ww);	  	  
	    FillHisto_floatFromMap("track_dz",      flav, isGluonSplit ,Track_dz[itrk]      ,ww);	  
	    FillHisto_intFromMap(  "track_isfromSV",flav, isGluonSplit ,Track_isfromSV[itrk] ,ww);	  
	    FillHisto_floatFromMap("track_len",     flav, isGluonSplit ,Track_length[itrk]  , ww);
	    FillHisto_floatFromMap("track_dist",    flav, isGluonSplit ,fabs(Track_dist[itrk])    , ww);
	    FillHisto_floatFromMap("track_pt",      flav, isGluonSplit ,Track_pt[itrk]      , ww);	  
	    //           cout << "track_chi2 :" << Track_chi2[itrk]<<endl;
	    //           cout << "track_HStrip :" << Track_nHitStrip[itrk]<<endl;
	    //           cout << "track_HPix :" << Track_nHitPixel[itrk]<<endl;
	    //           cout << "track_IP :" << Track_IP[itrk]<<endl;	
	    
	  //Tracks sorted by IP
	
	    Float_t sig  =Track_IP[itrk]/Track_IPerr[itrk];
	    Float_t sig2D=Track_IP2D[itrk]/Track_IP2Derr[itrk];
            if (sig>sig1_ip) {
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
	    else if (sig>sig2_ip) {
	      sig3_ip=sig2_ip;
	      sig2_ip=sig;	      
	      sig32D_ip=sig22D_ip;
	      sig22D_ip=sig2D;
	      n3_ip=n2_ip;
	      n2_ip=itrk;
	    }
	    else if (sig>sig3_ip) {
	      sig3_ip=sig;
	      sig32D_ip=sig2D;
	      n3_ip=itrk;
            }	      

	  }//end selected tracks

	}//end tracks loop
	
	if (n1_ip>-1)  {
	  FillHisto_floatFromMap("track_IPs1tr",    flav, isGluonSplit ,sig1_ip               , ww);
	  FillHisto_floatFromMap("track_IP1tr",     flav, isGluonSplit ,Track_IP[n1_ip]       , ww);
	  FillHisto_floatFromMap("track_IPerr1tr",  flav, isGluonSplit ,Track_IPerr[n1_ip]    , ww);
	  FillHisto_floatFromMap("track_IP2Ds1tr",  flav, isGluonSplit ,sig12D_ip             , ww);
	  FillHisto_floatFromMap("track_IP2D1tr",   flav, isGluonSplit ,Track_IP2D[n1_ip]     , ww);
	  FillHisto_floatFromMap("track_IP2Derr1tr",flav, isGluonSplit ,Track_IP2Derr[n1_ip]  , ww);		  
	}

	if (n2_ip>-1) {
	  FillHisto_floatFromMap("track_IPs2tr",    flav, isGluonSplit ,sig2_ip               , ww);
	  FillHisto_floatFromMap("track_IP2tr",     flav, isGluonSplit ,Track_IP[n2_ip]       , ww);
	  FillHisto_floatFromMap("track_IPerr2tr",  flav, isGluonSplit ,Track_IPerr[n2_ip]    , ww);	
	  FillHisto_floatFromMap("track_IP2Ds2tr",  flav, isGluonSplit ,sig22D_ip             , ww);
	  FillHisto_floatFromMap("track_IP2D2tr",   flav, isGluonSplit ,Track_IP2D[n2_ip]     , ww);
	  FillHisto_floatFromMap("track_IP2Derr2tr",flav, isGluonSplit ,Track_IP2Derr[n2_ip]  , ww);
	}
                
	if (n3_ip>-1) {
	  FillHisto_floatFromMap("track_IPs3tr",    flav, isGluonSplit ,sig3_ip               , ww);
	  FillHisto_floatFromMap("track_IP3tr",     flav, isGluonSplit ,Track_IP[n3_ip]       , ww);
	  FillHisto_floatFromMap("track_IPerr3tr",  flav, isGluonSplit ,Track_IPerr[n3_ip]    , ww);	
	  FillHisto_floatFromMap("track_IP2Ds3tr",  flav, isGluonSplit ,sig32D_ip             , ww);
	  FillHisto_floatFromMap("track_IP2D3tr",   flav, isGluonSplit ,Track_IP2D[n3_ip]     , ww);
	  FillHisto_floatFromMap("track_IP2Derr3tr",flav, isGluonSplit ,Track_IP2Derr[n3_ip]  , ww);
	}

	
	FillHisto_intFromMap(        "trk_multi_sel",     flav, isGluonSplit ,ntracksel	         , ww);  
	FillHisto2D_int_floatFromMap("seltrack_vs_jetpt", flav, isGluonSplit ,ptjet ,  ntracksel , ww);
	
	//---------------------------------
	//fill information related to SV
	//---------------------------------
	
	  n_sv           = Jet_SV_multi[ijet];	  
	  FillHisto_intFromMap(  "sv_multi_0",      flav, isGluonSplit ,n_sv 	 ,         ww);

	if (n_sv>0){  
	  chi2norm_sv    = SV_chi2[Jet_nFirstSV[ijet]]/SV_ndf[Jet_nFirstSV[ijet]];
	  flightSig_sv   = SV_flight[Jet_nFirstSV[ijet]]/SV_flightErr[Jet_nFirstSV[ijet]];
	  flight2DSig_sv = SV_flight2D[Jet_nFirstSV[ijet]]/SV_flight2DErr[Jet_nFirstSV[ijet]];
	  mass_sv        =Jet_SvxMass[Jet_nFirstSV[ijet]];
	  sv_dR_jet      =SV_deltaR_jet[Jet_nFirstSV[ijet]];
	  sv_dR_dir_sum  =SV_deltaR_sum_dir[Jet_nFirstSV[ijet]];
	  sv_dR_jet_sum  =SV_deltaR_sum_jet[Jet_nFirstSV[ijet]];
	  sv_en_rat      =SV_energy_ratio[Jet_nFirstSV[ijet]];
	  sv_abovC       =SV_aboveC[Jet_nFirstSV[ijet]];	    
	  sv_pt          =SV_vtx_pt[Jet_nFirstSV[ijet]];
	  sveta          =SV_vtx_eta[Jet_nFirstSV[ijet]];
	  svphi          =SV_vtx_phi[Jet_nFirstSV[ijet]];
	  
	  
	  sv_flight3D    =SV_flight[Jet_nFirstSV[ijet]] ;  
	  sv_flight3Derr =SV_flight[Jet_nFirstSV[ijet]]; 
	  sv_flight2D    =SV_flight2D[Jet_nFirstSV[ijet]];    
	  sv_flight2Derr =SV_flight2DErr[Jet_nFirstSV[ijet]];    
	  sv_totchar     =SV_totCharge[Jet_nFirstSV[ijet]] ;
	     
	  sv_nTrk        =SV_nTrk[Jet_nFirstSV[ijet]] ;  
	  sv_1st_nTrk    =SV_nTrk_firstVxt[Jet_nFirstSV[ijet]];
	 //-------------------------//
	//-----SV histograms-----//
        //-------------------------//   
	  FillHisto_intFromMap(  "sv_multi",     flav, isGluonSplit ,n_sv ,  ww);
	  FillHisto_floatFromMap("sv_chi2norm",     flav, isGluonSplit ,chi2norm_sv        , ww);
	  FillHisto_floatFromMap("sv_mass",         flav, isGluonSplit ,mass_sv,             ww);
	  FillHisto_floatFromMap("sv_deltaR_jet",   flav, isGluonSplit ,sv_dR_jet,           ww);
	  FillHisto_floatFromMap("sv_deltaR_sumJet",flav, isGluonSplit ,sv_dR_dir_sum,       ww);
	  FillHisto_floatFromMap("sv_deltaR_sumDir",flav, isGluonSplit ,sv_dR_jet_sum,       ww);
	  FillHisto_floatFromMap("sv_en_ratio",     flav, isGluonSplit ,sv_en_rat,           ww);
	  FillHisto_floatFromMap("sv_aboveC",       flav, isGluonSplit ,sv_abovC,            ww);
	  FillHisto_floatFromMap("sv_pt",           flav, isGluonSplit ,sv_pt,               ww);
	  FillHisto_floatFromMap("sv_flight2D",     flav, isGluonSplit ,sv_flight2D,         ww);
	  FillHisto_floatFromMap("sv_flight2Derr",  flav, isGluonSplit ,sv_flight2Derr,      ww);
	  FillHisto_floatFromMap("sv_flightSig2D",  flav, isGluonSplit ,flight2DSig_sv,      ww);
	  FillHisto_intFromMap("sv_tot_charge",     flav, isGluonSplit ,sv_totchar,          ww);
	  FillHisto_intFromMap(  "svnTrk",          flav, isGluonSplit ,sv_nTrk,             ww);
	  FillHisto_intFromMap(  "svnTrk_firstVxt", flav, isGluonSplit ,sv_1st_nTrk,         ww);
	  FillHisto_floatFromMap("sv_eta",          flav, isGluonSplit ,sveta,               ww);
	  FillHisto_floatFromMap("sv_phi",          flav, isGluonSplit ,svphi,               ww);	
	  FillHisto_floatFromMap("sv_flight3D",     flav, isGluonSplit ,sv_flight3D,         ww);
	  FillHisto_floatFromMap("sv_flight3Derr",  flav, isGluonSplit ,sv_flight3Derr,      ww);
	  FillHisto_floatFromMap("sv_flight3DSig",  flav, isGluonSplit ,flightSig_sv,        ww);
	
	  if (sv_nTrk >2)FillHisto_floatFromMap("sv_mass_3trk", flav, isGluonSplit ,mass_sv,ww);
	
	  FillHisto2D_float_floatFromMap("sv_mass_vs_flightDist3D"     ,flav,isGluonSplit ,sv_flight3D,mass_sv,ww);	
	  FillHisto2D_float_floatFromMap("avg_sv_mass_vs_jetpt"        ,flav,isGluonSplit ,ptjet,mass_sv,ww);
	  FillHisto2D_float_floatFromMap("sv_deltar_jet_vs_jetpt"      ,flav,isGluonSplit ,ptjet,sv_dR_jet,ww);
	  FillHisto2D_float_floatFromMap("sv_deltar_sum_jet_vs_jetpt"  ,flav,isGluonSplit ,ptjet,sv_dR_dir_sum,ww);
	  FillHisto2D_float_floatFromMap("sv_deltar_sum_dir_vs_jetpt"  ,flav,isGluonSplit ,ptjet,sv_dR_jet_sum,ww);	    
	    
	    
	}

      }//end produce jetProbaTree
	    

	
	//FillHisto_floatFromMap("sv_dist_jet_axis",flav, isGluonSplit ,SV_vtxDistJetAxis[0],ww);
	
	//         cout << " sv_multi:" << n_sv<<endl;
	//         cout << " sv_chi2norm:" << chi2norm_sv<<endl;
	//         cout << " sv_mass:" << Jet_SvxMass[0]<<endl;
	//         cout << " sv_flightSig:" << flightSig_sv<<endl;		      
	
            
      //Taggers
      FillHisto_floatFromMap("TCHE",  flav, isGluonSplit, tche	  ,   ww);
      FillHisto_floatFromMap("TCHP",  flav, isGluonSplit, tchp	  ,   ww);
      FillHisto_floatFromMap("JP",    flav, isGluonSplit, jetproba  , ww);
      FillHisto_floatFromMap("JBP",   flav, isGluonSplit, jetbproba , ww);
      FillHisto_floatFromMap("SSV",   flav, isGluonSplit, ssvhe	  ,   ww);
      FillHisto_floatFromMap("SSVHP", flav, isGluonSplit, ssvhp	  ,   ww);
      FillHisto_floatFromMap("CSV",   flav, isGluonSplit, csv	  ,   ww);
      //       cout << " TCHE:" << tche<<endl;
      //       cout << " TCHP:" << tchp<<endl;
      //       cout << " JP:" << jetproba<<endl;
      //       cout << " JBP:" << jetbproba<<endl;      
      //       cout << " SSV:" << ssvhe<<endl;
      //       cout << " SSVHP:" << ssvhp<<endl;
      //       cout << " CSV:" << csv<<endl;
      
      
      FillHisto_floatFromMap("TCHE_extended1",  flav, isGluonSplit, tche  , ww);
      FillHisto_floatFromMap("TCHP_extended1",  flav, isGluonSplit, tchp  , ww);
      
      FillHisto_floatFromMap("TCHE_extended2",  flav, isGluonSplit, tche  , ww);
      FillHisto_floatFromMap("TCHP_extended2",  flav, isGluonSplit, tchp  , ww);
      
      FillHisto_floatFromMap("discri_ssche0",   flav, isGluonSplit, ssvhe , ww);
      FillHisto_floatFromMap("discri_sschp0",   flav, isGluonSplit, ssvhp , ww);      
      
      //---------------------------------
      //fill information related to muons
      //---------------------------------
      int nselmuon = 0;
      int nmu = 0;
      if (nMuon>0){
        for (int imu=0; imu< nMuon; imu++){
          if (Muon_IdxJet[imu]==ijet ){
	    nmu++;
	    if (passMuonSelection(imu, ijet)){
	      if(nselmuon == 0) {
	        idxFirstMuon = imu;
	        nselmuon++;
	      }
	    }
          }
	}
      }
         
      
      FillHisto_intFromMap(  "muon_multi_sel",  flav, isGluonSplit , nselmuon   ,ww);
      FillHisto_intFromMap(  "muon_multi",      flav, isGluonSplit , nmu        ,ww);
      
      
      if(idxFirstMuon > -1){
        FillHisto_floatFromMap("mu_ptrel",    flav, isGluonSplit ,Muon_ptrel[idxFirstMuon] ,ww);
        FillHisto_floatFromMap("mu_chi2",     flav, isGluonSplit ,Muon_chi2[idxFirstMuon]  ,ww);
        FillHisto_floatFromMap("muon_Pt",     flav, isGluonSplit, Muon_pt[idxFirstMuon] ,     ww);
	FillHisto_floatFromMap("muon_eta",    flav, isGluonSplit, Muon_eta[idxFirstMuon] ,    ww);
	FillHisto_floatFromMap("muon_phi",    flav, isGluonSplit, Muon_phi[idxFirstMuon] ,    ww);
        FillHisto_floatFromMap("muon_Ip3d",   flav, isGluonSplit, Muon_IP[idxFirstMuon] ,     ww);     
        FillHisto_floatFromMap("muon_Ip2d",   flav, isGluonSplit, Muon_IP2D[idxFirstMuon] ,   ww);     
        FillHisto_floatFromMap("muon_Sip3d",  flav, isGluonSplit, Muon_IPsig[idxFirstMuon] ,  ww);     
        FillHisto_floatFromMap("muon_Sip2d",  flav, isGluonSplit, Muon_IP2Dsig[idxFirstMuon] ,ww); 
	//         cout << "muon_Pt :" <<Muon_pt[idxFirstMuon] <<endl;
	// 	cout << "muon_eta :" <<Muon_eta[idxFirstMuon] <<endl;
	// 	cout << "muon_phi :" <<Muon_phi[idxFirstMuon] <<endl;
// 	
	
	
        TLorentzVector themuon, thejet;
	
        thejet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], 0);
        themuon.SetPtEtaPhiM(Muon_pt[idxFirstMuon], Muon_eta[idxFirstMuon], Muon_phi[idxFirstMuon], 0);
	
        FillHisto_floatFromMap("muon_DeltaR",         flav, isGluonSplit, themuon.DeltaR(thejet) , ww);
	//cout << "muon_DeltaR :" << themuon.DeltaR(thejet) <<endl;
	
	
	FillHisto2D_float_floatFromMap("muon_ptrel_vs_jetpt", flav, isGluonSplit,ptjet,Muon_ptrel[idxFirstMuon],ww);
        FillHisto2D_float_floatFromMap("muon_DeltaR_vs_jetpt",flav, isGluonSplit,ptjet,themuon.DeltaR(thejet),ww);
      }

      
    }
    //----------------------------------
    //End Loop on jets 
    //-----------------------------------
    
    //---------------------------------
    //fill jet multiplicity
    //---------------------------------
    //if(njet_mc > 0 || njet_data > 0){ 
    if(isData) HistoBtag[4]->Fill(njet_data);
    else{
      HistoBtag[0]->Fill(float(njet_b)     /float(njet_mc) ,ww);
      HistoBtag[1]->Fill(float(njet_bfromg)/float(njet_mc) ,ww);
      HistoBtag[2]->Fill(float(njet_c)     /float(njet_mc) ,ww);
      HistoBtag[3]->Fill(float(njet_l)     /float(njet_mc) ,ww);
    }
    
    
    
    //cout << "---------------------End event---------------------------" << endl;
    
  }//----------------------------------------End loop events ----------------------------------------------------------------------------------------------//  
  
  
  cout << "----------------Configuration-------------------" <<endl;
  cout <<" Run over "<< N_event_data_before_sel<<" data events and " << N_event_mc_before_sel<< " mc events."<<endl;
  cout <<" After selection -> "<< N_event_data_after_sel<<" data events and "<< N_event_mc_after_sel<< "mc events left" <<endl;   
  // //   
  if  (use15_30)  cout << "Run with QCD15_30   sample with " << n15_30  << " events" <<endl; 
  if  (use30_50)  cout << "Run with QCD30_50   sample with " << n30_50  << " events" <<endl;
  if  (use50_80)  cout << "Run with QCD50_80   sample with " << n50_80  << " events" <<endl;
  if  (use80_120) cout << "Run with QCD80_120  sample with " << n80_120 << " events" <<endl;
  if  (use120_170)cout << "Run with QCD120_170 sample with " << n120_170<< " events" <<endl;
  if  (use170_300)cout << "Run with QCD170_300 sample with " << n170_300<< " events" <<endl;
  if  (use300_470)cout << "Run with QCD300_470 sample with " << n300_470<< " events" <<endl;
  if  (use470_600)cout << "Run with QCD470_600 sample with " << n470_600<< " events" <<endl;
  if  (use600_800)cout << "Run with QCD600_800 sample with " << n600_800<< " events" <<endl;

  
  myfile->cd();
  for (unsigned int i=0; i<HistoBtag.size(); i++) {
    HistoBtag[i]->Write();
  }
  for (unsigned int i=0; i<HistoBtag2D.size(); i++) {
    HistoBtag2D[i]->Write();
  }  

}

void CommPlotProducer::AddHisto(TString name, TString title, int nbins, float min, float max)  {
  
  TH1F* h_b      = new TH1F(name+"_b",title+"_b",nbins,min,max);
  TH1F* h_bfromg = new TH1F(name+"_bfromg",title+"_bfromg",nbins,min,max);  
  TH1F* h_c      = new TH1F(name+"_c",title+"_c",nbins,min,max);  
  TH1F* h_l      = new TH1F(name+"_l",title+"_l",nbins,min,max);
  TH1F* h_data   = new TH1F(name+"_data",title+"_data",nbins,min,max);
  
  
  h_b        ->Sumw2();
  h_bfromg   ->Sumw2();  
  h_c        ->Sumw2();  
  h_l        ->Sumw2(); 
  h_data     ->Sumw2();
  
  HistoBtag.push_back(h_b);
  HistoBtag.push_back(h_bfromg);  
  HistoBtag.push_back(h_c);  
  HistoBtag.push_back(h_l);  
  HistoBtag.push_back(h_data);  
  HistoBtag_map[name] = numb_histo;
  
  numb_histo++;
  
}


void CommPlotProducer::FillHisto_float(int flavour, bool isGS, int number, float value, float weight)  {
  
  if (!isData){
    if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
    
  }  
  else                                              HistoBtag[number*5 +4]->Fill(value);
  
  
}
void CommPlotProducer::FillHisto_int(int flavour, bool isGS, int number, int value, float weight)  {
  
  if (!isData){
    if (fabs(flavour)==5 && !isGS)             HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
  }  
  else                                         HistoBtag[number*5 +4]->Fill(value);
  
}


void CommPlotProducer::FillHisto_floatFromMap(TString name, int flavour, bool isGS, float value, float weight)  {
  
  
  int number = HistoBtag_map[name] ;
  if (!isData){
    if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
    
  }  
  else                                              HistoBtag[number*5 +4]->Fill(value);
  
   
}


void CommPlotProducer::FillHisto_intFromMap(TString name, int flavour, bool isGS, int value, float weight)  {
  
  int number = HistoBtag_map[name] ;
  if (!isData){
    if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*5 +0]->Fill(value,weight);
    else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
    else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
  }  
  else                                         HistoBtag[number*5 +4]->Fill(value);
  
}

//-----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------2D PLOTS-----------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//


void CommPlotProducer::AddHisto2D(TString name, TString title, int nbins, float min, float max, int nbins2, float
min2, float max2)  {
  
  TH2F* h_b      = new TH2F(name+"_b",title+"_b",nbins,min,max,nbins2,min2,max2);
  TH2F* h_bfromg = new TH2F(name+"_bfromg",title+"_bfromg",nbins,min,max,nbins2,min2,max2);  
  TH2F* h_c      = new TH2F(name+"_c",title+"_c",nbins,min,max,nbins2,min2,max2);  
  TH2F* h_l      = new TH2F(name+"_l",title+"_l",nbins,min,max,nbins2,min2,max2);
  TH2F* h_data   = new TH2F(name+"_data",title+"_data",nbins,min,max,nbins2,min2,max2);
  
  
  h_b        ->Sumw2();
  h_bfromg   ->Sumw2();  
  h_c        ->Sumw2();  
  h_l        ->Sumw2(); 
  h_data     ->Sumw2();
  
  HistoBtag2D.push_back(h_b);
  HistoBtag2D.push_back(h_bfromg);  
  HistoBtag2D.push_back(h_c);  
  HistoBtag2D.push_back(h_l);  
  HistoBtag2D.push_back(h_data);  
  HistoBtag2D_map[name] = numb_histo2D;
  numb_histo2D++;
  
}



void CommPlotProducer::FillHisto2D_int_floatFromMap(TString name, int flavour, bool isGS, int value, float value2, float weight)  {
  
  
  int number = HistoBtag2D_map[name] ;
  if (!isData){
    if (fabs(flavour)==5 && !isGS)                  HistoBtag2D[number*5 +0]->Fill(value,value2,weight);
    else if (fabs(flavour)==5 && isGS)              HistoBtag2D[number*5 +1]->Fill(value,value2,weight);
    else if (fabs(flavour)==4)                      HistoBtag2D[number*5 +2]->Fill(value,value2,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*5 +3]->Fill(value,value2,weight);
    
  }  
  else                                              HistoBtag2D[number*5 +4]->Fill(value,value2);
  
   
}

void CommPlotProducer::FillHisto2D_float_floatFromMap(TString name, int flavour, bool isGS, float value, float value2, float weight)  {
  
  
  int number = HistoBtag2D_map[name] ;
  if (!isData){
    if (fabs(flavour)==5 && !isGS)                  HistoBtag2D[number*5 +0]->Fill(value,value2,weight);
    else if (fabs(flavour)==5 && isGS)              HistoBtag2D[number*5 +1]->Fill(value,value2,weight);
    else if (fabs(flavour)==4)                      HistoBtag2D[number*5 +2]->Fill(value,value2,weight); 
    else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*5 +3]->Fill(value,value2,weight);
    
  }  
  else                                              HistoBtag2D[number*5 +4]->Fill(value,value2);
  
   
}


//-----------------------------------------------------------------------------------------------------------------//

bool CommPlotProducer::passMuonSelection(int muidx, int ijet){
  
  
  TLorentzVector muon, jet;
  
  
  jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], 0);
  muon.SetPtEtaPhiM(Muon_pt[muidx], Muon_eta[muidx], Muon_phi[muidx], 0);
  
  //double deltaRMuJet = 0;
  
  
  bool cut_mu_pass=false;
  if (Muon_pt[muidx]> 5	&& TMath::Abs(Muon_eta[muidx]) < 2.4 && Muon_isGlobal[muidx] == 1     &&
      Muon_nMuHit[muidx]> 0 && Muon_nMatched[muidx]>1	 && Muon_nTkHit[muidx]>10 &&
      Muon_nPixHit[muidx]> 1	&& Muon_nOutHit[muidx]<3  && Muon_chi2Tk[muidx]< 10    &&
      Muon_chi2[muidx]< 10	//&& Muon_vz[muidx]< 2  
      && 	jet.DeltaR(muon) < 0.4
      && TMath::Abs(Muon_vz[muidx]-PV_z[0]) < 2.) 
    cut_mu_pass=true;
  
  //cout << "jet.DeltaR(muon)  " << jet.DeltaR(muon)  << endl;
  
  //if(!cut_mu_pass) cout << " not found a good muon " << endl;
  
  
  return cut_mu_pass;
  
}


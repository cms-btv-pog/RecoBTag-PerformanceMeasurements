#define CommPlotProducer_cxx
#include "CommPlotProducer.h"

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
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>

using namespace std;


//---------------------------------------------------------------------------------------//
//-------------Fill number of events for each sample--------------------------------------//
//----------------------------------------------------------------------------------------//

void CommPlotProducer::Fill_nevent(double n15,double n20,double n30,double n50,double n80,double n120,double n170,double n300,double n470,double n600, double n800, double n1000){

 
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
    nmc_evt_vect[10]=n800;
    nmc_evt_vect[11]=n1000;
  
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
        if      ( pthat >=  15. && pthat <  20. ){
            use15_20=true;
            n15_20++;
        } 
        else if ( pthat >=  20. && pthat <  30. ){
            use20_30=true;    
            n20_30++;
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
        else if ( pthat >= 800. && pthat <= 1000. ){
            use800_1000=true;    
            n800_1000++;
        }
        else if ( pthat >= 1000. ){
            use1000_inf=true;    
            n1000_inf++;
        }

        if      ( pthat >=  15. && pthat <  30. ){
            use15_30=true;
            n15_30++;
        } 
        else if ( pthat >= 120. && pthat < 150. ){
            use120_150=true;    
            n120_150++;
        }
        else if ( pthat >= 150. ){
            use150_inf=true;    
            n150_inf++;
        }

    }  
  
    if  (use15_20)  cout << "Run with QCD15_20   sample with " << n15_20  << " events" <<endl; 
    if  (use20_30)  cout << "Run with QCD20_30   sample with " << n20_30  << " events" <<endl; 
    if  (use30_50)  cout << "Run with QCD30_50   sample with " << n30_50  << " events" <<endl;
    if  (use50_80)  cout << "Run with QCD50_80   sample with " << n50_80  << " events" <<endl;
    if  (use80_120) cout << "Run with QCD80_120  sample with " << n80_120 << " events" <<endl;
    if  (use120_170)cout << "Run with QCD120_170 sample with " << n120_170<< " events" <<endl;
    if  (use170_300)cout << "Run with QCD170_300 sample with " << n170_300<< " events" <<endl;
    if  (use300_470)cout << "Run with QCD300_470 sample with " << n300_470<< " events" <<endl;
    if  (use470_600)cout << "Run with QCD470_600 sample with " << n470_600<< " events" <<endl;
    if  (use600_800)cout << "Run with QCD600_800 sample with " << n600_800<< " events" <<endl;
    if  (use800_1000)cout << "Run with QCD800_1000 sample with " << n800_1000<< " events" <<endl;
    if  (use1000_inf)cout << "Run with QCD1000 sample with " << n1000_inf<< " events" <<endl;
    cout << endl;
    if (use15_30) cout << "Run with QCD15_30   sample with " << n15_30  << " events" <<endl;
    if (use150_inf) cout << "Run with QCD150   sample with " << n150_inf  << " events" <<endl;

    cout << endl;
    cout << endl;

    /*
      pythia Inclusive  8TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000        --> 10 fichiers
      pythia MuEnriched 8TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 12 fichies
      herwig Inclusive  8TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 11 fichiers
      herwig MuEnriched 8TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 12 fichies

      pythia Inclusive  7TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000;       --> 10 fichiers
      pythia MuEnriched 7TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-150, n150-plus                                                  --> 7 fichiers
    */

    cout << " To write in runCode.C " << endl;
    if (qcdtype==0) { // inclusive qcd 
        cout << " double   n15    = 0. ; " <<endl;
        cout << " double   n20    = "<< n15_30 << "; " << endl;
    }
    else { // MuEnriched qcd
        cout << " double   n15    = "<< n15_20 << "; " << endl;
        cout << " double   n20    = "<< n20_30 << "; " << endl;
    }
    cout << " double   n30    = "<< n30_50 << "; " << endl;
    cout << " double   n50    = "<< n50_80 << "; " << endl;
    cout << " double   n80    = "<< n80_120 << "; " << endl;
    if (sqrtstev!=7) { // 8 TeV
        cout << " double   n120  = "<< n120_170 << "; " << endl;
        cout << " double   n170  = "<< n170_300 << "; " << endl;
        cout << " double   n300  = "<< n300_470 << "; " << endl;
        cout << " double   n470  = "<< n470_600 << "; " << endl;
        cout << " double   n600  = "<< n600_800 << "; " << endl;
        cout << " double   n800  = "<< n800_1000 << "; " << endl;
        cout << " double   n1000 = "<< n1000_inf << "; " << endl;
    }
    else if (qcdtype==1) { // 7 TeV MuEnriched qcd
        cout << " double   n120  = "<< n120_150 << "; " << endl;
        cout << " double   n170  = "<< n150_inf << "; " << endl;
        cout << " double   n300  = 0.; " << endl;
        cout << " double   n470  = 0.; " << endl; 
        cout << " double   n600  = 0.; " << endl;
        cout << " double   n800  = 0.; " << endl;
        cout << " double   n1000 = 0.; " << endl;
    }
    cout << " m.Fill_nevent(n15,n20,n30,n50,n80,n120,n170,n300,n470,n600,n800,n1000);" << endl;


}

//------------------------------------------------------------------------------//
//-------------Set PU--------------------------------------------------------//
//-----------------------------------------------------------------------------//
void CommPlotProducer::SetPU(vector<float> PUvector, TString PUdataFile){


    vector<float> mc_vect;
    vector<float> data_vect;
  
  
    TFile *filepuest = new TFile(PUdataFile,"READ");
    TH1F* npu_dat= (TH1F*) filepuest->Get("pileup");
    
    for( int i=0; i<60; ++i) {
        mc_vect.push_back(PUvector[i]);
    }
    for (int i=0; i<(int)PUvector.size(); i++){
        data_vect.push_back(npu_dat->GetBinContent(i+1));
    }
  
    LumiWeights = reweight::LumiReWeighting(mc_vect,data_vect);
  
}

void CommPlotProducer::SetPU2012_S7(TString PUdataFile){
 
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
  
}

void CommPlotProducer::SetPU2012_S10(TString PUdataFile){
 
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
  

}
void CommPlotProducer::SetPV(){


    vector<float> mc_vect;
    vector<float> data_vect;
    puweight=false;
  
    //TFile *filepuest = new TFile(filePUstring);
    // Note, in this version SetPV is dummy ;-)
    TFile *filepuest = new TFile("/home/fynu/bfrancois/bTag/CMSSW_7_4_5/src/RecoBTag/PerformanceMeasurements/test/BTagAnalyzerMacros/inclusiveQCD_ptBin15to300_star1_TrigPFjet40_Jet60_250_notPUReweighted_goldenJSON/inclusiveQCD_ptBin15to300_star1_TrigPFjet40_Jet60_250_notPUReweighted_goldenJSON.root");
    //TFile *filepuest = new TFile("/home/fynu/bfrancois/bTag/CMSSW_7_4_5/src/RecoBTag/PerformanceMeasurements/test/BTagAnalyzerMacros/muEnQCD_ptBin20to300_star1_TrigDijet20_Jet60_250_notPUReweighted_goldenJSON/muEnQCD_ptBin20to300_star1_TrigDijet20_Jet60_250_notPUReweighted_goldenJSON.root");
    TH1D* npu_mc= (TH1D*) filepuest->Get("nPV_mc_unw");
    TH1D* npu_dat= (TH1D*) filepuest->Get("nPV_data");
    
    float n_integral_mc = npu_mc->Integral();
    float n_integral_da = npu_dat->Integral();
    npu_mc->Scale(1./n_integral_mc);
    npu_dat->Scale(1./n_integral_da);

    for (int i=0; i<60; i++){
        mc_vect.push_back(npu_mc->GetBinContent(i+1));
        data_vect.push_back(npu_dat->GetBinContent(i+1));
    }
  
    LumiWeights = reweight::LumiReWeighting(mc_vect,data_vect);
  
}
//------------------------------------------------------------------------------//
//-------------Set cross sections--------------------------------------------------------//
//-----------------------------------------------------------------------------//

void CommPlotProducer::SetInfo(TString generator, bool qcd, int TeV){
    gentype = generator;
    qcdtype = qcd;
    sqrtstev = TeV;
}

void CommPlotProducer::SetXS(){

    if (gentype=="pythia" || gentype=="herwig") {
        SetXS(gentype,qcdtype,sqrtstev);  // use what has been defined with SetInfo()
    } 
    else {
        if (sqrtstev!=-1) {
            cout << " MC Info not recognized! " << endl;
        }
        cout << " Set Default values to Pythia inclusive 13 TeV " << endl;
        SetInfo("pythia",0,13);
        SetXS("pythia",0,13);   
    }

}

void CommPlotProducer::SetXS(TString generator, bool MuEnriched, int TeV){
    /*
      pythia Inclusive  8TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000        --> 10 fichiers
      pythia MuEnriched 8TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 12 fichies
      herwig Inclusive  8TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 11 fichiers
      herwig MuEnriched 8TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 12 fichies

      pythia Inclusive  7TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000;       --> 10 fichiers
      pythia MuEnriched 7TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-150, n150-plus                                                  --> 7 fichiers
    */


    // http://cms.cern.ch/iCMS/jsp/mcprod/admin/requestmanagement.jsp?pwg=BTV&campid=Summer12_DR53X
    //                              0  &     15-30
    //                          15-20        20-30         30-50         50-80      80-120       120-170      170-300      300-470     470-600    600-800   800-1000   10000         
    double pythia_xs8  [12]={         0.,  9.8828742E8,  66285328.0,   8148778.0,  1033680.0,    156293.3,   34138.15,    1759.549,   113.8791,      26.9921, 3.550036,     0.};  // PTHAT 0, 15-30
    double pythia_xs8MU[12]={2.73858e+06,   1.8655e+06,   806298.0,    176187.6,      40448,     7463.94,    2299.752,    151.8048,    11.79648,    2.690196,   0.368781,     0.0849078};
  
    double herwig_xs8  [12]={         0.,        7.9E8,     5.32E7,   6545700.0,   833630.0,    126490.0,    27935.0,      1461.0,      95.25,     22.73,     2.997,     0.665}; // PTHAT 0, 15-30
    double herwig_xs8MU[12]={1.56419e+06,      783598.,     334139,     81821.2,   13754.9,     3162.25,    751.452,     53.4726,    3.22898,   1.00467,     0.128871,  0.025935};
  
    double pythia_xs7  [12]={         0.,   815900000.0, 53120000.0,    6359000.0,  784300.0,    115100.0,     24260.0,       1168.0,     70.22,      15.55,   1.844,    0 };
    double pythia_xs7MU[12]={  1471168.0,     1224034.0,   578463.0,    144421.74,   29048.7,   4440.2215,  2837.6712,          0,        0.,          0.,     0.,     0.};  
    //remark for pythia_xs7MU, it is for PTHAT 120-150 & 150-inf (and not 120-170 and 170-300)

//    double pythia_xs13[12]   = { 0., 2.237E9, 1.615E8, 2.211E7, 3000114.3, 493200., 120300., 7475., 587.1, 167., 28.25, 0. };
//    double pythia_xs13MU[12] = { 1.576E9*0.0039, 6.753E8*0.0065, 1.643E8*0.00816, 2.181E7*0.01522, 2.999E6*0.02424, 493200.*0.0473, 120300.*0.0676, 7475.*0.0864, 587.1*0.1024, 167.*0.0996, 28.25*0.1033, 8.975*0.1097 };
    //RunIISpring15
    double pythia_xs13[12]   = { 0., 1.83741E9, 1.40932E8, 1.92043E7, 2762530., 471100., 117276., 7823., 648.2, 186.9, 32.293, 9.4183 };
    double pythia_xs13MU[12] = { 1.27319E9*0.003, 5.58528E8*0.0053, 1.39803E8*0.01182, 1.92225E7*0.02276, 2.758420E6*0.03844, 469797.*0.05362, 117989.*0.07335, 7820.25*0.10196, 645.528*0.12242, 187.109*0.13412, 32.3486*0.14552, 10.4305*0.15544 };

    cout << generator << " " << MuEnriched << " " << TeV << endl;
 
    if      (generator=="pythia" && !MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs8[i]; }        }
    else if (generator=="pythia" &&  MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs8MU[i]; }        }

    else if (generator=="herwig" && !MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=herwig_xs8[i]; }        }
    else if (generator=="herwig" &&  MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=herwig_xs8MU[i]; }        }

    else if (generator=="pythia" && !MuEnriched && TeV==7){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs7[i]; }        }
    else if (generator=="pythia" &&  MuEnriched && TeV==7){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs7MU[i]; }        }
  
    else if (generator=="pythia" && !MuEnriched && TeV==13){      for (int i=0; i<12; i++){ x_section[i]=pythia_xs13[i]; }      }
    else if (generator=="pythia" &&  MuEnriched && TeV==13){      for (int i=0; i<12; i++){ x_section[i]=pythia_xs13MU[i]; }      }

}

void CommPlotProducer::SetSumXS(){
    sum_xs=0.0;
  
    for (int i=0; i<12; i++) {
        if (nmc_evt_vect[i]>0)  sum_xs+=x_section[i];
    }

}


//-------------------------------------------------------------------------------------------//
//-------------Compute the event weight given the pthat-------------------------------------------//
//------------------------------------------------------------------------------------------//

double CommPlotProducer::GetEvtWeight(){

    WeightXS=0.0;
    float  nevt  =0.0;
    float  xs    =0.0;
  
  
    if (qcdtype==0) { // inclusive
        if ( pthat >=  15. && pthat <  30. ){
            nevt  =nmc_evt_vect[1];
            xs=x_section[1];
        } 
    }
    else {
        if ( pthat >=  15. && pthat <  20. ){
            nevt  =nmc_evt_vect[0];
            xs=x_section[0];
        } 
        else if ( pthat >=  20. && pthat <  30. ){
            nevt  =nmc_evt_vect[1];
            xs=x_section[1];
        } 
    }
     
    if ( pthat >=  30. && pthat <  50. ){
        nevt  =nmc_evt_vect[2];
        xs=x_section[2];
    }
    else if ( pthat >=  50. && pthat <  80. ){
        nevt  =nmc_evt_vect[3];
        xs=x_section[3];
    }     
    else if ( pthat >=  80. && pthat <  120. ){
        nevt  =nmc_evt_vect[4];
        xs=x_section[4];
    } 
    if (gentype=="pythia" && qcdtype==1 && sqrtstev==7) {
        if ( pthat >=  120. && pthat <  150. ){
            nevt  =nmc_evt_vect[5];
            xs=x_section[5];
        } 
        else if ( pthat >=  150. ){
            nevt  =nmc_evt_vect[6];
            xs=x_section[6];
        }   
        else {
            nevt=0;
            xs=0;
        }
    }
    else {
        if ( pthat >=  120. && pthat <  170. ){
            nevt  =nmc_evt_vect[5];
            xs=x_section[5];
        }
        else if ( pthat >=  170. && pthat <  300. ){
            nevt  =nmc_evt_vect[6];
            xs=x_section[6];
        }    
        else if ( pthat >=  300. && pthat <  470. ){
            nevt  =nmc_evt_vect[7];
            xs=x_section[7];
        } 
        else if ( pthat >=  470. && pthat <  600. ){
            nevt  =nmc_evt_vect[8];
            xs=x_section[8];
        }
        else if ( pthat >=  600. && pthat <  800. ){
            nevt  =nmc_evt_vect[9];
            xs=x_section[9];
        }
        else if ( pthat >=  800. && pthat <  1000. ){
            nevt  =nmc_evt_vect[10];
            xs=x_section[10];
        }
        else if ( pthat >=  1000. ){
            nevt  =nmc_evt_vect[11];
            xs=x_section[11];
        }
    }     
  
    if (nevt>0) { WeightXS =xs/(sum_xs*nevt);  }
    else { WeightXS=0.; }
  
    if ( pthat <1. ) WeightXS=1. ;  // data

    return  WeightXS;
   
}

void CommPlotProducer::Loop(TString trigname, int trigger, float PtMin_Cut, float PtMax_Cut, TString output_name, TString puRewFileName, bool officalRecipe) { 

    // Prepare for PU reweighting 
    // if  PUreweightingFileName == "" --> no PU reweighting
    int nBinForPuRew = 0;
    TH1F* pu_weight_th1; 
    if (!isData && puRewFileName != ""){
        cout << "Start PU reweighting" << endl;
        if (officalRecipe){

            TFile *filepuest = new TFile(puRewFileName,"READ");
            TH1F* npu_dat = (TH1F*) filepuest->Get("pileup");
            nBinForPuRew = npu_dat->GetNbinsX(); 
            float n_integral_da = npu_dat->Integral();
            npu_dat->Scale(1./n_integral_da);

            static const Double_t arr_Spring2016[] = { 
                0.000829312873542, 0.00124276120498, 0.00339329181587, 0.00408224735376, 0.00383036590008, 0.00659159288946, 0.00816022734493, 0.00943640833116, 0.0137777376066, 0.017059392038, 
                0.0213193035468, 0.0247343174676, 0.0280848773878, 0.0323308476564, 0.0370394341409, 0.0456917721191, 0.0558762890594, 0.0576956187107, 0.0625325287017, 
                0.0591603758776, 0.0656650815128, 0.0678329011676, 0.0625142146389, 0.0548068448797, 0.0503893295063, 0.040209818868, 0.0374446988111, 0.0299661572042, 
                0.0272024759921, 0.0219328403791, 0.0179586571619, 0.0142926728247, 0.00839941654725, 0.00522366397213, 0.00224457976761, 0.000779274977993, 0.000197066585944,
                7.16031761328e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            vector<Double_t> Spring2016 (arr_Spring2016, arr_Spring2016 + sizeof(arr_Spring2016)/sizeof(arr_Spring2016[0]));


            TH1F* npu_mc = new TH1F("npu_mc", "npu_mc", nBinForPuRew, 0, nBinForPuRew);
            for( int i=0; i<nBinForPuRew; ++i) {
                npu_mc->SetBinContent(i+1, Spring2016.at(i));
            }
            float n_integral_mc = npu_mc->Integral();
            npu_mc->Scale(1./n_integral_mc);
            pu_weight_th1 = new TH1F("pu_weight_th1", "pu_weight_th1", nBinForPuRew, 0, nBinForPuRew);
            pu_weight_th1 = npu_dat;
            pu_weight_th1->Divide(npu_mc);

        }
        else {
            vector<float> mc_vect;
            vector<float> data_vect;
            puweight=false;
          
            TFile * puFile = new TFile(puRewFileName);

            TH1D* npu_mc= (TH1D*) puFile->Get("nPV_mc_unw");
            TH1D* npu_dat= (TH1D*) puFile->Get("nPV_data");
            cout << "Got file PU histo" << endl;
            
            float n_integral_mc = npu_mc->Integral();
            float n_integral_da = npu_dat->Integral();
            npu_mc->Scale(1./n_integral_mc);
            npu_dat->Scale(1./n_integral_da);

            for (int i=0; i<60; i++){
                mc_vect.push_back(npu_mc->GetBinContent(i+1));
                data_vect.push_back(npu_dat->GetBinContent(i+1));
            }
          
            LumiWeights = reweight::LumiReWeighting(mc_vect,data_vect);
        }
    }   

    //---------------Configuration-----------------------------------------// 
    //produceJetProbaTree=true;
    //produceRecoMuon = false;
    //produceTagVarCSVTree = true;
    //produceTagVarTree = true;
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
  
    int nHitAllCut, nHitPixCut;
    bool oldTrackSel = false;
    if (oldTrackSel) {
        nHitAllCut = 8;
        nHitPixCut = 2;
    }
    else {
        nHitAllCut = 0;
        nHitPixCut = 1;
    }


  
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
    TH1D* nPU_mc                  = new TH1D("nPU_mc",                "nPU_mc",                50,0,50);
    TH1D* nPV_data                = new TH1D("nPV_data",              "nPV_data",              50,0,50);
    TH1D* nPV_mc                  = new TH1D("nPV_mc",                "nPV_mc",                50,0,50);
    TH1D* nPV_mc_unw              = new TH1D("nPV_mc_unw",            "nPV_mc_unw",                50,0,50);
    TH1D* pt_hat                  = new TH1D("pt_hat",                "pt_hat",                100,0,1000);
    TH1D* pt_hat_fin              = new TH1D("pt_hat_fin",            "pt_hat_fin",             100,0,1000);
    TH1D* jet_pt_mc               = new TH1D("jet_pt_mc",  	    "jet_pt_mc", 	     PtMax/10,0,PtMax );
    TH1D* nJet_data               = new TH1D("nJet_data",             "nJet_data",             10,0,10    );
    TH1D* nJet_mc                 = new TH1D("nJet_mc",               "nJet_mc",               10,0,10    );
    TH1D* nJet_mc_inc             = new TH1D("nJet_mc_inc",           "nJet_mc_inc",               10,0,10    );
    TH1D* PVz_mc                  = new TH1D("PVz_mc",                "PVz_mc",                120,-30.,30.);
    TH1D* dPVz_mc                 = new TH1D("dPVz_mc",                "dPVz_mc",              50,0.,0.005);
    TH1D* nJet_csvLoose_data               = new TH1D("nJet_csvLoose_data",             "nJet_csvLoose_data",             10,0,10    );
    TH1D* nJet_csvMedium_data               = new TH1D("nJet_csvMedium_data",             "nJet_csvMedium_data",             10,0,10    );
    TH1D* nJet_csvTight_data               = new TH1D("nJet_csvTight_data",             "nJet_csvTight_data",             10,0,10    );
    TH1D* nJet_csvLoose_mc               = new TH1D("nJet_csvLoose_mc",             "nJet_csvLoose_mc",             10,0,10    );
    TH1D* nJet_csvMedium_mc               = new TH1D("nJet_csvMedium_mc",             "nJet_csvMedium_mc",             10,0,10    );
    TH1D* nJet_csvTight_mc               = new TH1D("nJet_csvTight_mc",             "nJet_csvTight_mc",             10,0,10    );


    nPU_mc->Sumw2();
    nPV_data->Sumw2();
    nPV_mc->Sumw2();
    nPV_mc_unw->Sumw2();
    pt_hat->Sumw2();
    pt_hat_fin->Sumw2();
    jet_pt_mc->Sumw2();
    nJet_data->Sumw2();
    nJet_mc->Sumw2();
    nJet_mc_inc->Sumw2();
     PVz_mc->Sumw2();
    dPVz_mc->Sumw2();
    nJet_csvLoose_data->Sumw2();
    nJet_csvMedium_data->Sumw2();
    nJet_csvTight_data->Sumw2();
    nJet_csvLoose_mc->Sumw2();
    nJet_csvMedium_mc->Sumw2();
    nJet_csvTight_mc->Sumw2();
  
    // --------------------------------------Histograms declaration -----------------------------------------//
    //  jets before cuts
    AddHisto("incjet_ptgen"  ,"pTgen of jets",		     100,0,1000 );
    AddHisto("incjet_pt"  ,"pT of jets",		     100,0,1000 );
    AddHisto("incjet_eta"  ,"eta of all jets",		     20,-2.5,2.5);
    AddHisto("incjet_phi"   ,"phi of all jets",                20,-1.*pi,pi);
    AddHisto("incjet_diffpt"  ,"abs(pTgen-pTreco) of jets",		     100,0,100 );
    AddHisto("incjet_diffrel"  ,"(pTgen-pTreco)/pTreco of jets",		     80,-2.,2. );
  
    
    // jets after cuts
    AddHisto("jet_ptgen"	  ,"pTgen of jets",		     PtMax/10,0,PtMax );
    AddHisto("jet_pt_all"	  ,"pT of jets",		     PtMax/10,0,PtMax );
    AddHisto("jet_pt_sv"	  ,"pT of jets containing a SV",     PtMax/10,0,PtMax);
    AddHisto("jet_eta"	  ,"eta of all jets",		     20,-2.5,2.5);
    AddHisto("jet_phi"    ,"phi of all jets",                20,-1.*pi,pi);
    AddHisto("jet_diffpt"  ,"abs(pTgen-pTreco) of jets",		     100,0,100 );
    AddHisto("jet_ratiopt"  ,"pTreco/pTgen of jets",		     150,0,3. );
    AddHisto("jet_diffrel"  ,"(pTgen-pTreco)/pTreco of jets",		     80,-2,2 );
    AddHisto("jet_flav"  ,"flavour of jets",		     28,-5.5, 22.5);
    AddHisto("jet_phi_sv"	  ,"phi of jets containing a SV",     50,-1*pi,pi);
    AddHisto("jet_eta_sv"	  ,"eta of jets containing a SV",     50,-1*pi,pi);
    AddHisto("jet_phi_CSVv2Loose"	  ,"phi of CSVv2 Loose jets",     50,-1*pi,pi);
    AddHisto("jet_eta_CSVv2Loose"	  ,"eta of CSVv2 Loose jets",     50,-1*pi,pi);
    AddHisto("jet_phi_CSVv2Medium"	  ,"phi of CSVv2 Medium jets",     50,-1*pi,pi);
    AddHisto("jet_eta_CSVv2Medium"	  ,"eta of CSVv2 Medium jets",     50,-1*pi,pi);
    AddHisto("jet_phi_CSVv2Tight"	  ,"phi of CSVv2 Tight jets",     50,-1*pi,pi);
    AddHisto("jet_eta_CSVv2Tight"	  ,"eta of CSVv2 Tight jets",     50,-1*pi,pi);

    AddHisto("jet_pt_csvv2_l"       ,"pT of jets",              PtMax/10,0,PtMax );
    AddHisto("jet_pt_csvv2_m"       ,"pT of jets",              PtMax/10,0,PtMax );
    AddHisto("jet_pt_csvv2_t"       ,"pT of jets",              PtMax/10,0,PtMax );
    AddHisto("jet_pt_csvl"       ,"pTgen of jets",                  PtMax/10,0,PtMax );
    AddHisto("jet_pt_csvm"       ,"pTgen of jets",                  PtMax/10,0,PtMax );
    AddHisto("jet_pt_csvivfl"    ,"pTgen of jets",                  PtMax/10,0,PtMax );
    AddHisto("jet_pt_csvivfm"    ,"pTgen of jets",                  PtMax/10,0,PtMax );
    AddHisto("jet_pt_jpl"        ,"pTgen of jets",                  PtMax/10,0,PtMax );
    AddHisto("jet_pt_tchpt"      ,"pTgen of jets",                  PtMax/10,0,PtMax );


  
    int nMubins=20;
    if (produceRecoMuon) {
      AddHisto("muon_multi"   ,      "number of muons",	   7,-0.5,6.5    );
      AddHisto("muon_ptrel"	  ,      "pT rel. of the muon",	   nMubins,0,5        );
      AddHisto("muon_chi2"	  ,      "norm. chi2 of the muon", nMubins,0,10       );  
      AddHisto("muon_Pt",		 "Muon p_{T}",  	   50, 0, 100	 );
      AddHisto("muon_eta",		 "Muon #eta",  	           nMubins, -2.5, 2.5 );  
      AddHisto("muon_phi",                 "Muon #phi",              nMubins, -1.*pi,pi);
      AddHisto("muon_Ip3d", 	 "Muon 3D IP",  	   nMubins, -0.1, 0.1 );
      AddHisto("muon_Ip2d", 	 "Muon 2D IP",  	   nMubins, -0.1, 0.1 );
      AddHisto("muon_Sip3d",	 "Muon 3D IP significance",nMubins, -35, 35   );
      AddHisto("muon_Sip2d",	 "Muon 2D IP significance",nMubins, -35, 35   );
      AddHisto("muon_DeltaR",	 "Muon deltaR",nMubins,0,0.4); 
    }
  
    int nSVbins = 20; 
    AddHisto("sv_multi_0"	  ,      "number of secondary vertex",                          6,-0.5,5.5   );
    AddHisto("sv_multi"	  ,      "number of secondary vertex",                          6,-0.5,5.5   );
    if (produceSVinfo) {
    AddHisto("sv_deltaR_jet",      "sv_deltaR_jet",                                   nSVbins,0.,0.4    );
    AddHisto("sv_deltaR_sumJet",   "SVvtxSumJetDeltaR",                                   nSVbins,0.,0.4    );
    AddHisto("sv_deltaR_sumDir",   "SVvtxSumVtxDirDeltaR",                                nSVbins,0.,0.4    );
    AddHisto("sv_pt",              "Vtx p_{T}",                                           nSVbins,0.,100.   );
    AddHisto("sv_eta",             "Vtx #eta",                                            nSVbins, -2.5, 2.5);
    AddHisto("sv_phi",             "Vtx #phi",                                            nSVbins, -1.*pi,pi);
    AddHisto("sv_flightSig2D",     "Flight distance significance 2D",                     nSVbins,0.,80.    );
    AddHisto("sv_flight2D",        "Flight distance 2D",                                  nSVbins,0.,2.5    );
    AddHisto("sv_flight2D_3trk",   "Flight distance 2D with >=3 tracks",                  nSVbins,0.,2.5    );
    AddHisto("sv_flight2DSig_3trk","Flight distance significance 2D with >=3 tracks",     nSVbins,0.,80.    );
    AddHisto("sv_flight3D",        "Flight distance 3D",                                  nSVbins,0.,15.    );  
    AddHisto("sv_flight3DSig" ,    "flight distance significance 3D",	                nSVbins,0.,80.    );
    AddHisto("sv_mass"	  ,      "invariant mass of the secondary vertex",              nSVbins,0.,8.     );
    AddHisto("sv_chi2norm"  ,      "normalized chi2 of the secondary vertex",             nSVbins,0.,10.    );
    AddHisto("svnTrk",	         "Track multiplicity : SVnVertexTracks (centered)",     13,-0.5,12.5 );
    AddHisto("sv_flight3Derr",     "Flight distance error 3D",                            nSVbins,0.,0.5);
    AddHisto("sv_flight2Derr",     "Flight distance error 2D",                            nSVbins,0.,0.2);
    AddHisto("sv_mass_3trk"	,"invariant mass of the secondary vertex with at least 3 SV tracks",  nSVbins,0.,8.     );
    }
  
    int nTrackbins = 100; 
    AddHisto("track_multi"  ,      "number of tracks in the jets",                40,-0.5,39.5 );
    AddHisto("trk_multi_sel"  ,    "number of selected tracks in the jets",       30,-0.5,29.5 );
    AddHisto("trk_multi_sel03"  ,    "number of selected tracks in the jets",       30,-0.5,29.5 );

    if (produceJetProbaTree) {
    AddHisto("track_multi04"  ,      "number of tracks in the jets04",                40,-0.5,39.5 );
    AddHisto("track_dr"   ,      "Delta R (jet, track)",               40,0.,0.4   );
    AddHisto("track_chi2"   ,      "normalized chi2 of the tracks",               nTrackbins,0.,5.   );
    AddHisto("track_nHit" ,        "number of hits ",               35,0, 35);
    AddHisto("track_HPix"   ,      "number of hits in the Pixel",                 10,0,10  );
    AddHisto("track_nHitStrip"   ,      "number of hits in the Strip",                 27, 0, 27  );
    AddHisto("track_nHitTIB"   ,      "number of hits in the TIB",                 12,0, 12  );
    AddHisto("track_nHitTID"   ,      "number of hits in the TID",                 10,0, 10  );
    AddHisto("track_nHitTOB"   ,      "number of hits in the TOB",                 16,0, 16  );
    AddHisto("track_nHitTEC"   ,      "number of hits in the TEC",                 21,0, 21  );
    AddHisto("track_nHitPXB"   ,      "number of hits in the PXB",                 7,0, 7  );
    AddHisto("track_nHitPXF"   ,      "number of hits in the PXF",                 7,0, 7  );
  
    AddHisto("track_IPs"    ,      "3D IP significance of all tracks",	        nTrackbins,-35.,35.  );
    AddHisto("track_IPs1tr" ,      "3D IP significance of the first track",       nTrackbins,-35.,35.  );
    AddHisto("track_IPs2tr" ,      "3D IP significance of the second track",      nTrackbins,-35.,35.  );
    AddHisto("track_IP"     ,      "3D IP of all tracks",	                        nTrackbins,-0.1,0.1);
    AddHisto("track_IP1tr"  ,      "3D IP of the first track",                    nTrackbins,-0.1,0.1);
    AddHisto("track_IP2tr"  ,      "3D IP of the second track",                   nTrackbins,-0.1,0.1); 
    AddHisto("track_IP2Ds"    ,    "2D IP significance of all tracks",	        nTrackbins,-35.,35.  );
    AddHisto("track_IP2Ds1tr" ,    "2D IP significance of the first track",       nTrackbins,-35.,35.  );
    AddHisto("track_IP2Ds2tr" ,    "2D IP significance of the second track",      nTrackbins,-35.,35.  );
    AddHisto("track_IP2D"    ,     "2D IP of all tracks",	                        nTrackbins,-0.1,0.1);
    AddHisto("track_IP2D1tr" ,     "2D IP of the first track",                    nTrackbins,-0.1,0.1);
    AddHisto("track_IP2D2tr" ,     "2D IP of the second track",                   nTrackbins,-0.1,0.1);
    AddHisto("track_IP2Derr1tr" ,  "2D IP error of the first track",              nTrackbins,0,0.1);    
    AddHisto("track_IPerr1tr"   ,  "3D IP error of the first track",              nTrackbins,0,0.1);  
    AddHisto("track_IP2Derr2tr" ,  "2D IP error of the second track",             nTrackbins,0,0.1);    
    AddHisto("track_IPerr2tr"   ,  "3D IP error of the second track",             nTrackbins,0,0.1); 
    AddHisto("track_IP2Derr" ,     "2D IP error",                                 nTrackbins,0,0.1);    
    AddHisto("track_IPerr"   ,     "3D IP error",                                 nTrackbins,0,0.1); 
    AddHisto("track_IPs3tr" ,      "3D IP significance of the third track",       nTrackbins,-35.,35.  );
    AddHisto("track_IP3tr"  ,      "3D IP of the third track",                    nTrackbins,-0.1,0.1);
    AddHisto("track_IPerr3tr"   ,  "3D IP error of the third track",              nTrackbins,0,0.1); 
    AddHisto("track_IP2Ds3tr" ,    "2D IP significance of the second track",      nTrackbins,-35.,35.  );
    AddHisto("track_IP2D3tr" ,     "2D IP of the third track",                    nTrackbins,-0.1,0.1);
    AddHisto("track_IP2Derr3tr" ,  "2D IP error of the third track",              nTrackbins,0,0.1);   
   
    AddHisto("track_len"     ,     "decay length",		                nTrackbins,0,5.     );
    AddHisto("track_dist"    ,     "distance to the jet axis",                    nTrackbins,0.,0.08    );
    AddHisto("track_dz"     ,     "transverse IP",                               nTrackbins,-3,3  );  
    AddHisto("track_isfromSV",     "Track is from SV",                            2,-0.5, 1.5   );  
    AddHisto("track_pt"	  ,      "pT of all the tracks",	                80,0.,200.    );
    AddHisto("track_pt15"   ,      "pT of all the tracks",                        150,0.,15.    );
    AddHisto("track_chi2_cut"     ,"normalized chi2 ",  	                        nTrackbins,0.,20.    );
    AddHisto("track_nHit_cut"   ,"number of hits  ",               35,-0.5, 34.5 );
    AddHisto("track_HPix_cut"     ,"number of hits in the Pixel ",                 10,-0.5, 9.5  );
    AddHisto("track_len_cut"      ,"decay length ",		                nTrackbins,0,25.     );
    AddHisto("track_dist_cut"     ,"distance to the jet axis ",                    nTrackbins,0.,0.3    );
    AddHisto("track_dz_cut"      ,"transverse IP ",		                nTrackbins,-20., 20.  );  
    AddHisto("track_pt_cut"       ,"pT ",	                                        80,0.,200.);
    AddHisto("track_pt15_cut"     , "pT of all the tracks",                       150,0.,15.    );
    AddHisto("track_IP2D_cut"     ,"IP2D ",	                                nTrackbins,-1.,1.);
    AddHisto2D("track_nHit_vs_eta"     ,"number of hits  vs eta",               35,-0.5, 34.5 ,50, -2.5, 2.5);
    AddHisto2D("track_HPix_vs_eta"     ,"number of hits in the Pixel vs eta",                 10,-0.5, 9.5, 50,  -1.*pi,pi);
    AddHisto2D("track_nHit_vs_phi"     ,"number of hits  vs phi",               35,-0.5, 34.5 ,50, -2.5, 2.5);
    AddHisto2D("track_HPix_vs_phi"     ,"number of hits in the Pixel vs phi",                 10,-0.5, 9.5, 50,  -1.*pi,pi);
    AddHisto2D("track_eta_vs_phi"     ,"eta vs phi",              50, -2.5, 2.5, 50,  -1.*pi,pi);

    AddHisto("track_len_sel_zoom"     ,     "decay length",		                nTrackbins,0,0.25     );
    AddHisto("track_len_all_zoom"     ,     "decay length",		                nTrackbins,0,0.25     );
    }

    int nDiscbins = 50; 
   
    AddHisto("TCHE"	  ,"TCHE",				     nDiscbins,0.,30. );
    AddHisto("TCHP"	  ,"TCHP",				     nDiscbins,0.,30. );  
    AddHisto("JP" 	  ,"JP",				     nDiscbins,0.,1.5 );
    AddHisto("JBP"	  ,"JBP",				     nDiscbins,0.,4.  );
    AddHisto("SSV"	  ,"SSVHE",				     70,0.,7.  );
    AddHisto("SSVHP"  ,"SSVHP",				             70,0.,7.  ); 
    AddHisto("CSV"	  ,"CSV",				     nDiscbins,0.,1.  );
    AddHisto("CSVIVF" ,"CSVIVF",                 nDiscbins,0.,1.  );
    AddHisto("cMVA" ,"cMVA",                 nDiscbins,0.,1.  );
    AddHisto("cMVAv2" ,"cMVAv2",                 nDiscbins,-1.,1.);
    AddHisto("cTagCvsB" ,"cTagCvsB",                 nDiscbins,-1.,1.);
    AddHisto("cTagCvsBN" ,"cTagCvsBN",                 nDiscbins,-1.,1.);
    AddHisto("cTagCvsBP" ,"cTagCvsBP",                 nDiscbins,-1.,1.);
    AddHisto("cTagCvsL" ,"cTagCvsL",                 nDiscbins,-1.,1.);
    AddHisto("cTagCvsLN" ,"cTagCvsLN",                 nDiscbins,-1.,1.);
    AddHisto("cTagCvsLP" ,"cTagCvsLP",                 nDiscbins,-1.,1.);

    AddHisto("TCHE_extended1"	  ,"TCHE_extended1",				     60, -30.,30. );
    AddHisto("TCHP_extended1"	  ,"TCHP_extended1",				     60, -30.,30. );
    AddHisto("TCHE_extended2"	  ,"TCHE_extended2",				     100,-30.,30. );
    AddHisto("TCHP_extended2"	  ,"TCHP_extended2",				     100,-30.,30. );
    AddHisto("discri_ssche0",	 "SSVHE Discriminator",    80, -1., 7.   ); 
    AddHisto("discri_sschp0",	 "SSVHP Discriminator",    80, -1., 7.   ); 

    AddHisto("SoftMu"	  ,"SoftMu",				     nDiscbins,0.,1.  );
    AddHisto("SoftEl"	  ,"SoftEl",				     nDiscbins,0.,1.  );
  
    AddHisto("pfmuon_multi",      "number of pfmuons",	   7,-0.5,6.5    );
    AddHisto("pfmuon_goodquality","quality of pfmuons",3,-0.5,2.5   );
    AddHisto("pfmuon_pt",		"pfmuon p_{T}",  	   50, 0, 100	 );
    AddHisto("pfmuon_eta",   	"pfmuon #eta",  	  nMubins, -2.5, 2.5 );  
    AddHisto("pfmuon_phi",        "pfmuon #phi",              nMubins, -1.*pi,pi);
    AddHisto("pfmuon_ip",	"3D IP pfmuon",50, -50., 50.   );
    AddHisto("pfmuon_ptrel",      "pT rel. of the muon",	   nMubins,0,5        );
    AddHisto("pfmuon_ratio",      "ratio of pfmuon", nMubins,0,1       );  
    AddHisto("pfmuon_ratiorel",   "ratioRel of pfmuon", nMubins,0,0.03       );  
    AddHisto("pfmuon_deltar",	"#DeltaR(pfmuon,jet)",nMubins,0,0.4);
  
    AddHisto("pfelectron_multi",      "number of pfelectron",	   7,-0.5,6.5    );
    AddHisto("pfelectron_pt",		"pfelectron p_{T}",  	   50, 0, 100	 );
    AddHisto("pfelectron_eta",   	"pfelectron #eta",  	           nMubins, -2.5, 2.5 );  
    AddHisto("pfelectron_phi",        "pfelectron #phi",              nMubins, -1.*pi,pi);
    AddHisto("pfelectron_ip",	"3D IP of pfelectron",50, -50., 50.   );
    AddHisto("pfelectron_ptrel",      "pT rel. of the pfelectron",	   nMubins,0,5        );
    AddHisto("pfelectron_ratio",      "ratio of pfelectron", nMubins,0,1       );  
    AddHisto("pfelectron_ratiorel",   "ratioRel of pfelectron", nMubins,0,0.03      );  
    AddHisto("pfelectron_deltar",	"#DeltaR(pfelectron,jet)",nMubins,0,0.5);  

    AddHisto("pfelectron_eta_1",  "pfelectron #eta",                 nMubins, -2.5, 2.5 );
    AddHisto("pfelectron_jeta_1", "jet #eta",                 nMubins, -2.5, 2.5 );
    AddHisto("pfelectron_ip_1",	 "3D IP of pfelectron",140, -350., 350.   );
    AddHisto("pfelectron_ip2d_1", "2D IP of pfelectron",140, -350., 350.   );	
    AddHisto("pfelectron_ptratio_1", "el pT / jet pT", 55, 0., 1.1); 
    AddHisto("pfelectron_ratio_1",    "ratio of pfelectron", nMubins,0,1       );
    AddHisto("pfelectron_ratiorel_1",  "ratioRel of pfelectron", nMubins,0,0.03      );
    AddHisto("pfelectron_deltar_1", "#DeltaR(pfelectron,jet)",nMubins,0,0.5);
    AddHisto("pfelectron_eta_2",  "pfelectron #eta",                 nMubins, -2.5, 2.5 );
    AddHisto("pfelectron_jeta_2", "jet #eta",                 nMubins, -2.5, 2.5 );
    AddHisto("pfelectron_ip_2",	 "3D IP of pfelectron",140, -350., 350.   );
    AddHisto("pfelectron_ip2d_2", "2D IP of pfelectron",140, -350., 350.   );	
    AddHisto("pfelectron_ptratio_2", "el pT / jet pT", 55, 0., 1.1); 
    AddHisto("pfelectron_ratio_2",    "ratio of pfelectron", nMubins,0,1       );
    AddHisto("pfelectron_ratiorel_2",  "ratioRel of pfelectron", nMubins,0,0.03      );
    AddHisto("pfelectron_deltar_2", "#DeltaR(pfelectron,jet)",nMubins,0,0.5);

    if (produceTagVarTree) {
      AddHisto("tagvar_vertexNTracks",  "# SV tracks",  13,-0.5,12.5 );
      AddHisto("tagvar_vertexmass",     "SV mass", nSVbins,0.,8.     );
      AddHisto("tagvar_vertexmass_cat0",     "SV mass", nSVbins,0.,8.     );
      AddHisto("tagvar_vertexmass3trk", "SV mass (at least 3 SV tracks)",nSVbins,0.,8.     );
      AddHisto("tagvar_vertexJetDeltaR","DeltaR(SV,jet) ",               nSVbins,0.,0.4    );
    }
    if (produceTagVarCSVTree) { 
      AddHisto("tagvarCSV_vertexCategory","vertex category",4,-0.5, 3.5);
      AddHisto("tagvarCSV_Sig2dAboveCharm",   "IP significance 2D charm",                            nSVbins,-35.,35.  );
      AddHisto("tagvarCSV_trackEtaRel",   "Track etaRel",100,0.,8.);
      AddHisto("tagvarCSV_vertexmass",    "SV mass", nSVbins,0.,8.     );
      AddHisto("tagvarCSV_vertexmass_cat0",    "SV mass", nSVbins,0.,8.     );
      AddHisto("tagvarCSV_vertexmass3trk",    "SV mass (at least 3 SV tracks)", nSVbins,0.,8.     );
      AddHisto("tagvarCSV_vertexmass3trk_cat0",    "SV mass (at least 3 SV tracks)", nSVbins,0.,8.     );
      AddHisto("tagvarCSV_vertexNTracks",  "# SV tracks",  13,-0.5,12.5 );
      AddHisto("tagvarCSV_vertexNTracks_cat0",  "# SV tracks",  13,-0.5,12.5 );
      AddHisto("tagvarCSV_energyratio","Fractional energy",                                   nSVbins,0.,1.     ); 
      AddHisto("tagvarCSV_trackSip3dSig","3D IP significance",       nTrackbins,-35.,35.  );
      AddHisto("tagvarCSV_2DsigFlightDist",     "Flight distance significance 2D",                     nSVbins,0.,80.    );
      AddHisto("tagvarCSV_2DsigFlightDist_cat0",     "Flight distance significance 2D",                     nSVbins,0.,80.    );
      AddHisto("tagvarCSV_vertexJetDeltaR","DeltaR(SV,jet) ",               nSVbins,0.,0.4    );
      AddHisto("tagvarCSV_vertexJetDeltaR_cat0","DeltaR(SV,jet) ",               nSVbins,0.,0.4    );
    }
  
 
    if (produceJetProbaTree) { 
    AddHisto2D("seltrack_vs_jetpt", "sel track multiplicity vs jet pt",         PtMax/20,0,PtMax, 100,-0.5,99.5);
    }

    if (produceSVinfo) {
    AddHisto2D("sv_mass_vs_flightDist3D", " SVMass vs SV 3D flight distance ",  50,0, 10,60,0,6);			
    AddHisto2D("avg_sv_mass_vs_jetpt","Avg SVMass vs jet pt",                   PtMax/20,0,PtMax, 100,0,6);
    AddHisto2D("sv_deltar_jet_vs_jetpt","SVJetDeltaR vs jet pt",                PtMax/20,0,PtMax, 40,0.,0.4);  
    AddHisto2D("sv_deltar_sum_jet_vs_jetpt","SVvtxSumJetDeltaR vs jet pt",      PtMax/20,0,PtMax, 40,0.,0.4);
    AddHisto2D("sv_deltar_sum_dir_vs_jetpt","SVvtxSumVtxDirDeltaR vs jet pt",   PtMax/20,0,PtMax, 40,0.,0.4); 
    }
    if (produceRecoMuon) {
    AddHisto2D("muon_ptrel_vs_jetpt","Muon_p{T}^{rel} vs jet pt",               PtMax/20,0,PtMax,50,0,5);  
    AddHisto2D("muon_DeltaR_vs_jetpt","Muon1 DeltaR vs jet pt",                 PtMax/20,0,PtMax,40,0,0.4);
    }
    
    Nevent = 0;
    int Ntrigevent = 0;
    if (fChain == 0){
        cout << "fchain = 0 why?\n";
        return;
    }
  
    Long64_t nentries = fChain->GetEntries();
  
    Long64_t nbytes = 0, nb = 0;

    //------------------------------------------------------------------------------------------------------------------//  
    //----------------------------------------EVENT LOOP ---------------------------------------------------------------// 
    //------------------------------------------------------------------------------------------------------------------//  
    
    bool printinfo = false;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        if (printinfo) cout << " loading tree\n";
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if (printinfo) cout << " getting entry\n";
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //-----------------------------------
        //is data or MC ?
        //-----------------------------------
        if(pthat == -1) isData=true;
        else            isData=false;
    
        if (isData) N_event_data_before_sel++;
        else        N_event_mc_before_sel++;
        //-----------------------------------
        //initialaze the weight at 1
        //-----------------------------------
        double ww=1; 

       
        //-----------------------------------//
        //-------xs reweighting--------------//
        //-----------------------------------//
        if (printinfo) cout << " getting event weight\n";
        ww=GetEvtWeight();
        if (sqrtstev<0) ww=1.;
        if (jentry<10) cout << " Evt weight for cross section : " << ww << endl;
        

        double ww_woPU=ww;


    
        //--------------------------------------------//  
        //-------------pile-up reweighting------------//  
        //--------------------------------------------//  
        if (!isData && puRewFileName != ""){
            float WeightPU =1;
            if (officalRecipe) {
                if (nPUtrue > nBinForPuRew-1) WeightPU = 0;
                else WeightPU = pu_weight_th1->GetBinContent(pu_weight_th1->GetBin(nPUtrue));
            }
            else WeightPU  = LumiWeights.weight( nPV );
            ww*=WeightPU;
        }    
        if (!isData)  pt_hat         ->Fill(pthat,ww);
    
        //-----------------------------------
        //counter of events
        //-----------------------------------
        Nevent++;
        if(Nevent%100000 ==0 && Nevent!=0) cout << " number of processed events is " << Nevent <<  " = " << (100.*Nevent)/(1.*nentries) << "%" <<endl;
        int njet_c     =0;
        int njet_l     =0;  
        int njet_mc    =0; 
        int njet_data  =0; 
        int njet_csvLoose_data = 0;    
        int njet_csvMedium_data = 0;    
        int njet_csvTight_data = 0;    
        int njet_csvLoose_mc = 0;    
        int njet_csvMedium_mc = 0;    
        int njet_csvTight_mc = 0;
    
        //at least 1 jet in the event
        if (nJet<=0) continue;

        //select good PV:
        //if (abs(PVz-GenPVz)>0.005) continue;

        //-----------------------------------
        //Apply trigger selection
        //-----------------------------------
        if (printinfo) cout << " checking trigger " << trigger << endl;
        bool isTrigOK = passTrigger(trigname, trigger);
        if (!isTrigOK) continue;
        ++Ntrigevent;

        //-----------------------------------
        //Fill control plot
        //-----------------------------------
    
        if(isData){
            N_event_data_after_sel++;
            nPV_data       ->Fill(nPV);
        }
        else{
            N_event_mc_after_sel++;
            nPU_mc         ->Fill(nPUtrue,ww);
            nPV_mc         ->Fill(nPV,ww);
            nPV_mc_unw     ->Fill(nPV,ww_woPU);
            nJet_mc_inc->Fill(nJet, ww);
        }

      
        //-----------------------------------
        //Loop on jets 
        //-----------------------------------
        if (printinfo) cout << " loop over jets\n";
        bool finselect = false;
        for (int ijet = 0; ijet < nJet; ijet++) {
      
            float ptjet    = Jet_pt[ijet];
            float etajet   = Jet_eta[ijet];
            float phijet   = Jet_phi[ijet];      
            float ntrkjet  = Jet_ntracks[ijet];  
            int   flav     = Jet_flavour[ijet];
            if (Jet_genpt[ijet]<8) flav=22;

            FillHisto_floatFromMap("incjet_ptgen",   flav, 0 ,Jet_genpt[ijet]    ,ww);
            FillHisto_floatFromMap("incjet_pt",  flav, 0 ,ptjet    ,ww);
            FillHisto_floatFromMap("incjet_eta",     flav, 0 ,etajet   ,ww);
            FillHisto_floatFromMap("incjet_phi",     flav, 0 ,phijet   ,ww);
            FillHisto_floatFromMap("incjet_diffpt",   flav, 0 ,abs(Jet_genpt[ijet]-ptjet)    ,ww);
            FillHisto_floatFromMap("incjet_diffrel",   flav, 0 ,(Jet_genpt[ijet]-ptjet)/ptjet    ,ww);

          
            if (   ptjet  < PtMin  || ptjet  > PtMax   ) continue;
            if (   fabs(etajet) > EtaCut               ) continue;

            int   idxFirstMuon = -1;
            int nselmuon = 0;
            int nmu = 0;
            if (usePFMuon) {
             if (nPFMuon>0){
                for (int imu=0; imu< nPFMuon; imu++){
                    if (PFMuon_IdxJet[imu]==ijet ){
                        nmu++;
                        if (passPFMuonSelection(imu)){
                            if(nselmuon == 0) {
                                idxFirstMuon = imu;
                            }
                            nselmuon++;
                        }
                    }
                }
              }
            }
            else if (produceRecoMuon) {
             if (nMuon>0){
                for (int imu=0; imu< nMuon; imu++){
                    if (Muon_IdxJet[imu]==ijet ){
                        nmu++;
                        if (passMuonSelection(imu, ijet)){
                            if(nselmuon == 0) {
                                idxFirstMuon = imu;
                            }
                            nselmuon++;
                        }
                    }
                }
             }
            }

            if (trigname=="btag" && nselmuon==0) continue;  // look at jet with muon for btag selection!

            finselect = true;

            int n_sv             =0.;
            float mass_sv        =-1.;
            float chi2norm_sv    =-1.;
            float flightSig_sv   =-1.;    
            float flight2DSig_sv =-1.;    
            float sv_dR_jet      =-1.;
            float sv_dR_dir_sum  =-1.; 
            float sv_dR_jet_sum  =-1.;
            float sv_pt	   =-1.;      
            float sveta         =-1000.; 
            float svphi         =-1000.; 
            float sv_flight3D    =-1.;
            float sv_flight3Derr =-1.;
            float sv_flight2D    =-1.;
            float sv_flight2Derr =-1.;
            int   sv_totchar     =-1.;
            float sv_nTrk        =0.;
      
      
            float tche     = Jet_Ip2P[ijet];
            float tchp     = Jet_Ip3P[ijet];
            float jetproba = Jet_Proba[ijet];
            float jetbproba= Jet_Bprob[ijet];
            float ssvhe    = Jet_Svx[ijet] ;
            float ssvhp    = Jet_SvxHP[ijet];
            float csv      = Jet_CombSvx[ijet];
            float csvivf   = Jet_CombIVF[ijet];
            float csvmva   = Jet_cMVA[ijet];
            float csvmvav2   = Jet_cMVAv2[ijet];
            float ctagCvsB = CTag_Jet_CvsB[ijet];
            float ctagCvsBN = CTag_Jet_CvsBN[ijet];
            float ctagCvsBP = CTag_Jet_CvsBP[ijet];
            float ctagCvsL = CTag_Jet_CvsL[ijet];
            float ctagCvsLN = CTag_Jet_CvsLN[ijet];
            float ctagCvsLP = CTag_Jet_CvsLP[ijet];
      
      
            bool isGluonSplit=false;
      
      
            //fill jet multiplicity
            if (!isData) {
                if (fabs(flav)==4)                  njet_c++;
                if (fabs(flav)<4 || fabs(flav)==21) njet_l++;
                njet_mc++;
                jet_pt_mc   ->Fill(ptjet,ww);
            }
            else njet_data++;
       
            // Gluon Splitting from Daniel
           int nsplit = 0;
           double dRqj;
           isGluonSplit = false;

           if ( fabs(flav) == 4 ) {
              for (int k = 0; k < nDHadrons; k++) {
                dRqj = DeltaR( Jet_eta[ijet], Jet_phi[ijet], DHadron_eta[k], DHadron_phi[k] );
                if ( dRqj < 0.4 ) nsplit++;
              }
           }
           if ( fabs(flav) == 5 ) {
               for (int k = 0; k < nBHadrons; k++) {
                 if ( BHadron_hasBdaughter[k] ) continue;
                 dRqj = DeltaR( Jet_eta[ijet], Jet_phi[ijet], BHadron_eta[k], BHadron_phi[k] );
                 if ( dRqj < 0.4 ) nsplit++;
               }
           }
           if ( nsplit > 1 ) isGluonSplit = true; // then you can fill or not the discriminator plots with isGluonSplit true

      
            FillHisto_floatFromMap("jet_ptgen",                 flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);
            FillHisto_floatFromMap("jet_pt_all",                 flav, isGluonSplit ,ptjet    ,ww);
            FillHisto_floatFromMap("jet_diffpt",   flav, isGluonSplit ,abs(Jet_genpt[ijet]-ptjet)    ,ww);
            if (Jet_genpt[ijet]>0) FillHisto_floatFromMap("jet_ratiopt",   flav, isGluonSplit ,ptjet/Jet_genpt[ijet]    ,ww);
            FillHisto_floatFromMap("jet_diffrel",   flav, isGluonSplit ,(Jet_genpt[ijet]-ptjet)/ptjet   ,ww);
            if (nSV > 0)
            {
                FillHisto_floatFromMap("jet_pt_sv",      flav, isGluonSplit ,ptjet    ,ww);
                FillHisto_floatFromMap("jet_phi_sv",      flav, isGluonSplit ,phijet    ,ww);
                FillHisto_floatFromMap("jet_eta_sv",      flav, isGluonSplit ,etajet    ,ww);
            }
      
            FillHisto_floatFromMap("jet_eta",     flav, isGluonSplit ,etajet   ,ww);
            FillHisto_floatFromMap("jet_phi",     flav, isGluonSplit ,phijet   ,ww);
            FillHisto_intFromMap(  "jet_flav", flav, isGluonSplit ,flav  ,ww);
            if (csvivf>0.460) {
                FillHisto_floatFromMap(  "jet_phi_CSVv2Loose", flav, isGluonSplit, phijet         ,   ww);
                FillHisto_floatFromMap(  "jet_eta_CSVv2Loose", flav, isGluonSplit, etajet         ,   ww);
                FillHisto_floatFromMap("jet_pt_csvv2_l",                 flav, isGluonSplit ,ptjet    ,ww);
                if(isData) njet_csvLoose_data++;
                else njet_csvLoose_mc++;
            }
            if (csvivf>0.800) {
                FillHisto_floatFromMap(  "jet_phi_CSVv2Medium", flav, isGluonSplit, phijet         ,   ww);
                FillHisto_floatFromMap(  "jet_eta_CSVv2Medium", flav, isGluonSplit, etajet         ,   ww);
                FillHisto_floatFromMap("jet_pt_csvv2_m",                 flav, isGluonSplit ,ptjet    ,ww);
                if(isData) njet_csvMedium_data++;
                else njet_csvMedium_mc++;
            }
            if (csvivf>0.935) {
                FillHisto_floatFromMap(  "jet_phi_CSVv2Tight", flav, isGluonSplit, phijet         ,   ww);
                FillHisto_floatFromMap(  "jet_eta_CSVv2Tight", flav, isGluonSplit, etajet         ,   ww);
                FillHisto_floatFromMap("jet_pt_csvv2_t",                 flav, isGluonSplit ,ptjet    ,ww);
                if(isData) njet_csvTight_data++;
                else njet_csvTight_mc++;
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
            int ntracksel03  =0;   
            int ntrkjet04 =0;
            ntrkjet=0;
         
            if ( produceJetProbaTree ) {
	
                for (int itrk=Jet_nFirstTrack[ijet]; itrk<Jet_nLastTrack[ijet] ; itrk++){
                    // temporary fix for different MiniAOD versions
                    // if (Track_pt[itrk]<0.95) continue;
                    ntrkjet++;

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
	    
                    if (Track_nHitAll[itrk] >= nHitAllCut)   passNhit=true;
                    if (Track_nHitPixel[itrk] >= nHitPixCut)   passPix= true;
                    if (fabs(Track_dz[itrk]) < 17)   passIPz=true;
                    if (Track_pt[itrk] > 1)           passPt=true;
                    if (Track_chi2[itrk] < 5)         passnormchi2=true;
                    if (fabs(Track_dist[itrk]) < 0.07)passtrkdist=true;
                    if (Track_length[itrk] < 5)       passtrklen=true;
                    if (fabs(Track_IP2D[itrk]) < 0.2)   passTrackIP2D=true;
	  
                    float dr_track= DeltaR( Jet_eta[ijet], Jet_phi[ijet], Track_eta[itrk], Track_phi[itrk] );
                    if (dr_track<0.4) ntrkjet04++;
                    if (!use_selected_tracks){
	    
                        FillHisto_floatFromMap("track_len_all_zoom",          flav, isGluonSplit ,Track_length[itrk] , ww);
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
                            FillHisto_floatFromMap("track_pt15_cut",           flav, isGluonSplit ,Track_pt[itrk]     , ww);
                            FillHisto2D_float_floatFromMap("track_nHit_vs_eta"     ,flav,isGluonSplit ,Track_nHitAll[itrk], Track_eta[itrk],ww);	
                            FillHisto2D_float_floatFromMap("track_nHit_vs_phi"     ,flav,isGluonSplit ,Track_nHitAll[itrk], Track_phi[itrk],ww);	
                            FillHisto2D_float_floatFromMap("track_HPix_vs_eta"     ,flav,isGluonSplit ,Track_nHitPixel[itrk], Track_eta[itrk],ww);	
                            FillHisto2D_float_floatFromMap("track_HPix_vs_phi"     ,flav,isGluonSplit ,Track_nHitPixel[itrk], Track_phi[itrk],ww);	
                            FillHisto2D_float_floatFromMap("track_eta_vs_phi"     ,flav,isGluonSplit ,Track_eta[itrk], Track_phi[itrk],ww);	
                        }	    
                        if (passNhit && passPix && passPt && passnormchi2 && passtrkdist && passtrklen){
                            FillHisto_floatFromMap("track_dz_cut",          flav, isGluonSplit ,Track_dz[itrk]      ,ww);
                        }
                        if (passNhit && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
                            FillHisto_intFromMap(  "track_HPix_cut",         flav, isGluonSplit ,Track_nHitPixel[itrk],ww);
                        }	
                        if (passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
                            FillHisto_intFromMap(  "track_nHit_cut",       flav, isGluonSplit ,Track_nHitAll[itrk],ww);	  
                            if (Track_nHitAll[itrk]<3) cout << "probl nHit " << Run << " " << Evt << " ijet " << ijet 
                                                            << " itrk " << itrk  << " nHit " << Track_nHitAll[itrk] <<
                                                           " pix " << Track_nHitPixel[itrk] <<  endl;
                        }
                        if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen ){
                            FillHisto_floatFromMap(  "track_IP2D_cut",         flav, isGluonSplit ,Track_IP2D[itrk],ww);	  
                        }	    
                    }
                    if (passNhit && passPix && passIPz && passPt && passnormchi2 && passtrkdist && passtrklen && passTrackIP2D){
                        ntracksel++;
                        if (dr_track<0.3) ntracksel03++;
	    
                        FillHisto_floatFromMap("track_dr",    flav, isGluonSplit ,dr_track	   ,ww);
                        FillHisto_floatFromMap("track_chi2",    flav, isGluonSplit ,Track_chi2[itrk]	   ,ww);
                        FillHisto_intFromMap(  "track_nHit",  flav, isGluonSplit ,Track_nHitAll[itrk],ww);
                        FillHisto_intFromMap(  "track_HPix",    flav, isGluonSplit ,Track_nHitPixel[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitStrip",    flav, isGluonSplit ,Track_nHitStrip[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitTIB",    flav, isGluonSplit ,Track_nHitTIB[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitTID",    flav, isGluonSplit ,Track_nHitTID[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitTOB",    flav, isGluonSplit ,Track_nHitTOB[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitTEC",    flav, isGluonSplit ,Track_nHitTEC[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitPXB",    flav, isGluonSplit ,Track_nHitPXB[itrk],ww);
                        FillHisto_intFromMap(  "track_nHitPXF",    flav, isGluonSplit ,Track_nHitPXF[itrk],ww);
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
                        FillHisto_floatFromMap("track_pt15",      flav, isGluonSplit ,Track_pt[itrk]      , ww);
                        FillHisto_floatFromMap("track_len_sel_zoom",  flav, isGluonSplit ,Track_length[itrk] , ww);
	    
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
                FillHisto_intFromMap(        "trk_multi_sel03",     flav, isGluonSplit ,ntracksel03	         , ww);  
                FillHisto2D_int_floatFromMap("seltrack_vs_jetpt", flav, isGluonSplit ,ptjet ,  ntracksel , ww);
                FillHisto_intFromMap(  "track_multi04", flav, isGluonSplit ,ntrkjet04  ,ww);
                FillHisto_intFromMap(  "track_multi", flav, isGluonSplit ,ntrkjet  ,ww);
	

            }//end produce jetProbaTree
            else {
                FillHisto_intFromMap(        "trk_multi_sel03",     flav, isGluonSplit ,Jet_nseltracks[ijet]	         , ww);  
            }

            //---------------------------------
            //fill information related to SV
            //---------------------------------
	
            n_sv           = Jet_SV_multi[ijet];	  
            FillHisto_intFromMap(  "sv_multi_0",      flav, isGluonSplit ,n_sv 	 ,         ww);

            if (n_sv>0){  
              FillHisto_intFromMap(  "sv_multi",     flav, isGluonSplit ,n_sv ,  ww);
              if (produceSVinfo) {
                chi2norm_sv    = SV_chi2[Jet_nFirstSV[ijet]]/SV_ndf[Jet_nFirstSV[ijet]];
                flightSig_sv   = SV_flight[Jet_nFirstSV[ijet]]/SV_flightErr[Jet_nFirstSV[ijet]];
                flight2DSig_sv = SV_flight2D[Jet_nFirstSV[ijet]]/SV_flight2DErr[Jet_nFirstSV[ijet]];
                mass_sv        =SV_mass[Jet_nFirstSV[ijet]];
                sv_dR_jet      =SV_deltaR_jet[Jet_nFirstSV[ijet]];
                sv_dR_dir_sum  =SV_deltaR_sum_dir[Jet_nFirstSV[ijet]];
                sv_dR_jet_sum  =SV_deltaR_sum_jet[Jet_nFirstSV[ijet]];
                sv_pt          =SV_vtx_pt[Jet_nFirstSV[ijet]];
                sveta          =SV_vtx_eta[Jet_nFirstSV[ijet]];
                svphi          =SV_vtx_phi[Jet_nFirstSV[ijet]];
	  
	  
                sv_flight3D    =SV_flight[Jet_nFirstSV[ijet]] ;  
                sv_flight3Derr =SV_flightErr[Jet_nFirstSV[ijet]];
                sv_flight2D    =SV_flight2D[Jet_nFirstSV[ijet]];    
                sv_flight2Derr =SV_flight2DErr[Jet_nFirstSV[ijet]];    
                sv_totchar     =SV_totCharge[Jet_nFirstSV[ijet]] ;
	     
                sv_nTrk        =SV_nTrk[Jet_nFirstSV[ijet]] ;  

                //-------------------------//
                //-----SV histograms-----//
                //-------------------------//   
                FillHisto_floatFromMap("sv_chi2norm",     flav, isGluonSplit ,chi2norm_sv        , ww);
                FillHisto_floatFromMap("sv_mass",         flav, isGluonSplit ,mass_sv,             ww);
                FillHisto_floatFromMap("sv_deltaR_jet",   flav, isGluonSplit ,sv_dR_jet,           ww);
                FillHisto_floatFromMap("sv_deltaR_sumJet",flav, isGluonSplit ,sv_dR_dir_sum,       ww);
                FillHisto_floatFromMap("sv_deltaR_sumDir",flav, isGluonSplit ,sv_dR_jet_sum,       ww);
                FillHisto_floatFromMap("sv_pt",           flav, isGluonSplit ,sv_pt,               ww);
                FillHisto_floatFromMap("sv_flight2D",     flav, isGluonSplit ,sv_flight2D,         ww);
                FillHisto_floatFromMap("sv_flight2Derr",  flav, isGluonSplit ,sv_flight2Derr,      ww);
                FillHisto_floatFromMap("sv_flightSig2D",  flav, isGluonSplit ,flight2DSig_sv,      ww);
                FillHisto_intFromMap(  "svnTrk",          flav, isGluonSplit ,sv_nTrk,             ww);
                FillHisto_floatFromMap("sv_eta",          flav, isGluonSplit ,sveta,               ww);
                FillHisto_floatFromMap("sv_phi",          flav, isGluonSplit ,svphi,               ww);	
                FillHisto_floatFromMap("sv_flight3D",     flav, isGluonSplit ,sv_flight3D,         ww);
                FillHisto_floatFromMap("sv_flight3Derr",  flav, isGluonSplit ,sv_flight3Derr,      ww);
                FillHisto_floatFromMap("sv_flight3DSig",  flav, isGluonSplit ,flightSig_sv,        ww);
	
                if (sv_nTrk >2) {
                    FillHisto_floatFromMap("sv_mass_3trk", flav, isGluonSplit ,mass_sv,ww);
                    FillHisto_floatFromMap("sv_flight2D_3trk",  flav, isGluonSplit ,sv_flight2D,        ww);
                    FillHisto_floatFromMap("sv_flight2DSig_3trk",  flav, isGluonSplit ,flight2DSig_sv,        ww);
                }

                FillHisto2D_float_floatFromMap("sv_mass_vs_flightDist3D"     ,flav,isGluonSplit ,sv_flight3D,mass_sv,ww);	
                FillHisto2D_float_floatFromMap("avg_sv_mass_vs_jetpt"        ,flav,isGluonSplit ,ptjet,mass_sv,ww);
                FillHisto2D_float_floatFromMap("sv_deltar_jet_vs_jetpt"      ,flav,isGluonSplit ,ptjet,sv_dR_jet,ww);
                FillHisto2D_float_floatFromMap("sv_deltar_sum_jet_vs_jetpt"  ,flav,isGluonSplit ,ptjet,sv_dR_dir_sum,ww);
                FillHisto2D_float_floatFromMap("sv_deltar_sum_dir_vs_jetpt"  ,flav,isGluonSplit ,ptjet,sv_dR_jet_sum,ww);	    
                }
	    
            }

            //-------------------------//
            //-----Tagger histograms-----//
            //-------------------------// 
            
            FillHisto_floatFromMap("TCHE",  flav, isGluonSplit, tche	  ,   ww);
            FillHisto_floatFromMap("TCHP",  flav, isGluonSplit, tchp	  ,   ww);
            FillHisto_floatFromMap("JP",    flav, isGluonSplit, jetproba  , ww);
            FillHisto_floatFromMap("JBP",   flav, isGluonSplit, jetbproba , ww);
            FillHisto_floatFromMap("SSV",   flav, isGluonSplit, ssvhe	  ,   ww);
            FillHisto_floatFromMap("SSVHP", flav, isGluonSplit, ssvhp	  ,   ww);
            FillHisto_floatFromMap("CSV",   flav, isGluonSplit, csv	  ,   ww);
            FillHisto_floatFromMap("CSVIVF",   flav, isGluonSplit, csvivf         ,   ww);
            FillHisto_floatFromMap("cMVA",   flav, isGluonSplit, csvmva         ,   ww);
            FillHisto_floatFromMap("cMVAv2",   flav, isGluonSplit, csvmvav2         ,   ww);
            FillHisto_floatFromMap("cTagCvsB",   flav, isGluonSplit, ctagCvsB         ,   ww);
            FillHisto_floatFromMap("cTagCvsBN",   flav, isGluonSplit, ctagCvsBN         ,   ww);
            FillHisto_floatFromMap("cTagCvsBP",   flav, isGluonSplit, ctagCvsBP         ,   ww);
            FillHisto_floatFromMap("cTagCvsL",   flav, isGluonSplit, ctagCvsL         ,   ww);
            FillHisto_floatFromMap("cTagCvsLN",   flav, isGluonSplit, ctagCvsLN         ,   ww);
            FillHisto_floatFromMap("cTagCvsLP",   flav, isGluonSplit, ctagCvsLP         ,   ww);
      
            FillHisto_floatFromMap("TCHE_extended1",  flav, isGluonSplit, tche  , ww);
            FillHisto_floatFromMap("TCHP_extended1",  flav, isGluonSplit, tchp  , ww);
      
            FillHisto_floatFromMap("TCHE_extended2",  flav, isGluonSplit, tche  , ww);
            FillHisto_floatFromMap("TCHP_extended2",  flav, isGluonSplit, tchp  , ww);
      
            FillHisto_floatFromMap("discri_ssche0",   flav, isGluonSplit, ssvhe , ww);
            FillHisto_floatFromMap("discri_sschp0",   flav, isGluonSplit, ssvhp , ww);      

            if (csv>0.244) FillHisto_floatFromMap("jet_pt_csvl",                 flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);
            if (csv>0.679) FillHisto_floatFromMap("jet_pt_csvm",                 flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);
            if (csvivf>0.423) FillHisto_floatFromMap("jet_pt_csvivfl",           flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);
            if (csvivf>0.814) FillHisto_floatFromMap("jet_pt_csvivfm",           flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);
            if (jetproba>0.275) FillHisto_floatFromMap("jet_pt_jpl",             flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);
            if (tchp>3.41 ) FillHisto_floatFromMap("jet_pt_tchpt",               flav, isGluonSplit ,Jet_genpt[ijet]    ,ww);

            float softmu   = Jet_SoftMu[ijet]  ;
            float solfel   = Jet_SoftEl[ijet];
            FillHisto_floatFromMap("SoftMu",      flav, isGluonSplit, softmu  ,   ww);
            FillHisto_floatFromMap("SoftEl",      flav, isGluonSplit, solfel  ,   ww);
      
            //---------------------------------
            //fill information related to muons
            //---------------------------------
         
            if (produceRecoMuon) {
             int nrecomu=0;
             int indrecomu=-1;
             if (!usePFMuon) {
                nrecomu=nselmuon;
                indrecomu=idxFirstMuon;
             }
             else {
               if (nMuon>0){
                for (int imu=0; imu< nMuon; imu++){
                    if (Muon_IdxJet[imu]==ijet ){
                        if (passMuonSelection(imu, ijet)){
                            if(nrecomu == 0) {
                                indrecomu = imu;
                            }
                            nrecomu++;
                        }
                    }
                }
               }
             }
             FillHisto_intFromMap(  "muon_multi",  flav, isGluonSplit , nrecomu   ,ww);
             if(indrecomu > -1){
                FillHisto_floatFromMap("muon_ptrel",    flav, isGluonSplit ,Muon_ptrel[indrecomu] ,ww);
                FillHisto_floatFromMap("muon_chi2",     flav, isGluonSplit ,Muon_chi2[indrecomu]  ,ww);
                FillHisto_floatFromMap("muon_Pt",     flav, isGluonSplit, Muon_pt[indrecomu] ,     ww);
                FillHisto_floatFromMap("muon_eta",    flav, isGluonSplit, Muon_eta[indrecomu] ,    ww);
                FillHisto_floatFromMap("muon_phi",    flav, isGluonSplit, Muon_phi[indrecomu] ,    ww);
                FillHisto_floatFromMap("muon_Ip3d",   flav, isGluonSplit, Muon_IP[indrecomu] ,     ww);     
                FillHisto_floatFromMap("muon_Ip2d",   flav, isGluonSplit, Muon_IP2D[indrecomu] ,   ww);     
                FillHisto_floatFromMap("muon_Sip3d",  flav, isGluonSplit, Muon_IPsig[indrecomu] ,  ww);     
                FillHisto_floatFromMap("muon_Sip2d",  flav, isGluonSplit, Muon_IP2Dsig[indrecomu] ,ww); 
	
                TLorentzVector themuon, thejet;
                thejet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], 0);
                themuon.SetPtEtaPhiM(Muon_pt[indrecomu], Muon_eta[indrecomu], Muon_phi[indrecomu], 0);
                FillHisto_floatFromMap("muon_DeltaR",         flav, isGluonSplit, themuon.DeltaR(thejet) , ww);

                FillHisto2D_float_floatFromMap("muon_ptrel_vs_jetpt", flav, isGluonSplit,ptjet,Muon_ptrel[indrecomu],ww);
                FillHisto2D_float_floatFromMap("muon_DeltaR_vs_jetpt",flav, isGluonSplit,ptjet,themuon.DeltaR(thejet),ww);
             }
           } // end of produceRecoMuon if
      

            // PFMuon
            int npfmu=0;
            int indpfmu=-1;
            if (usePFMuon) {
                npfmu=nselmuon;
                indpfmu= idxFirstMuon;
            }
            else {
             float minpf=5.;
             for (int im=0; im<nPFMuon; im++) {
                    if (PFMuon_IdxJet[im]==ijet && PFMuon_GoodQuality[im]>0) {
                        if (PFMuon_pt[im]> minpf) {
                            indpfmu=im;
                            minpf=PFMuon_pt[im];
                        }
                        npfmu++;
                    }
             }
            }  


            FillHisto_intFromMap(  "pfmuon_multi",  flav, isGluonSplit , npfmu   ,ww);
            if (indpfmu>-1) {
                    FillHisto_intFromMap(  "pfmuon_goodquality", flav, isGluonSplit, PFMuon_GoodQuality[indpfmu]        ,ww);
                    FillHisto_floatFromMap("pfmuon_pt",          flav, isGluonSplit, PFMuon_pt[indpfmu],                 ww);
                    FillHisto_floatFromMap("pfmuon_eta",   	 flav, isGluonSplit, PFMuon_eta[indpfmu],                ww);
                    FillHisto_floatFromMap("pfmuon_phi",         flav, isGluonSplit, PFMuon_phi[indpfmu],                ww);
                    FillHisto_floatFromMap("pfmuon_ip",	         flav, isGluonSplit, PFMuon_IP[indpfmu],              ww);
                    FillHisto_floatFromMap("pfmuon_ptrel",       flav, isGluonSplit, PFMuon_ptrel[indpfmu],              ww);
                    FillHisto_floatFromMap("pfmuon_ratio",       flav, isGluonSplit, PFMuon_ratio[indpfmu],              ww);
                    FillHisto_floatFromMap("pfmuon_ratiorel",    flav, isGluonSplit, PFMuon_ratioRel[indpfmu],           ww);
                    FillHisto_floatFromMap("pfmuon_deltar",	 flav, isGluonSplit, PFMuon_deltaR[indpfmu],             ww);
            }
      
            //--------------------------------------
            //fill information related to electrons
            //---------------------------------------

            //PFElectron
            int npfel=0;
            int indpfel=-1;
            float minpfe=0;


            for (int im=0; im<nPFElectron; im++) {
                    if (PFElectron_IdxJet[im]==ijet && PFElectron_pt[im]>2.) {
                        if (PFElectron_pt[im]> minpfe) {  
                            indpfel=im;
                            minpfe=PFElectron_pt[im];
                        }
                        npfel++;
                    }
            }


            FillHisto_intFromMap(  "pfelectron_multi",  flav, isGluonSplit , npfel   ,ww);
            if (indpfel>-1) {
                    FillHisto_floatFromMap("pfelectron_pt",          flav, isGluonSplit, PFElectron_pt[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_eta",   	 flav, isGluonSplit, PFElectron_eta[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_phi",         flav, isGluonSplit, PFElectron_phi[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ip",	 flav, isGluonSplit, PFElectron_IP[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ptrel",       flav, isGluonSplit, PFElectron_ptrel[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ratio",       flav, isGluonSplit, PFElectron_ratio[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ratiorel",    flav, isGluonSplit, PFElectron_ratioRel[indpfel],             ww);
                    FillHisto_floatFromMap("pfelectron_deltar",	 flav, isGluonSplit, PFElectron_deltaR[indpfel],                ww);

                    if (PFElectron_pt[indpfel]<50) {
		    FillHisto_floatFromMap("pfelectron_eta_1",   	 flav, isGluonSplit, PFElectron_eta[indpfel],                ww);
		    FillHisto_floatFromMap("pfelectron_jeta_1",   	 flav, isGluonSplit, etajet,                ww);
                    FillHisto_floatFromMap("pfelectron_ip_1",	 flav, isGluonSplit, PFElectron_IP[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ip2d_1",	 flav, isGluonSplit, PFElectron_IP2D[indpfel],                ww);
		    FillHisto_floatFromMap("pfelectron_ptratio_1",   	 flav, isGluonSplit, PFElectron_pt[indpfel]/ptjet,                ww);
                    FillHisto_floatFromMap("pfelectron_ratio_1",       flav, isGluonSplit, PFElectron_ratio[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ratiorel_1",    flav, isGluonSplit, PFElectron_ratioRel[indpfel],             ww);
                    FillHisto_floatFromMap("pfelectron_deltar_1",	 flav, isGluonSplit, PFElectron_deltaR[indpfel],                ww);
                    }
                    else {
		    FillHisto_floatFromMap("pfelectron_eta_2",   	 flav, isGluonSplit, PFElectron_eta[indpfel],                ww);
		    FillHisto_floatFromMap("pfelectron_jeta_2",   	 flav, isGluonSplit, etajet,                ww);
                    FillHisto_floatFromMap("pfelectron_ip_2",	 flav, isGluonSplit, PFElectron_IP[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ip2d_2",	 flav, isGluonSplit, PFElectron_IP2D[indpfel],                ww);
		    FillHisto_floatFromMap("pfelectron_ptratio_2",   	 flav, isGluonSplit, PFElectron_pt[indpfel]/ptjet,                ww);
                    FillHisto_floatFromMap("pfelectron_ratio_2",       flav, isGluonSplit, PFElectron_ratio[indpfel],                ww);
                    FillHisto_floatFromMap("pfelectron_ratiorel_2",    flav, isGluonSplit, PFElectron_ratioRel[indpfel],             ww);
                    FillHisto_floatFromMap("pfelectron_deltar_2",	 flav, isGluonSplit, PFElectron_deltaR[indpfel],                ww);

                    }
            }

            //------------------------------------
            //fill information related to TagVar 
            //------------------------------------

            if (produceTagVarTree) {

                if (Jet_nFirstSVTagVar[ijet]>-1) {
                // look at first vertex info only (possible to access all SV info for vertexNTracks & vertexmass
                FillHisto_floatFromMap("tagvar_vertexNTracks",      flav, isGluonSplit, TagVar_vertexNTracks[Jet_nFirstSVTagVar[ijet]],   ww);
                FillHisto_floatFromMap("tagvar_vertexmass",         flav, isGluonSplit, TagVar_vertexMass[Jet_nFirstSVTagVar[ijet]],   ww);
                if (TagVar_vertexNTracks[ijet]>=3) FillHisto_floatFromMap("tagvar_vertexmass3trk",     flav, isGluonSplit, TagVar_vertexMass[Jet_nFirstSVTagVar[ijet]],   ww);
                FillHisto_floatFromMap("tagvar_vertexJetDeltaR",    flav, isGluonSplit, TagVar_vertexJetDeltaR[Jet_nFirstSVTagVar[ijet]],   ww);
                }

            }

            //--------------------------------------
            //fill information related to TagVarCSV
            //--------------------------------------

            if (produceTagVarCSVTree) { 
                FillHisto_floatFromMap("tagvarCSV_vertexCategory", flav, isGluonSplit, TagVarCSV_vertexCategory[ijet],   ww);
                FillHisto_floatFromMap("tagvarCSV_Sig2dAboveCharm", flav, isGluonSplit, TagVarCSV_trackSip2dSigAboveCharm[ijet],   ww);
                for (int inrel=Jet_nFirstTrkEtaRelTagVarCSV[ijet]; inrel<Jet_nLastTrkEtaRelTagVarCSV[ijet]; inrel++) {
                   FillHisto_floatFromMap("tagvarCSV_trackEtaRel", flav, isGluonSplit, TagVarCSV_trackEtaRel[inrel],ww);
                }
                FillHisto_floatFromMap("tagvarCSV_vertexmass", flav, isGluonSplit, TagVarCSV_vertexMass[ijet],ww);
                if (TagVarCSV_vertexNTracks[ijet]>=3) FillHisto_floatFromMap("tagvarCSV_vertexmass3trk",     flav, isGluonSplit, TagVarCSV_vertexMass[ijet],   ww);
                FillHisto_floatFromMap("tagvarCSV_vertexNTracks", flav, isGluonSplit, TagVarCSV_vertexNTracks[ijet],ww);
                FillHisto_floatFromMap("tagvarCSV_energyratio", flav, isGluonSplit, TagVarCSV_vertexEnergyRatio[ijet],   ww);
                for (int inrel=Jet_nFirstTrkTagVarCSV[ijet]; inrel<Jet_nLastTrkTagVarCSV[ijet]; inrel++) {
                   FillHisto_floatFromMap("tagvarCSV_trackSip3dSig", flav, isGluonSplit, TagVarCSV_trackSip3dSig[inrel],ww);
                }
                FillHisto_floatFromMap("tagvarCSV_2DsigFlightDist", flav, isGluonSplit, TagVarCSV_flightDistance2dSig[ijet],   ww);
                FillHisto_floatFromMap("tagvarCSV_vertexJetDeltaR",    flav, isGluonSplit, TagVarCSV_vertexJetDeltaR[ijet],   ww);
                if (TagVarCSV_vertexCategory[ijet]==0) {
                  FillHisto_floatFromMap("tagvarCSV_vertexmass_cat0", flav, isGluonSplit, TagVarCSV_vertexMass[ijet],ww);
                  if (TagVarCSV_vertexNTracks[ijet]>=3) FillHisto_floatFromMap("tagvarCSV_vertexmass3trk_cat0",     flav, isGluonSplit, TagVarCSV_vertexMass[ijet],   ww);
                  FillHisto_floatFromMap("tagvarCSV_vertexNTracks_cat0", flav, isGluonSplit, TagVarCSV_vertexNTracks[ijet],ww);
                  FillHisto_floatFromMap("tagvarCSV_2DsigFlightDist_cat0", flav, isGluonSplit, TagVarCSV_flightDistance2dSig[ijet],   ww);
                  FillHisto_floatFromMap("tagvarCSV_vertexJetDeltaR_cat0",    flav, isGluonSplit, TagVarCSV_vertexJetDeltaR[ijet],   ww);
                }
            }
 

        }
        //----------------------------------
        //End Loop on jets 
        //-----------------------------------
    
        //---------------------------------
        //fill jet multiplicity
        //---------------------------------
        if(isData)
        {
            nJet_data->Fill(njet_data);
            nJet_csvLoose_data->Fill(njet_csvLoose_data);
            nJet_csvMedium_data->Fill(njet_csvMedium_data);
            nJet_csvTight_data->Fill(njet_csvTight_data);
        }
        else {
            nJet_mc->Fill(njet_mc, ww);
            nJet_csvLoose_mc->Fill(njet_csvLoose_mc);
            nJet_csvMedium_mc->Fill(njet_csvMedium_mc);
            nJet_csvTight_mc->Fill(njet_csvTight_mc);
            if (finselect) {
              pt_hat_fin->Fill(pthat,ww);
              PVz_mc->Fill(PVz,ww);
              dPVz_mc->Fill(abs(PVz-GenPVz),ww);
            }
        }
    
    
    
        //cout << "---------------------End event---------------------------" << endl;
    
    }//----------------------------------------End loop events ----------------------------------------------------------------------------------------------//  
  
  
    cout << "----------------Configuration-------------------" <<endl;
    cout <<" Run over "<< N_event_data_before_sel<<" data events and " << N_event_mc_before_sel<< " mc events."<<endl;
    cout <<" After selection -> "<< N_event_data_after_sel<<" data events and "<< N_event_mc_after_sel<< "mc events left" <<endl;   
    cout << endl;
  
    myfile->cd();
    nPU_mc->Write();
    nPV_data->Write();
    nPV_mc->Write();
    nPV_mc_unw->Write();
    pt_hat->Write();
    pt_hat_fin->Write();
    jet_pt_mc->Write();
    nJet_data->Write();
    nJet_mc->Write();
    nJet_mc_inc->Write();
    PVz_mc->Write();
    dPVz_mc->Write();
    nJet_csvLoose_data->Write();
    nJet_csvMedium_data->Write();
    nJet_csvTight_data->Write();
    nJet_csvLoose_mc->Write();
    nJet_csvMedium_mc->Write();
    nJet_csvTight_mc->Write();
  

    for (unsigned int i=0; i<HistoBtag.size(); i++) {
        HistoBtag[i]->Write();
    }
    for (unsigned int i=0; i<HistoBtag2D.size(); i++) {
        HistoBtag2D[i]->Write();
    }  
    myfile->Close();


}

void CommPlotProducer::AddHisto(TString name, TString title, int nbins, float min, float max)  {
  
    TH1D* h_b      = new TH1D(name+"_b",title+"_b",nbins,min,max);
    TH1D* h_bfromg = new TH1D(name+"_bfromg",title+"_bfromg",nbins,min,max);  
    TH1D* h_c      = new TH1D(name+"_c",title+"_c",nbins,min,max);  
    TH1D* h_l      = new TH1D(name+"_l",title+"_l",nbins,min,max);
    TH1D* h_g      = new TH1D(name+"_g",title+"_g",nbins,min,max);
    TH1D* h_data   = new TH1D(name+"_data",title+"_data",nbins,min,max);
    TH1D* h_cfromg = new TH1D(name+"_cfromg",title+"_cfromg",nbins,min,max);  
    TH1D* h_puu = new TH1D(name+"_puu",title+"_puu",nbins,min,max);  
  
  
    h_b        ->Sumw2();
    h_bfromg   ->Sumw2();  
    h_c        ->Sumw2();  
    h_cfromg   ->Sumw2();  
    h_l        ->Sumw2(); 
    h_data     ->Sumw2();
    h_puu     ->Sumw2();
    h_g     ->Sumw2();
  
    HistoBtag.push_back(h_b);
    HistoBtag.push_back(h_bfromg);  
    HistoBtag.push_back(h_c);  
    HistoBtag.push_back(h_l);  
    HistoBtag.push_back(h_data);  
    HistoBtag.push_back(h_cfromg);  
    HistoBtag.push_back(h_g);  
    HistoBtag.push_back(h_puu);  
    HistoBtag_map[name.Data()] = numb_histo;
  
    numb_histo++;
  
}


void CommPlotProducer::FillHisto_float(int flavour, bool isGS, int number, float value, double weight)  {
  
    if (!isData){
        if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*8 +0]->Fill(value,weight);
        else if (fabs(flavour)==5 && isGS)              HistoBtag[number*8 +1]->Fill(value,weight);
        else if (fabs(flavour)==4 && !isGS)             HistoBtag[number*8 +2]->Fill(value,weight); 
        else if (fabs(flavour)==4 && isGS)              HistoBtag[number*8 +5]->Fill(value,weight); 
        else if (fabs(flavour)< 4)                      HistoBtag[number*8 +3]->Fill(value,weight);
        else if (fabs(flavour)==21)                     HistoBtag[number*8 +6]->Fill(value,weight);
        else                                            HistoBtag[number*8 +7]->Fill(value,weight); 
    
    }  
    else                                              HistoBtag[number*8 +4]->Fill(value);
  
  
}
void CommPlotProducer::FillHisto_int(int flavour, bool isGS, int number, int value, double weight)  {
  
    if (!isData){
        if (fabs(flavour)==5 && !isGS)             HistoBtag[number*8 +0]->Fill(value,weight);
        else if (fabs(flavour)==5 && isGS)              HistoBtag[number*8 +1]->Fill(value,weight);
        else if (fabs(flavour)==4 && !isGS)             HistoBtag[number*8 +2]->Fill(value,weight); 
        else if (fabs(flavour)==4 && isGS)              HistoBtag[number*8 +5]->Fill(value,weight); 
        else if (fabs(flavour)< 4)                      HistoBtag[number*8 +3]->Fill(value,weight);
        else if (fabs(flavour)==21)                     HistoBtag[number*8 +6]->Fill(value,weight);
        else                                            HistoBtag[number*8 +7]->Fill(value,weight);

    }  
    else                                         HistoBtag[number*8 +4]->Fill(value);
  
}


void CommPlotProducer::FillHisto_floatFromMap(TString name, int flavour, bool isGS, float value, double weight)  {
  
  
    int number = HistoBtag_map[name.Data()] ;
    if (!isData){
        if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*8 +0]->Fill(value,weight);
        else if (fabs(flavour)==5 && isGS)              HistoBtag[number*8 +1]->Fill(value,weight);
        else if (fabs(flavour)==4 && !isGS)             HistoBtag[number*8 +2]->Fill(value,weight); 
        else if (fabs(flavour)==4 && isGS)              HistoBtag[number*8 +5]->Fill(value,weight); 
        else if (fabs(flavour)< 4)                      HistoBtag[number*8 +3]->Fill(value,weight);
        else if (fabs(flavour)==21)                     HistoBtag[number*8 +6]->Fill(value,weight);
        else                                            HistoBtag[number*8 +7]->Fill(value,weight);

    
    }  
    else                                              HistoBtag[number*8 +4]->Fill(value);
  
   
}


void CommPlotProducer::FillHisto_intFromMap(TString name, int flavour, bool isGS, int value, double weight)  {
  
    int number = HistoBtag_map[name.Data()] ;
    if (!isData){
        if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*8 +0]->Fill(value,weight);
        else if (fabs(flavour)==5 && isGS)              HistoBtag[number*8 +1]->Fill(value,weight);
        else if (fabs(flavour)==4 && !isGS)             HistoBtag[number*8 +2]->Fill(value,weight); 
        else if (fabs(flavour)==4 && isGS)              HistoBtag[number*8 +5]->Fill(value,weight); 
        else if (fabs(flavour)< 4)                      HistoBtag[number*8 +3]->Fill(value,weight);
        else if (fabs(flavour)==21)                     HistoBtag[number*8 +6]->Fill(value,weight);
        else                                            HistoBtag[number*8 +7]->Fill(value,weight);

    }  
    else                                         HistoBtag[number*8 +4]->Fill(value);
  
}

//-----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------2D PLOTS-----------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//


void CommPlotProducer::AddHisto2D(TString name, TString title, int nbins, float min, float max, int nbins2, float
                                  min2, float max2)  {
  
    TH2D* h_b      = new TH2D(name+"_b",title+"_b",nbins,min,max,nbins2,min2,max2);
    TH2D* h_bfromg = new TH2D(name+"_bfromg",title+"_bfromg",nbins,min,max,nbins2,min2,max2);  
    TH2D* h_c      = new TH2D(name+"_c",title+"_c",nbins,min,max,nbins2,min2,max2);  
    TH2D* h_cfromg = new TH2D(name+"_cfromg",title+"_cfromg",nbins,min,max,nbins2,min2,max2);  
    TH2D* h_l      = new TH2D(name+"_l",title+"_l",nbins,min,max,nbins2,min2,max2);
    TH2D* h_data   = new TH2D(name+"_data",title+"_data",nbins,min,max,nbins2,min2,max2);
    TH2D* h_puu = new TH2D(name+"_puu",title+"_puu",nbins,min,max,nbins2,min2,max2);  
  
  
    h_b        ->Sumw2();
    h_bfromg   ->Sumw2();  
    h_c        ->Sumw2();  
    h_cfromg   ->Sumw2();  
    h_l        ->Sumw2(); 
    h_data     ->Sumw2();
    h_puu     ->Sumw2();
  
    HistoBtag2D.push_back(h_b);
    HistoBtag2D.push_back(h_bfromg);  
    HistoBtag2D.push_back(h_c);  
    HistoBtag2D.push_back(h_l);  
    HistoBtag2D.push_back(h_data);  
    HistoBtag2D.push_back(h_cfromg);  
    HistoBtag2D.push_back(h_puu);  
    HistoBtag2D_map[name.Data()] = numb_histo2D;
    numb_histo2D++;
  
}



void CommPlotProducer::FillHisto2D_int_floatFromMap(TString name, int flavour, bool isGS, int value, float value2, double weight)  {
  
  
    int number = HistoBtag2D_map[name.Data()] ;
    if (!isData){
        if (fabs(flavour)==5 && !isGS)                  HistoBtag2D[number*7 +0]->Fill(value,value2,weight);
        else if (fabs(flavour)==5 && isGS)              HistoBtag2D[number*7 +1]->Fill(value,value2,weight);
        else if (fabs(flavour)==4 && !isGS)             HistoBtag2D[number*7 +2]->Fill(value,value2,weight); 
        else if (fabs(flavour)==4 && isGS)              HistoBtag2D[number*7 +5]->Fill(value,value2,weight); 
        else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*7 +3]->Fill(value,value2,weight);
        else                                            HistoBtag2D[number*7 +6]->Fill(value,value2,weight); 
    
    }  
    else                                              HistoBtag2D[number*7 +4]->Fill(value,value2);
  
   
}

void CommPlotProducer::FillHisto2D_float_floatFromMap(TString name, int flavour, bool isGS, float value, float value2, double weight)  {
  
  
    int number = HistoBtag2D_map[name.Data()] ;
    if (!isData){
        if (fabs(flavour)==5 && !isGS)                  HistoBtag2D[number*7 +0]->Fill(value,value2,weight);
        else if (fabs(flavour)==5 && isGS)              HistoBtag2D[number*7 +1]->Fill(value,value2,weight);
        else if (fabs(flavour)==4 && !isGS)              HistoBtag2D[number*7 +2]->Fill(value,value2,weight); 
        else if (fabs(flavour)==4 && isGS)               HistoBtag2D[number*7 +5]->Fill(value,value2,weight); 
        else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*7 +3]->Fill(value,value2,weight);
        else                                            HistoBtag2D[number*7 +6]->Fill(value,value2,weight); 
    
    }  
    else                                              HistoBtag2D[number*7 +4]->Fill(value,value2);
  
   
}


//-----------------------------------------------------------------------------------------------------------------//

bool CommPlotProducer::passMuonSelection(int muidx, int ijet){
  
  
    TLorentzVector muon, jet;
  
  
    jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], 0);
    muon.SetPtEtaPhiM(Muon_pt[muidx], Muon_eta[muidx], Muon_phi[muidx], 0);
  
    bool cut_mu_pass=false;
    if (Muon_pt[muidx]> 5	&& TMath::Abs(Muon_eta[muidx]) < 2.4 && Muon_isGlobal[muidx] == 1     &&
        Muon_nMuHit[muidx]> 0 && Muon_nMatched[muidx]>1	 && Muon_nTkHit[muidx]>10 &&
        Muon_nPixHit[muidx]> 1	&& Muon_nOutHit[muidx]<3  && Muon_chi2Tk[muidx]< 10    &&
        Muon_chi2[muidx]< 10	//&& Muon_vz[muidx]< 2  
        && 	jet.DeltaR(muon) < 0.4
        && TMath::Abs(Muon_vz[muidx]-PV_z[0]) < 2.) 
        cut_mu_pass=true;
  
  
    return cut_mu_pass;
  
}

bool CommPlotProducer::passPFMuonSelection(int muidx){
  
    bool cut_mu_pass=false;
    if (PFMuon_pt[muidx]> 5	&& TMath::Abs(PFMuon_eta[muidx]) < 2.4 && PFMuon_GoodQuality[muidx] >0    
        && PFMuon_deltaR[muidx] < 0.4 )
        cut_mu_pass=true;
  
    return cut_mu_pass;
  
}
//--------------------------------------------------------------------------------------------------------------------//
bool CommPlotProducer::passTrigger(TString trigger, int pttrig) {


    // FOR 2012 Trigger ! Not valid for 2011...

    bool passTrig=false;
    //bool Jet30  = false, Jet60  = false, Jet150 = false, Jet190 = false, Jet240 = false;
    bool Jet40  = false, Jet60=false,  Jet80  = false, Jet140 = false;
    bool Jet200 = false, Jet260 = false, Jet320 = false;
    bool Jet400 = false, Jet450 = false, Jet500 = false;
    bool Jet20  = false, Jet70  = false, Jet110 = false, Jet170 = false; 
    int triggerIdx = 0, bitIdx = 0;

    if ( trigger=="jet") {
        triggerIdx = 0;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;

        triggerIdx = 1;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet60  = true;
   
        triggerIdx = 2;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet80  = true;
   
        triggerIdx = 3;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet140 = true;
   
        triggerIdx = 4;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet200 = true;
   
        triggerIdx = 5;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet260 = true;
   
        triggerIdx = 6;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet320 = true;

        triggerIdx = 7;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet400 = true;

        triggerIdx = 8;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet450 = true;

        triggerIdx = 9;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet500 = true;


        if ( pttrig ==  40 && Jet40 )  passTrig=true;
        if ( pttrig ==  60 && Jet60 )  passTrig=true;
        if ( pttrig ==  80 && Jet80 )  passTrig=true;
        if ( pttrig == 140 && Jet140 ) passTrig=true;
        if ( pttrig == 200 && Jet200 ) passTrig=true;
        if ( pttrig == 260 && Jet260 ) passTrig=true;
        if ( pttrig == 320 && Jet320 ) passTrig=true;
        if ( pttrig == 400 && Jet400 ) passTrig=true;
        if ( pttrig == 450 && Jet450 ) passTrig=true;
        if ( pttrig == 500 && Jet500 ) passTrig=true;

        if (!passTrig) { return false; }

        //-----------------------------------
        //Determine if there is at least 
        //one jet which pass the trigger
        //in the event => away from the TO
        //-----------------------------------

        bool JetPtCut = false;
        for (int ijet=0; ijet<nJet ; ijet++) {
            float ptjet = Jet_pt[ijet];
            float etajet = fabs(Jet_eta[ijet]);
//            if (      pttrig ==  40 && ptjet >  60. && etajet < 2.4 ) JetPtCut = true;
            if (      pttrig ==  40 && ptjet >  50. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig ==  60 && ptjet >  70. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig ==  80 && ptjet > 100. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 140 && ptjet > 160. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 200 && ptjet > 220. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 260 && ptjet > 300. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 320 && ptjet > 360. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 400 && ptjet > 425. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 450 && ptjet > 475. && etajet < 2.4 ) JetPtCut = true;
            else if ( pttrig == 500 && ptjet > 525. && etajet < 2.4 ) JetPtCut = true;
        }
        if (passTrig && JetPtCut) {return true;}
        else {return false;}
    }
   
    else if ( trigger=="btag" ) {
        triggerIdx = 30;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet20  = true;
   
        triggerIdx = 31;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;
   
        triggerIdx = 32;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet70  = true;
   
        triggerIdx = 33;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet110 = true;
   
        triggerIdx = 34;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet170 = true;

        if ( pttrig ==  20 && Jet20 )  passTrig=true;
        if ( pttrig ==  40 && Jet40 )  passTrig=true;
        if ( pttrig ==  70 && Jet70 )  passTrig=true;
        if ( pttrig == 110 && Jet110 ) passTrig=true;
        if ( pttrig == 170 && Jet170 ) passTrig=true;
    
    
        if (!passTrig) { return false; }
    
        //-----------------------------------
        //Determine if there is at least 
        //two jets which pass the trigger
        //in the event => away from the TO
        //-----------------------------------

        int njtrig=0;
        if (pttrig ==300) njtrig+=1;
        for (int ijet = 0; ijet < nJet; ijet++) {
            float ptjet = Jet_pt[ijet];
            float etajet = fabs(Jet_eta[ijet]);
//            if ( pttrig ==20  &&  ptjet > 60. && etajet < 2.4 )  njtrig++;
            if ( pttrig == 20  &&  ptjet > 30. && etajet < 2.4 )  njtrig++;
            if ( pttrig == 40  &&  ptjet > 50. && etajet < 2.4 )  njtrig++;
            if ( pttrig == 70  &&  ptjet > 80. && etajet < 2.4 )  njtrig++;
            if ( pttrig == 110  &&  ptjet > 120. && etajet < 2.4 ) njtrig++;
            if ( pttrig == 170  &&  ptjet > 190. && etajet < 2.4 ) njtrig++;
        }
        if (passTrig && njtrig>1) {return true;}
        else {return false;}

    }
    return false;
}

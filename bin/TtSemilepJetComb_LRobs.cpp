#include <iostream>
#include <cassert>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TString.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TKey.h>
#include <vector>
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
// #include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"
#include "TopQuarkAnalysis/TopTools/interface/LRHelpFunctions.h"

using namespace std;



///////////////////////
// Constants         //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//input files
const TString dir = "/msa1/jmmaes/CMSSW_1_6_9_NOPAT/src/TopQuarkAnalysis/Examples/crab/SanityCheck/rootfiles/";
const  TString  inputFile    = "plots_obs_all_semimu.root";
const  TString  outputFile   = "final_obs_plots_all_semimu.root";
const  TString  outputPSfile = "LRJetCombAllObs_all_semimu.ps";

//const  TString  inputFile    = "plots_obs_tt0j_semimu.root";
//const  TString  outputFile   = "final_obs_plots_tt0j_semimu.root";
//const  TString  outputPSfile = "LRJetCombAllObs_tt0j_semimu.ps";

//observable histogram variables
const  int      nrJetCombObs  		= 66;
//const  int      nrJetCombObs  		= 1;
//const  int      JetCombObs[nrJetCombObs] 	= {62};
const  int      JetCombObs[nrJetCombObs] 	= {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66};

//definitions to make life more easy
TFormula gauss("gauss", "gaus");
TFormula symgauss("symgauss", "[0]*(exp(-0.5*(x/[1])**2))");
TFormula dblgauss("dblgauss", "[0]*(exp(-0.5*((x-[1])/[2])**2)+exp(-0.5*((x+[3])/[4])**2))");
TFormula symdblgauss("symdblgauss", "[0]*(exp(-0.5*((x-[1])/[2])**2)+exp(-0.5*((x+[1])/[2])**2))");
TFormula sigm("sigm", "[0]/(1 + 1/exp([1]*([2] - x)))");
TFormula sigmc("sigmc", "[0]/(1 + 1/exp([1]*([2] - x)))+[3]");
TFormula dblsigm("dblsigm", "[0]/(1 + 1/exp([1]**2*([2] - x)))/(1 + 1/exp([3]**2*(x - [4])))");
TFormula symdblsigm("symdblsigm", "[0]/(1 + 1/exp([1]**2*([2] - x)))/(1 + 1/exp([1]**2*([2] + x)))");

//observable fit functions
/*const char*     JetCombObsFits[nrJetCombObs]= {           
   "pol4" //obs15
   };*/
const char*     JetCombObsFits[nrJetCombObs] 	= { 
  "[0]/(1 + 1/exp([1]*([2] - x)))",  //obs1	
  "[0]/(1 + 1/exp([1]*([2] - x)))",  //obs2	
  "gaus",  //obs3
  "gaus", //obs4
  "gaus", //obs5
  "([0]+[3]*abs(x)/x)*(1-exp([1]*(abs(x)-[2])))",  //obs6	
  "[0]/(1 + 1/exp([1]*([2] - x)))",  //obs7
  "sigmc", //obs8
  "gauss",//obs9
  "pol4",//obs10
  "sigmc",//obs11
  "gauss",//obs12
  "pol4",//obs13
  "sigmc",//obs14
  "gauss",//obs15
  "pol4",//obs16
  "sigm+sigmc",//obs17
  "pol4",//obs18
  "pol8",//obs19
  "sigmc",//obs20 
  "gauss",//obs21
  "pol4",//obs22
  "sigmc",//obs23
  "gauss",//obs24
  "gauss",//obs25
  "sigmc",//obs26
  "gauss",//obs27
  "gauss",//obs28
  "sigmc",//obs29
  "gauss",//obs30
  "pol4",//obs31
  "sigmc",//obs32
  "pol2",//obs33
  "pol2",//obs34
  "sigmc",//obs35
  "pol0",//obs36
  "pol0",//obs37
  "pol2",//obs38
  "gauss",//obs39
  "gauss",//obs40 
  //"",//obs40
  "pol5",//obs41
  "pol3",//obs42
  "pol3",//obs43
  "pol3",//obs44
  //"",//obs44
  "pol5",//obs45
  "gaus+gaus(3)",//obs46
  "gaus",//obs47
  "sigmc",//obs48
  //"",//obs48
  "pol5",//obs49
  "pol3",//obs50
  "gauss",//obs51
  "gauss",//obs52
  //"",//obs52
  "pol5",//obs53
  "gaus+gaus(3)",//obs54
  "sigm",//obs55
  "sigmc",//obs56
  //"",//obs56
  "pol5",//obs57
  "pol3",//obs58
  "pol2",//obs59
  "pol2",//obs60
  //"",//obs60
  "pol5",//obs61
  "pol7",//obs62
  "gauss+gauss+gauss",//obs63
  "pol5",//obs64
  "gauss+gauss",//obs65
  "pol5"//obs66
  };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//
// Global variables
//
LRHelpFunctions *myLRhelper;
vector<int> obsNrs;
vector<double> obsMin,obsMax;
vector<const char*> obsFits;

//
// Main analysis
//

int main() { 
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  vector<string> xlabels;

  xlabels.push_back("min (|#theta_{b1}-#pi/2|, (|#theta_{b2}-#pi/2|)"); //obs1
  xlabels.push_back("max (|#theta_{b1}-#pi/2|, (|#theta_{b2}-#pi/2|)"); //obs2
  xlabels.push_back("min(p_{T,b1}, p_{T,b2}) [GeV/c]"); //obs3
  xlabels.push_back("max(p_{T,b1}, p_{T,b2}) [GeV/c]"); //obs4
  xlabels.push_back("min(p_{T,l1}, p_{T,l2}) [GeV/c]"); //obs5
  xlabels.push_back("max(p_{T,l1}, p_{T,l2}) [GeV/c]"); //obs6
  xlabels.push_back("|#Delta#phi(b_{1}, b_{2})|"); //obs7
  xlabels.push_back("|#Delta#theta(b_{1}, b_{2})|"); //obs8
  xlabels.push_back("min(|#Delta#phi(b_{1}, l_{1})|,|#Delta#phi(b_{2}, l_{2})|)"); //obs9
  xlabels.push_back("max(|#Delta#phi(b_{1}, l_{1})|,|#Delta#phi(b_{2}, l_{2})|)"); //obs10
  xlabels.push_back("min(|#Delta#theta(b_{1}, l_{1})|,|#Delta#theta(b_{2}, l_{2})|)"); //obs11
  xlabels.push_back("max(|#Delta#theta(b_{1}, l_{1})|,|#Delta#theta(b_{2}, l_{2})|)"); //obs12
  xlabels.push_back("m_{l_{1},l_{2}} [GeV/c^{2}]"); //obs13
  xlabels.push_back("E_{T,jet3} / E_{T,jet1}"); //obs14
  xlabels.push_back("E_{T,jet3} / E_{T,jet2}"); //obs15

  // define all histograms & fit functions
  //to replace with something more elegant
  for(int j = 0; j < nrJetCombObs; j++){
    obsNrs.push_back(JetCombObs[j]);
    obsFits.push_back(JetCombObsFits[j]);
  }
  myLRhelper = new LRHelpFunctions();
  //load plots
  myLRhelper->readObsHistsAndFits(dir+inputFile, obsNrs, false);
  cout << "Histos loaded\n";
  myLRhelper->recreateFitFct(obsNrs, obsFits);
  cout << "fit functions loaded\n";

  // manually set some initial values for fit function parameters
  vector<double> parsFobs1; parsFobs1.push_back(0.65); parsFobs1.push_back(-0.1); parsFobs1.push_back(30);
  myLRhelper -> setObsFitParameters(3,parsFobs1);
    myLRhelper -> setObsFitParameters(4,parsFobs1);

  vector<double> parsFobs2; parsFobs2.push_back(0.9); parsFobs2.push_back(-0.1); parsFobs2.push_back(15);
  myLRhelper -> setObsFitParameters(6,parsFobs2);
  
  vector<double> parsFobs8; 
  parsFobs8.push_back(.5);
  parsFobs8.push_back(-.03);
  parsFobs8.push_back(200);
  parsFobs8.push_back(.4);
  myLRhelper -> setObsFitParameters(8,parsFobs8);
  
  vector<double> parsFobs9; 
  parsFobs9.push_back(.3);
  parsFobs9.push_back(-.5);
  parsFobs9.push_back(4);
  myLRhelper -> setObsFitParameters(9,parsFobs9);  

  vector<double> parsFobs11; 
  parsFobs11.push_back(.3);
  parsFobs11.push_back(-.04);
  parsFobs11.push_back(100);
  parsFobs11.push_back(.4);
  myLRhelper -> setObsFitParameters(11,parsFobs11);

  vector<double> parsFobs12; 
  parsFobs12.push_back(.4);
  parsFobs12.push_back(.1);
  parsFobs12.push_back(4);
  myLRhelper -> setObsFitParameters(12,parsFobs12); 

  vector<double> parsFobs14; 
  parsFobs14.push_back(-10000000);
  parsFobs14.push_back(-.2);
  parsFobs14.push_back(-100);
  parsFobs14.push_back(54);
  myLRhelper -> setObsFitParameters(14,parsFobs14);

  vector<double> parsFobs15; 
  parsFobs15.push_back(.5);
  parsFobs15.push_back(.1);
  parsFobs15.push_back(4);
  myLRhelper -> setObsFitParameters(15,parsFobs15);

  vector<double> parsFobs17; 
  parsFobs17.push_back(2);
  parsFobs17.push_back(-.06);
  parsFobs17.push_back(-10);
  parsFobs17.push_back(7);
  parsFobs17.push_back(.03);
  parsFobs17.push_back(-30);
  parsFobs17.push_back(-2);
  myLRhelper -> setObsFitParameters(17,parsFobs17); 

  vector<double> parsFobs20; 
  parsFobs20.push_back(2000);
  parsFobs20.push_back(-.06);
  parsFobs20.push_back(-100);
  parsFobs20.push_back(-2000);
  myLRhelper -> setObsFitParameters(20,parsFobs20);
  
  vector<double> parsFobs21; 
  parsFobs21.push_back(.4);
  parsFobs21.push_back(.3);
  parsFobs21.push_back(2);
  myLRhelper -> setObsFitParameters(21,parsFobs21); 

  vector<double> parsFobs23; 
  parsFobs23.push_back(.6);
  parsFobs23.push_back(-.03);
  parsFobs23.push_back(-100);
  parsFobs23.push_back(.2);
  myLRhelper -> setObsFitParameters(23,parsFobs23);

  vector<double> parsFobs24; 
  parsFobs24.push_back(.6);
  parsFobs24.push_back(.1);
  parsFobs24.push_back(2);
  myLRhelper -> setObsFitParameters(24,parsFobs24);
 
  vector<double> parsFobs25; 
  parsFobs25.push_back(.6);
  parsFobs25.push_back(2);
  parsFobs25.push_back(1);
  myLRhelper -> setObsFitParameters(25,parsFobs25);

  vector<double> parsFobs26; 
  parsFobs26.push_back(.09);
  parsFobs26.push_back(-.04);
  parsFobs26.push_back(10);
  parsFobs26.push_back(.06);
  myLRhelper -> setObsFitParameters(26,parsFobs26);
 
  vector<double> parsFobs27; 
  parsFobs27.push_back(.2);
  parsFobs27.push_back(-1);
  parsFobs27.push_back(1);
  myLRhelper -> setObsFitParameters(27,parsFobs27);

  vector<double> parsFobs28; 
  parsFobs28.push_back(.1);
  parsFobs28.push_back(2);
  parsFobs28.push_back(.6);
  myLRhelper -> setObsFitParameters(28,parsFobs28);

  vector<double> parsFobs29; 
  parsFobs29.push_back(-10);
  parsFobs29.push_back(.02);
  parsFobs29.push_back(-400);
  parsFobs29.push_back(.9);
  myLRhelper -> setObsFitParameters(29,parsFobs29);

  vector<double> parsFobs30; 
  parsFobs30.push_back(.5);
  parsFobs30.push_back(-.4);
  parsFobs30.push_back(4);
  myLRhelper -> setObsFitParameters(30,parsFobs30);
  
  vector<double> parsFobs32; 
  parsFobs32.push_back(-.02);
  parsFobs32.push_back(.3);
  parsFobs32.push_back(-30);
  parsFobs32.push_back(1);
  myLRhelper -> setObsFitParameters(32,parsFobs32);

  vector<double> parsFobs39; 
  parsFobs39.push_back(.05);
  parsFobs39.push_back(.2);
  parsFobs39.push_back(1);
  myLRhelper -> setObsFitParameters(39,parsFobs39);

  vector<double> parsFobs40; 
  parsFobs40.push_back(.05);
  parsFobs40.push_back(.3);
  parsFobs40.push_back(1);
  myLRhelper -> setObsFitParameters(40,parsFobs40);

  vector<double> parsFobs46; 
  parsFobs46.push_back(.5);
  parsFobs46.push_back(1);
  parsFobs46.push_back(2);  
  parsFobs46.push_back(.1);
  parsFobs46.push_back(1.5);
  parsFobs46.push_back(1);
  myLRhelper -> setObsFitParameters(46,parsFobs46);
 
  vector<double> parsFobs47; 
  parsFobs47.push_back(.5);
  parsFobs47.push_back(1);
  parsFobs47.push_back(2);  
  myLRhelper -> setObsFitParameters(47,parsFobs47);
  
  vector<double> parsFobs48; 
  parsFobs48.push_back(.7);
  parsFobs48.push_back(2);
  parsFobs48.push_back(2);  
  parsFobs48.push_back(-.1);  
  myLRhelper -> setObsFitParameters(48,parsFobs48);
 
  vector<double> parsFobs51; 
  parsFobs51.push_back(.06);
  parsFobs51.push_back(.6);
  parsFobs51.push_back(1);  
  myLRhelper -> setObsFitParameters(51,parsFobs51);
  
  vector<double> parsFobs52; 
  parsFobs52.push_back(.06);
  parsFobs52.push_back(.7);
  parsFobs52.push_back(.7);  
  myLRhelper -> setObsFitParameters(52,parsFobs52);

  vector<double> parsFobs54; 
  parsFobs54.push_back(.5);
  parsFobs54.push_back(1);
  parsFobs54.push_back(2);    
  parsFobs54.push_back(.1);
  parsFobs54.push_back(1.5);
  parsFobs54.push_back(1); 
  myLRhelper -> setObsFitParameters(54,parsFobs54);

  vector<double> parsFobs55; 
  parsFobs55.push_back(.6);
  parsFobs55.push_back(2.4);
  parsFobs55.push_back(2);  
  myLRhelper -> setObsFitParameters(55,parsFobs55);
  
  vector<double> parsFobs56; 
  parsFobs56.push_back(.6);
  parsFobs56.push_back(3);
  parsFobs56.push_back(2);  
  parsFobs56.push_back(-.06);
  myLRhelper -> setObsFitParameters(56,parsFobs56);

  vector<double> parsFobs63; 
  parsFobs63.push_back(.3);
  parsFobs63.push_back(.1);
  parsFobs63.push_back(.5);    
  parsFobs63.push_back(.04);
  parsFobs63.push_back(.9);
  parsFobs63.push_back(.07);   
  parsFobs63.push_back(.08);
  parsFobs63.push_back(.3);
  parsFobs63.push_back(.04);
  myLRhelper -> setObsFitParameters(63,parsFobs63);

  vector<double> parsFobs65; 
  parsFobs65.push_back(.4);
  parsFobs65.push_back(40);
  parsFobs65.push_back(30);    
  parsFobs65.push_back(.7);
  parsFobs65.push_back(-20);
  parsFobs65.push_back(-33);   
  myLRhelper -> setObsFitParameters(65,parsFobs65);



  // normalize the S and B histograms to construct the pdf's
  //myLRhelper -> normalizeSandBhists();
  
  // produce and fit the S/S+N histograms
  myLRhelper -> makeAndFitSoverSplusBHists();
  myLRhelper->setXlabels(xlabels);

  // store histograms and fits in root-file
  myLRhelper -> storeToROOTfile(dir+outputFile);
       cout << "store\n";

  // make some control plots and put them in a .ps file
  myLRhelper -> storeControlPlots(dir+outputPSfile);
//   pair<double, double> match = myLRhelper->getBnumbers();
//   cout << "Matched B jets   : " <<match.first <<endl;
//   cout << "Non-matched jets : "<<match.second<<endl;
//   cout << "Purity           : "<< match.first/(match.first+match.second)<<endl;
//   for(int j = 1; j < nrJetCombObs+1; ++j){
//     TString name = "Obs"; name += j;
//     myLRhelper->singlePlot(name,j,"eps");
//   }
}


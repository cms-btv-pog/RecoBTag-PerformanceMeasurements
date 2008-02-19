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
const TString dir = "/msa1/jmmaes/CMSSW_1_6_9/src/TopQuarkAnalysis/Examples/crab/test/";
const  TString  inputFile    = "plots.root";
const  TString  outputFile   = "final_obs_plots.root";
const  TString  outputPSfile = "LRSignalSelAllObs.ps";

//observable histogram variables
const  int      nrSignalSelObs  		= 15;
const  int      SignalSelObs[nrSignalSelObs] 	= {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

//observable fit functions
const char*     SignalSelObsFits[nrSignalSelObs]= {           
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs1
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs2
   "([0]+[3]*x+[4]*x*x)/(1 + 1/exp([1]*([2] - x)))", //obs3
   "([0]+[3]*x+[4]*x*x)/(1 + 1/exp([1]*([2] - x)))", //obs4
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs5
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs6
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs7
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs8
   "pol4", //obs9
   "pol4", //obs10
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs11
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs12
   "[0]/(1 + 1/exp([1]*([2] - x)))", //obs13
   "pol4", //obs14
   "pol4" //obs15
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
  for(int j = 0; j < nrSignalSelObs; j++){
    obsNrs.push_back(SignalSelObs[j]);
    obsFits.push_back(SignalSelObsFits[j]);
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
//   for(int j = 1; j < nrSignalSelObs+1; ++j){
//     TString name = "Obs"; name += j;
//     myLRhelper->singlePlot(name,j,"eps");
//   }
}


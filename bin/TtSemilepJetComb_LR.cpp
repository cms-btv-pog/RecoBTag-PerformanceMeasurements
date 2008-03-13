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
#include <TF1.h>
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
const TString dir = "/msa1/jmmaes/CMSSW_1_6_9_NOPAT/src/TopQuarkAnalysis/Examples/crab/test/";
const  TString  inputFile    = "plots_CLR_test9.root";
const  TString  outputFile   = "final_lr_plots_test9.root";
const  TString  outputPSfile = "LRplots_test9.ps";

//observable histogram variables
const  int      nrJetCombObs  		= 8;
//const  int      JetCombObs[nrJetCombObs] 	= {1,2,3,4,5,6,7,8,9,10,11,12,13};
const  int      JetCombObs[nrJetCombObs] 	= {8,29,42,46,48,54,56,65};
//const  int      JetCombObs[nrJetCombObs] 	= {1,2,3,4,7,8,9,10,11,12};

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
  
  
  // define all histograms & fit functions
  //to replace with something more elegant
  for(int j = 0; j < nrJetCombObs; j++){
    obsNrs.push_back(JetCombObs[j]);
  }
  myLRhelper = new LRHelpFunctions();
  //load plots
  myLRhelper->readObsHistsAndFits(dir+inputFile, obsNrs, true);
//  myLRhelper->recreateFitFct(obsNrs, obsFits);
//  cout << "fit functions loaded\n";
  
  // produce and fit the S/S+N histograms
  myLRhelper -> makeAndFitPurityHists();
  myLRhelper -> makeAndFitSoverSplusBHists();
     cout << "make\n";

  // store histograms and fits in root-file
  myLRhelper -> storeToROOTfile(dir+outputFile);

  // make some control plots and put them in a .ps file
  myLRhelper -> storeControlPlots(dir+outputPSfile);
//  myLRhelper->purityPlot("lh","eps");
//   pair<double, double> match = myLRhelper->getBnumbers();
//   cout << "Matched B jets   : " <<match.first <<endl;
//   cout << "Non-matched jets : "<<match.second<<endl;
//   cout << "Purity           : "<< match.first/(match.first+match.second)<<endl;
}


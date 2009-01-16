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
TString dir = "";
TString  inputFile    = "raw_lr_plots.root";
TString  outputFile   = "final_lr_plots.root";
TString  outputPSfile = "final_lr_plots.ps";

//observable histogram variables
const  int      nrSignalSelObs  		= 14;
const  int      SignalSelObs[nrSignalSelObs] 	= {1,2,3,4,5,6,7,8,9,10,11,12,13,15};
//const  int      SignalSelObs[nrSignalSelObs] 	= {1,2,3,4,7,8,9,10,11,12};

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

int main(int argc, const char* argv[]) { 
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  inputFile = dir+inputFile;
  outputPSfile = dir+outputPSfile;
  outputFile = dir + outputFile;

  // Filenames
  if (argc>=2) {
    inputFile = TString(argv[1]);
  }
  if (argc>=3) {
    outputFile = TString(argv[2])+TString(".root");
    outputPSfile = TString(argv[2])+TString(".ps");
  }

  
  // define all histograms & fit functions
  //to replace with something more elegant
  for(int j = 0; j < nrSignalSelObs; j++){
    obsNrs.push_back(SignalSelObs[j]);
  }
  myLRhelper = new LRHelpFunctions();
  //load plots
  myLRhelper->readObsHistsAndFits(inputFile, obsNrs, true);
//  myLRhelper->recreateFitFct(obsNrs, obsFits);
//  cout << "fit functions loaded\n";
  
  // produce and fit the S/S+N histograms
  myLRhelper -> makeAndFitPurityHists();
  myLRhelper -> makeAndFitSoverSplusBHists();
     cout << "make\n";

  // store histograms and fits in root-file
  myLRhelper -> storeToROOTfile(outputFile);

  // make some control plots and put them in a .ps file
  myLRhelper -> storeControlPlots(outputPSfile);
  myLRhelper->purityPlot("lh","eps");
  myLRhelper->getBnumbersFromLRC(0.8);
}


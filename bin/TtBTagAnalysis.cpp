#include <iostream>
#include <cassert>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TLegend.h>
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
#include "RecoBTag/PerformanceMeasurements/interface/TtEtEtaHistoCollector.h"
#include "TopQuarkAnalysis/TopTools/test/tdrstyle.C"

using namespace std;



///////////////////////
// Constants         //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//input files
//const TString dir = "/localscratch/s/speer/top_eff/data/new_skim/plots_xsect_syst/";
//  const TString dir = "/localscratch/s/speer/top_eff/data/new_skim/plots_smear_minus_3/";
 const TString dir = "/localscratch/s/speer/top_eff/data/new_skim/plots_obs/";

// const  TString  mcFile    = "/localscratch/s/speer/top_eff/data/basic_skim_all/ttj0/plots.root";
// const  TString  dataFile    = "/localscratch/s/speer/top_eff/data/basic_skim_all/ttj0/plots.root";

// const  TString  mcFile    = "/localscratch/s/speer/top_eff/data/basic_skim_all/plots_obs/raw_btag_plots.root";
// const  TString  dataFile    = "/localscratch/s/speer/top_eff/data/basic_skim_all/plots_obs/raw_btag_plots.root";
TString  mcFile    = "/localscratch/s/speer/top_eff/data/new_skim/plots_obs/raw_btag_plots.root";
TString  dataFile    = "raw_btag_plots.root";

TString  outputFile   = "final_BTag_plots.root";
TString  outputPSfile = "LRplots.ps";

bool doEps = false;
TString epsType = "eps";

bool compare = false;
TString  referenceFile    = "/localscratch/s/speer/top_eff/data/new_skim/plots_obs/final_BTag_plots.root";

bool useSystFiles = true;
const int syst = 1;
const char*     systFiles[syst]= { "/localscratch/s/speer/top_eff/data/new_skim/plots_xsect_syst/final_BTag_plots.root" };


//observable histogram variables
const  int      nrSignalSelObs  		= 13;
const  int      SignalSelObs[nrSignalSelObs] 	= {1,2,3,4,5,6,7,8,9,10,11,12,13};
//const  int      SignalSelObs[nrSignalSelObs] 	= {1,2,3,4,7,8,9,10,11,12};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const bool debug = true;
const int wp = 3;
const char*     points[wp]= { "trackCountingJetTags_2.3" , "trackCountingJetTags_5.3", "trackCountingJetTags3th_4.8"};

const int lrBins   = 20;
const double lrMin   = 0.;
const double lrMax   = 0.7;




//
// Global variables
//
LRHelpFunctions *myLRhelper;
vector<int> obsNrs;
vector<double> obsMin,obsMax;
vector<const char*> obsFits;


void compute(TH1F *b_effH, TH1F *tagged_jetsH, TH1F *b_fractionH, TH1F *eff_nonBH);	
void computeStatError(TH1F *b_eff_errH, TH1F *tagged_jetsH, TH1F *b_fractionH);
void computeSystError(TH1F *b_eff_errH, TH1F *tagged_jetsH, TH1F *b_fractionH, TH1F *eff_nonBH);
void addSystError(TH1F *b_eff_errH, TH1F *new_systH);

void totalError(TH1F *b_effH, TH1F* b_eff_totH, TH1F *b_eff_statH, TH1F *b_eff_systH);
void eps(TH1F * histo, TH1F * histo2, TH1F * histo3, TLegend* leg);
void epsEff(vector<TH1F *>, TLegend* leg);

int main( int argc, const char* argv[]) {


  setTDRStyle();

// Filenames
  if (argc>=2) {
    mcFile = dataFile    = argv[1];
  }
  if (argc>=3) {
    outputFile = argv[2];
  }

  cout << "Input file  : " << mcFile<<endl;
  cout << "Output file : " << dir<<outputFile<<endl;

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();


//   // define all histograms & fit functions
//   myLRhelper = new LRHelpFunctions(lrBins, lrMin, lrMax, lrFits);
//   myLRhelper->readObsHistsAndFits(obsFileName, obsNrs, false);
//   if (debug) cout << "LHR observables read" << endl;


  // eta and pt ranges (bins for differential plots)
  vector<double> etaRanges;
  vector<double> etRanges;
  etaRanges.push_back(0.);
  etaRanges.push_back(1.4);
  etaRanges.push_back(2.4);
  etRanges.push_back(30.);
  etRanges.push_back(50.);
  etRanges.push_back(80.);
  etRanges.push_back(120.);  
  etRanges.push_back(500.);  

//   etRanges.push_back(160.);
//   etRanges.push_back(250.);


  if (debug) cout << "Load Data Plots from file " << dir<<dataFile<<endl;


  TFile *theFile = new TFile (dir+dataFile) ;
  TH1::AddDirectory(kFALSE);

  vector<TtEtEtaHistoCollector *> taggedJetsHistos;
  vector<TString> bTagger;

  for (int iWP = 0; iWP != wp; ++iWP) {
    bTagger.push_back( TString(points[iWP]) );

    if (debug) cout << "Load B-tagged plots for Working Point " << bTagger[iWP]<<endl;
    taggedJetsHistos.push_back(new TtEtEtaHistoCollector("TaggedJets_"+bTagger[iWP], etRanges, etaRanges,
	lrMin, lrMax, lrBins, true) ) ;
  }

  if (debug) cout << "Load plots for AllEvents" << endl;

  TtEtEtaHistoCollector * allJetsHistos = new TtEtEtaHistoCollector(TString("AllJets"), etRanges, etaRanges,
      lrMin, lrMax, lrBins, true);

  theFile->Close();

  if (debug) cout << "Load MC Plots from file " << mcFile<<endl;
  theFile = new TFile (mcFile) ;

  vector<TtEtEtaHistoCollector*> taggedLightJetsHistos, taggedBJetsHistos;
  for (int iWP = 0; iWP != wp; ++iWP) {
    if (debug) cout << "Load MC plots for Working Point " << bTagger[iWP]<<endl;
    taggedLightJetsHistos.push_back(new TtEtEtaHistoCollector("TaggedLightJets_"+bTagger[iWP], etRanges, etaRanges,
      lrMin, lrMax, lrBins, true) ) ;
    taggedBJetsHistos.push_back(new TtEtEtaHistoCollector("TaggedBJets_"+bTagger[iWP], etRanges, etaRanges,
      lrMin, lrMax, lrBins, true) ) ;
  }

  TtEtEtaHistoCollector * allJetsMCHistos = new TtEtEtaHistoCollector(TString("AllJets"), etRanges, etaRanges,
      lrMin, lrMax, lrBins, true);

  TtEtEtaHistoCollector * lrSHistos = new TtEtEtaHistoCollector(TString("LRsignal"), etRanges, etaRanges,
      lrMin, lrMax, lrBins, true);

  TtEtEtaHistoCollector * lrBHistos = new TtEtEtaHistoCollector(TString("LRbackground"), etRanges, etaRanges,
      lrMin, lrMax, lrBins, true);

  TtEtEtaHistoCollector * allLightHistos = new TtEtEtaHistoCollector(TString("All_LightJets"), etRanges, etaRanges,
      lrMin, lrMax, lrBins, true);

  theFile->Close();

  if (debug) cout << "All histos loaded! " << endl;
delete theFile;

//////////// Process

  TtEtEtaHistoCollector* intAllJetsHistos = allJetsHistos->buildIntegratedCollector();
  TtEtEtaHistoCollector* intAllJetsMCHistos = allJetsMCHistos->buildIntegratedCollector();

  TtEtEtaHistoCollector* intLrSHistos  =  lrSHistos->buildIntegratedCollector();
  TtEtEtaHistoCollector* intLrBHistos  =  lrBHistos->buildIntegratedCollector();
  TtEtEtaHistoCollector* intAllLightHistos = allLightHistos->buildIntegratedCollector();

  TtEtEtaHistoCollector* b_fraction = intLrSHistos->buildRatioCollector(TString("B_Purity"),
			*intAllJetsMCHistos);  	// w_b
  if (debug) cout << "b_fraction built " << endl;

  vector<TtEtEtaHistoCollector*> eff_nonB;		// eff_0
  vector<TtEtEtaHistoCollector*> tagged_fraction;  // x_tag
  vector<TtEtEtaHistoCollector*> eff_B;		// eff_B
  vector<TtEtEtaHistoCollector*> eff_B_stat;	// Stat error on eff_B
  vector<TtEtEtaHistoCollector*> eff_B_syst;	// syst error on eff_B
  vector<TtEtEtaHistoCollector*> eff_B_tot;	// total error on eff_B
  vector<TtEtEtaHistoCollector*> eff_B_MC;	// eff_B from MC

  vector<TtEtEtaHistoCollector*> intTaggedJetsHistos;	 // nbr of tagged jets
  vector<TtEtEtaHistoCollector*> intTaggedLightJetsHistos; // nbr of tagged light jets
  vector<TtEtEtaHistoCollector*> intTaggedBJetsHistos;	 // nbr of tagged B jets

  for (int iWP = 0; iWP != wp; ++iWP) {
    intTaggedJetsHistos.push_back(taggedJetsHistos[iWP]->buildIntegratedCollector());
    intTaggedLightJetsHistos.push_back(taggedLightJetsHistos[iWP]->buildIntegratedCollector());
    intTaggedBJetsHistos.push_back(taggedBJetsHistos[iWP]->buildIntegratedCollector());

    if (debug) cout << "Compute tagged_fraction for Working Point " << bTagger[iWP]<<endl;
    tagged_fraction.push_back( intTaggedJetsHistos[iWP]->buildRatioCollector(
		TString("tagged_fraction_")+bTagger[iWP], *intAllJetsHistos) ) ;

    if (debug) cout << "Compute eff_nonB for Working Point " << bTagger[iWP]<<endl;
    eff_nonB.push_back(intTaggedLightJetsHistos[iWP]->buildRatioCollector(
		TString("eff_nonB_")+bTagger[iWP], *intAllLightHistos) ) ;

    if (debug) cout << "Compute eff_B_MC for Working Point " << bTagger[iWP]<<endl;
    eff_B_MC.push_back(intTaggedBJetsHistos[iWP]->buildRatioCollector(
		TString("eff_B_MC_")+bTagger[iWP], *intLrSHistos) ) ;

    if (debug) cout << "Book eff_B for Working Point " << bTagger[iWP]<<endl;
    eff_B.push_back(new TtEtEtaHistoCollector("eff_B_"+bTagger[iWP], etRanges, etaRanges,
      lrMin, lrMax, lrBins, false) ) ;
    eff_B_stat.push_back(new TtEtEtaHistoCollector("eff_B_statErr_"+bTagger[iWP], etRanges, etaRanges,
      lrMin, lrMax, lrBins, false) ) ;
    eff_B_syst.push_back(new TtEtEtaHistoCollector("eff_B_systErr_"+bTagger[iWP], etRanges, etaRanges,
      lrMin, lrMax, lrBins, false) ) ;
    eff_B_tot.push_back(new TtEtEtaHistoCollector("eff_B_totErr_"+bTagger[iWP], etRanges, etaRanges,
      lrMin, lrMax, lrBins, false) ) ;

    if (debug) cout << "Compute eff_B for Working Point " << bTagger[iWP]<<endl;
    for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B[iWP]->histoBegin();
	bin!=eff_B[iWP]->histoEnd(); ++bin) {
      compute(bin->second, tagged_fraction[iWP]->getHisto(*bin->first),
	b_fraction->getHisto(*bin->first), eff_nonB[iWP]->getHisto(*bin->first) );
    }
    for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B_stat[iWP]->histoBegin();
	bin!=eff_B_stat[iWP]->histoEnd(); ++bin) {
      computeStatError(bin->second, tagged_fraction[iWP]->getHisto(*bin->first),
	b_fraction->getHisto(*bin->first));
    }
    for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B_syst[iWP]->histoBegin();
	bin!=eff_B_syst[iWP]->histoEnd(); ++bin) {
      computeSystError(bin->second, tagged_fraction[iWP]->getHisto(*bin->first),
	b_fraction->getHisto(*bin->first), eff_nonB[iWP]->getHisto(*bin->first) );
    }

  }
////////////

  vector<TtEtEtaHistoCollector*> difference;	// Difference
  if (compare) {
    vector<TtEtEtaHistoCollector*> eff_B_ref;	// eff_B of reference
    if (debug) cout << "Comparison with reference (for systematics): " << referenceFile << endl;
    theFile = new TFile (referenceFile) ;
    for (int iWP = 0; iWP != wp; ++iWP) {
      eff_B_ref.push_back(new TtEtEtaHistoCollector("eff_B_"+bTagger[iWP], etRanges, etaRanges,
	lrMin, lrMax, lrBins, true) ) ;
      difference.push_back(eff_B_ref[iWP]->buildDiffCollector(
		TString("difference_")+bTagger[iWP], *eff_B[iWP]) ) ;
    }
    theFile->Close();
    delete theFile;
  }
////////////
  if (useSystFiles) {
    for (int isyst = 0; isyst != syst; ++isyst) {
      if (debug) cout << "Systematics from file: " << systFiles[isyst] << endl;
      theFile = new TFile (systFiles[isyst]) ;
      for (int iWP = 0; iWP != wp; ++iWP) {
	TtEtEtaHistoCollector * diff = new TtEtEtaHistoCollector("difference_"+bTagger[iWP], etRanges, etaRanges,
	  lrMin, lrMax, lrBins, true);
	for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B_syst[iWP]->histoBegin();
	    bin!=eff_B_syst[iWP]->histoEnd(); ++bin) {
          addSystError(bin->second, diff->getHisto(*bin->first));
	}
	delete diff;
      }
      theFile->Close();
      delete theFile;

    }
  }

////////////
  for (int iWP = 0; iWP != wp; ++iWP) {
    for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B[iWP]->histoBegin();
	bin!=eff_B[iWP]->histoEnd(); ++bin) {
      totalError(bin->second, eff_B_tot[iWP]->getHisto(*bin->first),
	eff_B_stat[iWP]->getHisto(*bin->first), eff_B_syst[iWP]->getHisto(*bin->first));
    }

    eff_B[iWP]->buildDifferentialPlots(7);
    cout << "MC diffPlot\n";
    eff_B_MC[iWP]->buildDifferentialPlots(1);

  }

//////////// Save all histos.

  TFile * theFile2 = new TFile (dir+TString (outputFile) , "RECREATE" ) ;
  theFile2->cd();

  if (debug) cout << "write Histos" << endl;
  allJetsHistos->write();
  intAllJetsMCHistos->write();
  lrSHistos->write();
  lrBHistos->write();
  allLightHistos->write();

  intAllJetsHistos->write();
  intLrSHistos->write();
  intLrBHistos->write();
  intAllLightHistos->write();
  b_fraction->write();


  for (unsigned int iWP = 0; iWP != bTagger.size(); ++iWP) {
    eff_nonB[iWP]->setLabels(TString("Cut on Combined LR"), TString("Mistag rate"));
    eff_B[iWP]->setLabels(TString("Cut on Combined LR"), TString("b-tagging efficiency"));
    eff_B_stat[iWP]->setLabels(TString("Cut on Combined LR"), TString("Stat. uncert. b-tag eff."));
    eff_B_syst[iWP]->setLabels(TString("Cut on Combined LR"), TString("Syst. uncert. b-tag eff."));
    eff_B_tot[iWP]->setLabels(TString("Cut on Combined LR"), TString("Tot. uncert. b-tag eff."));

    taggedJetsHistos[iWP]->write();
    taggedLightJetsHistos[iWP]->write();
    intTaggedJetsHistos[iWP]->write();
    intTaggedLightJetsHistos[iWP]->write();
    taggedBJetsHistos[iWP]->write();
    intTaggedBJetsHistos[iWP]->write();

    tagged_fraction[iWP]->write();
    eff_nonB[iWP]->write();
    eff_B[iWP]->write();
    eff_B_stat[iWP]->write();
    eff_B_syst[iWP]->write();
    eff_B_tot[iWP]->write();
    eff_B_MC[iWP]->write();
    if (compare) difference[iWP]->write();
    if (doEps) {
      eff_nonB[iWP]->setLabels(TString("Cut on Combined LR"), TString("Mistag rate"));
      eff_nonB[iWP]->eps(epsType);

      eff_B[iWP]->setRange(0.0,1.1);
      
   eff_B_stat[iWP]->setLabels(TString("Cut on Combined LR"), TString("Uncertainty"));
   eff_B_syst[iWP]->setLabels(TString("Cut on Combined LR"), TString("Uncertainty"));
   eff_B_tot[iWP]->setLabels(TString("Cut on Combined LR"), TString("Uncertainty"));
   eff_B_tot[iWP]->setRange(0., 0.08);

	for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B_tot[iWP]->histoBegin();
	    bin!=eff_B_tot[iWP]->histoEnd(); ++bin) {
  TLegend* leg = new TLegend(0.7,0.8,0.90,0.93);
  leg->AddEntry(bin->second,"Total","p");
  leg->AddEntry(eff_B_syst[iWP]->getHisto(*bin->first),"Systematic","p");
  leg->AddEntry(eff_B_stat[iWP]->getHisto(*bin->first),"Statistic","p");
           eps(bin->second, eff_B_syst[iWP]->getHisto(*bin->first), eff_B_stat[iWP]->getHisto(*bin->first), leg);
	 }
   eff_B[iWP]->setLabels(TString("Cut on Combined LR"), TString("b-tagging efficiency"));
    eff_B[iWP]->eps(epsType);
//       eff_B_stat[iWP]->eps(epsType);
//       eff_B_syst[iWP]->eps(epsType);
    }
  }
// 	for (TtEtEtaHistoCollector::BinHistoMapIter bin = eff_B[0]->histoBegin();
// 	    bin!=eff_B[0]->histoEnd(); ++bin) {
//            epsEff(bin->second, eff_B[1]->getHisto(*bin->first), eff_B[2]->getHisto(*bin->first));
// 	 }

  TtEtEtaHistoCollector::diffPlotVect diffA = eff_B[0]->getDifferentialPlots();
  TtEtEtaHistoCollector::diffPlotVect diffB = eff_B[1]->getDifferentialPlots();
  TtEtEtaHistoCollector::diffPlotVect diffC = eff_B[2]->getDifferentialPlots();
  TtEtEtaHistoCollector::diffPlotVect diffD = eff_B_MC[0]->getDifferentialPlots();
  TtEtEtaHistoCollector::diffPlotVect diffE = eff_B_MC[1]->getDifferentialPlots();
  TtEtEtaHistoCollector::diffPlotVect diffF = eff_B_MC[2]->getDifferentialPlots();
  
  for (unsigned int i = 0; i != diffA.size(); ++i) {
  TLegend* leg = new TLegend(0.7,0.8,0.90,0.93);
  leg->AddEntry(diffA[i]->getHisto(),"Loose","p");
  leg->AddEntry(diffB[i]->getHisto(),"Medium","p");
  leg->AddEntry(diffC[i]->getHisto(),"Tight","p");
  
  leg->AddEntry(diffD[i]->getHisto(),"Loose MC","p");
  leg->AddEntry(diffE[i]->getHisto(),"Medium MC","p");
  leg->AddEntry(diffF[i]->getHisto(),"Tight MC","p");
  vector<TH1F *> histoVec;
  histoVec.push_back(diffA[i]->getHisto());
  histoVec.push_back(diffB[i]->getHisto());
  histoVec.push_back(diffC[i]->getHisto());
  histoVec.push_back(diffD[i]->getHisto());
  histoVec.push_back(diffE[i]->getHisto());
  histoVec.push_back(diffF[i]->getHisto());
    epsEff(histoVec, leg);
  }

  if (debug) cout << "done " << endl;


  theFile2->Close();

//Plots:
  
//   eff_B_stat->setlabels(TString("Cut on Combined LR"), TString(""));
//   eff_B_syst->setlabels(TString("Cut on Combined LR"), TString(""));
//   ->setlabels(TString("Cut on Combined LR"), TString(""));
//   ->setlabels(TString("Cut on Combined LR"), TString(""));

  
  
}

void compute(TH1F *b_effH, TH1F *tagged_jetsH, TH1F *b_fractionH, TH1F *eff_nonBH)
{
//  cout <<"Enter compute "<<b_effH->GetName()<<endl;
//  cout <<"Enter compute "<<tagged_jetsH->GetName()<<endl;
//  cout <<"Enter compute "<<b_fractionH->GetName()<<endl;
//  cout <<"Enter compute "<<eff_nonBH->GetName()<<endl;

  int bins = b_effH->GetNbinsX();
// cout <<"Bins "<<bins<<endl;
  double eff, tagged_jets, b_fraction, eff_nonB;
  for(int bin=0; bin<=bins+1; ++bin){
    tagged_jets = tagged_jetsH->GetBinContent(bin);
// cout <<"tagged_jets "<<tagged_jets<<endl;
    b_fraction = b_fractionH->GetBinContent(bin);
    eff_nonB = eff_nonBH->GetBinContent(bin);
// cout <<"eff_nonB "<<eff_nonB<<endl;
//  cout <<"bin "<< bin <<" "<< tagged_jets<<" "<<eff_nonB <<" "<< b_fraction <<" "<<
//  ( tagged_jets - eff_nonB * (1-b_fraction) ) /b_fraction<<endl;
    if (b_fraction!=0.) {
      eff = ( tagged_jets - eff_nonB * (1-b_fraction) ) /b_fraction;
      b_effH->SetBinContent(bin, eff );
    }

// cout <<"eff "<<eff<<endl;
  }
}

void computeStatError(TH1F *b_eff_errH, TH1F *tagged_jetsH, TH1F *b_fractionH)
{

  int bins = b_eff_errH->GetNbinsX();
  double b_fraction, tagged_jets_err;
  for(int bin=0; bin<=bins+1; ++bin){
    tagged_jets_err = tagged_jetsH->GetBinError(bin);
    b_fraction = b_fractionH->GetBinContent(bin);
    if (b_fraction!=0.) 
      b_eff_errH->SetBinContent(bin, tagged_jets_err/ b_fraction);
//     cout << bin <<" "<<tagged_jets_err<<" "<< b_fraction<<endl;
  }
}

void computeSystError(TH1F *b_eff_errH, TH1F *tagged_jetsH, TH1F *b_fractionH, TH1F *eff_nonBH)
{
 cout <<"Enter computeSyst "<<b_eff_errH->GetName()<<endl;
 cout <<"Enter computeSyst "<<tagged_jetsH->GetName()<<endl;
 cout <<"Enter computeSyst "<<b_fractionH->GetName()<<endl;
 cout <<"Enter computeSyst "<<eff_nonBH->GetName()<<endl;
  int bins = b_eff_errH->GetNbinsX();
  double b_fraction, eff_nonB;
  double eff_nonB_err;
  double systErr;
  for(int bin=0; bin<=bins+1; ++bin){
    b_fraction = b_fractionH->GetBinContent(bin);
    eff_nonB = eff_nonBH->GetBinContent(bin);

    systErr = 0.;
    eff_nonB_err = 0.2* eff_nonB;
    if (b_fraction!=0.) {
      systErr+= eff_nonB_err * (1/b_fraction-1);
      b_eff_errH->SetBinContent(bin, systErr);
    }
      cout << bin <<" "<<eff_nonB<<" "<<systErr<<" "<< b_fraction<<endl;
  }
}

void addSystError(TH1F *b_eff_errH, TH1F *new_systH)
{
  cout <<"Enter addSystError "<<b_eff_errH->GetName()<<endl;
  cout <<"Enter addSystError "<<new_systH->GetName()<<endl;
  int bins = b_eff_errH->GetNbinsX();
  double error;
  for(int bin=0; bin<=bins+1; ++bin){
    error = sqrt(b_eff_errH->GetBinContent(bin) * b_eff_errH->GetBinContent(bin)
	+ new_systH->GetBinContent(bin) * new_systH->GetBinContent(bin) );
    b_eff_errH->SetBinContent(bin, error);
//       cout << bin <<" "<<eff_nonB<<" "<<systErr<<" "<< b_fraction<<endl;
  }
}

void totalError(TH1F *b_effH, TH1F* b_eff_totH, TH1F *b_eff_statH, TH1F *b_eff_systH)
{
  int bins = b_effH->GetNbinsX();
  double error;
  for(int bin=0; bin<=bins+1; ++bin){
    error = sqrt(b_eff_statH->GetBinContent(bin) * b_eff_statH->GetBinContent(bin)
	+ b_eff_systH->GetBinContent(bin) * b_eff_systH->GetBinContent(bin) );
    b_effH->SetBinError(bin, error);
    b_eff_totH->SetBinContent(bin, error);
  }
}

void eps(TH1F * histo1, TH1F * histo2, TH1F * histo3, TLegend* leg)
{
  TStyle *tdrStyle = gROOT->GetStyle("tdrStyle");
  tdrStyle->SetOptFit(0);
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptTitle(0);

  TCanvas c2("c2","",600,600);
  histo1->SetMarkerColor(1);
histo2->SetMarkerColor(2);
histo3->SetMarkerColor(3);
  histo1->Draw("P");
  
  histo2->Draw("Psame");
  histo3->Draw("Psame");
  leg->Draw();
  cout << histo1->GetName()+TString(".")+epsType<<endl;
  c2.Print(histo1->GetName()+TString(".")+epsType);
}
void epsEff(vector<TH1F *> histV, TLegend* leg)
{
  TStyle *tdrStyle = gROOT->GetStyle("tdrStyle");
  tdrStyle->SetOptFit(0);
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptTitle(0);
  TCanvas c2("c2","",600,600);

  for (vector<TH1F *>::iterator i = histV.begin(); i!=histV.end(); ++i){
  (*i)->GetYaxis()->SetRangeUser(0.2,1.01);
  (*i)->SetMarkerColor(1+ (i - histV.begin()));
  (*i)->SetMarkerStyle(20+ (i - histV.begin()));
  (*i)->SetMarkerSize(1.);
  if (i == histV.begin()) (*i)->Draw("P");
    else (*i)->Draw("Psame");
  }
  leg->Draw();
  c2.Print(histV[0]->GetName()+TString(".")+epsType);
}

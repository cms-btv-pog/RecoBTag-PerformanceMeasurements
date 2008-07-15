#include "RecoBTag/PerformanceMeasurements/interface/TtDifferentialPlot.h"
#include "RecoBTag/PerformanceMeasurements/interface/TtEtEtaHistoCollector.h"

// #include "RecoBTag/Analysis/interface/TtDifferentialPlot.h"
// #include "RecoBTag/Analysis/interface/EffPurFromHistos.h"
// #include "RecoBTag/Analysis/interface/Tools.h"
#include "FWCore/Utilities/interface/CodedException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TF1.h"
#include "TStyle.h"

#include <iostream>
using namespace std ;



TtDifferentialPlot::TtDifferentialPlot (TtEtEtaHistoCollector * aCollector,
	ConstVarType constVariable, const TString & tagName) :
	processed(false), theCollector(aCollector), constVar(constVariable), 
	constVariableName( "" ), diffVariableName( "" ),
	constVariableValue( 999.9 , 999.9 ), commonName( "Beff_"+tagName),
	theHisto(0) {}


TtDifferentialPlot::~TtDifferentialPlot () {
  if (processed) {
    delete theHisto;
  }
}



void TtDifferentialPlot::write() {
  if (processed) {
//   gFile->cd();
//   gFile->mkdir(commonName);
//   gFile->cd(commonName);
    theHisto->Write() ;
  }
}



// void TtDifferentialPlot::plot (TCanvas & thePlotCanvas ) {
// 
// //   thePlotCanvas = new TCanvas(  commonName ,
// // 				commonName ,
// // 				btppXCanvas , btppYCanvas ) ;
// //
// //   if ( !btppTitle ) gStyle->SetOptTitle ( 0 ) ;
// 
//   if (!processed) return;
// //fixme:
//   bool btppNI = false;
//   bool btppColour = true;
// 
//   thePlotCanvas.SetFillColor ( 0 ) ;
//   thePlotCanvas.cd ( 1 ) ;
//   gPad->SetLogy  ( 1 ) ;
//   gPad->SetGridx ( 1 ) ;
//   gPad->SetGridy ( 1 ) ;
// 
//   int col_b   ;
//   int col_c   ;
//   int col_g   ;
//   int col_dus ;
//   int col_ni  ;
// 
//   int mStyle_b   ;
//   int mStyle_c   ;
//   int mStyle_g   ;
//   int mStyle_dus ;
//   int mStyle_ni  ;
// 
//   // marker size (same for all)
//   float mSize = 1.5 ;
// 
//   if ( btppColour ) {
//     col_b    = 2 ;
//     col_c    = 6 ;
//     col_g    = 3 ;
//     col_dus  = 4 ;
//     col_ni   = 5 ;
//     mStyle_b   = 20 ;
//     mStyle_c   = 20 ;
//     mStyle_g   = 20 ;
//     mStyle_dus = 20 ;
//     mStyle_ni  = 20 ;
//   }
//   else {
//     col_b    = 1 ;
//     col_c    = 1 ;
//     col_g    = 1 ;
//     col_dus  = 1 ;
//     col_ni   = 1 ;
//     mStyle_b   = 12 ;
//     mStyle_c   = 22 ;
//     mStyle_g   = 29 ;
//     mStyle_dus = 20 ;
//     mStyle_ni  = 27 ;
//   }
// 
//   // for the moment: plot b (to see what the constant b-efficiency is), c, g, uds
//   // b in red
//   // No, do not plot b (because only visible for the soft leptons)
//   // theDifferentialHistoB_b   -> GetXaxis()->SetTitle ( diffVariableName ) ;
//   // theDifferentialHistoB_b   -> GetYaxis()->SetTitle ( "non b-jet efficiency" ) ;
//   // theDifferentialHistoB_b   -> GetYaxis()->SetTitleOffset ( 1.25 ) ;
//   // theDifferentialHistoB_b   -> SetMaximum ( 0.4 )  ;
//   // theDifferentialHistoB_b   -> SetMinimum ( 1.e-4 )  ;
//   // theDifferentialHistoB_b   -> SetMarkerColor ( col_b ) ;
//   // theDifferentialHistoB_b   -> SetLineColor   ( col_b ) ;
//   // theDifferentialHistoB_b   -> SetMarkerSize  ( mSize ) ;
//   // theDifferentialHistoB_b   -> SetMarkerStyle ( mStyle_b ) ;
//   // theDifferentialHistoB_b   -> SetStats ( false ) ;
//   // theDifferentialHistoB_b   -> Draw ( "pe" ) ;
// 
// }
// 
// void TtDifferentialPlot::epsPlot(const TString & name)
// {
//   plot(name, ".eps");
// }
// 
// void TtDifferentialPlot::psPlot(const TString & name)
// {
//   plot(name, ".ps");
// }
// 
// void TtDifferentialPlot::plot(const TString & name, const TString & ext)
// {
//   if (!processed) return;
//    TCanvas tc(commonName, commonName);
//    plot(tc);
//    tc.Print(TString(name + commonName + ext));
// }


void TtDifferentialPlot::process (int bin) {
  setVariableName () ;
  bookHisto () ;
  fillHisto (bin) ;
  processed = true;
}


void TtDifferentialPlot::setVariableName ()
{
  if ( constVar==constETA ) {
    constVariableName  = "eta" ;
    diffVariableName   = "et"  ;
    constVariableValue = make_pair ( theBins[0]->getEtaMin() , theBins[0]->getEtaMax() ) ;
  }
  if ( constVar==constET  ) {
    constVariableName = "et"  ;
    diffVariableName  = "eta" ;
    constVariableValue = make_pair ( theBins[0]->getEtMin() , theBins[0]->getEtMax() ) ;
  }

  //edm::LogInfo("Info") 
   cout << "====>>>> TtDifferentialPlot::setVariableName() : set const/diffVariableName to : "
     << constVariableName << " / " << diffVariableName << endl
     << "====>>>>                                            constant value interval : "
     << constVariableValue.first  << " - " << constVariableValue.second << endl ;
}



void TtDifferentialPlot::bookHisto () {

  // vector with ranges
  vector<double> variableRanges ;

  for ( unsigned int iP = 0 ; iP < theBins.size() ; iP++ ) {
    if ( diffVariableName == "eta" ) {
      // only active bins in the variable on x-axis
      if ( theBins[iP]->getEtaActive() ) {
	variableRanges.push_back ( theBins[iP]->getEtaMin() ) ;
	// also max if last one
	if ( iP == theBins.size()-1 ) variableRanges.push_back ( theBins[iP]->getEtaMax() ) ;
      }
    }
    if ( diffVariableName == "et" ) {
      // only active bins in the variable on x-axis
      if ( theBins[iP]->getEtActive() ) {
	variableRanges.push_back ( theBins[iP]->getEtMin() ) ;
	// also max if last one
	if ( iP == theBins.size()-1 ) variableRanges.push_back ( theBins[iP]->getEtMax() ) ;
      }
    }
  }

  // to book histo with variable binning -> put into array
  int      nBins    = variableRanges.size() - 1 ;
  double * binArray = new double [nBins+1] ;

  for ( int i = 0 ; i < nBins + 1 ; i++ ) {
    binArray[i] = variableRanges[i] ;
  }


  commonName += constVariableName ;
  commonName += "_" ;
  commonName += constVariableValue.first ;
  commonName += "-" ;
  commonName += constVariableValue.second ;
  commonName += "_Vs_" ;
  commonName += diffVariableName ;
  commonName.ReplaceAll ( " " , "" ) ;

  theHisto = new TH1F (commonName , commonName , nBins , binArray ) ;
}


void TtDifferentialPlot::fillHisto (int bin) {
  // loop over bins and find corresponding misid. in the MisIdVs..... histo
  for ( unsigned int iP = 0 ; iP < theBins.size() ; iP++ ) {
//     EffPurFromHistos * currentEffPurFromHistos = theBinPlotters[iP]->getEffPurFromHistos() ;
    //
    bool   isActive   = true ;
    double valueXAxis = -999.99 ;
    // find right bin based on middle of the interval
    if ( diffVariableName == "eta" ) {
      isActive = theBins[iP]->getEtaActive() ;
      valueXAxis = 0.5 * ( theBins[iP]->getEtaMin() + theBins[iP]->getEtaMax() ) ;
    } else if ( diffVariableName == "et"  ) {
      isActive = theBins[iP]->getEtActive() ;
      valueXAxis = 0.5 * ( theBins[iP]->getEtMin() + theBins[iP]->getEtMax() ) ;
    } else {
      throw cms::Exception("Configuration")
	<< "====>>>> TtDifferentialPlot::fillHisto() : illegal diffVariableName = " << diffVariableName << endl;
    }

    // for the moment: ignore inactive bins
    // (maybe later: if a Bin is inactive -> set value to fill well below left edge of histogram to have it in the underflow)

    if ( !isActive ) continue ;

    TH1F* scanHisto = theCollector->getHisto(*theBins[iP]);
    Measurement1D bEff(scanHisto->GetBinContent(bin),scanHisto->GetBinError(bin));
    if (bin == -1) bEff = getBestValue(scanHisto);

    int iBinSet = theHisto->FindBin(valueXAxis) ;
    theHisto->SetBinContent(iBinSet, bEff.value());
    theHisto->SetBinError(iBinSet, bEff.error());
  }
}

Measurement1D TtDifferentialPlot::getBestValue(const TH1F * scanHisto)
{
  int bins = scanHisto->GetNbinsX();
  int binMin = -1;
  float errorMin = 100000, error;
  for(int bin=1; bin<=bins; ++bin){
    error = scanHisto->GetBinError(bin);
    if ((error>0.0) && (error < errorMin)) {
      errorMin = error;
      binMin = bin;
    }
  }
  cout << "binMin " <<binMin << " "<<scanHisto->GetBinError(binMin)<<endl;
  return Measurement1D(scanHisto->GetBinContent(binMin),scanHisto->GetBinError(binMin));
}


void TtDifferentialPlot::eps(const TString & ext)
{
  TStyle *tdrStyle = gROOT->GetStyle("tdrStyle");
  tdrStyle->SetOptFit(0);
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptTitle(0);
//   tdrStyle->SetPadTopMargin(0.01);
//   tdrStyle->SetPadBottomMargin(0.01);
//   tdrStyle->SetPadLeftMargin(0.01);
//   tdrStyle->SetPadRightMargin(0.01);

  TCanvas c2("c2","",600,600);
//   c2.SetTopMargin(0.01);
//   c2.SetBottomMargin(0.01);
//   c2.SetLeftMargin(0.01);
//   c2.SetRightMargin(0.01);
  theHisto->Draw("P");
  cout << theHisto->GetName()+TString(".")+ext<<endl;
  c2.Print(theHisto->GetName()+TString(".")+ext);
}

void TtDifferentialPlot::setLabels(const TString& yLabel)
{
  theHisto->GetXaxis()->SetTitle(diffVariableName);
  theHisto->GetYaxis()->SetTitle(yLabel);
  theHisto->GetYaxis()->SetRangeUser(0.,1.1);
}

void TtDifferentialPlot::setRange(const float min, const float max)
{
  theHisto->GetYaxis()->SetRangeUser(min, max);
}

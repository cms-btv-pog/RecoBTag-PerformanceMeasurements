#ifndef TtDifferentialPlot_H
#define TtDifferentialPlot_H

#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"

#include <vector>

#include "RecoBTag/PerformanceMeasurements/interface/EtaEtBin.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

class TtEtEtaHistoCollector;

class TtDifferentialPlot {
  
 public:

  typedef std::vector <EtaEtBin*> BinVector;
  typedef BinVector::iterator BinVectorIter;
  typedef BinVector::const_iterator BinVectorConstIter;

  enum ConstVarType {constET, constETA };

  TtDifferentialPlot (TtEtEtaHistoCollector * aCollector, ConstVarType constVariable,
	const TString & tagName);

  ~TtDifferentialPlot () ;


  void addBin ( const EtaEtBin * bin ) { theBins.push_back ( bin->clone() ) ; }

  void process (int bin = -1) ;

  void write () ; 
// 
//   void epsPlot(const TString & name);
// 
//   void psPlot(const TString & name);
// 
//   void plot (TCanvas & theCanvas) ;
// 
//   void plot(const TString & name, const TString & ext);

// 
//   void print () const ;

  TH1F * getHisto() {return theHisto;}

  void eps(const TString & ext);
  void setLabels(const TString& yLabel);
  void setRange(const float min, const float max);


 private:
  
  void setVariableName () ;
  
  void bookHisto () ;

  void fillHisto (int bin = -1);
  Measurement1D getBestValue(const TH1F * scanHisto);

  bool processed;

  TtEtEtaHistoCollector * theCollector;

  ConstVarType constVar;
  // the name for the variable with constant value
  TString constVariableName ;
  // the name of the variable to be plotted on the x-axis (e.g. "eta", "pt")
  TString diffVariableName ;

  // value of the constant variable (lower/upper edge of interval)
  std::pair<double,double> constVariableValue ;

  // the common name to describe histograms
  TString commonName ;

  BinVector theBins ;

  TH1F * theHisto    ;

  // the plot Canvas
//   TCanvas * thePlotCanvas ;
  
} ;


#endif

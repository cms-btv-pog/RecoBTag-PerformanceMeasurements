#ifndef TtEtEtaHistoCollector_H
#define TtEtEtaHistoCollector_H

#include "AnalysisDataFormats/TopObjects/interface/TopJet.h" 
#include "RecoBTag/PerformanceMeasurements/interface/EtaEtBin.h"
#include "RecoBTag/PerformanceMeasurements/interface/TtDifferentialPlot.h"

#include <string>
#include <vector>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TSystem.h>

/** \class TtEtEtaHistoCollector
 *
 *  Top level steering routine for b tag performance analysis.
 *
 */

class TtEtEtaHistoCollector {
public:

  typedef std::map <EtaEtBin* , TH1F*> BinHistoMap;
  typedef BinHistoMap::iterator BinHistoMapIter;
  typedef BinHistoMap::const_iterator BinHistoMapConstIter;
  typedef std::vector<TtDifferentialPlot*> diffPlotVect;

  TtEtEtaHistoCollector(){}
  TtEtEtaHistoCollector(const TString & baseName, const std::vector<double> & eTRanges, 
    const std::vector<double> & etaRanges, const double lowerBound,
    const double upperBound, const int bins, const bool update);
//   TtEtEtaHistoCollector(const TtEtEtaHistoCollector & other);

  ~TtEtEtaHistoCollector();
  
//   TtEtEtaHistoCollector & operator= (TtEtEtaHistoCollector & other);

  /**
   *   Calls Sumw2 for each histo. The error per bin will be computed as 
   *   sqrt(sum of squares of weight) for each bin.
   */
  void Sumw2();
  
  void analyze(const TopJet & jet, double lhr, const double weight);

  void write();

  TString baseName() const {return baseName_;}
  void baseName(const TString & newName) {baseName_= newName;}

  void setHisto(const EtaEtBin& bin, TH1F * histo);
  TH1F * getHisto(const EtaEtBin& bin);

  BinHistoMapIter histoBegin() {return binHistoMap.begin();}
  BinHistoMapIter histoEnd() {return binHistoMap.end();}

  BinHistoMapConstIter histoBegin() const {return binHistoMap.begin();}
  BinHistoMapConstIter histoEnd() const {return binHistoMap.end();}

//   void integratedHisto(const EtaEtBin& bin, TH1F *histo);
  void integratedHisto(const TH1F * const histo, TH1F * newHisto);
//   TH1F* integratedHisto(const TH1F * const histo, TString newName);

  TtEtEtaHistoCollector * buildIntegratedCollector();

  TtEtEtaHistoCollector * buildRatioCollector(const TString& name,
			TtEtEtaHistoCollector & otherColl);
  TtEtEtaHistoCollector * buildDiffCollector(const TString& name,
			TtEtEtaHistoCollector & otherColl);

  void buildDifferentialPlots(int bin = -1);

  void eps( TH1F * histo, const TString & ext);
  void eps(const TString & ext);
  void setLabels(const TString& xLabel, const TString& yLabel);
  void setRange(const float min, const float max);

  const diffPlotVect & getDifferentialPlots() const { return differentialPlots;}


   private:

  // Contains plots for each bin of rapidity and et.
  diffPlotVect differentialPlots;
  bool diffPlotOK;

  std::vector<double> etRanges_, etaRanges_;
  double lowerBound_, upperBound_;
  int bins_;

  BinHistoMap binHistoMap;
  //std::vector<BinHisto*>   binPlotters;
  TString baseName_;



};


#endif

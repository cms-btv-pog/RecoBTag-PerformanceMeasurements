#ifndef TtBTagAnalysis_H
#define TtBTagAnalysis_H


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h" 
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h" 
#include "RecoBTag/PerformanceMeasurements/interface/TtEtEtaHistoCollector.h" 

#include "TopQuarkAnalysis/TopEventSelection/interface/TtDilepLRSignalSelObservables.h"


#include "TopQuarkAnalysis/TopTools/interface/LRHelpFunctions.h"

// #include "TtBTagAnalysisHisto.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility> 


class TtBTagAnalysis : public edm::EDAnalyzer {

  public:

    /// Constructor
    TtBTagAnalysis(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~TtBTagAnalysis();
    virtual void endJob();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);


  private: 

    void writePlotter(TtEtEtaHistoCollector* plotter);
    void analyzeJet(const TopJet& jet, double lr, bool matchB);

    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug, mcMode;

    std::string rootFileName, obsFileName;
    double weight;
    

    edm::InputTag evtsols;
    LRHelpFunctions *myLRhelper;
    std::vector<int> obsNrs;
    int lrBins;
    double lrMin,lrMax;
    const char* lrFits;
    vector<double> etaRanges, etRanges;

    //observable histogram variables
    int nrSignalSelObs;

    TtDilepLRSignalSelObservables * myLRSignalSelObservables;

  TtEtEtaHistoCollector * allJetsHistos, *lrSHistos, *lrBHistos, *allLightHistos;
  std::vector< TtEtEtaHistoCollector *> taggedJetsHistos;
  std::vector< TtEtEtaHistoCollector *> taggedLightJetsHistos, taggedBJetsHistos;
  typedef std::pair<std::string, double> StringIntPair;
  typedef vector< StringIntPair > StringIntPairVect;
  StringIntPairVect bTagger;
//   TopJetSelection jetSelector;
};


#endif

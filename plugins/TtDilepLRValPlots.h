#ifndef TtDilepLRValPlots_H
#define TtDilepLRValPlots_H


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h" 
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h" 

#include "TopQuarkAnalysis/TopEventSelection/interface/TtDilepLRSignalSelObservables.h"


#include "TopQuarkAnalysis/TopTools/interface/LRHelpFunctions.h"

// #include "TtDilepLRValPlotsHisto.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility> 


class TtDilepLRValPlots : public edm::EDAnalyzer {

  public:

    /// Constructor
    TtDilepLRValPlots(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~TtDilepLRValPlots();
    virtual void endJob();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);


  private: 


    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug;

    std::string rootFileName, obsFileName;
    bool csa;
    double weight;
    

    edm::InputTag evtsols;
    LRHelpFunctions *myLRhelper;
    std::vector<int> obsNrs;
    int lrBins;
    double lrMin,lrMax;
    const char* lrFits;

    //observable histogram variables
    int nrSignalSelObs;

    TtDilepLRSignalSelObservables * myLRSignalSelObservables;

};


#endif

#ifndef TtDilepLRObsPlots_H
#define TtDilepLRObsPlots_H


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "AnalysisDataFormats/TopObjects/interface/TopLepton.h" 
#include "AnalysisDataFormats/TopObjects/interface/TopObject.h" 
#include "AnalysisDataFormats/TopObjects/interface/TopParticle.h" 
#include "AnalysisDataFormats/TopObjects/interface/TopMET.h" 
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h" 
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h" 
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h" 
#include "TopQuarkAnalysis/TopEventSelection/interface/TtDilepLRSignalSelObservables.h"


#include "TopQuarkAnalysis/TopTools/interface/LRHelpFunctions.h"

// #include "TtDilepLRObsPlotsHisto.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility> 


class TtDilepLRObsPlots : public edm::EDAnalyzer {

  public:

    /// Constructor
    TtDilepLRObsPlots(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~TtDilepLRObsPlots();
    virtual void endJob();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);


  private: 

    bool allMatch(const TtDilepEvtSolution & sol) ;

    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug;

    std::string rootFileName;
    std::string leptonFlavour;
    double weight;
    int signal, background, goodSolution, allSolution, B, nonB, tau;
    int bestSol;
    

    edm::InputTag evtsols, jetSource_;
    LRHelpFunctions *myLRhelper;
    std::vector<int> obsNrs;
    std::vector<double> obsMin,obsMax;
    std::vector<const char*> obsFits;
    std::vector<std::string> obsFitsStr;

    int nrSignalSelObs, nrSignalSelHistBins;
    TtDilepLRSignalSelObservables * myLRSignalSelObservables;


};


#endif

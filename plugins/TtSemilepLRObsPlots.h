#ifndef TtSemilepLRObsPlots_H
#define TtSemilepLRObsPlots_H


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
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h" 
#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiLRJetCombObservables.h"


#include "TopQuarkAnalysis/TopTools/interface/LRHelpFunctions.h"

// #include "TtSemilepLRObsPlotsHisto.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility> 


class TtSemilepLRObsPlots : public edm::EDAnalyzer {

  public:

    /// Constructor
    TtSemilepLRObsPlots(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~TtSemilepLRObsPlots();
    virtual void endJob();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);


  private: 

    bool allMatch(const TtSemiEvtSolution & sol) ;

    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug;

    std::string rootFileName;
    std::string leptonFlavour; 
    bool csa;
    double weight;
    float signal, background, allSolution, B, nonB, tau, goodSolution, goodSolution1, goodSolution2;
    int bestSol;
    std::string bTagCutLabel;
    double bCut;

    edm::InputTag evtsols, jetSource_;
    LRHelpFunctions *myLRhelper;
    std::vector<int> obsNrs;
    std::vector<double> obsMin,obsMax;
    std::vector<const char*> obsFits;
    std::vector<std::string> obsFitsStr;

    int nrJetCombObs, nrJetCombHistBins;
    TtSemiLRJetCombObservables * myLRJetCombObservables;


};


#endif

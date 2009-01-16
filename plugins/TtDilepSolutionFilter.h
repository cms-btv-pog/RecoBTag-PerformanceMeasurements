#ifndef TtDilepSolutionFilter_H
#define TtDilepSolutionFilter_H


#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h" 
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h" 

#include <vector>
#include <map>
#include <string>
#include <utility> 
#include "TFile.h"
#include "TH1F.h"


class TtDilepSolutionFilter : public edm::EDFilter {

  public:

    /// Constructor
    TtDilepSolutionFilter(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~TtDilepSolutionFilter();
    virtual void endJob();

    /// Perform the real analysis
    virtual bool filter(edm::Event & iEvent, edm::EventSetup const& iSetup);


  private: 

    bool checkElectron (const pat::Electron & electron) const;

    // Switch for debug output
    bool debug;
    TString rootFileName;

    // The file which will store the histos, if any
    TFile *theFile;
    TH1F *procA, *procB, *procBee, *procBem, *procBmm;
    double jetPtCut, leptonPtCut, leptonTriggerPtCut, weight, jetEtCut;
    double signal, background, goodSolution, allSolution, B, nonB, tau;
    int bestSol,exception;
    edm::InputTag evtsols;
};


#endif

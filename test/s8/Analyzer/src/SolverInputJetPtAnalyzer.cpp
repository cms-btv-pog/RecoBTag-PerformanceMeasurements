/**
 * SolverInputJetPtAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include <TDirectory.h>
#include <TH1F.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8TreeInfo.h"
#include "S8Tree/interface/S8TriggerCenter.h"
#include "Utility/interface/LeptonInJetDiscriminator.h"
#include "Utility/interface/TwoMuonsInJet.h"
#include "Utility/interface/SolverLeptonInJetPlots.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/SolverInputJetPtAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;

using s8::SolverInputJetPtAnalyzer;

SolverInputJetPtAnalyzer::SolverInputJetPtAnalyzer() throw()
{
}

SolverInputJetPtAnalyzer::~SolverInputJetPtAnalyzer() throw()
{
}



// Analyzer interface
//
void SolverInputJetPtAnalyzer::init()
{
    _muonInJetTagger.reset(new TaggerOperatingPoint());
    _awayJetTagger.reset(new TaggerOperatingPoint());

    _muonInJet.reset(new TwoMuonsInJet());
    _muonInJet->setDelegate(this);
    _muonInJet->setAwayJetTaggerOperatingPoint(_awayJetTagger.get());

    _leptonInJetDiscriminator.reset(new LeptonInJetPtDiscriminator());

    _solverNPlots.reset(new SolverLeptonInJetPlots("n"));
    _solverPPlots.reset(new SolverLeptonInJetPlots("p"));

    _solverNPlots->setDiscriminator(_leptonInJetDiscriminator.get());
    _solverPPlots->setDiscriminator(_leptonInJetDiscriminator.get());

    _solverNPlots->setTaggerOperatingPoint(_muonInJetTagger.get());
    _solverPPlots->setTaggerOperatingPoint(_muonInJetTagger.get());

    _solverNPlots->init();
    _solverPPlots->init();
}

void SolverInputJetPtAnalyzer::eventDidLoad(const Event *event)
{
    // call muon in jet
    //
    (*_muonInJet)(event);
}

void SolverInputJetPtAnalyzer::print(std::ostream &) const
{
}

void SolverInputJetPtAnalyzer::save(TDirectory *directory) const
{
    // Save (n) and (p) plots
    //
    _solverNPlots->save(directory);
    _solverPPlots->save(directory);
}



// SolverInputOptionsDelegate interface
//
void SolverInputJetPtAnalyzer::optionDataIsSet(const bool &value)
{
    if (value)
    {
        _solverNPlots->setIsFlavoured(false);
        _solverPPlots->setIsFlavoured(false);
    }
}

// MuonInJetOptionsDelegate interface
//
void SolverInputJetPtAnalyzer::optionTagIsSet(const std::string &tag)
{
    *_muonInJetTagger << tag;
}

void SolverInputJetPtAnalyzer::optionAwayTagIsSet(const std::string &tag)
{
    *_awayJetTagger << tag;
}

// MuonInJetDelegate interface
//
void SolverInputJetPtAnalyzer::muonIsInJetPlusAwayJet(const Lepton *,
                                                      const Lepton *lepton,
                                                      const Jet *jet)
{
    _solverNPlots->fill(lepton, jet);
}

void SolverInputJetPtAnalyzer::muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                            const Lepton *lepton,
                                                            const Jet *jet)
{
    _solverPPlots->fill(lepton, jet);
}

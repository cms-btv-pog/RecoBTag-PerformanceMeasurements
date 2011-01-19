/**
 * SolverInputAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <stdexcept>

#include <TDirectory.h>
#include <TH1F.h>
#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8Fwd.h"
#include "S8Tree/interface/S8GenEvent.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8TreeInfo.h"
#include "S8Tree/interface/S8TriggerCenter.h"
#include "Utility/interface/LeptonInJetDiscriminator.h"
#include "Utility/interface/MuonInJet.h"
#include "Utility/interface/SolverLeptonInJetPlots.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/SolverInputAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;

using s8::SolverInputAnalyzer;

SolverInputAnalyzer::SolverInputAnalyzer() throw():
    _gluonSplitting(KEEP),
    _minPtHat(0),
    _maxPtHat(0),
    _otherEvent(false)
{
}

SolverInputAnalyzer::~SolverInputAnalyzer() throw()
{
}



// Analyzer interface
//
void SolverInputAnalyzer::init()
{
    _gluonSplitting = KEEP;
    _minPtHat = 0;
    _maxPtHat = 0;
    _otherEvent = false;

    _muonInJetTagger.reset(new TaggerOperatingPoint());
    _awayJetTagger.reset(new TaggerOperatingPoint());

    _muonInJet.reset(new MuonInJet());
    _muonInJet->setDelegate(this);
    _muonInJet->setAwayJetTaggerOperatingPoint(_awayJetTagger.get());

    _leptonInJetDiscriminator.reset(new LeptonInJetPtRelDiscriminator());

    _solverNPlots.reset(new SolverLeptonInJetPlots("n"));
    _solverPPlots.reset(new SolverLeptonInJetPlots("p"));

    _solverNPlots->setDiscriminator(_leptonInJetDiscriminator.get());
    _solverPPlots->setDiscriminator(_leptonInJetDiscriminator.get());

    _solverNPlots->setTaggerOperatingPoint(_muonInJetTagger.get());
    _solverPPlots->setTaggerOperatingPoint(_muonInJetTagger.get());

    _solverNPlots->init();
    _solverPPlots->init();
}

void SolverInputAnalyzer::eventDidLoad(const Event *event)
{
    using s8::Jets;

    // PtHat cut
    //
    if (_minPtHat &&
        event->gen()->ptHat() < _minPtHat ||
        
        _maxPtHat &&
        event->gen()->ptHat() > _maxPtHat)

        return;

    // Gluon Splitting
    //
    switch(_gluonSplitting)
    {
        case KEEP: break;

        case REMOVE:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB))
                    return;

                break;
            }

        case ONLY:
            {
                if (!event->gen()->isGluonSplitting(s8::GenEvent::BB))
                    return;

                break;
            }

        case ENHANCE:
            {
                if (!event->gen()->isGluonSplitting(s8::GenEvent::BB))
                {
                    _otherEvent = !_otherEvent;

                    if (_otherEvent)
                        return;
                }

                break;
            }

        case SIMULATE:
            {
                const Jets *jets = event->jets();

                if (2 != jets->size())
                    return;

                for(Jets::const_iterator jet = jets->begin();
                    jets->end() != jet;
                    ++jet)
                {
                    if (30 >= (*jet)->p4()->Pt())
                        return;
                }

                const Jet *jet1 = (*jets)[0];
                const Jet *jet2 = (*jets)[1];

                if (1.5 <= jet1->p4()->DeltaR(*(jet2->p4())))
                    return;
            }
    }


    // call muon in jet
    //
    (*_muonInJet)(event);
}

void SolverInputAnalyzer::save(TDirectory *directory) const
{
    // Save (n) and (p) plots
    //
    _solverNPlots->save(directory);
    _solverPPlots->save(directory);
}



// SolverInputOptionsDelegate interface
//
void SolverInputAnalyzer::optionDataIsSet(const bool &value)
{
    if (value)
    {
        _solverNPlots->setIsFlavoured(false);
        _solverPPlots->setIsFlavoured(false);
    }
}

void SolverInputAnalyzer::optionGluonSplittingIsSet(const GluonSplitting &value)
{
    _gluonSplitting = value;
}

void SolverInputAnalyzer::optionPtHatIsSet(const int &min, const int &max)
{
    using std::runtime_error;

    if (0 > min || 0 > max)
        throw runtime_error("Negative value of the pThat cut range is supplied");

    _minPtHat = min;
    _maxPtHat = max;
}

// MuonInJetOptionsDelegate interface
//
void SolverInputAnalyzer::optionTagIsSet(const std::string &tag)
{
    *_muonInJetTagger << tag;
}

void SolverInputAnalyzer::optionAwayTagIsSet(const std::string &tag)
{
    *_awayJetTagger << tag;
}

void SolverInputAnalyzer::optionMuonPtIsSet(const double &pt)
{
    _muonInJet->setMuonMinimumPtCut(pt);
}

// MuonInJetDelegate interface
//
void SolverInputAnalyzer::muonIsInJetPlusAwayJet(const Lepton *lepton,
                                                 const Jet *jet)
{
    _solverNPlots->fill(lepton, jet);
}

void SolverInputAnalyzer::muonIsInJetPlusTaggedAwayJet(const Lepton *lepton,
                                                       const Jet *jet)
{
    _solverPPlots->fill(lepton, jet);
}

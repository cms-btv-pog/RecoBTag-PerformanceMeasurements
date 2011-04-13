/**
 * SolverInputAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <iostream>
#include <stdexcept>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8Fwd.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Utility/interface/LeptonInJetDiscriminator.h"
#include "Utility/interface/MuonInJet.h"
#include "Selector/interface/S8Selector.h"
#include "Utility/interface/SolverLeptonInJetPlots.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/SolverInputAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;
using std::runtime_error;

using s8::SolverInputAnalyzer;

SolverInputAnalyzer::SolverInputAnalyzer() throw()
{
    _s8_selector = 0;
    _reweights = 0;
}

SolverInputAnalyzer::~SolverInputAnalyzer() throw()
{
    if (_s8_selector)
        delete _s8_selector;
}



// Analyzer interface
//
void SolverInputAnalyzer::init()
{
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

    _s8_selector = new S8Selector();
    _n_primary_vertices = 0;
}

void SolverInputAnalyzer::treeDidLoad(const TreeInfo *,
                                      const TriggerCenter *trigger_center)
{
    _s8_selector->treeDidLoad(trigger_center);
}


void SolverInputAnalyzer::eventDidLoad(const Event *event)
{
    const Event *modified_event = (*_s8_selector)(event);
    if (!modified_event)
        return;

    event = modified_event;
    _n_primary_vertices = event->primaryVertices()->size();

    // call muon in jet
    //
    (*_muonInJet)(event);
}

void SolverInputAnalyzer::print(std::ostream &out) const
{
    _s8_selector->print(out);
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

void SolverInputAnalyzer::optionMuonPtIsSet(const Range &value)
{
    _n_muon_pt = value;
}

void SolverInputAnalyzer::optionJetPtIsSet(const Range &value)
{
    _n_jet_pt = value;
}

void SolverInputAnalyzer::optionJetEtaIsSet(const Range &value)
{
    _n_jet_eta = value;
}

// PythiaOptionsDelegate interface
//
void SolverInputAnalyzer::optionGluonSplittingIsSet(const GluonSplitting &value)
{
    _s8_selector->optionGluonSplittingIsSet(value);
}

void SolverInputAnalyzer::optionPtHatIsSet(const Range &value)
{
    _s8_selector->optionPtHatIsSet(value);
}

// Trigger options
//
void SolverInputAnalyzer::optionTriggerIsSet(const Trigger &trigger)
{
    _s8_selector->optionTriggerIsSet(trigger);
}

void SolverInputAnalyzer::optionSimulateTriggerIsSet(const bool &value)
{
    _s8_selector->optionSimulateTriggerIsSet(value);
}

void SolverInputAnalyzer::optionReweightTriggerIsSet(const std::string &filename)
{
    _reweight_file.reset(new TFile(filename.c_str(), "read"));
    if (!_reweight_file->IsOpen())
        throw runtime_error("Failed to open reweights file: " + filename);

    _reweights = dynamic_cast<TH3 *>(_reweight_file->Get("n_pt_eta_pv"));
    if (!_reweights)
        throw runtime_error("Failed to extract reweights from: " + filename);
}

// Misc Options
//
void SolverInputAnalyzer::optionPrimaryVerticesIsSet(const Range &primary_vertices)
{
    _s8_selector->optionPrimaryVerticesIsSet(primary_vertices);
}

// MuonInJetDelegate interface
//
bool SolverInputAnalyzer::shouldSkipMuonInJetPlusAwayJet(const Lepton *muon,
                                                         const Jet *jet)
{
    return !isValueInRange(muon->p4()->Pt(), _n_muon_pt) ||
           !isValueInRange(jet->p4()->Pt(), _n_jet_pt) ||
           !isValueInRange(fabs(jet->p4()->Eta()), _n_jet_eta);
}

void SolverInputAnalyzer::muonIsInJetPlusAwayJet(const Lepton *lepton,
                                                 const Jet *jet)
{
    double weight = 1;
    const double jet_pt = jet->p4()->Pt();

    if (_reweights &&
        isValueInRange(jet_pt, _n_jet_pt))

        weight = _reweights->GetBinContent(_reweights->FindBin(jet->p4()->Pt(),
                                                               jet->p4()->Eta(),
                                                               _n_primary_vertices));

    _solverNPlots->fill(lepton, jet, weight);
}

void SolverInputAnalyzer::muonIsInJetPlusTaggedAwayJet(const Lepton *lepton,
                                                       const Jet *jet)
{
    double weight = 1;
    const double jet_pt = jet->p4()->Pt();

    if (_reweights &&
        isValueInRange(jet_pt, _n_jet_pt))

        weight = _reweights->GetBinContent(_reweights->FindBin(jet->p4()->Pt(),
                                                               jet->p4()->Eta(),
                                                               _n_primary_vertices));

    _solverPPlots->fill(lepton, jet, weight);
}

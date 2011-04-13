/**
 * MonitorTriggerAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Feb 25, 2011
 * Copyright 2011, All rights reserved
 */

#include <cmath>

#include <iostream>
#include <stdexcept>

#include <TDirectory.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8GenEvent.h"
#include "S8Tree/interface/S8Fwd.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Selector/interface/S8Selector.h"
#include "Utility/interface/MonitorPlots.h"
#include "Utility/interface/MuonInJet.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/MonitorTriggerAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;

using s8::MonitorTriggerAnalyzer;

MonitorTriggerAnalyzer::MonitorTriggerAnalyzer() throw():
    _use_trigger_prescale(false),
    _prescale(0)
{
    _s8_selector = 0;
}

MonitorTriggerAnalyzer::~MonitorTriggerAnalyzer() throw()
{
    if (_s8_selector)
        delete _s8_selector;
}



// Analyzer interface
//
void MonitorTriggerAnalyzer::init()
{
    _muonInJetTagger.reset(new TaggerOperatingPoint());
    _awayJetTagger.reset(new TaggerOperatingPoint());

    _muonInJet.reset(new MuonInJet());
    _muonInJet->setDelegate(this);
    _muonInJet->setAwayJetTaggerOperatingPoint(_awayJetTagger.get());

    _muons.reset(new TH1D("muons", "N muons",
                            10, 0, 10));
    // Force Manual memory management
    //
    _muons->SetDirectory(0);
    _muons->GetXaxis()->SetTitle("N_{#mu}");

    _pthat.reset(new TH1D("pthat", "pT hat", 100, 0, 200));
    // Force Manual memory management
    //
    _pthat->SetDirectory(0);
    _pthat->GetXaxis()->SetTitle("#hat{p_{T}}");

    _jets.reset(new TH1D("jets", "N jets",
                            10, 0, 10));
    // Force Manual memory management
    //
    _jets->SetDirectory(0);
    _jets->GetXaxis()->SetTitle("N_{jet}");

    // Generic
    //
    _monitorMuons.reset(new MonitorLepton("allmus"));
    _monitorJets.reset(new MonitorJet("alljets"));
    _monitorJets->plot(MonitorBase::PT)->SetBins(230, 0, 230); 
    _monitorLeadingJet.reset(new MonitorJet("leadingjet"));

    // (n)
    //
    _monitorNMuons.reset(new MonitorLepton("nmu"));
    _monitorNJets.reset(new MonitorJet("njet"));
    _monitorNJets->plot(MonitorBase::PT)->SetBins(230, 0, 230); 
    _monitorNDelta.reset(new MonitorDelta("nmujet"));

    // (p)
    //
    _monitorPMuons.reset(new MonitorLepton("pmu"));
    _monitorPJets.reset(new MonitorJet("pjet"));
    _monitorPJets->plot(MonitorBase::PT)->SetBins(230, 0, 230); 
    _monitorPDelta.reset(new MonitorDelta("pmujet"));

    // change binning in the Delta plots
    //
    MonitorDelta *deltas[] = { _monitorNDelta.get(), _monitorPDelta.get() };
    for(int i = 0; 2 > i; ++i)
    {
        deltas[i]->plot(MonitorDelta::PTREL)->SetBins(70, 0, 7);
    }

    _s8_selector = new S8Selector();
}

void MonitorTriggerAnalyzer::treeDidLoad(const TreeInfo *,
                                  const TriggerCenter *trigger_center)
{
    _s8_selector->treeDidLoad(trigger_center);
}

void MonitorTriggerAnalyzer::eventDidLoad(const Event *event)
{
    using s8::Jets;
    using s8::Leptons;

    // check if event passes the selection
    //
    const Event *modified_event = (*_s8_selector)(event);
    if (!modified_event)
        return;

    event = modified_event;

    if (_use_trigger_prescale)
    {
        // Extract prescale of the trigger
        //
        // Search for user defined triggers among those that passed the event
        //
        s8::Triggers::const_iterator found_trigger =
            find_if(event->triggers()->begin(),
                    event->triggers()->end(),
                    _trigger_predicator);

        // Check if trigger was found and passed the event
        //
        if (event->triggers()->end() != found_trigger)
            _prescale = (*found_trigger)->prescale();
    }

    _pthat->Fill(event->gen()->ptHat());
    _muons->Fill(event->muons()->size());
    _jets->Fill(event->jets()->size());

    for(Leptons::const_iterator muon = event->muons()->begin();
        event->muons()->end() != muon;
        ++muon)

        _monitorMuons->fill(_prescale ? _prescale : 1, *muon);

    const Jet *leadingJet = 0;
    for(Jets::const_iterator jet = event->jets()->begin();
        event->jets()->end() != jet;
        ++jet)
    {
        _monitorJets->fill(_prescale ? _prescale : 1, *jet);

        if (leadingJet &&
            leadingJet->p4()->Pt() >= (*jet)->p4()->Pt())

            continue;

        leadingJet = *jet;
    }

    if (leadingJet)
        _monitorLeadingJet->fill(_prescale ? _prescale : 1, leadingJet);

    (*_muonInJet)(event);

    _prescale = 0;
}

void MonitorTriggerAnalyzer::print(std::ostream &out) const
{
    _s8_selector->print(out);
}

void MonitorTriggerAnalyzer::save(TDirectory *directory) const
{
    TDirectory *output = directory->mkdir("MonitorTriggerAnalyzer");
    if (!output)
    {
        cerr << "MonitorTriggerAnalyzer: Failed to create output folder: "
            << "no output is saved" << endl;

        return;
    }

    output->cd();

    saveGenericPlots(output);
    saveNPlots(output);
    savePPlots(output);
}

// MuonInJetOptionsDelegate interface
//
void MonitorTriggerAnalyzer::optionTagIsSet(const std::string &tag)
{
    *_muonInJetTagger << tag;
}

void MonitorTriggerAnalyzer::optionAwayTagIsSet(const std::string &tag)
{
    *_awayJetTagger << tag;
}

void MonitorTriggerAnalyzer::optionMuonPtIsSet(const Range &value)
{
    _n_muon_pt = value;
}

void MonitorTriggerAnalyzer::optionJetPtIsSet(const Range &value)
{
    _n_jet_pt = value;
}

void MonitorTriggerAnalyzer::optionJetEtaIsSet(const Range &value)
{
    _n_jet_eta = value;
}

// MuonInJetDelegate interface
//
bool MonitorTriggerAnalyzer::shouldSkipMuonInJetPlusAwayJet(const Lepton *muon,
                                                     const Jet *jet)
{
    return !isValueInRange(muon->p4()->Pt(), _n_muon_pt) ||
           !isValueInRange(jet->p4()->Pt(), _n_jet_pt) ||
           !isValueInRange(fabs(jet->p4()->Eta()), _n_jet_eta);
}

void MonitorTriggerAnalyzer::muonIsInJetPlusAwayJet(const Lepton *muon,
                                             const Jet *jet)
{
    _monitorNMuons->fill(_prescale ? _prescale : 1, muon);
    _monitorNJets->fill(_prescale ? _prescale : 1, jet);
    _monitorNDelta->fill(_prescale ? _prescale : 1, muon->p4(), jet->p4());
}

void MonitorTriggerAnalyzer::muonIsInJetPlusTaggedAwayJet(const Lepton *muon,
                                                   const Jet *jet)
{
    _monitorPMuons->fill(_prescale ? _prescale : 1, muon);
    _monitorPJets->fill(_prescale ? _prescale : 1, jet);
    _monitorPDelta->fill(_prescale ? _prescale : 1, muon->p4(), jet->p4());
}

// PythiaOptionsDelegate interface
//
void MonitorTriggerAnalyzer::optionGluonSplittingIsSet(const GluonSplitting &value)
{
    _s8_selector->optionGluonSplittingIsSet(value);
}

void MonitorTriggerAnalyzer::optionPtHatIsSet(const Range &value)
{
    _s8_selector->optionPtHatIsSet(value);
}

// Trigger options
//
void MonitorTriggerAnalyzer::optionTriggerIsSet(const Trigger &trigger)
{
    _trigger = trigger;
    _trigger_predicator.setSearchTrigger(_trigger);

    _s8_selector->optionTriggerIsSet(trigger);
}

void MonitorTriggerAnalyzer::optionUseTriggerPrescaleIsSet(const bool &value)
{
    _use_trigger_prescale = value;
}



void MonitorTriggerAnalyzer::saveGenericPlots(TDirectory *directory) const
{
    TDirectory *subdir = directory->mkdir("generic");
    if (!subdir)
    {
        cerr << "Failed to create 'generic' TDirectory. Plots are not saved"
            << endl;

        return;
    }

    subdir->cd();

    _muons->Write();
    _jets->Write();
    _pthat->Write();

    // Generic
    //
    _monitorMuons->save(subdir);
    _monitorJets->save(subdir);
    _monitorLeadingJet->save(subdir);

    directory->cd();
}

void MonitorTriggerAnalyzer::saveNPlots(TDirectory *directory) const
{
    TDirectory *subdir = directory->mkdir("n");
    if (!subdir)
    {
        cerr << "Failed to create 'n' TDirectory. Plots are not saved"
            << endl;

        return;
    }

    subdir->cd();

    // (n)
    //
    _monitorNMuons->save(subdir);
    _monitorNJets->save(subdir);
    _monitorNDelta->save(subdir);

    directory->cd();
}

void MonitorTriggerAnalyzer::savePPlots(TDirectory *directory) const
{
    TDirectory *subdir = directory->mkdir("p");
    if (!subdir)
    {
        cerr << "Failed to create 'p' TDirectory. Plots are not saved"
            << endl;

        return;
    }

    subdir->cd();

    // (p)
    //
    _monitorPMuons->save(subdir);
    _monitorPJets->save(subdir);
    _monitorPDelta->save(subdir);

    directory->cd();
}

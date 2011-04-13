/**
 * MonitorAnalyzer
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
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8GenEvent.h"
#include "S8Tree/interface/S8Fwd.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "S8Tree/interface/S8PrimaryVertex.h"
#include "Selector/interface/S8Selector.h"
#include "Utility/interface/MonitorPlots.h"
#include "Utility/interface/MuonInJet.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/MonitorAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;
using std::runtime_error;

using s8::MonitorAnalyzer;

MonitorAnalyzer::MonitorAnalyzer() throw()
{
    _s8_selector = 0;
    _reweights = 0;
    _event = 0;
}

MonitorAnalyzer::~MonitorAnalyzer() throw()
{
    if (_s8_selector)
        delete _s8_selector;
}



// Analyzer interface
//
void MonitorAnalyzer::init()
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
    _muons->Sumw2();

    _pthat.reset(new TH1D("pthat", "pT hat", 100, 0, 200));
    // Force Manual memory management
    //
    _pthat->SetDirectory(0);
    _pthat->GetXaxis()->SetTitle("#hat{p_{T}}");
    _pthat->Sumw2();

    _jets.reset(new TH1D("jets", "N jets",
                            10, 0, 10));
    // Force Manual memory management
    //
    _jets->SetDirectory(0);
    _jets->GetXaxis()->SetTitle("N_{jet}");
    _jets->Sumw2();

    _primary_vertices.reset(new TH1D("pvs", "Primary Vertices",
                                     10, 0, 10));
    // Forse Manual memory management
    //
    _primary_vertices->SetDirectory(0);
    _primary_vertices->GetXaxis()->SetTitle("N_{PV}");
    _primary_vertices->Sumw2();

    // Generic
    //
    _monitorMuons.reset(new MonitorLepton("allmus"));
    _monitorJets.reset(new MonitorJet("alljets"));
    _monitorJets->plot(MonitorBase::PT)->SetBins(230, 0, 230); 
    _monitorLeadingJet.reset(new MonitorJet("leadingjet"));

    // (n)
    //
    _monitorNMuonInJet.reset(new MonitorMuonInJet("n"));
    _n_primary_vertices.reset(new TH1D("npvs", "(n) Primary Vertices",
                                        10, 0, 10));
    // Forse Manual memory management
    //
    _n_primary_vertices->SetDirectory(0);
    _n_primary_vertices->GetXaxis()->SetTitle("(n) N_{PV}");
    _n_primary_vertices->Sumw2();

    _n_primary_vertices_z.reset(new TH1D("npvz", "(n) Primary Vertices Z",
                                        400, -20, 20));
    // Forse Manual memory management
    //
    _n_primary_vertices_z->SetDirectory(0);
    _n_primary_vertices_z->GetXaxis()->SetTitle("cm");
    _n_primary_vertices_z->Sumw2();

    // (p)
    //
    _monitorPMuonInJet.reset(new MonitorMuonInJet("p"));

    // change binning in the Delta plots
    //
    MonitorDelta *deltas[] = { _monitorNMuonInJet->delta(),
                               _monitorPMuonInJet->delta() };
    for(int i = 0; 2 > i; ++i)
    {
        deltas[i]->plot(MonitorDelta::PTREL)->SetBins(70, 0, 7);
    }

    _s8_selector = new S8Selector();
    _are_n_plots_filled = false;
}

void MonitorAnalyzer::treeDidLoad(const TreeInfo *,
                                  const TriggerCenter *trigger_center)
{
    _s8_selector->treeDidLoad(trigger_center);
}

void MonitorAnalyzer::eventDidLoad(const Event *event)
{
    using s8::Jets;
    using s8::Leptons;

    // check if event passes the selection
    //
    const Event *modified_event = (*_s8_selector)(event);
    if (!modified_event)
        return;

    event = modified_event;

    _pthat->Fill(event->gen()->ptHat());
    _muons->Fill(event->muons()->size());
    _jets->Fill(event->jets()->size());
    _primary_vertices->Fill(event->primaryVertices()->size());

    for(Leptons::const_iterator muon = event->muons()->begin();
        event->muons()->end() != muon;
        ++muon)

        _monitorMuons->fill(1, *muon);

    const Jet *leadingJet = 0;
    for(Jets::const_iterator jet = event->jets()->begin();
        event->jets()->end() != jet;
        ++jet)
    {
        _monitorJets->fill(1, *jet);

        if (leadingJet &&
            leadingJet->p4()->Pt() >= (*jet)->p4()->Pt())

            continue;

        leadingJet = *jet;
    }

    if (leadingJet)
        _monitorLeadingJet->fill(1, leadingJet);

    _are_n_plots_filled = false;
    _event = event;

    (*_muonInJet)(event);

    _event = 0;
}

void MonitorAnalyzer::print(std::ostream &out) const
{
    _s8_selector->print(out);
}

void MonitorAnalyzer::save(TDirectory *directory) const
{
    TDirectory *output = directory->mkdir("MonitorAnalyzer");
    if (!output)
    {
        cerr << "MonitorAnalyzer: Failed to create output folder: "
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
void MonitorAnalyzer::optionTagIsSet(const std::string &tag)
{
    *_muonInJetTagger << tag;
}

void MonitorAnalyzer::optionAwayTagIsSet(const std::string &tag)
{
    *_awayJetTagger << tag;
}

void MonitorAnalyzer::optionMuonPtIsSet(const Range &value)
{
    _n_muon_pt = value;
}

void MonitorAnalyzer::optionJetPtIsSet(const Range &value)
{
    _n_jet_pt = value;
}

void MonitorAnalyzer::optionJetEtaIsSet(const Range &value)
{
    _n_jet_eta = value;
}

// MuonInJetDelegate interface
//
bool MonitorAnalyzer::shouldSkipMuonInJetPlusAwayJet(const Lepton *muon,
                                                     const Jet *jet)
{
    return !isValueInRange(muon->p4()->Pt(), _n_muon_pt) ||
           !isValueInRange(jet->p4()->Pt(), _n_jet_pt) ||
           !isValueInRange(fabs(jet->p4()->Eta()), _n_jet_eta);
}

void MonitorAnalyzer::muonIsInJetPlusAwayJet(const Lepton *muon,
                                             const Jet *jet)
{
    double weight = 1;
    const double jet_pt = jet->p4()->Pt();

    if (_reweights &&
        isValueInRange(jet_pt, _n_jet_pt))

        weight = _reweights->GetBinContent(_reweights->FindBin(jet->p4()->Pt(),
                                                               jet->p4()->Eta(),
                                                               _event->primaryVertices()->size()));

    _monitorNMuonInJet->fill(weight, muon, jet, _event->primaryVertices()->size());
    _monitorNMuonInJet->jet()->fillDiscriminator(weight,
            jet->btag(_muonInJetTagger->btag()));

    if (!_are_n_plots_filled)
    {
        _n_primary_vertices->Fill(_event->primaryVertices()->size(), weight);

        for(PrimaryVertices::const_iterator pv = _event->primaryVertices()->begin();
            _event->primaryVertices()->end() != pv;
            ++pv)
        {
            _n_primary_vertices_z->Fill((*pv)->vertex()->Z(), weight);
        }

        _are_n_plots_filled = true;
    }
}

void MonitorAnalyzer::muonIsInJetPlusTaggedAwayJet(const Lepton *muon,
                                                   const Jet *jet)
{
    double weight = 1;
    const double jet_pt = jet->p4()->Pt();

    if (_reweights &&
        isValueInRange(jet_pt, _n_jet_pt))

        weight = _reweights->GetBinContent(_reweights->FindBin(jet->p4()->Pt(),
                                                               jet->p4()->Eta(),
                                                               _event->primaryVertices()->size()));

    _monitorPMuonInJet->fill(weight, muon, jet, _event->primaryVertices()->size());
    _monitorPMuonInJet->jet()->fillDiscriminator(weight,
            jet->btag(_muonInJetTagger->btag()));
}

// PythiaOptionsDelegate interface
//
void MonitorAnalyzer::optionGluonSplittingIsSet(const GluonSplitting &value)
{
    _s8_selector->optionGluonSplittingIsSet(value);
}

void MonitorAnalyzer::optionPtHatIsSet(const Range &value)
{
    _s8_selector->optionPtHatIsSet(value);
}

// Trigger options
//
void MonitorAnalyzer::optionTriggerIsSet(const Trigger &trigger)
{
    _s8_selector->optionTriggerIsSet(trigger);
}

void MonitorAnalyzer::optionSimulateTriggerIsSet(const bool &value)
{
    _s8_selector->optionSimulateTriggerIsSet(value);
}

void MonitorAnalyzer::optionReweightTriggerIsSet(const std::string &filename)
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
void MonitorAnalyzer::optionPrimaryVerticesIsSet(const Range &primary_vertices)
{
    _s8_selector->optionPrimaryVerticesIsSet(primary_vertices);
}



void MonitorAnalyzer::saveGenericPlots(TDirectory *directory) const
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
    _primary_vertices->Write();

    // Generic
    //
    _monitorMuons->save(subdir);
    _monitorJets->save(subdir);
    _monitorLeadingJet->save(subdir);

    directory->cd();
}

void MonitorAnalyzer::saveNPlots(TDirectory *directory) const
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
    _monitorNMuonInJet->save(subdir);
    _n_primary_vertices->Write();
    _n_primary_vertices_z->Write();

    directory->cd();
}

void MonitorAnalyzer::savePPlots(TDirectory *directory) const
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
    _monitorPMuonInJet->save(subdir);

    directory->cd();
}

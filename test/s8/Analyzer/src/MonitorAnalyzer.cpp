/**
 * MonitorAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <stdexcept>

#include <TDirectory.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8GenEvent.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Utility/interface/MonitorPlots.h"
#include "Utility/interface/MuonInJet.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/MonitorAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;

using s8::MonitorAnalyzer;

MonitorAnalyzer::MonitorAnalyzer() throw():
    _minPtHat(0),
    _maxPtHat(0)
{
}

MonitorAnalyzer::~MonitorAnalyzer() throw()
{
}



// Analyzer interface
//
void MonitorAnalyzer::init()
{
    _minPtHat = 0;
    _maxPtHat = 0;

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
    _monitorJets->plot(MonitorBase::PT)->SetBins(46, 0, 230); 
    _monitorLeadingJet.reset(new MonitorJet("leadingjet"));

    // (n)
    //
    _monitorNMuons.reset(new MonitorLepton("nmu"));
    _monitorNJets.reset(new MonitorJet("njet"));
    _monitorNJets->plot(MonitorBase::PT)->SetBins(46, 0, 230); 
    _monitorNDelta.reset(new MonitorDelta("nmujet"));

    // (p)
    //
    _monitorPMuons.reset(new MonitorLepton("pmu"));
    _monitorPJets.reset(new MonitorJet("pjet"));
    _monitorPJets->plot(MonitorBase::PT)->SetBins(46, 0, 230); 
    _monitorPDelta.reset(new MonitorDelta("pmujet"));

    // change binning in the Delta plots
    //
    MonitorDelta *deltas[] = { _monitorNDelta.get(), _monitorPDelta.get() };
    for(int i = 0; 2 > i; ++i)
    {
        deltas[i]->plot(MonitorDelta::PTREL)->SetBins(70, 0, 7);
    }
}

void MonitorAnalyzer::eventDidLoad(const Event *event)
{
    // PtHat cut
    //
    if (_minPtHat &&
        event->gen()->ptHat() < _minPtHat ||
        
        _maxPtHat &&
        event->gen()->ptHat() > _maxPtHat)

        return;

    _pthat->Fill(event->gen()->ptHat());
    _muons->Fill(event->muons()->size());
    _jets->Fill(event->jets()->size());

    for(Leptons::const_iterator muon = event->muons()->begin();
        event->muons()->end() != muon;
        ++muon)

        _monitorMuons->fill(*muon);

    const Jet *leadingJet = 0;
    for(Jets::const_iterator jet = event->jets()->begin();
        event->jets()->end() != jet;
        ++jet)
    {
        _monitorJets->fill(*jet);

        if (leadingJet &&
            leadingJet->p4()->Pt() >= (*jet)->p4()->Pt())

            continue;

        leadingJet = *jet;
    }

    if (leadingJet)
        _monitorLeadingJet->fill(leadingJet);

    (*_muonInJet)(event);
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

// MuonInJetDelegate interface
//
void MonitorAnalyzer::muonIsInJetPlusAwayJet(const Lepton *muon,
                                             const Jet *jet)
{
    _monitorNMuons->fill(muon);
    _monitorNJets->fill(jet);
    _monitorNDelta->fill(muon->p4(), jet->p4());
}

void MonitorAnalyzer::muonIsInJetPlusTaggedAwayJet(const Lepton *muon,
                                                   const Jet *jet)
{
    _monitorPMuons->fill(muon);
    _monitorPJets->fill(jet);
    _monitorPDelta->fill(muon->p4(), jet->p4());
}

// PythiaOptionsDelegate interface
//
void MonitorAnalyzer::optionGluonSplittingIsSet(const GluonSplitting &)
{
    cout << "Gluon Splitting is not supported by the Monitor Analyzer at the moment" << endl;
}

void MonitorAnalyzer::optionPtHatIsSet(const int &min, const int &max)
{
    using std::runtime_error;

    if (0 > min || 0 > max)
        throw runtime_error("Negative value of the pThat cut range is supplied");

    _minPtHat = min;
    _maxPtHat = max;
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
    _monitorNMuons->save(subdir);
    _monitorNJets->save(subdir);
    _monitorNDelta->save(subdir);

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
    _monitorPMuons->save(subdir);
    _monitorPJets->save(subdir);
    _monitorPDelta->save(subdir);

    directory->cd();
}

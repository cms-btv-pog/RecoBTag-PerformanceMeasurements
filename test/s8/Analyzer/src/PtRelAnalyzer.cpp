/**
 * PtRelAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 14, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <iostream>

#include <TDirectory.h>
#include <TH2F.h>
#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "S8Tree/interface/S8TreeInfo.h"
#include "S8Tree/interface/S8TriggerCenter.h"
#include "Utility/interface/MuonInJet.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Analyzer/interface/PtRelAnalyzer.h"

using std::cerr;
using std::clog;
using std::endl;

using s8::PtRelAnalyzer;

PtRelAnalyzer::PtRelAnalyzer() throw():
    _isData(false)
{
}

PtRelAnalyzer::~PtRelAnalyzer() throw()
{
}



// Analyzer interface
//
void PtRelAnalyzer::init()
{
    _awayJetTagger.reset(new TaggerOperatingPoint());

    _muonInJet.reset(new MuonInJet());
    _muonInJet->setDelegate(this);
    _muonInJet->setAwayJetTaggerOperatingPoint(_awayJetTagger.get());

    // Pt bins
    //
    const int nPtBins = 7;
    const double ptBins[] = {30, 40, 50, 60, 70, 80, 120, 230};

    _ptrel.reset(new TH2F("ptrel", "All Jets",
                            50, 0, 5,
                            nPtBins, ptBins));

    _ptrel_b.reset(new TH2F("ptrel_b", "b-Jets",
                            50, 0, 5,
                            nPtBins, ptBins));

    _ptrel_c.reset(new TH2F("ptrel_c", "c-Jets",
                            50, 0, 5,
							nPtBins, ptBins));

    _ptrel_cl.reset(new TH2F("ptrel_cl", "cl-Jets",
                            50, 0, 5,
							nPtBins, ptBins));

    _ptrel_l.reset(new TH2F("ptrel_l", "l-Jets",
                            50, 0, 5,
							nPtBins, ptBins));

    _ptrel_g.reset(new TH2F("ptrel_g", "g splitting",
                            50, 0, 5,
							nPtBins, ptBins));

    // Force Manual memory management
    //
    TH2 *plots[] = { _ptrel.get(),
                     _ptrel_b.get(),
                     _ptrel_c.get(), _ptrel_cl.get(),
                     _ptrel_l.get(),
                     _ptrel_g.get() };

    for(int i = 0; 6 > i; ++i)
    {
        TH2 *plot = *(plots + i);
        plot->SetDirectory(0);
        plot->GetXaxis()->SetTitle("p_{T}^{rel}(#mu-in-jet,jet)");
        plot->GetYaxis()->SetTitle("Jet p_{T} [GeV/c]");
    }
}

void PtRelAnalyzer::eventDidLoad(const Event *event)
{
    (*_muonInJet)(event);
}

void PtRelAnalyzer::print(std::ostream &) const
{
}

void PtRelAnalyzer::save(TDirectory *directory) const
{
    TDirectory *output = directory->mkdir("PtRelAnalyzer");
    if (!output)
    {
        cerr << "PtRelAnalyzer: Failed to create output folder: "
            << "no output is saved" << endl;

        return;
    }

    output->cd();

    _ptrel->Write();
    _ptrel_b->Write();
    _ptrel_c->Write();
    _ptrel_cl->Write();
    _ptrel_l->Write();
    _ptrel_g->Write();
}



// SolverInputOptionsDelegate interface
//
void PtRelAnalyzer::optionAwayTagIsSet(const std::string &tag)
{
    *_awayJetTagger << tag;
}

void PtRelAnalyzer::optionDataIsSet(const bool &value)
{
    _isData = value;
}

// MuonInJetDelegate interface
//
void PtRelAnalyzer::muonIsInJetPlusAwayJet(const Lepton *muon, const Jet *jet)
{
    const double ptrel = muon->p4()->Vect().Perp(jet->p4()->Vect());
    const double pt = jet->p4()->Pt();
    _ptrel->Fill(ptrel, pt);

    if (_isData)
        return;

    switch(abs(jet->flavour()))
    {
        case 5:  _ptrel_b->Fill(ptrel, pt);
                 break;

        case 4: _ptrel_c->Fill(ptrel, pt);
                _ptrel_cl->Fill(ptrel, pt);
                break;

        case 3: // Fall through
        case 2: // Fall through
        case 1: _ptrel_cl->Fill(ptrel, pt);
                _ptrel_l->Fill(ptrel, pt);
                break;

        case 21: _ptrel_g->Fill(ptrel, pt);
                 _ptrel_l->Fill(ptrel, pt);
                 break;
    }
}

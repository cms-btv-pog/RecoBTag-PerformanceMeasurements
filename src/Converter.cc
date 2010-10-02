/**
 * Converter
 * 
 *
 * Created by Samvel Khalatian on Oct 1, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <memory>
#include <stdexcept>

#include <TChain.h>
#include <TH2F.h>
#include <TFile.h>

#include "RecoBTag/PerformanceMeasurements/interface/Event.h"
#include "RecoBTag/PerformanceMeasurements/interface/OperatingPoint.h"

#include "RecoBTag/PerformanceMeasurements/interface/Converter.h"

using std::auto_ptr;
using std::cerr;
using std::endl;
using std::string;
using std::runtime_error;

using s8::Event;
using s8::Jet;
using s8::Muon;
using s8::Converter;

Converter::Plots::Plots(const string &prefix,
                        const string &suffix)
try
    :_isInitializationFailed(false),
     _operatingPoint(0)
{
    // Pt bins
    //
    const int nbins = 3;
    const double bins[] = {30, 50, 80, 230};

    all = new TH2F((prefix + "_pT_" + suffix).c_str(),
                    (prefix + " p_{T}^rel vs p_{T} " +
                     suffix).c_str(),
                    nbins, bins,
                    50, 0, 5);

    tag = new TH2F((prefix + "tag_pT_" + suffix).c_str(),
                   (prefix + " tag p_{T}^rel vs p_{T} " +
                    suffix).c_str(),
                   nbins, bins,
                   50, 0, 5);
}
catch(const std::exception &error)
{
    cerr << "[error] " << error.what() << endl;

    if (all)
        delete all;

    if (tag)
        delete tag;

    _isInitializationFailed = true;
}

Converter::Plots::~Plots()
{
    delete tag;
    delete all;
}

void Converter::Plots::setOperatingPoint(const OperatingPoint &op)
{
    _operatingPoint = op;
}

void Converter::Plots::fill(const Muon *muon, const Jet *jet)
{
    if (_isInitializationFailed)
        throw runtime_error("Plots initialization failed.");

    all->Fill(jet->p4().Pt(),
              muon->p4().Vect().Perp(jet->p4().Vect()));

    if (_operatingPoint < jet->btag(Jet::TCHE))
        tag->Fill(jet->p4().Pt(),
                  muon->p4().Vect().Perp(jet->p4().Vect()));
}

void Converter::Plots::save() const
{
    if (_isInitializationFailed)
        throw runtime_error("Plots initialization failed.");

    all->Write();
    tag->Write();
}

Converter::FlavouredPlots::FlavouredPlots(const string &prefix):
    b(prefix, "b"),
    cl(prefix, "cl")
{
}

void Converter::FlavouredPlots::setOperatingPoint(const OperatingPoint &op)
{
    b.setOperatingPoint(op);
    cl.setOperatingPoint(op);
}

void Converter::FlavouredPlots::fill(const Muon *muon, const Jet *jet)
{
    switch(jet->flavour())
    {
        case 5:  b.fill(muon, jet);
                 break;

        case 4: // Fall through
        case 3: // Fall through
        case 2: // Fall through
        case 1: // Fall through
        case 21: cl.fill(muon, jet);
                 break;
    }
}

void Converter::FlavouredPlots::save() const
{
    b.save();
    cl.save();
}

Converter::Converter() throw():
    _n("n"),
    _p("p")
{
}

Converter::~Converter() throw()
{
}

bool Converter::run(const int argc, char *argv[])
{
    if (3 > argc)
        throw runtime_error("syntax: executable [TCHEM,TCHPL,etc.] tree.root");

    OperatingPoint op;
    op << argv[1];

    _n.setOperatingPoint(op);
    _p.setOperatingPoint(op);

    process(argv[2]);

    return true;
}

void Converter::process(const string &file)
{
    auto_ptr<Event> event(new Event());
    auto_ptr<TChain> chain(new TChain());

    chain->Add((file + "/TreeMaker/s8").c_str());

    chain->SetBranchAddress("event", &event);

    const int entries = chain->GetEntries();
    for(int entry = 0; entries > entry; ++entry)
    {
        if (1000 < entry)
            break;

        chain->GetEntry(entry);

        analyze(event.get());
    }

    auto_ptr<TFile> output(new TFile("s8input.root", "RECREATE"));
    if (!output->IsOpen())
        throw runtime_error("failed to open output file. Results are not saved.");

    _n.save();
    _p.save();
}

void Converter::analyze(const Event *event)
{
    // Skip event if number of jets is less than 2 or there are no muons
    //
    if (2 > event->jets().size() ||
        !event->muons().size())

        return;

    const PrimaryVertex *primaryVertex = &*event->primaryVertices().begin();

    for(JetCollection::const_iterator jet = event->jets().begin();
        event->jets().end() != jet;
        ++jet)
    {
        // find muon (with highest Pt) inside the jet
        //
        const Muon *muonInJet = 0;
        for(MuonCollection::const_iterator muon = event->muons().begin();
            event->muons().end() != muon;
            ++muon)
        {
            if (2 <= abs(muon->vertex().z() - primaryVertex->vertex().z()))
                continue;

            const double deltaR = muon->p4().DeltaR(jet->p4());
            if (0.01 > deltaR ||
                .4 <= deltaR ||
                -1 >= muon->p4().Vect().Perp(jet->p4().Vect()))

                continue;

            if (muonInJet &&
                muonInJet->p4().Pt() >= muon->p4().Pt())
                continue;

            muonInJet = &*muon;
        }

        // Test if muon in jet was found
        //
        if (!muonInJet)
            continue;

        // Find Away jet with highest Pt
        //
        const Jet *awayJet = 0;
        for(JetCollection::const_iterator jetIter = event->jets().begin();
            event->jets().end() != jetIter;
            ++jetIter)
        {
            if (jetIter == jet)
                continue;
            
            if (awayJet &&
                awayJet->p4().Pt() >= jetIter->p4().Pt())
                continue;

            awayJet = &*jetIter;
        }

        // Test if away jet was found
        //
        if (!awayJet)
            continue;

        // (n) muon-jet + away-jet
        //
        _n.fill(muonInJet, &*jet);

        // Test if away jet is tagged
        //
        if (1.19 >= awayJet->btag(Jet::TCHP))
            continue;

        // (p) muon-jet + tagged-away-jet
        //
        _p.fill(muonInJet, &*jet);
    }
}

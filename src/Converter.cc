/**
 * Converter
 * 
 *
 * Created by Samvel Khalatian on Oct 1, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <boost/program_options.hpp>

#include <TChain.h>
#include <TH2F.h>
#include <TFile.h>

#include "RecoBTag/PerformanceMeasurements/interface/Event.h"
#include "RecoBTag/PerformanceMeasurements/interface/OperatingPoint.h"
#include "RecoBTag/PerformanceMeasurements/interface/Plots.h"

#include "RecoBTag/PerformanceMeasurements/interface/Converter.h"

using std::auto_ptr;
using std::cerr;
using std::cout;
using std::endl;
using std::left;
using std::setw;
using std::string;
using std::runtime_error;
using std::vector;

using s8::Event;
using s8::Jet;
using s8::Muon;
using s8::Converter;

namespace po = boost::program_options;

Converter::Converter() throw()
{
    _n = 0;
    _p = 0;
}

Converter::~Converter() throw()
{
    if (_n)
        delete _n;

    if (_p)
        delete _p;
}

bool Converter::run(const int &argc, char **argv)
{
    if (!parseArguments(argc, argv))
        return true;

    OperatingPoint op;
    op << _config.tag;

    if (_config.isData)
    {
        _n = new NonFlavouredPlots("n");
        _p = new NonFlavouredPlots("p");
    }
    else
    {
        _n = new CombinedPlots("n");
        _p = new CombinedPlots("p");
    }

    _n->setOperatingPoint(op);
    _p->setOperatingPoint(op);

    process();

    return true;
}

bool Converter::parseArguments(const int &argc, char **argv)
{
    po::options_description generic("Generic Options");
    generic.add_options()
        ("help,h", "produce help message.")

        ("events,e",
            po::value<int>(&(_config.events))->default_value(0),
            "Maximum number of events to be processed. 0 - all")

        ("data", "Data Input flag. Default: Monte-Carlo")

        ("tag",
            po::value<string>(&(_config.tag))->default_value("TCHEM"),
            "b-Tagger. Default: TCHEM")

        ("output,o",
            po::value<string>(&(_config.output))->default_value("s8input.root"),
            "Output file.")
    ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input",
            po::value<string>(&(_config.input)),
             "input file.")
    ;

    po::positional_options_description positional;
    positional.add("input", 1);

    // Group options
    //
    po::options_description options("Command Line Options");
    options.add(generic).add(hidden);

    // Parse arguments
    //
    po::variables_map arguments;
    po::store(po::command_line_parser(argc, argv).
              options(options).positional(positional).run(), arguments);
    po::notify(arguments);

    if (arguments.count("help"))
    {
        cout << "Usage: <executable> [Options] input" << endl
            << endl;
        cout << generic << endl;

        return false;
    }

    if (0 > _config.events)
        throw runtime_error("Wrong number of events is specified.");

    _config.isData = arguments.count("data") ? true : false;

    if (_config.output.empty())
        cout << "output file is not specified. Results are not saved." << endl;

    if (_config.input.empty())
        throw runtime_error("Input is not specified.");

    cout << "Arguments" << endl;
    cout << setw(25) << left << " [+] Events" << _config.events << endl;
    cout << setw(25) << left << " [+] Input Type"
        << (_config.isData ? "Data" : "Monte-Carlo") << endl;
    cout << setw(25) << left << " [+] Tag" << _config.tag << endl;
    cout << setw(25) << left << " [+] Output" << _config.output << endl;
    cout << setw(25) << left << " [+] Input" << _config.input << endl;
    cout << endl;

    return true;
}

void Converter::process()
{
    auto_ptr<Event> event(new Event());
    auto_ptr<TChain> chain(new TChain());

    chain->Add((_config.input + "/TreeMaker/s8").c_str());

    chain->SetBranchAddress("event", &event);

    int entries = chain->GetEntries();
    if (_config.events &&
        entries > _config.events)

        entries = _config.events;

    const int fractions = 10;
    int fraction = 1;
    int entriesFraction = entries * fraction/ fractions;
    cout << entries << " entries to be processed." << endl
        << endl;

    for(int entry = 0; entries > entry; ++entry)
    {
        if (entriesFraction < entry)
        {
            cout << (100 * fraction / fractions) << "% processed" << endl;

            ++fraction;
            entriesFraction = entries * fraction / fractions;
        }

        chain->GetEntry(entry);

        analyze(event.get());
    }
    cout << "100% processed" << endl
        << endl;

    cout << "saving output." << endl;
    auto_ptr<TFile> output(new TFile(_config.output.c_str(), "RECREATE"));
    if (!output->IsOpen())
        throw runtime_error("failed to open output file. Results are not saved.");

    _n->save(output.get());
    _p->save(output.get());
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
        _n->fill(muonInJet, &*jet);

        // Test if away jet is tagged
        //
        if (1.19 >= awayJet->btag(Jet::TCHP))
            continue;

        // (p) muon-jet + tagged-away-jet
        //
        _p->fill(muonInJet, &*jet);
    }
}

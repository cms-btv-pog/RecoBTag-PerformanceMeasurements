/**
 * Plots
 * s8 
 *
 * Created by Samvel Khalatian on Oct 3, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <stdexcept>

#include <TDirectory.h>
#include <TH2F.h>

#include "RecoBTag/PerformanceMeasurements/interface/Jet.h"
#include "RecoBTag/PerformanceMeasurements/interface/Muon.h"
#include "RecoBTag/PerformanceMeasurements/interface/OperatingPoint.h"

#include "RecoBTag/PerformanceMeasurements/interface/Plots.h"

using std::cerr;
using std::endl;
using std::runtime_error;
using std::string;

using s8::Plots;
using s8::PlotGroup;
using s8::NonFlavouredPlots;
using s8::FlavouredPlots;
using s8::CombinedPlots;

Plots::~Plots()
{
}



PlotGroup::PlotGroup(const string &prefix,
                     const string &suffix)
try
    :_isInitializationFailed(false),
     _operatingPoint(0)
{
    _all = 0;
    _tag = 0;

    // Pt bins
    //
    const int nbins = 3;
    const double bins[] = {30, 50, 80, 230};

    const string newSuffix = suffix.empty() ?
        "_pT" :
        "_pT_" + suffix;

    _all = new TH2F((prefix + newSuffix).c_str(),
                    (prefix + " p_{T}^rel vs p_{T} " +
                        suffix).c_str(),
                    nbins, bins,
                    50, 0, 5);

    _tag = new TH2F((prefix + "tag" + newSuffix).c_str(),
                    (prefix + " tag p_{T}^rel vs p_{T} " +
                        suffix).c_str(),
                    nbins, bins,
                    50, 0, 5);
}
catch(const std::exception &error)
{
    cerr << "[error] " << error.what() << endl;

    if (_all)
    {
        delete _all;
        _all = 0;
    }

    if (_tag)
    {
        delete _tag;
        _tag = 0;
    }

    _isInitializationFailed = true;
}

PlotGroup::~PlotGroup()
{
    if (_isInitializationFailed)
        return;

    delete _tag;
    delete _all;
}

void PlotGroup::setOperatingPoint(const OperatingPoint &op)
{
    _operatingPoint = op;
}

void PlotGroup::fill(const Muon *muon, const Jet *jet)
{
    if (_isInitializationFailed)
        throw runtime_error("Plots initialization failed.");

    _all->Fill(jet->p4().Pt(),
               muon->p4().Vect().Perp(jet->p4().Vect()));

    if (_operatingPoint < jet->btag(Jet::TCHE))
        _tag->Fill(jet->p4().Pt(),
                   muon->p4().Vect().Perp(jet->p4().Vect()));
}

void PlotGroup::save(TDirectory *dir) const
{
    if (_isInitializationFailed)
        throw runtime_error("Plots initialization failed.");

    _all->Write();
    _tag->Write();
}



NonFlavouredPlots::NonFlavouredPlots(const string &prefix):
    _plots(prefix)
{
}

void NonFlavouredPlots::setOperatingPoint(const OperatingPoint &op)
{
    _plots.setOperatingPoint(op);
}

void NonFlavouredPlots::fill(const Muon *muon, const Jet *jet)
{
    _plots.fill(muon, jet);
}

void NonFlavouredPlots::save(TDirectory *dir) const
{
    const std::string subdirName = "muon_in_jet";

    if (dir->FindKey(subdirName.c_str()))
        dir->cd(subdirName.c_str());
    else
    {
        TDirectory *subdir = dir->mkdir(subdirName.c_str());
        if (!subdir)
            throw runtime_error("Failed to create subdirectory: " + subdirName);

        subdir->cd();
    }

    _plots.save(dir);

    dir->cd();
}



FlavouredPlots::FlavouredPlots(const string &prefix):
    _b(prefix, "b"),
    _cl(prefix, "cl")
{
}

void FlavouredPlots::setOperatingPoint(const OperatingPoint &op)
{
    _b.setOperatingPoint(op);
    _cl.setOperatingPoint(op);
}

void FlavouredPlots::fill(const Muon *muon, const Jet *jet)
{
    switch(jet->flavour())
    {
        case 5:  _b.fill(muon, jet);
                 break;

        case 4: // Fall through
        case 3: // Fall through
        case 2: // Fall through
        case 1: // Fall through
        case 21: _cl.fill(muon, jet);
                 break;
    }
}

void FlavouredPlots::save(TDirectory *dir) const
{
    const std::string subdirName = "MCTruth";

    if (dir->FindKey(subdirName.c_str()))
        dir->cd(subdirName.c_str());
    else
    {
        TDirectory *subdir = dir->mkdir(subdirName.c_str());
        if (!subdir)
            throw runtime_error("Failed to create subdirectory: " + subdirName);

        subdir->cd();
    }

    _b.save(dir);
    _cl.save(dir);

    dir->cd();
}



CombinedPlots::CombinedPlots(const string &prefix):
    _flavoured(prefix),
    _nonFlavoured(prefix)
{
}

void CombinedPlots::setOperatingPoint(const OperatingPoint &op)
{
    _flavoured.setOperatingPoint(op);
    _nonFlavoured.setOperatingPoint(op);
}

void CombinedPlots::fill(const s8::Muon *muon, const s8::Jet *jet)
{
    _flavoured.fill(muon, jet);
    _nonFlavoured.fill(muon, jet);
}

void CombinedPlots::save(TDirectory *dir) const
{
    _flavoured.save(dir);
    _nonFlavoured.save(dir);
}

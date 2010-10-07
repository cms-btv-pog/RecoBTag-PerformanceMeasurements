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
    _allPt = 0;
    _tagPt = 0;

    _allEta = 0;
    _tagEta = 0;

    // Pt bins
    //
    const int nPtBins = 3;
    const double ptBins[] = {30, 50, 80, 230};

    // Eta bins
    //
    const int nEtaBins = 4;
    const double etaBins[] = {0.0, 0.5, 1.0, 1.5, 2.5};

    string newSuffix = suffix.empty() ?
        "_pT" :
        "_pT_" + suffix;

    _allPt = new TH2F((prefix + newSuffix).c_str(),
                        (prefix + " p_{T}^rel vs p_{T} " +
                            suffix).c_str(),
                        nPtBins, ptBins,
                        50, 0, 5);

    _tagPt = new TH2F((prefix + "tag" + newSuffix).c_str(),
                        (prefix + " tag p_{T}^rel vs p_{T} " +
                            suffix).c_str(),
                        nPtBins, ptBins,
                        50, 0, 5);

    newSuffix = suffix.empty() ?
        "_eta" :
        "_eta_" + suffix;

    _allEta = new TH2F((prefix + newSuffix).c_str(),
                        (prefix + " p_{T}^rel vs #eta " +
                            suffix).c_str(),
                        nEtaBins, etaBins,
                        50, 0, 5);

    _tagEta = new TH2F((prefix + "tag" + newSuffix).c_str(),
                        (prefix + " tag p_{T}^rel vs #eta " +
                            suffix).c_str(),
                        nEtaBins, etaBins,
                        50, 0, 5);
}
catch(const std::exception &error)
{
    cerr << "[error] " << error.what() << endl;

    if (_allPt)
    {
        delete _allPt;
        _allPt = 0;
    }

    if (_tagPt)
    {
        delete _tagPt;
        _tagPt = 0;
    }

    if (_allEta)
    {
        delete _allEta;
        _allEta = 0;
    }

    if (_tagEta)
    {
        delete _tagEta;
        _tagEta = 0;
    }

    _isInitializationFailed = true;
}

PlotGroup::~PlotGroup()
{
    if (_isInitializationFailed)
        return;

    delete _tagEta;
    delete _allEta;

    delete _tagPt;
    delete _allPt;
}

void PlotGroup::setOperatingPoint(const OperatingPoint &op)
{
    _operatingPoint = op;
}

void PlotGroup::fill(const Muon *muon, const Jet *jet)
{
    if (_isInitializationFailed)
        throw runtime_error("Plots initialization failed.");

    _allPt->Fill(jet->p4().Pt(),
               muon->p4().Vect().Perp(jet->p4().Vect() + muon->p4().Vect()));

    _allEta->Fill(jet->p4().Eta(),
                  muon->p4().Vect().Perp(jet->p4().Vect() + muon->p4().Vect()));

    if (_operatingPoint < jet->btag(Jet::TCHE))
    {
        _tagPt->Fill(jet->p4().Pt(),
                   muon->p4().Vect().Perp(jet->p4().Vect() +
                                          muon->p4().Vect()));

        _tagEta->Fill(jet->p4().Eta(),
                   muon->p4().Vect().Perp(jet->p4().Vect() +
                                          muon->p4().Vect()));
    }
}

void PlotGroup::save(TDirectory *dir) const
{
    if (_isInitializationFailed)
        throw runtime_error("Plots initialization failed.");

    _allPt->Write();
    _tagPt->Write();

    _allEta->Write();
    _tagEta->Write();
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

/**
 * SolverLeptonInJetPlots
 * s8 
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <TDirectory.h>
#include <TH2F.h>
#include <TLorentzVector.h>

#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Utility/interface/TaggerOperatingPoint.h"
#include "Utility/interface/LeptonInJetDiscriminator.h"

#include "Utility/interface/SolverLeptonInJetPlots.h"

using std::cerr;
using std::endl;
using std::string;
using std::runtime_error;

using s8::Jet;
using s8::Lepton;

using s8::plot::LeptonInJetPlots;
using s8::plot::LeptonInJetGroup;
using s8::plot::LeptonInJetPlotCategory;
using s8::plot::LeptonInJetPlotGroup;
using s8::plot::NonFlavouredLeptonInJetPlotGroup;
using s8::plot::FlavouredLeptonInJetPlotGroup;
using s8::SolverLeptonInJetPlots;

LeptonInJetPlots::~LeptonInJetPlots()
{
}



LeptonInJetGroup::~LeptonInJetGroup()
{
}



LeptonInJetPlotCategory::LeptonInJetPlotCategory(const string &prefix,
                                                 const string &suffix):
    _didNotInitialize(true),
    _prefix(prefix),
    _suffix(suffix)
{
    _discriminator = 0;

}

LeptonInJetPlotCategory::~LeptonInJetPlotCategory()
{
}

void LeptonInJetPlotCategory::init()
{
    if (!_discriminator)
        throw runtime_error("Can not init: discriminator is not set");

    // Check if object already initialized
    //
    if (!_didNotInitialize)
        return;

    // Pt bins
    //
    const int nPtBins = 7;
    const double ptBins[] = {20, 30, 40, 50, 60, 80, 140, 230};

    // Eta bins
    //
    const int nEtaBins = 2;
    const double etaBins[] = {0.0, 1.4, 2.5};

    // Phi bins
    //
    const int nPhiBins = 4;
    const double phiBins[] = {0.0, 0.785, 1.57, 2.355, 3.14};

    // PtRel bins
    //
    BinGroup discriminatorBins = _discriminator->bins();

    string newSuffix = _suffix.empty() ?
        "_pT" :
        "_pT_" + _suffix;

    _pt.reset(new TH2F((_prefix + newSuffix).c_str(),
                       (_prefix + " p_{T} " +
                           _suffix).c_str(),
                       nPtBins, ptBins,
                       discriminatorBins.bins,
                       discriminatorBins.min,
                       discriminatorBins.max));
    _pt->Sumw2();

    newSuffix = _suffix.empty() ?
        "_eta" :
        "_eta_" + _suffix;

    _eta.reset(new TH2F((_prefix + newSuffix).c_str(),
                        (_prefix + " #eta " +
                            _suffix).c_str(),
                        nEtaBins, etaBins,
                        discriminatorBins.bins,
                        discriminatorBins.min,
                        discriminatorBins.max));
    _eta->Sumw2();
    
    newSuffix = _suffix.empty() ?
        "_phi" :
        "_phi_" + _suffix;

    _phi.reset(new TH2F((_prefix + newSuffix).c_str(),
                        (_prefix + " #phi " +
                            _suffix).c_str(),
                        nPhiBins, phiBins,
                        discriminatorBins.bins,
                        discriminatorBins.min,
                        discriminatorBins.max));
    _phi->Sumw2();

    _didNotInitialize = false;
}

void LeptonInJetPlotCategory::setDiscriminator(LeptonInJetDiscriminator *discriminator)
{
    _discriminator = discriminator;
}

void LeptonInJetPlotCategory::fill(const Lepton *lepton, const Jet *jet,
                                   const double &weight)
{
    if (_didNotInitialize)
        throw runtime_error("LeptonInJetPlotCategory did not initialize");

    if (20 > jet->p4()->Pt())
        return;

    const double discriminant = _discriminator->operator()(lepton, jet);

    _pt->Fill(jet->p4()->Pt(), discriminant, weight);
    _eta->Fill(fabs(jet->p4()->Eta()), discriminant, weight);
    _phi->Fill(fabs(jet->p4()->Phi()), discriminant, weight);
}

void LeptonInJetPlotCategory::save(TDirectory *) const
{
    if (_didNotInitialize)
    {
        cerr << "LeptonInJetPlotCategory did not initialize" << endl;

        return;
    }

    _pt->Write();
    _eta->Write();
    _phi->Write();
}



LeptonInJetPlotGroup::LeptonInJetPlotGroup(const string &prefix,
                     const string &suffix)
{
    _operatingPoint = 0;

    _all.reset(new LeptonInJetPlotCategory(prefix, suffix));
    _tag.reset(new LeptonInJetPlotCategory(prefix + "tag", suffix));
}

LeptonInJetPlotGroup::~LeptonInJetPlotGroup()
{
}

void LeptonInJetPlotGroup::init()
{
    _all->init();
    _tag->init();
}

void LeptonInJetPlotGroup::setDiscriminator(LeptonInJetDiscriminator *discriminator)
{
    _all->setDiscriminator(discriminator);
    _tag->setDiscriminator(discriminator);
}

void LeptonInJetPlotGroup::setTaggerOperatingPoint(const TaggerOperatingPoint
                                                *operatingPoint)
{
    _operatingPoint = operatingPoint;
}

void LeptonInJetPlotGroup::fill(const Lepton *lepton, const Jet *jet,
                                const double &weight)
{
    if (!_operatingPoint)
        throw runtime_error("LeptonInJetPlotGroup: tagger operating point is not set");

    _all->fill(lepton, jet, weight);

    if (*_operatingPoint < jet->btag(_operatingPoint->btag()))
        _tag->fill(lepton, jet, weight);
}

void LeptonInJetPlotGroup::save(TDirectory *dir) const
{
    _all->save(dir);
    _tag->save(dir);
}



NonFlavouredLeptonInJetPlotGroup::NonFlavouredLeptonInJetPlotGroup(const string
                                                                    &prefix)
{
    _plots.reset(new LeptonInJetPlotGroup(prefix));
}

void NonFlavouredLeptonInJetPlotGroup::init()
{
    _plots->init();
}

void NonFlavouredLeptonInJetPlotGroup::setDiscriminator(LeptonInJetDiscriminator *discriminator)
{
    _plots->setDiscriminator(discriminator);
}

void NonFlavouredLeptonInJetPlotGroup::setTaggerOperatingPoint(const TaggerOperatingPoint
                                                            *operatingPoint)
{
    _plots->setTaggerOperatingPoint(operatingPoint);
}

void NonFlavouredLeptonInJetPlotGroup::fill(const Lepton *lepton,
                                            const Jet *jet,
                                            const double &weight)
{
    _plots->fill(lepton, jet, weight);
}

void NonFlavouredLeptonInJetPlotGroup::save(TDirectory *directory) const
{
    const std::string subdirName = "lepton_in_jet";

    // Check if subdir already exists
    //
    TDirectory *subdir = 0;
    TObject *object = directory->FindObject(subdirName.c_str());
    if (object)
    {
        subdir = dynamic_cast<TDirectory *>(object);
        subdir->cd();
    }
    else
    {
        // Subdir does not exit. Create one.
        //
        subdir = directory->mkdir(subdirName.c_str());
        if (!subdir)
        {
            cerr << "Failed to create TDirectory: " << subdirName
                << endl;

            return;
        }

        subdir->cd();
    }

    _plots->save(subdir);

    directory->cd();
}



FlavouredLeptonInJetPlotGroup::FlavouredLeptonInJetPlotGroup(const string &prefix)
{
    _b.reset(new LeptonInJetPlotGroup(prefix, "b"));
    _cl.reset(new LeptonInJetPlotGroup(prefix, "cl"));
}

void FlavouredLeptonInJetPlotGroup::init()
{
    _b->init();
    _cl->init();
}

void FlavouredLeptonInJetPlotGroup::setDiscriminator(LeptonInJetDiscriminator *discriminator)
{
    _b->setDiscriminator(discriminator);
    _cl->setDiscriminator(discriminator);
}

void FlavouredLeptonInJetPlotGroup::setTaggerOperatingPoint(const TaggerOperatingPoint *operatingPoint)
{
    _b->setTaggerOperatingPoint(operatingPoint);
    _cl->setTaggerOperatingPoint(operatingPoint);
}

void FlavouredLeptonInJetPlotGroup::fill(const Lepton *lepton, const Jet *jet,
                                         const double &weight)
{
    // PDG IDs: http://pdg.lbl.gov/mc_particle_id_contents.html
    //
    // QUARKS:
    //      d   1 
    //      u   2 
    //      s   3 
    //      c   4 
    //      b   5 
    //
    //  GAUGE BOSONS
    //      g   21
    //
    // cl is a reference to c-quark + light-quarks where the latter one
    // includes uds-quarks + gluon
    //
    switch(abs(jet->flavour()))
    {
        case 5:  _b->fill(lepton, jet, weight);
                 break;

        case 4: // Fall through
        case 3: // Fall through
        case 2: // Fall through
        case 1: // Fall through
        case 21: _cl->fill(lepton, jet, weight);
                 break;
    }
}

void FlavouredLeptonInJetPlotGroup::save(TDirectory *directory) const
{
    const std::string subdirName = "MCTruth";

    // Check if folder already exists in the output directory
    //
    TDirectory *subdir = 0;
    TObject *object = directory->FindObject(subdirName.c_str());
    if (object)
    {
        subdir = dynamic_cast<TDirectory *>(object);
        subdir->cd();
    }
    else
    {
        // There is no folder: create one
        //
        subdir = directory->mkdir(subdirName.c_str());
        if (!subdir)
        {
            cerr << "Failed to create TDirectory: " << subdirName
                << endl;

            return;
        }

        subdir->cd();
    }

    _b->save(subdir);
    _cl->save(subdir);

    directory->cd();
}



SolverLeptonInJetPlots::SolverLeptonInJetPlots(const string &prefix):
    _isFlavoured(true)
{
    _flavoured.reset(new plot::FlavouredLeptonInJetPlotGroup(prefix));
    _nonFlavoured.reset(new plot::NonFlavouredLeptonInJetPlotGroup(prefix));
}

void SolverLeptonInJetPlots::init()
{
    _flavoured->init();
    _nonFlavoured->init();
}

bool SolverLeptonInJetPlots::isFlavoured() const
{
    return _isFlavoured;
}

void SolverLeptonInJetPlots::setIsFlavoured(const bool &flag)
{
    _isFlavoured = flag;
}

void SolverLeptonInJetPlots::setDiscriminator(LeptonInJetDiscriminator *discriminator)
{
    _flavoured->setDiscriminator(discriminator);
    _nonFlavoured->setDiscriminator(discriminator);
}

void SolverLeptonInJetPlots::setTaggerOperatingPoint(const TaggerOperatingPoint *operatingPoint)
{
    _flavoured->setTaggerOperatingPoint(operatingPoint);
    _nonFlavoured->setTaggerOperatingPoint(operatingPoint);
}

void SolverLeptonInJetPlots::fill(const Lepton *lepton, const Jet *jet,
                                  const double &weight)
{
    if (_isFlavoured)
        _flavoured->fill(lepton, jet, weight);

    _nonFlavoured->fill(lepton, jet, weight);
}

void SolverLeptonInJetPlots::save(TDirectory *dir) const
{
    if (_isFlavoured)
        _flavoured->save(dir);

    _nonFlavoured->save(dir);
}

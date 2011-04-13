/**
 * TopTreeMaker
 * 
 *
 * Created by Samvel Khalatian on August 20, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include <TTree.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "EDModule/Analyzer/interface/Tools.h"
#include "Tree/Top/interface/TopElectron.h"
#include "Tree/Top/interface/TopElectronIsolation.h"
#include "Tree/Top/interface/TopEventID.h"
#include "Tree/Top/interface/TopJet.h"
#include "Tree/Top/interface/TopJetEnergy.h"
#include "Tree/Top/interface/TopMuon.h"
#include "Tree/Top/interface/TopMuonIsolation.h"

#include "EDModule/Analyzer/interface/TopTreeMaker.h"

using std::cout;
using std::endl;
using std::string;

using cms::Exception;

using edm::errors::ErrorCodes;
using edm::Handle;
using edm::InputTag;
using edm::LogWarning;
using edm::LogInfo;
using edm::ParameterSet;

using reco::BeamSpot;

TopTreeMaker::TopTreeMaker(const edm::ParameterSet &config)
{
    _beamSpots = config.getParameter<string>("beamSpots");
    _electrons = config.getParameter<string>("electrons");
    _genParticles = config.getParameter<string>("genParticles");
    _jets = config.getParameter<string>("jets");
    _mets = config.getParameter<string>("mets");
    _muons = config.getParameter<string>("muons");
}

TopTreeMaker::~TopTreeMaker()
{
}

void TopTreeMaker::beginJob()
{
    edm::Service<TFileService> fileService;

    _tree = fileService->make<TTree>("top", "Top ttmuj tree.");

    _event.reset(new top::Event());
    _tree->Branch("event", _event.get(), 32000, 0);
}

void TopTreeMaker::endJob()
{
    if (!_event.get())
        return;

    // Note: Event should be destroyed after ROOT file is written and closed.
    //
    _event.reset();
}

void TopTreeMaker::analyze(const edm::Event &event, const edm::EventSetup &)
{
    if (!_event.get())
        return;

    using namespace top::tools;

    using pat::ElectronCollection;
    using pat::JetCollection;
    using pat::METCollection;
    using pat::MuonCollection;

    // Extract BeamSpot
    //
    Handle<BeamSpot> beamSpot;
    event.getByLabel(InputTag(_beamSpots), beamSpot);

    if (!beamSpot.isValid())
        throw Exception("NotFound")
            << "failed to extract BeamSpot." << endl;

    // Extract MET
    //
    Handle<METCollection> mets;
    event.getByLabel(InputTag(_mets), mets);

    if (!mets.isValid())
        throw Exception("NotFound")
            << "failed to extract METs." << endl;

    // Extract Muons
    //
    Handle<MuonCollection> muons;
    event.getByLabel(InputTag(_muons), muons);

    if (!muons.isValid())
        throw Exception("NotFound")
            << "failed to extract muons." << endl;

    // Extract Jets
    //
    Handle<JetCollection> jets;
    event.getByLabel(InputTag(_jets), jets);

    if (!jets.isValid())
        throw Exception("NotFound")
            << "failed to extract jets." << endl;

    // Extract Electrons
    //
    Handle<ElectronCollection> electrons;
    event.getByLabel(InputTag(_electrons), electrons);

    if (!electrons.isValid())
        throw Exception("NotFound")
            << "failed to extract electrons." << endl;

    // Extract Monte-Carlo GenParticles
    //
    /*
    Handle<GenParticleCollection> genParticles;
    event.getByLabel(InputTag(_genParticles), genParticles);

    if (!genParticles.isValid())
    {
        LogWarning("TopTreeMaker")
            << "Failed to extract genParticles.";

        return;
    }
    */

    _event->reset();

    {
        top::EventID id;
        id.setRun(event.id().run());
        id.setLumiBlock(event.id().luminosityBlock());
        id.setEvent(event.id().event());

        _event->setID(id);
    }

    setP4(_event->met().p4(), mets->begin()->p4());

    // Process all Muons
    //
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        // Only Muons with basic cuts are saved:
        if (muon->isGlobalMuon() &&
            10 <= muon->pt() &&
            2.5 >= fabs(muon->eta()))
        {
            top::Muon topMuon;
            topMuon.setIsGlobal(muon->isGlobalMuon());
            topMuon.setIsTracker(muon->isTrackerMuon());

            topMuon.setEta(muon->eta());
            topMuon.setPhi(muon->phi());

            topMuon.setMatches(muon->numberOfMatches());

            setP4(topMuon.p4(), muon->p4());

            setIsolation(topMuon, top::Muon::R03, muon->isolationR03());
            setIsolation(topMuon, top::Muon::R05, muon->isolationR05());

            // Inner Track is only available for the Tracker Muons
            if (muon->isTrackerMuon())
            {
                top::ImpactParameter ip;

                ip.setValue(muon->innerTrack()->dxy(beamSpot->position()));

                ip.setError(sqrt(muon->innerTrack()->d0Error() *
                                 muon->innerTrack()->d0Error() +
                                 0.5 *
                                 beamSpot->BeamWidthX() *
                                 beamSpot->BeamWidthX() +
                                 0.5 *
                                 beamSpot->BeamWidthY() *
                                 beamSpot->BeamWidthY()));

                topMuon.setImpactParameter(top::Muon::BS2D, ip);

                topMuon.setInnerValidHits(muon->innerTrack()->numberOfValidHits());
            }

            // Next properties are only defined for the Global Muon
            topMuon.setOuterValidHits(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
            topMuon.setChi2(muon->globalTrack()->chi2());
            topMuon.setNdof(muon->globalTrack()->ndof());

            _event->muons().push_back(topMuon);
        }
    } // End loop over muons

    // Process All Jets
    //
    for(JetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        // Select only energetic jets
        if (30 <= jet->pt() &&
            2.4 >= fabs(jet->eta()))
        {
            top::Jet topJet;

            setP4(topJet.p4(), jet->p4());
            setEnergy(topJet, jet->caloSpecific());

            topJet.setEta(jet->eta());
            topJet.setPhi(jet->phi());

            topJet.setHits90(jet->jetID().n90Hits);
            topJet.setHpd(jet->jetID().fHPD);

            _event->jets().push_back(topJet);
        }
    }

    // Process All Electrons
    for(ElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        // Select electrons to be stored
        if (15 <= electron->et() &&
            2.5 >= fabs(electron->eta()))
        {
            top::Electron topElectron;

            setP4(topElectron.p4(), electron->p4());

            topElectron.setEta(electron->eta());
            topElectron.setPhi(electron->phi());

            setIsolation(topElectron,
                         top::Electron::R03,
                         electron->isolationVariables03());

            setIsolation(topElectron,
                         top::Electron::R04,
                         electron->isolationVariables04());

            _event->electrons().push_back(topElectron);
        }
    }

    _tree->Fill();
}

DEFINE_FWK_MODULE(TopTreeMaker);

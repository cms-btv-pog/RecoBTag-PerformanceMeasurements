/**
 * TreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <vector>

#include <TTree.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoBTag/PerformanceMeasurements/interface/EventID.h"
#include "RecoBTag/PerformanceMeasurements/interface/GenParticle.h"
#include "RecoBTag/PerformanceMeasurements/interface/Jet.h"
#include "RecoBTag/PerformanceMeasurements/interface/Muon.h"
#include "RecoBTag/PerformanceMeasurements/interface/Tools.h"

#include "RecoBTag/PerformanceMeasurements/interface/TreeMaker.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using edm::errors::ErrorCodes;
using edm::Handle;
using edm::InputTag;
using edm::LogWarning;
using edm::LogInfo;
using edm::ParameterSet;

using reco::Vertex;

using namespace s8::tools;

TreeMaker::TreeMaker(const edm::ParameterSet &config)
{
    _primaryVertices = config.getParameter<string>("primaryVertices");
    _jets = config.getParameter<string>("jets");
    _muons = config.getParameter<string>("muons");
}

TreeMaker::~TreeMaker()
{
}

void TreeMaker::beginJob()
{
    edm::Service<TFileService> fileService;

    _tree = fileService->make<TTree>("s8", "System8 tree.");

    _event.reset(new s8::Event());
    _tree->Branch("event", _event.get(), 32000, 0);
}

void TreeMaker::endJob()
{
    if (!_event.get())
        return;

    // Note: Event should be destroyed after ROOT file is written and closed.
    _event.reset();
}

void TreeMaker::analyze(const edm::Event &event, const edm::EventSetup &)
{
    if (!_event.get())
        return;

    using pat::JetCollection;
    using pat::MuonCollection;

    // Extract Muons
    //
    Handle<MuonCollection> muons;
    event.getByLabel(InputTag(_muons), muons);

    if (!muons.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Muons.";

        return;
    }

    // Extract Jets
    //
    Handle<JetCollection> jets;
    event.getByLabel(InputTag(_jets), jets);

    if (!jets.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Jets.";

        return;
    }

    // Primary Vertices
    //
    typedef vector<Vertex> PVCollection;

    Handle<PVCollection> primaryVertices;
    event.getByLabel(InputTag(_primaryVertices), primaryVertices);

    if (!primaryVertices.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Primary Vertices.";

        return;
    }

    if (primaryVertices->empty())
    {
        LogWarning("TreeMaker")
            << "primary vertices collection is empty.";

        return;
    }

    _event->reset();

    {
        s8::EventID &id = _event->id();

        id.setRun(event.id().run());
        id.setLumiBlock(event.id().luminosityBlock());
        id.setEvent(event.id().event());
    }

    // Process Primary Vertices
    //
    for(PVCollection::const_iterator vertex = primaryVertices->begin();
        primaryVertices->end() != vertex;
        ++vertex)
    {
        if (!isGoodPrimaryVertex(*vertex, event.isRealData()))
            continue;

        s8::PrimaryVertex s8Vertex;

        setVertex(s8Vertex.vertex(), vertex->position());
        s8Vertex.setNdof(vertex->ndof());
        s8Vertex.setRho(vertex->position().Rho());

        _event->primaryVertices().push_back(s8Vertex);
    }

    // Process all Muons
    //
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        // Only Muons with basic cuts are saved:
        //
        if (1 >= muon->numberOfMatches())
            continue;

        s8::Muon s8Muon;

        setP4(s8Muon.p4(), muon->p4());
        setVertex(s8Muon.vertex(), muon->vertex());

        s8Muon.impactParameter().first = muon->dB();
        s8Muon.impactParameter().second = muon->edB();

        if (muon->genLepton())
        {
            s8::GenParticle &s8GenParticle = s8Muon.genParticle();

            setP4(s8GenParticle.p4(), muon->genLepton()->p4());
            setVertex(s8GenParticle.vertex(), muon->genLepton()->vertex());

            s8GenParticle.setId(muon->genLepton()->pdgId());
            if (muon->genLepton()->mother())
                s8GenParticle.setParentId(muon->genLepton()->mother()->pdgId());
        }

        _event->muons().push_back(s8Muon);
    } // End loop over muons

    // Process All Jets
    //
    for(JetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        using s8::Jet;
        Jet s8Jet;

        setP4(s8Jet.p4(), jet->p4());

        s8Jet.setFlavour(jet->partonFlavour());
        s8Jet.setTracks(jet->associatedTracks().size());

        // Save b-taggers
        //
        s8Jet.setBTag(Jet::TCHE,
                      jet->bDiscriminator("trackCountingHighEffBJetTags"));

        s8Jet.setBTag(Jet::TCHP,
                      jet->bDiscriminator("trackCountingHighPurBJetTags"));

        s8Jet.setBTag(Jet::JP,
                      jet->bDiscriminator("jetProbabilityBJetTags"));

        s8Jet.setBTag(Jet::SSV,
                      jet->bDiscriminator("simpleSecondaryVertexBJetTags"));

        s8Jet.setBTag(Jet::SSVHE,
                      jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));

        s8Jet.setBTag(Jet::SSVHP,
                      jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));

        _event->jets().push_back(s8Jet);
    }

    _tree->Fill();
}

bool TreeMaker::isGoodPrimaryVertex(const Vertex &vertex,
                                    const bool &isData)
{
    return !vertex.isFake() &&
            4 <= vertex.ndof() &&
            (isData ? 24 : 15) >= fabs(vertex.z()) &&
            2 >= fabs(vertex.position().Rho());
}

DEFINE_FWK_MODULE(TreeMaker);

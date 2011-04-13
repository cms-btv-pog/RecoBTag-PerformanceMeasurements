/**
 * S8TreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/regex.hpp>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "Tree/System8/interface/S8TreeInfo.h"
#include "Tree/System8/interface/S8EventID.h"
#include "Tree/System8/interface/S8GenEvent.h"
#include "Tree/System8/interface/S8GenParticle.h"
#include "Tree/System8/interface/S8Jet.h"
#include "Tree/System8/interface/S8Lepton.h"
#include "Tree/System8/interface/S8PrimaryVertex.h"
#include "Tree/System8/interface/S8TreeInfo.h"
#include "Tree/System8/interface/S8Trigger.h"
#include "Tree/System8/interface/S8TriggerCenter.h"
#include "EDModule/Analyzer/interface/Tools.h"

#include "EDModule/Analyzer/interface/S8TreeMaker.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using edm::errors::ErrorCodes;
using edm::Handle;
using edm::InputTag;
using edm::LogWarning;
using edm::LogInfo;
using edm::ParameterSet;
using edm::TriggerResults;
using edm::TriggerNames;

using reco::Vertex;

using s8::TreeInfo;
using s8::TriggerCenter;

using namespace top::tools;

void fillSplit(const int &status, boost::tribool &split)
{
    switch(status)
    {
        case 3: split = false;
                break;

        case 2: if (boost::indeterminate(split))
                    split = true;

                break;
    }
}

S8TreeMaker::S8TreeMaker(const edm::ParameterSet &config):
    _didInitializeHltConfigProvider(false)
{
    _primaryVertices = config.getParameter<string>("primaryVertices");
    _jets = config.getParameter<string>("jets");
    _muons = config.getParameter<string>("muons");
    _electrons = config.getParameter<string>("electrons");
    _triggers = config.getParameter<string>("triggers");

    _jetSelector = config.getParameter<ParameterSet>("jetSelector");

    _isPythia = config.getParameter<bool>("isPythia");
    _saveTriggers = config.getParameter<bool>("saveTriggers");
}

S8TreeMaker::~S8TreeMaker()
{
}

void S8TreeMaker::beginJob()
{
    cout << "Start BeginJob" << endl;

    edm::Service<TFileService> fileService;

    _treeInfo.reset(new TreeInfo);
    _triggerCenter.reset(new TriggerCenter);

    _tree = fileService->make<TTree>("s8", "System8 tree.");

    // Prepare Branches
    //
    _eventID.reset(new s8::EventID());
    _tree->Branch("eventID", _eventID.get());

    _genEvent.reset(new s8::GenEvent());
    _tree->Branch("genEvent", _genEvent.get());

    _s8Electrons.reset(new s8::Leptons());
    _tree->Branch("electrons.", _s8Electrons.get());

    _s8Jets.reset(new s8::Jets());
    _tree->Branch("jets", _s8Jets.get());

    _s8Muons.reset(new s8::Leptons());
    _tree->Branch("muons.", _s8Muons.get());

    _s8PrimaryVertices.reset(new s8::PrimaryVertices());
    _tree->Branch("primaryVertices", _s8PrimaryVertices.get());

    _s8Triggers.reset(new s8::Triggers());
    _tree->Branch("triggers", _s8Triggers.get());

    _didInitializeHltConfigProvider = false;

    cout << "BeginJob is over" << endl;
}

void S8TreeMaker::endJob()
{
    if (!_tree)
        return;

    // Tree Info is disabled for the moment until Hadd is fixed
    //
    edm::Service<TFileService> fileService;
    TDirectory *dir = fileService->cd();
    dir->WriteObject(_treeInfo.get(), "s8info");
    dir->WriteObject(_triggerCenter.get(), "s8triggers");
}

void S8TreeMaker::beginRun(const edm::Run &run,
                           const edm::EventSetup &eventSetup)
{
    using s8::Trigger;

    // Initialize HLT Config Provider for new Run
    //
    bool didChange = true;
    if (!_hltConfigProvider.init(run,
                                 eventSetup,
                                 InputTag(_triggers).process(),
                                 didChange))

        throw cms::Exception("S8TreeMaker")
            << "Failed to initialize HLTConfig for: " << _triggers;

    // Test if Trigger Menu has changed
    //
    if (!didChange &&
        _didInitializeHltConfigProvider)

        return;

    _hlts.clear();

    // Triggers are found and changed
    //
    typedef std::vector<std::string> Triggers;

    const Triggers &triggerNames = _hltConfigProvider.triggerNames();

    using s8::tools::make_hash;

    TriggerCenter::TriggerMap &triggers = _triggerCenter->triggers();

    for(Triggers::const_iterator trigger = triggerNames.begin();
        triggerNames.end() != trigger;
        ++trigger)
    {
        smatch matches;
        if (!regex_match(*trigger, matches,
                         regex("^(\\w+?)(?:_[vV](\\d+))?$")))
        {
            cout << "Do not understand Trigger Name: " << *trigger
                << endl;

            continue;
        }

        string trigger_name = matches[1];
        boost::to_lower(trigger_name);
        triggers.insert(make_pair(make_hash(trigger_name), matches[1]));

        // Found Trigger of the interest
        //
        HLT foundHLT;

        foundHLT.hash = make_hash(trigger_name);
        foundHLT.id = distance(triggerNames.begin(), trigger);
        foundHLT.version = matches[2].matched
            ? lexical_cast<int>(matches[2])
            : 0;

        _hlts[*trigger] = foundHLT;
    }

    if (_hlts.empty())
        LogWarning("S8TreeMaker")
            << "None of the searched HLT Triggers is found" << endl;

    _didInitializeHltConfigProvider = true;
}

void S8TreeMaker::analyze(const edm::Event &event,
                          const edm::EventSetup &eventSetup)
{
    if (!_tree)
        throw cms::Exception("S8TreeMaker")
            << "Tree does not exist";

    processEventID(event);
    processGenEvent(event);
    processElectrons(event);
    processJets(event);
    processMuons(event);
    processPrimaryVertices(event);

    if (_saveTriggers)
        processTriggers(event, eventSetup);

    // Write Tree entry
    //
    _tree->Fill();

    // Reset event
    //
    _eventID->reset();
    _genEvent->reset();

    for(s8::Jets::iterator jet = _s8Jets->begin();
        _s8Jets->end() != jet;
        ++jet)
    {
        delete *jet;
    }
    _s8Jets->clear();

    for(s8::Leptons::iterator electron = _s8Electrons->begin();
        _s8Electrons->end() != electron;
        ++electron)
    {
        delete *electron;
    }
    _s8Electrons->clear();

    for(s8::Leptons::iterator muon = _s8Muons->begin();
        _s8Muons->end() != muon;
        ++muon)
    {
        delete *muon;
    }
    _s8Muons->clear();

    for(s8::PrimaryVertices::iterator primaryVertex = _s8PrimaryVertices->begin();
        _s8PrimaryVertices->end() != primaryVertex;
        ++primaryVertex)
    {
        delete *primaryVertex;
    }
    _s8PrimaryVertices->clear();

    for(s8::Triggers::iterator trigger = _s8Triggers->begin();
        _s8Triggers->end() != trigger;
        ++trigger)
    {
        delete *trigger;
    }
    _s8Triggers->clear();
}

void S8TreeMaker::processEventID(const edm::Event &event)
{
    _eventID->setRun(event.id().run());
    _eventID->setLumiBlock(event.id().luminosityBlock());
    _eventID->setEvent(event.id().event());
}

void S8TreeMaker::processGenEvent(const edm::Event &event)
{
    using reco::GenParticleCollection;

    if (event.isRealData())
        return;

    // check if Event is Pythia
    //
    if (_isPythia)
    {
        Handle<GenEventInfoProduct> generator;
        event.getByLabel(InputTag("generator"), generator);

        if (!generator.isValid())
        {
            LogWarning("S8TreeMaker")
                << "failed to extract Generator";

            return;
        }

        _genEvent->setPtHat(generator->qScale());
    }

    Handle<GenParticleCollection> genParticles;
    event.getByLabel(InputTag("genParticles"), genParticles);

    boost::tribool bsplit = boost::indeterminate;
    boost::tribool csplit = boost::indeterminate;
    for(GenParticleCollection::const_iterator particle = genParticles->begin();
        genParticles->end() != particle;
        ++particle)
    {
        const int pdgID = abs(particle->pdgId());

        switch(pdgID)
        {
            case 5: fillSplit(particle->status(), bsplit);
                    break;

            case 4: fillSplit(particle->status(), csplit);
                    break;

            default: break;
        }
    }

    if (bsplit)
        _genEvent->setGluonSplitting(s8::GenEvent::BB, true);

    if (csplit)
        _genEvent->setGluonSplitting(s8::GenEvent::CC, true);
}

void S8TreeMaker::processElectrons(const edm::Event &event)
{
    using pat::ElectronCollection;

    // Extract Electrons
    //
    Handle<ElectronCollection> electrons;
    event.getByLabel(InputTag(_electrons), electrons);

    if (!electrons.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Electrons.";

        return;
    }

    // Process all Electrons
    //
    for(ElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        s8::Lepton *s8Electron = new s8::Lepton();;

        setP4(s8Electron->p4(), electron->p4());
        setVertex(s8Electron->vertex(), electron->vertex());

        s8Electron->impactParameter()->first = electron->dB();
        s8Electron->impactParameter()->second = electron->edB();

        // Extract GenParticle information
        //
        if (electron->genLepton())
        {
            s8::GenParticle *s8GenParticle = s8Electron->genParticle();

            setP4(s8GenParticle->p4(), electron->genLepton()->p4());
            setVertex(s8GenParticle->vertex(), electron->genLepton()->vertex());

            s8GenParticle->setId(electron->genLepton()->pdgId());
            s8GenParticle->setStatus(electron->genLepton()->status());
            if (electron->genLepton()->mother())
                s8GenParticle->setParentId(electron->genLepton()->mother()->pdgId());
        }

        _s8Electrons->push_back(s8Electron);
    } // End loop over electrons
}

void S8TreeMaker::processJets(const edm::Event &event)
{
    using pat::JetCollection;

    // Extract Jets
    //
    Handle<JetCollection> jets;
    event.getByLabel(InputTag(_jets), jets);

    if (!jets.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Jets.";

        return;
    }

    // Process All Jets
    //
    pat::strbitset jetBitset = _jetSelector.getBitTemplate();
    for(JetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        if (!_jetSelector(*jet, jetBitset))
            continue;

        using s8::Jet;
        Jet *s8Jet = new Jet();

        setP4(s8Jet->p4(), jet->p4());

        s8Jet->setFlavour(jet->partonFlavour());
        s8Jet->setTracks(jet->associatedTracks().size());

        if (jet->genParton())
        {
            s8::GenParticle *s8GenParticle = s8Jet->genParticle();

            setP4(s8GenParticle->p4(), jet->genParton()->p4());
            setVertex(s8GenParticle->vertex(), jet->genParton()->vertex());

            s8GenParticle->setId(jet->genParton()->pdgId());
            s8GenParticle->setStatus(jet->genParton()->status());
            if (jet->genParton()->mother())
                s8GenParticle->setParentId(jet->genParton()->mother()->pdgId());
        }

        // Save b-taggers
        //
        s8Jet->setBTag(Jet::TCHE,
                       jet->bDiscriminator("trackCountingHighEffBJetTags"));

        s8Jet->setBTag(Jet::TCHP,
                       jet->bDiscriminator("trackCountingHighPurBJetTags"));

        s8Jet->setBTag(Jet::JP,
                       jet->bDiscriminator("jetProbabilityBJetTags"));

        s8Jet->setBTag(Jet::SSV,
                       jet->bDiscriminator("simpleSecondaryVertexBJetTags"));

        s8Jet->setBTag(Jet::SSVHE,
                       jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));

        s8Jet->setBTag(Jet::SSVHP,
                       jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));

        s8Jet->setBTag(Jet::CSV,
                       jet->bDiscriminator("combinedSecondaryVertexBJetTags"));

        s8Jet->setBTag(Jet::CSV_MVA,
                       jet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));

        s8Jet->setBTag(Jet::JBP,
                       jet->bDiscriminator("jetBProbabilityBJetTags"));

        _s8Jets->push_back(s8Jet);
    }
}

void S8TreeMaker::processMuons(const edm::Event &event)
{
    using pat::MuonCollection;

    // Extract Muons
    //
    Handle<MuonCollection> muons;
    event.getByLabel(InputTag(_muons), muons);

    if (!muons.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Muons.";

        return;
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

        s8::Lepton *s8Muon = new s8::Lepton();

        setP4(s8Muon->p4(), muon->p4());
        setVertex(s8Muon->vertex(), muon->vertex());

        s8Muon->impactParameter()->first = muon->dB();
        s8Muon->impactParameter()->second = muon->edB();

        // Extract GenParticle information
        //
        if (muon->genLepton())
        {
            s8::GenParticle *s8GenParticle = s8Muon->genParticle();

            setP4(s8GenParticle->p4(), muon->genLepton()->p4());
            setVertex(s8GenParticle->vertex(), muon->genLepton()->vertex());

            s8GenParticle->setId(muon->genLepton()->pdgId());
            s8GenParticle->setStatus(muon->genLepton()->status());
            if (muon->genLepton()->mother())
                s8GenParticle->setParentId(muon->genLepton()->mother()->pdgId());
        }

        _s8Muons->push_back(s8Muon);
    } // End loop over muons
}

void S8TreeMaker::processPrimaryVertices(const edm::Event &event)
{
    // Primary Vertices
    //
    typedef vector<Vertex> PVCollection;

    Handle<PVCollection> primaryVertices;
    event.getByLabel(InputTag(_primaryVertices), primaryVertices);

    if (!primaryVertices.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Primary Vertices.";

        return;
    }

    if (primaryVertices->empty())
    {
        LogWarning("S8TreeMaker")
            << "primary vertices collection is empty.";

        return;
    }

    // Process Primary Vertices
    //
    for(PVCollection::const_iterator vertex = primaryVertices->begin();
        primaryVertices->end() != vertex;
        ++vertex)
    {
        if (!isGoodPrimaryVertex(*vertex, event.isRealData()))
            continue;

        s8::PrimaryVertex *s8Vertex = new s8::PrimaryVertex();

        setVertex(s8Vertex->vertex(), vertex->position());
        s8Vertex->setNdof(vertex->ndof());
        s8Vertex->setRho(vertex->position().Rho());

        _s8PrimaryVertices->push_back(s8Vertex);
    }
}

void S8TreeMaker::processTriggers(const edm::Event &event,
                                  const edm::EventSetup &eventSetup)
{
    if (_hlts.empty())
        return;

    // Triggers
    //
    Handle<TriggerResults> triggers;
    event.getByLabel(InputTag(_triggers), triggers);

    if (!triggers.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Triggers";

        return;
    }

    // Process only found HLTs
    //
    for(HLTs::const_iterator hlt = _hlts.begin();
        _hlts.end() != hlt;
        ++hlt)
    {
        s8::Trigger *s8Trigger = new s8::Trigger();
        s8Trigger->setHash(hlt->second.hash);
        if (hlt->second.version)
            s8Trigger->setVersion(hlt->second.version);

        s8Trigger->setIsPass(triggers->accept(hlt->second.id));
        s8Trigger->setPrescale(_hltConfigProvider.prescaleValue(event,
            eventSetup, hlt->first));


        _s8Triggers->push_back(s8Trigger);
    }
}

bool S8TreeMaker::isGoodPrimaryVertex(const Vertex &vertex,
                                    const bool &isData)
{
    return !vertex.isFake() &&
            4 <= vertex.ndof() &&
            (isData ? 24 : 15) >= fabs(vertex.z()) &&
            2 >= fabs(vertex.position().Rho());
}

DEFINE_FWK_MODULE(S8TreeMaker);

/**
 * S8TreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TREEMAKER
#define S8_TREEMAKER

#include <memory>
#include <string>

#include "TFile.h"

#include "Tree/System8/interface/S8Tools.h"
#include "Tree/System8/interface/S8Fwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

class TTree;

namespace reco
{
    class Vertex;
}

namespace s8
{
    class TreeInfo;
    class TriggerCenter;
}

class S8TreeMaker : public edm::EDAnalyzer
{
    /*
     * Produce S8 ROOT Tree
     */
    public:
        S8TreeMaker(const edm::ParameterSet &);
        virtual ~S8TreeMaker();

    private:
        virtual void beginJob();
        virtual void beginRun(const edm::Run &, const edm::EventSetup &);
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        void processEventID(const edm::Event &);
        void processGenEvent(const edm::Event &);
        void processElectrons(const edm::Event &);
        void processJets(const edm::Event &);
        void processMuons(const edm::Event &);
        void processPrimaryVertices(const edm::Event &);
        void processTriggers(const edm::Event &, const edm::EventSetup &);

        bool isGoodPrimaryVertex(const reco::Vertex &, const bool & = false); 

        std::auto_ptr<s8::EventID>         _eventID;
        std::auto_ptr<s8::GenEvent>        _genEvent;
        std::auto_ptr<s8::Jets>            _s8Jets;
        std::auto_ptr<s8::Leptons>         _s8Electrons;
        std::auto_ptr<s8::Leptons>         _s8Muons;
        std::auto_ptr<s8::PrimaryVertices> _s8PrimaryVertices;
        std::auto_ptr<s8::Triggers>        _s8Triggers;

        std::auto_ptr<s8::TreeInfo>       _treeInfo;
        std::auto_ptr<s8::TriggerCenter>  _triggerCenter;

        TTree                            *_tree;
        PFJetIDSelectionFunctor           _jetSelector;
        HLTConfigProvider                 _hltConfigProvider;

        struct HLT
        {
            s8::tools::Hash hash;
            int             id;
            int             version;
        };

        typedef std::map<std::string, HLT>  HLTs;

        HLTs _hlts;

        std::string _primaryVertices;
        std::string _jets;
        std::string _muons;
        std::string _electrons;
        std::string _triggers;
        bool        _isPythia;
        bool        _didInitializeHltConfigProvider;
        bool        _saveTriggers;
};

#endif

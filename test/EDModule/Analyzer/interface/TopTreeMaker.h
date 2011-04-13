/**
 * TopTreeMaker
 * 
 *
 * Created by Samvel Khalatian on August 20, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_TREEMAKER
#define TOP_TREEMAKER

#include <string>

#include "Tree/Top/interface/TopEvent.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

class TopTreeMaker : public edm::EDAnalyzer
{
    /*
     * Produce Top ROOT Tree
     */
    public:
        TopTreeMaker(const edm::ParameterSet &);
        virtual ~TopTreeMaker();

    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        std::auto_ptr<top::Event>  _event;
        TTree                     *_tree;

        std::string _beamSpots;
        std::string _electrons;
        std::string _genParticles;
        std::string _jets;
        std::string _mets;
        std::string _muons;
};

#endif

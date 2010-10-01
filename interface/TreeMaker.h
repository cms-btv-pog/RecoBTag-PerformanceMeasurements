/**
 * TreeMaker
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

#include "RecoBTag/PerformanceMeasurements/interface/Event.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

class TTree;

namespace reco
{
    class Vertex;
}

class TreeMaker : public edm::EDAnalyzer
{
    /*
     * Produce S8 ROOT Tree
     */
    public:
    TreeMaker(const edm::ParameterSet &);
    virtual ~TreeMaker();

    private:
    virtual void beginJob();
    virtual void analyze(const edm::Event &, const edm::EventSetup &);
    virtual void endJob();

    bool isGoodPrimaryVertex(const reco::Vertex &, const bool & = false); 

    std::auto_ptr<s8::Event>  _event;
    TTree                    *_tree;

    std::string _primaryVertices;
    std::string _jets;
    std::string _muons;
};

#endif

/**
 * PVFilter
 * top
 *
 * Created by Samvel Khalatian on August 19, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_TREE_PVFILTER
#define TOP_TREE_PVFILTER

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

namespace top
{
    /*
     * Filter out events passing Primary Vertex cuts.
     *
     * Input: RECO
     * Doc  : Muon Synchronization Exercise (ttmuj group)
     */
    class PVFilter : public edm::EDFilter
    {
        public:
            PVFilter(const edm::ParameterSet &);
            virtual ~PVFilter();

        private:
            virtual bool filter(edm::Event &, const edm::EventSetup &);

            std::string _pvTag;
            bool        _isDataInput;
    };
}

#endif

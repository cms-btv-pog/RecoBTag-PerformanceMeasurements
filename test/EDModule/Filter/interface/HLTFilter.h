/**
 * HLTFilter
 * top
 *
 * Created by Samvel Khalatian on August 19, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_TREE_HLTFILTER
#define TOP_TREE_HLTFILTER

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

namespace top
{
    /*
     * Filter out events passing specific High Level Trigger.
     *
     * Input: RECO
     * Doc  : Muon Synchronization Exercise (ttmuj group)
     */
    class HLTFilter : public edm::EDFilter
    {
        public:
            HLTFilter(const edm::ParameterSet &);
            virtual ~HLTFilter();

        private:
            virtual bool filter(edm::Event &, const edm::EventSetup &);

            std::string _tag;
            std::string _hlt;
    };
}

#endif

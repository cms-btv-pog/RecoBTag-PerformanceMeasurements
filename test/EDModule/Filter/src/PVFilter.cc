/**
 * PVFilter
 * top
 *
 * Created by Samvel Khalatian on August 19, 2010
 * Copyright 2010, All rights reserved
 */

#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "EDModule/Filter/interface/PVFilter.h"

using top::PVFilter;

PVFilter::PVFilter(const edm::ParameterSet &parameters)
{
    using std::string;

    _pvTag  = parameters.getParameter<string>("pvTag");
    _isDataInput = parameters.getParameter<bool>("isDataInput");
}

PVFilter::~PVFilter()
{
}

bool PVFilter::filter(edm::Event &event, const edm::EventSetup &)
{
    using edm::Handle;
    using reco::Vertex;

    typedef std::vector<Vertex> PVCollection;

    // Extract Trigger Restuls
    Handle<PVCollection> primaryVertices;
    event.getByLabel(edm::InputTag(_pvTag), primaryVertices);

    if (!primaryVertices.isValid() ||
        !primaryVertices->size())
        return false;

    // Use the first available Primary Vertex
    const Vertex &primaryVertex = primaryVertices->at(0);

    return !primaryVertex.isFake() &&
            4 <= primaryVertex.ndof() &&
            (_isDataInput ? 24 : 15) >= fabs(primaryVertex.z()) &&
            2 >= fabs(primaryVertex.position().Rho());
}

DEFINE_FWK_MODULE(PVFilter);

/**
 * PrimaryVertex
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include "RecoBTag/PerformanceMeasurements/interface/PrimaryVertex.h"

using std::runtime_error;

using s8::PrimaryVertex;

void PrimaryVertex::setNdof(const int &ndof)
{
    if (0 > ndof)
        throw runtime_error("[PrimaryVertex] Wrong number of NDOF supplied.");

    _ndof = ndof;
}

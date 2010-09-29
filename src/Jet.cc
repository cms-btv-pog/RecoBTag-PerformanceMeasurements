/**
 * Jet
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include "RecoBTag/PerformanceMeasurements/interface/Jet.h"

using std::runtime_error;

using s8::Jet;

Jet::Jet() throw():
    _flavour(0),
    _tracks(0),
    _btag(Jet::BTAGS)
{
}

void Jet::setTracks(const int &tracks)
{
    if (0 > tracks)
        throw runtime_error("[Jet] Wrong number of tracks supplied.");

    _tracks = tracks;
}

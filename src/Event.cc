/**
 * Event
 * s8
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#include "RecoBTag/PerformanceMeasurements/interface/EventID.h"

#include "RecoBTag/PerformanceMeasurements/interface/Event.h"

using s8::Event;

void Event::reset()
{
    _id = EventID();

    _jets.clear();
    _muons.clear();
    _primaryVertices.clear();
}

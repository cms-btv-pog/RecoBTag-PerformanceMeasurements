/**
 * GenEvent
 * s8 
 *
 * Created by Samvel Khalatian on Oct 11, 2010
 * Copyright 2010, All rights reserved
 */

#include <algorithm>
#include <stdexcept>

#include "interface/S8GenEvent.h"

using std::find;
using std::runtime_error;

using s8::GenEvent;

GenEvent::GenEvent() throw():
    _ptHat(0)
{
    for(int gluonSplitting = 0;
        GLUON_SPLITTINGS > gluonSplitting;
        ++gluonSplitting)
    {
        *(_gluonSplittings + gluonSplitting) = false;
    }
}

GenEvent::GenEvent(const GenEvent &event):
    _ptHat(event.ptHat())
{
    for(int gluonSplitting = 0;
        GLUON_SPLITTINGS > gluonSplitting;
        ++gluonSplitting)
    {
        *(_gluonSplittings + gluonSplitting) =
            event.isGluonSplitting(GluonSplitting(gluonSplitting));
    }
}

GenEvent &GenEvent::operator =(const GenEvent &event)
{
    _ptHat = event.ptHat();

    for(int gluonSplitting = 0;
        GLUON_SPLITTINGS > gluonSplitting;
        ++gluonSplitting)
    {
        *(_gluonSplittings + gluonSplitting) =
            event.isGluonSplitting(GluonSplitting(gluonSplitting));
    }

    return *this;
}

void GenEvent::reset()
{
    _ptHat = 0;

    for(int gluonSplitting = 0;
        GLUON_SPLITTINGS > gluonSplitting;
        ++gluonSplitting)
    {
        *(_gluonSplittings + gluonSplitting) = false;
    }
}

bool GenEvent::isGluonSplitting() const
{
    return isGluonSplitting(BB) ||
           isGluonSplitting(CC);
}

bool GenEvent::isGluonSplitting(const GluonSplitting &gluonSplitting) const
{
    if (GLUON_SPLITTINGS <= gluonSplitting)
        throw runtime_error("[GenEvent] Unsupported GluonSplitting value used");

    return _gluonSplittings[gluonSplitting];
}

double GenEvent::ptHat() const
{
    return _ptHat;
}

void GenEvent::setGluonSplitting(const GluonSplitting &gluonSplitting,
                                 const bool &value)
{
    switch(gluonSplitting)
    {
        case BB:   // Fall through
        case CC:   _gluonSplittings[gluonSplitting] = value; 
                   break;

        default:
            throw runtime_error("[GenEvent] Unsupported Gluon Splitting value supplied");
    }
}

void GenEvent::setPtHat(const double &ptHat)
{
    if (ptHat < 0)
        throw runtime_error("[GenEvent] Negative PtHat supplied");

    _ptHat = ptHat;
}

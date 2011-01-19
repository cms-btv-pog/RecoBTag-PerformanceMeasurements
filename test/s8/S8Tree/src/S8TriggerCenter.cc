/**
 * TriggerCenter
 * s8 
 *
 * Created by Samvel Khalatian on Nov 5, 2010
 * Copyright 2010, All rights reserved
 */

#include "interface/S8TriggerCenter.h"

using s8::TriggerCenter;

TriggerCenter::TriggerCenter()
{
}

TriggerCenter::TriggerMap &TriggerCenter::triggers()
{
    return _triggers;
}

const TriggerCenter::TriggerMap &TriggerCenter::triggers() const
{
    return _triggers;
}

void TriggerCenter::merge(const TriggerCenter &triggerCenter)
{
    for(TriggerMap::const_iterator trigger = triggerCenter._triggers.begin();
        triggerCenter._triggers.end() != trigger;
        ++trigger)
    {
        // Map will automatically take care of elements with the same hash
        //
        _triggers.insert(*trigger);
    }
}

ClassImp(TriggerCenter)

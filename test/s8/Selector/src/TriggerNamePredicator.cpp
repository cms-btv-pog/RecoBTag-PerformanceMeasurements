/**
 * TriggerNamePredicator
 * s8
 *
 * Created by Samvel Khalatian on Feb 8, 2011
 * Copyright 2010, All rights reserved
 */

#include "S8Tree/interface/S8Trigger.h"
#include "Selector/interface/TriggerNamePredicator.h"

using s8::TriggerNamePredicator;

TriggerNamePredicator::TriggerNamePredicator() throw()
{
    _search_trigger = 0;
}

void TriggerNamePredicator::setSearchTrigger(const Trigger &trigger)
{
    _search_trigger = &trigger;
}

bool TriggerNamePredicator::operator()(const Trigger *trigger) const
{
    if (!_search_trigger)
        return false;

    return _search_trigger->hash() == trigger->hash();
}

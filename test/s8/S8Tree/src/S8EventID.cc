/**
 * EventID
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include "interface/S8EventID.h"

using std::runtime_error;

using s8::EventID;

EventID::EventID() throw():
    _run(0),
    _lumiBlock(0),
    _event(0)
{
}

EventID::~EventID() throw()
{
}

void EventID::reset()
{
    _run = 0;
    _lumiBlock = 0;
    _event = 0;
}

int EventID::run() const
{
    return _run;
}

int EventID::lumiBlock() const
{
    return _lumiBlock;
}

int EventID::event() const
{
    return _event;
}

void EventID::setRun(const int &run)
{
    if (run < 0)
        throw runtime_error("[EventID] Bad run number supplied.");

    _run = run;
}

void EventID::setLumiBlock(const int &lumiBlock)
{
    if (lumiBlock < 0)
        throw runtime_error("[EventID] Bad LumiBlock number supplied.");

    _lumiBlock = lumiBlock;
}

void EventID::setEvent(const int &event)
{
    if (event < 0)
        throw runtime_error("[EventID] Bad Event number supplied.");

    _event = event;
}

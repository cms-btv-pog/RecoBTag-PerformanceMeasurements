/**
 * Event
 * s8
 *
 * Created by Samvel Khalatian on Nov 21, 2010
 * Copyright 2010, All rights reserved
 */

#include "S8Tree/interface/S8EventID.h"
#include "S8Tree/interface/S8GenEvent.h"

#include "IO/interface/Event.h"

using s8::Event;

Event::Event() throw()
{
    _id = 0;
    _gen = 0;

    _jets = 0;
    _electrons = 0;
    _muons = 0;
    _primaryVertices = 0;
    _triggers = 0;
}

const s8::EventID *Event::id() const
{
    return _id;
}

const s8::GenEvent *Event::gen() const
{
    return _gen;
}

const s8::Jets *Event::jets() const
{
    return _jets;
}

const s8::Leptons *Event::electrons() const
{
    return _electrons;
}

const s8::Leptons *Event::muons() const
{
    return _muons;
}

const s8::PrimaryVertices *Event::primaryVertices() const
{
    return _primaryVertices;
}

const s8::Triggers *Event::triggers() const
{
    return _triggers;
}



// Setters
//
void Event::setID(const EventID *id)
{
    _id = id;
}

void Event::setGen(const GenEvent *gen)
{
    _gen = gen;
}

void Event::setJets(const Jets *jets)
{
    _jets = jets;
}

void Event::setElectrons(const Leptons *electrons)
{
    _electrons = electrons;
}

void Event::setMuons(const Leptons *muons)
{
    _muons = muons;
}

void Event::setPrimaryVertices(const PrimaryVertices *primaryVertices)
{
    _primaryVertices = primaryVertices;
}

void Event::setTriggers(const Triggers *triggers)
{
    _triggers = triggers;
}

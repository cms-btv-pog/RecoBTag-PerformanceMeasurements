/**
 * GenParticle
 * s8 
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>
#include <TVector3.h>

#include "interface/S8GenParticle.h"

using s8::GenParticle;

GenParticle::GenParticle() throw():
    _id(0),
    _parentId(0),
    _status(0)
{
    _p4 = new TLorentzVector();
    _vertex = new TVector3();
}

GenParticle::~GenParticle() throw()
{
    delete _vertex;
    delete _p4;
}

GenParticle::GenParticle(const GenParticle &particle):
    _id(particle.id()),
    _parentId(particle.parentId()),
    _status(particle.status())
{
    _p4 = new TLorentzVector(*particle.p4());
    _vertex = new TVector3(*particle.vertex());
}

GenParticle &GenParticle::operator =(const GenParticle &particle)
{
    _id = particle.id();
    _parentId = particle.parentId();
    _status = particle.status();

    *_p4 = *particle.p4();
    *_vertex = *particle.vertex();

    return *this;
}

void GenParticle::reset()
{
    _id = 0;
    _parentId = 0;
    _status = 0;

    _p4->SetPxPyPzE(0, 0, 0, 0);
    _vertex->SetXYZ(0, 0, 0);
}

int GenParticle::id() const
{
    return _id;
}

int GenParticle::parentId() const
{
    return _parentId;
}

int GenParticle::status() const
{
    return _status;
}

TLorentzVector *GenParticle::p4()
{
    return _p4;
}

const TLorentzVector *GenParticle::p4() const
{
    return _p4;
}

TVector3 *GenParticle::vertex()
{
    return _vertex;
}

const TVector3 *GenParticle::vertex() const
{
    return _vertex;
}




void GenParticle::setId(const int &id)
{
    _id = id;
}

void GenParticle::setParentId(const int &id)
{
    _parentId = id;
}

void GenParticle::setStatus(const int &status)
{
    _status = status;
}

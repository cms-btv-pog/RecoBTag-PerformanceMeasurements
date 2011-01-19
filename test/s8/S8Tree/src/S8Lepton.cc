/**
 * Lepton
 * s8 
 *
 * Created by Samvel Khalatian on Nov 15, 2010
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>
#include <TVector3.h>

#include "interface/S8GenParticle.h"

#include "interface/S8Lepton.h"

using s8::Lepton;

Lepton::Lepton() throw()
{
    _impactParameter = new ImpactParameter();
    _p4 = new TLorentzVector();
    _vertex = new TVector3();
    _genParticle = new GenParticle();
}

Lepton::~Lepton() throw()
{
    delete _genParticle;
    delete _vertex;
    delete _p4;
    delete _impactParameter;
}

Lepton::Lepton(const Lepton &lepton)
{
    _impactParameter = new ImpactParameter(*lepton.impactParameter());
    _p4 = new TLorentzVector(*lepton.p4());
    _vertex = new TVector3(*lepton.vertex());
    _genParticle = new GenParticle(*lepton.genParticle());
}

Lepton &Lepton::operator =(const Lepton &lepton)
{
    *_impactParameter = *lepton.impactParameter();
    *_p4 = *lepton.p4();
    *_vertex = *lepton.vertex();
    *_genParticle = *lepton.genParticle();

    return *this;
}

void Lepton::reset()
{
    _impactParameter->first = 0;
    _impactParameter->second = 0;

    _p4->SetPxPyPzE(0, 0, 0, 0);
    _vertex->SetXYZ(0, 0, 0);
    _genParticle->reset();
}

Lepton::ImpactParameter *Lepton::impactParameter()
{
    return _impactParameter;
}

const Lepton::ImpactParameter *Lepton::impactParameter() const
{
    return _impactParameter;
}

TLorentzVector *Lepton::p4()
{
    return _p4;
}

const TLorentzVector *Lepton::p4() const
{
    return _p4;
}

TVector3 *Lepton::vertex()
{
    return _vertex;
}

const TVector3 *Lepton::vertex() const
{
    return _vertex;
}

s8::GenParticle *Lepton::genParticle()
{
    return _genParticle;
}

const s8::GenParticle *Lepton::genParticle() const
{
    return _genParticle;
}

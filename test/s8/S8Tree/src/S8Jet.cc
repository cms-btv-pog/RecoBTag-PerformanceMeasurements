/**
 * Jet
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include <TLorentzVector.h>

#include "interface/S8GenParticle.h"

#include "interface/S8Jet.h"

using std::runtime_error;
using std::out_of_range;

using s8::Jet;

Jet::Jet() throw():
    _flavour(0),
    _tracks(0)
{
    _p4 = new TLorentzVector();
    _genParticle = new GenParticle();

    for(int btag = 0; BTAGS > btag; ++btag)
    {
        *(_btag + btag) = 0;
    }
}

Jet::~Jet() throw()
{
    delete _genParticle;
    delete _p4;
}

Jet::Jet(const Jet &jet):
    _flavour(jet.flavour()),
    _tracks(jet.tracks())
{
    _p4 = new TLorentzVector(*jet.p4());
    _genParticle = new GenParticle(*jet.genParticle());

    for(int btag = 0; BTAGS > btag; ++btag)
    {
        *(_btag + btag) = *(jet._btag + btag);
    }
}

Jet &Jet::operator =(const Jet &jet)
{
    _flavour = jet.flavour();
    _tracks = jet.tracks();

    *_p4 = *jet.p4();
    *_genParticle = *jet.genParticle();

    for(int btag = 0; BTAGS > btag; ++btag)
    {
        *(_btag + btag) = *(jet._btag + btag);
    }

    return *this;
}

void Jet::reset()
{
    _flavour = 0;
    _tracks = 0;

    _p4->SetPxPyPzE(0, 0, 0, 0);
    _genParticle->reset();

    for(int btag = 0; BTAGS > btag; ++btag)
    {
        *(_btag + btag) = 0;
    }
}

int Jet::flavour() const
{
    return _flavour;
}

int Jet::tracks() const
{
    return _tracks;
}

TLorentzVector *Jet::p4()
{
    return _p4;
}

const TLorentzVector *Jet::p4() const
{
    return _p4;
}

s8::GenParticle *Jet::genParticle()
{
    return _genParticle;
}

const s8::GenParticle *Jet::genParticle() const
{
    return _genParticle;
}

double Jet::btag(const BTag &tag) const
{
    if (TCHE > tag || SSVHP < tag)
        throw out_of_range("BTag id is out of range");

    return *(_btag + tag);
}



void Jet::setFlavour(const int &flavour)
{
    _flavour = flavour;
}

void Jet::setTracks(const int &tracks)
{
    if (0 > tracks)
        throw runtime_error("[Jet] Wrong number of tracks supplied.");

    _tracks = tracks;
}

void Jet::setBTag(const BTag &tag, const double &value)
{
    if (TCHE > tag || SSVHP < tag)
        throw out_of_range("BTag id is out of range");

    *(_btag + tag) = value;
}

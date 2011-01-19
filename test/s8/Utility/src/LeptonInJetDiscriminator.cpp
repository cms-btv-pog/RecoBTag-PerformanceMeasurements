/**
 * LeptonInJetDiscriminator
 * s8 
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>

#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"

#include "Utility/interface/LeptonInJetDiscriminator.h"

using s8::Jet;
using s8::Lepton;

using s8::LeptonInJetDiscriminator;
using s8::LeptonInJetPtRelDiscriminator;
using s8::LeptonInJetPtDiscriminator;

LeptonInJetDiscriminator::~LeptonInJetDiscriminator()
{
}



LeptonInJetPtRelDiscriminator::LeptonInJetPtRelDiscriminator()
{
    _binGroup.bins = 700;
    _binGroup.min = 0;
    _binGroup.max = 7;
}

LeptonInJetPtRelDiscriminator::~LeptonInJetPtRelDiscriminator()
{
}

double LeptonInJetPtRelDiscriminator::operator()(const Lepton *lepton,
                                                 const Jet *jet)
{
    return lepton->p4()->Vect().Perp(jet->p4()->Vect());
}

s8::BinGroup LeptonInJetPtRelDiscriminator::bins() const
{
    return _binGroup;
}



LeptonInJetPtDiscriminator::LeptonInJetPtDiscriminator()
{
    _binGroup.bins = 510;
    _binGroup.min = -1;
    _binGroup.max = 50;
}

LeptonInJetPtDiscriminator::~LeptonInJetPtDiscriminator()
{
}

double LeptonInJetPtDiscriminator::operator()(const Lepton *lepton,
                                              const Jet *jet)
{
    return lepton ? lepton->p4()->Pt() : -1;
}

s8::BinGroup LeptonInJetPtDiscriminator::bins() const
{
    return _binGroup;
}

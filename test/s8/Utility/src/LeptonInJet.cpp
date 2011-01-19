/**
 * LeptonInJet
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "S8Tree/interface/S8Lepton.h"
#include "S8Tree/interface/S8PrimaryVertex.h"
#include "Utility/interface/LeptonInJetDelegate.h"

#include "Utility/interface/LeptonInJet.h"

using s8::LeptonInJet;
using s8::LeptonInJetDelegate;

LeptonInJet::LeptonInJet() throw()
{
    _delegate = 0;
}

LeptonInJet::~LeptonInJet() throw()
{
}

LeptonInJetDelegate *LeptonInJet::delegate() const
{
    return _delegate;
}

void LeptonInJet::setDelegate(LeptonInJetDelegate *delegate)
{
    _delegate = delegate;
}

void LeptonInJet::operator()(const PrimaryVertex *primaryVertex,
                             const Leptons *leptons,
                             const Jet *jet)
{
    // There is nothing to do if Delegate is not set
    //
    if (!_delegate)
        return;

    for(Leptons::const_iterator lepton = leptons->begin();
        leptons->end() != lepton;
        ++lepton)
    {
        const Lepton *particle = *lepton;

        if (1 <= fabs(particle->vertex()->z() - primaryVertex->vertex()->z()))
            continue;

        const double deltaR = particle->p4()->DeltaR(*jet->p4());
        if (0.01 > deltaR ||
            .4 <= deltaR ||
            -1 >= particle->p4()->Vect().Perp(jet->p4()->Vect()))

            continue;

        // found Lepton in Jet
        //
        _delegate->leptonIsInJet(particle, jet);

        if (_delegate->shouldLookForMoreLeptons())
            continue;

        break;
    }
}

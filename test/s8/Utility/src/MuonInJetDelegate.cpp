/**
 * MuonInJetDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#include "Utility/interface/MuonInJetDelegate.h"

using s8::MuonInJetDelegate;

MuonInJetDelegate::~MuonInJetDelegate() throw()
{
}

bool MuonInJetDelegate::muonInJetShouldProcessJet(const Jet *)
{
    return true;
}

bool MuonInJetDelegate::shouldSkipMuonInJetPlusAwayJet(const Lepton *,
                                                       const Jet *)
{
    return false;
}

void MuonInJetDelegate::muonIsInJetPlusAwayJet(const Lepton *, const Jet *)
{
}

void MuonInJetDelegate::muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                     const Jet *)
{
}

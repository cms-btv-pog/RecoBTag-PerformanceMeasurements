/**
 * TwoMuonsInJetDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#include "Utility/interface/TwoMuonsInJetDelegate.h"

using s8::TwoMuonsInJetDelegate;

TwoMuonsInJetDelegate::~TwoMuonsInJetDelegate() throw()
{
}

bool TwoMuonsInJetDelegate::muonInJetShouldProcessJet(const Jet *)
{
    return true;
}

void TwoMuonsInJetDelegate::muonIsInJetPlusAwayJet(const Lepton *,
                                                   const Lepton *,
                                                   const Jet *)
{
}

void TwoMuonsInJetDelegate::muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                         const Lepton *,
                                                          const Jet *)
{
}

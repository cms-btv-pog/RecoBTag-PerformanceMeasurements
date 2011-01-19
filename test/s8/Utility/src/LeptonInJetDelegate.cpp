/**
 * LeptonInJetDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#include "Utility/interface/LeptonInJetDelegate.h"

using s8::LeptonInJetDelegate;

LeptonInJetDelegate::~LeptonInJetDelegate() throw()
{
}

void LeptonInJetDelegate::leptonIsInJet(const Lepton *, const Jet *)
{
}

bool LeptonInJetDelegate::shouldLookForMoreLeptons()
{
    return false;
}

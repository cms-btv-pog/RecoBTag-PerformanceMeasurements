/**
 * TriggerOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include "Option/interface/TriggerOptionsDelegate.h"

using s8::TriggerOptionsDelegate;

TriggerOptionsDelegate::TriggerOptionsDelegate() throw()
{
}

TriggerOptionsDelegate::~TriggerOptionsDelegate() throw()
{
}

void TriggerOptionsDelegate::optionTriggerIsSet(const Trigger &)
{
}

void TriggerOptionsDelegate::optionUseTriggerPrescaleIsSet(const bool &)
{
}

void TriggerOptionsDelegate::optionSimulateTriggerIsSet(const bool &)
{
}

void TriggerOptionsDelegate::optionReweightTriggerIsSet(const std::string &)
{
}

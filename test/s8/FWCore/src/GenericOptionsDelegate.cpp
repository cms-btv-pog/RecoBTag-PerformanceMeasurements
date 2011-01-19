/**
 * GenericOptionsDelegate
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include "interface/GenericOptionsDelegate.h"

using std::string;

using core::GenericOptionsDelegate;

GenericOptionsDelegate::~GenericOptionsDelegate() throw()
{
}

void GenericOptionsDelegate::optionDebugIsSet(const string &)
{
}

void GenericOptionsDelegate::optionOutputIsSet(const string &)
{
}

void GenericOptionsDelegate::optionEventsIsSet(const int &)
{
}

void GenericOptionsDelegate::optionInputIsSet(const string &)
{
}

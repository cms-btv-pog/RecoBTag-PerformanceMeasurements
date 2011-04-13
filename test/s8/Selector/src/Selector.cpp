/**
 * Selector
 * s8
 *
 * Created by Samvel Khalatian on Feb 3, 2011
 * Copyright 2010, All rights reserved
 */

#include "Selector/interface/Selector.h"

using s8::Selector;

Selector::Selector() throw()
{
}

Selector::~Selector() throw()
{
}



// Helpers
//
std::ostream &s8::operator<<(std::ostream &out, const Selector &selector)
{
    selector.print(out);

    return out;
}

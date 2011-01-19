/**
 * Options
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <ostream>

#include <boost/program_options.hpp>

#include "interface/Options.h"

using core::Options;

Options::Options() throw()
{
}

Options::~Options() throw()
{
}



std::ostream &core::operator<<(std::ostream &out, const Options &options)
{
    options.print(out);

    return out;
}

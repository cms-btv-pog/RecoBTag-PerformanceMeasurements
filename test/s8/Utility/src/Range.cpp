/**
 * Range
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#include <ostream>

#include "Utility/interface/Range.h"

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using s8::Range;

Range::Range() throw():
    _minimum(0),
    _maximum(0)
{
}

Range::~Range() throw()
{
}

double Range::minimum() const throw()
{
    return _minimum;
}

double Range::maximum() const throw()
{
    return _maximum;
}

void Range::setMinimum(const double &minimum)
{
    _minimum = minimum;
}

void Range::setMaximum(const double &maximum)
{
    _maximum = maximum;
}



// Helpers
//
bool s8::isValueInRange(const double &value, const Range &range)
{
    return (range.minimum() ? value >= range.minimum() : true) &&
        (range.maximum() ? value < range.maximum() : true);
}

void s8::parse(Range &range, const std::string &value)
{
    smatch matches;
    regex pattern("^(?:(\\d+(?:\\.\\d*)?)\\.\\.(\\d+(?:\\.\\d*)?)|(\\d+(?:\\.\\d*)?)\\.\\.|\\.\\.(\\d+(?:\\.\\d*)?))$");
    if (!regex_match(value, matches, pattern))
    {
        range.setMinimum(0);
        range.setMaximum(0);

        return;
    }

    if (matches[1].matched &&
        matches[2].matched)
    {
        range.setMinimum(lexical_cast<double>(matches[1]));
        range.setMaximum(lexical_cast<double>(matches[2]));
    }
    else if(matches[3].matched)
    {
        range.setMinimum(lexical_cast<double>(matches[3]));
    }
    else
    {
        range.setMaximum(lexical_cast<double>(matches[4]));
    }
}

std::ostream &s8::operator<<(std::ostream &out, const Range &range)
{
    if (range.minimum() ||
        range.maximum())
    {
        if (range.minimum())
            out << range.minimum();
        
        out << "..";
        
        if (range.maximum())
            out << range.maximum();
    }
    else
        out << "---";

    return out;
}

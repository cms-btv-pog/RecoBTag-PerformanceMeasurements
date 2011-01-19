/**
 * PythiaOptions
 * s8
 *
 * Created by Samvel Khalatian on Dec 16, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <stdexcept>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "Option/interface/PythiaOptionsDelegate.h"

#include "Option/interface/PythiaOptions.h"

using std::runtime_error;
using std::string;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using s8::PythiaOptions;
using s8::PythiaOptionsDelegate;

PythiaOptions::PythiaOptions() throw():
    _minPtHat(0),
    _maxPtHat(0)
{
    _delegate = 0;
}

PythiaOptions::~PythiaOptions() throw()
{
}

PythiaOptionsDelegate *PythiaOptions::delegate() const
{
    return _delegate;
}

void PythiaOptions::setDelegate(PythiaOptionsDelegate *delegate)
{
    _delegate = delegate;
}

void PythiaOptions::init()
{
    _description.reset(new po::options_description("Pythia Options"));
    _description->add_options()
        ("gluon-splitting",
            po::value<string>()->default_value("keep")->notifier(
                boost::bind(&PythiaOptions::optionGluonSplittingIsSet,
                            this, _1)),
            "Use gluon Splitting (bbbar only): keep, remove, only, enhance, simulate (2 jets events with delta R and pT cut)")

        ("pthat",
            po::value<string>()->default_value("")->notifier(
                boost::bind(&PythiaOptions::optionPtHatIsSet,
                            this, _1)),
            "Cut Pt Hat. Only integers are supported. Format: --pthat=MIN..MAX")
    ;
}

po::options_description *PythiaOptions::description() const
{
    return _description.get();
}

void PythiaOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "Pythia Options" << endl;
    out << setw(25) << left << " [+] Gluon Splitting"
        << _gluonSplitting << endl;

    out << setw(25) << left << " [+] Pt Hat Cut";
    if (_minPtHat || _maxPtHat)
    {
        if (_minPtHat)
            out << _minPtHat;

        out << "..";

        if (_maxPtHat)
            out << _maxPtHat;

        out << endl;
    }
    else
        out << "undefined" << endl;
}



void PythiaOptions::optionGluonSplittingIsSet(const string &value)
{
    string option = value;
    boost::to_lower(option);

    if ("keep" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::KEEP);
    }
    else if ("remove" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::REMOVE);
    }
    else if ("only" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::ONLY);
    }
    else if ("simulate" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::SIMULATE);
    }
    else if ("enhance" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::ENHANCE);
    }
    else
        throw runtime_error("Unsupported Gluon Splitting option used: " + value);

    _gluonSplitting = value;
}

void PythiaOptions::optionPtHatIsSet(const string &value)
{
    if (value.empty())
        return;

    smatch matches;
    regex pattern("^(?:(\\d+)..(\\d+)|(\\d+)..|..(\\d+))$");
    if (!regex_match(value, matches, pattern))
        throw runtime_error("Unsuppored pthat format. Use: min..max");

    if (matches[4].matched)
    {
        _maxPtHat = lexical_cast<int>(matches[4]);
    }
    else if(matches[3].matched)
    {
        _minPtHat = lexical_cast<int>(matches[3]);
    }
    else
    {
        _minPtHat = lexical_cast<int>(matches[1]);
        _maxPtHat = lexical_cast<int>(matches[2]);
    }

    if (_delegate)
        _delegate->optionPtHatIsSet(_minPtHat, _maxPtHat);
}

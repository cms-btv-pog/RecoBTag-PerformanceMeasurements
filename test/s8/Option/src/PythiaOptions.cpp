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
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "Option/interface/PythiaOptionsDelegate.h"

#include "Option/interface/PythiaOptions.h"

using std::runtime_error;
using std::string;

using boost::regex;
using boost::smatch;

using s8::PythiaOptions;
using s8::PythiaOptionsDelegate;

PythiaOptions::PythiaOptions() throw():
    _gluon_splitting("---")
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
        ("g-split",
            po::value<string>()->notifier(
                boost::bind(&PythiaOptions::optionGluonSplittingIsSet,
                            this, _1)),
            "Use gluon splitting. Values:\nBB, CC, ONLY (BB or CC)\nNO_BB, NO_CC, REMOVE (BB and CC)\nADD_BB, ADD_CC, ENHANCE (BB and CC)")

        ("pthat",
            po::value<string>()->notifier(
                boost::bind(&PythiaOptions::optionPtHatIsSet,
                            this, _1)),
            "Cut Pt Hat. Format: --pthat=MIN..MAX")
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
        << _gluon_splitting << endl;
    out << setw(25) << left << " [+] Pt Hat Cut" << _pt_hat << endl;
}



void PythiaOptions::optionGluonSplittingIsSet(const string &value)
{
    string option = value;
    boost::to_lower(option);

    // [see PythiaOptionsDelegate for description]
    //
    if ("bb" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::BB);
    }
    else if ("cc" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::CC);
    }
    else if ("keep" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::KEEP);
    }

    else if ("no_bb" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::NO_BB);
    }
    else if ("no_cc" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::NO_CC);
    }
    else if ("remove" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::REMOVE);
    }

    else if ("add_bb" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::ADD_BB);
    }
    else if ("add_cc" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::ADD_CC);
    }
    else if ("enhance" == option)
    {
        if (_delegate)
            _delegate->optionGluonSplittingIsSet(PythiaOptionsDelegate::ENHANCE);
    }
    else
        throw runtime_error("Unsupported Gluon Splitting option used: " + value);

    _gluon_splitting = value;
}

void PythiaOptions::optionPtHatIsSet(const string &value)
{
    if (value.empty())
        return;

    parse(_pt_hat, value);

    if (_delegate)
        _delegate->optionPtHatIsSet(_pt_hat);
}

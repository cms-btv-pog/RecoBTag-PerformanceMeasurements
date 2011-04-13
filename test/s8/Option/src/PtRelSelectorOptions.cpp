/**
 * PtRelSelectorOptions
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <stdexcept>
#include <string>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include "Option/interface/PtRelSelectorOptionsDelegate.h"

#include "Option/interface/PtRelSelectorOptions.h"

using std::string;
using std::runtime_error;

using s8::PtRelSelectorOptions;
using s8::PtRelSelectorOptionsDelegate;

PtRelSelectorOptions::PtRelSelectorOptions() throw():
    _muonPt(0)
{
    _delegate = 0;
}

PtRelSelectorOptions::~PtRelSelectorOptions() throw()
{
}

PtRelSelectorOptionsDelegate *PtRelSelectorOptions::delegate() const
{
    return _delegate;
}

void PtRelSelectorOptions::setDelegate(PtRelSelectorOptionsDelegate *delegate)
{
    _delegate = delegate;
}

void PtRelSelectorOptions::init()
{
    _description.reset(new po::options_description("PtRel Selector Options"));
    _description->add_options()
        ("muon-pt",
            po::value<double>()->default_value(5.0)->notifier(
                boost::bind(&PtRelSelectorOptions::optionMuonPtIsSet, this, _1)),
            "Muon pT cut")
    ;
}

po::options_description *PtRelSelectorOptions::description() const
{
    return _description.get();
}

void PtRelSelectorOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "PtRel Selector Options" << endl;
    out << setw(25) << left << " [+] Muon pT" << _muonPt << endl;
}



void PtRelSelectorOptions::optionMuonPtIsSet(const double &pt)
{
    if (0 > pt)
        throw runtime_error("Negative muon pT value");

    _muonPt = pt;

    if (_delegate)
        _delegate->optionMuonPtIsSet(pt);
}

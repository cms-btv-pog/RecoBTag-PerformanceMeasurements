/**
 * MuonInJetOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <stdexcept>
#include <string>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include "Option/interface/MuonInJetOptionsDelegate.h"

#include "Option/interface/MuonInJetOptions.h"

using std::string;
using std::runtime_error;

using s8::MuonInJetOptions;
using s8::MuonInJetOptionsDelegate;

MuonInJetOptions::MuonInJetOptions() throw():
    _muonPt(0)
{
    _delegate = 0;
}

MuonInJetOptions::~MuonInJetOptions() throw()
{
}

MuonInJetOptionsDelegate *MuonInJetOptions::delegate() const
{
    return _delegate;
}

void MuonInJetOptions::setDelegate(MuonInJetOptionsDelegate *delegate)
{
    _delegate = delegate;
}

void MuonInJetOptions::init()
{
    _description.reset(new po::options_description("Muon-in-Jet Options"));
    _description->add_options()
        ("tag",
            po::value<string>()->default_value("TCHEL")->notifier(
                boost::bind(&MuonInJetOptions::optionTagIsSet, this, _1)),
            "muon-in-jet tagger Operating Point")

        ("away-tag",
            po::value<string>()->default_value("TCHPL")->notifier(
                boost::bind(&MuonInJetOptions::optionAwayTagIsSet, this, _1)),
            "away jet tagger Operating Point")

        ("muon-pt",
            po::value<double>()->default_value(5.0)->notifier(
                boost::bind(&MuonInJetOptions::optionMuonPtIsSet, this, _1)),
            "Muon pT cut")
    ;
}

po::options_description *MuonInJetOptions::description() const
{
    return _description.get();
}

void MuonInJetOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "Muon-In-Jet Options" << endl;
    out << setw(25) << left << " [+] Jet Tag" << _tag << endl;
    out << setw(25) << left << " [+] Away Jet Tag" << _awayTag << endl;
    out << setw(25) << left << " [+] Muon pT" << _muonPt << endl;
}



void MuonInJetOptions::optionTagIsSet(const string &tag)
{
    if (tag.empty())
        throw runtime_error("Tag is empty");

    _tag = tag;

    if (_delegate)
        _delegate->optionTagIsSet(tag);
}

void MuonInJetOptions::optionAwayTagIsSet(const string &tag)
{
    if (tag.empty())
        throw runtime_error("Away Tag is empty");

    _awayTag = tag;

    if (_delegate)
        _delegate->optionAwayTagIsSet(tag);
}

void MuonInJetOptions::optionMuonPtIsSet(const double &pt)
{
    if (0 > pt)
        throw runtime_error("Negative muon pT value");

    _muonPt = pt;

    if (_delegate)
        _delegate->optionMuonPtIsSet(pt);
}

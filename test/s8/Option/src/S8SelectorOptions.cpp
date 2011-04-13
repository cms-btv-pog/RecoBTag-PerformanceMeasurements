/**
 * S8SelectorOptions
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <stdexcept>
#include <string>

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "Option/interface/S8SelectorOptionsDelegate.h"

#include "Option/interface/S8SelectorOptions.h"

using std::string;
using std::runtime_error;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using s8::S8SelectorOptions;
using s8::S8SelectorOptionsDelegate;

S8SelectorOptions::S8SelectorOptions() throw():
    _min_jet_pt(0),
    _min_jet_eta(0),
    _max_jet_eta(0)
{
    _delegate = 0;
}

S8SelectorOptions::~S8SelectorOptions() throw()
{
}

S8SelectorOptionsDelegate *S8SelectorOptions::delegate() const
{
    return _delegate;
}

void S8SelectorOptions::setDelegate(S8SelectorOptionsDelegate *delegate)
{
    _delegate = delegate;
}

void S8SelectorOptions::init()
{
    _description.reset(new po::options_description("Muon-in-Jet Options"));
    _description->add_options()
        ("jet-pt",
            po::value<double>()->default_value(0)->notifier(
                boost::bind(&S8SelectorOptions::optionJetPtIsSet,
                            this, _1)),
            "Cut all jets by pT")

        ("jet-eta",
            po::value<string>()->default_value("")->notifier(
                boost::bind(&S8SelectorOptions::optionJetEtaIsSet,
                            this, _1)),
            "Cut all jets by eta. Format: --jet-eta=MIN..MAX")
    ;
}

po::options_description *S8SelectorOptions::description() const
{
    return _description.get();
}

void S8SelectorOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "S8 Selector Options" << endl;
    out << setw(25) << left << " [+] Jet pT";
    if (_min_jet_pt)
        out << _min_jet_pt;
    else
        out << "---";
    out << endl;

    out << setw(25) << left << " [+] Jet Eta";
    if (_min_jet_eta ||
        _max_jet_eta)
    {
        if (_min_jet_eta)
            out << _min_jet_eta;

        out << "..";

        if (_max_jet_eta)
            out << _max_jet_eta;

        out << endl;
    }
    else
        out << "---" << endl;
}



void S8SelectorOptions::optionJetPtIsSet(const double &pt)
{
    if (0 > pt)
        throw runtime_error("Negative jet pT value");

    _min_jet_pt = pt;

    if (_delegate)
        _delegate->optionJetPtIsSet(pt);
}

void S8SelectorOptions::optionJetEtaIsSet(const string &value)
{
    if (value.empty())
        return;

    smatch matches;
    regex pattern("^(?:(\\d+(?:\\.\\d*))..(\\d+(?:\\.\\d*))|(\\d+(?:\\.\\d*))..|..(\\d+(?:\\.\\d*)))$");
    if (!regex_match(value, matches, pattern))
        throw runtime_error("Unsuppored jet eta format (only possitive values are supported). Use: MIN..MAX");

    if (matches[4].matched)
    {
        _max_jet_eta = lexical_cast<double>(matches[4]);
    }
    else if(matches[3].matched)
    {
        _min_jet_eta = lexical_cast<double>(matches[3]);
    }
    else
    {
        _min_jet_eta = lexical_cast<double>(matches[1]);
        _max_jet_eta = lexical_cast<double>(matches[2]);
    }

    if (_delegate &&
        (_min_jet_eta ||
         _max_jet_eta))

        _delegate->optionJetEtaIsSet(_min_jet_eta, _max_jet_eta);
}

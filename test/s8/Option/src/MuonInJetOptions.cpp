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
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "Option/interface/MuonInJetOptionsDelegate.h"

#include "Option/interface/MuonInJetOptions.h"

using std::string;
using std::runtime_error;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using s8::MuonInJetOptions;
using s8::MuonInJetOptionsDelegate;

MuonInJetOptions::MuonInJetOptions() throw()
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
            po::value<string>()->default_value("")->notifier(
                boost::bind(&MuonInJetOptions::optionMuonPtIsSet,
                            this, _1)),
            "Cut muon by pt. Format: --muon-pt=MIN..MAX")

        ("jet-pt",
            po::value<string>()->default_value("")->notifier(
                boost::bind(&MuonInJetOptions::optionJetPtIsSet,
                            this, _1)),
            "Cut jet by pt. Format: --jet-pt=MIN..MAX")

        ("jet-eta",
            po::value<string>()->default_value("")->notifier(
                boost::bind(&MuonInJetOptions::optionJetEtaIsSet,
                            this, _1)),
            "Cut jet by eta. Format: --jet-eta=MIN..MAX")
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
    out << setw(25) << left << " [+] Muon pT" << _muon_pt << endl;
    out << setw(25) << left << " [+] Jet pT" << _jet_pt << endl;
    out << setw(25) << left << " [+] Jet Eta" << _jet_eta << endl;
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

void MuonInJetOptions::optionMuonPtIsSet(const std::string &value)
{
    if (value.empty())
        return;

    parse(_muon_pt, value);

    if (_delegate)
        _delegate->optionMuonPtIsSet(_muon_pt);
}

void MuonInJetOptions::optionJetPtIsSet(const std::string &value)
{
    if (value.empty())
        return;

    parse(_jet_pt, value);

    if (_delegate)
        _delegate->optionJetPtIsSet(_jet_pt);
}

void MuonInJetOptions::optionJetEtaIsSet(const string &value)
{
    if (value.empty())
        return;

    parse(_jet_eta, value);

    if (_delegate)
        _delegate->optionJetEtaIsSet(_jet_eta);
}

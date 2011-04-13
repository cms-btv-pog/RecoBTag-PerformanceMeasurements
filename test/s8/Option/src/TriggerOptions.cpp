/**
 * TriggerOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "Option/interface/TriggerOptionsDelegate.h"
#include "S8Tree/interface/S8Trigger.h"
#include "S8Tree/interface/S8Tools.h"

#include "Option/interface/TriggerOptions.h"

using std::clog;
using std::endl;
using std::string;
using std::runtime_error;

namespace fs = boost::filesystem;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using s8::TriggerOptions;
using s8::TriggerOptionsDelegate;

TriggerOptions::TriggerOptions() throw():
    _use_trigger_prescale(false),
    _simulate_trigger(false)
{
    _delegate = 0;
}

TriggerOptions::~TriggerOptions() throw()
{
}

TriggerOptionsDelegate *TriggerOptions::delegate() const
{
    return _delegate;
}

void TriggerOptions::setDelegate(TriggerOptionsDelegate *delegate)
{
    _delegate = delegate;
}

void TriggerOptions::init()
{
    _description.reset(new po::options_description("Trigger Options"));
    _description->add_options()
        ("trigger",
            po::value<Triggers>()->multitoken()->notifier(
                boost::bind(&TriggerOptions::optionTriggerIsSet, this, _1)),
            "Trigger(s) to be used (space separated in case of many). Syntax: hlt_xxx[:Prescale]")

        ("use-trigger-prescale",
            po::value<bool>()->implicit_value(true)->notifier(
                boost::bind(&TriggerOptions::optionUseTriggerPrescaleIsSet,
                            this, _1)),
            "Use trigger prescale in the analysis")

        ("simulate-trigger",
            po::value<bool>()->implicit_value(true)->notifier(
                boost::bind(&TriggerOptions::optionSimulateTriggerIsSet,
                            this, _1)),
            "Simulate trigger if it is missing")

        ("reweight-trigger",
            po::value<string>()->notifier(
                boost::bind(&TriggerOptions::optionReweightTriggerIsSet,
                            this, _1)),
            "Reweight trigger")
    ;
}

po::options_description *TriggerOptions::description() const
{
    return _description.get();
}

void TriggerOptions::print(std::ostream &out) const
{
    using std::left;
    using std::setw;

    out << "Trigger Options" << endl;
    out << setw(25) << left << " [+] Triggers";
    if (_trigger_map.empty())
    {
       out << "---" << endl;
    }
    else
    {
        out << "(" << _trigger_map.size() << ")" << endl;
        out << endl;
        out << "      Name" << setw(25) << left << " "
            << " Hash" << endl;
        out << "---------------------------------------------------"
            << endl;
        for(TriggerMap::const_iterator trigger = _trigger_map.begin();
            _trigger_map.end() != trigger;
            ++trigger)
        {
            out << " [+] " << setw(30) << left << trigger->second
                << trigger->first << endl;
        }
        out << "---------------------------------------------------"
            << endl;
    }
    out << setw(25) << left << " [+] Use Prescale" <<
                (_use_trigger_prescale ? "YES" : "---") << endl;
    out << setw(25) << left << " [+] Simulate" <<
                (_simulate_trigger ? "YES" : "---") << endl;
    out << setw(25) << left << " [+] Reweight" <<
                (_reweight_trigger.empty() ? "---" : _reweight_trigger)
                << endl;
}



void TriggerOptions::optionTriggerIsSet(const Triggers &triggers)
{
    using std::cerr;
    using s8::tools::make_hash;

    if (triggers.empty())
        throw runtime_error("No trigger is specified");

    _trigger_map.clear();

    if (!_delegate)
    {
        clog << "Triggers are specified but not used: delegate is not set"
            << endl;

        return;
    }

    for(Triggers::const_iterator trigger = triggers.begin();
        triggers.end() != trigger;
        ++trigger)
    {
        string trigger_name = *trigger;
        boost::to_lower(trigger_name);

        smatch matches;
        regex pattern("^(hlt_\\w+?)(?:_[v](\\d+))?(?:#(\\d+))?$");
        if (!regex_match(trigger_name, matches, pattern))
        {
            cerr << "Didn't understand trigger: " << *trigger << endl;
            cerr << "Syntax: hlt_xxx[:prescale]" << endl;
            cerr << endl;

            continue;
        }

        Trigger s8_trigger;
        s8_trigger.setHash(make_hash(matches[1]));

        if (matches[2].matched)
            s8_trigger.setVersion(lexical_cast<int>(matches[2]));

        if (matches[3].matched)
            s8_trigger.setPrescale(lexical_cast<int>(matches[3]));

        _delegate->optionTriggerIsSet(s8_trigger);

        _trigger_map[s8_trigger.hash()] = *trigger;
    }
}

void TriggerOptions::optionUseTriggerPrescaleIsSet(const bool &value)
{
    _use_trigger_prescale = value;

    if (_delegate)
        _delegate->optionUseTriggerPrescaleIsSet(value);
}

void TriggerOptions::optionSimulateTriggerIsSet(const bool &value)
{
    _simulate_trigger = value;

    if (_delegate)
        _delegate->optionSimulateTriggerIsSet(value);
}

void TriggerOptions::optionReweightTriggerIsSet(const string &filename)
{
    // Test if file exists
    //
    fs::path file(filename);
    if (!fs::exists(file))

        throw runtime_error("Reweights file does not exist: " + filename);

    if (!fs::is_regular_file(file) &&
        !fs::is_symlink(file))

        throw runtime_error("It seems reweights is neither regular file nor symlink: " + filename);

    _reweight_trigger = filename;

    if (_delegate)
        _delegate->optionReweightTriggerIsSet(filename);
}

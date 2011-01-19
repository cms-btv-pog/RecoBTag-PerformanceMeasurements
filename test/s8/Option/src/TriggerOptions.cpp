/**
 * TriggerOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <stdexcept>
#include <string>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include "Option/interface/TriggerOptionsDelegate.h"

#include "Option/interface/TriggerOptions.h"

using std::string;
using std::runtime_error;

using s8::TriggerOptions;
using s8::TriggerOptionsDelegate;

TriggerOptions::TriggerOptions() throw()
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
            po::value<Triggers>()->notifier(
                boost::bind(&TriggerOptions::optionTriggerIsSet, this, _1)),
            "Trigger(s) to be used (wildcard is accepted)")
    ;
}

po::options_description *TriggerOptions::description() const
{
    return _description.get();
}

void TriggerOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "Trigger Options" << endl;
    out << setw(25) << left << " [+] Triggers";
    if (_triggers.empty())
    {
       out << "not used" << endl;
    }
    else
    {
        for(Triggers::const_iterator trigger = _triggers.begin();
            _triggers.end() != trigger;
            ++trigger)
        {
            out << "       " << *trigger << endl;
        }
    }
}



void TriggerOptions::optionTriggerIsSet(const Triggers &triggers)
{
    if (triggers.empty())
        throw runtime_error("No trigger is specified");

    _triggers = triggers;

    if (_delegate)
    {
        for(Triggers::const_iterator trigger = triggers.begin();
            triggers.end() != trigger;
            ++trigger)
        {
            _delegate->optionTriggerIsSet(*trigger);
        }
    }
}

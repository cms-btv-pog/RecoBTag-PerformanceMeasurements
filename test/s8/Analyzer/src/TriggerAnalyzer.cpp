/**
 * TriggerAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <TDirectory.h>
#include <TH1I.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8EventID.h"
#include "S8Tree/interface/S8TriggerCenter.h"
#include "S8Tree/interface/S8Fwd.h"

#include "Analyzer/interface/TriggerAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;
using std::find_if;
using std::left;
using std::ostringstream;
using std::setw;

using s8::Trigger;
using s8::TriggerAnalyzer;

class TriggerPredicator
{
    public:
        TriggerPredicator()
        {
            _search_trigger = 0;
        }

        void setSearchTrigger(const Trigger &trigger)
        {
            _search_trigger = &trigger;
        }

        bool operator()(const Trigger *trigger) const
        {
            if (!_search_trigger)
                return false;

            return _search_trigger->hash() == trigger->hash();
        }

    private:
        const Trigger *_search_trigger;
};



TriggerAnalyzer::TriggerCounter::TriggerCounter(const Trigger &hlt):
    trigger(hlt),
    counter(0),
    name("")
{
    ostringstream name;
    name << "ps_" << trigger.hash();

    prescale.reset(new TH1I(name.str().c_str(), "Prescale vs Run",
                            3000, 147000, 150000));
}



TriggerAnalyzer::TriggerAnalyzer() throw()
{
}

TriggerAnalyzer::~TriggerAnalyzer() throw()
{
}



// Analyzer interface
//
void TriggerAnalyzer::init()
{
}

void TriggerAnalyzer::treeDidLoad(const TreeInfo *,
                                  const TriggerCenter *triggerCenter)
{
    using std::setw;
    using std::left;

    cout << "Triggers List is saved in the log file" << endl
        << endl;

    clog << "Triggers:" << endl;
    const TriggerCenter::TriggerMap &triggers = triggerCenter->triggers();
    for(TriggerCenter::TriggerMap::const_iterator trigger = triggers.begin();
        triggers.end() != trigger;
        ++trigger)
    {
        clog << " > " << setw(20) << left << trigger->first << " "
            << trigger->second << endl;
    }
    clog << endl;

    for(Triggers::iterator trigger = _triggers.begin();
        _triggers.end() != trigger;
        ++trigger)
    {
        TriggerCenter::TriggerMap::const_iterator found_trigger =
            triggers.find(trigger->trigger.hash());

        if (triggers.end() == found_trigger)
            continue;

        trigger->name = found_trigger->second;
    }
}

void TriggerAnalyzer::eventDidLoad(const Event *event)
{
    if (_triggers.empty())
        return;

    TriggerPredicator trigger_predicator;

    for(Triggers::iterator trigger = _triggers.begin();
        _triggers.end() != trigger;
        ++trigger)
    {
        trigger_predicator.setSearchTrigger(trigger->trigger);

        s8::Triggers::const_iterator found_trigger =
            find_if(event->triggers()->begin(),
                    event->triggers()->end(),
                    trigger_predicator);

        if (event->triggers()->end() == found_trigger ||
            !(**found_trigger))

            continue;

        // Trigger is found
        //
        ++(trigger->counter);
        const int bin = trigger->prescale->FindBin(event->id()->run());
        trigger->prescale->SetBinContent(bin,
                (*found_trigger)->prescale());
    }
}

void TriggerAnalyzer::print(std::ostream &out) const
{
    out << "Trigger cutflow (events passed trigger)" << endl;
    out << "      Name" << setw(25) << left << " "
        << " Hash" << endl;
    out << "---------------------------------------------------"
        << endl;
    for(Triggers::const_iterator trigger = _triggers.begin();
        _triggers.end() != trigger;
        ++trigger)
    {
        out << " [+] " << setw(30) << left;
        if (trigger->name.empty())
            out << trigger->trigger.hash();
        else
            out << trigger->name;

        out << trigger->counter << endl;
    }
    out << "---------------------------------------------------"
        << endl;
}

void TriggerAnalyzer::save(TDirectory *directory) const
{
    TDirectory *output = directory->mkdir("TriggerAnalyzer");
    if (!output)
    {
        cerr << "TriggerAnalyzer: Failed to create output folder: "
            << "no output is saved" << endl;

        return;
    }

    output->cd();

    for(Triggers::const_iterator trigger = _triggers.begin();
        _triggers.end() != trigger;
        ++trigger)
    {
        ostringstream title;
        title << trigger->name << ": Prescale vs Run";

        trigger->prescale->SetTitle(title.str().c_str());
        trigger->prescale->Write();
    }
}

// Trigger options
//
void TriggerAnalyzer::optionTriggerIsSet(const Trigger &trigger)
{
    _triggers.push_back(TriggerCounter(trigger));
}

/**
 * S8Selector
 * s8
 *
 * Created by Samvel Khalatian on Feb 3, 2011
 * Copyright 2010, All rights reserved
 */

#include <algorithm>
#include <iomanip>
#include <ostream>

#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8GenEvent.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8TriggerCenter.h"
#include "S8Tree/interface/S8Tools.h"
#include "Selector/interface/S8Selector.h"

using std::endl;
using std::find_if;
using std::left;
using std::setw;

using s8::S8Selector;

S8Selector::S8Selector() throw():
    _gluon_splitting(KEEP),
    _is_first_event(true),
    _simulate_trigger(false)
{
    using s8::tools::make_hash;

    _hlt_jet20u.setHash(make_hash("hlt_btagmu_jet20u"));
    _hlt_dijet20u.setHash(make_hash("hlt_btagmu_dijet20u"));
    _hlt_dijet30u.setHash(make_hash("hlt_btagmu_dijet30u"));
}

S8Selector::~S8Selector() throw()
{
}

void S8Selector::optionGluonSplittingIsSet(const GluonSplitting &value)
{
    _gluon_splitting = value;
}

void S8Selector::optionPtHatIsSet(const Range &value)
{
    _pt_hat = value;
}

void S8Selector::optionTriggerIsSet(const Trigger &trigger)
{
    _triggers.push_back(TriggerCounter(trigger));
}

void S8Selector::optionSimulateTriggerIsSet(const bool &value)
{
    _simulate_trigger = value;
}

void S8Selector::optionPrimaryVerticesIsSet(const Range &primary_vertices)
{
    _primary_vertices = primary_vertices;
}

void S8Selector::treeDidLoad(const TriggerCenter *trigger_center)
{
    const TriggerCenter::TriggerMap &triggers = trigger_center->triggers();

    // Search for trigger in the Trigger Menu
    //
    for(Triggers::iterator trigger = _triggers.begin();
        _triggers.end() != trigger;
        ++trigger)
    {
        trigger->disable = true;
        trigger->simulate = false;

        TriggerCenter::TriggerMap::const_iterator found_trigger =
            triggers.find(trigger->trigger.hash());

        // check if trigger is found
        //
        if (triggers.end() == found_trigger)
        {
            // Trigger is not found
            //
            if (_simulate_trigger &&
                    
                (trigger->trigger.hash() != _hlt_dijet20u.hash() ||
                 trigger->trigger.hash() != _hlt_dijet30u.hash()))

                trigger->simulate = true;
            else
                continue;
        }
        else
            trigger->name = found_trigger->second;

        trigger->disable = false;
    }
}

const s8::Event *S8Selector::operator()(const Event *event)
{
    if (!event)
        return 0;

    ++_events.total;

    if (!isValueInRange(event->primaryVertices()->size(), _primary_vertices))
        return 0;

    ++_events.primary_vertices;

    if (!isValueInRange(event->gen()->ptHat(), _pt_hat))
        return 0;

    ++_events.pt_hat;

    // Triggers
    //
    if (!_triggers.empty() &&
        !processTriggers(event))

        return 0;

    ++_events.trigger;

    // Gluon Splitting
    // [see PythiaOptionsDelegate for description]
    //
    switch(_gluon_splitting)
    {
        case KEEP:
            break;

        case BB:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB))
                    break;

                return 0;
            }

        case CC:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::CC))
                    break;

                return 0;
            }

        case ONLY:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB) ||
                    event->gen()->isGluonSplitting(s8::GenEvent::CC))

                    break;

                return 0;
            }

        case NO_BB:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB))
                    return 0;

                break;
            }

        case NO_CC:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::CC))
                    return 0;

                break;
            }

        case REMOVE:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB) ||
                    event->gen()->isGluonSplitting(s8::GenEvent::CC))

                    return 0;

                break;
            }

        case ADD_BB:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB))
                    break;

                _is_first_event = !_is_first_event;

                if (_is_first_event)
                    return 0;

                break;
            }

        case ADD_CC:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::CC))
                    break;

                _is_first_event = !_is_first_event;

                if (_is_first_event)
                    return 0;

                break;
            }

        case ENHANCE:
            {
                if (event->gen()->isGluonSplitting(s8::GenEvent::BB) ||
                    event->gen()->isGluonSplitting(s8::GenEvent::CC))

                    break;

                _is_first_event = !_is_first_event;

                if (_is_first_event)
                    return 0;

                break;
            }
    }

    ++_events.gluon_splitting;

    return event;
}

void S8Selector::print(std::ostream &out) const
{
    out << "S8Selector cutflow (events passed)" << endl;
    out << " " << setw(20) << left << "Total" << _events.total << endl;
    out << " " << setw(20) << left << "Primary Vertices"
        << _events.primary_vertices << endl;
    out << " " << setw(20) << left << "Pt Hat" << _events.pt_hat << endl;
    out << " " << setw(20) << left << "Triggers" << _events.trigger << endl;
    out << " " << setw(20) << left << "Gluon Splitting"
        << _events.gluon_splitting <<endl;
    out << endl;
    out << "Triggers Cutflow" << endl;
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



bool S8Selector::processTriggers(const Event *event)
{
    bool triggers_pass = false;

    // Search for user defined triggers among those that passed the event
    //
    for(Triggers::iterator trigger = _triggers.begin();
        _triggers.end() != trigger;
        ++trigger)
    {
        if (trigger->disable)
            continue;

        if (trigger->simulate)
        {
            if (!simulateTrigger(trigger->trigger, event))
                continue;
        }
        else
        {
            if (!didTriggerPass(trigger->trigger, event))
                continue;
        }

        ++(trigger->counter);

        triggers_pass = true;
    }

    return triggers_pass;
}

bool S8Selector::didTriggerPass(const Trigger &trigger,
                                const Event *event)
{
    _trigger_predicator.setSearchTrigger(trigger);

    s8::Triggers::const_iterator found_trigger =
        find_if(event->triggers()->begin(),
                event->triggers()->end(),
                _trigger_predicator);

    // Check if trigger was found and passed the event
    //
    return event->triggers()->end() != found_trigger &&
        (**found_trigger);
}

bool S8Selector::simulateTrigger(const Trigger &trigger, const Event *event)
{
    if (!didTriggerPass(_hlt_jet20u, event))
        return false;

    double max_pt = 0;
    if (_hlt_dijet20u.hash() == trigger.hash())
        max_pt = 60;
    else if (_hlt_dijet30u.hash() == trigger.hash())
        max_pt = 80;
    else
        return false;

    // Apply offline jet pT cut
    //
    int trigger_jets = 0;
    for(Jets::const_iterator jet = event->jets()->begin();
        event->jets()->end() != jet &&
            2 > trigger_jets;
        ++jet)
    {
        if (max_pt < (*jet)->p4()->Pt() &&
            3.0 > fabs((*jet)->p4()->Eta()))

            ++trigger_jets;
    }

    return 1 < trigger_jets;
}

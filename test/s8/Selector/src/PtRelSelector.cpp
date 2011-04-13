/**
 * PtRelSelector
 * s8
 *
 * Created by Samvel Khalatian on Feb 3, 2011
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Selector/interface/PtRelSelector.h"

using s8::Event;
using s8::PtRelSelector;

PtRelSelector::PtRelSelector() throw():
    _min_muon_pt(0)
    
{
    _modified_event = new Event();
    _modified_muons = new Leptons();
}

PtRelSelector::~PtRelSelector() throw()
{
    delete _modified_muons;
    delete _modified_event;
}

void PtRelSelector::treeDidLoad(const TriggerCenter *)
{
}

const Event *PtRelSelector::operator()(const Event *event)
{
    // Event should contain one muon with pT above XXX GeV/c
    //
    _modified_muons->clear();
    for(Leptons::const_iterator muon = event->muons()->begin();
            event->muons()->end() != muon;
            ++muon)
    {
        if (_min_muon_pt > (*muon)->p4()->Pt())
            continue;

        _modified_muons->push_back(*muon);
    }

    if (1 != _modified_muons->size())
        return 0;

    if (2 != event->jets()->size())
        return 0;

    // Event has passed the selection. Copying is cheap b/c only pointers
    // to collections are copied.
    //
    *_modified_event = *event;

    // Substitute muons and jets
    //
    _modified_event->setMuons(_modified_muons);

    return _modified_event;
}

void PtRelSelector::optionMuonPtIsSet(const double &value)
{
    _min_muon_pt = value;
}

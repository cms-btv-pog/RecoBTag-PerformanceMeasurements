/**
 * TwoMuonsInJet
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Utility/interface/LeptonInJet.h"
#include "Utility/interface/TaggerOperatingPoint.h"
#include "Utility/interface/TwoMuonsInJetDelegate.h"

#include "Utility/interface/TwoMuonsInJet.h"

using s8::TwoMuonsInJet;
using s8::TwoMuonsInJetDelegate;

TwoMuonsInJet::TwoMuonsInJet() throw()
{
    _delegate = 0;

    _leptonInJet.reset(new LeptonInJet());
    _leptonInJet->setDelegate(this);

    _awayJetTaggerOperatingPoint = 0;

    _firstLepton = 0;
    _secondLepton = 0;
    _jet = 0;
}

TwoMuonsInJet::~TwoMuonsInJet() throw()
{
}

TwoMuonsInJetDelegate *TwoMuonsInJet::delegate() const
{
    return _delegate;
}

void TwoMuonsInJet::setDelegate(TwoMuonsInJetDelegate *delegate)
{
    _delegate = delegate;
}

void TwoMuonsInJet::setAwayJetTaggerOperatingPoint(const TaggerOperatingPoint *tag)
{
    _awayJetTaggerOperatingPoint = tag;
}

void TwoMuonsInJet::operator()(const Event *event)
{
    if (!_delegate ||
        !_awayJetTaggerOperatingPoint)
        return;

    // One or more primary vertices should exist
    //
    if (1 > event->primaryVertices()->size())
        return;

    // Event should have at least 2 jets
    //
    if (2 > event->jets()->size())
        return;

    // There should be at lease one muon
    //
    if (1 > event->muons()->size())
        return;

    const PrimaryVertex *primaryVertex = *(event->primaryVertices()->begin());
    for(Jets::const_iterator jet = event->jets()->begin();
        event->jets()->end() != jet;
        ++jet)
    {
        _firstLepton = 0;
        _secondLepton = 0;
        _jet = 0;

        if (!_delegate->muonInJetShouldProcessJet(*jet))
            continue;

        // away jet exists b/c number of jets in the event was tested before
        // (at least 2 jets should exist in the event)
        //
        _leptonInJet->operator()(primaryVertex, event->muons(), *jet);

        // There should be at least first lepton
        //
        if (!_firstLepton)
            continue;

        _delegate->muonIsInJetPlusAwayJet(_firstLepton,
                                          _secondLepton,
                                          _jet);

        // check if there is tagged away jet
        //
        bool awayJetIsNotTagged = true;
        for(Jets::const_iterator awayJet = event->jets()->begin();
            event->jets()->end() != awayJet;
            ++awayJet)
        {
            // Skip muon-in-jet
            //
            if (jet == awayJet)
                continue;

            if ((*awayJet)->btag(_awayJetTaggerOperatingPoint->btag())
                    <= *_awayJetTaggerOperatingPoint)

                continue;

            awayJetIsNotTagged = false;

            break;
        }
        
        if (awayJetIsNotTagged)
            continue;

        _delegate->muonIsInJetPlusTaggedAwayJet(_firstLepton,
                                                _secondLepton,
                                                _jet);
    }

    _firstLepton = 0;
    _secondLepton = 0;
    _jet = 0;
}

void TwoMuonsInJet::leptonIsInJet(const Lepton *lepton, const Jet *jet)
{
    _jet = jet;

    // Check if first lepton has lower pT than argument
    //
    if (_firstLepton &&
        _firstLepton->p4()->Pt() < lepton->p4()->Pt())
    {
        // Move first lepton into second
        //
        _secondLepton = _firstLepton;
        _firstLepton = lepton;

        return;
    }

    // Is there first lepton saved?
    //
    if (!_firstLepton)
    {
        _firstLepton = lepton;

        return;
    }

    // Does second lepton has lower pT than argument?
    //
    if (_secondLepton &&
        _secondLepton->p4()->Pt() >= lepton->p4()->Pt())
    {
        // No
        //
        return;
    }

    // Store second lepton
    //
    _secondLepton = lepton;
}

// Go over all muons in order to pick the one with highest pT
//
bool TwoMuonsInJet::shouldLookForMoreLeptons()
{
    return true;
}

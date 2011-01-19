/**
 * MuonInJet
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include <TLorentzVector.h>

#include "IO/interface/Event.h"
#include "S8Tree/interface/S8Fwd.h"
#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"
#include "Utility/interface/LeptonInJet.h"
#include "Utility/interface/MuonInJetDelegate.h"
#include "Utility/interface/TaggerOperatingPoint.h"

#include "Utility/interface/MuonInJet.h"

using s8::MuonInJet;
using s8::MuonInJetDelegate;

MuonInJet::MuonInJet() throw()
{
    _delegate = 0;

    _leptonInJet.reset(new LeptonInJet());
    _leptonInJet->setDelegate(this);

    _awayJetTaggerOperatingPoint = 0;

    _muonMinimumPtCut = 0;

    _lepton = 0;
    _jet = 0;
}

MuonInJet::~MuonInJet() throw()
{
}

MuonInJetDelegate *MuonInJet::delegate() const
{
    return _delegate;
}

void MuonInJet::setDelegate(MuonInJetDelegate *delegate)
{
    _delegate = delegate;
}

void MuonInJet::setAwayJetTaggerOperatingPoint(const TaggerOperatingPoint *tag)
{
    _awayJetTaggerOperatingPoint = tag;
}

void MuonInJet::setMuonMinimumPtCut(const double &pt)
{
    using std::runtime_error;

    if (0 > pt)
        throw runtime_error("Negative pT supplied");

    _muonMinimumPtCut = pt;
}

void MuonInJet::operator()(const Event *event)
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
        _lepton = 0;
        _jet = 0;

        if (!_delegate->muonInJetShouldProcessJet(*jet))
            continue;

        // away jet exists b/c number of jets in the event was tested before
        // (at least 2 jets should exist in the event)
        //
        _leptonInJet->operator()(primaryVertex, event->muons(), *jet);

        if (!_lepton)
            continue;

        _delegate->muonIsInJetPlusAwayJet(_lepton, _jet);

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

        _delegate->muonIsInJetPlusTaggedAwayJet(_lepton, _jet);
    }

    _lepton = 0;
    _jet = 0;
}

void MuonInJet::leptonIsInJet(const Lepton *lepton, const Jet *jet)
{
    if (_muonMinimumPtCut &&
        lepton->p4()->Pt() < _muonMinimumPtCut)

        return;

    if (_lepton &&
        _lepton->p4()->Pt() >= lepton->p4()->Pt())

        return;

    _lepton = lepton;
    _jet = jet;
}

// Go over all muons in order to pick the one with highest pT
//
bool MuonInJet::shouldLookForMoreLeptons()
{
    return true;
}

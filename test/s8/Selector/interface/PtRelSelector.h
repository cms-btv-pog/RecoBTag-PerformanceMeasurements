/**
 * PtRelSelector
 * s8
 *
 * Created by Samvel Khalatian on Feb 3, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PTREL_SELECTOR
#define S8_PTREL_SELECTOR

#include "Option/interface/MuonInJetOptionsDelegate.h"
#include "Selector/interface/Selector.h"

namespace s8
{
    /*
     * The selection is defined on Feb 3, 2011
     *
     * 1. Select events with only one muon with pT > 6 GeV/c
     * 2. Select events with only two jets with quality cuts passed:
     *      pT > 30 GeV/c, |eta| < 2.4, JetID, etc.
     * 3. Apply additional cut on the muon-in-jet jet: pT > X GeV/c
     *    (where X is the offline jet pT cut depending on the trigger used).
     *
     * Note: jets quality cuts are applied elsewhere, pre-selection for example.
     *       cut on the muon-in-jet jet is used in the MuonInJet
     */
    class PtRelSelector : public Selector,
                          public MuonInJetOptionsDelegate
    {
        public:
            PtRelSelector() throw();
            virtual ~PtRelSelector() throw();

            // Muon In Jet options
            //
            virtual void optionMuonPtIsSet(const double &);

            virtual void treeDidLoad(const TriggerCenter *);

            virtual const Event *operator()(const Event *);

        private:
            Event   *_modified_event;
            Leptons *_modified_muons;

            double _min_muon_pt;
    };
}

#endif

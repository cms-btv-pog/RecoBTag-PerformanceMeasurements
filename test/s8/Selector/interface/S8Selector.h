/**
 * S8Selector
 * s8
 *
 * Created by Samvel Khalatian on Feb 3, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_S8SELECTOR
#define S8_S8SELECTOR

#include <vector>

#include "Option/interface/MiscOptionsDelegate.h"
#include "Option/interface/PythiaOptionsDelegate.h"
#include "Option/interface/TriggerOptionsDelegate.h"
#include "Selector/interface/Selector.h"
#include "S8Tree/interface/S8Trigger.h"
#include "Selector/interface/TriggerNamePredicator.h"
#include "Utility/interface/Range.h"

namespace s8
{
    class TriggerNamePredicator;

    class S8Selector : public Selector,
                       public PythiaOptionsDelegate,
                       public TriggerOptionsDelegate
    {
        public:
            S8Selector() throw();
            virtual ~S8Selector() throw();

            // Pythia options
            //
            virtual void optionGluonSplittingIsSet(const GluonSplitting &);
            virtual void optionPtHatIsSet(const Range &);

            // Trigger options
            //
            virtual void optionTriggerIsSet(const Trigger &);
            virtual void optionSimulateTriggerIsSet(const bool &);

            // Misc options
            //
            virtual void optionPrimaryVerticesIsSet(const Range &);

            virtual void treeDidLoad(const TriggerCenter *);

            virtual const Event *operator()(const Event *);

            virtual void print(std::ostream &) const;

        private:
            bool processTriggers(const Event *);

            bool simulateTrigger(const Trigger &, const Event *);
            bool didTriggerPass(const Trigger &, const Event *);

            GluonSplitting _gluon_splitting;
            Range          _pt_hat;
            bool           _is_first_event;

            struct TriggerCounter
            {
                TriggerCounter(const Trigger &hlt):
                    trigger(hlt),
                    counter(0),
                    name(""),
                    disable(true),
                    simulate(false)
                {}

                Trigger     trigger;
                int         counter;
                std::string name;
                bool        disable;
                bool        simulate;
            };

            typedef std::vector<TriggerCounter> Triggers;

            Triggers _triggers;

            struct EventCounter
            {
                EventCounter():
                    total(0),
                    primary_vertices(0),
                    pt_hat(0),
                    gluon_splitting(0)
                {}

                int total;
                int primary_vertices;
                int pt_hat;
                int trigger;
                int gluon_splitting;
            };

            EventCounter _events;

            Trigger      _hlt_jet20u;
            Trigger      _hlt_dijet20u;
            Trigger      _hlt_dijet30u;

            bool _simulate_trigger;

            TriggerNamePredicator _trigger_predicator;

            Range _primary_vertices;
    };
}

#endif

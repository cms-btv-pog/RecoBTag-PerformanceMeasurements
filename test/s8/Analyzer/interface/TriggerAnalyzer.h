/**
 * TriggerAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_ANALYZER
#define S8_TRIGGER_ANALYZER

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/TriggerOptionsDelegate.h"
#include "S8Tree/interface/S8Trigger.h"

class TH1;

namespace s8
{
    class TaggerOperatingPoint;

    class TriggerAnalyzer: public Analyzer,
                           public TriggerOptionsDelegate
    {
        public:
            TriggerAnalyzer() throw();
            virtual ~TriggerAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void treeDidLoad(const TreeInfo *,
                                     const TriggerCenter *);

            virtual void eventDidLoad(const Event *);
            virtual void print(std::ostream &) const;
            virtual void save(TDirectory *) const;

            // Trigger options
            //
            void optionTriggerIsSet(const Trigger &);

        private:
            struct TriggerCounter
            {
                typedef boost::shared_ptr<TH1> TH1Ptr;

                TriggerCounter(const Trigger &hlt);

                Trigger     trigger;
                int         counter;
                std::string name;
                TH1Ptr      prescale;
            };

            typedef std::vector<TriggerCounter> Triggers;

            Triggers   _triggers;
    };
}

#endif

/**
 * TriggerAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_ANALYZER
#define S8_TRIGGER_ANALYZER

#include "Analyzer/interface/Analyzer.h"

class TH1;

namespace s8
{
    class TaggerOperatingPoint;

    class TriggerAnalyzer: public Analyzer
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

            virtual void save(TDirectory *) const;
    };
}

#endif

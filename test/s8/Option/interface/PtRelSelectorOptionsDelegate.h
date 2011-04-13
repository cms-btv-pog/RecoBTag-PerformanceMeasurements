/**
 * PtRelSelectorOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PTREL_SELECTOR_OPTIONS_DELEGATE
#define S8_PTREL_SELECTOR_OPTIONS_DELEGATE

#include <string>

namespace s8
{
    class PtRelSelectorOptionsDelegate
    {
        public:
            PtRelSelectorOptionsDelegate() throw();
            virtual ~PtRelSelectorOptionsDelegate() throw();

            virtual void optionMuonPtIsSet(const double &);
    };
}

#endif

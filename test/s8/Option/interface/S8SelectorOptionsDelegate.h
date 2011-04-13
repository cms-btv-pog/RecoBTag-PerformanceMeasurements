/**
 * S8SelectorOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_S8SELECTOR_OPTIONS_DELEGATE
#define S8_S8SELECTOR_OPTIONS_DELEGATE

#include <string>

namespace s8
{
    class S8SelectorOptionsDelegate
    {
        public:
            S8SelectorOptionsDelegate() throw();
            virtual ~S8SelectorOptionsDelegate() throw();

            virtual void optionJetPtIsSet(const double &);
            virtual void optionJetEtaIsSet(const double &, const double &);
    };
}

#endif

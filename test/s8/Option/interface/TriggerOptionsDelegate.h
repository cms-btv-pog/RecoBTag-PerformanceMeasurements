/**
 * TriggerOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_OPTIONS_DELEGATE
#define S8_TRIGGER_OPTIONS_DELEGATE

#include <string>

namespace s8
{
    class TriggerOptionsDelegate
    {
        public:
            TriggerOptionsDelegate() throw();
            virtual ~TriggerOptionsDelegate() throw();

            virtual void optionTriggerIsSet(const std::string &);
    };
}

#endif

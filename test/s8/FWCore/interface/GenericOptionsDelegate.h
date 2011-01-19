/**
 * GenericOptionsDelegate
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_GENERIC_OPTIONS_DELEGATE
#define CORE_GENERIC_OPTIONS_DELEGATE

#include <string>

namespace core
{
    class GenericOptionsDelegate
    {
        public:
            virtual ~GenericOptionsDelegate() throw();

            virtual void optionDebugIsSet(const std::string &);
            virtual void optionOutputIsSet(const std::string &);
            virtual void optionEventsIsSet(const int &);
            virtual void optionInputIsSet(const std::string &);
    };
}

#endif

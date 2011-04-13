/**
 * MiscOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Mar 9, 2011
 * Copyright 2011, All rights reserved
 */

#ifndef S8_MISC_OPTIONS_DELEGATE
#define S8_MISC_OPTIONS_DELEGATE

namespace s8
{
    class Range;

    class MiscOptionsDelegate
    {
        public:
            MiscOptionsDelegate() throw();
            virtual ~MiscOptionsDelegate() throw();

            virtual void optionPrimaryVerticesIsSet(const Range &);
    };
}

#endif

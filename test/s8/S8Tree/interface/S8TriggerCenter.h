/**
 * TriggerCenter
 * s8 
 *
 * Created by Samvel Khalatian on Nov 5, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGERCENTER
#define S8_TRIGGERCENTER

#include <map>

#include <TObject.h>
#include <Rtypes.h>

namespace s8
{
    class TriggerCenter: public TObject
    {
        public:
            typedef unsigned int Hash;

            typedef std::map<Hash, std::string> TriggerMap;

            TriggerCenter();

            TriggerMap &triggers();
            const TriggerMap &triggers() const;

            void merge(const TriggerCenter &);

        private:
            TriggerMap _triggers;

            ClassDef(TriggerCenter, 1);
    };
}

#endif

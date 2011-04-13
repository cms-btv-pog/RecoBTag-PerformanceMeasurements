/**
 * TriggerNamePredicator
 * s8
 *
 * Created by Samvel Khalatian on Feb 8, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_NAME_PREDICATOR
#define S8_TRIGGER_NAME_PREDICATOR

namespace s8
{
    class Trigger;

    class TriggerNamePredicator
    {
        public:
            TriggerNamePredicator() throw();

            void setSearchTrigger(const Trigger &trigger);

            bool operator()(const Trigger *trigger) const;

        private:
            const Trigger *_search_trigger;
    };
}

#endif

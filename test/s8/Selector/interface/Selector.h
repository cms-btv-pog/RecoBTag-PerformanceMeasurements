/**
 * Selector
 * s8
 *
 * Created by Samvel Khalatian on Feb 3, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SELECTOR
#define S8_SELECTOR

#include <iosfwd>

namespace s8
{
    class Event;
    class TriggerCenter;

    /*
     * Apply select events matching custom requirements
     */
    class Selector
    {
        public:
            Selector() throw();
            virtual ~Selector() throw();

            // Pass Trigger Menu
            //
            virtual void treeDidLoad(const TriggerCenter *) = 0;

            // Always check pointer before use
            //
            virtual const Event *operator()(const Event *) = 0;

            virtual void print(std::ostream &) const = 0;
    };

    std::ostream &operator<<(std::ostream &, const Selector &);
}

#endif

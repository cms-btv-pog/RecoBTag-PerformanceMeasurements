/**
 * TriggerOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_OPTIONS
#define S8_TRIGGER_OPTIONS

#include <string>
#include <vector>

#include "FWCore/interface/Options.h"

namespace s8
{
    class TriggerOptionsDelegate;

    class TriggerOptions: public core::Options
    {
        public:
            TriggerOptions() throw();
            virtual ~TriggerOptions() throw();

            TriggerOptionsDelegate *delegate() const;
            void setDelegate(TriggerOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            typedef std::vector<std::string> Triggers;

            void optionTriggerIsSet(const Triggers &);

            TriggerOptionsDelegate                 *_delegate;
            std::auto_ptr<po::options_description>  _description;

            Triggers _triggers;
    };
}

#endif

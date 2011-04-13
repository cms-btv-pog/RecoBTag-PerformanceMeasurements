/**
 * TriggerOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_OPTIONS
#define S8_TRIGGER_OPTIONS

#include <map>
#include <memory>
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
            typedef unsigned int                Hash;
            typedef std::vector<std::string>    Triggers;
            typedef std::map<Hash, std::string> TriggerMap;

            void optionTriggerIsSet(const Triggers &);
            void optionUseTriggerPrescaleIsSet(const bool &);
            void optionSimulateTriggerIsSet(const bool &);
            void optionReweightTriggerIsSet(const std::string &);

            TriggerOptionsDelegate                 *_delegate;
            std::auto_ptr<po::options_description>  _description;

            TriggerMap _trigger_map;
            bool _use_trigger_prescale;
            bool _simulate_trigger;
            std::string _reweight_trigger;
    };
}

#endif

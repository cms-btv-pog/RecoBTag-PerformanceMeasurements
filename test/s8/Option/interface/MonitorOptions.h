/**
 * MonitorOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_OPTIONS
#define S8_MONITOR_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"

namespace s8
{
    class MonitorOptionsDelegate;
    class MuonInJetOptions;
    class PythiaOptions;
    class TriggerOptions;

    class MonitorOptions: public core::Options
    {
        public:
            MonitorOptions() throw();
            virtual ~MonitorOptions() throw();

            MonitorOptionsDelegate *delegate() const;
            void setDelegate(MonitorOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionDataIsSet(const bool &);

            MonitorOptionsDelegate                 *_delegate;
            std::auto_ptr<po::options_description>  _description;
            std::auto_ptr<po::options_description>  _hiddenDescription;

            std::auto_ptr<MuonInJetOptions> _muonInJetOptions;
            std::auto_ptr<PythiaOptions>    _pythiaOptions;
            std::auto_ptr<TriggerOptions>   _triggerOptions;
    };
}

#endif

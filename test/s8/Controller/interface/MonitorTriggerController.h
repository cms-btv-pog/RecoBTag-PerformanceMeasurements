/**
 * MonitorTriggerController
 * s8
 *
 * Created by Samvel Khalatian on Feb 25, 2011
 * Copyright 2011, All rights reserved
 */

#ifndef S8_MONITOR_TRIGGER_CONTROLLER
#define S8_MONITOR_TRIGGER_CONTROLLER

#include <memory>

#include "Controller/interface/AppController.h"

namespace s8
{
    class MonitorTriggerAnalyzer;
    class MonitorOptions;

    class MonitorTriggerController: public AppController
    {
        public:
            MonitorTriggerController() throw();
            virtual ~MonitorTriggerController() throw();

        protected:
            // core::AppController interface
            //
            virtual void init();

        private:
            // s8::AppController interface
            //
            virtual Analyzer *createAnalyzer();

            // core::AppController interface
            //
            virtual core::Options *createOptions();

            // InputFileDelegate interface
            //
            virtual bool inputFileShouldLoadTriggers();

            std::auto_ptr<MonitorTriggerAnalyzer> _analyzer;
            std::auto_ptr<MonitorOptions>  _options;
    };
}

#endif

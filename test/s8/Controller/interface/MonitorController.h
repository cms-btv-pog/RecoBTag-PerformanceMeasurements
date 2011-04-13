/**
 * MonitorController
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_CONTROLLER
#define S8_MONITOR_CONTROLLER

#include <memory>

#include "Controller/interface/AppController.h"

namespace s8
{
    class MonitorAnalyzer;
    class MonitorOptions;

    class MonitorController: public AppController
    {
        public:
            MonitorController() throw();
            virtual ~MonitorController() throw();

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

            std::auto_ptr<MonitorAnalyzer> _analyzer;
            std::auto_ptr<MonitorOptions>  _options;
    };
}

#endif

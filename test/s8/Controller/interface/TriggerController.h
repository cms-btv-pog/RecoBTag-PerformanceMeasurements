/**
 * TriggerController
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER_CONTROLLER
#define S8_TRIGGER_CONTROLLER

#include <memory>

#include "Controller/interface/AppController.h"

namespace s8
{
    class TriggerAnalyzer;
    class TriggerOptions;

    class TriggerController: public AppController
    {
        public:
            TriggerController() throw();
            virtual ~TriggerController() throw();

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

            std::auto_ptr<TriggerAnalyzer> _analyzer;
            std::auto_ptr<TriggerOptions>  _options;
    };
}

#endif

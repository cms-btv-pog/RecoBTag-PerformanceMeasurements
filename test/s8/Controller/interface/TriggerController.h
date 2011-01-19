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
    class MuonInJetOptions;

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

            // InputFileDelegate interface
            //
            virtual bool inputFileShouldContinue();

            std::auto_ptr<TriggerAnalyzer>    _analyzer;
    };
}

#endif

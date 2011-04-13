/**
 * MonitorController + SolverInputController
 * s8
 *
 * Created by Samvel Khalatian on Feb 14, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_SOLVER_INPUT_CONTROLLER
#define S8_MONITOR_SOLVER_INPUT_CONTROLLER

#include <memory>

#include "Controller/interface/AppController.h"

namespace s8
{
    class MonitorSolverInputAnalyzer;
    class SolverInputOptions;

    class MonitorSolverInputController: public AppController
    {
        public:
            MonitorSolverInputController() throw();
            virtual ~MonitorSolverInputController() throw();

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

            std::auto_ptr<MonitorSolverInputAnalyzer> _analyzer;
            std::auto_ptr<SolverInputOptions>  _options;
    };
}

#endif

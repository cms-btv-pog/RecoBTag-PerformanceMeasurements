/**
 * PtRelController
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PTREL_CONTROLLER
#define S8_PTREL_CONTROLLER

#include <memory>

#include "Controller/interface/AppController.h"
#include "Option/interface/SolverInputOptionsDelegate.h"

namespace s8
{
    class PtRelAnalyzer;
    class SolverInputOptions;

    class PtRelController: public AppController,
                           public SolverInputOptionsDelegate
    {
        public:
            PtRelController() throw();
            virtual ~PtRelController() throw();

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

            std::auto_ptr<PtRelAnalyzer>      _analyzer;
            std::auto_ptr<SolverInputOptions> _options;
    };
}

#endif

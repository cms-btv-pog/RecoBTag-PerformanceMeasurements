/**
 * SolverInputJetPtController
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PTREL_CONTROLLER
#define S8_PTREL_CONTROLLER

#include <memory>

#include "Controller/interface/AppController.h"

namespace s8
{
    class SolverInputJetPtAnalyzer;
    class SolverInputOptions;

    class SolverInputJetPtController: public AppController
    {
        public:
            SolverInputJetPtController() throw();
            virtual ~SolverInputJetPtController() throw();

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

            std::auto_ptr<SolverInputJetPtAnalyzer>    _analyzer;
            std::auto_ptr<SolverInputOptions> _options;
    };
}

#endif

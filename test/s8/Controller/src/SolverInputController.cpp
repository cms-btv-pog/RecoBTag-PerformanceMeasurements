/**
 * SolverInputController
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/SolverInputAnalyzer.h"
#include "Option/interface/SolverInputOptions.h"

#include "Controller/interface/SolverInputController.h"

using s8::SolverInputController;

SolverInputController::SolverInputController() throw()
{
}

SolverInputController::~SolverInputController() throw()
{
}

// Override parent init()
//
void SolverInputController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new SolverInputAnalyzer());
    _analyzer->init();

    _options.reset(new SolverInputOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

s8::Analyzer *SolverInputController::createAnalyzer()
{
    return _analyzer.get();
}

core::Options *SolverInputController::createOptions()
{
    return _options.get();
}

/**
 * MonitorSolverInputController
 * s8
 *
 * Created by Samvel Khalatian on Feb 14, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/MonitorSolverInputAnalyzer.h"
#include "Option/interface/SolverInputOptions.h"

#include "Controller/interface/MonitorSolverInputController.h"

using s8::MonitorSolverInputController;

MonitorSolverInputController::MonitorSolverInputController() throw()
{
}

MonitorSolverInputController::~MonitorSolverInputController() throw()
{
}

// Override parent init()
//
void MonitorSolverInputController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new MonitorSolverInputAnalyzer());
    _analyzer->init();

    _options.reset(new SolverInputOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

s8::Analyzer *MonitorSolverInputController::createAnalyzer()
{
    return _analyzer.get();
}

core::Options *MonitorSolverInputController::createOptions()
{
    return _options.get();
}

bool MonitorSolverInputController::inputFileShouldLoadTriggers()
{
    return true;
}

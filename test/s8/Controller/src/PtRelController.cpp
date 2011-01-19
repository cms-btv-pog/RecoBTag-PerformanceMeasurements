/**
 * PtRelController
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/PtRelAnalyzer.h"
#include "Option/interface/SolverInputOptions.h"

#include "Controller/interface/PtRelController.h"

using s8::PtRelController;

PtRelController::PtRelController() throw()
{
}

PtRelController::~PtRelController() throw()
{
}

s8::Analyzer *PtRelController::createAnalyzer()
{
    return _analyzer.get();
}

// Override parent init()
//
void PtRelController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new PtRelAnalyzer());
    _analyzer->init();

    _options.reset(new SolverInputOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

core::Options *PtRelController::createOptions()
{
    return _options.get();
}

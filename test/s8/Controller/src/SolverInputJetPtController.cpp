/**
 * SolverInputJetPtController
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/SolverInputJetPtAnalyzer.h"
#include "Option/interface/SolverInputOptions.h"

#include "Controller/interface/SolverInputJetPtController.h"

using s8::SolverInputJetPtController;

SolverInputJetPtController::SolverInputJetPtController() throw()
{
}

SolverInputJetPtController::~SolverInputJetPtController() throw()
{
}

// Override parent init()
//
void SolverInputJetPtController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new SolverInputJetPtAnalyzer());
    _analyzer->init();

    _options.reset(new SolverInputOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

s8::Analyzer *SolverInputJetPtController::createAnalyzer()
{
    return _analyzer.get();
}

core::Options *SolverInputJetPtController::createOptions()
{
    return _options.get();
}

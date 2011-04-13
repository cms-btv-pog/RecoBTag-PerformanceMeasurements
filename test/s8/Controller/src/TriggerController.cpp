/**
 * TriggerController
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/TriggerAnalyzer.h"
#include "Option/interface/TriggerOptions.h"

#include "Controller/interface/TriggerController.h"

using s8::TriggerController;

TriggerController::TriggerController() throw()
{
}

TriggerController::~TriggerController() throw()
{
}

// Override parent init()
//
void TriggerController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new TriggerAnalyzer());
    _analyzer->init();

    _options.reset(new TriggerOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

s8::Analyzer *TriggerController::createAnalyzer()
{
    return _analyzer.get();
}

core::Options *TriggerController::createOptions()
{
    return _options.get();
}

bool TriggerController::inputFileShouldLoadTriggers()
{
    return true;
}

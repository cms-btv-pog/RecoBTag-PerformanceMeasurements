/**
 * TriggerController
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/TriggerAnalyzer.h"

#include "Controller/interface/TriggerController.h"

using s8::TriggerController;

TriggerController::TriggerController() throw()
{
}

TriggerController::~TriggerController() throw()
{
}

s8::Analyzer *TriggerController::createAnalyzer()
{
    return _analyzer.get();
}

bool TriggerController::inputFileShouldContinue()
{
    return false;
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
}

/**
 * MonitorTriggerController
 * s8
 *
 * Created by Samvel Khalatian on Feb 25, 2011
 * Copyright 2011, All rights reserved
 */

#include "Analyzer/interface/MonitorTriggerAnalyzer.h"
#include "Option/interface/MonitorOptions.h"

#include "Controller/interface/MonitorTriggerController.h"

using s8::MonitorTriggerController;

MonitorTriggerController::MonitorTriggerController() throw()
{
}

MonitorTriggerController::~MonitorTriggerController() throw()
{
}

// Override parent init()
//
void MonitorTriggerController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new MonitorTriggerAnalyzer());
    _analyzer->init();

    _options.reset(new MonitorOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

s8::Analyzer *MonitorTriggerController::createAnalyzer()
{
    return _analyzer.get();
}

core::Options *MonitorTriggerController::createOptions()
{
    return _options.get();
}

bool MonitorTriggerController::inputFileShouldLoadTriggers()
{
    return true;
}

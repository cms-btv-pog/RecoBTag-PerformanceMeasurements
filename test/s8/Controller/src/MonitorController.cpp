/**
 * MonitorController
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include "Analyzer/interface/MonitorAnalyzer.h"
#include "Option/interface/MonitorOptions.h"

#include "Controller/interface/MonitorController.h"

using s8::MonitorController;

MonitorController::MonitorController() throw()
{
}

MonitorController::~MonitorController() throw()
{
}

// Override parent init()
//
void MonitorController::init()
{
    // Let parent to initialize
    //
    AppController::init();

    _analyzer.reset(new MonitorAnalyzer());
    _analyzer->init();

    _options.reset(new MonitorOptions());
    _options->init();
    _options->setDelegate(_analyzer.get());
}

s8::Analyzer *MonitorController::createAnalyzer()
{
    return _analyzer.get();
}

core::Options *MonitorController::createOptions()
{
    return _options.get();
}

bool MonitorController::inputFileShouldLoadTriggers()
{
    return true;
}

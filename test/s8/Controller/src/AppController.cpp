/**
 * AppController
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <TFile.h>

#include "Analyzer/interface/Analyzer.h"
#include "IO/interface/File.h"
#include "S8Tree/interface/S8TreeInfo.h"
#include "S8Tree/interface/S8TriggerCenter.h"
#include "Redirector/interface/Debug.h"

#include "Controller/interface/AppController.h"

using std::clog;
using std::cout;
using std::endl;

using s8::AppController;

AppController::AppController() throw():
    _maxEvents(0),
    _processedEvents(0)
{
    _analyzer = 0;
}

AppController::~AppController() throw()
{
}

void AppController::init()
{
    _inputFile.reset(new InputFile());
    _inputFile->setDelegate(this);

    _outputFile.reset(new OutputFile());
    _outputFile->setDelegate(this);

    _debug.reset(new core::Debug());
}

void AppController::applicationWillRun()
{
    _startClock = clock();

    _analyzer = createAnalyzer();
}

void AppController::applicationDidRun()
{
    using std::setw;
    using std::left;

    clock_t stopClock = clock();

    cout << "Stats" << endl;
    cout << setw(25) << left << " - Run time"
        << (static_cast<double>(stopClock) - _startClock) / CLOCKS_PER_SEC
        << " (s)" << endl;

    cout << setw(25) << left << "     (per event)"
        << (static_cast<double>(stopClock) - _startClock) /
            _processedEvents / CLOCKS_PER_SEC
        << " (s)" << endl;

    cout << setw(25) << left << " - Processed events"
        << _processedEvents << endl;
}

void AppController::optionDebugIsSet(const std::string &fileName)
{
    _debug->init(fileName);
}

void AppController::optionEventsIsSet(const int &value)
{
    _maxEvents = value;
}


// No Options by default
//
core::Options *AppController::createOptions()
{
    return 0;
}

core::InputFile *AppController::createInputFile()
{
    return _inputFile.get();
}

core::OutputFile *AppController::createOutputFile()
{
    return _outputFile.get();
}



bool AppController::inputFileShouldOpen(const std::string &name)
{
    // Do not open files if analyzer is not created or max events are processed
    //
    return _analyzer ?
        (!_maxEvents ||
            _processedEvents < _maxEvents) :
        false;
}

void AppController::inputFileDidOpen(TFile *file)
{
    using std::runtime_error;

    clog << "<< " << file->GetName() << endl;

    TreeInfo *info = 0;
    file->GetObject("S8TreeMaker/s8info", info);

    if (!info)
        throw runtime_error("Failed to extract System8 Tree info");

    TriggerCenter *triggerCenter = 0;
    file->GetObject("S8TreeMaker/s8triggers", triggerCenter);

    if (!triggerCenter)
        throw runtime_error("Failed to extract System8 Trigger center");

    clog << " Tree " << info->version() << " with "
        << triggerCenter->triggers().size() << " triggers." << endl
        << endl;

    // call analyzer with tree info and trigger center
    //
    _analyzer->treeDidLoad(info, triggerCenter);
}

bool AppController::inputFileShouldLoadJets()
{
	return true;
}

bool AppController::inputFileShouldLoadMuons()
{
	return true;
}

bool AppController::inputFileShouldLoadPrimaryVertices()
{
	return true;
}

void AppController::inputFileDidLoadEvent(const s8::Event *event)
{
    ++_processedEvents;

    // pass event to analyzer
    //
    _analyzer->eventDidLoad(event);
}

bool AppController::inputFileShouldContinue()
{
    return !_maxEvents ||
            _processedEvents < _maxEvents;
}

void AppController::inputFileWillClose(TFile *file)
{
    // Do nothing
    //
}

void AppController::outputFileDidOpen(TFile *file)
{
    // Do nothing
    //
}

void AppController::outputFileWillClose(TFile *file)
{
    // let analyzer save results
    //
    if (file)
        _analyzer->save(file);
}

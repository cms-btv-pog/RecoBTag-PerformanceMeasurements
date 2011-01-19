/**
 * AppController
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>

#include "interface/File.h"
#include "interface/GenericOptions.h"

#include "interface/AppController.h"

using std::auto_ptr;
using std::cerr;
using std::cout;
using std::clog;
using std::endl;
using std::exception;
using std::string;

using core::AppController;

AppController::AppController() throw()
{
}

AppController::~AppController() throw()
{
}

bool AppController::init(int argc, char *argv[]) throw()
try
{
    // Let descendants initialize: an exception may be thrown
    //
    init();

    // Generic Options are always present
    //
    auto_ptr<GenericOptions> genericOptions(new GenericOptions());
    genericOptions->init();
    genericOptions->setDelegate(this);

    // Prepare combination of the Generic and Custom options
    //
    auto_ptr<po::options_description>
        programOptions(new po::options_description());
    programOptions->add(*genericOptions->description());

    // Create Custom options and add to the Generic Ones
    //
    Options *customOptions = createOptions();
    if (customOptions)
        programOptions->add(*customOptions->description());

    // Parse options
    //
    auto_ptr<po::variables_map> arguments(new po::variables_map());
    po::store(po::command_line_parser(argc, argv).
            options(*programOptions).positional(*genericOptions->positional()).
            run(),
            *arguments);

    // config file is temporarily broken
    //
    if (false &&
        arguments->count("config"))
    {
        const std::string config =
            arguments->operator[]("config").as<string>();

        /*
        po::store(po::parse_config_file(config.c_str(), *programOptions),
                  *arguments);
        */
    }

    // Exit if no arguments is supplied or --help is used
    if (1 == argc ||
        arguments->count("help"))
    {
        cout << "Usage: " << argv[0] << " [Options] inputs" << endl;
        cout << *programOptions << endl;

        return false;
    }

    // Notify observers with processed arguments
    //
    po::notify(*arguments);

    // Print options
    //
    cout << *genericOptions << endl;

    if (customOptions)
        cout << *customOptions << endl;
    
    return true;
}
catch(const std::exception &error)
{
    cerr << "[exception] " << error.what() << endl;

    return false;
}
catch(...)
{
    cerr << "[exception] unknown exception was raised" << endl;

    return false;
}

bool AppController::run() throw()
try
{
    if (_inputFileNames.empty())
    {
        cout << "No input is specified: nothing to be done" << endl;

        return true;
    }

    applicationWillRun();

    OutputFile *output = _outputFileName.empty() ? 0 : createOutputFile();
    if (output)
        output->setName(_outputFileName);

    if (output)
        output->open();

    const bool result = processInput();

    if (output)
        output->close();

    applicationDidRun();

    return result;
}
catch(const std::exception &error)
{
    cerr << "[exception] " << error.what() << endl;

    return false;
}
catch(...)
{
    cerr << "[exception] unknown exception was raised" << endl;

    return false;
}

void AppController::optionOutputIsSet(const std::string &fileName)
{
    if (!fileName.empty())
        _outputFileName = fileName;
    else
        cerr << "Output FileName is empty: no output is created" << endl;
}

void AppController::optionInputIsSet(const std::string &fileName)
{
    if (fileName.empty())
        cerr << "Input FileName is empty: skip" << endl;
    else
        _inputFileNames.push_back(fileName);
}

bool AppController::processInput() throw()
try
{
    for(InputFileNames::const_iterator filename = _inputFileNames.begin();
        _inputFileNames.end() != filename;
        ++filename)
    {
        InputFile *input = createInputFile();
        if (!input)
            // Nothing to be done
            //
            return false;

        input->setName(*filename);

        input->open();
        if (!input->isOpen())
            continue;

        input->process();

        input->close();
    }

    return true;
}
catch(const std::exception &error)
{
    cerr << "[exception] " << error.what() << endl;

    return false;
}
catch(...)
{
    cerr << "[exception] unknown exception was raised" << endl;

    return false;
}

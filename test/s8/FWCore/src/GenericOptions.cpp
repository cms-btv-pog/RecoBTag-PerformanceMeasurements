/**
 * GenericOptions
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <iostream>
#include <fstream>
#include <ostream>
#include <stdexcept>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "interface/GenericOptionsDelegate.h"

#include "interface/GenericOptions.h"

namespace sys = boost::filesystem;

using std::cerr;
using std::endl;
using std::runtime_error;
using std::string;
using std::vector;

using boost::regex;

using core::GenericOptions;
using core::GenericOptionsDelegate;

GenericOptions::GenericOptions() throw()
{
    _delegate = 0;
}

GenericOptions::~GenericOptions() throw()
{
}

void GenericOptions::init()
{
    _description.reset(new po::options_description("Generic Options"));
    _description->add_options()
        ("help,h",
             "Produce help message")

        ("config,c",
            po::value<string>(),
            "Config file with arguments (used in conjuction with command line args)")

        ("debug,d",
            po::value<string>()->default_value("debug.log")->notifier(
                boost::bind(&GenericOptions::optionDebugIsSet, this, _1)),
            "Save debug information in the output file")

        ("output,o",
            po::value<string>()->default_value("output.root")->notifier(
                boost::bind(&GenericOptions::optionOutputIsSet, this, _1)),
            "Output file")

        ("events,e",
            po::value<int>()->default_value(0)->notifier(
                boost::bind(&GenericOptions::optionEventsIsSet, this, _1)),
            "Maximum number of events to be processed: 0 - all")

        ("input,i",
            po::value<vector<string> >()->notifier(
                boost::bind(&GenericOptions::optionInputIsSet, this, _1)),
            "Input file(s)")
    ;

    _positional.reset(new po::positional_options_description());
    _positional->add("input", -1);
}

GenericOptionsDelegate *GenericOptions::delegate() const
{
    return _delegate;
}

void GenericOptions::setDelegate(GenericOptionsDelegate *delegate)
{
    _delegate = delegate;
}

po::options_description *GenericOptions::description() const
{
    return _description.get();
}

po::positional_options_description *GenericOptions::positional() const
{
    return _positional.get();
}

void GenericOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "Generic Options" << endl;
    out << setw(25) << left << " [+] Debug" <<
        (_debug.empty() ? "---" : _debug) << endl;

    out << setw(25) << left << " [+] Output" <<
        (_output.empty() ? "---" : _output) << endl;

    out << setw(25) << left << " [+] Events";
    if (_events)
        out << _events;
    else
        out << "all";
    out << endl;

    out << setw(25) << left << " [+] Inputs (" << _inputs.size()
        << "):" << endl;
    for(Files::const_iterator input = _inputs.begin();
        _inputs.end() != input;
        ++input)
    {
        out << "      " << (input->second ? " " : "x")
            << " " << input->first << endl;
    }
}



void GenericOptions::optionDebugIsSet(const string &value)
{
    _debug = value;

    // Forward message
    //
    if (_delegate)
        _delegate->optionDebugIsSet(value);
}

void GenericOptions::optionOutputIsSet(const string &value)
{
    _output = value;

    // Forward message
    //
    if (_delegate)
        _delegate->optionOutputIsSet(value);
}

void GenericOptions::optionEventsIsSet(const int &value)
{
    if (0 > value)
        throw runtime_error("Negative number of events is specified");

    _events = value;

    if (_delegate)
        _delegate->optionEventsIsSet(value);
}

void GenericOptions::optionInputIsSet(const vector<string> &value)
{
    using std::make_pair;

    if (value.empty())
        return;

    bool didPrintErrorMessage = false;

    for(vector<string>::const_iterator input = value.begin();
        value.end() != input;
        ++input)
    {
        // Restrict input files to ROOT or TXT (with list of input files)
        //
        if (!regex_search(*input, regex("\\.(txt|root)$")))
        {
            cerr << "Input '" << *input << "' is not TXT or ROOT" << endl;

            didPrintErrorMessage = true;

            _inputs.push_back(make_pair(*input, false));

            continue;
        }

        // Check input file for existance
        //
        if (!sys::exists(input->c_str()))
        {
            cerr << "Input '" << *input << "' does not exist" << endl;

            didPrintErrorMessage = true;

            _inputs.push_back(make_pair(*input, false));

            continue;
        }

        // Test if input is a TXT file with list of ROOTs
        //
        if (regex_search(*input, regex("\\.txt$")))
        {
            // Process input file
            //
            readInputFilesFromTXT(*input);
            
            continue;
        }

        _inputs.push_back(make_pair(*input, true));

        if (_delegate)
            _delegate->optionInputIsSet(*input);
    }

    if (didPrintErrorMessage)
        cerr << endl;
}

void GenericOptions::readInputFilesFromTXT(const std::string &filename)
{
    std::ifstream input(filename.c_str());
    if (!input.is_open())
    {
        cerr << "Failed to open '" << filename << "' to read: skip" << endl;

        return;
    }

    // Create vector of input files and pass it to optionInputIsSet as if
    // user has specified input files explicitly
    //
    vector<string> files;
    for(string file; input >> file; )
    {
        files.push_back(file);
    }

    optionInputIsSet(files);
}

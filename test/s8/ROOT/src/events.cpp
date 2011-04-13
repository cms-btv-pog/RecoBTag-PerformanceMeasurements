/**
 * Extract total number of events from tree(s)
 * s8
 *
 * Created by Samvel Khalatian on Mar 02, 2011
 * Copyright 2011, All rights reserved
 */

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>

#include <TFile.h>
#include <TTree.h>

using namespace std;

using boost::regex;
using boost::shared_ptr;
using boost::smatch;

namespace fs = boost::filesystem;

int extractEventsFromInput(const string &);
int extractEventsFromTXTFile(const string &filename);
int extractEventsFromFile(const string &filename);

int main(int argc, char *argv[])
try
{
    if (2 > argc)
    {
        cerr << "Usage: " << argv[0] << " input" << endl;

        return 1;
    }

    int total_events = 0;
    for(int i = 1; argc > i; ++i)
        total_events += extractEventsFromInput(argv[i]);

    cout << "Total Events: " << total_events << endl;
    
    return 0;
}
catch(const std::exception &error)
{
    cerr << "Error: " << error.what() << endl;

    return 1;
}
catch(...)
{
    cerr << "Unknown error" << endl;

    return 1;
}

int extractEventsFromInput(const string &filename)
{
    smatch matches;
    if (!regex_match(filename, matches, regex(".*\\.(?:(txt)|(root))$")))
    {
        cerr << "Input is not supported: " << filename << endl;

        return 0;
    }

    if (!fs::exists(filename.c_str()))
    {
        cerr << "Input does not exist: " << filename << endl;

        return 0;
    }

    // TXT File
    //
    if (matches[1].matched)
        return extractEventsFromTXTFile(filename);

    // ROOT file
    //
    else if (matches[2].matched)
        return extractEventsFromFile(filename);

    // Code should not get here
    //
    else
    {
        cerr << "Usupported file: " << filename << endl;

        return 0;
    }
}

int extractEventsFromTXTFile(const string &filename)
{
    ifstream txt_input(filename.c_str());
    if (!txt_input.is_open())
    {
        cerr << "Failed to open file: " << filename << endl;

        return 0;
    }

    int total_events = 0;
    for(string input; txt_input >> input; )
        total_events += extractEventsFromInput(input);

    return total_events;
}

int extractEventsFromFile(const string &filename)
{
    shared_ptr<TFile> input(new TFile(filename.c_str(), "read"));
    if (!input->IsOpen())
    {
        cerr << "Failed to open input file: " << filename << endl;

        return 0;
    }

    TTree *tree = dynamic_cast<TTree *>(input->Get("S8TreeMaker/s8"));
    if (!tree)
    {
        cerr << "Failed to exract Tree from file: " << filename << endl;
        
        return 0;
    }

    return tree->GetEntries();
}

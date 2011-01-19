/**
 * Generic Scale Factors calculator
 * s8
 *
 * Created by Samvel Khalatian on Dec 30, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/regex.hpp>

#include <TFile.h>
#include <TGraphErrors.h>

#include "Stat/interface/Measurement.h"

using std::auto_ptr;
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::runtime_error;
using std::string;

using boost::regex;
using boost::smatch;

enum Type { BTAG, CLTAG };

string getGraphName(const Type &type);

ostream &operator<<(ostream &, const TGraphErrors &);

void extractScaleFactors(const TGraphErrors *, const TGraphErrors *);

class Input
{
    public:
        Input() throw();

        void init(const string &);
        void read();

        TGraphErrors *graph() const;

    private:
        string  _filename;
        string  _path;
        string  _graphname;

        TFile        *_file;
        TGraphErrors *_graph;
};

int main(int argc, char *argv[])
try
{
    if (3 != argc)
        throw runtime_error("usage: " + string(argv[0]) +
                            " numerator.root:path/plot denominator.root:path/plot");

    auto_ptr<Input> numerator(new Input());
    numerator->init(argv[1]);
    numerator->read();

    auto_ptr<Input> denominator(new Input());
    denominator->init(argv[2]);
    denominator->read();

    cout << "Numerator:" << endl;
    cout << *(numerator->graph()) << endl;

    cout << "Denominator:" << endl;
    cout << *(denominator->graph()) << endl;

    extractScaleFactors(numerator->graph(), denominator->graph());

    return 0;
}
catch(const std::exception &error)
{
    cerr << "Error" << endl;
    cerr << error.what() << endl;

    return 1;
}
catch(...)
{
    cerr << "Unknown error" << endl;

    return 1;
}

Input::Input() throw()
{
    _file = 0;
    _graph = 0;
}

void Input::init(const string &input)
{
    using std::setw;
    using std::left;

    if (_file)
    {
        _graph = 0;

        _file->Close();
        _file = 0;
    }

    // Parse input
    //
    smatch matches;
    regex pattern("^((?:\\w+/)*\\w+\\.root):(?:((?:\\w+/)+)(\\w+)|(\\w+))$");
    if (!regex_match(input, matches, pattern))
        throw runtime_error("Supported Input format: file.root:path/graph");

    if (matches[2].matched)
    {
        _path = matches[2];
        _graphname = matches[3];
    }
    else if (matches[4].matched)
    {
        _path = "";
        _graphname = matches[4];
    }

    _filename = matches[1];

    cout << "Input is parsed" << endl;
    cout << setw(25) << left << " [+] File" << _filename << endl;
    cout << setw(25) << left << " [+] Path" << (_path.empty() ? "---" : _path) << endl;
    cout << setw(25) << left << " [+] Graph" << _graphname << endl;
    cout << endl;
}

void Input::read()
{
    if (_filename.empty())
        throw runtime_error("Input is not initialized");

    if (!_file)
    {
        _file = TFile::Open(_filename.c_str());
        if (!_file->IsOpen())
            throw runtime_error("Failed to open input: " + _filename);
    }

    string object = (_path.empty() ? "" : _path) + _graphname;

    _graph = dynamic_cast<TGraphErrors *>(_file->Get(object.c_str()));
    if (!_graph)
        throw runtime_error("Failed to get object [" + object +
                            "] from input: " + _filename);
}

TGraphErrors *Input::graph() const
{
    return _graph;
}

string getGraphName(const Type &type)
{
    switch(type)
    {
        case Type(BTAG):  return "eff_tag_b";
        case Type(CLTAG): return "eff_tag_cl";
    }

    throw runtime_error("Unsupported Graph type supplied");
}

void extractScaleFactors(const TGraphErrors *s8, const TGraphErrors *mc)
{
    using s8::Measurement;

    for(int s8Point = 0, mcPoint = 0; s8->GetN() > s8Point; ++s8Point)
    {
        double s8x;
        double s8y;

        s8->GetPoint(s8Point, s8x, s8y);

        Measurement s8Efficiency(s8y, pow(s8->GetErrorY(s8Point), 2));

        // Find Systematic
        //
        double mcX = 0;
        double mcY;

        while(mc->GetN() > mcPoint)
        {
            mc->GetPoint(mcPoint, mcX, mcY);

            if (mcX < s8x)
            {
                ++mcPoint;

                continue;
            }

            if (mcX > s8x)
                mcX = -1;

            break;
        }

        // Check if more points are left in the MC Graph
        //
        if (!mcX)
            break;

        // Is point missing?
        //
        if (-1 == mcX)
            continue;

        Measurement mcEfficiency(mcY, pow(mc->GetErrorY(mcPoint), 2));
        
        Measurement scaleFactor = s8Efficiency % mcEfficiency;

        cout << "bin " << (s8Point + 1) << ": " << scaleFactor << endl;
    }
}

ostream &operator<<(ostream &out, const TGraphErrors &graph)
{
    using s8::Measurement;

    for(int point = 0; graph.GetN() > point; ++point)
    {
        double x;
        double y;

        graph.GetPoint(point, x, y);

        Measurement efficiency(y, pow(graph.GetErrorY(point), 2));

        out << "bin " << (point + 1) << ": " << efficiency << endl;
    }

    return out;
}

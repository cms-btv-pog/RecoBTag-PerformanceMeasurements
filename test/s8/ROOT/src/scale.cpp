/**
 * Calculate Scale Factors
 * s8
 *
 * Created by Samvel Khalatian on Dec 20, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TGraphErrors.h>

#include "Stat/interface/Measurement.h"

using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using std::string;

enum Type { BTAG, CLTAG };

string getGraphName(const Type &type);

void processEfficiency(TFile *input, const Type &type);

void extractEfficiencies(const TGraphErrors *s8);
void extractScaleFactors(const TGraphErrors *, const TGraphErrors *);

int main(int argc, char *argv[])
try
{
    TFile *input = TFile::Open(argv[1]);
    if (!input->IsOpen())
        throw runtime_error("Failed to open input file");

    processEfficiency(input, Type(BTAG));
    processEfficiency(input, Type(CLTAG));

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

string getGraphName(const Type &type)
{
    switch(type)
    {
        case Type(BTAG):  return "eff_tag_b";
        case Type(CLTAG): return "eff_tag_cl";
    }

    throw runtime_error("Unsupported Graph type supplied");
}

void processEfficiency(TFile *input, const Type &type)
{
    string graphName = getGraphName(type);

    TGraphErrors *s8efficiency =
        dynamic_cast<TGraphErrors *>(input->Get(("s8efficiency/" +
                                                graphName).c_str()));
    if (!s8efficiency)
        throw runtime_error("Failed to extract System8 efficiency: " +
                            graphName);

    TGraphErrors *mcefficiency =
        dynamic_cast<TGraphErrors *>(input->Get(("mcefficiency/" +
                                                graphName).c_str()));
    if (!mcefficiency)
        throw runtime_error("Failed to extract Monte-Carlo efficiency: " +
                            graphName);

    string typeName = "b-Tag";
    if (Type(CLTAG) == type)
        typeName = "cl-Tag";

    cout << typeName << " efficiencies:" << endl;
    extractEfficiencies(s8efficiency);
    cout << endl;

    cout << typeName << " scale factors:" << endl;
    extractScaleFactors(s8efficiency, mcefficiency);
    cout << endl;
}

void extractEfficiencies(const TGraphErrors *s8)
{
    using s8::Measurement;

    for(int s8Point = 0; s8->GetN() > s8Point; ++s8Point)
    {
        double s8x;
        double s8y;

        s8->GetPoint(s8Point, s8x, s8y);

        Measurement s8Efficiency(s8y, pow(s8->GetErrorY(s8Point), 2));

        cout << "bin " << (s8Point + 1) << ": " << s8Efficiency << endl;
    }
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
        
        Measurement scaleFactor = s8Efficiency / mcEfficiency;

        cout << "bin " << (s8Point + 1) << ": " << scaleFactor << endl;
    }
}

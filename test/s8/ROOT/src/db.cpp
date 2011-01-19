/**
 * Generate DB ASCII file
 * s8
 *
 * Created by Samvel Khalatian on Jan 04, 2011
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TGraphErrors.h>

#include "S8Tree/interface/S8Jet.h"
#include "Stat/interface/Measurement.h"
#include "Utility/interface/TaggerOperatingPoint.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;
using std::map;
using std::ostream;
using std::runtime_error;
using std::string;
using std::vector;

using s8::Jet;
using s8::Measurement;
using s8::TaggerOperatingPoint;

struct Result
{
    Measurement point;
    Measurement eff_tag_b;
    Measurement sf_eff_tag_b;
};

enum ResultType
{
    //
    // BTAG
    //
    BTAGBEFF=1001, BTAGBERR=1002, BTAGCEFF=1003, 
    BTAGCERR=1004, BTAGLEFF=1005, BTAGLERR=1006, BTAGNBEFF=1007, BTAGNBERR=1008,
    //
    // add corrections in case the table is for weights and not efficiencies
    //
    BTAGBEFFCORR=1009, BTAGBERRCORR=1010, BTAGCEFFCORR=1011, 
    BTAGCERRCORR=1012, BTAGLEFFCORR=1013, BTAGLERRCORR=1014, BTAGNBEFFCORR=1015, BTAGNBERRCORR=1016,
    //
    // MUONS
    //
    MUEFF=2001, MUERR=2002, MUFAKE=2003, MUEFAKE=2004
};

enum  BinningVariablesType
{
    // Jets
    JetEta=1, JetEt=2, JetPhi=3, JetNTracks=4, JetAbsEta=5,

    // Muons 
    MuonPt=1001, MuonCharge=1002,MuonEta=1003, MuonPhi=1004
};

typedef map<double, Measurement> BinnedSystematics;

void processEfficiency(TFile *input, const TaggerOperatingPoint &);
vector<Result> extractScaleFactors(const TGraphErrors *, const TGraphErrors *,
                                   BinnedSystematics &);

BinnedSystematics getSystematics(const TaggerOperatingPoint &op);

ostream &operator<<(ostream &, const Jet::BTag &);
ostream &operator<<(ostream &, const Measurement &);

ostream &printMinMax(ostream &out, const Measurement &m, const double & = 0);

template<typename T>
    ostream &operator<<(ostream &out, const vector<T> &array)
    {
        typename vector<T>::const_iterator element = array.begin();

        while(true)
        {
            out << *element;

            if (array.end() == ++element)
                break;

            out << " ";
        }

        return out;
    }

int main(int argc, char *argv[])
try
{
    if (3 != argc)
        throw runtime_error("Usage: " + string(argv[0]) +
                            " OperatingPoint s8results.root");

    TaggerOperatingPoint operatingPoint;
    operatingPoint << argv[1];

    TFile *input = TFile::Open(argv[2]);
    if (!input->IsOpen())
        throw runtime_error("Failed to open input file: " + string(argv[2]));

    vector<ResultType> results;
    results.push_back(ResultType(1001));
    results.push_back(ResultType(1002));
    results.push_back(ResultType(1009));
    results.push_back(ResultType(1010));

    vector<BinningVariablesType> types;
    types.push_back(BinningVariablesType(5));
    types.push_back(BinningVariablesType(2));

    cout << operatingPoint.btag() << endl;
    cout << operatingPoint << endl;
    cout << "PerformancePayloadFromTable" << endl;
    cout << results.size() << endl;
    cout << types.size() << endl;
    cout << results << endl;
    cout << types << endl;

    processEfficiency(input, operatingPoint);

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

ostream &operator <<(ostream &out, const Jet::BTag &btag)
{
    switch(btag)
    {
        case Jet::TCHE:  out << "TrackCountingHighEfficiency";
                         break;

        case Jet::TCHP:  out << "TrackCountingHighPurity";
                         break;

        case Jet::JP:    out << "JetProbability";
                         break;

        case Jet::SSVHE: out << "SimpleSecondaryVertexHighEfficiency";
                         break;

        case Jet::SSVHP: out << "SimpleSecondaryVertexHighPurity";
                         break;

        default:         throw runtime_error("Unsupported BTag is used. Only TCHE(P), JP, SSVHE(P) can be used");
    }

    return out;
}

ostream &operator <<(ostream &out, const Measurement &m)
{
    return out << m.value() << " " << sqrt(m.variance());
}

ostream &printMinMax(ostream &out, const Measurement &m, const double &max)
{
    double sigma = sqrt(m.variance());

    return out << (m.value() - sigma) << " " << ( 0 == max ? m.value() + sigma : max);
}

void processEfficiency(TFile *input, const TaggerOperatingPoint &operatingPoint)
{
    string graphName = "eff_tag_b";

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

    BinnedSystematics systematics = getSystematics(operatingPoint);
    
    typedef vector<Result> Results;
    Results dbResults = extractScaleFactors(s8efficiency, mcefficiency,
                                            systematics);

    for(Results::iterator result = dbResults.begin();
        dbResults.end() != result;
        ++result)
    {
        cout << "0 2.4 ";

        if (1 == distance(result, dbResults.end()))
            printMinMax(cout, result->point, 9999) << " ";
        else
            printMinMax(cout, result->point) << " ";

        ::operator<<(cout, result->eff_tag_b) << " ";
        ::operator<<(cout, result->sf_eff_tag_b) << endl; 
    }
}

vector<Result> extractScaleFactors(const TGraphErrors *s8,
                                   const TGraphErrors *mc,
                                   BinnedSystematics &systematics)
{
    vector<Result> results;

    for(int s8Point = 0, mcPoint = 0; s8->GetN() > s8Point; ++s8Point)
    {
        double s8x;
        double s8y;

        s8->GetPoint(s8Point, s8x, s8y);

        Measurement s8point(s8x, pow(s8->GetErrorX(s8Point), 2));
        Measurement s8Efficiency(s8y, pow(s8->GetErrorY(s8Point), 2));
        
        // Add systematics
        //
        Measurement s8systematics = systematics[s8x];

        // Convert relative systematic error to absolute
        //
        s8systematics.setVariance(s8systematics.variance() * pow(s8Efficiency.value(), 2));

        // Add systematic error to the statistical in the efficiency
        //
        s8Efficiency += s8systematics;

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

        Result result;
        result.point = s8point;
        result.eff_tag_b = s8Efficiency;
        result.sf_eff_tag_b = scaleFactor;

        results.push_back(result);
    }

    return results;
}

// Note: errors are relative (not absolute)
//
BinnedSystematics getSystematics(const TaggerOperatingPoint &op)
{
    BinnedSystematics systematics;

    if (TaggerOperatingPoint::TCHEL == op)
    {
        systematics[25] = Measurement(0, pow(0.083, 2));
        systematics[35] = Measurement(0, pow(0.038, 2));
        systematics[45] = Measurement(0, pow(0.053, 2));
        systematics[55] = Measurement(0, pow(0.051, 2));
        systematics[70] = Measurement(0, pow(0.056, 2));
        systematics[110] = Measurement(0, pow(0.075, 2));
        systematics[185] = Measurement(0, pow(0.175, 2));
    }
    else if (TaggerOperatingPoint::TCHEM == op)
    {
        systematics[25] = Measurement(0, pow(0.092, 2));
        systematics[35] = Measurement(0, pow(0.050, 2));
        systematics[45] = Measurement(0, pow(0.055, 2));
        systematics[55] = Measurement(0, pow(0.059, 2));
        systematics[70] = Measurement(0, pow(0.101, 2));
        systematics[110] = Measurement(0, pow(0.144, 2));
        systematics[185] = Measurement(0, pow(0.328, 2));
    }
    else if (TaggerOperatingPoint::SSVHEM == op)
    {
        systematics[25] = Measurement(0, pow(0.194, 2));
        systematics[35] = Measurement(0, pow(0.049, 2));
        systematics[45] = Measurement(0, pow(0.052, 2));
        systematics[55] = Measurement(0, pow(0.063, 2));
        systematics[70] = Measurement(0, pow(0.091, 2));
        systematics[110] = Measurement(0, pow(0.153, 2));
        systematics[185] = Measurement(0, pow(0.516, 2));
    }
    else
    {
        systematics[25] = Measurement(0, 0);
        systematics[35] = Measurement(0, 0);
        systematics[45] = Measurement(0, 0);
        systematics[55] = Measurement(0, 0);
        systematics[70] = Measurement(0, 0);
        systematics[110] = Measurement(0, 0);
        systematics[185] = Measurement(0, 0);
    }

    return systematics;
}

/**
 * Reweight triggered evetns
 * s8
 *
 * Created by Samvel Khalatian on Mar 01, 2011
 * Copyright 2011, All rights reserved
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include <TFile.h>
#include <TH3.h>

#include "Utility/interface/Range.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::runtime_error;
using std::string;

using s8::Range;

typedef pair<int, int> Bins;

struct BinsGroup
{
    BinsGroup(const TH3 *hist):
        x(hist->GetXaxis()->GetFirst(), hist->GetXaxis()->GetLast()),
        y(hist->GetYaxis()->GetFirst(), hist->GetYaxis()->GetLast()),
        z(hist->GetZaxis()->GetFirst(), hist->GetZaxis()->GetLast())
    {}

    Bins x;
    Bins y;
    Bins z;
};

void calculateWeights(const Range &, TFile *out, TDirectory *data, TDirectory *mc);

int main(int argc, char *argv[])
try
{
    if (5 > argc)
    {
        cout << "Usage: " << argv[0] << " MIN..MAX out.root data_monitor.root"
            << " mc_monitor.root" << endl;

        return 1;
    }

    Range range;
    parse(range, argv[1]);

    TFile *out = TFile::Open(argv[2], "recreate");
    if (!out->IsOpen())
        throw runtime_error("Failed to open output file");

    TFile *data = TFile::Open(argv[3], "read");
    if (!data->IsOpen())
        throw runtime_error("Failed to open Data");

    TFile *mc = TFile::Open(argv[4], "read");
    if (!mc->IsOpen())
        throw runtime_error("Failed to open Monte-Carlo");

    calculateWeights(range, out, data, mc);

    out->Close();

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

void calculateWeights(const Range &range, TFile *out, TDirectory *data, TDirectory *mc)
{
    const string path = "MonitorAnalyzer/n/n_pt_eta_pv";

    TH3 *h_data = dynamic_cast<TH3 *>(data->Get(path.c_str()));
    if (!h_data)
        throw runtime_error("Failed to extract " + path + " from Data");

    TH3 *h_mc = dynamic_cast<TH3 *>(mc->Get(path.c_str()));
    if (!h_mc)
        throw runtime_error("Failed to extract " + path + " from Monte-Carlo");

    BinsGroup data_bins(h_data);
    BinsGroup mc_bins(h_mc);

    if (range.minimum())
    {
        data_bins.x.first = h_data->GetXaxis()->FindBin(range.minimum());
        mc_bins.x.first = h_mc->GetXaxis()->FindBin(range.minimum());
    }

    if (range.maximum())
    {
        data_bins.x.second = h_data->GetXaxis()->FindBin(range.maximum());
        mc_bins.x.second = h_mc->GetXaxis()->FindBin(range.maximum());
    }

    h_mc->Scale(h_data->Integral(data_bins.x.first, data_bins.x.second,
                                 data_bins.y.first, data_bins.y.second,
                                 data_bins.z.first, data_bins.z.second) /
                h_mc->Integral(mc_bins.x.first, mc_bins.x.second,
                               mc_bins.y.first, mc_bins.y.second,
                               mc_bins.z.first, mc_bins.z.second));

    TH3 *h_weights = dynamic_cast<TH3 *>(h_data->Clone());
    h_weights->Divide(h_mc);

    out->WriteTObject(h_weights);
}

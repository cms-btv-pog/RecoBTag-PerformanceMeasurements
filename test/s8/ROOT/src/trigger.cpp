/**
 * Compare MC and Data processed with different triggers
 * s8
 *
 * Created by Samvel Khalatian on Feb 07, 2011
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TRint.h>

using namespace std;

typedef pair<TDirectory *, string> Input;
typedef vector<Input> Inputs;
typedef stack<TObject *> Heap;

Heap heaps;

void compare(TDirectory *mc, const Inputs &data);

int main(int argc, char *argv[])
try
{
    if (3 > argc)
    {
        cerr << "Usage: " << argv[0] << " mc.root data.root [data.root]"
            << endl;

        return 1;
    }

    std::auto_ptr<TRint> application(new TRint("trigger", 0, 0));

    TFile *mc = TFile::Open(argv[1]);
    if (!mc->IsOpen())
        throw runtime_error("Failed to open MC file");

    Inputs data;
    for(int i = 2; argc > i; ++i)
    {
        TFile *input = TFile::Open(argv[i]);
        if (!input->IsOpen())
        {
            cerr << "Failed to open Data file: " << argv[i] << endl;

            continue;
        }

        data.push_back(make_pair(input, argv[i]));
    }

    if (data.empty())
    {
        cerr << "No input data is available: exit" << endl;

        return 0;
    }

    compare(mc, data);

    application->Run(true);

    cout << endl;

    while(!heaps.empty())
    {
        delete heaps.top();

        heaps.pop();
    }
    
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

void compare(TDirectory *mc, const Inputs &data)
{
    string plot_name = "MonitorAnalyzer/n/njet_pt";
    //string plot_name = "MonitorAnalyzer/generic/leadingjet_pt";

    // Extract and normalize MC plot
    //
    int color = kRed + 1;

    TH1 *plot_mc = dynamic_cast<TH1 *>(mc->Get(plot_name.c_str())->Clone());
    if (!plot_mc)
    {
        cerr << "Failed to extract MC plot: " << plot_name << endl;

        return;
    }
    plot_mc->SetLineColor(color);
    plot_mc->SetMarkerColor(color);
    plot_mc->GetYaxis()->SetTitle("Events");
    plot_mc->GetYaxis()->SetTitleOffset(1.5);

    double plot_mc_integral = plot_mc->Integral();

    // Process each Data input
    //
    for(Inputs::const_iterator input = data.begin();
        data.end() != input;
        ++input)
    {
        TH1 *plot_data = dynamic_cast<TH1 *>(input->first->Get(plot_name.c_str())->Clone());
        if (!plot_data)
        {
            cerr << "Failed to extract Data plot: " << plot_name << endl;
            cerr << "Input: " << input->second << endl;

            continue;
        }

        TCanvas *canvas = new TCanvas(input->second.c_str(),
            input->second.c_str(), 640, 480);
        heaps.push(canvas);

        canvas->SetGrid();

        TH1 *plot_mc_clone = dynamic_cast<TH1 *>(plot_mc->Clone());
        plot_mc_clone->Scale(plot_data->Integral() / plot_mc_integral);

        plot_mc_clone->Draw("h p e1");
        plot_data->Draw("same ap e1");

        canvas->RedrawAxis();
        canvas->Update();

        std::ostringstream filename;
        filename << input->second << ".png";

        canvas->SaveAs(filename.str().c_str());
    }
}

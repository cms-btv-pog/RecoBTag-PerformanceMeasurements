#include <cmath>

#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include <TCanvas.h>
#include <TClass.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <THStack.h>
#include <TKey.h>
#include <TLatex.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TObject.h>
#include <TRint.h>
#include <TLegend.h>

namespace sys = boost::filesystem;

using std::cerr;
using std::clog;
using std::cout;
using std::endl;
using std::runtime_error;
using std::string;

double luminosity = 33.3;

enum Type { Normal, Closure, Gluon };

typedef std::pair<double, double> Measurement;

std::ostream &operator<<(std::ostream &out,
                         const Measurement &measurement)
{
    const double sigma = sqrt(measurement.second);

    return out << std::setprecision(4) << std::fixed
        << measurement.first << " +/- " << sigma
        << " ("  << 100 * sigma / measurement.first << "%)";
}

typedef std::stack<TObject *> Heap;

Heap heaps;

void systematics(TDirectory *target, TList *sourcelist, const Type &type = Type(Normal) );
std::string style(TGraphErrors *graph, TFile *file);
TLegend *createLegend(const string &name = "");

void style(TLegend *legend)
{ 
    legend->SetBorderSize(1);
    legend->SetLineStyle(1);
    legend->SetTextFont(43);
    legend->SetTextSizePixels(24);
    legend->SetFillColor(0);
}

TLatex *createLabel(const double &luminosity = 33.3)
{
    std::ostringstream title;
    title << "#splitline{CMS Preliminary 2010}{";

    if (luminosity)
        title << static_cast<unsigned int>(ceil(luminosity)) << " pb^{-1}";
    else
        title << "Simulation";

    title << " at #sqrt{s} = 7 TeV}";

    TLatex *label = new TLatex(3.570061, 23.08044, title.str().c_str());
    label->SetNDC();
    label->SetTextAlign(13);
    label->SetX(0.4);
    label->SetY(0.908);

    return label;
}

int main(int argc, char *argv[])
try
{
    // Test if sufficient number of arguments is specified.
    if (4 > argc)
        throw std::invalid_argument("usage: merge out.root central.root systematics1.root [systematics2.root]");

    std::auto_ptr<TRint> application(new TRint("histInMemory", 0, 0));

    TFile *output = 0;
    TList *inputs = 0;
    try
    {
        if (sys::exists(argv[1]))
            throw runtime_error("Output file exists");

        // Create output file
        output = new TFile(argv[1], "RECREATE");
        if (!output->IsOpen())
            throw std::runtime_error("Failed to open output file");

        // open specified input files to the List
        //
        inputs = new TList();
        cout << "Inputs" << endl;
        for( int i = 2; argc > i; ++i)
        {
            cout << " [+] " << argv[i] << endl;
            inputs->Add(TFile::Open(argv[i]));
        }

        // Call merge
        //
        if ("closure.root" == string(argv[1]))
            systematics(output, inputs, Type(Closure));
        else if ("gsplit.root" == string(argv[1]))
            systematics(output, inputs, Type(Gluon));
        else
            systematics(output, inputs);

        // memory cleanup
        delete inputs;
        delete output;
    }
    catch(const std::exception &error)
    {
        // memory cleanup in case of error
        if (inputs)
            delete inputs;

        if (output)
            delete output;

        throw;
    }

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
    cerr << error.what() << endl;

    return 1;
}

void systematics(TDirectory *target, TList *sourcelist, const Type &type)
{
    using std::make_pair;

    TCanvas *canvas = new TCanvas();
    heaps.push(canvas);

    canvas->SetGrid();

    TMultiGraph *graph = new TMultiGraph();
    heaps.push(graph);

    TLegend *legend = new TLegend(0.57,0.22,0.97,0.38);
    heaps.push(legend);

    TLatex *label = (Type(Normal) == type ? createLabel() : createLabel(0));
    heaps.push(label);

    style(legend);

    TFile *first_source = (TFile*)sourcelist->First();
    TGraphErrors *central = dynamic_cast<TGraphErrors *>(first_source->Get("s8efficiency/eff_tag_b")->Clone());
    legend->AddEntry(central, style(central, first_source).c_str(), "p");

    graph->Add(central, "p");

    if (!central)
    {
        cerr << "Failed to extract Eff_tag_b central value" << endl;

        return;
    }

    // Graph, point counter
    //
    typedef std::pair<TGraphErrors *, int> Graph;
    typedef std::vector<Graph> Graphs;

    Graphs graphs;

    if (Type(Closure) == type)
    {
        TGraphErrors *mc = dynamic_cast<TGraphErrors *>(first_source->Get("mcefficiency/eff_tag_b")->Clone());
        legend->AddEntry(mc, "MCTrue", "p");

        graph->Add(mc, "p");
        graphs.push_back(make_pair(mc, 0));
    }
    else
    {
        for(TFile *nextsource = (TFile*)sourcelist->After(first_source);
            nextsource;
            nextsource = (TFile*)sourcelist->After(nextsource))
        {
            TGraphErrors *eff = dynamic_cast<TGraphErrors *>(nextsource->Get("s8efficiency/eff_tag_b")->Clone());
            legend->AddEntry(eff, style(eff, nextsource).c_str(), "p");

            graph->Add(eff, "p");
            graphs.push_back(make_pair(eff, 0));
        }
    }

    for(int point = 0; central->GetN() > point; ++point)
    {
        double centralX;
        double y;

        central->GetPoint(point, centralX, y);

        Measurement measurement = make_pair(y, 0);

        // Find Systematic
        //
        for(Graphs::iterator iter = graphs.begin();
            graphs.end() != iter;
            ++iter)
        {
            double x = 0;

            while(iter->first->GetN() > iter->second)
            {
                iter->first->GetPoint(iter->second, x, y);

                if (x < centralX)
                {
                    ++(iter->second);

                    x = 0;

                    continue;
                }

                if (x > centralX)
                    x = -1;

                break;
            }

            if (!x)
            {
                clog << "No more points in graph" << endl;

                continue;
            }

            if (-1 == x)
            {
                clog << "missing point for Central x = " << centralX
                    << " (x = " << x << ")" << endl;

                continue;
            }
            
            const double sigma = pow(y - measurement.first, 2);
            if (sigma > measurement.second)
                measurement.second = sigma;
        }

        cout << "bin " << (point + 1) << ": " << measurement << endl;
    }

    graph->Draw("a");
    graph->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
    graph->GetYaxis()->SetTitle("#epsilon_{b}^{tag}");
    graph->GetXaxis()->SetNdivisions(6);
    graph->GetYaxis()->SetNdivisions(6);
    graph->Draw("a");

    legend->Draw();

    label->Draw();
}

std::string style(TGraphErrors *graph, TFile *file)
{
    const std::string filename = file->GetName();
    std::string label = "";
    int color = 0;
    
    if ("ttbar.root" == filename)
    {
        color = kBlack;
        label = "TTbarJets";
    }

    else if ("qcd.root" == filename)
    {
        color = kGreen + 1;
        label = "QCD";
    }

    else if ("pt50.root" == filename)
    {
        color = kRed + 1;
        label = "Pt50";
    }

    else if ("pt30.root" == filename)
    {
        color = kGreen + 1;
        label = "Pt30";
    }

    else if ("ptrel05.root" == filename)
    {
        color = kOrange - 3;
        label = "p_{T}^{rel} > 0.5";
    }

    else if ("ptrel08.root" == filename)
    {
        color = kRed;
        label = "p_{T}^{rel} > 0.8";
    }

    else if ("ptrel12.root" == filename)
    {
        color = kRed + 3;
        label = "p_{T}^{rel} > 1.2";
    }

    else if ("gluon.root" == filename)
    {
        color = kBlue + 1;
        label = "Gluon-Splitting";
    }

    else if ("awayTCHEL.root" == filename)
    {
        color = kAzure + 2;
        label = "Away TCHEL";
    }

    else if ("awayTCHPM.root" == filename)
    {
        color = kAzure + 6;
        label = "Away TCHPM";
    }

    else if ("mu5.root" == filename)
    {
        color = kOrange + 1;
        label = "#mu p_{T} > 5";
    }

    else if ("mu7.root" == filename)
    {
        color = kOrange + 2;
        label = "#mu p_{T} > 7";
    }

    else if ("mu10.root" == filename)
    {
        color = kOrange + 3;
        label = "#mu p_{T} > 10";
    }

    else if ("mc.root" == filename)
    {
        color = kGreen + 1;
        label = "System8";
    }

    else
    {
        cerr << "Can not aply style to plot: unrecognized filename - " << filename << endl;

        return "";
    }

    graph->SetLineColor(color);
    graph->SetMarkerColor(color);

    return label;
}

TLegend *createLegend(const string &name)
{
    TLegend *legend = new TLegend( .6, .6, .8, .9);
    if (!name.empty())
        legend->SetHeader(name.c_str());

    legend->SetMargin(0.12);
    legend->SetTextSize(0.035);
    legend->SetFillColor(10);
    legend->SetBorderSize(0);

    return legend;
}


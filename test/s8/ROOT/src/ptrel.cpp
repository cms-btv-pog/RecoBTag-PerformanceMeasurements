#include <cstdlib>

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TClass.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TKey.h>
#include <TList.h>
#include <TRint.h>
#include <TLegend.h>

using std::auto_ptr;
using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using std::string;
using std::vector;
using std::invalid_argument;
using std::exception;
using std::ostringstream;

double luminosity = 31.6;
vector<TObject *> heap;

class Input
{
    public:
        enum Type { Data, ppMuX, TTbarJets, Pt15, Pt30, Pt50, Pt150 };

        Input()
        {
            file = 0;
            type = Data;
        }

        TFile *file;
        Type   type;
};

void cleanUpMemory();
void merge(Input &, Input &);
TH1 *getProjection(TH2 *hist, const int &bin);
void printStats(const string &name, TH1 *plot);

std::string style(TFile *file, TH1 *);
TLegend *createLegend(const string &name = "");

double ptrelCut = 0.8;

int main(int argc, char *argv[])
try
{
    // Test if sufficient number of arguments is specified.
    if (4 != argc)
        throw invalid_argument("usage: " + string(argv[0]) + " PtRelCut data.root mc.root");

    auto_ptr<TRint> application(new TRint("ptrel", 0, 0));

    try
    {
        ::ptrelCut = atof(argv[1]);
        // Open DATA
        //
        Input data;
        data.file = TFile::Open(argv[2]);
        if (!data.file)
            throw runtime_error("Failed to open Data file");

        // Open MC
        //
        Input mc;
        mc.file = TFile::Open(argv[3]);
        if (!mc.file)
            throw runtime_error("Failed to open MC file");

        // Call merge
        //
        merge(data, mc);

        application->Run(true);
    }
    catch(const exception &error)
    {
        cleanUpMemory();

        throw;
    }

    // Explicit memory clean up
    //
    cleanUpMemory();

    return 0;
}
catch(const exception &error)
{
    cerr << "Error" << endl;
    cerr << error.what() << endl;

    return 1;
}

void cleanUpMemory()
{
    for(vector<TObject *>::iterator obj = heap.begin();
        heap.end() != obj;
        ++obj)
    {
        delete *obj;
    }
}

void merge(Input &data, Input &mc)
{
    TH2 *ptrelData = dynamic_cast<TH2 *>(data.file->Get("PtRelAnalyzer/ptrel"));
    TH2 *ptrel_b = dynamic_cast<TH2 *>(mc.file->Get("PtRelAnalyzer/ptrel_b"));
    TH2 *ptrel_cl = dynamic_cast<TH2 *>(mc.file->Get("PtRelAnalyzer/ptrel_cl"));

    cout << "PtRel cut at: " << ::ptrelCut << endl;

    TCanvas *canvas = new TCanvas("merge", "Merge", 1200, 768);
    heap.push_back(canvas);
    canvas->Divide(4, 2);

    for(int bin = 1, bins = ptrelData->GetNbinsY();
        bins >= bin;
        ++bin)
    {
        canvas->cd(bin);

        ostringstream title;
        title << "pT: " << ptrelData->GetYaxis()->GetBinLowEdge(bin)
            << ".." << ptrelData->GetYaxis()->GetBinUpEdge(bin);
        TLegend *legend = createLegend(title.str().c_str());
        heap.push_back(legend);

        TH1 *plotData = getProjection(ptrelData, bin);
        heap.push_back(plotData);
        legend->AddEntry(plotData, "Data");
        plotData->SetLineWidth(2);

        const double dataIntegral = plotData->Integral();

        TH1 *plot_b = getProjection(ptrel_b, bin);
        heap.push_back(plot_b);
        legend->AddEntry(plot_b, "b-jets");
        plot_b->SetLineWidth(2);
        plot_b->SetLineColor(kRed);

        TH1 *plot_cl = getProjection(ptrel_cl, bin);
        heap.push_back(plot_cl);
        legend->AddEntry(plot_cl, "cl-jets");
        plot_cl->SetLineColor(kViolet);
        plot_cl->SetLineWidth(2);

        const double denominator = plot_b->Integral() + plot_cl->Integral();

        plot_b->Scale(dataIntegral / denominator);
        plot_cl->Scale(dataIntegral / denominator);

        THStack *stack = new THStack();
        heap.push_back(stack);

        stack->Add(plotData);
        stack->Add(plot_b);
        stack->Add(plot_cl);

        stack->Draw("nostack hist");
        legend->Draw();

        cout << title.str() << endl;
        printStats("data   ", plotData);
        printStats("b-jets ", plot_b);
        printStats("cl-jets", plot_cl);
        cout << endl;
    }
}

void printStats(const string &name, TH1 *plot)
{
    const int binOfTheCut = plot->GetXaxis()->FindBin(::ptrelCut);
    double subIntegral = plot->Integral(1, binOfTheCut);
    cout << name << ": " << subIntegral
        << " (" << static_cast<int>(subIntegral * 100./ plot->Integral()) << "%)" << endl;
}

TH1 *getProjection(TH2 *hist, const int &bin)
{
    TH1 *plot = dynamic_cast<TH1 *>(hist->ProjectionX("proj", bin, bin)->Clone());
    plot->SetDirectory(0);

    return plot;
}

std::string style(TFile *file, TH1 *hist)
{
    std::string filename = file->GetName();

    if ("data.root" == filename)
    {
        hist->SetLineWidth(2);

        return "data";
    }

    if ("ppmux.root" == filename)
    {
        hist->SetLineWidth(4);
        hist->SetLineColor(kAzure + 1);

        hist->Scale(48440000 * 1.76 / 886239 * luminosity);

        return "ppmux";
    }

    if ("ttbar.root" == filename)
    {
        hist->SetLineWidth(4);
        hist->SetLineColor(kRed);

        hist->Scale(157.5 / 996859 * luminosity);

        return "ttbar";
    }

    if ("pt15.root" == filename)
    {
        hist->SetLineColor(kGray+1);
        hist->SetLineWidth(2);

        hist->Scale(874100000 * 0.0039 / 2986966 * luminosity);

        return "pt15";
    }

    if ("pt30.root" == filename)
    {
        hist->SetLineColor(kGreen + 1);
        hist->SetLineWidth(2);

        hist->Scale(61160000 * 0.0126 / 7628669 * luminosity);

        return "pt30";
    }

    if ("pt50.root" == filename)
    {
        hist->SetLineColor(kViolet);
        hist->SetLineWidth(2);

        hist->Scale(7290000 * 0.0246 / 4262850 * luminosity);

        return "pt150";
    }

    if ("pt150.root" == filename)
    {
        hist->SetLineColor(kBlue);
        hist->SetLineWidth(2);

        hist->Scale(48100 * 0.056 / 489271 * luminosity);

        return "pt150";
    }

    cerr << "Didn't understand filename: " << filename
        << " Input is not used" << endl;

    return "";
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


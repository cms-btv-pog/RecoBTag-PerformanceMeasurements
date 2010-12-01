/**
 * S8GraphGroup
 * 
 *
 * Created by Samvel Khalatian on Nov 23, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>
#include <iostream>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>

#include "S8GraphGroup.h"

using std::cerr;
using std::cout;
using std::endl;

void fill(TGraph *graph, const int &point,
                      const Measurement &x, const Measurement &y)
{
    graph->SetPoint(point, x.first, y.first);
}

void fill(TGraphErrors *graph, const int &point,
                      const Measurement &x, const Measurement &y)
{
    fill(dynamic_cast<TGraph *>(graph), point, x, y);

    graph->SetPointError(point, sqrt(x.second), sqrt(y.second));
}


FlavouredEffGraphGroup::FlavouredEffGraphGroup(
    const BinnedSolution &binnedSolution)
{
    init(binnedSolution.size());
}

FlavouredEffGraphGroup::FlavouredEffGraphGroup(
    const BinnedNumericInputGroup &binnedInput)
{
    init(binnedInput.size());
}

void FlavouredEffGraphGroup::init(const int &size)
{
    b.reset(new TGraphErrors(size));
    b->SetMarkerStyle(8);
    b->SetMarkerColor(1);
    b->SetLineColor(1); 
    b->SetMarkerSize(1.5); 

    cl.reset(new TGraphErrors(size));
    cl->SetMarkerStyle(8);
    cl->SetMarkerColor(1);
    cl->SetLineColor(1); 
    cl->SetMarkerSize(1.5); 
}



EffGraphGroup::EffGraphGroup(const std::string &prefix,
                             const BinnedSolution &binnedSolution):
    mu(binnedSolution),
    tag(binnedSolution),
    _prefix(prefix)
{
    int point = 0;
    for(BinnedSolution::const_iterator solutionInBin = binnedSolution.begin();
        binnedSolution.end() != solutionInBin;
        ++solutionInBin, ++point)
    {
        fill(mu.b.get(), point,
             solutionInBin->bin,
             solutionInBin->solution.at("eff_mu_b"));

        fill(mu.cl.get(), point,
             solutionInBin->bin,
             solutionInBin->solution.at("eff_mu_cl"));

        fill(tag.b.get(), point,
             solutionInBin->bin,
             solutionInBin->solution.at("eff_tag_b"));

        fill(tag.cl.get(), point,
             solutionInBin->bin,
             solutionInBin->solution.at("eff_tag_cl"));
    }
}

EffGraphGroup::EffGraphGroup(const std::string &prefix,
                             const BinnedNumericInputGroup &binnedInput):
    mu(binnedInput),
    tag(binnedInput),
    _prefix(prefix)
{
    int point = 0;
    for(BinnedNumericInputGroup::const_iterator inputGroup =
            binnedInput.begin();
        binnedInput.end() != inputGroup;
        ++inputGroup, ++point)
    {
        fill(mu.b.get(), point,
             inputGroup->bin, inputGroup->efficiency.mu.b);

        fill(mu.cl.get(), point,
             inputGroup->bin, inputGroup->efficiency.mu.cl);

        fill(tag.b.get(), point,
             inputGroup->bin, inputGroup->efficiency.tag.b);

        fill(tag.cl.get(), point,
             inputGroup->bin, inputGroup->efficiency.tag.cl);
    }
}

EffGraphGroup::~EffGraphGroup()
{
    while(!_heaps.empty())
    {
        delete _heaps.top();

        _heaps.pop();
    }

    cout << "Eff Graph group is destroyed" << endl;
}

void EffGraphGroup::draw()
{
    TCanvas *canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle((_prefix + "Eff Mu B").c_str());
    canvas->SetGrid();
    mu.b->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle((_prefix + "Eff Mu CL").c_str());
    canvas->SetGrid();
    mu.cl->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle((_prefix + "Eff Tag B").c_str());
    canvas->SetGrid();
    tag.b->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle((_prefix + "Eff Tag CL").c_str());
    canvas->SetGrid();
    tag.cl->Draw("ap");
}

void EffGraphGroup::save(TDirectory *folder)
{
    TDirectory *subdir = folder->mkdir((_prefix + "efficiency").c_str());
    if (!subdir)
    {
        cerr << "failed to create efficiency subdir in ROOT file" << endl;

        return;
    }

    subdir->cd();

    mu.b->Write("eff_mu_b");
    mu.cl->Write("eff_mu_cl");

    tag.b->Write("eff_tag_b");
    tag.cl->Write("eff_tag_cl");

    folder->cd();
}



InputGraph::InputGraph(const int &size)
{
    all.reset(new TGraphErrors(size));
    all->SetMarkerStyle(8);
    all->SetMarkerColor(1);
    all->SetLineColor(1); 
    all->SetMarkerSize(1.2); 

    mu.reset(new TGraphErrors(size));
    mu->SetMarkerStyle(8);
    mu->SetMarkerColor(kRed);
    mu->SetLineColor(kRed); 
    mu->SetMarkerSize(1.2); 

    tag.reset(new TGraphErrors(size));
    tag->SetMarkerStyle(8);
    tag->SetMarkerColor(kGreen);
    tag->SetLineColor(kGreen); 
    tag->SetMarkerSize(1.2); 

    muTag.reset(new TGraphErrors(size));
    muTag->SetMarkerStyle(8);
    muTag->SetMarkerColor(kBlue);
    muTag->SetLineColor(kBlue); 
    muTag->SetMarkerSize(1.2); 
}



InputGraphGroup::InputGraphGroup(const BinnedNumericInputGroup &binnedInput):
    n(binnedInput.size()),
    p(binnedInput.size())
{
    int point = 0;
    for(BinnedNumericInputGroup::const_iterator inputGroup =
            binnedInput.begin();
        binnedInput.end() != inputGroup;
        ++inputGroup, ++point)
    {
        fill(n.all.get(), point,
             inputGroup->bin, inputGroup->input.n.all);

        fill(p.all.get(), point,
             inputGroup->bin, inputGroup->input.p.all);

        fill(n.mu.get(), point,
             inputGroup->bin, inputGroup->input.n.mu);

        fill(p.mu.get(), point,
             inputGroup->bin, inputGroup->input.p.mu);

        fill(n.tag.get(), point,
             inputGroup->bin, inputGroup->input.n.tag);

        fill(p.tag.get(), point,
             inputGroup->bin, inputGroup->input.p.tag);

        fill(n.muTag.get(), point,
             inputGroup->bin, inputGroup->input.n.muTag);

        fill(p.muTag.get(), point,
             inputGroup->bin, inputGroup->input.p.muTag);
    }
}

InputGraphGroup::~InputGraphGroup()
{
    while(!_heaps.empty())
    {
        delete _heaps.top();

        _heaps.pop();
    }

    cout << "Input Graph Group is destroyed" << endl;
}

void InputGraphGroup::draw()
{
    TCanvas *canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("(n) Inputs");
    canvas->SetGrid();

    TMultiGraph *graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) n.all->Clone(), "lp");
    graph->Add((TGraphErrors *) n.mu->Clone(), "lp");
    graph->Add((TGraphErrors *) n.tag->Clone(), "lp");
    graph->Add((TGraphErrors *) n.muTag->Clone(), "lp");
    graph->Draw("a");

    TLegend *legend = new TLegend(0.57,0.22,0.87,0.38, "Input");
    _heaps.push(legend);
    legend->SetMargin(0.12);
    legend->SetTextSize(0.027);
    legend->SetFillColor(10);
    legend->AddEntry(n.all.get(), "n", "p");
    legend->AddEntry(n.mu.get(), "n mu", "p");
    legend->AddEntry(n.tag.get(), "n tag", "p");
    legend->AddEntry(n.muTag.get(), "n muTag", "p");
    legend->Draw();

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("(p) Inputs");
    canvas->SetGrid();

    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) p.all->Clone(), "lp");
    graph->Add((TGraphErrors *) p.mu->Clone(), "lp");
    graph->Add((TGraphErrors *) p.tag->Clone(), "lp");
    graph->Add((TGraphErrors *) p.muTag->Clone(), "lp");
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Input");
    _heaps.push(legend);
    legend->SetMargin(0.12);
    legend->SetTextSize(0.027);
    legend->SetFillColor(10);
    legend->AddEntry(p.all.get(), "p", "p");
    legend->AddEntry(p.mu.get(), "p mu", "p");
    legend->AddEntry(p.tag.get(), "p tag", "p");
    legend->AddEntry(p.muTag.get(), "p muTag", "p");
    legend->Draw();
}

void InputGraphGroup::save(TDirectory *folder)
{
    TDirectory *subdir = folder->mkdir("input");
    if (!subdir)
    {
        cerr << "failed to create input subdir in ROOT file" << endl;

        return;
    }

    subdir->cd();

    n.all->Write("n");
    p.all->Write("p");

    n.mu->Write("n_mu");
    p.mu->Write("p_mu");

    n.tag->Write("n_tag");
    p.tag->Write("p_tag");

    n.muTag->Write("n_muTag");
    p.muTag->Write("p_muTag");

    folder->cd();
}


GraphGroup::GraphGroup(const BinnedNumericInputGroup &binnedInput,
                       const BinnedSolution &binnedSolution):
    mcEfficiency("mc", binnedInput),
    s8Efficiency("s8", binnedSolution),
    input(binnedInput)
{
    alpha.reset(new TGraphErrors(binnedInput.size()));
    alpha->SetMarkerStyle(23);
    alpha->SetMarkerColor(4);
    alpha->SetLineColor(4);
    alpha->SetMinimum(0.7);
    alpha->SetMaximum(1.3);

    beta.reset(new TGraphErrors(binnedInput.size()));
    beta->SetMarkerStyle(22); 
    beta->SetMarkerColor(3); 
    beta->SetLineColor(3); 
    beta->SetMinimum(0.7);
    beta->SetMaximum(1.3);

    gamma.reset(new TGraphErrors(binnedInput.size()));
    gamma->SetMarkerStyle(23);
    gamma->SetMarkerColor(2);
    gamma->SetLineColor(2);
    gamma->SetMinimum(0.7);
    gamma->SetMaximum(1.3);

    delta.reset(new TGraphErrors(binnedInput.size()));
    delta->SetMarkerStyle(22);
    delta->SetMarkerColor(1);
    delta->SetLineColor(1);
    delta->SetMinimum(0.7);
    delta->SetMaximum(1.3);

    kappaB.reset(new TGraphErrors(binnedInput.size()));
    kappaB->SetMarkerStyle(20); 
    kappaB->SetMarkerColor(1); 
    kappaB->SetLineColor(1); 
    kappaB->SetMinimum(0.7);
    kappaB->SetMaximum(1.3);

    kappaCL.reset(new TGraphErrors(binnedInput.size()));
    kappaCL->SetMarkerStyle(21); 
    kappaCL->SetMarkerColor(2); 
    kappaCL->SetLineColor(2); 
    kappaCL->SetMinimum(0.7);
    kappaCL->SetMaximum(1.3);

    int point = 0;
    for(BinnedNumericInputGroup::const_iterator inputGroup =
            binnedInput.begin();
        binnedInput.end() != inputGroup;
        ++inputGroup, ++point)
    {
        fill(alpha.get(), point,
             inputGroup->bin, inputGroup->coefficients.alpha);

        fill(beta.get(), point,
             inputGroup->bin, inputGroup->coefficients.beta);
        
        fill(gamma.get(), point,
             inputGroup->bin, inputGroup->coefficients.gamma);

        fill(delta.get(), point,
             inputGroup->bin, inputGroup->coefficients.delta);

        fill(kappaB.get(), point,
             inputGroup->bin, inputGroup->coefficients.kappaB);

        fill(kappaCL.get(), point,
             inputGroup->bin, inputGroup->coefficients.kappaCL);
    }
}

GraphGroup::~GraphGroup()
{
    while(!_heaps.empty())
    {
        delete _heaps.top();

        _heaps.pop();
    }

    cout << "Graph group is destroyed" << endl;
}

void GraphGroup::save(TDirectory *folder)
{
    TDirectory *subdir = folder->mkdir("coefficient");
    if (!subdir)
    {
        cerr << "failed to create subdir in ROOT file" << endl;

        return;
    }

    subdir->cd();

    alpha->Write("alpha");
    beta->Write("beta");
    gamma->Write("gamma");
    delta->Write("delta");
    kappaB->Write("kappaB");
    kappaCL->Write("kappaCL");

    folder->cd();

    mcEfficiency.save(folder);
    s8Efficiency.save(folder);
    input.save(folder);
}

void GraphGroup::draw()
{
    TCanvas *canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("Alpha");
    canvas->SetGrid();
    alpha->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("Beta");
    canvas->SetGrid();
    beta->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("Gamma");
    canvas->SetGrid();
    gamma->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("Delta");
    canvas->SetGrid();
    delta->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("KappaB");
    canvas->SetGrid();
    kappaB->Draw("ap");

    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetTitle("KappaCL");
    canvas->SetGrid();
    kappaCL->Draw("ap");

    mcEfficiency.draw();
    s8Efficiency.draw();

    input.draw();
}

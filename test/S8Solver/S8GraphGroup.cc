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

void fill(TGraphErrors *scale, const TGraphErrors *s8, const TGraphErrors *mc)
{
    using std::make_pair;

    for(int point = 0; s8->GetN() > point; ++point)
    {
        double s8X;
        double s8Y;
        s8->GetPoint(point, s8X, s8Y);

        double mcX;
        double mcY;
        mc->GetPoint(point, mcX, mcY);

        Measurement x = make_pair(s8X, pow(s8->GetErrorX(point), 2));
        Measurement y = make_pair(s8Y, pow(s8->GetErrorY(point), 2)) /
            make_pair(mcY, pow(mc->GetErrorY(point), 2));

        fill(scale, point, x, y);
    }
}

void fill(EffGraphGroup &scale,
          const EffGraphGroup &s8,
          const EffGraphGroup &mc)
{
    fill(scale.mu.b.get(), s8.mu.b.get(), mc.mu.b.get());
    fill(scale.mu.cl.get(), s8.mu.cl.get(), mc.mu.cl.get());
    fill(scale.tag.b.get(), s8.tag.b.get(), mc.tag.b.get());
    fill(scale.tag.cl.get(), s8.tag.cl.get(), mc.tag.cl.get());
}


FlavouredEffGraphGroup::FlavouredEffGraphGroup(const int &size)
{
    init(size, true);
}

FlavouredEffGraphGroup::FlavouredEffGraphGroup(
    const BinnedSolution &binnedSolution)
{
    init(binnedSolution.size());
}

FlavouredEffGraphGroup::FlavouredEffGraphGroup(
    const BinnedNumericInputGroup &binnedInput)
{
    init(binnedInput.size(), true);
}

void FlavouredEffGraphGroup::init(const int &size, const bool &isMC)
{
    b.reset(new TGraphErrors(size));
    b->SetMarkerStyle(8);
    b->SetMarkerColor(isMC ? 1 : kGreen + 1);
    b->SetLineColor(1); 
    b->SetMarkerSize(1.5); 

    cl.reset(new TGraphErrors(size));
    cl->SetMarkerStyle(8);
    cl->SetMarkerColor(isMC ? 1 : kRed + 1);
    cl->SetLineColor(isMC ? 1 : kRed + 1); 
    cl->SetMarkerSize(1.5); 
}



EffGraphGroup::EffGraphGroup(const int &size):
    mu(size),
    tag(size)
{
}

EffGraphGroup::EffGraphGroup(const BinnedSolution &binnedSolution):
    mu(binnedSolution),
    tag(binnedSolution)
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

EffGraphGroup::EffGraphGroup(const BinnedNumericInputGroup &binnedInput):
    mu(binnedInput),
    tag(binnedInput)
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

void EffGraphGroup::save(TDirectory *)
{
    mu.b->Write("eff_mu_b");
    mu.cl->Write("eff_mu_cl");

    tag.b->Write("eff_tag_b");
    tag.cl->Write("eff_tag_cl");
}



EffGraph::EffGraph(const BinnedNumericInputGroup &binnedInput,
                   const BinnedSolution &binnedSolution):
    mc(binnedInput),
    s8(binnedSolution),
    scale(binnedSolution.size())
{
    fill(scale, s8, mc);
}

EffGraph::~EffGraph()
{
    while(!_heaps.empty())
    {
        delete _heaps.top();

        _heaps.pop();
    }

    cout << "Efficiency Group is destroyed" << endl;
}

void EffGraph::draw()
{
    TCanvas *canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetWindowSize(1024, 768);
    canvas->Divide(2, 2);
    canvas->SetTitle("Eff");
    canvas->SetGrid();

    canvas->cd(1)->SetGrid();
    TMultiGraph *graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.mu.b->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.mu.b->Clone(), "lp");
    graph->SetMinimum(0.4);
    graph->SetMaximum(1.0);
    graph->Draw("a");

    TLegend *legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Mu b");
    _heaps.push(legend);
    legend->AddEntry(mc.mu.b.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.mu.b.get(), "System8", "p");
    legend->Draw();

    canvas->cd(2)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.mu.cl->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.mu.cl->Clone(), "lp");
    graph->SetMinimum(0.1);
    graph->SetMaximum(0.5);
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Mu cl");
    _heaps.push(legend);
    legend->AddEntry(mc.mu.cl.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.mu.cl.get(), "System8", "p");
    legend->Draw();

    canvas->cd(3)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.tag.b->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.tag.b->Clone(), "lp");
    graph->SetMinimum(0.4);
    graph->SetMaximum(1.0);
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Tag b");
    _heaps.push(legend);
    legend->AddEntry(mc.tag.b.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.tag.b.get(), "System8", "p");
    legend->Draw();

    canvas->cd(4)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.tag.cl->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.tag.cl->Clone(), "lp");
    graph->SetMinimum(0.1);
    graph->SetMaximum(0.5);
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Tag cl");
    _heaps.push(legend);
    legend->AddEntry(mc.tag.cl.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.tag.cl.get(), "System8", "p");
    legend->Draw();

    // Scales
    //
    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetWindowSize(1024, 768);
    canvas->Divide(2, 2);
    canvas->SetTitle("Eff");
    canvas->SetGrid();

    canvas->cd(1)->SetGrid();
    scale.mu.b->Draw("ap");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Mu b");
    _heaps.push(legend);
    legend->AddEntry(scale.mu.b.get(), "Scale Factor", "p");
    legend->Draw();

    canvas->cd(2)->SetGrid();
    scale.mu.cl->Draw("ap");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Mu cl");
    _heaps.push(legend);
    legend->AddEntry(scale.mu.cl.get(), "Scale Factor", "p");
    legend->Draw();

    canvas->cd(3)->SetGrid();
    scale.tag.b->Draw("ap");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Tag b");
    _heaps.push(legend);
    legend->AddEntry(scale.tag.b.get(), "Scale Factor", "p");
    legend->Draw();

    canvas->cd(4)->SetGrid();
    scale.tag.cl->Draw("ap");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Eff Tag cl");
    _heaps.push(legend);
    legend->AddEntry(scale.tag.cl.get(), "Scale Factor", "p");
    legend->Draw();
}

void EffGraph::save(TDirectory *folder)
{
    TDirectory *subdir = folder->mkdir("mcefficiency");
    if (!subdir)
    {
        cerr << "failed to create mcefficiency subdir in ROOT file" << endl;

        return;
    }

    subdir->cd();
    mc.save(subdir);
    folder->cd();

    subdir = folder->mkdir("s8efficiency");
    if (!subdir)
    {
        cerr << "failed to create s8efficiency subdir in ROOT file" << endl;

        return;
    }

    subdir->cd();
    s8.save(subdir);
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
    canvas->SetWindowSize(1024, 640);
    canvas->Divide(2, 1);
    canvas->SetTitle("Inputs");

    // (n)
    //
    canvas->cd(1)->SetGrid();
    TMultiGraph *graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) n.all->Clone(), "lp");
    graph->Add((TGraphErrors *) n.mu->Clone(), "lp");
    graph->Add((TGraphErrors *) n.tag->Clone(), "lp");
    graph->Add((TGraphErrors *) n.muTag->Clone(), "lp");
    graph->Draw("a");

    TLegend *legend = new TLegend(0.57,0.22,0.87,0.38, "Input");
    _heaps.push(legend);
    legend->SetTextSize(0.03);
    legend->AddEntry(n.all.get(), "n", "p");
    legend->AddEntry(n.mu.get(), "n mu", "p");
    legend->AddEntry(n.tag.get(), "n tag", "p");
    legend->AddEntry(n.muTag.get(), "n muTag", "p");
    legend->Draw();

    // (p)
    //
    canvas->cd(2)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) p.all->Clone(), "lp");
    graph->Add((TGraphErrors *) p.mu->Clone(), "lp");
    graph->Add((TGraphErrors *) p.tag->Clone(), "lp");
    graph->Add((TGraphErrors *) p.muTag->Clone(), "lp");
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Input");
    _heaps.push(legend);
    legend->SetTextSize(0.03);
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
    efficiency(binnedInput, binnedSolution),
    input(binnedInput)
{
    alpha.reset(new TGraphErrors(binnedInput.size()));
    alpha->SetMarkerStyle(23);
    alpha->SetMarkerColor(kRed + 1);
    alpha->SetLineColor(kRed + 1);
    alpha->SetMinimum(0.7);
    alpha->SetMaximum(1.3);

    beta.reset(new TGraphErrors(binnedInput.size()));
    beta->SetMarkerStyle(22); 
    beta->SetMarkerColor(kGreen + 1); 
    beta->SetLineColor(kGreen + 1); 
    beta->SetMinimum(0.7);
    beta->SetMaximum(1.3);

    gamma.reset(new TGraphErrors(binnedInput.size()));
    gamma->SetMarkerStyle(23);
    gamma->SetMarkerColor(kRed + 1);
    gamma->SetLineColor(kRed + 1);
    gamma->SetMinimum(0.7);
    gamma->SetMaximum(1.3);

    delta.reset(new TGraphErrors(binnedInput.size()));
    delta->SetMarkerStyle(22);
    delta->SetMarkerColor(kGreen + 1);
    delta->SetLineColor(kGreen + 1);
    delta->SetMinimum(0.7);
    delta->SetMaximum(1.3);

    kappaB.reset(new TGraphErrors(binnedInput.size()));
    kappaB->SetMarkerStyle(20); 
    kappaB->SetMarkerColor(kGreen + 1); 
    kappaB->SetLineColor(kGreen + 1); 
    kappaB->SetMinimum(0.7);
    kappaB->SetMaximum(1.3);

    kappaCL.reset(new TGraphErrors(binnedInput.size()));
    kappaCL->SetMarkerStyle(21); 
    kappaCL->SetMarkerColor(kRed + 1); 
    kappaCL->SetLineColor(kRed + 1); 
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

    efficiency.save(folder);
    input.save(folder);
}

void GraphGroup::draw()
{
    TCanvas *canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetWindowSize(1200,640);
    canvas->Divide(3, 1);
    canvas->SetTitle("Coefficients");

    // mu
    //
    canvas->cd(1)->SetGrid();
    TMultiGraph *graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) gamma->Clone(), "lp");
    graph->Add((TGraphErrors *) delta->Clone(), "lp");
    graph->SetMinimum(0.8);
    graph->SetMaximum(1.2);
    graph->Draw("a");

    TLegend *legend = new TLegend(0.57,0.22,0.87,0.38, "Coeffcients");
    _heaps.push(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(gamma.get(), "gamma", "p");
    legend->AddEntry(delta.get(), "delta", "p");
    legend->Draw();

    // tag
    //
    canvas->cd(2)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) alpha->Clone(), "lp");
    graph->Add((TGraphErrors *) beta->Clone(), "lp");
    graph->SetMinimum(0.8);
    graph->SetMaximum(1.2);
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Coeffcients");
    _heaps.push(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(alpha.get(), "alpha", "p");
    legend->AddEntry(beta.get(), "beta", "p");
    legend->Draw();

    // muTag
    //
    canvas->cd(3)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) kappaCL->Clone(), "lp");
    graph->Add((TGraphErrors *) kappaB->Clone(), "lp");
    graph->SetMinimum(0.8);
    graph->SetMaximum(1.2);
    graph->Draw("a");

    legend = new TLegend(0.57,0.22,0.87,0.38, "Coeffcients");
    _heaps.push(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(kappaCL.get(), "kappaCL", "p");
    legend->AddEntry(kappaB.get(), "kappaB", "p");
    legend->Draw();

    efficiency.draw();

    input.draw();
}

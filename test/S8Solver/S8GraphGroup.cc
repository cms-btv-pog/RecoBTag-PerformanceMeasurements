/**
 * S8GraphGroup
 * 
 *
 * Created by Samvel Khalatian on Nov 23, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <TAxis.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMultiGraph.h>

#include "S8GraphGroup.h"

using std::cerr;
using std::cout;
using std::endl;

TLatex *createLabel(const double &luminosity = 31.6)
{
    std::ostringstream title;
    title << "#splitline{CMS Preliminary 2010}{"
        << std::fixed << std::setprecision(2) << luminosity
        << " pb^{-1} at #sqrt{s} = 7 TeV}";

    TLatex *label = new TLatex(3.570061, 23.08044, title.str().c_str());
    label->SetNDC();
    label->SetTextAlign(13);
    label->SetX(0.4);
    label->SetY(0.908);

    return label;
}

void style(TLegend *legend)
{ 
    legend->SetBorderSize(1);
    legend->SetLineStyle(1);
    legend->SetTextFont(43);
    legend->SetTextSizePixels(24);
    legend->SetFillColor(0);
}

std::string getXTitle(const Graph::Type &type)
{
    switch(type)
    {
        case Graph::PT:  return "jet p_{T} [GeV/c]";
        case Graph::ETA: return "jet #eta";
        case Graph::PHI: return "jet #phi";
    }

    return "";
}

void setXtitle(const Graph::Type &type, TGraphErrors *graph)
{
    graph->GetXaxis()->SetTitle(getXTitle(type).c_str());
}

void setXtitle(const Graph::Type &type, TMultiGraph *graph)
{
    graph->GetXaxis()->SetTitle(getXTitle(type).c_str());
}

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
    b->SetLineColor(isMC ? 1 : kGreen + 1); 
    b->SetMarkerSize(1.5); 

    cl.reset(new TGraphErrors(size));
    cl->SetMarkerStyle(8);
    cl->SetMarkerColor(isMC ? 1 : kRed + 1);
    cl->SetLineColor(isMC ? 1 : kRed + 1); 
    cl->SetMarkerSize(1.5); 
}



EffGraphGroup::EffGraphGroup(const Graph::Type &type, const int &size):
    mu(size),
    tag(size)
{
}

EffGraphGroup::EffGraphGroup(const Graph::Type &type,
                             const BinnedSolution &binnedSolution):
    mu(binnedSolution),
    tag(binnedSolution)
{
    mu.b->GetYaxis()->SetTitle("Eff Mu b");
    mu.cl->GetYaxis()->SetTitle("Eff Mu cl");

    tag.b->GetYaxis()->SetTitle("Eff Mu b");
    tag.cl->GetYaxis()->SetTitle("Eff Mu cl");

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

EffGraphGroup::EffGraphGroup(const Graph::Type &type,
                             const BinnedNumericInputGroup &binnedInput):
    mu(binnedInput),
    tag(binnedInput)
{
    mu.b->GetYaxis()->SetTitle("Eff Mu b");
    mu.cl->GetYaxis()->SetTitle("Eff Mu cl");

    setXtitle(type, mu.b.get());
    setXtitle(type, mu.cl.get());

    tag.b->GetYaxis()->SetTitle("Eff Mu b");
    tag.cl->GetYaxis()->SetTitle("Eff Mu cl");

    setXtitle(type, tag.b.get());
    setXtitle(type, tag.cl.get());

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



EffGraph::EffGraph(const Graph::Type &type,
                   const BinnedNumericInputGroup &binnedInput,
                   const BinnedSolution &binnedSolution):
    mc(type, binnedInput),
    s8(type, binnedSolution),
    scale(type, binnedSolution.size()),
    _type(type)
{
    scale.mu.b->GetYaxis()->SetTitle("Scale Factor");
    scale.mu.cl->GetYaxis()->SetTitle("Scale Factor");

    scale.tag.b->GetYaxis()->SetTitle("Scale Factor");
    scale.tag.cl->GetYaxis()->SetTitle("Scale Factor");

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
    graph->SetMaximum(1.1);
    graph->Draw("a");
    graph->GetYaxis()->SetTitle("#epsilon_{b}^{#mu}");
    setXtitle(_type, graph);

    TLegend *legend = new TLegend(0.57,0.22,0.91,0.38);
    _heaps.push(legend);
    legend->AddEntry(mc.mu.b.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.mu.b.get(), "Data", "p");
    style(legend);
    legend->Draw();

    TLatex *label = createLabel();
    _heaps.push(label);
    label->Draw();

    canvas->cd(2)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.mu.cl->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.mu.cl->Clone(), "lp");
    graph->SetMinimum(0.1);
    graph->SetMaximum(0.6);
    graph->Draw("a");
    graph->GetYaxis()->SetTitle("#epsilon_{cl}^{#mu}");
    setXtitle(_type, graph);

    legend = new TLegend(0.57,0.22,0.91,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(mc.mu.cl.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.mu.cl.get(), "Data", "p");
    legend->Draw();

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    canvas->cd(3)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.tag.b->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.tag.b->Clone(), "lp");
    graph->SetMinimum(0.4);
    graph->SetMaximum(1.1);
    graph->Draw("a");
    graph->GetYaxis()->SetTitle("#epsilon_{b}^{tag}");
    setXtitle(_type, graph);

    legend = new TLegend(0.57,0.22,0.91,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(mc.tag.b.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.tag.b.get(), "Data", "p");
    legend->Draw();

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    canvas->cd(4)->SetGrid();
    graph = new TMultiGraph();
    _heaps.push(graph);
    graph->Add((TGraphErrors *) mc.tag.cl->Clone(), "lp");
    graph->Add((TGraphErrors *) s8.tag.cl->Clone(), "lp");
    graph->SetMinimum(0.1);
    graph->SetMaximum(0.6);
    graph->Draw("a");
    graph->GetYaxis()->SetTitle("#epsilon_{cl}^{tag}");
    setXtitle(_type, graph);

    legend = new TLegend(0.57,0.22,0.91,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(mc.tag.cl.get(), "Monte-Carlo", "p");
    legend->AddEntry(s8.tag.cl.get(), "Data", "p");
    legend->Draw();

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    // Scales
    //
    canvas = new TCanvas();
    _heaps.push(canvas);
    canvas->SetWindowSize(1024, 768);
    canvas->Divide(2, 2);
    canvas->SetTitle("Eff");
    canvas->SetGrid();

    canvas->cd(1)->SetGrid();
    scale.mu.b->SetMinimum(.6);
    scale.mu.b->SetMaximum(1.4);
    scale.mu.b->Draw("ap");
    scale.mu.b->GetYaxis()->SetTitle("SF_{b}^{#mu}");
    setXtitle(_type, scale.mu.b.get());

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    canvas->cd(2)->SetGrid();
    scale.mu.cl->SetMinimum(.6);
    scale.mu.cl->SetMaximum(1.4);
    scale.mu.cl->Draw("ap");
    scale.mu.cl->GetYaxis()->SetTitle("SF_{cl}^{#mu}");
    setXtitle(_type, scale.mu.cl.get());

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    canvas->cd(3)->SetGrid();
    scale.tag.b->SetMinimum(.6);
    scale.tag.b->SetMaximum(1.4);
    scale.tag.b->Draw("ap");
    scale.tag.b->GetYaxis()->SetTitle("SF_{b}^{tag}");
    setXtitle(_type, scale.tag.b.get());

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    canvas->cd(4)->SetGrid();
    scale.tag.cl->SetMinimum(.6);
    scale.tag.cl->SetMaximum(1.4);
    scale.tag.cl->Draw("ap");
    scale.tag.cl->GetYaxis()->SetTitle("SF_{cl}^{tag}");
    setXtitle(_type, scale.tag.cl.get());

    label = createLabel();
    _heaps.push(label);
    label->Draw();
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



InputGraph::InputGraph(const Graph::Type &type, const int &size)
{
    all.reset(new TGraphErrors(size));
    all->SetMarkerStyle(8);
    all->SetMarkerColor(1);
    all->SetLineColor(1); 
    all->SetMarkerSize(1.2); 
    all->GetYaxis()->SetTitle("entries");
    setXtitle(type, all.get());

    mu.reset(new TGraphErrors(size));
    mu->SetMarkerStyle(8);
    mu->SetMarkerColor(kRed);
    mu->SetLineColor(kRed); 
    mu->SetMarkerSize(1.2); 
    mu->GetYaxis()->SetTitle("entries");
    setXtitle(type, mu.get());

    tag.reset(new TGraphErrors(size));
    tag->SetMarkerStyle(8);
    tag->SetMarkerColor(kGreen);
    tag->SetLineColor(kGreen); 
    tag->SetMarkerSize(1.2); 
    tag->GetYaxis()->SetTitle("entries");
    setXtitle(type, tag.get());

    muTag.reset(new TGraphErrors(size));
    muTag->SetMarkerStyle(8);
    muTag->SetMarkerColor(kBlue);
    muTag->SetLineColor(kBlue); 
    muTag->SetMarkerSize(1.2); 
    muTag->GetYaxis()->SetTitle("entries");
    setXtitle(type, muTag.get());
}



InputGraphGroup::InputGraphGroup(const Graph::Type &type,
                                 const BinnedNumericInputGroup &binnedInput):
    n(type, binnedInput.size()),
    p(type, binnedInput.size()),
    _type(type)
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
    graph->GetYaxis()->SetTitle("entries");
    setXtitle(_type, graph);

    TLegend *legend = new TLegend(0.57,0.22,0.87,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(n.all.get(), "n", "p");
    legend->AddEntry(n.mu.get(), "n mu", "p");
    legend->AddEntry(n.tag.get(), "n tag", "p");
    legend->AddEntry(n.muTag.get(), "n muTag", "p");
    legend->Draw();

    TLatex *label = createLabel();
    _heaps.push(label);
    label->Draw();

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
    graph->GetYaxis()->SetTitle("entries");
    setXtitle(_type, graph);

    legend = new TLegend(0.57,0.22,0.87,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(p.all.get(), "p", "p");
    legend->AddEntry(p.mu.get(), "p mu", "p");
    legend->AddEntry(p.tag.get(), "p tag", "p");
    legend->AddEntry(p.muTag.get(), "p muTag", "p");
    legend->Draw();

    label = createLabel();
    _heaps.push(label);
    label->Draw();
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
                       const BinnedSolution &binnedSolution,
                       const Graph::Type &type):
    efficiency(type, binnedInput, binnedSolution),
    input(type, binnedInput),
    _type(type)
{
    alpha.reset(new TGraphErrors(binnedInput.size()));
    alpha->SetMarkerStyle(23);
    alpha->SetMarkerColor(kRed + 1);
    alpha->SetLineColor(kRed + 1);
    alpha->SetMinimum(0.7);
    alpha->SetMaximum(1.3);
    alpha->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, alpha.get());

    beta.reset(new TGraphErrors(binnedInput.size()));
    beta->SetMarkerStyle(22); 
    beta->SetMarkerColor(kGreen + 1); 
    beta->SetLineColor(kGreen + 1); 
    beta->SetMinimum(0.7);
    beta->SetMaximum(1.3);
    beta->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, beta.get());

    gamma.reset(new TGraphErrors(binnedInput.size()));
    gamma->SetMarkerStyle(23);
    gamma->SetMarkerColor(kRed + 1);
    gamma->SetLineColor(kRed + 1);
    gamma->SetMinimum(0.7);
    gamma->SetMaximum(1.3);
    gamma->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, gamma.get());

    delta.reset(new TGraphErrors(binnedInput.size()));
    delta->SetMarkerStyle(22);
    delta->SetMarkerColor(kGreen + 1);
    delta->SetLineColor(kGreen + 1);
    delta->SetMinimum(0.7);
    delta->SetMaximum(1.3);
    delta->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, delta.get());

    kappaB.reset(new TGraphErrors(binnedInput.size()));
    kappaB->SetMarkerStyle(20); 
    kappaB->SetMarkerColor(kGreen + 1); 
    kappaB->SetLineColor(kGreen + 1); 
    kappaB->SetMinimum(0.7);
    kappaB->SetMaximum(1.3);
    kappaB->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, kappaB.get());

    kappaCL.reset(new TGraphErrors(binnedInput.size()));
    kappaCL->SetMarkerStyle(21); 
    kappaCL->SetMarkerColor(kRed + 1); 
    kappaCL->SetLineColor(kRed + 1); 
    kappaCL->SetMinimum(0.7);
    kappaCL->SetMaximum(1.3);
    kappaCL->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, kappaCL.get());

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
    graph->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, graph);

    TLegend *legend = new TLegend(0.57,0.22,0.97,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(gamma.get(), "gamma (cl)", "p");
    legend->AddEntry(delta.get(), "delta (b)", "p");
    legend->Draw();

    TLatex *label = createLabel();
    _heaps.push(label);
    label->Draw();

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
    graph->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, graph);

    legend = new TLegend(0.57,0.22,0.97,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(alpha.get(), "alpha (cl)", "p");
    legend->AddEntry(beta.get(), "beta (b)", "p");
    legend->Draw();

    label = createLabel();
    _heaps.push(label);
    label->Draw();

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
    graph->GetYaxis()->SetTitle("a.u.");
    setXtitle(_type, graph);

    legend = new TLegend(0.57,0.22,0.97,0.38);
    _heaps.push(legend);
    style(legend);
    legend->AddEntry(kappaCL.get(), "kappaCL", "p");
    legend->AddEntry(kappaB.get(), "kappaB", "p");
    legend->Draw();

    label = createLabel();
    _heaps.push(label);
    label->Draw();

    efficiency.draw();

    input.draw();
}

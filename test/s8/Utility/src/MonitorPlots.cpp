/**
 * MonitorPlots
 * s8
 *
 * Created by Samvel Khalatian on Oct 15, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include <TDirectory.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TLorentzVector.h>

#include "S8Tree/interface/S8Jet.h"
#include "S8Tree/interface/S8Lepton.h"

#include "Utility/interface/MonitorPlots.h"

using std::runtime_error;

using s8::Monitor;
using s8::MonitorBase;
using s8::MonitorDelta;
using s8::MonitorLepton;
using s8::MonitorJet;
using s8::MonitorMuonInJet;

Monitor::Monitor()
{
}

Monitor::~Monitor() throw()
{
}



MonitorBase::MonitorBase(const std::string &prefix)
{
    _pt.reset(new TH1F((prefix + "_pt").c_str(), "Pt", 230, 0, 230));
    _pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    _pt->Sumw2();

    _eta.reset(new TH1F((prefix + "_eta").c_str(), "Eta", 60, -3, 3));
    _eta->GetXaxis()->SetTitle("#eta");
    _eta->Sumw2();

    _phi.reset(new TH1F((prefix + "_phi").c_str(), "Phi", 70, -3.5, 3.5));
    _phi->GetXaxis()->SetTitle("#phi");
    _phi->Sumw2();
}

MonitorBase::~MonitorBase() throw()
{
}

void MonitorBase::fill(const double &weight, const Plot &plot, const double &value)
{
    TH1 *histogram = 0;

    switch(plot)
    {
        case PT:  histogram = _pt.get();
                  break;

        case ETA: histogram = _eta.get();
                  break;

        case PHI: histogram = _phi.get();
                  break;

        default: throw runtime_error("Unsuppored Plot: " + plot);
    }

    histogram->Fill(value, weight);
}

void MonitorBase::save(TDirectory *)
{
    _pt->Write();
    _eta->Write();
    _phi->Write();
}

TH1 *MonitorBase::plot(const Plot &plot)
{
    TH1 *histogram;

    switch(plot)
    {
        case PT:  histogram = _pt.get();
                  break;

        case ETA: histogram = _eta.get();
                  break;

        case PHI: histogram = _phi.get();
                  break;

        default: throw runtime_error("Unsupported Plot");
    }

    return histogram;
}



MonitorDelta::MonitorDelta(const std::string &prefix)
{
    _ptrel.reset(new TH1F((prefix + "_ptrel").c_str(), "Ptrel", 50, 0, 5));
    _ptrel->GetXaxis()->SetTitle("p_{T}^{rel} [GeV/c]");
    _ptrel->Sumw2();

    _deltaR.reset(new TH1F((prefix + "_dr").c_str(), "Delta R", 60, 0, 0.55));
    _deltaR->GetXaxis()->SetTitle("#DeltaR");
    _deltaR->Sumw2();

    _deltaPhi.reset(new TH1F((prefix + "_dphi").c_str(), "Delta Phi", 120, -0.55, 0.55));
    _deltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    _deltaPhi->Sumw2();

    _deltaEta.reset(new TH1F((prefix + "_deta").c_str(), "Delta Eta", 120, -0.55, 0.55));
    _deltaEta->GetXaxis()->SetTitle("#Delta#eta");
    _deltaEta->Sumw2();
}

MonitorDelta::~MonitorDelta() throw()
{
}

MonitorJet *MonitorMuonInJet::jet()
{
    return &_jet;
}

MonitorDelta *MonitorMuonInJet::delta()
{
    return &_delta;
}

void MonitorDelta::fill(const double &weight, const TLorentzVector *p4_1, const TLorentzVector *p4_2)
{
    _ptrel->Fill(p4_1->Vect().Perp(p4_2->Vect()), weight);
    _deltaR->Fill(p4_1->DeltaR(*p4_2), weight);
    _deltaPhi->Fill(p4_1->DeltaPhi(*p4_2), weight);
    _deltaEta->Fill(p4_1->Eta() - p4_2->Eta(), weight);
}

void MonitorDelta::save(TDirectory *)
{
    _ptrel->Write();
    _deltaR->Write();
    _deltaPhi->Write();
    _deltaEta->Write();
}

TH1 *MonitorDelta::plot(const Plot &plot)
{
    TH1 *histogram;

    switch(plot)
    {
        case PTREL:     histogram = _ptrel.get();
                        break;

        case DELTA_R:   histogram = _deltaR.get();
                        break;

        case DELTA_PHI: histogram = _deltaPhi.get();
                        break;

        case DELTA_ETA: histogram = _deltaEta.get();
                        break;

        default: throw runtime_error("Unsupported Plot");
    }

    return histogram;
}



MonitorLepton::MonitorLepton(const std::string &prefix):
    MonitorBase(prefix)
{
}

void MonitorLepton::fill(const double &weight, const Lepton *lepton)
{
    MonitorBase::fill(weight, PT,  lepton->p4()->Pt());
    MonitorBase::fill(weight, ETA, lepton->p4()->Eta());
    MonitorBase::fill(weight, PHI, lepton->p4()->Phi());
}



MonitorJet::MonitorJet(const std::string &prefix):
    MonitorBase(prefix)
{
    _discriminator.reset(new TH1F((prefix + "_disc").c_str(), "Discriminator", 100, 0, 10));
    _discriminator->GetXaxis()->SetTitle("b-Tag discriminator");
    _discriminator->Sumw2();
}

MonitorJet::~MonitorJet() throw()
{
}

void MonitorJet::fill(const double &weight, const Jet *jet)
{
    MonitorBase::fill(weight, PT,  jet->p4()->Pt());
    MonitorBase::fill(weight, ETA, jet->p4()->Eta());
    MonitorBase::fill(weight, PHI, jet->p4()->Phi());
}

void MonitorJet::fillDiscriminator(const double &weight, const double &value)
{
    _discriminator->Fill(value, weight);
}

void MonitorJet::save(TDirectory *dir)
{
    MonitorBase::save(dir);

    _discriminator->Write();
}



MonitorMuonInJet::MonitorMuonInJet(const std::string &prefix):
    _muon(prefix + "mu"),
    _jet(prefix + "jet"),
    _delta(prefix + "mujet")
{
    TAxis *pt_axis = _jet.plot(MonitorBase::PT)->GetXaxis();
    TAxis *eta_axis = _jet.plot(MonitorBase::ETA)->GetXaxis();

    _pt_eta_pv.reset(new TH3F((prefix + "_pt_eta_pv").c_str(), "Pt Eta Npv",
                    pt_axis->GetNbins(), pt_axis->GetXmin(), pt_axis->GetXmax(),
                    eta_axis->GetNbins(), eta_axis->GetXmin(), eta_axis->GetXmax(),
                    10, 0, 10));
    _pt_eta_pv->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    _pt_eta_pv->GetYaxis()->SetTitle("#eta");
    _pt_eta_pv->GetZaxis()->SetTitle("N_{PV}");
    _pt_eta_pv->Sumw2();
}

MonitorMuonInJet::~MonitorMuonInJet() throw()
{
}

void MonitorMuonInJet::fill(const double &weight,
                            const Lepton *muon,
                            const Jet *jet,
                            const int &npv)
{
    _muon.fill(weight, muon);
    _jet.fill(weight, jet);
    _delta.fill(weight, muon->p4(), jet->p4());
    _pt_eta_pv->Fill(jet->p4()->Pt(), jet->p4()->Eta(), npv, weight);
}

void MonitorMuonInJet::save(TDirectory *dir)
{
    _muon.save(dir);
    _jet.save(dir);
    _delta.save(dir);
    _pt_eta_pv->Write();
}

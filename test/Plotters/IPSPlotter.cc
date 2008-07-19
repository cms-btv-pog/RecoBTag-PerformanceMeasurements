
#include <assert.h>
#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TTree.h"

#include "RecoBTag/PerformanceMeasurements/test/Plotters/IPSPlotter.h"

#include "SimTracker/TrackHistory/interface/TrackCategories.h"

ClassImp(IPSPlotter)

void IPSPlotter::Book()
{
    // Histograms
    ips_ = new TH1F("ips","Positive and negative tag IPS distribution", 60, -10., 10.);
    ipsB_ = new TH1F("ipsB","Positive and negative tag IPS distribution for jet with B tracks", 60, -10., 10.);
    ipsC_ = new TH1F("ipsC","Positive and negative tag IPS distribution for jet with C tracks", 60, -10., 10.);
    ipsKs_ = new TH1F("ipsKs","Positive and negative tag IPS distribution for jet with Ks tracks", 60, -10., 10.);
    ipsElse_ = new TH1F("ipsElse","Positive and negative tag IPS distribution for jet with ELSE tracks", 60, -10., 10.);
    ipsLight_ = new TH1F("ipslight","Positive and negative tag IPS distribution for jet with light tracks", 60, -10., 10.);
    ipsLambda_ = new TH1F("ipsLambda","Positive and negative tag IPS distribution for jet with Lambda tracks", 60, -10., 10.);
    ipsInteraction_ = new TH1F("ipsInteraction","Positive and negative tag IPS distribution for jet with Interaction tracks", 60, -10., 10.);
    ipsConversion_ = new TH1F("ipsConversion","Positive and negative tag IPS distribution for jet with conversion tracks", 60, -10., 10.);
}

void IPSPlotter::Fill(BTagEvent * event)
{
    std::cout << event->njets << " " << event->tracks.size() << std::endl;

    for (std::size_t i = 0; i<event->tracks.size(); ++i)
        for (std::size_t j = 0; j<event->tracks[i].ip3d.size(); ++j)
        {
            TrackCategories::Flags const & is = event->tracks[i].is[j];

            Double_t value = event->tracks[i].ip3d[j];
            Double_t error = event->tracks[i].ip3dSigma[j];
            
            Double_t ips3d = 0;
            
            if(error != 0.0) ips3d = value/error;

            ips_->Fill(ips3d);

            if ( is[TrackCategories::BWeakDecay] )
                ipsB_->Fill(ips3d);
            else if ( is[TrackCategories::KsDecay] )
                ipsKs_->Fill(ips3d);
            else if ( is[TrackCategories::Conversion] )
                ipsConversion_->Fill(ips3d);
            else if ( is[TrackCategories::LambdaDecay] )
                ipsLambda_->Fill(ips3d);
            else if ( is[TrackCategories::Interaction] )
                ipsInteraction_->Fill(ips3d);
            else if ( is[TrackCategories::CWeakDecay] )
                ipsC_->Fill(ips3d);
            else if (
                !is[TrackCategories::BWeakDecay] &&
                !is[TrackCategories::CWeakDecay]
            )
                ipsLight_->Fill(ips3d);
            else
                ipsElse_->Fill(ips3d);
        }
}

void IPSPlotter::Write()
{
    std::cout << "Saving Histograms" << std::endl;

    TFile file("IPSPlotter.root", "RECREATE");

    // Writting the individual histograms
    ips_->Write();
    ipsB_->Write();
    ipsC_->Write();
    ipsKs_->Write();
    ipsElse_->Write();
    ipsLight_->Write();
    ipsLambda_->Write();
    ipsInteraction_->Write();
    ipsConversion_->Write();

    // Creating a stack plot with all the contributions
    THStack * stack = new THStack("stack", "Positive and negative tag IPS distribution.");

    // Adding the contributions to the stack
    ipsB_->SetLineColor(kBlue);
    ipsB_->Scale(1./ipsB_->Integral());
    stack->Add(ipsB_);

    ipsKs_->SetLineColor(kMagenta);
    ipsKs_->Scale(1./ipsKs_->Integral());
    stack->Add(ipsKs_);

    ipsLambda_->SetLineColor(kMagenta);
    ipsLambda_->Scale(1./ipsLambda_->Integral());
    stack->Add(ipsLambda_);

    ipsConversion_->SetLineColor(kMagenta);
    ipsConversion_->Scale(1./ipsConversion_->Integral());
    stack->Add(ipsConversion_);

    ipsInteraction_->SetLineColor(kRed);
    ipsInteraction_->Scale(1./ipsInteraction_->Integral());
    stack->Add(ipsInteraction_);

    ipsC_->SetLineColor(kGreen);
    ipsC_->Scale(1./ipsC_->Integral());
    stack->Add(ipsC_);

    ipsLight_->SetLineColor(kBlack);
    ipsLight_->Scale(1./ipsLight_->Integral());
    stack->Add(ipsLight_);

    // ipsElse_->SetLineColor(kMagenta);
    // ipsElse_->Scale(1./ipsElse_->Integral());
    // stack->Add(ipsElse_);

    // Writing the stack
    stack->Write();
}

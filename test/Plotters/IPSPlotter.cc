
#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TTree.h"

#include "RecoBTag/PerformanceMeasurements/test/Plotters/IPSPlotter.h"

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
  ipsDisplaced_ = new TH1F("ipsDisplaced","Positive and negative tag IPS distribution for jet with displaced tracks", 60, -10., 10.);
  ipsConversion_ = new TH1F("ipsConversion","Positive and negative tag IPS distribution for jet with conversion tracks", 60, -10., 10.);
}

void IPSPlotter::Fill(BTagEvent * event)
{	
  // Loop over the jets
  for (Int_t i = 0; i<event->njets; ++i)
  {
    if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
      ips_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
    if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
      ips_->Fill(event->btag_NegTag_disc3D_2trk[i]);
  	  	
    if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::Bottom] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsB_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsB_->Fill(event->btag_NegTag_disc3D_2trk[i]);
    } 
    else if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::Ks] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsKs_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsKs_->Fill(event->btag_NegTag_disc3D_2trk[i]);
    }
    else if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::PhotonConversion] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsConversion_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsConversion_->Fill(event->btag_NegTag_disc3D_2trk[i]);      
    }
    else if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::Lambda] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsLambda_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsLambda_->Fill(event->btag_NegTag_disc3D_2trk[i]);      
    }
    else if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::Displaced] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsDisplaced_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsDisplaced_->Fill(event->btag_NegTag_disc3D_2trk[i]);      
    }
    else if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::Charm] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsC_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsC_->Fill(event->btag_NegTag_disc3D_2trk[i]);      
    }
    else if ( event->btag_TrkCounting_disc3D_2trk_is[i][BTagTrackEvent::Light] )
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsLight_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsLight_->Fill(event->btag_NegTag_disc3D_2trk[i]);      
    }
    else
    {
      if ( event->btag_TrkCounting_disc3D_2trk[i] >= 0 )
        ipsElse_->Fill(event->btag_TrkCounting_disc3D_2trk[i]);
      if ( event->btag_NegTag_disc3D_2trk[i] < 0 )
        ipsElse_->Fill(event->btag_NegTag_disc3D_2trk[i]);      
    }    
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
  ipsDisplaced_->Write();
  ipsConversion_->Write();

  // Creating a stack plot with all the contributions
  THStack * stack = new THStack("stack", "Positive and negative tag IPS distribution.");
  
  // Adding the contributions to the stack
  ipsB_->SetLineColor(kBlue);
  ipsB_->Scale(1./ipsB_->Integral());  stack->Add(ipsB_);

  ipsKs_->SetLineColor(kMagenta);
  ipsKs_->Scale(1./ipsKs_->Integral());
  stack->Add(ipsKs_);
  
  ipsLambda_->SetLineColor(kMagenta);
  ipsLambda_->Scale(1./ipsLambda_->Integral());
  stack->Add(ipsLambda_);
  
  ipsConversion_->SetLineColor(kMagenta);
  ipsConversion_->Scale(1./ipsConversion_->Integral());  
  stack->Add(ipsConversion_);

  ipsDisplaced_->SetLineColor(kRed);
  ipsDisplaced_->Scale(1./ipsDisplaced_->Integral());
  stack->Add(ipsDisplaced_);

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

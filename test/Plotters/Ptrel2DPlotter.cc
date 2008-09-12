
#include <iostream>

#include "TFile.h"
#include "THStack.h"
#include "TTree.h"

#include "Ptrel2DPlotter.h"

#include <cmath>

ClassImp(Ptrel2DPlotter)


void Ptrel2DPlotter::Book() 
{

  pTB_ = new TH1F(
    "pTB",
    "Muon pt",
    50, 0., 100.
  );    
  pTC_ = new TH1F(
    "pTC",
    "Muon pt", 
    50, 0., 100.
  );
  pTLight_ = new TH1F(
    "pTLight",
    "Muon pt",
    50, 0., 100.
  );

  deltaRB_ = new TH1F(
    "deltaRB",
    "Moun delta R",
    50, 0., 0.5
  );    
  deltaRC_ = new TH1F(
    "detalRC",
    "Moun delta R", 
    50, 0., 0.5
  );
  deltaRLight_ = new TH1F(
    "deltaRLight",
    "Moun delta R",
    50, 0., 0.5
  );    

  ntrkB_ = new TH1F(
    "ntrkB",
    "Number of tracks",
    41, -0.5, 40.5
  );
  ntrkC_ = new TH1F(
    "ntrkC",
    "Number of tracks",
    41, -0.5, 40.5
  );
  ntrkLight_ = new TH1F(
    "ntrkLight",
    "Number of tracks",
    41, -0.5, 40.5
  );
 
  pTrelVspTB_ = new TH2F(
    "pTrelVspTB",
    "Ptrel vs moun pt",
    40, 0., 100., 40, 0., 4.
  );    
  pTrelVspTC_ = new TH2F(
    "pTrelVspTC",
    "Ptrel vs moun pt", 
    40, 0., 100., 40, 0., 4.
  );
  pTrelVspTLight_ = new TH2F(
    "pTrelVspTLight",
    "Ptrel vs moun pt",
    40, 0., 100., 40, 0., 4.
  );    

  pTrelVsDeltaRB_ = new TH2F(
    "pTrelVsDeltaRB",
    "Ptrel vs muon delta R",
    50, 0., 0.5, 40, 0., 4.
  );    
  pTrelVsDeltaRC_ = new TH2F( 
    "pTrelVsDeltaRC",
    "Ptrel vs muon delta R", 
    50, 0., 0.5, 40, 0., 4.
  );
  pTrelVsDeltaRLight_ = new TH2F(
    "pTrelVsDeltaRLight",
    "Ptrel vs muon delta R",
    50, 0., 0.5, 40, 0., 4.
  );    

  pTrelVsNtrkB_ = new TH2F(
    "pTrelVsNtrkB",
    "Ptrel vs number of tracks",
    41, -0.5, 40.5, 40, 0., 4.
  );    
  pTrelVsNtrkC_ = new TH2F( 
    "pTrelVsNtrkC",
    "Ptrel vs number of tracks", 
    41, -0.5, 40.5, 40, 0., 4.
  );
  pTrelVsNtrkLight_ = new TH2F(
    "pTrelVsNtrkLight",
    "Ptrel vs number of tracks",
    41, -0.5, 40.5, 40, 0., 4.
  );    
 
  pTrelVsIps2trkB_ = new TH2F(
    "pTrelVsIps2trkB",
    "Ptrel vs ips of 2 trk",
    40, -10., 10., 40, 0., 4.
  );  
  pTrelVsIps2trkC_ = new TH2F(
    "pTrelVsIps2trkC",
    "Ptrel vs ips of 2 trk", 
    40, -10., 10., 40, 0., 4.
  );
  pTrelVsIps2trkLight_ = new TH2F(
    "pTrelVsIps2trkLight",
    "Ptrel vs ips of 2 trk",
    40, -10., 10., 40, 0., 4.
  );    

  pTrelVsAbsIps2trkB_ = new TH2F(
    "pTrelVsAbsIps2trkB",
    "Ptrel vs ips of 2 trk",
    20, 0., 10., 40, 0., 4.
  );  
  pTrelVsAbsIps2trkC_ = new TH2F(
    "pTrelVsAbsIps2trkC",
    "Ptrel vs ips of 2 trk", 
    20, 0., 10., 40, 0., 4.
  );
  pTrelVsAbsIps2trkLight_ = new TH2F(
    "pTrelVsAbsIps2trkLight",
    "Ptrel vs ips of 2 trk",
    20, 0., 10., 40, 0., 4.
  );    

}

void Ptrel2DPlotter::Fill(BTagEvent * event)
{	
  // Loop over the jets
  for (Int_t jetIndex = 0; jetIndex < event->njets; ++jetIndex)
  {
  	float maxPt = 0.;
  	std::size_t muonIndex = 0;
  	std::size_t maxPtIndex = 0;
  	
    if (!event->lepton[jetIndex].pt.empty())
    {
      // Look for the muon with maximum pt.
      for (; muonIndex < event->lepton[jetIndex].pt.size(); ++muonIndex)
      {
        double pt = event->lepton[jetIndex].pt[muonIndex];
        if (pt > maxPt)
        {
          maxPt = pt;
          maxPtIndex = muonIndex; 
        }
      }
      	
      if ( event->jet_flavour[jetIndex] == 5 )
      {
        pTB_->Fill(
          event->lepton[jetIndex].pt[maxPtIndex]
        );        
        deltaRB_->Fill(
          event->lepton[jetIndex].jet_deltaR[maxPtIndex]
        );
        ntrkB_->Fill(
          event->jet_ntrks[jetIndex]
        );
        pTrelVsDeltaRB_->Fill(
          event->lepton[jetIndex].jet_deltaR[maxPtIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]        
        );
        pTrelVspTB_->Fill(
          event->lepton[jetIndex].pt[maxPtIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );        
        pTrelVsNtrkB_->Fill(
          event->jet_ntrks[jetIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );        
        pTrelVsIps2trkB_->Fill(
          event->btag_TrkCounting_disc3D_2trk[jetIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );
        pTrelVsAbsIps2trkB_->Fill(
          event->btag_TrkCounting_disc3D_2trk[jetIndex],
          fabs(event->lepton[jetIndex].jet_ptrel[maxPtIndex])
        );
      }
      else if ( event->jet_flavour[jetIndex] == 4 )
      {
        pTC_->Fill(
          event->lepton[jetIndex].pt[maxPtIndex]
        );
        deltaRC_->Fill(
          event->lepton[jetIndex].jet_deltaR[maxPtIndex]
        );
        ntrkC_->Fill(
          event->jet_ntrks[jetIndex]
        );
        pTrelVsDeltaRC_->Fill(
          event->lepton[jetIndex].jet_deltaR[maxPtIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]        
        );
        pTrelVspTC_->Fill(
          event->lepton[jetIndex].pt[maxPtIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );        
        pTrelVsNtrkC_->Fill(
          event->jet_ntrks[jetIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );        
        pTrelVsIps2trkC_->Fill(
          event->btag_TrkCounting_disc3D_2trk[jetIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );
        pTrelVsAbsIps2trkC_->Fill(
          event->btag_TrkCounting_disc3D_2trk[jetIndex],
          fabs(event->lepton[jetIndex].jet_ptrel[maxPtIndex])
        );
      }
      else
      {
        pTLight_->Fill(
          event->lepton[jetIndex].pt[maxPtIndex]
        );
        deltaRLight_->Fill(
          event->lepton[jetIndex].jet_deltaR[maxPtIndex]
        );   
        ntrkLight_->Fill(
          event->jet_ntrks[jetIndex]
        );
        pTrelVsDeltaRLight_->Fill(
          event->lepton[jetIndex].jet_deltaR[maxPtIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]        
        );
      	pTrelVspTLight_->Fill(
          event->lepton[jetIndex].pt[maxPtIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );        
        pTrelVsNtrkLight_->Fill(
          event->jet_ntrks[jetIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );        
        pTrelVsIps2trkLight_->Fill(
          event->btag_TrkCounting_disc3D_2trk[jetIndex],
          event->lepton[jetIndex].jet_ptrel[maxPtIndex]
        );      	
        pTrelVsAbsIps2trkLight_->Fill(
          event->btag_TrkCounting_disc3D_2trk[jetIndex],
          fabs(event->lepton[jetIndex].jet_ptrel[maxPtIndex])
        );
      }
    }        
  }
}

void Ptrel2DPlotter::Write()
{
  std::cout << "Saving Histograms" << std::endl;

  TFile file("Ptrel2DPlotter.root", "RECREATE");

  pTB_->Write();
  pTC_->Write();
  pTLight_->Write();

  deltaRB_->Write();
  deltaRC_->Write();
  deltaRLight_->Write();

  ntrkB_->Write();
  ntrkC_->Write();
  ntrkLight_->Write();

  pTrelVsDeltaRB_->Write();
  pTrelVsDeltaRC_->Write();
  pTrelVsDeltaRLight_->Write();

  pTrelVspTB_->Write();
  pTrelVspTC_->Write();
  pTrelVspTLight_->Write();

  pTrelVsNtrkB_->Write();
  pTrelVsNtrkC_->Write();
  pTrelVsNtrkLight_->Write();
  
  pTrelVsIps2trkB_->Write();
  pTrelVsIps2trkC_->Write();
  pTrelVsIps2trkLight_->Write();
  
  pTrelVsAbsIps2trkB_->Write();
  pTrelVsAbsIps2trkC_->Write();
  pTrelVsAbsIps2trkLight_->Write();  
}

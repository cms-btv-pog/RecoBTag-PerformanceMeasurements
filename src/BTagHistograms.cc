/**_________________________________________________________________
   class:   BTagHistograms.cc


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagHistograms.cc,v 1.7 2008/08/11 06:06:22 bazterra Exp $

________________________________________________________________**/


#ifdef NOSCRAMV
#include "BTagHistograms.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagHistograms.h"
#endif

#include "TF1.h"

#include<iostream>

//_______________________________________________________________
BTagHistograms::BTagHistograms()
{

}

//_______________________________________________________________
BTagHistograms::~BTagHistograms()
{

    this->DeleteHisto();

}

//_______________________________________________________________
void BTagHistograms::Init(TString type, TString suffix1, TString suffix2)
{

    const int nptarray = 11;
    const int netaarray = 10;
    Double_t jetptbins[nptarray] = {30., 40., 50., 60., 70., 80, 90., 100., 120., 140., 230.};
    Double_t jetetabins[netaarray] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5};

    if (suffix1 != "") suffix1 = "_" + suffix1;
    if (suffix2 != "") suffix1 += "_" + suffix2;


    if ( type == "efficiencies" )
    {

        h1["jet_pt"+suffix1] = new TH1D("jet_pt"+suffix1,"Jet p_{T} [GeV/c]",nptarray-1,jetptbins);
        h1["jet_pt_b"+suffix1] = new TH1D("jet_pt_b"+suffix1,"Jet p_{T} [GeV/c]",nptarray-1,jetptbins);
        h1["jet_pt_c"+suffix1] = new TH1D("jet_pt_c"+suffix1,"Jet p_{T} [GeV/c]",nptarray-1,jetptbins);
        h1["jet_pt_udsg"+suffix1] = new TH1D("jet_pt_udsg"+suffix1,"Jet p_{T} [GeV/c]",nptarray-1,jetptbins);

        h1["jet_eta"+suffix1] = new TH1D("jet_eta"+suffix1,"Jet |#eta|",netaarray-1,jetetabins);
        h1["jet_eta_b"+suffix1] = new TH1D("jet_eta_b"+suffix1,"Jet |#eta|",netaarray-1,jetetabins);
        h1["jet_eta_c"+suffix1] = new TH1D("jet_eta_c"+suffix1,"Jet |#eta|",netaarray-1,jetetabins);
        h1["jet_eta_udsg"+suffix1] = new TH1D("jet_eta_udsg"+suffix1,"Jet |#eta|",netaarray-1,jetetabins);

        // WHAT THE HECK ??
        // for (std::map<TString,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih)
        // {
        //     TH1 *htemp = ih->second;
        // }
    }
    if ( type == "ptrel")
    {

        h1["jet_ptrel"+suffix1] = new TH1D("jet_ptrel"+suffix1,"p_{Trel} [GeV/c]",50,0.,5.);
        h1["jet_ptrel_b"+suffix1] = new TH1D("jet_ptrel_b"+suffix1,"p_{Trel} [GeV/c]",50,0.,5.);
        h1["jet_ptrel_c"+suffix1] = new TH1D("jet_ptrel_c"+suffix1,"p_{Trel} [GeV/c]",50,0.,5.);
        h1["jet_ptrel_udsg"+suffix1] = new TH1D("jet_ptrel_udsg"+suffix1,"p_{Trel} [GeV/c]",50,0.,5.);

        for (std::map<TString,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih)
        {
            TH1 *htemp = ih->second;
            try
            {
                htemp->Sumw2();
            }
            catch (...) {}
        }
    }
    else if ( type == "n")
    {

        h2["n_pT"+suffix1] = new TH2D("n_pT"+suffix1,"muon-jet+away-jet pT vs pTrel",nptarray-1,jetptbins,50,0.,5.);
        h2["n_eta"+suffix1] = new TH2D("n_eta"+suffix1,"muon-jet+away-jet eta vs ptrel",netaarray-1,jetetabins,50,0.,5.);

        for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
        {
            TH2 *htemp = ih->second;
            try
            {
                htemp->Sumw2();
            }
            catch (...) {}

        }
    }
    else if ( type == "ntag")
    {
        h2["ntag_pT"+suffix1] = new TH2D("ntag_pT"+suffix1,"tagged muon-jet+away-jet pT vs pTrel",nptarray-1,jetptbins,50,0.,5.);
        h2["ntag_eta"+suffix1] = new TH2D("ntag_eta"+suffix1,"tagged muon-jet+away-jet eta vs ptrel",netaarray-1,jetetabins,50,0.,5.);

        for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
        {
            TH2 *htemp = ih->second;
            try
            {
                htemp->Sumw2();
            }
            catch (...) {}

        }
    }
    else if ( type == "p")
    {

        h2["p_pT"+suffix1] = new TH2D("p_pT"+suffix1,"muon-jet+tagged-away-jet pT vs pTrel",nptarray-1,jetptbins,50,0.,5.);
        h2["p_eta"+suffix1] = new TH2D("p_eta"+suffix1,"muon-jet+tagged-away-jet eta vs ptrel",netaarray-1,jetetabins,50,0.,5.);

        for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
        {
            TH2 *htemp = ih->second;
            try
            {
                htemp->Sumw2();
            }
            catch (...) {}

        }
    }
    else if ( type == "ptag")
    {
        h2["ptag_pT"+suffix1] = new TH2D("ptag_pT"+suffix1,"tagged muon-jet+tagged-away-jet pT vs pTrel",nptarray-1,jetptbins,50,0.,5.);
        h2["ptag_eta"+suffix1] = new TH2D("ptag_eta"+suffix1,"tagged muon-jet+tagged-away-jet eta vs pTrel",netaarray-1,jetetabins,50,0.,5.);

        for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
        {
            TH2 *htemp = ih->second;
            try
            {
                htemp->Sumw2();
            }
            catch (...) {}
        }
    }

    for (std::map<TString,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih)
    {
        TH1 *htemp = ih->second;
        htemp->SetXTitle( htemp->GetTitle() );
    }
    for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
    {
        TH2 *htemp = ih->second;
        htemp->SetXTitle( htemp->GetTitle() );
    }
}

//_______________________________________________________________
void BTagHistograms::Fill1d(TString name, Double_t x, Double_t weight)
{

    h1[name]->Fill(x,weight);
}

//_______________________________________________________________
void BTagHistograms::Fill2d(TString name, Double_t x, Double_t y, Double_t weight)
{

    h2[name]->Fill(x,y,weight);

}

//_______________________________________________________________
void BTagHistograms::Fit(TString name, Double_t mean)
{

}
//_______________________________________________________________
void BTagHistograms::Save()
{

    for (std::map<TString,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih)
    {
        //std::cout << "get histo: " << std::endl;
        TH1D *htemp = (TH1D*) ih->second;
        //std::cout << htemp->Print() << std::endl;
        if (htemp->GetEntries() > 0 ) htemp->Write();
    }
    for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
    {
        //std::cout << "get histo: " << std::endl;
        TH2D *htemp = (TH2D*) ih->second;
        //std::cout << htemp->Print() << std::endl;
        if (htemp->GetEntries() > 0 ) htemp->Write();
    }
}

//_______________________________________________________________
void BTagHistograms::SaveToFile(TString filename)
{

    foutfile = new TFile(filename,"RECREATE");
    for (std::map<TString,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih)
    {
        //std::cout << "get histo: " << std::endl;
        TH1D *htemp = (TH1D*) ih->second;
        //std::cout << htemp->Print() << std::endl;
        htemp->Write();
    }
    for (std::map<TString,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih)
    {
        //std::cout << "get histo: " << std::endl;
        TH2D *htemp = (TH2D*) ih->second;
        //std::cout << htemp->Print() << std::endl;
        htemp->Write();
    }

    foutfile->Write();
    foutfile->Close();

}

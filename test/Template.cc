#define Template_cxx
#include "Template.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>
#include "TLegend.h"

#include <iostream>
#include <iomanip>


void Template::Loop()
{

    TH1D *h1 = new TH1D("h1","jet E_{T} [GeV]",20,20.,180.);
    TH1D *h2 = new TH1D("h2","muon p_{T} [GeV]",20,0.,80.);

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    std::cout << " Total entries = " << fChain->GetEntries() << std::endl;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        if ( jentry%100 == 0 ) std::cout << " processing entry: " << jentry << std::endl;


        std::vector<float> jet_pt_vec = fS8evt->jet_pt;
        int vec_size = jet_pt_vec.size();

        for ( int ijet =0; ijet != vec_size; ++ijet)
        {

            h1->Fill( fS8evt->jet_et[ijet] );

            if ( fS8evt->jet_hasLepton[ijet] == 1 )
            {

                BTagLeptonEvent Muons = fS8evt->lepton[ijet];

                int mu_size = Muons.pt.size();

                for ( int imu = 0; imu != mu_size; ++imu )
                {

                    h2->Fill( Muons.pt[imu] );

                }
            }
        }// end jet loop

    } // end loop



    //______________________________________________________

    TCanvas *cv_pt = new TCanvas("jet_pt","jet_pt",700,700);
    h1->Draw();
    cv_pt->Update();

    //______________________________________________________

    TCanvas *cv_mupt = new TCanvas("muon_pt","muon_pt",700,700);
    h2->Draw();
    cv_mupt->Update();

}

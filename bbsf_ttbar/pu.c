#include <TH1D.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <TFile.h>

//https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py

void pu_addgen()
{
   std::cout << "adding gen hist" << std::endl;
  
   TH1D * h_gen = new TH1D("h_gen", ";probFunctionVariable;probValue", 100, 0., 100.);
   ifstream myfile;
   myfile.open("truepu.dat");
   for (int i = 0; i < 99; ++i) {
      double w;
      myfile >> w;
      h_gen->SetBinContent(i+1, w);
      h_gen->SetBinError(i+1, 0.);
   }
   myfile.close();
  
   TFile * f = new TFile("pu.root", "RECREATE"); 
   h_gen->Write(); 
   f->Close();
}

void pu_addmc(const TString label)
{
   std::cout << "adding " + label + " hists using getTrueNumInteractions()..." << std::endl;

   TH1D * h_denom = new TH1D("h_"+label, ";getTrueNumInteractions();events / 1", 100, 0., 100.);

   ifstream infile("./filelistdir/list_"+label);
   std::string instring;
   const TString dir = "/eos/uscms/store/user/fojensen/";
   int n = 0;
   while (infile >> instring) {
      TFile * f = TFile::Open(dir + instring);
      TH1D * h = (TH1D*)f->Get("plotPUTrue/h_getTrueNumInteractions");
      h_denom->Add(h);
      f->Close();
      ++n;
   }
   infile.close();
   std::cout << n << " hists added" << std::endl;

   TFile * f = TFile::Open("pu.root", "UPDATE");
   h_denom->Write();
   f->Close();
}

void pu_addsf(const TString label)
{
   std::cout << "adding sfs to " + label + "..." << std::endl;

   const TString filenames[3] = {"RA2DownPileupHistogram", "RA2CentralPileupHistogram", "RA2UpPileupHistogram"};
   const TString filestyle[3] = {"down", "central", "up"};
   TH1D *h_pu[3], *h_num[3];
   for (int i = 0; i < 3; ++i) {
      TFile * tempfile = TFile::Open("./RA2PileupHistograms/"+filenames[i]+".root");
      h_pu[i] = (TH1D*)tempfile->Get("pileup");
      h_num[i] = (TH1D*)h_pu[i]->Clone("h_num_"+filestyle[i]);
      h_num[i]->Scale(1./h_num[i]->Integral());
   }
   
   TFile * f = TFile::Open("pu.root", "UPDATE");
   TH1D * h_denom =  (TH1D*)f->Get("h_"+label);
   h_denom->Scale(1./h_denom->Integral());
   
   TH1D * sf[3];
   for (int i = 0; i < 3; ++i) {
      sf[i] = new TH1D("sf_"+filestyle[i]+"_"+label, ";;pile-up weights", 100, 0., 100.);
   }

   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 100; ++j) {
         const double dataval = h_num[i]->GetBinContent(j+1);
         const double mcval = h_denom->GetBinContent(j+1);
         if (dataval>0. && mcval>0.) {
            const double sfval = dataval/mcval;
            const double dataerr = h_num[i]->GetBinError(j+1);
            const double mcerr = h_denom->GetBinError(j+1);
            const double sferrval = sfval * sqrt(pow(dataerr/dataval, 2)+pow(mcerr/mcval, 2));
            sf[i]->SetBinContent(j+1, sfval);
            sf[i]->SetBinError(j+1, sferrval);
         } else {
            sf[i]->SetBinContent(j+1, 0.);
            sf[i]->SetBinError(j+1, 0.);
         }
      }
   }

   sf[0]->Write();
   sf[1]->Write();
   sf[2]->Write();
   f->Close();
}

void pu()
{
   pu_addgen();

   pu_addmc("TTToSemiLeptonic");
   pu_addsf("TTToSemiLeptonic");

   pu_addmc("TTToHadronic");
   pu_addsf("TTToHadronic");

   pu_addmc("TTTo2L2Nu");
   pu_addsf("TTTo2L2Nu");

   pu_addmc("W4JetsToLNu");
   pu_addsf("W4JetsToLNu");

   pu_addmc("W3JetsToLNu");
   pu_addsf("W3JetsToLNu");

   pu_addmc("WWToLNuQQ_NNPDF31");
   pu_addsf("WWToLNuQQ_NNPDF31");

   pu_addmc("ST_tW_antitop_5f");
   pu_addsf("ST_tW_antitop_5f");

   pu_addmc("ST_tW_top_5f");
   pu_addsf("ST_tW_top_5f");
}


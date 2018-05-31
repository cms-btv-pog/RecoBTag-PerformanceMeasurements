#include <TPaveText.h>
#include <TLatex.h>
#include <TExec.h>
#include <TCanvas.h>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>

void drawtext(const TString name)
{  
   TH1D * h = (TH1D*)gPad->GetListOfPrimitives()->FindObject(name);
   const int n = h->GetNbinsX();
   TLatex *l[n];
   for (int i = 0; i < n; i++) {
      const double x = h->GetBinCenter(i+1);
      const double y = h->GetBinContent(i+1);
      const double err = h->GetBinError(i+1);
      l[i] = new TLatex(x-0.5, y+0.1, Form("%.3f #pm %.3f",y, err));
      l[i]->SetTextSize(0.03);
      l[i]->SetTextAlign(11);
      l[i]->SetTextAngle(30);
      l[i]->Paint();
   }
}

void sf_addsf()
{
   const TString labels[4] = {"0p3", "0p6", "0p8", "0p9"};
   
   TH1D * h[4];
   for (int i = 0; i < 4; ++i) {
      TFile * f = TFile::Open("./sf_rootfiles/sf_"+labels[i]+"_SingleMuon_puWeight_noTopPtWeight_bkgNormNominal.root");
      h[i] = (TH1D*)f->Get("h_sf");
      h[i]->SetName("h_sf_"+labels[i]);
   }

   const int nbins = 3;
   const double errsys[3] = {0.02, 0.04, 0.06};   

   for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 3; ++j) {
         const double errstat = h[i]->GetBinError(j+1);
         const double err = sqrt(errstat*errstat + errsys[j]*errsys[j]);
         h[i]->SetBinError(j+1, err);
      }
   }

   TFile * f = new TFile("./sf_rootfiles/sf.root", "RECREATE");
   for (int i = 0; i < 4; ++i) {
      h[i]->SetStats(0);
      h[i]->Write("sf_"+labels[i]);
   }
   f->Close();

   const TString titles[4] = {"0.3", "0.6", "0.8", "0.9"};

   TExec *dth0 = new TExec("dth0", "drawtext(\"h_sf_0p3\");");
   h[0]->GetListOfFunctions()->Add(dth0);

   TExec *dth1 = new TExec("dth1", "drawtext(\"h_sf_0p6\");");
   h[1]->GetListOfFunctions()->Add(dth1);

   TExec *dth2 = new TExec("dth2", "drawtext(\"h_sf_0p8\");");
   h[2]->GetListOfFunctions()->Add(dth2);

   TExec *dth3 = new TExec("dth3", "drawtext(\"h_sf_0p9\");");
   h[3]->GetListOfFunctions()->Add(dth3);

   TCanvas * c = new TCanvas("c", "", 800, 800);
   c->Divide(2, 2);
   for (int i = 0; i < 4; ++i) {
      c->cd(i+1);
      h[i]->SetStats(0);
      h[i]->SetTitle(titles[i]);
      h[i]->SetMarkerStyle(20);
      h[i]->SetLineColor(1);
      h[i]->SetMarkerColor(1);
      h[i]->Draw("PE");
      h[i]->SetMaximum(1.4);
      h[i]->SetMinimum(0.6);
   }
   c->SaveAs("./plots/sf.pdf");
   c->SaveAs("./plots/sf.png");
}


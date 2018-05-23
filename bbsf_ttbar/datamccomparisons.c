#include "labelCMS.h"
#include <TChain.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TH1D.h>
#include <TLegend.h>
#include "addOverflow.h"
#include <THStack.h>

void drawcomparisons(const TH1D * h, const TString var, const bool dolog=false)
{
   std::cout << "running on " << TString(h->GetName()) << std::endl;
   const bool mcscale = true;

   const TString baseline = "1>0";

   const int n = 9;
   TH1D * hist[n];
   TString names[n], labels[n];
 
   names[0] = "SingleMuon";        labels[0] = "Single-Muon";
   names[1] = "TTToSemiLeptonic";  labels[1] = "t#bar{t}#rightarrow2ql#nu";
   names[2] = "TTToHadronic";      labels[2] = "t#bar{t}#rightarrow4q";
   names[3] = "TTTo2L2Nu";         labels[3] = "t#bar{t}#rightarrow2l2#nu";
   names[4] = "W3JetsToLNu";       labels[4] = "W+3Jets";
   names[5] = "W4JetsToLNu";       labels[5] = "W+4Jets";
   names[6] = "WWToLNuQQ_NNPDF31"; labels[6] = "WW";
   names[7] = "ST_tW_antitop_5f";  labels[7] = "single-top #bar{t}W^{+}";
   names[8] = "ST_tW_top_5f";      labels[8] = "single-top tW^{-}";

   for (int i = 0; i < n; ++i) {
      hist[i] = (TH1D*)h->Clone(TString(h->GetName()) + "_" + TString::Itoa(i, 10));
      TChain chain4("treeAK4");
      chain4.Add("./skims/" + names[i] + ".root"); 
      TChain chain8("treeAK8");
      chain8.Add("./skims/" + names[i] + ".root");
      chain4.AddFriend(&chain8);
      char buffer[1000];
      if (i==0) {
         sprintf(buffer, "(%s)",  baseline.Data());
      } else if (i==1||i==2||i==3) {
         sprintf(buffer, "41860. * xsWeight * puWeight * ttbar_ptweight * (%s)", baseline.Data());
         //sprintf(buffer, "41860. * xsWeight * ttbar_ptweight * (%s)", baseline.Data());
      } else {
         sprintf(buffer, "41860. * xsWeight * puWeight * (%s)", baseline.Data());
         //sprintf(buffer, "41860. * xsWeight * (%s)", baseline.Data());
      }
      chain4.Project(hist[i]->GetName(), var, buffer);
      addOverflow(hist[i]);
   }
  
   TLegend l(0.25, 0.75, 0.85, 0.875);
   l.SetBorderSize(0);
   l.SetNColumns(2);
   hist[0]->SetMarkerStyle(20);
   l.AddEntry(hist[0], labels[0], "P");

   std::cout << "data integral: " << hist[0]->Integral() << std::endl;

   double mcint = 0;
   for (int i = 1; i < n; ++i) {
      hist[i]->SetFillColor(i+1);
      l.AddEntry(hist[i], labels[i], "F");
      std::cout << labels[i] << " " << hist[i]->Integral() << std::endl;
      mcint += hist[i]->Integral();
   }
   
   if (mcscale) {
      const double scale = hist[0]->Integral()/mcint;
      std::cout << "scaling mc by: " << scale << std::endl;
      for (int i = 1; i < n; ++i) {
         hist[i]->Scale(scale);
      }
   }

   TH1D * mcsum = (TH1D*)hist[1]->Clone("mcsum");
   for (int i = 2; i < n; ++i) mcsum->Add(hist[i]);
   
   TH1D * ratio = (TH1D*)hist[0]->Clone("ratio");
   ratio->Divide(mcsum);
   ratio->SetMarkerStyle(20);
   ratio->GetYaxis()->SetTitle("data / mc");
   ratio->SetStats(0);

   THStack stack;
   for (int i = 1; i < n; ++i) stack.Add(hist[i]);
   char buffer[1000];
   sprintf(buffer, ";%s;%s", hist[0]->GetXaxis()->GetTitle(), hist[0]->GetYaxis()->GetTitle());
   stack.SetTitle(buffer);

   TCanvas c("c", "", 800, 400);
   c.Divide(2, 1);
   TPad * pad1 = (TPad*)c.cd(1);
   stack.Draw("HIST");
   hist[0]->Draw("P, SAME, E");   
   l.Draw();
   stack.SetMinimum(0.);
   stack.SetMaximum(1.5 * TMath::Max(stack.GetMaximum(), hist[0]->GetMaximum()));
   if (dolog) {
      stack.SetMinimum(10.);
      stack.SetMaximum(10. * TMath::Max(stack.GetMaximum(), hist[0]->GetMaximum()));
      pad1->SetLogy(); 
   }
   labelData(true);
   c.cd(2);
   ratio->Draw("PE");
   ratio->SetMinimum(0.);
   ratio->SetMaximum(2.);
   TLine line(hist[0]->GetBinLowEdge(1), 1., hist[0]->GetBinLowEdge(hist[0]->GetNbinsX()+1), 1.);
   line.SetLineColor(2);
   line.SetLineStyle(2);
   line.Draw();
   labelData(true);
   const TString savetag = "plots/datamccomparisons_"+TString(h->GetName());
   c.SaveAs(savetag + ".pdf");
   c.SaveAs(savetag + ".png");
   
   std::cout << "" << std::endl;
}

void datamccomparisons()
{
   TH1D * muon_pt = new TH1D("muon_pt", ";muon p_{T} [GeV];events / 25 GeV", 20, 50., 550.);
   drawcomparisons(muon_pt, "PatMuon_pt[idx_muon]", true);

   TH1D * hadJetAK8_pt = new TH1D("hadJetAK8_pt", ";AK8 jet p_{T} [GeV];events / 25 GeV", 20, 250., 750.);
   drawcomparisons(hadJetAK8_pt, "FatJetInfo.Jet_pt[idx_ak8]", true);

   TH1D * hadJetAK8_mass = new TH1D("hadJetAK8_mass", ";AK8 jet mass [GeV];events / 7.5 GeV", 20, 50., 200.);
   drawcomparisons(hadJetAK8_mass, "FatJetInfo.Jet_massSoftDrop[idx_ak8]");

   TH1D * hadJetAK8_t2t1 = new TH1D("hadJetAK8_t2t1", ";AK8 jet #tau_{2}/#tau_{1};events / 0.05", 20, 0., 1.);
   drawcomparisons(hadJetAK8_t2t1, "FatJetInfo.Jet_tau2[idx_ak8]/FatJetInfo.Jet_tau1[idx_ak8]");

   TH1D * hadJetAK8_disc = new TH1D("hadJetAK8_disc", ";AK8 jet bb-discriminator;events / 0.1", 20, -1., 1.);
   drawcomparisons(hadJetAK8_disc, "FatJetInfo.Jet_DoubleSV[idx_ak8]");

   TH1D * hadJetAK4_pt = new TH1D("hadJetAK4_pt", ";had-side AK4 jet p_{T} [GeV];events / 20 GeV", 20, 30., 430.);
   drawcomparisons(hadJetAK4_pt, "Jet_pt[idx_hadak4]", true);

   TH1D * hadJetAK4_disc = new TH1D("hadJetAK4_disc", ";had-side AK4 jet b disc;events / 0.05", 20, 0., 1.);
   drawcomparisons(hadJetAK4_disc, "Jet_DeepCSVb[idx_hadak4]");

   TH1D * lepJetAK4_pt = new TH1D("lepJetAK4_pt", ";lep-side AK4 jet p_{T} [GeV];events / 20 GeV", 20, 30., 430.);
   drawcomparisons(lepJetAK4_pt, "Jet_pt[idx_lepak4]", true);

   TH1D * lepJetAK4_disc = new TH1D("lepJetAK4_disc", ";lep-side AK4 jet b disc;events / 0.05", 20, 0., 1.);
   drawcomparisons(lepJetAK4_disc, "Jet_DeepCSVb[idx_lepak4]");
}


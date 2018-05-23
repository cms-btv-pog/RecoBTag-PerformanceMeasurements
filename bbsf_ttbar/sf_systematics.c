#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void sf_systematics_toppt()
{
   TFile * f0 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuon_puWeight_noTopPtWeight.root");
   TH1D * h0 = (TH1D*)f0->Get("h_sf");
   h0->SetLineWidth(2);
   h0->SetLineColor(8);
   h0->SetMarkerColor(8);
   h0->SetMarkerStyle(20);
   
   TFile * f1 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuon_puWeight_topPtWeight.root");
   TH1D * h1 = (TH1D*)f1->Get("h_sf");
   h1->SetLineWidth(2);
   h1->SetLineColor(9);
   h1->SetMarkerColor(9);
   h1->SetMarkerStyle(20);
   
   TCanvas * c = new TCanvas("c", "", 400, 400);
   h0->Draw("PE");
   h1->Draw("PE, SAME");
   h0->SetStats(0);
   h0->SetTitle(";AK8 jet p_{T} [GeV];scale factor");
   h0->SetMinimum(0.85);
   h0->SetMaximum(1.01);
  
   TLegend * l = new TLegend(0.4, 0.2, 0.8, 0.35);
   l->SetBorderSize(0);
   l->SetHeader("top p_{T} weights");
   l->AddEntry(h1, "with", "P");
   l->AddEntry(h0, "without", "P");
   l->Draw();

   c->SaveAs("plots/sf_toppt_0p3.pdf");
   c->SaveAs("plots/sf_toppt_0p3.png");

   /*c->cd(2);
   r->Draw("PE");
   r->SetStats(0);*/
}

void sf_systematics_pu()
{
   TFile * f0 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuon_puWeight_noTopPtWeight.root");
   TH1D * h0 = (TH1D*)f0->Get("h_sf");
   h0->SetLineWidth(2);
   h0->SetLineColor(9);
   h0->SetMarkerColor(9);
   h0->SetMarkerStyle(20);
   
   TFile * f1 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuon_puWeight_up_noTopPtWeight.root");
   TH1D * h1 = (TH1D*)f1->Get("h_sf");
   h1->SetLineWidth(2);
   h1->SetLineColor(6);
   h1->SetMarkerColor(6);
   h1->SetMarkerStyle(20);

   TFile * f2 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuon_puWeight_down_noTopPtWeight.root");
   TH1D * h2 = (TH1D*)f2->Get("h_sf");
   h2->SetLineWidth(2);
   h2->SetLineColor(7);
   h2->SetMarkerColor(7);
   h2->SetMarkerStyle(20);

   TFile * f3 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuon_noPuWeight_noTopPtWeight.root");
   TH1D * h3 = (TH1D*)f3->Get("h_sf");
   h3->SetLineWidth(2);
   h3->SetLineColor(8);
   h3->SetMarkerColor(8);
   h3->SetMarkerStyle(20);

   TCanvas * c = new TCanvas("c", "", 400, 400);
   h0->Draw("PE");
   h1->Draw("PE, SAME");
   h2->Draw("PE, SAME");
   h3->Draw("PE, SAME");
   h0->SetStats(0);
   h0->SetTitle(";AK8 jet p_{T} [GeV];scale factor");
   h0->SetMinimum(0.85);
   h0->SetMaximum(1.01);

   TLegend * l = new TLegend(0.4, 0.2, 0.8, 0.35);
   l->SetBorderSize(0);
   l->SetNColumns(2);
   l->SetHeader("pile-up weights");
   l->AddEntry(h0, "nominal", "P");
   l->AddEntry(h1, "+#sigma", "P");
   l->AddEntry(h2, "-#sigma", "P");
   l->AddEntry(h3, "none", "P");
   l->Draw();

   c->SaveAs("plots/sf_pu_0p3.pdf");
   c->SaveAs("plots/sf_pu_0p3.png");
}

void sf_systematics_dataset()
{
   TFile * f0 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuonB_puWeight_noTopPtWeight.root");
   TH1D * h0 = (TH1D*)f0->Get("h_sf");
   h0->SetLineWidth(2);
   h0->SetLineColor(2);
   h0->SetMarkerColor(2);
   h0->SetMarkerStyle(20);
   
   TFile * f1 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuonC_puWeight_noTopPtWeight.root");
   TH1D * h1 = (TH1D*)f1->Get("h_sf");
   h1->SetLineWidth(2);
   h1->SetLineColor(3);
   h1->SetMarkerColor(3);
   h1->SetMarkerStyle(20);

   TFile * f2 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuonD_puWeight_noTopPtWeight.root");
   TH1D * h2 = (TH1D*)f2->Get("h_sf");
   h2->SetLineWidth(2);
   h2->SetLineColor(4);
   h2->SetMarkerColor(4);
   h2->SetMarkerStyle(20);

   TFile * f3 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuonE_puWeight_noTopPtWeight.root");
   TH1D * h3 = (TH1D*)f3->Get("h_sf");
   h3->SetLineWidth(2);
   h3->SetLineColor(6);
   h3->SetMarkerColor(6);
   h3->SetMarkerStyle(20);

   TFile * f4 = TFile::Open("sf_rootfiles/sf_0p3_SingleMuonF_puWeight_noTopPtWeight.root");
   TH1D * h4 = (TH1D*)f4->Get("h_sf");
   h4->SetLineWidth(2);
   h4->SetLineColor(7);
   h4->SetMarkerColor(7);
   h4->SetMarkerStyle(20);

   TCanvas * c = new TCanvas("c", "", 400, 400);
   h0->Draw("PE");
   h1->Draw("PE, SAME");
   h2->Draw("PE, SAME");
   h3->Draw("PE, SAME");
   h4->Draw("PE, SAME");
   h0->SetStats(0);
   h0->SetTitle(";AK8 jet p_{T} [GeV];scale factor");
   h0->SetMinimum(0.8);
   h0->SetMaximum(1.05);

   TLegend * l = new TLegend(0.4, 0.2, 0.8, 0.35);
   l->SetBorderSize(0);
   l->SetNColumns(2);
   l->SetHeader("SingleMuon");
   l->AddEntry(h0, "B", "P");
   l->AddEntry(h1, "C", "P");
   l->AddEntry(h2, "D", "P");
   l->AddEntry(h3, "E", "P");
   l->AddEntry(h4, "F", "P");
   l->Draw();
   
   c->SaveAs("plots/sf_dataset_0p3.pdf"); 
   c->SaveAs("plots/sf_dataset_0p3.png");
}

void sf_systematics()
{
   sf_systematics_toppt();
   sf_systematics_pu();
   sf_systematics_dataset();
}


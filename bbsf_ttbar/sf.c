#include <TFile.h>
#include <TExec.h>
#include <TPaveText.h>
#include <TLatex.h>
#include "addOverflow.h"
#include <TH1D.h>
#include <TChain.h>
#include <TCanvas.h>
#include "labelCMS.h"

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

TH1D* makesf(TH1D *mistag_data, TH1D *mistag_mc, const TString name, const TH1D * h)
{
   TH1D* h_sf = (TH1D*)h->Clone(name);
   h_sf->SetTitle(";AK8 jet p_{T} [GeV];data / mc");
   for (int i = 1; i <= h_sf->GetNbinsX(); ++i) {
      const double x = h_sf->GetBinCenter(i);
      const double y_data = mistag_data->GetBinContent(i);
      const double yerr_data = mistag_data->GetBinError(i);      
      const double y_mc = mistag_mc->GetBinContent(i);
      const double yerr_mc = mistag_mc->GetBinError(i);
      const double sf = y_data/y_mc;
      h_sf->SetBinContent(i, sf);
      double err = sqrt(pow(yerr_data/y_data, 2) + pow(yerr_mc/y_mc, 2));
      err = sf * err;
      h_sf->SetBinError(i, err);
   }
   return h_sf;
}

void runPoint(const double wp, const TString wptag, TH1D * h)
{
   // pu weights
   TString puweight, puweighttag;
   puweight = "puWeight"; puweighttag = puweight;
   //puweight = "1.";  puweighttag = "noPuWeight";
   //puweight = "puWeight_up"; puweighttag = puweight;
   //puweight = "puWeight_down"; puweighttag = puweight;

   // top pt weights (ttbar only)
   TString topweight, topweighttag;
   //topweight = "ttbar_ptweight"; topweighttag = "topPtWeight";
   topweight = "1."; topweighttag = "noTopPtWeight";

   // dataset
   TString datatag = "";
   TString datatitle = "";
   double lumi = 0.;
   //datatag = "SingleMuon"; lumi = 41860.;  datatitle = "Single-Muon";
   //datatag = "SingleMuonB"; lumi = 4792.; datatitle = "Single-Muon B";
   //datatag = "SingleMuonC"; lumi = 9755.; datatitle = "Single-Muon C";
   //datatag = "SingleMuonD"; lumi = 4319.; datatitle = "Single-Muon D";
   //datatag = "SingleMuonE"; lumi = 9424.; datatitle = "Single-Muon E";
   datatag = "SingleMuonF"; lumi = 13570.; datatitle = "Single-Muon F";

   //////
   // fill data hists
   //////
 
   char filebuffer[10000];
   sprintf(filebuffer, "./sf_rootfiles/sf_%s_%s_%s_%s.root", wptag.Data(), datatag.Data(), puweighttag.Data(), topweighttag.Data());
   TFile * f = new TFile(filebuffer, "RECREATE");
 
   TChain * chain4_data = new TChain("treeAK4");
   TChain * chain8_data = new TChain("treeAK8");
   chain4_data->Add("skims/"+datatag+".root");
   chain8_data->Add("skims/"+datatag+".root");
   chain4_data->AddFriend(chain8_data);

   TH1D * h_denom_data = (TH1D*)h->Clone("h_denom_data");
   chain4_data->Project("h_denom_data", "FatJetInfo.Jet_pt[idx_ak8]");
   addOverflow(h_denom_data);
   
   TH1D * h_num_data = (TH1D*)h->Clone("h_num_data");
   char buffer_data[1000];
   sprintf(buffer_data, "FatJetInfo.Jet_DoubleSV[idx_ak8]>=%f", wp);
   chain4_data->Project("h_num_data", "FatJetInfo.Jet_pt[idx_ak8]", buffer_data);
   addOverflow(h_num_data);
   
   //////
   // fill signal hists
   //////

   TChain * chain4_sig[3], *chain8_sig[3];
   for (int i = 0; i < 3; ++i) {
      chain4_sig[i] = new TChain("treeAK4");
      chain8_sig[i] = new TChain("treeAK8");
   }
   chain4_sig[0]->Add("./skims/TTToSemiLeptonic.root"); chain8_sig[0]->Add("./skims/TTToSemiLeptonic.root");
   chain4_sig[1]->Add("./skims/TTToHadronic.root");     chain8_sig[1]->Add("./skims/TTToHadronic.root");
   chain4_sig[2]->Add("./skims/TTTo2L2Nu.root");        chain8_sig[2]->Add("./skims/TTTo2L2Nu.root");
   for (int i = 0; i < 3; ++i) chain4_sig[i]->AddFriend(chain8_sig[i]);
   
   TH1D * h_num_sig = (TH1D*)h->Clone("h_num_sig");
   TH1D * h_denom_sig = (TH1D*)h->Clone("h_denom_sig");
   for (int i = 0; i < 3; ++i) {
      char buffer_denom[1000], buffer_num[1000];
      sprintf(buffer_denom, "%f * xsWeight * %s * %s * (1>0)", lumi, puweight.Data(), topweight.Data());
      sprintf(buffer_num, "%f * xsWeight * %s * %s * (FatJetInfo.Jet_DoubleSV[idx_ak8]>=%f)", lumi, puweight.Data(), topweight.Data(), wp);
      chain4_sig[i]->Project("+h_denom_sig", "FatJetInfo.Jet_pt[idx_ak8]", buffer_denom);
      chain4_sig[i]->Project("+h_num_sig", "FatJetInfo.Jet_pt[idx_ak8]", buffer_num);
   }
   addOverflow(h_denom_sig);
   addOverflow(h_num_sig);

   //////
   // fill mc background hists
   //////

   TChain * chain4_bkg[5], *chain8_bkg[5];
   double scale_bkg[5];
   for (int i = 0; i < 5; ++i) {
      chain4_bkg[i] = new TChain("treeAK4");
      chain8_bkg[i] = new TChain("treeAK8");
   } 
   chain4_bkg[0]->Add("skims/W3JetsToLNu.root");       scale_bkg[0] = 1.; chain8_bkg[0]->Add("skims/W3JetsToLNu.root");
   chain4_bkg[1]->Add("skims/W4JetsToLNu.root");       scale_bkg[1] = 1.; chain8_bkg[1]->Add("skims/W4JetsToLNu.root");
   chain4_bkg[2]->Add("skims/ST_tW_antitop_5f.root");  scale_bkg[2] = 1.; chain8_bkg[2]->Add("skims/ST_tW_antitop_5f.root");
   chain4_bkg[3]->Add("skims/ST_tW_top_5f.root");      scale_bkg[3] = 1.; chain8_bkg[3]->Add("skims/ST_tW_top_5f.root");
   chain4_bkg[4]->Add("skims/WWToLNuQQ_NNPDF31.root"); scale_bkg[4] = 1.; chain8_bkg[4]->Add("skims/WWToLNuQQ_NNPDF31.root");
   for (int i = 0; i < 5; ++i) chain4_bkg[i]->AddFriend(chain8_bkg[i]);

   TH1D * h_num_bkg = (TH1D*)h->Clone("h_num_bkg");
   TH1D * h_denom_bkg = (TH1D*)h->Clone("h_denom_bkg");
   for (int i = 0; i < 5; ++i) {
      char buffer_denom[1000], buffer_num[1000];
      sprintf(buffer_denom, "%f * %f * xsWeight * %s * (1>0)", scale_bkg[i], lumi, puweight.Data());
      sprintf(buffer_num, "%f * %f * xsWeight * %s * (FatJetInfo.Jet_DoubleSV[idx_ak8]>=%f)", scale_bkg[i], lumi, puweight.Data(), wp);
      chain4_bkg[i]->Project("+h_denom_bkg", "FatJetInfo.Jet_pt[idx_ak8]", buffer_denom);
      chain4_bkg[i]->Project("+h_num_bkg", "FatJetInfo.Jet_pt[idx_ak8]", buffer_num);
   }
   addOverflow(h_denom_bkg);
   addOverflow(h_num_bkg);

   //////
   // make mistag rates and sfs
   //////

   // subtract background from data
   TH1D * h_num_data_bkgsub = (TH1D*)h_num_data->Clone("h_num_data_bkgsub");
   h_num_data_bkgsub->Add(h_num_bkg, -1.);
   TH1D * h_denom_data_bkgsub = (TH1D*)h_denom_data->Clone("h_denom_data_bkgsub");
   h_denom_data_bkgsub->Add(h_denom_bkg, -1.);
   // make mistag rate for background subtracted data
   TH1D * h_mistag_data_bkgsub = (TH1D*)h_num_data_bkgsub->Clone("h_mistag_data_bkgsub");
   h_mistag_data_bkgsub->Divide(h_denom_data_bkgsub);
   h_mistag_data_bkgsub->SetTitle(";AK8 jet p_{T} [GeV];mistag rate");

   // make mistag rate for signal
   TH1D * h_mistag_sig = (TH1D*)h_num_sig->Clone("h_mistag_sig");
   h_mistag_sig->Divide(h_denom_sig);
   h_mistag_sig->SetTitle(";AK8 jet p_{T} [GeV];mistag rate");

   // make sf
   TH1D * h_sf = makesf(h_mistag_data_bkgsub, h_mistag_sig, "h_sf", h);

   f->Write();
   f->Close();
}

void sf()
{
   double x_3[4] = {250., 350., 430., 700.};
   TH1D * h_3 = new TH1D("h_3", ";AK8 jet p_{T} [GeV];events / bin", 3, x_3);
   runPoint(0.3, "0p3", h_3);

   TH1D * h_6 = new TH1D("h_6", ";AK8 jet p_{T} [GeV];events / bin", 3, x_3);
   runPoint(0.6, "0p6", h_6);

   TH1D * h_8 = new TH1D("h_8", ";AK8 jet p_{T} [GeV];events / bin", 3, x_3);
   runPoint(0.8, "0p8", h_8);

   //double x_9[3] = {250., 350., 700.};
   TH1D * h_9 = new TH1D("h_9", ";AK8 jet p_{T} [GeV];events / bin", 3, x_3);
   runPoint(0.9, "0p9", h_9);
}

void sf_plot(const TString filetag)
{
   TFile * f = TFile::Open("./sf_rootfiles/" + filetag + ".root");
   TH1D * h_sf = (TH1D*)f->Get("h_sf");
   TH1D * h_mistag_sig = (TH1D*)f->Get("h_mistag_sig");
   TH1D * h_mistag_data_bkgsub = (TH1D*)f->Get("h_mistag_data_bkgsub");
   TH1D * h_num_data_bkgsub = (TH1D*)f->Get("h_num_data_bkgsub");
   TH1D * h_denom_data_bkgsub = (TH1D*)f->Get("h_denom_data_bkgsub");
   TH1D * h_num_sig = (TH1D*)f->Get("h_num_sig");
   TH1D * h_denom_sig = (TH1D*)f->Get("h_denom_sig");

   TExec *drawtext_mistag_sig = new TExec("drawtext_mistag_sig", "drawtext(\"h_mistag_sig\");");
   h_mistag_sig->GetListOfFunctions()->Add(drawtext_mistag_sig);

   TExec *drawtext_mistag_data_bkgsub = new TExec("drawtext_mistag_databkgsub", "drawtext(\"h_mistag_data_bkgsub\");");
   h_mistag_data_bkgsub->GetListOfFunctions()->Add(drawtext_mistag_data_bkgsub);

   TExec *drawtext_sf = new TExec("drawtext_sf", "drawtext(\"h_sf\");");
   h_sf->GetListOfFunctions()->Add(drawtext_sf);

   // calculate inclusive mistag data
   const double datanumint = h_num_data_bkgsub->Integral();
   const double datadenomint = h_denom_data_bkgsub->Integral();
   const double data_eff = datanumint/datadenomint;
   const double data_efferr = sqrt(datanumint*(1.-data_eff))/datadenomint;
   TPaveText *pt_1 = new TPaveText(.5, .8, .86, .875, "NDC");
   pt_1->SetBorderSize(0);
   pt_1->SetFillColor(0);
   pt_1->SetTextAlign(33);
   char buffer_1[1000];
   sprintf(buffer_1, "inclusive: %.3f #pm %.3f", data_eff, data_efferr);
   pt_1->AddText(buffer_1);
   
   // calculate inclusive mistag mc
   const double mcnumint = h_num_sig->Integral();
   const double mcdenomint = h_denom_sig->Integral();
   const double mc_eff = mcnumint/mcdenomint;
   const double mc_efferr = sqrt(mcnumint*(1.-mc_eff))/mcdenomint;
   TPaveText *pt_2 = new TPaveText(.5, .8, .86, .875, "NDC");
   pt_2->SetBorderSize(0);
   pt_2->SetFillColor(0);
   pt_2->SetTextAlign(33);
   char buffer_2[1000];
   sprintf(buffer_2, "inclusive: %.3f #pm %.3f", mc_eff, mc_efferr);
   pt_2->AddText(buffer_2);
 
   // calculate inclusive sf
   const double inceff = data_eff / mc_eff;
   double inceff_err = sqrt(pow(data_efferr/data_eff, 2) + pow(mc_efferr/mc_eff, 2));
   inceff_err = inceff * inceff_err;
   TPaveText *pt_3 = new TPaveText(.5, .8, .86, .875, "NDC");
   pt_3->SetBorderSize(0);
   pt_3->SetFillColor(0);
   pt_3->SetTextAlign(33);
   char buffer_3[1000];
   sprintf(buffer_3, "inclusive: %.3f #pm %.3f", inceff, inceff_err);
   pt_3->AddText(buffer_3);

   TCanvas * canvas = new TCanvas("canvas", "", 1200, 400);
   canvas->Divide(3, 1);

   canvas->cd(1);
   h_mistag_data_bkgsub->SetMarkerStyle(20);
   h_mistag_data_bkgsub->SetLineColor(1);
   h_mistag_data_bkgsub->SetStats(0);
   h_mistag_data_bkgsub->Draw("PE");
   h_mistag_data_bkgsub->SetMinimum(0.);
   h_mistag_data_bkgsub->SetMaximum(0.8);
   pt_1->Draw();

   canvas->cd(2);
   h_mistag_sig->SetMarkerStyle(20);
   h_mistag_sig->SetLineColor(1);
   h_mistag_sig->SetStats(0);
   h_mistag_sig->Draw("PE");
   h_mistag_sig->SetMinimum(0.);
   h_mistag_sig->SetMaximum(0.8);
   pt_2->Draw();
   labelSim();

   canvas->cd(3);
   h_sf->SetMarkerStyle(20);
   h_sf->SetLineColor(1);
   h_sf->SetStats(0);
   h_sf->Draw("PE");
   h_sf->SetStats(0);
   h_sf->SetMinimum(0.6);
   h_sf->SetMaximum(1.4);
   pt_3->Draw();

   canvas->SaveAs("./plots/" + filetag + ".pdf");
   canvas->SaveAs("./plots/" + filetag + ".png");
}



#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include<map>
#include<iostream>

// InclName=
// systName=
// sample = "TCL"
void SemileptonicSL(TString InclName, TString SemiName, TString SystName, TString sample) {

  gROOT->SetStyle("CMS");

  TString path = "~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/";
	
	std::map<TString, TCanvas*> cv_map;
	
	TH1D *taggedjet_pt_b_I;
	TH1D *taggedjet_pt_c_I;
	TH1D *jet_pt_b_I;
	TH1D *jet_pt_c_I;
	TH1D *eff_pt_b_I;
	
	TH1D *taggedjet_pt_b_S;
	TH1D *taggedjet_pt_c_S;
	TH1D *jet_pt_b_S;
	TH1D *jet_pt_c_S;
	TH1D *eff_pt_b_S;

	TH1D *taggedjet_pt_b_syst;
	TH1D *taggedjet_pt_c_syst;
	TH1D *jet_pt_b_syst;
	TH1D *jet_pt_c_syst;
	TH1D *eff_pt_b_syst;
	
	TH1D *scale_pt;

	TH1D *taggedjet_eta_b_I;
	TH1D *taggedjet_eta_c_I;
	TH1D *jet_eta_b_I;
	TH1D *jet_eta_c_I;
	TH1D *eff_eta_b_I;
	
	TH1D *taggedjet_eta_b_S;
	TH1D *taggedjet_eta_c_S;
	TH1D *jet_eta_b_S;
	TH1D *jet_eta_c_S;
	TH1D *eff_eta_b_S;

	TH1D *taggedjet_eta_b_syst;
	TH1D *taggedjet_eta_c_syst;
	TH1D *jet_eta_b_syst;
	TH1D *jet_eta_c_syst;
	TH1D *eff_eta_b_syst;
		
	TH1D *scale_eta;
	
	TFile *inclfile = new TFile(path+InclName);
	TFile *semifile = new TFile(path+SemiName);
	TFile *systfile = new TFile(path+SystName);
	
	//__________________________________________________________

	inclfile->cd();
	
	taggedjet_pt_b_I = (TH1D*) gDirectory->Get("taggedjet_pt_b");
	taggedjet_pt_c_I = (TH1D*) gDirectory->Get("taggedjet_pt_c");
	jet_pt_b_I = (TH1D*) gDirectory->Get("jet_pt_b");
	jet_pt_c_I = (TH1D*) gDirectory->Get("jet_pt_c");
	//std::cout << "hola1" << std::endl;
	
	
	taggedjet_eta_b_I = (TH1D*) gDirectory->Get("taggedjet_eta_b");
	taggedjet_eta_c_I = (TH1D*) gDirectory->Get("taggedjet_eta_c");
	jet_eta_b_I = (TH1D*) gDirectory->Get("jet_eta_b");
	jet_eta_c_I = (TH1D*) gDirectory->Get("jet_eta_c");
	
	semifile->cd();
	
	taggedjet_pt_b_S = (TH1D*) gDirectory->Get("taggedjet_pt_b");
	taggedjet_pt_c_S = (TH1D*) gDirectory->Get("taggedjet_pt_c");
	jet_pt_b_S = (TH1D*) gDirectory->Get("jet_pt_b");
	jet_pt_c_S = (TH1D*) gDirectory->Get("jet_pt_c");
	//std::cout << "hola2" << std::endl;
	taggedjet_eta_b_S = (TH1D*) gDirectory->Get("taggedjet_eta_b");
	taggedjet_eta_c_S = (TH1D*) gDirectory->Get("taggedjet_eta_c");
	jet_eta_b_S = (TH1D*) gDirectory->Get("jet_eta_b");
	jet_eta_c_S = (TH1D*) gDirectory->Get("jet_eta_c");

	systfile->cd();
	
	taggedjet_pt_b_syst = (TH1D*) gDirectory->Get("taggedjet_pt_b");
	taggedjet_pt_c_syst = (TH1D*) gDirectory->Get("taggedjet_pt_c");
	jet_pt_b_syst = (TH1D*) gDirectory->Get("jet_pt_b");
	jet_pt_c_syst = (TH1D*) gDirectory->Get("jet_pt_c");
	taggedjet_eta_b_syst = (TH1D*) gDirectory->Get("taggedjet_eta_b");
	taggedjet_eta_c_syst = (TH1D*) gDirectory->Get("taggedjet_eta_c");
	jet_eta_b_syst = (TH1D*) gDirectory->Get("jet_eta_b");
	jet_eta_c_syst = (TH1D*) gDirectory->Get("jet_eta_c");

	
	eff_pt_b_I = (TH1D*) jet_pt_b_S->Clone("eff_pt_b_I");
	eff_pt_b_S = (TH1D*) jet_pt_b_S->Clone("eff_pt_b_S");
	eff_pt_b_syst = (TH1D*) jet_pt_b_S->Clone("eff_pt_b_syst");
	scale_pt = (TH1D*) jet_pt_b_S->Clone("scale_pt");
	eff_pt_b_I->Reset();
	eff_pt_b_S->Reset();
	eff_pt_b_syst->Reset();
	scale_pt->Reset();
	eff_pt_b_I->Sumw2();
	eff_pt_b_S->Sumw2();
	eff_pt_b_syst->Sumw2();
	scale_pt->Sumw2();
		
	eff_eta_b_I = (TH1D*) jet_eta_b_S->Clone("eff_eta_b_I");
	eff_eta_b_S = (TH1D*) jet_eta_b_S->Clone("eff_eta_b_S");
	eff_eta_b_syst = (TH1D*) jet_eta_b_S->Clone("eff_eta_b_syst"); 
	scale_eta = (TH1D*) jet_eta_b_S->Clone("scale_eta");
	eff_eta_b_I->Reset();
	eff_eta_b_S->Reset();
	eff_eta_b_syst->Reset();
	scale_eta->Reset();
	eff_eta_b_I->Sumw2();
	eff_eta_b_S->Sumw2();
	eff_eta_b_syst->Sumw2();
	scale_eta->Sumw2();

	//std::cout << "hola3" << std::endl;
	eff_pt_b_I->Divide(taggedjet_pt_b_I,jet_pt_b_I,1.,1.,"B");
	eff_pt_b_S->Divide(taggedjet_pt_b_S,jet_pt_b_S,1.,1.,"B");
	scale_pt->Divide(eff_pt_b_I,eff_pt_b_S);
	eff_pt_b_syst->Divide(taggedjet_pt_b_syst,jet_pt_b_syst,1.,1.,"B");
	
	//TH1D * delta_pt = (TH1D*) jet_pt_b_S->Clone("delta_pt");
	//delta_pt->Reset();
	//delta_pt->Sumw2();
	double minsyst = 999999.;
	double maxsyst = -1.;
	for (int i=1; i<=eff_pt_b_I->GetNbinsX(); ++i) {
	  double delta_pt = TMath::Abs(eff_pt_b_I->GetBinContent(i) - eff_pt_b_syst->GetBinContent(i))/(eff_pt_b_S->GetBinContent(i));
	  double error = scale_pt->GetBinError(i);
	  std::cout << " pt bin # " << i << " stat= " << error << " syst= " << delta_pt << std::endl;
	  scale_pt->SetBinError(i, TMath::Sqrt(error*error+delta_pt*delta_pt));
	  if (minsyst > delta_pt && i>1) minsyst = delta_pt;
	  if (maxsyst < delta_pt && i>1) maxsyst = delta_pt;
	}
	std::cout << " pt min syst = " << minsyst << std::endl;
	std::cout << " pt max syst = " << maxsyst << std::endl;
	
	eff_eta_b_I->Divide(taggedjet_eta_b_I,jet_eta_b_I,1.,1.,"B");
	eff_eta_b_S->Divide(taggedjet_eta_b_S,jet_eta_b_S,1.,1.,"B");
	scale_eta->Divide(eff_eta_b_I,eff_eta_b_S);
	eff_eta_b_syst->Divide(taggedjet_eta_b_syst,jet_eta_b_syst,1.,1.,"B");
	minsyst = 999999.;
	maxsyst = -1.;
	for (int i=1; i<=eff_eta_b_I->GetNbinsX(); ++i) { 
          double delta_pt = TMath::Abs(eff_eta_b_I->GetBinContent(i) - eff_eta_b_syst->GetBinContent(i))/(eff_eta_b_S->GetBinContent(i)); 
          double error = scale_eta->GetBinError(i); 
		  std::cout << "eta bin # " << i << " stat= " << error << " syst= " << delta_pt << std::endl; 
          scale_eta->SetBinError(i, TMath::Sqrt(error*error+delta_pt*delta_pt));
		  if (minsyst > delta_pt) minsyst = delta_pt;
		  if (maxsyst < delta_pt) maxsyst = delta_pt;
	} 
	std::cout << " eta min syst = " << minsyst << std::endl;
	std::cout << " eta max syst = " << maxsyst << std::endl;
	
	
	//std::cout << "hola4" << std::endl;
	
	cv_map["SF_pt_"+sample] = new TCanvas("SF_pt_"+sample,"SF_pt_"+sample,700,700);

	eff_pt_b_I->SetXTitle("jet P_{T} [GeV/c]");
	
	eff_pt_b_I->SetMaximum(1.2);
	eff_pt_b_I->SetMinimum(0.1);
	
	eff_pt_b_I->SetMarkerStyle(26);
	eff_pt_b_S->SetMarkerStyle(23);
	eff_pt_b_I->SetMarkerColor(2);
	eff_pt_b_S->SetMarkerColor(4);
	eff_pt_b_I->SetLineColor(2);
	eff_pt_b_S->SetLineColor(4);
	scale_pt->SetLineColor(1);
	
	eff_pt_b_I->Draw("PE1");
	eff_pt_b_S->Draw("PE1 same");
	scale_pt->Draw("PE1 same");
	scale_pt->Fit("pol0","0");
	//TF1 *fscale = scale->GetFunction("pol0");
	//fscale->SetLineColor(3);
	//fscale->Draw("same");
	
	TLegend *legpt = new TLegend(0.65,0.35,0.85,0.55,sample,"NDC");
	legpt->SetFillColor(10);
	legpt->AddEntry(eff_pt_b_I,"#epsilon_{b} incl","P");
	legpt->AddEntry(eff_pt_b_S,"#epsilon_{b#rightarrow#mu}","P");
	legpt->AddEntry(scale_pt,"SF","P");
	legpt->Draw();

	cv_map["SF_pt_"+sample]->Print("SF_pt_"+sample+".eps");

	
	cv_map["SF_eta_"+sample] = new TCanvas("SF_eta_"+sample,"SF_eta_"+sample,700,700);

	eff_eta_b_I->SetXTitle("jet |#eta|");
	
	eff_eta_b_I->SetMaximum(1.2);
	eff_eta_b_I->SetMinimum(0.1);
	
	eff_eta_b_I->SetMarkerStyle(26);
	eff_eta_b_S->SetMarkerStyle(23);
	eff_eta_b_I->SetMarkerColor(2);
	eff_eta_b_S->SetMarkerColor(4);
	eff_eta_b_I->SetLineColor(2);
	eff_eta_b_S->SetLineColor(4);
	scale_eta->SetLineColor(1);

	eff_eta_b_I->Draw("PE1");
	eff_eta_b_S->Draw("PE1 same");
	scale_eta->Draw("PE1 same");
	scale_eta->Fit("pol0","0");
	//TF1 *fscale = scale->GetFunction("pol0");
	//fscale->SetLineColor(3);
	//fscale->Draw("same");
	
	TLegend *legeta = new TLegend(0.65,0.35,0.85,0.55,sample,"NDC");
	legeta->SetFillColor(10);
	legeta->AddEntry(eff_eta_b_I,"#epsilon_{b} incl","P");
	legeta->AddEntry(eff_eta_b_S,"#epsilon_{b#rightarrow#mu}","P");
	legeta->AddEntry(scale_eta,"SF","P");
	//legeta->Draw();

	cv_map["SF_eta_"+sample]->Print("SF_eta_"+sample+".eps");
	

}

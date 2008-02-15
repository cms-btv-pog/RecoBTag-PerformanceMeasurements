
#include "S8Solver.h"
#include "S8AnalyticSolver.h"
#include "S8NumericSolver.h"
#include "S8FitSolver.h"

#include "TF1.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TArrayD.h"
#include "TLegend.h"

#include<iostream>

ClassImp(S8Solver)

//____________________________________________________________
S8Solver::S8Solver() {
	
//}
//____________________________________________________________
//S8Solver::S8Solver(std::string name) {

	fAlphaConst = true;
	fBetaConst = false;
	fcategory = "pT";
	fminPtrel = 0.8;
	fMaxPtrel = 3.0;//-1;
	fVerbose = false;
	fmethod = "numeric";
	frebin = false;
	frecalculateFactors = true;
	fthename = "";
	fAlphaf = 1.;
	fBetaf = 1.;
	fKappabf = 1.;
	fKappaclf = 1.;
	fKappabConst = false;
	fKappaclConst = true;
	fisCorrFile = false;
	fDeltaConst = false;
	fGammaConst = false;

}
//____________________________________________________________
void S8Solver::LoadHistos() {

	//if (finputCorrFile) {
	//	fisCorrFile = true;
	//	std::cout << "use different sample for correlation factors" << std::endl;
	//}
	
	finputFile->cd();
	
	
	fnHistoBase = (TH2F*) gDirectory->Get("Histograms/samples/n_"+fcategory);
	fpHistoBase = (TH2F*) gDirectory->Get("Histograms/samples/p_"+fcategory);
	fnSvxHistoBase = (TH2F*) gDirectory->Get("Histograms/samples/ntag_"+fcategory+"_"+fthename);
	fpSvxHistoBase = (TH2F*) gDirectory->Get("Histograms/samples/ptag_"+fcategory+"_"+fthename);

	//std::cout << " 2D histos loaded " << fcategory << std::endl;

	// get bin for pTrel cut
	std::cout << "entries fn =" << fnHistoBase->GetEntries() << std::endl;
	int ith_ptrel_bin = (int) fnHistoBase->GetYaxis()->FindBin(fminPtrel);
	int ith_max_bin = -1;
	if (fMaxPtrel != -1) ith_max_bin = (int) fnHistoBase->GetYaxis()->FindBin(fMaxPtrel);
	std::cout << " ptrel bin " << ith_max_bin << std::endl;
	
	fnHisto = (TH1D*) fnHistoBase->ProjectionX("fnHisto",-1,ith_max_bin,"e");
	fpHisto = (TH1D*) fpHistoBase->ProjectionX("fpHisto",-1,ith_max_bin,"e");
	fnHistoMu = (TH1D*) fnHistoBase->ProjectionX("fnHistoMu",ith_ptrel_bin,ith_max_bin,"e");
	fpHistoMu = (TH1D*) fpHistoBase->ProjectionX("fpHistoMu",ith_ptrel_bin,ith_max_bin,"e");
	fnHistoSvx = (TH1D*) fnSvxHistoBase->ProjectionX("fnHistoSvx", -1, ith_max_bin, "e");
	fpHistoSvx = (TH1D*) fpSvxHistoBase->ProjectionX("fpHistoSvx", -1, ith_max_bin, "e");
	fnHistoAll = (TH1D*) fnSvxHistoBase->ProjectionX("fnHistoAll",ith_ptrel_bin, ith_max_bin,"e");
	fpHistoAll = (TH1D*) fpSvxHistoBase->ProjectionX("fpHistoAll",ith_ptrel_bin, ith_max_bin,"e");
	
	std::cout << " got projections" << std::endl;
	// rebin correlation factors
	const int ncorrptarray = 6;
	const int ncorretaarray = 3;
	Double_t corrptbins[ncorrptarray] = {30.,40.,60.,80,120,230.};//{30.,40.,60.,80,230.};
	Double_t corretabins[ncorrptarray] = {0.,1.5,2.5};
	//if (fcategory=="Eta") {
	//  ncorrptarray = tmpncorretaarray;
	  
	// recalculate correlation factors
	if (frecalculateFactors) {
		std::cout << "recalculate " << std::endl;
		fh_alpha = (TH1D*) fnHisto->Clone("fh_alpha");
		fh_beta = (TH1D*) fnHisto->Clone("fh_beta");
		fh_kcl = (TH1D*) fnHisto->Clone("fh_kcl");
		fh_kb = (TH1D*) fnHisto->Clone("fh_kb");
		fh_delta = (TH1D*) fnHisto->Clone("fh_delta");
		fh_gamma = (TH1D*) fnHisto->Clone("fh_gamma");
		fh_alpha->Reset();
		fh_beta->Reset();
		fh_kcl->Reset();
		fh_kb->Reset();
		fh_delta->Reset();
		fh_gamma->Reset();
		//fh_alpha->Sumw2();
		//fh_beta->Sumw2();
		//fh_kcl->Sumw2();
		//fh_kb->Sumw2();
		
		std::cout << "reset1" << std::endl;
		
		std::map< TString, TH1*> h1;
		std::map< TString, TH2*> h2;
		
		if (fisCorrFile) {
			std::cout << "what" << std::endl;
			finputCorrFile->cd();
			std::cout << " another file" << std::endl;
		}
		finputFile->cd();
		
		h2["b_npT"] = (TH2F*) gDirectory->Get("Histograms/samples/n_"+fcategory+"_b");
		std::cout << "got one" << std::endl;
		h2["cl_npT"] = (TH2F*) gDirectory->Get("Histograms/samples/n_"+fcategory+"_cl");
		h2["b_ppT"] = (TH2F*) gDirectory->Get("Histograms/samples/p_"+fcategory+"_b");
		h2["cl_ppT"] = (TH2F*) gDirectory->Get("Histograms/samples/p_"+fcategory+"_cl");
		h2["b_ncmbpT"] = (TH2F*) gDirectory->Get("Histograms/samples/ntag_"+fcategory+"_b_"+fthename);
		h2["cl_ncmbpT"] = (TH2F*) gDirectory->Get("Histograms/samples/ntag_"+fcategory+"_cl_"+fthename);
		h2["b_pcmbpT"] = (TH2F*) gDirectory->Get("Histograms/samples/ptag_"+fcategory+"_b_"+fthename);
		h2["cl_pcmbpT"] = (TH2F*) gDirectory->Get("Histograms/samples/ptag_"+fcategory+"_cl_"+fthename);

		std::cout << " got initial truth dist." << std::endl;
		
		TH1D *halljets_b = h2["b_npT"]->ProjectionX("halljets_b", -1 , ith_max_bin,"e");
		TH1D *halljets_cl = h2["cl_npT"]->ProjectionX("halljets_cl", -1 , ith_max_bin,"e");
		TH1D *halloppjets_b = h2["b_ppT"]->ProjectionX("halloppjets_b", -1 , ith_max_bin,"e");
		TH1D *halloppjets_cl = h2["cl_ppT"]->ProjectionX("halloppjets_cl", -1 , ith_max_bin,"e");
		
		h1["eff_pTrel_b"] = (TH1D*) fnHisto->Clone("eff_pTrel_b");
		h1["eff_pTrel_cl"] = (TH1D*) fnHisto->Clone("eff_pTrel_cl");
		h1["eff_pTrel_TaggedJet_b"] = (TH1D*) fnHisto->Clone("eff_pTrel_TaggedJet_b");
		h1["eff_pTrel_TaggedJet_cl"] = (TH1D*) fnHisto->Clone("eff_pTrel_TaggedJet_cl");
		h1["eff_TaggedJet_b"] = (TH1D*) fnHisto->Clone("eff_TaggedJet_b");
		h1["eff_TaggedJet_cl"] = (TH1D*) fnHisto->Clone("eff_TaggedJet_cl");
		h1["eff_TaggedBothJets_b"] = (TH1D*) fnHisto->Clone("eff_TaggedBothJets_b");
		h1["eff_TaggedBothJets_cl"] = (TH1D*) fnHisto->Clone("eff_TaggedBothJets_cl");

		h1["eff_mu_taggedaway_b"] = (TH1D*) fnHisto->Clone("eff_mu_taggedaway_b");
		h1["eff_mu_taggedaway_cl"] = (TH1D*) fnHisto->Clone("eff_mu_taggedaway_cl");

		h1["eff_pTrel_b"]->Reset();
		h1["eff_pTrel_cl"]->Reset();
		h1["eff_pTrel_TaggedJet_b"]->Reset();
		h1["eff_pTrel_TaggedJet_cl"]->Reset(); 
		h1["eff_TaggedJet_b"]->Reset(); 
		h1["eff_TaggedJet_cl"]->Reset();
		h1["eff_TaggedBothJets_b"]->Reset();
		h1["eff_TaggedBothJets_cl"]->Reset();
		h1["eff_mu_taggedaway_b"]->Reset();
		h1["eff_mu_taggedaway_cl"]->Reset();

		//h1["eff_pTrel_b"]->Sumw2();
		//h1["eff_pTrel_cl"]->Sumw2();
		//h1["eff_pTrel_TaggedJet_b"]->Sumw2();
		//h1["eff_pTrel_TaggedJet_cl"]->Sumw2(); 
		//h1["eff_TaggedJet_b"]->Sumw2(); 
		//h1["eff_TaggedJet_cl"]->Sumw2();  
		//h1["eff_TaggedBothJets_b"]->Sumw2();
		//h1["eff_TaggedBothJets_cl"]->Sumw2();
		std::cout << "reset" << std::endl;
		
		h1["eff_pTrel_b"]->Divide( h2["b_npT"]->ProjectionX("b_halljets_ptrel",ith_ptrel_bin, ith_max_bin ,"e") , halljets_b ,1.,1.,"B");
		h1["eff_pTrel_cl"]->Divide( h2["cl_npT"]->ProjectionX("cl_halljets_ptrel", ith_ptrel_bin, ith_max_bin,"e") , halljets_cl ,1.,1.,"B");
		h1["eff_TaggedJet_b"]->Divide( h2["b_ncmbpT"]->ProjectionX("b_halljets_tagged", -1 , ith_max_bin,"e") , halljets_b ,1.,1.,"B");
		h1["eff_TaggedJet_cl"]->Divide( h2["cl_ncmbpT"]->ProjectionX("cl_halljets_tagged", -1 , ith_max_bin,"e") , halljets_cl ,1.,1.,"B");
		h1["eff_pTrel_TaggedJet_b"]->Divide( h2["b_ncmbpT"]->ProjectionX("b_halljets_ptreltagged",ith_ptrel_bin, ith_max_bin ,"e") , halljets_b ,1.,1.,"B");
		h1["eff_pTrel_TaggedJet_cl"]->Divide( h2["cl_ncmbpT"]->ProjectionX("cl_halljets_ptrelltagged", ith_ptrel_bin, ith_max_bin,"e") , halljets_cl ,1.,1.,"B");
		
		h1["eff_TaggedBothJets_b"]->Divide( h2["b_pcmbpT"]->ProjectionX("b_halloppjets_tagged", -1 , ith_max_bin,"e") , halloppjets_b ,1.,1.,"B");
		h1["eff_TaggedBothJets_cl"]->Divide( h2["cl_pcmbpT"]->ProjectionX("cl_halloppjets_tagged", -1 , ith_max_bin,"e") , halloppjets_cl ,1.,1.,"B");

		h1["eff_mu_taggedaway_b"]->Divide( h2["b_ppT"]->ProjectionX("b_hallopppjets_ptrel",ith_ptrel_bin, ith_max_bin ,"e") , halloppjets_b, 1.,1.,"B");
		h1["eff_mu_taggedaway_cl"]->Divide( h2["cl_ppT"]->ProjectionX("cl_halloppjets_ptrel",ith_ptrel_bin, ith_max_bin ,"e") , halloppjets_cl, 1.,1.,"B");

		feffTag_b = (TH1D*) h1["eff_TaggedJet_b"]->Clone("feffTag_b"); 
		feffTag_cl = (TH1D*) h1["eff_TaggedJet_cl"]->Clone("feffTag_cl");
		std::cout << "clonned" << std::endl;
		
		if (frebin) {
			std::cout << "rebinning" << std::endl;
			TH1D* tmpfh_alpha = (TH1D*) fnHisto->Clone("tmpfh_alpha");
			TH1D* tmpfh_beta = (TH1D*) fnHisto->Clone("tmpfh_beta");
			TH1D* tmpfh_kb = (TH1D*)  fnHisto->Clone("tmpfh_kb");
			TH1D* tmpfh_kcl = (TH1D*) fnHisto->Clone("tmpfh_kcl");
			TH1D* tmpfh_gamma = (TH1D*) fnHisto->Clone("tmpfh_gamma");
			TH1D* tmpfh_delta = (TH1D*) fnHisto->Clone("tmpfh_delta");
			tmpfh_alpha->Reset();
			tmpfh_beta->Reset();
			tmpfh_kb->Reset();
			tmpfh_kcl->Reset();
			tmpfh_gamma->Reset();
			tmpfh_delta->Reset();
			//tmpfh_alpha->Sumw2();
			//tmpfh_beta->Sumw2();
			//tmpfh_kb->Sumw2();
			//tmpfh_kcl->Sumw2();
			
			tmpfh_alpha->Divide( h1["eff_TaggedBothJets_cl"], h1["eff_TaggedJet_cl"]);
			tmpfh_beta->Divide( h1["eff_TaggedBothJets_b"], h1["eff_TaggedJet_b"]);

			tmpfh_delta->Divide( h1["eff_mu_taggedaway_b"], h1["eff_pTrel_b"]);
			tmpfh_gamma->Divide( h1["eff_mu_taggedaway_cl"], h1["eff_pTrel_cl"]);
			
			tmpfh_kb->Divide( h1["eff_pTrel_TaggedJet_b"], h1["eff_pTrel_b"]  );
			tmpfh_kb->Divide( h1["eff_TaggedJet_b"] );

			tmpfh_kcl->Divide( h1["eff_pTrel_TaggedJet_cl"], h1["eff_pTrel_cl"] );
			tmpfh_kcl->Divide(  h1["eff_TaggedJet_cl"] );

			// rebin only beta and kb!!
			fh_alpha = (TH1D*) tmpfh_alpha->Rebin(ncorrptarray-1,"fh_alpha",corrptbins);
			//fh_alpha->Divide( h1["eff_TaggedBothJets_cl"], h1["eff_TaggedJet_cl"]);
			fh_beta = (TH1D*)tmpfh_beta->Rebin(ncorrptarray-1,"fh_beta",corrptbins);
			fh_kb = (TH1D*)tmpfh_kb->Rebin(ncorrptarray-1,"fh_kb",corrptbins);
			fh_kcl = (TH1D*)tmpfh_kcl->Rebin(ncorrptarray-1,"fh_kcl",corrptbins);
			//fh_kcl->Divide( h1["eff_pTrel_TaggedJet_cl"], h1["eff_pTrel_cl"] );
			//fh_kcl->Divide(  h1["eff_TaggedJet_cl"] );
			
			std::cout << "rebinning done" << std::endl;
		} else {

			fh_alpha->Divide( h1["eff_TaggedBothJets_cl"], h1["eff_TaggedJet_cl"]);
			fh_beta->Divide( h1["eff_TaggedBothJets_b"], h1["eff_TaggedJet_b"]);
		
			fh_kb->Divide( h1["eff_pTrel_TaggedJet_b"], h1["eff_pTrel_b"]  );
			fh_kb->Divide( h1["eff_TaggedJet_b"] );

			fh_kcl->Divide( h1["eff_pTrel_TaggedJet_cl"], h1["eff_pTrel_cl"] );
			fh_kcl->Divide(  h1["eff_TaggedJet_cl"] );
			
			fh_delta->Divide( h1["eff_mu_taggedaway_b"], h1["eff_pTrel_b"] );
			fh_gamma->Divide( h1["eff_mu_taggedaway_cl"], h1["eff_pTrel_cl"] );

		}
			
		// fit to pol0
		fh_alpha->Fit("pol0","0");
		fh_beta->Fit("pol0","0");
		fh_kb->Fit("pol0","0");
		fh_kcl->Fit("pol0","0");
		fh_delta->Fit("pol0","0");
		fh_gamma->Fit("pol0","0");

	}
	else {

		// true efficiency
		feffTag_b = (TH1D*) gDirectory->Get("eff_TaggedJet_b");
		feffTag_cl = (TH1D*) gDirectory->Get("eff_TaggedJet_cl");
		feffmu_b = (TH1D*) gDirectory->Get("eff_pTrel_b");
		feffmu_cl = (TH1D*) gDirectory->Get("eff_pTrel_cl");
		feffTagmu_b = (TH1D*) gDirectory->Get("eff_pTrel_TaggedJet_b");
		feffTagmu_cl = (TH1D*) gDirectory->Get("eff_pTrel_TaggedJet_cl");
		
		if (fisCorrFile) finputCorrFile->cd();
		
		if (fcategory=="pT") {
			if (frebin) {
				
				TH1D* tmpfh_alpha = (TH1D*) gDirectory->Get("alpha");
				TH1D* tmpfh_beta = (TH1D*) gDirectory->Get("beta");
				TH1D* tmpfh_kb = (TH1D*) gDirectory->Get("kappa_b");
				TH1D* tmpfh_kcl = (TH1D*) gDirectory->Get("kappa_cl");
				//TH1D* tmpfh_gamma = (TH1D*) gDirectory->Get("gamma");
				//TH1D* tmpfh_delta = (TH1D*) gDirectory->Get("delta");
				
				fh_alpha = (TH1D*) tmpfh_alpha->Rebin(ncorrptarray-1,"fh_alpha",corrptbins);
				fh_beta = (TH1D*)tmpfh_beta->Rebin(ncorrptarray-1,"fh_beta",corrptbins);
				fh_kb = (TH1D*)tmpfh_kb->Rebin(ncorrptarray-1,"fh_kb",corrptbins);
				fh_kcl = (TH1D*)tmpfh_kcl->Rebin(ncorrptarray-1,"fh_kcl",corrptbins);
				// fit to pol0
				fh_alpha->Fit("pol0","0");
				fh_beta->Fit("pol0","0");
				fh_kb->Fit("pol0","0");
				fh_kcl->Fit("pol0","0");
			} else {
				fh_alpha = (TH1D*) gDirectory->Get("alpha");
				fh_beta = (TH1D*) gDirectory->Get("beta");
				fh_kb = (TH1D*) gDirectory->Get("kappa_b");
				fh_kcl = (TH1D*) gDirectory->Get("kappa_cl");
			}
		
		} else {
			if (frebin) {
				TH1D* tmpfh_alpha = (TH1D*) gDirectory->Get("alpha_eta");
				TH1D* tmpfh_beta = (TH1D*) gDirectory->Get("beta_eta");
				TH1D* tmpfh_kb = (TH1D*) gDirectory->Get("kappa_eta_b");
				TH1D* tmpfh_kcl = (TH1D*) gDirectory->Get("kappa_eta_cl");
				
				tmpfh_alpha->Rebin(ncorretaarray-1,"fh_alpha_eta",corretabins);
				tmpfh_beta->Rebin(ncorretaarray-1,"fh_beta_eta",corretabins);
				tmpfh_kb->Rebin(ncorretaarray-1,"fh_kb_eta",corretabins);
				tmpfh_kcl->Rebin(ncorretaarray-1,"fh_kcl_eta",corretabins);
				// fit to pol0
				fh_alpha->Fit("pol0","0");
				fh_beta->Fit("pol0","0");
				fh_kb->Fit("pol0","0");
				fh_kcl->Fit("pol0","0");
				
			} else {
				fh_alpha = (TH1D*) gDirectory->Get("alpha_eta");
				fh_beta = (TH1D*) gDirectory->Get("beta_eta");
				fh_kb = (TH1D*) gDirectory->Get("kappa_b_eta");
				fh_kcl = (TH1D*) gDirectory->Get("kappa_cl_eta");
			}
		}
	}
	//std::cout << " got correlations" << std::endl;
}


void S8Solver::GetInput() {

	this->LoadHistos();
	
	// integrated input
	TotalInput["n"] = fnHisto->Integral();
	TotalInput["nMu"] = fnHistoMu->Integral();
	TotalInput["p"] = fpHisto->Integral();
	TotalInput["pMu"] = fpHistoMu->Integral();
	TotalInput["nTag"] = fnHistoSvx->Integral();
	TotalInput["nMuTag"] = fnHistoAll->Integral();
	TotalInput["pTag"] = fpHistoSvx->Integral();
	TotalInput["pMuTag"] = fpHistoAll->Integral();

	// asumming parameters fitted to a constant
	TF1 *Fkb = fh_kb->GetFunction("pol0");
	TF1 *Fkcl = fh_kcl->GetFunction("pol0");
	TF1 *Falpha = fh_alpha->GetFunction("pol0");
	TF1 *Fbeta = fh_beta->GetFunction("pol0");
	TF1 *Fdelta = fh_delta->GetFunction("pol0");
	TF1 *Fgamma = fh_gamma->GetFunction("pol0");

	TotalInput["kappa_b"] = fKappabf * Fkb->GetParameter(0);
	TotalInput["kappa_cl"] = fKappaclf * Fkcl->GetParameter(0);
	TotalInput["alpha"] = fAlphaf * Falpha->GetParameter(0);
	TotalInput["beta"] = fBetaf * Fbeta->GetParameter(0);
	TotalInput["delta"] = Fdelta->GetParameter(0);
	TotalInput["gamma"] = Fgamma->GetParameter(0);

	//check
	/*
	TotalInput["n"] = 439977;
	TotalInput["nMu"] = 257449;
	TotalInput["p"] = 174425;
	TotalInput["pMu"] = 104614;
	TotalInput["nTag"] = 131090;
	TotalInput["nMuTag"] = 84172;
	TotalInput["pTag"] = 57112;
	TotalInput["pMuTag"] = 36576;
	TotalInput["kappa_b"] = 0.984;
	TotalInput["kappa_cl"] = 0.86;
	TotalInput["alpha"] = 0.98;
	TotalInput["beta"] = 0.973;
	*/
	
	// binned input base in the n samples
	std::vector< TH1* > HistoList;
	HistoList.push_back(fnHisto);
	HistoList.push_back(fpHisto);
	HistoList.push_back(fnHistoMu);
	HistoList.push_back(fpHistoMu);
	HistoList.push_back(fnHistoSvx);
	HistoList.push_back(fpHistoSvx);
	HistoList.push_back(fnHistoAll);
	HistoList.push_back(fpHistoAll);
	HistoList.push_back(fh_kb);
	HistoList.push_back(fh_kcl);
	HistoList.push_back(fh_alpha);
	HistoList.push_back(fh_beta);
	HistoList.push_back(fh_delta);
	HistoList.push_back(fh_gamma);

	std::vector< TString> name;
	name.push_back("n");
	name.push_back("nMu");
	name.push_back("p");
	name.push_back("pMu");
	name.push_back("nTag");
	name.push_back("nMuTag");
	name.push_back("pTag");
	name.push_back("pMuTag");
	name.push_back("kappa_b");
	name.push_back("kappa_cl");
	name.push_back("alpha");
	name.push_back("beta");
	name.push_back("delta");
	name.push_back("gamma");

	for (int ibin = 1; ibin<= fnHisto->GetNbinsX(); ++ibin) {

		std::map<TString, double> tmpmap;

		double pt = fnHisto->GetXaxis()->GetBinCenter(ibin);
		
		for ( size_t ihisto=0; ihisto!=HistoList.size(); ++ihisto) {
			
			TH1D *htemp = (TH1D*) HistoList[ihisto];

			if (name[ihisto]=="kappa_b"||name[ihisto]=="kappa_cl"||
			    name[ihisto]=="alpha"||name[ihisto]=="beta"||
			    name[ihisto]=="delta"||name[ihisto]=="gamma") {

				if(name[ihisto]=="beta") {
					if (!fBetaConst) {
						int abin = htemp->GetXaxis()->FindBin(pt);
						tmpmap[name[ihisto]] = fBetaf * htemp->GetBinContent(abin);
					} else {
						tmpmap[name[ihisto]] = fBetaf * TotalInput["beta"];
					}
				}
				else if(name[ihisto]=="alpha") { 
					if (!fAlphaConst) { 
						int abin = htemp->GetXaxis()->FindBin(pt); 
						tmpmap[name[ihisto]] = fAlphaf * htemp->GetBinContent(abin); 
					} else { 
						tmpmap[name[ihisto]] = fAlphaf * TotalInput["alpha"]; 
					} 
				}
				else if(name[ihisto]=="kappa_b") {
				  if (!fKappabConst) {
				    int abin = htemp->GetXaxis()->FindBin(pt);  
				    tmpmap[name[ihisto]] = fKappabf * htemp->GetBinContent(abin);  
				  } else {
				    tmpmap[name[ihisto]] = fKappabf * TotalInput["kappa_b"]; 
				  }
				}
				else if(name[ihisto]=="kappa_cl") {
				  if (!fKappaclConst) {
				    int abin = htemp->GetXaxis()->FindBin(pt);  
				    tmpmap[name[ihisto]] = fKappaclf * htemp->GetBinContent(abin);  
				  } else {
				    tmpmap[name[ihisto]] = fKappaclf * TotalInput["kappa_cl"]; 
				  }
				}
				else if(name[ihisto]=="delta") {
                                  if (!fDeltaConst) {
                                    int abin = htemp->GetXaxis()->FindBin(pt);
                                    tmpmap[name[ihisto]] = htemp->GetBinContent(abin);
                                  } else {
                                    tmpmap[name[ihisto]] = TotalInput["delta"];
                                  }
                                }
				else if(name[ihisto]=="gamma") {
                                  if (!fGammaConst) {
                                    int abin = htemp->GetXaxis()->FindBin(pt);
                                    tmpmap[name[ihisto]] = htemp->GetBinContent(abin);
                                  } else {
                                    tmpmap[name[ihisto]] = TotalInput["gamma"];
                                  }
                                }

			} else {
				tmpmap[name[ihisto]] = htemp->GetBinContent(ibin);
			}
		}

		BinnedInput[ibin] = tmpmap;		

	}
	

}

void S8Solver::Solve() {

	this->GetInput();
	
	if (fmethod=="analytic" || fmethod=="fit" ) {

		S8AnalyticSolver sol;
		
		// average solution
		
		sol.Solve(TotalInput);

		fTotalSolution = sol.GetSolution();
		fTotalSolutionErr = sol.GetSolutionErr();

		// now fit
		if (fmethod=="fit") {
			//std::cout << "fit total" << std::endl;
			S8FitSolver s8fit;
			s8fit.Init(fTotalSolution);
			//std::cout << "fit init done" << std::endl;
			s8fit.Solve(TotalInput);
			//std::cout << "fit done" << std::endl;
			fTotalSolution = s8fit.GetSolution();
			fTotalSolutionErr = s8fit.GetSolutionErr();
		}
		
			// binned solution
			for( std::map<int,std::map<TString,double> >::const_iterator ibin = BinnedInput.begin(); ibin!=BinnedInput.end(); ++ibin) {
				sol.Solve(ibin->second);
				fBinnedSolution[ibin->first] = sol.GetSolution();
				fBinnedSolutionErr[ibin->first] = sol.GetSolutionErr();
				
				// now fit
				if (fmethod=="fit") {
					S8FitSolver as8fit;
					as8fit.Init(fBinnedSolution[ibin->first]);
					as8fit.Solve(ibin->second);
					fBinnedSolution[ibin->first] = as8fit.GetSolution();
					fBinnedSolutionErr[ibin->first] = as8fit.GetSolutionErr();
				}
			}
		
	}
	if (fmethod=="numeric") {

		// average solution
		Double_t inputs[8];
		inputs[0] = TotalInput["n"];
		inputs[1] = TotalInput["nMu"];
		inputs[2] = TotalInput["nTag"];
		inputs[3] = TotalInput["p"];
		inputs[4] = TotalInput["nMuTag"];
		inputs[5] = TotalInput["pTag"];
		inputs[6] = TotalInput["pMu"];
		inputs[7] = TotalInput["pMuTag"];
		
		S8NumericSolver sol("sol",inputs);
		//sol.SetCorr(TotalInput["kappa_b"],1.,TotalInput["beta"],
		//	    TotalInput["kappa_cl"],1.,TotalInput["alpha"]);
		sol.SetCorr(TotalInput["kappa_b"],TotalInput["beta"],1.,
			    TotalInput["kappa_cl"],TotalInput["alpha"],1.);
		//sol.SetCorr(TotalInput["kappa_b"],TotalInput["beta"],TotalInput["delta"],
		//    TotalInput["kappa_cl"],TotalInput["alpha"],TotalInput["gamma"]);

		sol.SetCorrError(0.,0.,0.,0.,0.,0.);
		sol.SetError(2);
		sol.SetNbErrorIteration(200);
		sol.SetInitialOrder(1,1);
		bool converge = true;
		
		if(!sol.Solve()) {
			//converge = false;
			sol.SetInitialOrder(0,-1);   // OK en principe...
			if(!sol.Solve()) {
				//converge = false;
				sol.SetInitialOrder(3,1);   // OK en principe...
				if(!sol.Solve()) {
					//converge = false;
					sol.SetInitialOrder(2,1);   // OK en principe...
					if(!sol.Solve())
						//converge = false;
						std::cout << "Arrggg : Impossible de determiner une solution physique I" << std::endl;
				}
				
			}
		}

		if (converge) {
			fTotalSolution["n_b"]       = sol.GetResult("na")*TotalInput["n"];
			fTotalSolution["n_cl"]      = sol.GetResult("nb")*TotalInput["n"];
			fTotalSolution["effMu_b"]   = sol.GetResult("ea1");
			fTotalSolution["effMu_cl"]  = sol.GetResult("eb1");
			fTotalSolution["effTag_b"]  = sol.GetResult("ea2");
			fTotalSolution["effTag_cl"] = sol.GetResult("eb2");
			fTotalSolution["p_b"]       = sol.GetResult("ea3")*fTotalSolution["n_b"];
			fTotalSolution["p_cl"]      = sol.GetResult("eb3")*fTotalSolution["n_cl"];
		
			// FIX errors
			fTotalSolutionErr["n_b"]       = (sol.GetErrorSup("na")+sol.GetErrorInf("na"))/2.;
			fTotalSolutionErr["n_cl"]      = (sol.GetErrorSup("nb")+sol.GetErrorInf("nb"))/2.;
			fTotalSolutionErr["effMu_b"]   = (sol.GetErrorSup("ea1")+sol.GetErrorInf("ea1"))/2.;
			fTotalSolutionErr["effMu_cl"]  = (sol.GetErrorSup("eb1")+sol.GetErrorInf("eb1"))/2.;
			fTotalSolutionErr["effTag_b"]  = (sol.GetErrorSup("ea2")+sol.GetErrorInf("ea2"))/2.;
			fTotalSolutionErr["effTag_cl"] = (sol.GetErrorSup("eb2")+sol.GetErrorInf("eb2"))/2.;
			fTotalSolutionErr["p_b"]       = (sol.GetErrorSup("ea3")+sol.GetErrorInf("ea3"))/2.;
			fTotalSolutionErr["p_cl"]      =  (sol.GetErrorSup("eb3")+sol.GetErrorInf("eb3"))/2.;
		} else {
			for( std::map<TString,double>::const_iterator ii=fTotalSolution.begin(); ii!=fTotalSolution.end(); ++ii) {
				fTotalSolution[ii->first] = 0.;
				fTotalSolutionErr[ii->first] = 0.;
			}
		}
		
		// binned solutions
		for( std::map<int,std::map<TString,double> >::const_iterator ibin = BinnedInput.begin(); ibin!=BinnedInput.end(); ++ibin) {

			converge = true;
			std::map<TString,double> tmpinput;
			tmpinput = ibin->second;
			
			inputs[0] = tmpinput["n"];
			inputs[1] = tmpinput["nMu"];
			inputs[2] = tmpinput["nTag"];
			inputs[3] = tmpinput["p"];
			inputs[4] = tmpinput["nMuTag"];
			inputs[5] = tmpinput["pTag"];
			inputs[6] = tmpinput["pMu"];
			inputs[7] = tmpinput["pMuTag"];
			
			S8NumericSolver solu("solu",inputs);
			//solu.SetCorr(tmpinput["kappa_b"],1.,tmpinput["beta"],
			//     tmpinput["kappa_cl"],1.,tmpinput["alpha"]);
			//solu.SetCorr(tmpinput["kappa_b"],tmpinput["beta"],1., 
			//tmpinput["kappa_cl"],tmpinput["alpha"],1.); 
			solu.SetCorr(tmpinput["kappa_b"],tmpinput["beta"],tmpinput["delta"],
						 tmpinput["kappa_cl"],tmpinput["alpha"],tmpinput["gamma"]);

			solu.SetCorrError(0.,0.,0.,0.,0.,0.);
			solu.SetError(2);
			solu.SetNbErrorIteration(200);
			solu.SetInitialOrder(1,1);
			if(!solu.Solve()) {
				//converge = false;
				solu.SetInitialOrder(0,0);   // OK en principe...
				if(!solu.Solve()) {
					//converge = false;
					solu.SetInitialOrder(3,1);   // OK en principe...
					if(!solu.Solve()) {
						//converge = false;
						solu.SetInitialOrder(2,1);   // OK en principe...
						if(!solu.Solve())
							//converge = false;
							std::cout << "Arrggg : Impossible de determiner une solution physique" << std::endl;
					}
				}
			}

			std::map<TString,double> tmpsolu;
			std::map<TString,double> tmpsoluerr;
				
			if (converge) {
				
				tmpsolu["n_b"]       = solu.GetResult("na")*tmpinput["n"];
				tmpsolu["n_cl"]      = solu.GetResult("nb")*tmpinput["n"];
				tmpsolu["effMu_b"]   = solu.GetResult("ea1");
				tmpsolu["effMu_cl"]  = solu.GetResult("eb1");
				tmpsolu["effTag_b"]  = solu.GetResult("ea2");
				tmpsolu["effTag_cl"] = solu.GetResult("eb2");
				tmpsolu["p_b"]       = solu.GetResult("ea3")*tmpsolu["n_b"];
				tmpsolu["p_cl"]      = solu.GetResult("eb3")*tmpsolu["n_cl"];
							
				tmpsoluerr["n_b"]       = (solu.GetErrorSup("na")+solu.GetErrorInf("na"))/2.;
				tmpsoluerr["n_cl"]      = (solu.GetErrorSup("nb")+solu.GetErrorInf("nb"))/2.;
				tmpsoluerr["effMu_b"]   = (solu.GetErrorSup("ea1")+solu.GetErrorInf("ea1"))/2.;
				tmpsoluerr["effMu_cl"]  = (solu.GetErrorSup("eb1")+solu.GetErrorInf("eb1"))/2.;
				tmpsoluerr["effTag_b"]  = (solu.GetErrorSup("ea2")+solu.GetErrorInf("ea2"))/2.;
				tmpsoluerr["effTag_cl"] = (solu.GetErrorSup("eb2")+solu.GetErrorInf("eb2"))/2.;
				tmpsoluerr["p_b"]       = (solu.GetErrorSup("ea3")+solu.GetErrorInf("ea3"))/2.;
				tmpsoluerr["p_cl"]      = (solu.GetErrorSup("eb3")+solu.GetErrorInf("eb3"))/2.;
							
			} else {
				for( std::map<TString,double>::const_iterator ii=fTotalSolution.begin(); ii!=fTotalSolution.end(); ++ii) {
					tmpsolu[ii->first] = 0.;
					tmpsoluerr[ii->first] = 0.;
				}
			}
			fBinnedSolution[ibin->first] = tmpsolu;
			fBinnedSolutionErr[ibin->first] = tmpsoluerr;

		}
	}

}

void S8Solver::PrintData(TString option) {

	std::cout << " Name: " << fthename << std::endl;
	std::cout << " Method: " << fmethod << std::endl;
	std::cout << " Category: " << fcategory << std::endl;
	std::cout << " alpha constant: " << fAlphaConst << std::endl;
	std::cout << " beta constant: " << fBetaConst << std::endl;
	std::cout << " k_b constant: " << fKappabConst << std::endl;
	std::cout << " k_cl constant: " << fKappaclConst << std::endl;
	std::cout << " max Ptrel: " << fMaxPtrel << std::endl;
	std::cout << " rebin correlation histograms: " << frebin << std::endl;
	std::cout << " recalculate correlation histograms: " << frecalculateFactors << std::endl;
	std::cout << " alpha scale factor: " << fAlphaf << std::endl;
	std::cout << " beta scale factor: " << fBetaf << std::endl;
	std::cout << " k_b scale factor: " << fKappabf << std::endl;
	std::cout << " k_cl scale factor: " << fKappaclf << std::endl;
	
	if (option=="input") {
		std::cout << " Input values" << std::endl;
		std::cout << " Average: " << std::endl;
		for( std::map<TString,double>::const_iterator i = TotalInput.begin(); i!=TotalInput.end(); ++i) {
			std::cout << i->first << " = " << i->second << std::endl;
		}
		std::cout << " Binned: " << std::endl;
		for( std::map<int,std::map< TString, double > >::const_iterator ibin = BinnedInput.begin(); ibin!=BinnedInput.end(); ++ibin) {

			std::cout << "### bin: " << ibin->first << std::endl;
			std::map<TString, double> tmpmap = ibin->second;
			
			for( std::map<TString,double>::const_iterator i = tmpmap.begin(); i!=tmpmap.end(); ++i) {
				std::cout << i->first << " = " << i->second << std::endl;
			}
		}
		
	}
	else {
		std::cout << " Solutions: " << std::endl;
		
		std::cout << " Average: " << std::endl;
		for( std::map<TString,double>::const_iterator i = fTotalSolution.begin(); i!=fTotalSolution.end(); ++i) {
			std::cout << i->first << " = " << i->second << " \\pm " << fTotalSolutionErr[i->first] << std::endl;
		}
		std::cout << " Binned: " << std::endl;
		for( std::map<int,std::map<TString,double> >::const_iterator ibin = fBinnedSolution.begin(); ibin!=fBinnedSolution.end(); ++ibin) {
			
			std::cout << "### bin: " << ibin->first << std::endl;
			std::map<TString, double> tmpmap = ibin->second;
			std::map<TString, double> tmpmaperr = fBinnedSolutionErr[ibin->first];
			for( std::map<TString,double>::const_iterator i = tmpmap.begin(); i!=tmpmap.end(); ++i) {
				std::cout << i->first << " = " << i->second << " \\pm " << tmpmaperr[i->first] << std::endl;
			}
		}
	}
	
}

void S8Solver::Draw(int maxNbins) {

	
	
  Int_t nxbins = fBinnedSolution.size();

  if (maxNbins != 0 ) nxbins = maxNbins;
  
  TArrayD ptarray(nxbins);
  TArrayD ptarrayErr(nxbins);
  TArrayD s8effTag_b(nxbins);
  TArrayD s8effTag_bErr(nxbins);

  
  TArrayD effTag_b(nxbins);
  TArrayD effTag_bErr(nxbins);
  TArrayD effTag_cl(nxbins);
  TArrayD effTag_clErr(nxbins);
  TArrayD effmu_b(nxbins);
  TArrayD effmu_bErr(nxbins);
  TArrayD effmu_cl(nxbins);
  TArrayD effmu_clErr(nxbins);
  TArrayD effTagmu_b(nxbins);
  TArrayD effTagmu_bErr(nxbins);
  TArrayD effTagmu_cl(nxbins);
  TArrayD effTagmu_clErr(nxbins);

  TArrayD input_n(nxbins);
  TArrayD input_p(nxbins);
  TArrayD input_ntag(nxbins);
  TArrayD input_ptag(nxbins);
  TArrayD input_nmu(nxbins);
  TArrayD input_pmu(nxbins);
  TArrayD input_ntagmu(nxbins);
  TArrayD input_ptagmu(nxbins);

  TArrayD input_nErr(nxbins);
  TArrayD input_pErr(nxbins);
  TArrayD input_ntagErr(nxbins);
  TArrayD input_ptagErr(nxbins);
  TArrayD input_nmuErr(nxbins);
  TArrayD input_pmuErr(nxbins);
  TArrayD input_ntagmuErr(nxbins);
  TArrayD input_ptagmuErr(nxbins);
  
	for (int ibin = 1; ibin<= nxbins; ++ibin) {

		ptarray[ibin-1] = fnHisto->GetXaxis()->GetBinCenter(ibin);
		ptarrayErr[ibin-1] = 0.5 * fnHisto->GetXaxis()->GetBinWidth(ibin);
		effTag_b[ibin-1] = feffTag_b->GetBinContent(ibin);
		effTag_bErr[ibin-1] = feffTag_b->GetBinError(ibin);
		
	}
	// fnHisto->GetXaxis()->GetXbins()->GetArray();

	// binned input
	int jjbin= 0;
	for( std::map<int,std::map<TString,double> >::const_iterator ibin = BinnedInput.begin(); ibin!=BinnedInput.end(); ++ibin) {

		std::map<TString,double> tmpinput;
		tmpinput = ibin->second;
			
		input_n[jjbin] = tmpinput["n"];
		input_nmu[jjbin] = tmpinput["nMu"];
		input_ntag[jjbin] = tmpinput["nTag"];
		input_p[jjbin] = tmpinput["p"];
		input_ntagmu[jjbin] = tmpinput["nMuTag"];
		input_ptag[jjbin] = tmpinput["pTag"];
		input_pmu[jjbin] = tmpinput["pMu"];
		input_ptagmu[jjbin] = tmpinput["pMuTag"];

		input_nErr[jjbin] = sqrt(tmpinput["n"]);
		input_nmuErr[jjbin] = sqrt(tmpinput["nMu"]);
		input_ntagErr[jjbin] = sqrt(tmpinput["nTag"]);
		input_pErr[jjbin] = sqrt(tmpinput["p"]);
		input_ntagmuErr[jjbin] = sqrt(tmpinput["nMuTag"]);
		input_ptagErr[jjbin] = sqrt( tmpinput["pTag"]);
		input_pmuErr[jjbin] = sqrt(tmpinput["pMu"]);
		input_ptagmuErr[jjbin] = sqrt(tmpinput["pMuTag"]);

		jjbin++;
	}
	
	for( std::map<int,std::map<TString,double> >::const_iterator ibin = fBinnedSolution.begin(); ibin!=fBinnedSolution.end(); ++ibin) {
			
		//std::cout << " bin # " << ibin->first << std::endl;
		std::map<TString, double> tmpmaps8 = ibin->second;
		std::map<TString, double> tmpmaps8err = fBinnedSolutionErr[ibin->first];

		if (ibin->first <= nxbins) {
			s8effTag_b[ibin->first -1] = tmpmaps8["effTag_b"];
			s8effTag_bErr[ibin->first -1] = tmpmaps8err["effTag_b"];
		}
		//for( std::map<TString,double>::const_iterator i = tmpmap.begin(); i!=tmpmap.end(); ++i) {
			//std::cout << i->first << " = " << i->second << " \\pm " << tmpmaperr[i->first] << std::endl;
			
		//}
	}
	
	geffTag_b = new TGraphErrors(nxbins,ptarray.GetArray(),effTag_b.GetArray(),ptarrayErr.GetArray(),effTag_bErr.GetArray());
	gS8effTag_b = new TGraphErrors(nxbins,ptarray.GetArray(),s8effTag_b.GetArray(),ptarrayErr.GetArray(),s8effTag_bErr.GetArray());

	geffTag_b->SetTitle("True b-efficiency");
	gS8effTag_b->SetTitle("System8 Results");
	geffTag_b->SetName("geffTag_b");
	gS8effTag_b->SetName("gS8effTag_b");
	
	geffTag_b->SetMarkerStyle(8);
	gS8effTag_b->SetMarkerStyle(8);
	geffTag_b->SetMarkerColor(1);
	gS8effTag_b->SetMarkerColor(2);
	geffTag_b->SetLineColor(1); 
	gS8effTag_b->SetLineColor(2); 
	geffTag_b->SetMarkerSize(1.5); 
	gS8effTag_b->SetMarkerSize(1.5); 

	TString prefix = "effTag_b";	
	cv_map[prefix+"_"+fthename] = new TCanvas(prefix+"_"+fthename,prefix+"_"+fthename,700,700); 

	TMultiGraph *multi_eff_b = new TMultiGraph();
	multi_eff_b->Add(geffTag_b,"p");
	multi_eff_b->Add(gS8effTag_b,"p");
	multi_eff_b->Draw("a");
	multi_eff_b->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
	multi_eff_b->GetYaxis()->SetTitle("b-jet Efficiency");
	
	TLegend *legb = new TLegend(0.57,0.22,0.87,0.38,"","NDC");
	legb->SetMargin(0.12);
	legb->SetTextSize(0.027);
	legb->SetFillColor(10);
	legb->AddEntry(geffTag_b,"True b-efficiency","P");
	legb->AddEntry(gS8effTag_b,"System8 Results","P");
	legb->Draw();
	gPad->SetGrid();

	TString acvname = "correlations_"+fthename;
	cv_map[acvname] = new TCanvas(acvname,acvname,700,700);
	fh_kb->SetMarkerColor(1); 
	fh_kcl->SetMarkerColor(2); 
	fh_beta->SetMarkerColor(3); 
	fh_alpha->SetMarkerColor(4);
	fh_kb->SetLineColor(1); 
	fh_kcl->SetLineColor(2); 
	fh_beta->SetLineColor(3); 
	fh_alpha->SetLineColor(4);
	fh_kb->SetMinimum(0.7);
	fh_kb->SetMaximum(1.3);
	fh_kb->Draw("P");
	fh_kcl->Draw("Psame");
	fh_beta->Draw("Psame");
	fh_alpha->Draw("Psame");

	TLegend *legc = new TLegend(0.77,0.22,0.87,0.38,"","NDC");
	//legc->SetMargin(0.12);
	legc->SetTextSize(0.029);
	legc->SetFillColor(10);
	legc->AddEntry(fh_kb,"#kappa_{b}","P");
	legc->AddEntry(fh_kcl,"#kappa_{cl}","P");
	legc->AddEntry(fh_alpha,"#alpha","P");
	legc->AddEntry(fh_beta,"#beta","P");
	legc->Draw();
	gPad->SetGrid();

	acvname = "xcorrelations_"+fthename;
        cv_map[acvname] = new TCanvas(acvname,acvname,700,700);
        fh_delta->SetMarkerColor(1);
        fh_gamma->SetMarkerColor(2);
        fh_delta->SetLineColor(1);
        fh_gamma->SetLineColor(2);
	fh_delta->SetMinimum(0.7);
        fh_delta->SetMaximum(1.3);
        fh_delta->Draw("P");
        fh_gamma->Draw("Psame");
        
        TLegend *legcc = new TLegend(0.77,0.22,0.87,0.38,"","NDC");
        //legc->SetMargin(0.12);
        legcc->SetTextSize(0.029);
        legcc->SetFillColor(10);
        legcc->AddEntry(fh_delta,"#delta","P");
        legcc->AddEntry(fh_gamma,"#gamma","P");
	legcc->Draw();
        gPad->SetGrid();

	ginput_n   = new TGraphErrors(nxbins,ptarray.GetArray(),input_n.GetArray(),ptarrayErr.GetArray(),0);
	ginput_nmu = new TGraphErrors(nxbins,ptarray.GetArray(),input_nmu.GetArray(),ptarrayErr.GetArray(),0);
	ginput_ntag   = new TGraphErrors(nxbins,ptarray.GetArray(),input_ntag.GetArray(),ptarrayErr.GetArray(),0);
	ginput_ntagmu = new TGraphErrors(nxbins,ptarray.GetArray(),input_ntagmu.GetArray(),ptarrayErr.GetArray(),0);
	ginput_p   = new TGraphErrors(nxbins,ptarray.GetArray(),input_p.GetArray(),ptarrayErr.GetArray(),0);
	ginput_pmu = new TGraphErrors(nxbins,ptarray.GetArray(),input_pmu.GetArray(),ptarrayErr.GetArray(),0);
	ginput_ptag   = new TGraphErrors(nxbins,ptarray.GetArray(),input_ptag.GetArray(),ptarrayErr.GetArray(),0);
	ginput_ptagmu = new TGraphErrors(nxbins,ptarray.GetArray(),input_ptagmu.GetArray(),ptarrayErr.GetArray(),0);

	ginput_n->SetMarkerStyle(8);
	ginput_nmu->SetMarkerStyle(8);
	ginput_ntag->SetMarkerStyle(8);
	ginput_ntagmu->SetMarkerStyle(8);
	ginput_p->SetMarkerStyle(8);
	ginput_pmu->SetMarkerStyle(8);
	ginput_ptag->SetMarkerStyle(8);
	ginput_ptagmu->SetMarkerStyle(8);

	ginput_n->SetMarkerColor(1);
	ginput_nmu->SetMarkerColor(2);
	ginput_ntag->SetMarkerColor(3);
	ginput_ntagmu->SetMarkerColor(4);
	ginput_p->SetMarkerColor(1);
	ginput_pmu->SetMarkerColor(2);
	ginput_ptag->SetMarkerColor(3);
	ginput_ptagmu->SetMarkerColor(4);

	ginput_n->SetLineColor(1);
	ginput_nmu->SetLineColor(2);
	ginput_ntag->SetLineColor(3);
	ginput_ntagmu->SetLineColor(4);
	ginput_p->SetLineColor(1);
	ginput_pmu->SetLineColor(2);
	ginput_ptag->SetLineColor(3);
	ginput_ptagmu->SetLineColor(4);

	acvname = "ninputs_"+fthename;
	cv_map[acvname] = new TCanvas(acvname,acvname,700,700);
	
	TMultiGraph *multi_n = new TMultiGraph();
	multi_n->Add(ginput_n,"p");
	multi_n->Add(ginput_nmu,"p");
	multi_n->Add(ginput_ntag,"p");
	multi_n->Add(ginput_ntagmu,"p");
	multi_n->Draw("a");
	multi_n->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
	multi_n->GetYaxis()->SetTitle("Events");
	gPad->SetLogy();
	TLegend *legd = new TLegend(0.79,0.69,0.89,0.86,"","NDC");
	//legd->SetMargin(0.12);
	legd->SetTextSize(0.029);
	legd->SetFillColor(10);
	legd->AddEntry(ginput_n,"n","P");
	legd->AddEntry(ginput_nmu,"n^{mu}","P");
	legd->AddEntry(ginput_ntag,"n^{tag}","P");
	legd->AddEntry(ginput_ntagmu,"n^{tag,mu}","P");
	legd->Draw();
			
	acvname = "pinputs_"+fthename;
	cv_map[acvname] = new TCanvas(acvname,acvname,700,700);
	
	TMultiGraph *multi_p = new TMultiGraph();
	multi_p->Add(ginput_p,"p");
	multi_p->Add(ginput_pmu,"p");
	multi_p->Add(ginput_ptag,"p");
	multi_p->Add(ginput_ptagmu,"p");
	multi_p->Draw("a");
	multi_p->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
	multi_p->GetYaxis()->SetTitle("Events");
	gPad->SetLogy();
	TLegend *lege = new TLegend(0.79,0.69,0.89,0.86,"","NDC");
	//lege->SetMargin(0.12);
	lege->SetTextSize(0.029);
	lege->SetFillColor(10);
	lege->AddEntry(ginput_p,"p","P");
	lege->AddEntry(ginput_pmu,"p^{mu}","P");
	lege->AddEntry(ginput_ptag,"p^{tag}","P");
	lege->AddEntry(ginput_ptagmu,"p^{tag,mu}","P");
	lege->Draw();
	
}

void S8Solver::Print(TString extension ) {

	for(std::map<TString,TCanvas*>::const_iterator icv=cv_map.begin(); icv!=cv_map.end(); ++icv){

		TString tmpname = icv->first;
		TCanvas *acv = icv->second;
		acv->Print(TString(tmpname+"."+extension));
	}

}

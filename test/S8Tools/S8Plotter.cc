#define S8Plotter_cxx
#include "S8Plotter.h"

#include "TArrayD.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include <iostream>
#include <iomanip>

ClassImp(S8Plotter)


//______________________________________________________________
bool SortByMaximum( TH1* lh, TH1* rh ) {

	return lh->GetMaximum() > rh->GetMaximum();

}

//______________________________________________________________
void S8Plotter::Book() {

	quark_label.push_back("");
	quark_label.push_back("b");
	quark_label.push_back("c");
	quark_label.push_back("uds");
	quark_label.push_back("g");

	quark_color[""] = 1;
	quark_color["b"] = 2;
	quark_color["c"] = 3;
	quark_color["uds"] = 4;
	quark_color["g"] = 6;

	cut_label["cut0"] = "";
	cut_label["cut1"] = "";
	cut_label["cut1"] = "";
	cut_label["cut1"] = "";
	
	std::string hname;
	std::string htitle;

	hname = "njets"; 
	htitle = "number of jets";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),40,0.,40.);
	
	hname = "nmuons"; 
	htitle = "number of muons";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),12,0.,12.);

	hname = "muon_pt"; 
	htitle = "muon p_{T} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,80.);
	
	hname = "muon_eta"; 
	htitle = "muon |#eta|";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,-5., 5.);

	htitle = "jet p_{T} [GeV/c]";
	for( size_t i = 0; i != quark_label.size(); ++i ) {
		hname = "jet_pt"; 
		if (i!=0) hname += "_"+quark_label[i];
		h1[hname] = new TH1D(hname.c_str(),htitle.c_str(), fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray() );
		h1[hname]->SetLineColor( quark_color[quark_label[i]] );
	}
	htitle = "jet p_{T} [GeV/c]";
	for( size_t i = 0; i != quark_label.size(); ++i ) {
		hname = "taggedjet_pt"; 
		if (i!=0) hname += "_"+quark_label[i];
		h1[hname] = new TH1D(hname.c_str(),htitle.c_str(), fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray() );
		h1[hname]->SetLineColor( quark_color[quark_label[i]] );
	}
	htitle = "jet |#eta|";
	for( size_t i = 0; i != quark_label.size(); ++i ) {
		hname = "jet_eta"; 
		if (i!=0) hname += "_"+quark_label[i];
		h1[hname] = new TH1D(hname.c_str(),htitle.c_str(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() );
		h1[hname]->SetLineColor( quark_color[quark_label[i]] );
	}
	htitle = "jet |#eta|";
	for( size_t i = 0; i != quark_label.size(); ++i ) {
		hname = "taggedjet_eta"; 
		if (i!=0) hname += "_"+quark_label[i];
		h1[hname] = new TH1D(hname.c_str(),htitle.c_str(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() );
		h1[hname]->SetLineColor( quark_color[quark_label[i]] );
	}
	
	hname = "jet_deltaR"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);

	hname = "jet_deltaR_b"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);
	h1[hname]->SetLineColor(quark_color["b"]);

	hname = "jet_deltaR_c"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);
	h1[hname]->SetLineColor(quark_color["c"]);

	hname = "jet_deltaR_udsg"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);
	h1[hname]->SetLineColor(quark_color["udsg"]);

	hname = "jet_deltaphi"; 
	htitle = "#Delta #phi";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,3.2);

	hname = "jet_deltaphi_b"; 
	htitle = "#Delta #phi";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,3.2);
	h1[hname]->SetLineColor(quark_color["b"]);
	
	hname = "jet_pTrel"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_b"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_c"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_udsg"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_pt50"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_pt100"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_pt150"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_b_pt50"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_b_pt100"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_b_pt150"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_c_pt50"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_c_pt100"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_c_pt150"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_udsg_pt50"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_udsg_pt100"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_udsg_pt150"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_mupt6"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_mupt10"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_mupt20"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_b_mupt6"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_b_mupt10"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_b_mupt20"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_c_mupt6"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_c_mupt10"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_c_mupt20"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_udsg_mupt6"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);
	
	hname = "jet_pTrel_udsg_mupt10"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_pTrel_udsg_mupt20"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),50,0.,5.);

	hname = "jet_ntrks"; 
	htitle = "jet_ntrks";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),130,0.,260.);
}

void S8Plotter::Loop()
{

	// setup tree
	Init();
	// book histograms
	Book();

	// book additional histograms
	int nptbins = fJetPtAxis.GetNbins();
	const Double_t *jetptbins = (fJetPtAxis.GetXbins())->GetArray();
	int netabins = fJetEtaAxis.GetNbins();
	const Double_t *jetetabins = (fJetEtaAxis.GetXbins())->GetArray();
	// int ncorrptbins = fCorrPtAxis.GetNbins();
	// const Double_t *corrptbins = (fCorrPtAxis.GetXbins())->GetArray();
	// int ncorretabins = fCorrEtaAxis.GetNbins();
	// const Double_t *corretabins = (fCorrEtaAxis.GetXbins())->GetArray();
	
	h2["npT"] = new TH2F("npT","MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["ppT"] = new TH2F("ppT","MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["ncmbpT"] = new TH2F("ncmbpT","opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["pcmbpT"] = new TH2F("pcmbpT","opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.);

	h2["qpT"] = new TH2F("qpT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["qcmbpT"] = new TH2F("qcmbpT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	
	
	h2["nEta"] = new TH2F("nEta","MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["pEta"] = new TH2F("pEta","MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["ncmbEta"] = new TH2F("ncmbEta","opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["pcmbEta"] = new TH2F("pcmbEta","opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.);

	h2["qEta"] = new TH2F("qEta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["qcmbEta"] = new TH2F("qcmbEta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.);
	
	h2["b_npT"] = new TH2F("b_npT","b MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["b_ppT"] = new TH2F("b_ppT","b MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["b_ncmbpT"] = new TH2F("b_ncmbpT","b opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["b_pcmbpT"] = new TH2F("b_pcmbpT","b opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.);

	h2["b_qpT"] = new TH2F("b_qpT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["b_qcmbpT"] = new TH2F("b_qcmbpT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	
	h2["b_nEta"] = new TH2F("b_nEta","b MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["b_pEta"] = new TH2F("b_pEta","b MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["b_ncmbEta"] = new TH2F("b_ncmbEta","b opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["b_pcmbEta"] = new TH2F("b_pcmbEta","b opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.);

	h2["b_qEta"] = new TH2F("b_qEta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["b_qcmbEta"] = new TH2F("b_qcmbEta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.);
	
	
	h2["cl_npT"] = new TH2F("cl_npT","cl MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["cl_ppT"] = new TH2F("cl_ppT","cl MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["cl_ncmbpT"] = new TH2F("cl_ncmbpT","cl opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["cl_pcmbpT"] = new TH2F("cl_pcmbpT","cl opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.);

	h2["cl_qpT"] = new TH2F("cl_qpT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["cl_qcmbpT"] = new TH2F("cl_qcmbpT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	
	h2["c_npT"] = new TH2F("c_npT","c MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["c_ncmbpT"] = new TH2F("c_ncmbpT","c opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["g_npT"] = new TH2F("g_npT","g MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["g_ncmbpT"] = new TH2F("g_ncmbpT","g opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["uds_npT"] = new TH2F("uds_npT","uds MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);
	h2["uds_ncmbpT"] = new TH2F("uds_ncmbpT","uds opp tag MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.);

	h2["cl_nEta"] = new TH2F("cl_nEta","cl MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["cl_pEta"] = new TH2F("cl_pEta","cl MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["cl_ncmbEta"] = new TH2F("cl_ncmbEta","cl opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["cl_pcmbEta"] = new TH2F("cl_pcmbEta","cl opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.);

	h2["cl_qEta"] = new TH2F("cl_qEta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.);
	h2["cl_qcmbEta"] = new TH2F("cl_qcmbEta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.);
	
	h2["ptVsEta"] = new TH2F("ptVsEta","pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() );
	h2["ptVsEta_b"] = new TH2F("ptVsEta_b","pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() );
	h2["taggedjet_ptVsEta"] = new TH2F("taggedjet_ptVsEta","taggedjet pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() );
	h2["taggedjet_ptVsEta_b"] = new TH2F("taggedjet_ptVsEta_b","taggedjet pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() );

	h1["alpha"] = new TH1D("alpha","alpha",nptbins,jetptbins);
	h1["beta"] = new TH1D("beta","beta",nptbins,jetptbins);
	h1["kappa_cl"] = new TH1D("kappa_cl","kappa_cl",nptbins,jetptbins);
	h1["kappa_b"] = new TH1D("kappa_b","kappa_b",nptbins,jetptbins);
	h1["delta"] = new TH1D("delta","delta",nptbins,jetptbins);
	h1["gamma"] = new TH1D("gamma","gamma",nptbins,jetptbins);

	h1["alpha_eta"] = new TH1D("alpha_eta","alpha_eta",netabins,jetetabins);
	h1["beta_eta"] = new TH1D("beta_eta","beta_eta",netabins,jetetabins);
	h1["kappa_eta_cl"] = new TH1D("kappa_eta_cl","kappa_eta_cl",netabins,jetetabins);
	h1["kappa_eta_b"] = new TH1D("kappa_eta_b","kappa_eta_b",netabins,jetetabins);
	h1["delta_eta"] = new TH1D("delta_eta","delta_eta",netabins,jetetabins);
	h1["gamma_eta"] = new TH1D("gamma_eta","gamma_eta",netabins,jetetabins);
/*
	h1["alpha"] = new TH1D("alpha","alpha",ncorrptbins,corrptbins);
	h1["beta"] = new TH1D("beta","beta",ncorrptbins,corrptbins);
	h1["kappa_cl"] = new TH1D("kappa_cl","kappa_cl",ncorrptbins,corrptbins);
	h1["kappa_b"] = new TH1D("kappa_b","kappa_b",ncorrptbins,corrptbins);

	h1["alpha_eta"] = new TH1D("alpha_eta","alpha_eta",ncorretabins,corretabins);
	h1["beta_eta"] = new TH1D("beta_eta","beta_eta",ncorretabins,corretabins);
	h1["kappa_eta_cl"] = new TH1D("kappa_eta_cl","kappa_eta_cl",ncorretabins,corretabins);
	h1["kappa_eta_b"] = new TH1D("kappa_eta_b","kappa_eta_b",ncorretabins,corretabins);
*/
	// enable errors
	for(std::map<std::string,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih){
		TH2 *htemp = ih->second;
		htemp->Sumw2();
	}
	for(std::map<std::string,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih){
		TH1 *htemp = ih->second;
		htemp->Sumw2();
	}

	// variable definitions
	std::map<std::string, double> totalTagged; // total tagged jets
	std::map<std::string, double> totalJets;  // total selected jets
	int count_multiple_mu = 0; // number of events with multiple muons in jet
	int count_mu_jets = 0;     // number of events with multiple jets with good muons

	
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	std::cout << " Total entries = " << fChain->GetEntries() << std::endl;
	TString tmpfilename = "";

	/////////////////////////////////////
	/////// LOOP OVER ENTRIES ///////////
	/////////////////////////////////////
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (fVerbose) {
			std::cout << "### processing entry: " << jentry << std::endl;
			TFile *tmpfile = fChain->GetFile();
			if (tmpfilename != TString(tmpfile->GetName()) ) {
					tmpfilename = TString(tmpfile->GetName());
					std::cout << " File: "<< tmpfilename << std::endl;
				}
		}
		else if ( jentry%10000 == 0 ) std::cout << "### processing entry: " << jentry << std::endl;
		
		
		if (fVerbose) std::cout << " fill njets and nmuons histograms" << std::endl;
		h1["njets"]->Fill(fS8evt->njets);
		h1["nmuons"]->Fill(fS8evt->nmuons);

			
		std::vector<float> jet_pt_vec = fS8evt->jet_pt;
		int vec_size = jet_pt_vec.size();

		int nmultiple_muons = 0;
		bool event_with_mult_mu = 0;
		int njets_with_lepton = 0;
		//int nopposite_jets = 0;
		
		double ptrel = 0.;
		//bool isTaggability = false;
		//bool isOppositeJetSample = false;
		//bool passGoodMuon = false;
		//bool passJetbTagger = false;
		//bool passOppJetbTagger = false;
		//bool passptrel = false;
		
		
		TLorentzVector p4Jet;
		TLorentzVector p4MuJet;
		TLorentzVector p4OppJet;
		int JetFlavor = -1;
		//int OppJetFlavor = -1;
				
		int ntagtracks = 0;
		////////// Loop over jets ////////////////////////
		for ( int ijet =0; ijet != vec_size; ++ijet) {

			ntagtracks = fS8evt->jet_ntrks[ijet]; //taggability
			h1["jet_ntrks"]->Fill(fS8evt->jet_ntrks[ijet]);
			double jetcorr = fS8evt->jetcorrection[ijet];
						
			//if ( fS8evt->jet_hasLepton[ijet] == 1 && ntagtracks>=1 ) { // no taggability for the moment
			
			//isMuonJetSample = true;
				
			p4Jet.SetPtEtaPhiE(fS8evt->jet_pt[ijet], fS8evt->jet_eta[ijet], fS8evt->jet_phi[ijet], fS8evt->jet_e[ijet] );
			p4Jet = jetcorr * p4Jet;

			//// Jet Selection /////
			if ( p4Jet.Pt() <= 20.|| TMath::Abs( p4Jet.Eta() ) >= 2.5 ) continue;
			// get MC flavor of jet
			JetFlavor = fS8evt->jet_flavour[ijet];

			nmultiple_muons = 0;
			bool passGoodMuon = false;
			if ( fS8evt->jet_hasLepton[ijet] == 1 && fS8evt->njets>0 ) {
			//if ( fS8evt->jet_hasLepton[ijet] == 1 ) {		
			  // get a muon
			  BTagLeptonEvent Muons = fS8evt->lepton[ijet];	
					
					if ( fVerbose ) std::cout << " Muons in jet = " << Muons.pt.size() << std::endl;
					int ith_mu_highest_pt = -1;
					double mu_highest_pt = 0;
				
					for ( size_t imu = 0; imu != Muons.pt.size(); ++imu ) {
						
					  //// Muon Selection /////
					  //if ( ( Muons.trkrechits[imu] >= 8 ) && ( Muons.chi2[imu]/Muons.ndof[imu] <5 ) ) {
							passGoodMuon = true;
							nmultiple_muons++;
							
							if (Muons.pt[imu] > mu_highest_pt ) { mu_highest_pt = Muons.pt[imu]; ith_mu_highest_pt = imu; }
						
							ptrel = Muons.jet_ptrel[imu];
							// double deltaR = Muons.jet_deltaR[imu];
							//std::cout << "ptrel and deltaR " << ptrel << " and " << deltaR <<std::endl;

                            //if ( ptrel > 0.8 ) {
							//	passptrel = true;
							//	}

							if ( fVerbose) std::cout << " muon " << imu << " pt= " << Muons.pt[imu] << " eta= " << Muons.eta[imu] << " chamber hits= " << Muons.SArechits[imu] << " chi2/ndof = " << Muons.chi2[imu]/Muons.ndof[imu] << " IPS= " << Muons.d0sigma[imu] << " mcpdgid= " << Muons.mc_pdgid[imu] << std::endl;
					//}
					}// select only one muon in jet, the one with the highest pt
					
							
					h1["muon_pt"]->Fill(Muons.pt[ith_mu_highest_pt]);
					h1["muon_eta"]->Fill(Muons.eta[ith_mu_highest_pt]);
					h1["jet_deltaR"]->Fill(Muons.jet_deltaR[ith_mu_highest_pt]);
					if ( JetFlavor == 5 ) {
						h1["jet_deltaR_b"]->Fill( Muons.jet_deltaR[ith_mu_highest_pt]);
					}
					if ( JetFlavor == 4 ) {
						h1["jet_deltaR_c"]->Fill( Muons.jet_deltaR[ith_mu_highest_pt]);
					}
					if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor==21 ) {
						h1["jet_deltaR_udsg"]->Fill( Muons.jet_deltaR[ith_mu_highest_pt]);
					}

					h1["jet_pTrel"]->Fill(Muons.jet_ptrel[ith_mu_highest_pt]);
					if ( JetFlavor == 5 ) {
						h1["jet_pTrel_b"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						if ( p4Jet.Pt() > 50 ) {
							h1["jet_pTrel_b_pt50"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( p4Jet.Pt() > 100 ) {
							h1["jet_pTrel_b_pt100"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( p4Jet.Pt() > 150 ) {
							h1["jet_pTrel_b_pt150"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( Muons.pt[ith_mu_highest_pt] > 6 ) {
							h1["jet_pTrel_b_mupt6"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( Muons.pt[ith_mu_highest_pt] > 10 ) {
							h1["jet_pTrel_b_mupt10"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if (Muons.pt[ith_mu_highest_pt] > 20 ) {
							h1["jet_pTrel_b_mupt20"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
					}
					if ( JetFlavor == 4 ) {
						h1["jet_pTrel_c"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						if ( p4Jet.Pt() > 50 ) {
							h1["jet_pTrel_c_pt50"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( p4Jet.Pt() > 100 ) {
							h1["jet_pTrel_c_pt100"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( p4Jet.Pt() > 150 ) {
							h1["jet_pTrel_c_pt150"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( Muons.pt[ith_mu_highest_pt] > 6 ) {
							h1["jet_pTrel_c_mupt6"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( Muons.pt[ith_mu_highest_pt] > 10 ) {
							h1["jet_pTrel_c_mupt10"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if (Muons.pt[ith_mu_highest_pt] > 20 ) {
							h1["jet_pTrel_c_mupt20"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
					}
					if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor==21 ) {
						h1["jet_pTrel_udsg"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						if ( p4Jet.Pt() > 50 ) {
							h1["jet_pTrel_udsg_pt50"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( p4Jet.Pt() > 100 ) {
							h1["jet_pTrel_udsg_pt100"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( p4Jet.Pt() > 150 ) {
							h1["jet_pTrel_udsg_pt150"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( Muons.pt[ith_mu_highest_pt] > 6 ) {
							h1["jet_pTrel_udsg_mupt6"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if ( Muons.pt[ith_mu_highest_pt] > 10 ) {
							h1["jet_pTrel_udsg_mupt10"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
						if (Muons.pt[ith_mu_highest_pt] > 20 ) {
							h1["jet_pTrel_udsg_mupt20"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
						}
					}
					if ( p4Jet.Pt() > 50 ) {
						h1["jet_pTrel_pt50"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
					}
					if ( p4Jet.Pt() > 100 ) {
						h1["jet_pTrel_pt100"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
					}
					if ( p4Jet.Pt() > 150 ) {
						h1["jet_pTrel_pt150"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
					}
					if ( Muons.pt[ith_mu_highest_pt] > 6 ) {
						h1["jet_pTrel_mupt6"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
					}
					if ( Muons.pt[ith_mu_highest_pt] > 10 ) {
						h1["jet_pTrel_mupt10"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
					}
					if (Muons.pt[ith_mu_highest_pt] > 20 ) {
						h1["jet_pTrel_mupt20"]->Fill( Muons.jet_ptrel[ith_mu_highest_pt]);
					}	
					//h1["muon_pt"]->Fill(Muons.pt[ith_mu_highest_pt]);
					//h1["muon_eta"]->Fill(Muons.eta[ith_mu_highest_pt]);
					//std::cout << "filled muon_pt and eta " << std::endl;
					if ( passGoodMuon ) {
					  p4MuJet.SetPtEtaPhiE(fS8evt->jet_pt[ijet]*fS8evt->jetcorrection[ijet],
							       fS8evt->jet_eta[ijet], fS8evt->jet_phi[ijet],
							       fS8evt->jet_e[ijet]*fS8evt->jetcorrection[ijet]);

					  TLorentzVector vmu, vtot;
					  vmu.SetPtEtaPhiE(Muons.pt[ith_mu_highest_pt], Muons.eta[ith_mu_highest_pt], Muons.phi[ith_mu_highest_pt], Muons.e[ith_mu_highest_pt]);

					  vtot = vmu + p4MuJet;
					  ptrel = ( vmu.Px() * vtot.Px()
						    + vmu.Py() * vtot.Py()
						    + vmu.Pz() * vtot.Pz() ) / vtot.P();
					  ptrel = TMath::Sqrt( vmu.P() * vmu.P() - ptrel * ptrel );

					  // check btagging
					  bool isbTaggedJet = false;
					  if ( ftagger == "TrackCounting" ) {
					    if ( flevel == "Tight" ) {
					      if ( fS8evt->btag_TrkCounting_disc3D_3trk[ijet] > fbTaggerMap["Tight"] ) isbTaggedJet = true;
					    } else {
					      if (fVerbose) std::cout << "discriminator= " << fbTaggerMap[flevel] << std::endl;
					      if ( fS8evt->btag_TrkCounting_disc3D_2trk[ijet] > fbTaggerMap[flevel] ) isbTaggedJet = true;

					    }
					  }
					  if ( ftagger == "TrackProbability" ) {
					    if ( fS8evt->btag_JetProb_disc3D[ijet] > fbTaggerMap[flevel] ) isbTaggedJet = true;
					  }

					  FillHistos("n",p4MuJet, ptrel, JetFlavor,isbTaggedJet);
					  
					  ///// find away jet /////
					  bool OtherTaggedJet = false;
					  bool OtherMuonJet = false;
					  for ( int kjet =0; kjet != vec_size; ++kjet) {
					    // skip selected mu+jet
					    if ( ijet == kjet ) continue;

					    TLorentzVector p4OppJet;
					    p4OppJet.SetPtEtaPhiE(fS8evt->jet_pt[kjet], fS8evt->jet_eta[kjet], fS8evt->jet_phi[kjet], fS8evt->jet_e[kjet] );
					    double ojetcorr = fS8evt->jetcorrection[kjet];
					    p4OppJet = ojetcorr * p4OppJet;
					    //int OppJetFlavor = fS8evt->jet_flavour_alg[kjet];
						if ( p4OppJet.Pt() < 20. || TMath::Abs( p4OppJet.Eta() ) >= 2.5 ) continue;
						
					    if ( !OtherTaggedJet ) {
					      
							bool isbTaggedOtherJet = false;
							if ( fAwaytagger == "TrackCounting" ) {
								if ( fAwaylevel == "Tight" ) {
									if ( fS8evt->btag_TrkCounting_disc3D_3trk[kjet] > fbAwayTaggerMap["Tight"] ) isbTaggedOtherJet = true;
								} else {
									if (fVerbose) std::cout << "discriminator= " << fbAwayTaggerMap[fAwaylevel] << std::endl;
									if ( fS8evt->btag_TrkCounting_disc3D_2trk[kjet] > fbAwayTaggerMap[fAwaylevel] ) isbTaggedOtherJet = true;

								}
							}
							if ( fAwaytagger == "TrackProbability" ) {
								if ( fS8evt->btag_JetProb_disc3D[kjet] > fbAwayTaggerMap[fAwaylevel] ) isbTaggedOtherJet = true;
							}
					      
							if (isbTaggedOtherJet) {

								FillHistos("p",p4MuJet, ptrel, JetFlavor, isbTaggedJet);
								OtherTaggedJet = true;
							}
					    }

					    if ( !OtherMuonJet ) {
					      // skip selected mu+jet
					      //if ( ijet == kjet ) continue;
					      // look for another muon-in-jet
					      if ( fS8evt->jet_hasLepton[kjet] == 1 ) {
							  // get a muon
							  BTagLeptonEvent oMuons = fS8evt->lepton[kjet];
							  bool pass2GoodMuon = false;
							  if ( fVerbose ) std::cout << " other Muons in jet = " << oMuons.pt.size() << std::endl;
							  int ith_omu_highest_pt = -1;
							  double omu_highest_pt = 0;
							  for ( size_t imu = 0; imu != oMuons.pt.size(); ++imu ) {

								  //// Muon Selection /////
								  //if ( ( oMuons.trkrechits[imu] >= 8 ) && ( oMuons.chi2[imu]/oMuons.ndof[imu] <5 ) ) {
									  pass2GoodMuon = true;
									  
									  if (oMuons.pt[imu] > omu_highest_pt ) { omu_highest_pt = oMuons.pt[imu]; ith_omu_highest_pt = imu; }
									  
									  //}
							  }
							  // select only one muon in jet, the one with the highest pt
							  if ( pass2GoodMuon ) {
								  OtherMuonJet = true;
								  TLorentzVector p4oMuJet;
								  p4oMuJet.SetPtEtaPhiE(fS8evt->jet_pt[kjet]*fS8evt->jetcorrection[kjet],
														fS8evt->jet_eta[kjet], fS8evt->jet_phi[kjet],
														fS8evt->jet_e[kjet]*fS8evt->jetcorrection[kjet]);
						  
								  //TLorentzVector vmu, vtot;
								  vmu.SetPtEtaPhiE(oMuons.pt[ith_omu_highest_pt], Muons.eta[ith_omu_highest_pt], Muons.phi[ith_omu_highest_pt], Muons.e[ith_omu_highest_pt]);
								  
								  vtot = vmu + p4MuJet;
								  double optrel = ( vmu.Px() * vtot.Px()
													+ vmu.Py() * vtot.Py()
													+ vmu.Pz() * vtot.Pz() ) / vtot.P();
								  optrel = TMath::Sqrt( vmu.P() * vmu.P() - ptrel * ptrel );
								  // do we need ptrel of this muon?
								  
								  // check btagging
								  bool isbTaggedOtherMuJet = false;
								  if ( fAwaytagger == "TrackCounting" ) {
									  if ( fAwaylevel == "Tight" ) {
										  if ( fS8evt->btag_TrkCounting_disc3D_3trk[kjet] > fbAwayTaggerMap["Tight"] ) isbTaggedOtherMuJet = true;
									  } else {
										  if (fVerbose) std::cout << "discriminator= " << fbAwayTaggerMap[fAwaylevel] << std::endl;
										  if ( fS8evt->btag_TrkCounting_disc3D_2trk[kjet] > fbAwayTaggerMap[fAwaylevel] ) isbTaggedOtherMuJet = true;
										  
									  }
								  }
								  if ( fAwaytagger == "TrackProbability" ) {
									  if ( fS8evt->btag_JetProb_disc3D[kjet] > fbAwayTaggerMap[fAwaylevel] ) isbTaggedOtherMuJet = true;
								  }

								  FillHistos("q",p4MuJet, ptrel, JetFlavor, isbTaggedOtherMuJet);
							  }
					      }
					    }// otherMuonJet
					  }// second loop over jets
					  
					  njets_with_lepton++;
					  if (nmultiple_muons>1) event_with_mult_mu = true;
					  
					} // first MuonJet
			}
	       		
			// Fill histograms for all jets
			h1["jet_pt"]->Fill( p4Jet.Pt() );
			h1["jet_eta"]->Fill( p4Jet.Eta() );
			h2["ptVsEta"]->Fill( p4Jet.Pt(), p4Jet.Eta());
			if ( JetFlavor == 5 ) {
			  totalJets["b"]++;
			  h1["jet_pt_b"]->Fill( p4Jet.Pt() );
			  h1["jet_eta_b"]->Fill( p4Jet.Eta() );
			  h2["ptVsEta_b"]->Fill( p4Jet.Pt(), p4Jet.Eta());
			}
			if ( JetFlavor == 4 ) {
			  totalJets["c"]++;
			  h1["jet_pt_c"]->Fill( p4Jet.Pt() );
			  h1["jet_eta_c"]->Fill( p4Jet.Eta() );
			}
			if ( JetFlavor>0 && JetFlavor<4 ) {
			  totalJets["uds"]++;
			  h1["jet_pt_uds"]->Fill( p4Jet.Pt() );
			  h1["jet_eta_uds"]->Fill( p4Jet.Eta() );
			}
			if ( JetFlavor==21 ) {
			  totalJets["g"]++;
			  h1["jet_pt_g"]->Fill( p4Jet.Pt() );
			  h1["jet_eta_g"]->Fill( p4Jet.Eta() );
			}
			bool ataggedjet = false;
			if ( ftagger == "TrackCounting" ) {
			  if ( flevel == "Tight" ) {
			    if ( fS8evt->btag_TrkCounting_disc3D_3trk[ijet] > fbTaggerMap["Tight"] ) ataggedjet = true;
			  } else {
			    if (fVerbose) std::cout << "discriminator= " << fbTaggerMap[flevel] << std::endl;
			    if ( fS8evt->btag_TrkCounting_disc3D_2trk[ijet] > fbTaggerMap[flevel] ) ataggedjet = true;
			    
			  }
			}
			if ( ftagger == "TrackProbability" ) {
			  if ( fS8evt->btag_JetProb_disc3D[ijet] > fbTaggerMap[flevel] ) ataggedjet = true;
			}
			if (ataggedjet) {
			  h1["taggedjet_pt"]->Fill( p4Jet.Pt() );
			  h1["taggedjet_eta"]->Fill( p4Jet.Eta() );
			  h2["taggedjet_ptVsEta"]->Fill( p4Jet.Pt(), p4Jet.Eta());
		  if ( JetFlavor == 5 ) {
			    totalTagged["b"]++;
			    h1["taggedjet_pt_b"]->Fill( p4Jet.Pt() );
				h1["taggedjet_eta_b"]->Fill( p4Jet.Eta() );
			    h2["taggedjet_ptVsEta_b"]->Fill( p4Jet.Pt(), p4Jet.Eta());
			  }
			  if ( JetFlavor == 4 ) {
			    totalTagged["c"]++;
			    h1["taggedjet_pt_c"]->Fill( p4Jet.Pt() );
				h1["taggedjet_eta_c"]->Fill( p4Jet.Eta() );
			  }
			  if ( JetFlavor>0 && JetFlavor<4 ) {
			    totalTagged["uds"]++;
			    h1["taggedjet_pt_uds"]->Fill( p4Jet.Pt() );
				h1["taggedjet_eta_uds"]->Fill( p4Jet.Eta() );
			  }
			  if ( JetFlavor==21 ) {
			    totalTagged["g"]++;
			    h1["taggedjet_pt_g"]->Fill( p4Jet.Pt() );
				h1["taggedjet_eta_g"]->Fill( p4Jet.Eta() );
			  }
			}
			
			// performance plots
			fperformanceTC2.Add(fS8evt->btag_TrkCounting_disc3D_2trk[ijet],JetFlavor);
			fperformanceTC3.Add(fS8evt->btag_TrkCounting_disc3D_3trk[ijet],JetFlavor);
			fperformanceTP.Add(fS8evt->btag_JetProb_disc3D[ijet],JetFlavor);
			
		}//// end loop over jets
		
		if (event_with_mult_mu) count_multiple_mu++;
		if (njets_with_lepton>1) count_mu_jets++;
			
	}// end loop over entries

	//recover b performace plots
	fperformanceTC2.Eval();
	fperformanceTC3.Eval();
	fperformanceTP.Eval();
	
	gTC2_b = new TGraphErrors(fperformanceTC2.GetN(),
											fperformanceTC2.GetArray("b").GetArray(),fperformanceTC2.GetArray("b").GetArray(),
											fperformanceTC2.GetArray("bErr").GetArray(),fperformanceTC2.GetArray("bErr").GetArray());
	
	gTC2_c = new TGraphErrors(fperformanceTC2.GetN(),
											fperformanceTC2.GetArray("b").GetArray(),fperformanceTC2.GetArray("c").GetArray(),
											fperformanceTC2.GetArray("bErr").GetArray(),fperformanceTC2.GetArray("cErr").GetArray());
	
	gTC2_udsg = new TGraphErrors(fperformanceTC2.GetN(),
											fperformanceTC2.GetArray("b").GetArray(),fperformanceTC2.GetArray("udsg").GetArray(),
											fperformanceTC2.GetArray("bErr").GetArray(),fperformanceTC2.GetArray("udsgErr").GetArray());
	
	gTC3_b = new TGraphErrors(fperformanceTC3.GetN(),
											fperformanceTC3.GetArray("b").GetArray(),fperformanceTC3.GetArray("b").GetArray(),
											fperformanceTC3.GetArray("bErr").GetArray(),fperformanceTC3.GetArray("bErr").GetArray());
	
	gTC3_c = new TGraphErrors(fperformanceTC3.GetN(),
											fperformanceTC3.GetArray("b").GetArray(),fperformanceTC3.GetArray("c").GetArray(),
											fperformanceTC3.GetArray("bErr").GetArray(),fperformanceTC3.GetArray("cErr").GetArray());
	
	gTC3_udsg = new TGraphErrors(fperformanceTC3.GetN(),
											fperformanceTC3.GetArray("b").GetArray(),fperformanceTC3.GetArray("udsg").GetArray(),
											fperformanceTC3.GetArray("bErr").GetArray(),fperformanceTC3.GetArray("udsgErr").GetArray());
	
	gTP_b = new TGraphErrors(fperformanceTP.GetN(),
											fperformanceTP.GetArray("b").GetArray(),fperformanceTP.GetArray("b").GetArray(),
											fperformanceTP.GetArray("bErr").GetArray(),fperformanceTP.GetArray("bErr").GetArray());
	
	gTP_c = new TGraphErrors(fperformanceTP.GetN(),
											fperformanceTP.GetArray("b").GetArray(),fperformanceTP.GetArray("c").GetArray(),
											fperformanceTP.GetArray("bErr").GetArray(),fperformanceTP.GetArray("cErr").GetArray());
	
	gTP_udsg = new TGraphErrors(fperformanceTP.GetN(),
											fperformanceTP.GetArray("b").GetArray(),fperformanceTP.GetArray("udsg").GetArray(),
											fperformanceTP.GetArray("bErr").GetArray(),fperformanceTP.GetArray("udsgErr").GetArray());
	
	
	//multiTC2 = new TMultiGraph();
	multiTC3 = new TMultiGraph();
	multiTP = new TMultiGraph();

	// put TC2 and TC3 in one plot
	
	//multiTC2->Add(gTC2_b,"p");
	multiTC3->Add(gTC2_c,"p");
	multiTC3->Add(gTC2_udsg,"p");
	//multiTC3->Add(gTC3_b,"p");
	multiTC3->Add(gTC3_c,"p");
	multiTC3->Add(gTC3_udsg,"p");
	//multiTP->Add(gTP_b,"p");
	multiTP->Add(gTP_c,"p");
	multiTP->Add(gTP_udsg,"p");
	gTC2_b->SetMarkerColor(quark_color["b"]);
	gTC2_c->SetMarkerColor(quark_color["c"]);
	gTC2_udsg->SetMarkerColor(quark_color["uds"]);
	gTC3_b->SetMarkerColor(quark_color["b"]);
	gTC3_c->SetMarkerColor(quark_color["c"]);
	gTC3_udsg->SetMarkerColor(quark_color["uds"]);
	gTP_b->SetMarkerColor(quark_color["b"]);
	gTP_c->SetMarkerColor(quark_color["c"]);
	gTP_udsg->SetMarkerColor(quark_color["uds"]);
	
	gTC2_b->SetMarkerStyle(20);
	gTC2_c->SetMarkerStyle(24);
	gTC2_udsg->SetMarkerStyle(26);
	gTC3_b->SetMarkerStyle(20);
	gTC3_c->SetMarkerStyle(22);
	gTC3_udsg->SetMarkerStyle(23);
	gTP_b->SetMarkerStyle(20);
	gTP_c->SetMarkerStyle(22);
	gTP_udsg->SetMarkerStyle(23);
	/*
	TLegend *legtc2 = new TLegend(0.65,0.65,0.85,0.85,"","NDC");
	legtc2->SetFillColor(10);
	//legtc2->AddEntry(gTC2_b,"b-jet","P");
	legtc2->AddEntry(gTC2_c,"c-jet","P");
	legtc2->AddEntry(gTC2_udsg,"udsg-jet","P");
	legtc2->SetHeader("TrkCounting");
	cv_map["bPerformanceTC2"] = new TCanvas("bPerformanceTC2","bPerformanceTC2",700,700);
	multiTC2->Draw("a");
	multiTC2->GetXaxis()->SetTitle("b-jet efficiency");
	multiTC2->GetYaxis()->SetTitle("non-b jet efficiency");
	multiTC2->SetMinimum(1.e-4);
	multiTC2->SetMaximum(1);
	gPad->SetLogy();
	gPad->SetGrid();
	legtc2->Draw();
	*/
	TLegend *legtc3 = new TLegend(0.65,0.65,0.85,0.85,"","NDC");
	legtc3->SetFillColor(10);
	//legtc3->AddEntry(gTC2_b,"b-jet","P");
	legtc3->AddEntry(gTC2_c,"TC2 c-jet","P");
	legtc3->AddEntry(gTC2_udsg,"TC2 udsg-jet","P");
	legtc3->AddEntry(gTC3_c,"TC3 c-jet","P");
	legtc3->AddEntry(gTC3_udsg,"TC3 udsg-jet","P");
	legtc3->SetHeader("TrkCounting");
	cv_map["bPerformanceTC3"] = new TCanvas("bPerformanceTC3","bPerformanceTC3",700,700);
	multiTC3->Draw("a");
	multiTC3->GetXaxis()->SetTitle("b-jet efficiency");
	multiTC3->GetYaxis()->SetTitle("non-b jet efficiency");
	multiTC3->SetMinimum(1.e-4);
	multiTC3->SetMaximum(1);
	gPad->SetLogy();
	gPad->SetGrid();
	legtc3->Draw();
	TLegend *legtp = new TLegend(0.65,0.65,0.85,0.85,"","NDC");
	legtp->SetFillColor(10);
	//legtp->AddEntry(gTC2_b,"b-jet","P");
	legtp->AddEntry(gTC2_c,"c-jet","P");
	legtp->AddEntry(gTC2_udsg,"udsg-jet","P");
	legtp->SetHeader("TrkProbability");
	cv_map["bPerformanceTP"] = new TCanvas("bPerformanceTP","bPerformanceTP",700,700);
	multiTP->Draw("a");
	multiTP->GetXaxis()->SetTitle("b-jet efficiency");
	multiTP->GetYaxis()->SetTitle("non-b jet efficiency");
	multiTP->SetMinimum(1.e-4);
	multiTP->SetMaximum(1);
	gPad->SetLogy();
	gPad->SetGrid();
	legtp->Draw();
	
	std::cout << std::setfill('#') << std::setw(100) << "#" << std::endl;
	std::cout << std::setfill(' ');
	std::cout << " Total entries = " << fChain->GetEntries() << std::endl;
	std::cout << " Events with multiple muon-in-jet = " << count_mu_jets << std::endl;
	std::cout << " Events with multiple muons in a jet = " << count_multiple_mu << std::endl;
	PrintInfo();
	//std::cout << "no cut" << std::setw(10) << totalJets["uds"] << std::setw(10) << totalJets["c"] << std::setw(10) << totalJets["b"] << std::setw(10) << totalJets["g"] << std::endl;
	//std::cout << "tagged" << std::setw(10) << totalTagged["uds"] << std::setw(10) << totalTagged["c"] << std::setw(10) << totalTagged["b"] << std::setw(10) << totalTagged["g"] << std::endl;
	std::cout << std::fixed;
	std::cout << " <b-efficiency> = " << std::setprecision(2) << std::setw(10) << 100.*totalTagged["b"]/totalJets["b"] <<" \\pm "<< 100.*EffErr(totalTagged["b"],totalJets["b"])<< std::endl;
	std::cout << " <c-mistag>     = " << std::setprecision(2) << std::setw(10) << 100.*totalTagged["c"]/totalJets["c"] <<" \\pm "<< 100.*EffErr(totalTagged["c"],totalJets["c"])<< std::endl;
	std::cout << " <uds-mistag>   = " << std::setprecision(2) << std::setw(10) << 100.*totalTagged["uds"]/totalJets["uds"] <<" \\pm "<< 100.*EffErr(totalTagged["uds"],totalJets["uds"])<< std::endl;
	std::cout << " <g-mistag>     = " << std::setprecision(2) << std::setw(10) << 100.*totalTagged["g"]/totalJets["g"] <<" \\pm "<< 100.*EffErr(totalTagged["g"],totalJets["g"])<< std::endl;
	std::cout << " <udsg-mistag>  = " << std::setprecision(2) << std::setw(10) << 100.*(totalTagged["g"]+totalTagged["uds"])/(totalJets["g"]+totalJets["uds"]) <<" \\pm "<< 100.*EffErr(totalTagged["g"]+totalTagged["uds"],totalJets["g"]+totalJets["uds"])<< std::endl;
	std::cout << std::setfill('#') << std::setw(100) << "#" << std::endl;
	std::cout << std::setfill(' ');


	//______________________________________________________

	cv_map["jet_pt"] = new TCanvas("jet_pt","jet_pt",700,700);
	h1["jet_pt"]->SetXTitle("jet p_{T} [GeV/c]");
	h1["jet_pt"]->Draw();

	//______________________________________________________

	cv_map["muon_pt"] = new TCanvas("muon_pt","muon_pt",700,700);
	h1["muon_pt"]->SetXTitle("muon p_{T} [GeV/c]");
	h1["muon_pt"]->Draw();

		
	//______________________________________________________
	cv_map["deltaphi"] = new TCanvas("deltaphi","deltaphi",700,700);
	h1["jet_deltaphi"]->SetXTitle( h1["jet_deltaphi"]->GetTitle() );
	h1["jet_deltaphi"]->Draw();
	h1["jet_deltaphi_b"]->Draw("same");
		
	//______________________________________________________

	h1["eff_pTrel"] = new TH1D("eff_pTrel","eff_pTrel",nptbins,jetptbins);
	h1["eff_pTrel_b"] = new TH1D("eff_pTrel_b","eff_pTrel_b",nptbins,jetptbins);
	h1["eff_pTrel_cl"] = new TH1D("eff_pTrel_cl","eff_pTrel_cl",nptbins,jetptbins);

	h1["eff_TaggedJet"] = new TH1D("eff_TaggedJet","eff_TaggedJet",nptbins,jetptbins);
	h1["eff_TaggedJet_b"] = new TH1D("eff_TaggedJet_b","eff_TaggedJet_b",nptbins,jetptbins);
	h1["eff_TaggedJet_cl"] = new TH1D("eff_TaggedJet_cl","eff_TaggedJet_cl",nptbins,jetptbins);

	h1["eff_TaggedJet_c"] = new TH1D("eff_TaggedJet_c","eff_TaggedJet_c",nptbins,jetptbins);
	h1["eff_TaggedJet_uds"] = new TH1D("eff_TaggedJet_uds","eff_TaggedJet_uds",nptbins,jetptbins);
	
	h1["eff_OppTaggedJet"] = new TH1D("eff_OppTaggedJet","eff_OppTaggedJet",nptbins,jetptbins);
	h1["eff_OppTaggedJet_b"] = new TH1D("eff_OppTaggedJet_b","eff_OppTaggedJet_b",nptbins,jetptbins);
	h1["eff_OppTaggedJet_cl"] = new TH1D("eff_OppTaggedJet_cl","eff_OppTaggedJet_cl",nptbins,jetptbins);

	h1["eff_TaggedBothJets"] = new TH1D("eff_TaggedBothJets","eff_TaggedBothJets",nptbins,jetptbins);
	h1["eff_TaggedBothJets_b"] = new TH1D("eff_TaggedBothJets_b","eff_TaggedBothJets_b",nptbins,jetptbins);
	h1["eff_TaggedBothJets_cl"] = new TH1D("eff_TaggedBothJets_cl","eff_TaggedBothJets_cl",nptbins,jetptbins);

	h1["eff_pTrel_TaggedJet"] = new TH1D("eff_pTrel_TaggedJet","eff_pTrel_TaggedJet",nptbins,jetptbins);
	h1["eff_pTrel_TaggedJet_b"] = new TH1D("eff_pTrel_TaggedJet_b","eff_pTrel_TaggedJet_b",nptbins,jetptbins);
	h1["eff_pTrel_TaggedJet_cl"] = new TH1D("eff_pTrel_TaggedJet_cl","eff_pTrel_TaggedJet_cl",nptbins,jetptbins);

	h1["eff_pTrel_OppTaggedJet"] = new TH1D("eff_pTrel_OppTaggedJet","eff_pTrel_OppTaggedJet",nptbins,jetptbins);
	h1["eff_pTrel_OppTaggedJet_b"] = new TH1D("eff_pTrel_OppTaggedJet_b","eff_pTrel_OppTaggedJet_b",nptbins,jetptbins);
	h1["eff_pTrel_OppTaggedJet_cl"] = new TH1D("eff_pTrel_OppTaggedJet_cl","eff_pTrel_OppTaggedJet_cl",nptbins,jetptbins);

	h1["eff_pTrel_TaggedBothJets"] = new TH1D("eff_pTrel_TaggedBothJets","eff_pTrel_TaggedBothJets",nptbins,jetptbins);
	h1["eff_pTrel_TaggedBothJets_b"] = new TH1D("eff_pTrel_TaggedBothJets_b","eff_pTrel_TaggedBothJets_b",nptbins,jetptbins);
	h1["eff_pTrel_TaggedBothJets_cl"] = new TH1D("eff_pTrel_TaggedBothJets_cl","eff_pTrel_TaggedBothJets_cl",nptbins,jetptbins);

	h1["eff_mu_taggedaway"]    = new TH1D("eff_mu_taggedaway",   "eff_mu_taggedaway",nptbins,jetptbins);
	h1["eff_mu_taggedaway_b"]  = new TH1D("eff_mu_taggedaway_b", "eff_mu_taggedaway_b",nptbins,jetptbins);
	h1["eff_mu_taggedaway_cl"] = new TH1D("eff_mu_taggedaway_cl","eff_mu_taggedaway_cl",nptbins,jetptbins);
	// Booking for eta dependencyjetetabins

        h1["eff_pTrel_eta"] = new TH1D("eff_pTrel_eta","eff_pTrel_eta",netabins,jetetabins);
        h1["eff_pTrel_eta_b"] = new TH1D("eff_pTrel_eta_b","eff_pTrel_eta_b",netabins,jetetabins);
        h1["eff_pTrel_eta_cl"] = new TH1D("eff_pTrel_eta_cl","eff_pTrel_eta_cl",netabins,jetetabins);

        h1["eff_TaggedJet_eta"] = new TH1D("eff_TaggedJet_eta","eff_TaggedJet_eta",netabins,jetetabins);
        h1["eff_TaggedJet_eta_b"] = new TH1D("eff_TaggedJet_eta_b","eff_TaggedJet_eta_b",netabins,jetetabins);
        h1["eff_TaggedJet_eta_cl"] = new TH1D("eff_TaggedJet_eta_cl","eff_TaggedJet_eta_cl",netabins,jetetabins);

        h1["eff_TaggedJet_eta_c"] = new TH1D("eff_TaggedJet_eta_c","eff_TaggedJet_eta_c",netabins,jetetabins);
        h1["eff_TaggedJet_eta_uds"] = new TH1D("eff_TaggedJet_eta_uds","eff_TaggedJet_eta_uds",netabins,jetetabins);

        h1["eff_OppTaggedJet_eta"] = new TH1D("eff_OppTaggedJet_eta","eff_OppTaggedJet_eta",netabins,jetetabins);
        h1["eff_OppTaggedJet_eta_b"] = new TH1D("eff_OppTaggedJet_eta_b","eff_OppTaggedJet_eta_b",netabins,jetetabins);
        h1["eff_OppTaggedJet_eta_cl"] = new TH1D("eff_OppTaggedJet_eta_cl","eff_OppTaggedJet_eta_cl",netabins,jetetabins);

        h1["eff_TaggedBothJets_eta"] = new TH1D("eff_TaggedBothJets_eta","eff_TaggedBothJets_eta",netabins,jetetabins);
        h1["eff_TaggedBothJets_eta_b"] = new TH1D("eff_TaggedBothJets_eta_b","eff_TaggedBothJets_eta_b",netabins,jetetabins);
        h1["eff_TaggedBothJets_eta_cl"] = new TH1D("eff_TaggedBothJets_eta_cl","eff_TaggedBothJets_eta_cl",netabins,jetetabins);

        h1["eff_pTrel_TaggedJet_eta"] = new TH1D("eff_pTrel_TaggedJet_eta","eff_pTrel_TaggedJet_eta",netabins,jetetabins);
        h1["eff_pTrel_TaggedJet_eta_b"] = new TH1D("eff_pTrel_TaggedJet_eta_b","eff_pTrel_TaggedJet_eta_b",netabins,jetetabins);
        h1["eff_pTrel_TaggedJet_eta_cl"] = new TH1D("eff_pTrel_TaggedJet_eta_cl","eff_pTrel_TaggedJet_eta_cl",netabins,jetetabins);

        h1["eff_pTrel_OppTaggedJet_eta"] = new TH1D("eff_pTrel_OppTaggedJet_eta","eff_pTrel_OppTaggedJet_eta",netabins,jetetabins);
        h1["eff_pTrel_OppTaggedJet_eta_b"] = new TH1D("eff_pTrel_OppTaggedJet_eta_b","eff_pTrel_OppTaggedJet_eta_b",netabins,jetetabins);
        h1["eff_pTrel_OppTaggedJet_eta_cl"] = new TH1D("eff_pTrel_OppTaggedJet_eta_cl","eff_pTrel_OppTaggedJet_eta_cl",netabins,jetetabins);

        h1["eff_pTrel_TaggedBothJets_eta"] = new TH1D("eff_pTrel_TaggedBothJets_eta","eff_pTrel_TaggedBothJets_eta",netabins,jetetabins);
        h1["eff_pTrel_TaggedBothJets_eta_b"] = new TH1D("eff_pTrel_TaggedBothJets_eta_b","eff_pTrel_TaggedBothJets_eta_b",netabins,jetetabins);
        h1["eff_pTrel_TaggedBothJets_eta_cl"] = new TH1D("eff_pTrel_TaggedBothJets_eta_cl","eff_pTrel_TaggedBothJets_eta_cl",netabins,jetetabins);

//	h1["eff_mu_taggedaway_eta"]    = new TH1D("eff_mu_taggedaway_eta",   "eff_mu_taggedaway_eta",netabins,jetetabins);
//	h1["eff_mu_taggedaway_eta_b"]  = new TH1D("eff_mu_taggedaway_eta_b", "eff_mu_taggedaway_eta_b",netabins,jetetabins);
//	h1["eff_mu_taggedawa_etay_cl"] = new TH1D("eff_mu_taggedaway_eta_cl","eff_mu_taggedaway_eta_cl",netabins,jetetabins);
	//TH1D *halljets = (TH1D*) h1["eff_pTrel"]->Clone("halljets");
	TH1D *halljets = h2["npT"]->ProjectionX("halljets", -1 , -1,"e");
	TH1D *halljets_b = h2["b_npT"]->ProjectionX("halljets_b", -1 , -1,"e");
	TH1D *halljets_cl = h2["cl_npT"]->ProjectionX("halljets_cl", -1 , -1,"e");
	TH1D *halljets_c = h2["c_npT"]->ProjectionX("halljets_c", -1 , -1,"e");
	TH1D *halljets_uds = h2["uds_npT"]->ProjectionX("halljets_uds", -1 , -1,"e");

	h1["eff_pTrel"]->Divide( h2["npT"]->ProjectionX("halljets_ptrel", 9 , -1,"e") , halljets ,1.,1.,"B");
	h1["eff_pTrel_b"]->Divide( h2["b_npT"]->ProjectionX("b_halljets_ptrel", 9 , -1,"e") , halljets_b ,1.,1.,"B");
	h1["eff_pTrel_cl"]->Divide( h2["cl_npT"]->ProjectionX("cl_halljets_ptrel", 9 , -1,"e") , halljets_cl ,1.,1.,"B");
	
	h1["eff_TaggedJet"]->Divide( h2["ncmbpT"]->ProjectionX("halljets_tagged", -1 , -1,"e") , halljets ,1.,1.,"B");
	h1["eff_TaggedJet_b"]->Divide( h2["b_ncmbpT"]->ProjectionX("b_halljets_tagged", -1 , -1,"e") , halljets_b ,1.,1.,"B");
	h1["eff_TaggedJet_cl"]->Divide( h2["cl_ncmbpT"]->ProjectionX("cl_halljets_tagged", -1 , -1,"e") , halljets_cl ,1.,1.,"B");


	h1["eff_TaggedJet_c"]->Divide( h2["c_ncmbpT"]->ProjectionX("c_halljets_tagged", -1 , -1,"e") , halljets_c ,1.,1.,"B");
	h1["eff_TaggedJet_uds"]->Divide( h2["uds_ncmbpT"]->ProjectionX("uds_halljets_tagged", -1 , -1,"e") , halljets_uds ,1.,1.,"B");
	
	
	h1["eff_pTrel_TaggedJet"]->Divide( h2["ncmbpT"]->ProjectionX("halljets_ptreltagged", 9 , -1,"e") , halljets ,1.,1.,"B");
	h1["eff_pTrel_TaggedJet_b"]->Divide( h2["b_ncmbpT"]->ProjectionX("b_halljets_ptreltagged", 9 , -1,"e") , halljets_b ,1.,1.,"B");
	h1["eff_pTrel_TaggedJet_cl"]->Divide( h2["cl_ncmbpT"]->ProjectionX("cl_halljets_ptrelltagged", 9 , -1,"e") , halljets_cl ,1.,1.,"B");

	// Eta dependency

        TH1D *halljets_eta = h2["nEta"]->ProjectionX("halljets_eta", -1 , -1,"e");
        TH1D *halljets_eta_b = h2["b_nEta"]->ProjectionX("halljets_eta_b", -1 , -1,"e");
        TH1D *halljets_eta_cl = h2["cl_nEta"]->ProjectionX("halljets_eta_cl", -1 , -1,"e");
        // TH1D *halljets_eta_c = h2["c_nEta"]->ProjectionX("halljetseta_c", -1 , -1,"e");
        // TH1D *halljets_eta_uds = h2["uds_nEta"]->ProjectionX("halljetseta_uds", -1 , -1,"e");

        h1["eff_pTrel_eta"]->Divide( h2["nEta"]->ProjectionX("halljetseta_ptrel", 9 , -1,"e") , halljets_eta ,1.,1.,"B");
        h1["eff_pTrel_eta_b"]->Divide( h2["b_nEta"]->ProjectionX("b_halljetseta_ptrel", 9 , -1,"e") , halljets_eta_b ,1.,1.,"B");
        h1["eff_pTrel_eta_cl"]->Divide( h2["cl_nEta"]->ProjectionX("cl_halljetseta_ptrel", 9 , -1,"e") , halljets_eta_cl ,1.,1.,"B");

        h1["eff_TaggedJet_eta"]->Divide( h2["ncmbEta"]->ProjectionX("halljetseta_tagged", -1 , -1,"e") , halljets_eta ,1.,1.,"B");
        h1["eff_TaggedJet_eta_b"]->Divide( h2["b_ncmbEta"]->ProjectionX("b_halljetseta_tagged", -1 , -1,"e") , halljets_eta_b ,1.,1.,"B");
        h1["eff_TaggedJet_eta_cl"]->Divide( h2["cl_ncmbEta"]->ProjectionX("cl_halljetseta_tagged", -1 , -1,"e") , halljets_eta_cl ,1.,1.,"B");

        // h1["eff_TaggedJet_eta_c"]->Divide( h2["c_ncmbEta"]->ProjectionX("c_halljetseta_tagged", -1 , -1,"e") , halljets_eta_c ,1.,1.,"B");
        // h1["eff_TaggedJet_eta_uds"]->Divide( h2["uds_ncmbEta"]->ProjectionX("uds_halljetseta_tagged", -1 , -1,"e") , halljets_eta_uds ,1.,1.,"B");

        h1["eff_pTrel_TaggedJet_eta"]->Divide( h2["ncmbEta"]->ProjectionX("halljetseta_ptreltagged", 9 , -1,"e") , halljets_eta ,1.,1.,"B");
        h1["eff_pTrel_TaggedJet_eta_b"]->Divide( h2["b_ncmbEta"]->ProjectionX("b_halljetseta_ptreltagged", 9 , -1,"e") , halljets_eta_b ,1.,1.,"B");
        h1["eff_pTrel_TaggedJet_eta_cl"]->Divide( h2["cl_ncmbEta"]->ProjectionX("cl_halljetseta_ptrelltagged", 9 , -1,"e") , halljets_eta_cl ,1.,1.,"B");

        TH1D *halloppjets_eta = h2["pEta"]->ProjectionX("halloppjets_eta", -1 , -1,"e");
        TH1D *halloppjets_eta_b = h2["b_pEta"]->ProjectionX("halloppjets_eta_b", -1 , -1,"e");
        TH1D *halloppjets_eta_cl = h2["cl_pEta"]->ProjectionX("halloppjets_eta_cl", -1 , -1,"e");

        h1["eff_TaggedBothJets_eta"]->Divide( h2["pcmbEta"]->ProjectionX("halloppjetseta_tagged", -1 , -1,"e") , halloppjets_eta ,1.,1.,"B");
        h1["eff_TaggedBothJets_eta_b"]->Divide( h2["b_pcmbEta"]->ProjectionX("b_halloppjetseta_tagged", -1 , -1,"e") , halloppjets_eta_b ,1.,1.,"B");
        h1["eff_TaggedBothJets_eta_cl"]->Divide( h2["cl_pcmbEta"]->ProjectionX("cl_halloppjetseta_tagged", -1 , -1,"e") , halloppjets_eta_cl ,1.,1.,"B");


	TH1D *halloppjets = h2["ppT"]->ProjectionX("halloppjets", -1 , -1,"e");
	TH1D *halloppjets_b = h2["b_ppT"]->ProjectionX("halloppjets_b", -1 , -1,"e");
	TH1D *halloppjets_cl = h2["cl_ppT"]->ProjectionX("halloppjets_cl", -1 , -1,"e");

	h1["eff_TaggedBothJets"]->Divide( h2["pcmbpT"]->ProjectionX("halloppjets_tagged", -1 , -1,"e") , halloppjets ,1.,1.,"B");
	h1["eff_TaggedBothJets_b"]->Divide( h2["b_pcmbpT"]->ProjectionX("b_halloppjets_tagged", -1 , -1,"e") , halloppjets_b ,1.,1.,"B");
	h1["eff_TaggedBothJets_cl"]->Divide( h2["cl_pcmbpT"]->ProjectionX("cl_halloppjets_tagged", -1 , -1,"e") , halloppjets_cl ,1.,1.,"B");
	h1["eff_mu_taggedaway"]->Divide( h2["ppT"]->ProjectionX("halloppjets_ptrel", 9 , -1,"e"), halloppjets, 1.,1.,"B");
	h1["eff_mu_taggedaway_b"]->Divide( h2["b_ppT"]->ProjectionX("b_halloppjets_ptrel", 9 , -1,"e"), halloppjets_b, 1.,1.,"B");
	h1["eff_mu_taggedaway_cl"]->Divide(h2["cl_ppT"]->ProjectionX("cl_halloppjets_ptrel", 9, -1,"e"), halloppjets_cl, 1.,1.,"B");

//i am here
//	h1["eff_mu_taggedaway_eta_b"]->Divide( h2["b_ppT"]->ProjectionX("b_halloppjetseta_ptrel", 9 , -1,"e"), halloppjets_eta_b, 1.,1.,"B");
//	h1["eff_mu_taggedaway_eta_cl"]->Divide(h2["cl_ppT"]->ProjectionX("cl_halloppjetseta_ptrel", 9, -1,"e"), halloppjets_eta_cl, 1.,1.,"B");
		
	// alpha
	h1["alpha"]->Divide( h1["eff_TaggedBothJets_cl"], h1["eff_TaggedJet_cl"]);
	//h1["alpha"]->Divide( halloppjets_cl, h1["eff_TaggedJet_cl"]);
	// beta
	h1["beta"]->Divide( h1["eff_TaggedBothJets_b"], h1["eff_TaggedJet_b"]);
	//h1["beta"]->Divide( halloppjets_b, h1["eff_TaggedJet_b"]);

	
	// kappa_b
	h1["kappa_b"]->Divide( h1["eff_pTrel_TaggedJet_b"], h1["eff_pTrel_b"]  );
	h1["kappa_b"]->Divide( h1["eff_TaggedJet_b"] );
	// kappa_cl
	h1["kappa_cl"]->Divide( h1["eff_pTrel_TaggedJet_cl"], h1["eff_pTrel_cl"] );
	h1["kappa_cl"]->Divide(  h1["eff_TaggedJet_cl"] );
	
	h1["delta"]->Divide( h1["eff_mu_taggedaway_b"], h1["eff_pTrel_b"] );
	h1["gamma"]->Divide( h1["eff_mu_taggedaway_cl"], h1["eff_pTrel_cl"] );
	        // eta dependency for the correlations

        // alpha
        h1["alpha_eta"]->Divide( h1["eff_TaggedBothJets_eta_cl"], h1["eff_TaggedJet_eta_cl"]);
        //h1["alpha"]->Divide( halloppjets_cl, h1["eff_TaggedJet_cl"]);
        // beta
        h1["beta_eta"]->Divide( h1["eff_TaggedBothJets_eta_b"], h1["eff_TaggedJet_eta_b"]);
        //h1["beta"]->Divide( halloppjets_b, h1["eff_TaggedJet_b"]);

        // kappa_b
        h1["kappa_eta_b"]->Divide( h1["eff_pTrel_TaggedJet_eta_b"], h1["eff_pTrel_eta_b"]  );
        h1["kappa_eta_b"]->Divide( h1["eff_TaggedJet_eta_b"] );
        // kappa_cl
        h1["kappa_eta_cl"]->Divide( h1["eff_pTrel_TaggedJet_eta_cl"], h1["eff_pTrel_eta_cl"] );
        h1["kappa_eta_cl"]->Divide(  h1["eff_TaggedJet_eta_cl"] );

//		h1["delta_eta"]->Divide( h1["eff_mu_taggedaway_eta_b"], h1["eff_pTrel_eta_b"] );
//		h1["gamma_eta"]->Divide( h1["eff_mu_taggedaway_eta_cl"], h1["eff_pTrel_eta_cl"] );
	
	//______________________________________________________

	cv_map["eff_pTrel"] = new TCanvas("eff_pTrel","eff_pTrel",700,700);
	h1["eff_pTrel"]->SetYTitle("Efficiency");
	h1["eff_pTrel"]->SetMarkerStyle(8);
	h1["eff_pTrel"]->SetMarkerSize(1.5);
	h1["eff_pTrel"]->Draw("PE1");
	gPad->SetGrid();
	//______________________________________________________
	
	cv_map["eff_TaggedJet"] = new TCanvas("eff_TaggedJet","eff_TaggedJet",700,700);
	h1["eff_TaggedJet"]->SetYTitle("Efficiency");
	h1["eff_TaggedJet"]->SetMarkerStyle(8);
	h1["eff_TaggedJet"]->SetMarkerSize(1.5);
	h1["eff_TaggedJet"]->Draw("PE1");
	gPad->SetGrid();
	//______________________________________________________
	
	cv_map["eff_pTrel_TaggedJet"] = new TCanvas("eff_pTrel_TaggedJet","eff_pTrel_TaggedJet",700,700);
	h1["eff_TaggedBothJets"]->SetYTitle("Efficiency");
	h1["eff_TaggedBothJets"]->SetMarkerStyle(8);
	h1["eff_TaggedBothJets"]->SetMarkerSize(1.5);
	h1["eff_TaggedBothJets"]->Draw("PE1");
	gPad->SetGrid();
	//______________________________________________________

	double maxYaxis = 1.1;
	double minYaxis = 0.1;
	cv_map["kappa_b"] = new TCanvas("kappa_b","kappa_b",700,700);
	h1["kappa_b"]->SetMarkerStyle(8);
	h1["kappa_b"]->SetMarkerSize(1.5);
	h1["kappa_b"]->SetLineColor(quark_color["b"]);
	h1["kappa_b"]->SetMarkerColor(quark_color["b"]);
	h1["kappa_b"]->SetXTitle("jet p_{T} [GeV/c]");
	h1["kappa_b"]->SetMaximum(maxYaxis);
	h1["kappa_b"]->SetMinimum(minYaxis);
	h1["kappa_b"]->Draw("PE1");
	std::cout << " kappa_b Fit" << std::endl;
	h1["kappa_b"]->Fit("pol0","0");
	TF1 *f1_kb = h1["kappa_b"]->GetFunction("pol0");
	f1_kb->SetLineColor(quark_color["b"]);
	f1_kb->Draw("same");
	h1["eff_TaggedJet_b"]->SetMarkerStyle(23);
	h1["eff_TaggedJet_b"]->Draw("PE1 same");
	h1["eff_pTrel_b"]->SetMarkerStyle(21);
	h1["eff_pTrel_b"]->Draw("PE1 same");
	h1["eff_pTrel_TaggedJet_b"]->SetMarkerStyle(26);
	h1["eff_pTrel_TaggedJet_b"]->Draw("PE1 same");
	TLegend *legkb = new TLegend(0.6,0.65,0.85,0.85,"","NDC");
	legkb->SetFillColor(10);
	legkb->AddEntry(h1["kappa_b"],"#kappa_{b}","P");
	TString strTag;
	if ( ftagger == "TrackCounting")  strTag = "TC";
	if ( ftagger == "TrackProbability") strTag = "TP";
	legkb->AddEntry(h1["eff_pTrel_TaggedJet_b"],"#epsilon(p_{T}^{Rel}#wedge"+strTag+"-"+flevel+")_{b}^{#mu-jet}","P");
	legkb->AddEntry(h1["eff_pTrel_b"],"#epsilon(p_{T}^{Rel})_{b}^{muon-jet}","P");
	legkb->AddEntry(h1["eff_TaggedJet_b"],"#epsilon("+strTag+"-"+flevel+")_{b}^{#mu-jet}","P");
	legkb->Draw();
	gPad->SetGrid();
    //______________________________________________________

	cv_map["kappa_cl"] = new TCanvas("kappa_cl","kappa_cl",700,700);
	h1["kappa_cl"]->SetMarkerStyle(8);
	h1["kappa_cl"]->SetMarkerSize(1.5);
	h1["kappa_cl"]->SetLineColor(quark_color["c"]);
	h1["kappa_cl"]->SetMarkerColor(quark_color["c"]);
	h1["kappa_cl"]->SetXTitle("jet p_{T} [GeV/c]");
	h1["kappa_cl"]->SetMaximum(maxYaxis);
	h1["kappa_cl"]->SetMinimum(minYaxis);
	h1["kappa_cl"]->Draw("PE1");
	std::cout << " kappa_cl Fit" << std::endl;
	h1["kappa_cl"]->Fit("pol0","0");
	TF1 *f1_kcl = h1["kappa_cl"]->GetFunction("pol0");
	f1_kcl->SetLineColor(quark_color["c"]);
	f1_kcl->Draw("same");
	h1["eff_TaggedJet_cl"]->SetMarkerStyle(23);
	h1["eff_TaggedJet_cl"]->Draw("PE1 same");
	h1["eff_pTrel_cl"]->SetMarkerStyle(21);
	h1["eff_pTrel_cl"]->Draw("PE1 same");
	h1["eff_pTrel_TaggedJet_cl"]->SetMarkerStyle(26);
        h1["eff_pTrel_TaggedJet_cl"]->Draw("PE1 same");
	TLegend *legkc = new TLegend(0.6,0.65,0.85,0.85,"","NDC");
	legkc->SetFillColor(10);
	legkc->AddEntry(h1["kappa_cl"],"#kappa_{cl}","P");
	legkc->AddEntry(h1["eff_pTrel_TaggedJet_cl"],"#epsilon(p_{T}^{Rel}#wedge"+strTag+"-"+flevel+")_{cl}^{#mu-jet}","P");
	legkc->AddEntry(h1["eff_pTrel_cl"],"#epsilon(p_{T}^{Rel})_{cl}^{#mu-jet}","P");
	legkc->AddEntry(h1["eff_TaggedJet_cl"],"#epsilon("+strTag+"-"+flevel+")_{cl}^{#mu-jet}","P");
	legkc->Draw();
	gPad->SetGrid();
    //______________________________________________________
	
	cv_map["alpha"] = new TCanvas("alpha","alpha",700,700);
	h1["alpha"]->SetMarkerStyle(8);
	h1["alpha"]->SetMarkerSize(1.5);
	h1["alpha"]->SetLineColor(quark_color["c"]);
	h1["alpha"]->SetMarkerColor(quark_color["c"]);
	h1["alpha"]->SetXTitle("jet p_{T} [GeV/c]");
	h1["alpha"]->SetMaximum(maxYaxis);
	h1["alpha"]->SetMinimum(minYaxis);
	h1["alpha"]->Draw("PE1");
	std::cout << " alpha Fit" << std::endl;
	h1["alpha"]->Fit("pol0","0");
	TF1 *f1_alpha = h1["alpha"]->GetFunction("pol0");
	f1_alpha->SetLineColor(quark_color["c"]);
	f1_alpha->Draw("same");
	h1["eff_TaggedBothJets_cl"]->SetMarkerStyle(26);
	h1["eff_TaggedBothJets_cl"]->Draw("PE1 same");
	h1["eff_TaggedJet_cl"]->SetMarkerStyle(23);
	h1["eff_TaggedJet_cl"]->Draw("PE1 same");
	TLegend *lega = new TLegend(0.6,0.65,0.85,0.85,"","NDC");
        lega->SetFillColor(10);
        lega->AddEntry(h1["alpha"],"#alpha","P");
        lega->AddEntry(h1["eff_TaggedJet_cl"],"#epsilon("+strTag+"-"+flevel+")_{cl}^{#mu-jet}","P");
        lega->AddEntry(h1["eff_TaggedBothJets_cl"],"#epsilon("+strTag+"-"+flevel+")_{cl}^{#mu-jet-away-jet}","P");
        lega->Draw();
	gPad->SetGrid();

    //______________________________________________________
	
	cv_map["beta"] = new TCanvas("beta","beta",700,700);
	h1["beta"]->SetMarkerStyle(8);
	h1["beta"]->SetMarkerSize(1.5);
	h1["beta"]->SetLineColor(quark_color["b"]);
	h1["beta"]->SetMarkerColor(quark_color["b"]);
	h1["beta"]->SetXTitle("jet p_{T} [GeV/c]");
	h1["beta"]->SetMaximum(maxYaxis);
	h1["beta"]->SetMinimum(minYaxis);
	h1["beta"]->Draw("PE1");
	std::cout << " beta Fit" << std::endl;
	h1["beta"]->Fit("pol0","0");
	TF1 *f1_beta = h1["beta"]->GetFunction("pol0");
	f1_beta->SetLineColor(quark_color["b"]);
	f1_beta->Draw("same");
	h1["eff_TaggedBothJets_b"]->SetMarkerStyle(26);
	h1["eff_TaggedBothJets_b"]->Draw("PE1 same");
	h1["eff_TaggedJet_b"]->SetMarkerStyle(23);
	h1["eff_TaggedJet_b"]->Draw("PE1 same");
	TLegend *legb = new TLegend(0.6,0.65,0.85,0.85,"","NDC");
        legb->SetFillColor(10);
        legb->AddEntry(h1["beta"],"#beta","P");
        legb->AddEntry(h1["eff_TaggedJet_b"],"#epsilon("+strTag+"-"+flevel+")_{b}^{#mu-jet}","P");
        legb->AddEntry(h1["eff_TaggedBothJets_b"],"#epsilon("+strTag+"-"+flevel+")_{b}^{#mu-jet-away-jet}","P");
        legb->Draw();
	gPad->SetGrid();
	
}

void S8Plotter::FillHistos(std::string type, TLorentzVector p4Jet, double ptrel, int flavor, bool tagged) {

  if ( ( p4Jet.Pt() >= fMinPt ) && ( TMath::Abs(p4Jet.Eta()) >= fMinEta ) && (TMath::Abs(p4Jet.Eta()) < fMaxEta ) ) {
			
		
	h2[type+"pT"]->Fill( p4Jet.Pt() , ptrel );
	h2[type+"Eta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
	if ( flavor == 5 ) {
		h2["b_"+type+"pT"]->Fill( p4Jet.Pt() , ptrel );
		h2["b_"+type+"Eta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
	}
	if ( (flavor>0 && flavor<5) || flavor==21 ) {
		h2["cl_"+type+"pT"]->Fill( p4Jet.Pt() , ptrel);
		h2["cl_"+type+"Eta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
	}
	if ( flavor == 4 ) { 
		//h2["c_"+type+"pT"]->Fill( p4Jet.Pt() , ptrel ); 
	}
	if ( flavor>0 && flavor<4 ) { 
		//h2["uds_"+type+"pT"]->Fill( p4Jet.Pt() , ptrel ); 
	}
	if ( flavor==21 ) {
		//h2["g_"+type+"pT"]->Fill( p4Jet.Pt() , ptrel );
	}
	
	if ( tagged ) {
		h2[type+"cmbpT"]->Fill( p4Jet.Pt() , ptrel );
		h2[type+"cmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		if ( flavor == 5 ) {
			h2["b_"+type+"cmbpT"]->Fill( p4Jet.Pt() , ptrel );
			h2["b_"+type+"cmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		}
		if ( (flavor>0 && flavor<5) || flavor==21 ) {
			h2["cl_"+type+"cmbpT"]->Fill( p4Jet.Pt() , ptrel);
			h2["cl_"+type+"cmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		}
		if ( flavor == 4 ) {
			//h2["c_"+type+"cmbpT"]->Fill( p4Jet.Pt() , ptrel );
		}
		if ( flavor>0 && flavor<4 ) {
			//h2["uds_"+type+"cmbpT"]->Fill( p4Jet.Pt() , ptrel );
		}
		if ( flavor==21 ) {
			//h2["g_"+type+"cmbpT"]->Fill( p4Jet.Pt() , ptrel );
		}
		
	}

	}

}


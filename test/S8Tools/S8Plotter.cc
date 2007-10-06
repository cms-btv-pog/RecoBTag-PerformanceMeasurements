#define S8Plotter_cxx
#include "S8Plotter.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>
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
	quark_label.push_back("_b");
	quark_label.push_back("_c");
	quark_label.push_back("_udsg");
	
	quark_color["b"] = 2;
	quark_color["c"] = 4;
	quark_color["udsg"] = 6;

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

	hname = "jet_pt"; 
	htitle = "jet p_{T} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),80,0.,300.);
	
	hname = "jet_deltar"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);

	hname = "jet_deltar_b"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);
	h1[hname]->SetLineColor(quark_color["b"]);

	hname = "jet_deltar_cl"; 
	htitle = "#Delta R";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),60,0.,0.55);
	h1[hname]->SetLineColor(quark_color["c"]);

	hname = "jet_deltar_udsg"; 
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
	
	hname = "jet_ptrel"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),30,0.,4.);
	
	hname = "jet_ptrel_b"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),30,0.,4.);
	
	hname = "jet_ptrel_cl"; 
	htitle = "p_{Trel} [GeV/c]";
	h1[hname] = new TH1D(hname.c_str(),htitle.c_str(),30,0.,4.);
	
	hname = "jet_et";
	htitle = "E_{T} [GeV]";
	
}

void S8Plotter::Loop()
{
  
	Book();
	
	const int nptbins = 14;
	Double_t jetptbins[nptbins] = {30., 40., 50., 60., 70., 80, 90., 100., 120., 140., 160., 180., 230., 300.};
	
	Double_t jetetabins[8] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,2.0};

	h2["npT"] = new TH2F("npT","MuTag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["ppT"] = new TH2F("ppT","MuTag && CMBtag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["ncmbpT"] = new TH2F("ncmbpT","opp tag: MuTag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["pcmbpT"] = new TH2F("pcmbpT","opp tag MuTag && CMBtag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);

	h2["nEta"] = new TH2F("nEta","MuTag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["pEta"] = new TH2F("pEta","MuTag && CMBtag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["ncmbEta"] = new TH2F("ncmbEta","opp tag: MuTag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["pcmbEta"] = new TH2F("pcmbEta","opp tag MuTag && CMBtag Eta vs pTrel",7,jetetabins,50,0.,5.);

	h2["b_npT"] = new TH2F("b_npT","b MuTag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["b_ppT"] = new TH2F("b_ppT","b MuTag && CMBtag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["b_ncmbpT"] = new TH2F("b_ncmbpT","b opp tag: MuTag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["b_pcmbpT"] = new TH2F("b_pcmbpT","b opp tag MuTag && CMBtag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);

	h2["b_nEta"] = new TH2F("b_nEta","b MuTag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["b_pEta"] = new TH2F("b_pEta","b MuTag && CMBtag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["b_ncmbEta"] = new TH2F("b_ncmbEta","b opp tag: MuTag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["b_pcmbEta"] = new TH2F("b_pcmbEta","b opp tag MuTag && CMBtag Eta vs pTrel",7,jetetabins,50,0.,5.);

	h2["cl_npT"] = new TH2F("cl_npT","cl MuTag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["cl_ppT"] = new TH2F("cl_ppT","cl MuTag && CMBtag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["cl_ncmbpT"] = new TH2F("cl_ncmbpT","cl opp tag: MuTag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);
	h2["cl_pcmbpT"] = new TH2F("cl_pcmbpT","cl opp tag MuTag && CMBtag pT vs pTrel",nptbins-1,jetptbins,50,0.,5.);

	h2["cl_nEta"] = new TH2F("cl_nEta","cl MuTag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["cl_pEta"] = new TH2F("cl_pEta","cl MuTag && CMBtag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["cl_ncmbEta"] = new TH2F("cl_ncmbEta","cl opp tag: MuTag Eta vs pTrel",7,jetetabins,50,0.,5.);
	h2["cl_pcmbEta"] = new TH2F("cl_pcmbEta","cl opp tag MuTag && CMBtag Eta vs pTrel",7,jetetabins,50,0.,5.);
	
	h1["alpha"] = new TH1D("alpha","alpha",nptbins-1,jetptbins);
	h1["beta"] = new TH1D("beta","beta",nptbins-1,jetptbins);
	h1["kappa_cl"] = new TH1D("kappa_cl","kappa_cl",nptbins-1,jetptbins);
	h1["kappa_b"] = new TH1D("kappa_b","kappa_b",nptbins-1,jetptbins);
	
	// errors
	h2["npT"]->Sumw2();
	h2["ppT"]->Sumw2();
	h2["ncmbpT"]->Sumw2();
	h2["pcmbpT"]->Sumw2();
	h2["nEta"]->Sumw2();
	h2["pEta"]->Sumw2();
	h2["ncmbEta"]->Sumw2();
	h2["pcmbEta"]->Sumw2();
	h2["b_npT"]->Sumw2();
	h2["b_ppT"]->Sumw2();
	h2["b_ncmbpT"]->Sumw2();
	h2["b_pcmbpT"]->Sumw2();
	h2["b_nEta"]->Sumw2();
	h2["b_pEta"]->Sumw2();
	h2["b_ncmbEta"]->Sumw2();
	h2["b_pcmbEta"]->Sumw2();
	h2["cl_npT"]->Sumw2();
	h2["cl_ppT"]->Sumw2();
	h2["cl_ncmbpT"]->Sumw2();
	h2["cl_pcmbpT"]->Sumw2();
	h2["cl_nEta"]->Sumw2();
	h2["cl_pEta"]->Sumw2();
	h2["cl_ncmbEta"]->Sumw2();
	h2["cl_pcmbEta"]->Sumw2();
	//h1["alpha"]->Sumw2(); 
	//h1["beta"]->Sumw2(); 
	//h1["kappa_cl"]->Sumw2();
	//h1["kappa_b"]->Sumw2(); 
	
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	std::cout << " Total entries = " << fChain->GetEntries() << std::endl;
	//int count_multiple_mu = 0;
		
	//double TCdiscriminator = 4.0;
	//double TCdiscriminatorLoose = 2.0;
	
	std::vector< int > counters;
	TString tmpfilename = "";
	
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
		
		//counters[0]++;

		if (fVerbose) std::cout << " fill njets and nmuons histograms" << std::endl;
		h1["njets"]->Fill(fS8evt->njets);
		h1["nmuons"]->Fill(fS8evt->nmuons);

			
		std::vector<float> jet_pt_vec = fS8evt->jet_pt;
		int vec_size = jet_pt_vec.size();

		int njets_with_lepton = 0;
		int nopposite_jets = 0;
		
		double ptrel = 0.;
		bool isMuonJetSample = false;
		//bool isOppositeJetSample = false;
		bool passGoodMuon = false;
		bool passJetbTagger = false;
		bool passOppJetbTagger = false;
		//bool passptrel = false;
		
		
		TLorentzVector p4Jet;
		TLorentzVector p4OppJet;
		int JetFlavor = -1;
		int OppJetFlavor = -1;
		
		// first check that we have a lepton in a jet
		if (fVerbose) std::cout << " loop over jets, collection size = " << vec_size << std::endl;
		for ( int ijet =0; ijet != vec_size; ++ijet) {
			
			if ( fS8evt->jet_hasLepton[ijet] == 1 ) isMuonJetSample = true;
		}
		
		// no lepton+jet in this event therefore continue
		if (!isMuonJetSample) continue;
		if (fVerbose) std::cout << " found event with lepton in jet" << std::endl;
		//counters[1]++;
		
		for ( int ijet =0; ijet != vec_size; ++ijet) {

			double jetcorr = fS8evt->jetcorrection[ijet];
			//double jet_e = sqrt( fS8evt->jet_et[ijet]*fS8evt->jet_et[ijet] - fS8evt->jet_pt[ijet]*fS8evt->jet_pt[ijet] + fS8evt->jet_p[ijet]*fS8evt->jet_p[ijet] );
			
			if ( fS8evt->jet_hasLepton[ijet] == 1 ) {

				//isMuonJetSample = true;
				
				p4Jet.SetPtEtaPhiE(fS8evt->jet_pt[ijet], fS8evt->jet_eta[ijet], fS8evt->jet_phi[ijet], fS8evt->jet_e[ijet] );
				p4Jet = jetcorr * p4Jet;
								
				h1["jet_pt"]->Fill( p4Jet.Pt() );

				// get MC flavor of jet
				JetFlavor = fS8evt->jet_flavour_alg[ijet];
			
				BTagLeptonEvent Muons = fS8evt->lepton[ijet];	
				int mu_size = Muons.pt.size();

				if ( fVerbose && mu_size > 1 ) std::cout << " Muons in jet = " << mu_size << std::endl;
				int ith_mu_highest_pt = -1;
				double mu_highest_pt = 0;
				
				for ( int imu = 0; imu != mu_size; ++imu ) {

					if ( Muons.trkrechits[imu] >= 8 ) {
						passGoodMuon = true;
						njets_with_lepton++;

						if (Muons.pt[imu] > mu_highest_pt ) { mu_highest_pt = Muons.pt[imu]; ith_mu_highest_pt = imu; }
						
						//ptrel = Muons.jet_ptrel[imu];
						//if ( ptrel > 0.8 ) {
						//	passptrel = true;
						//	}

						if ( fVerbose && mu_size>1 ) std::cout << " muon " << imu << " pt= " << Muons.pt[imu] << " eta= " << Muons.eta[imu] << " chamber hits= " << Muons.SArechits[imu] << " chi2/ndof = " << Muons.chi2[imu]/Muons.ndof[imu] << " IPS= " << Muons.d0sigma[imu] << " mcpdgid= " << Muons.mc_pdgid[imu] << std::endl;
						
					}
				}
				// select only one muon in jet, the one with the highest pt
				if (passGoodMuon) {
					TLorentzVector vmu, vjet, vtot;
					vmu.SetPtEtaPhiE(Muons.pt[ith_mu_highest_pt], Muons.eta[ith_mu_highest_pt], Muons.phi[ith_mu_highest_pt], Muons.e[ith_mu_highest_pt]);
					vjet.SetPtEtaPhiE(fS8evt->jet_pt[ijet]*fS8evt->jetcorrection[ijet],
									  fS8evt->jet_eta[ijet], fS8evt->jet_phi[ijet],
									  fS8evt->jet_e[ijet]*fS8evt->jetcorrection[ijet]);
					vtot = vjet + vmu;
					ptrel = ( vmu.Px() * vtot.Px() 
									  + vmu.Py() * vtot.Py()
									  + vmu.Pz() * vtot.Pz() ) / vtot.P();
					ptrel = TMath::Sqrt( vmu.P() * vmu.P() - ptrel * ptrel );

					//ptrel = Muons.jet_ptrel[ith_mu_highest_pt];
					
					h1["muon_pt"]->Fill( Muons.pt[ith_mu_highest_pt] );
					h1["jet_deltar"]->Fill( Muons.jet_deltaR[ith_mu_highest_pt] );
					h1["jet_ptrel"]->Fill( ptrel );
					if ( JetFlavor == 5 ) {
						h1["jet_deltar_b"]->Fill( Muons.jet_deltaR[ith_mu_highest_pt] );
						h1["jet_ptrel_b"]->Fill( ptrel );
					}
					if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
						h1["jet_deltar_cl"]->Fill( Muons.jet_deltaR[ith_mu_highest_pt] );
						h1["jet_ptrel_cl"]->Fill( ptrel );
					}
					
					if ( fS8evt->btag_TrkCounting_disc3D_3trk[ijet] > 4. ) passJetbTagger = true;
				}
			}//check lepton in jet
			
		}// end loop over jets

		if (!passGoodMuon) continue;
		
		// find opposite tagged jet
		int ith_oppjet_highest_pt = -1;
		double oppjet_highest_pt = 0;
		for ( int ijet =0; ijet != vec_size; ++ijet) {
			// select jet other than the jet with lepton
			if ( fS8evt->jet_hasLepton[ijet] == 0 ) {
				double jetcorr = fS8evt->jetcorrection[ijet];
				//double jet_e = sqrt( fS8evt->jet_et[ijet]*fS8evt->jet_et[ijet] - fS8evt->jet_pt[ijet]*fS8evt->jet_pt[ijet] + fS8evt->jet_p[ijet]*fS8evt->jet_p[ijet] );
				
				if ( fS8evt->btag_TrkCounting_disc3D_3trk[ijet] > 4. ) {
					nopposite_jets++;
					passOppJetbTagger = true;
					p4OppJet.SetPtEtaPhiE(fS8evt->jet_pt[ijet], fS8evt->jet_eta[ijet], fS8evt->jet_phi[ijet], fS8evt->jet_e[ijet] );
					p4OppJet = jetcorr * p4OppJet;

					if ( p4OppJet.Pt() > oppjet_highest_pt ) { ith_oppjet_highest_pt = ijet; oppjet_highest_pt = p4OppJet.Pt(); }
					
					
				}
			}
		}
		if (passOppJetbTagger) {
			p4OppJet.SetPtEtaPhiE(fS8evt->jet_pt[ith_oppjet_highest_pt], fS8evt->jet_eta[ith_oppjet_highest_pt], fS8evt->jet_phi[ith_oppjet_highest_pt], fS8evt->jet_e[ith_oppjet_highest_pt] );
			double jetcorr = fS8evt->jetcorrection[ith_oppjet_highest_pt];
			p4OppJet = jetcorr * p4OppJet;
			OppJetFlavor = fS8evt->jet_flavour_alg[ith_oppjet_highest_pt];
			
			double deltaphi = TMath::Abs(p4Jet.Phi() - p4OppJet.Phi() );
			if ( deltaphi > TMath::Pi() ) deltaphi = TMath::Abs(2.*TMath::Pi() - deltaphi);
			h1["jet_deltaphi"]->Fill( deltaphi );
			if ( JetFlavor==5 && OppJetFlavor==5 ) h1["jet_deltaphi_b"]->Fill( deltaphi );

			if(fVerbose) {
				if (deltaphi == 0 ) std::cout << " deltaphi is zero " << std::endl;
				std::cout << " lepton-jet: px= " << p4Jet.Px() << " py= " << p4Jet.Py() << " pz = " << p4Jet.Pz() << " e = " << p4Jet.E() << std::endl;
				std::cout << " opp.   jet: px= " << p4OppJet.Px() << " py= " << p4OppJet.Py() << " pz = " << p4OppJet.Pz() << " e = " << p4OppJet.E() << std::endl;
			}
		}
		
		if (fVerbose && nopposite_jets>1) std::cout << " number of opposite jets = " << nopposite_jets << std::endl;
		
		
							
		h2["npT"]->Fill( p4Jet.Pt() , ptrel );
		h2["nEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		if ( JetFlavor == 5 ) {
			h2["b_npT"]->Fill( p4Jet.Pt() , ptrel );
			h2["b_nEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		}
		if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
			h2["cl_npT"]->Fill( p4Jet.Pt() , ptrel);
			h2["cl_nEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		}

		if ( passJetbTagger ) {
			h2["ncmbpT"]->Fill( p4Jet.Pt() , ptrel );
			h2["ncmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			if ( JetFlavor == 5 ) {
				h2["b_ncmbpT"]->Fill( p4Jet.Pt() , ptrel );
				h2["b_ncmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			}
			if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
				h2["cl_ncmbpT"]->Fill( p4Jet.Pt() , ptrel);
				h2["cl_ncmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			}
		}
		if ( passOppJetbTagger ) {
			h2["ppT"]->Fill( p4Jet.Pt() , ptrel );
			h2["pEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			if ( JetFlavor == 5 ) {
				h2["b_ppT"]->Fill( p4Jet.Pt() , ptrel );
				h2["b_pEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			}
			if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
				h2["cl_ppT"]->Fill( p4Jet.Pt() , ptrel);
				h2["cl_pEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			}
		}
		
		if ( passJetbTagger && passOppJetbTagger ) {
			h2["pcmbpT"]->Fill( p4Jet.Pt() , ptrel);
			h2["pcmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			if ( JetFlavor == 5 ) {
				h2["b_pcmbpT"]->Fill( p4Jet.Pt() , ptrel);
				h2["b_pcmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			}
			if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
				h2["cl_pcmbpT"]->Fill( p4Jet.Pt() , ptrel);
				h2["cl_pcmbEta"]->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			}
		}
		
	}// end loop over entries

	std::cout << " Total entries = " << fChain->GetEntries() << std::endl;
	
	//std::cout << " Events with a jet with multiple muons " <<  << std::endl;
	

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

	h1["eff_pTrel"] = new TH1D("eff_pTrel","eff_pTrel",nptbins-1,jetptbins);
	h1["eff_pTrel_b"] = new TH1D("eff_pTrel_b","eff_pTrel_b",nptbins-1,jetptbins);
	h1["eff_pTrel_cl"] = new TH1D("eff_pTrel_cl","eff_pTrel_cl",nptbins-1,jetptbins);

	h1["eff_TaggedJet"] = new TH1D("eff_TaggedJet","eff_TaggedJet",nptbins-1,jetptbins);
	h1["eff_TaggedJet_b"] = new TH1D("eff_TaggedJet_b","eff_TaggedJet_b",nptbins-1,jetptbins);
	h1["eff_TaggedJet_cl"] = new TH1D("eff_TaggedJet_cl","eff_TaggedJet_cl",nptbins-1,jetptbins);

	h1["eff_OppTaggedJet"] = new TH1D("eff_OppTaggedJet","eff_OppTaggedJet",nptbins-1,jetptbins);
	h1["eff_OppTaggedJet_b"] = new TH1D("eff_OppTaggedJet_b","eff_OppTaggedJet_b",nptbins-1,jetptbins);
	h1["eff_OppTaggedJet_cl"] = new TH1D("eff_OppTaggedJet_cl","eff_OppTaggedJet_cl",nptbins-1,jetptbins);

	h1["eff_TaggedBothJets"] = new TH1D("eff_TaggedBothJets","eff_TaggedBothJets",nptbins-1,jetptbins);
	h1["eff_TaggedBothJets_b"] = new TH1D("eff_TaggedBothJets_b","eff_TaggedBothJets_b",nptbins-1,jetptbins);
	h1["eff_TaggedBothJets_cl"] = new TH1D("eff_TaggedBothJets_cl","eff_TaggedBothJets_cl",nptbins-1,jetptbins);

	h1["eff_pTrel_TaggedJet"] = new TH1D("eff_pTrel_TaggedJet","eff_pTrel_TaggedJet",nptbins-1,jetptbins);
	h1["eff_pTrel_TaggedJet_b"] = new TH1D("eff_pTrel_TaggedJet_b","eff_pTrel_TaggedJet_b",nptbins-1,jetptbins);
	h1["eff_pTrel_TaggedJet_cl"] = new TH1D("eff_pTrel_TaggedJet_cl","eff_pTrel_TaggedJet_cl",nptbins-1,jetptbins);

	h1["eff_pTrel_OppTaggedJet"] = new TH1D("eff_pTrel_OppTaggedJet","eff_pTrel_OppTaggedJet",nptbins-1,jetptbins);
	h1["eff_pTrel_OppTaggedJet_b"] = new TH1D("eff_pTrel_OppTaggedJet_b","eff_pTrel_OppTaggedJet_b",nptbins-1,jetptbins);
	h1["eff_pTrel_OppTaggedJet_cl"] = new TH1D("eff_pTrel_OppTaggedJet_cl","eff_pTrel_OppTaggedJet_cl",nptbins-1,jetptbins);

	h1["eff_pTrel_TaggedBothJets"] = new TH1D("eff_pTrel_TaggedBothJets","eff_pTrel_TaggedBothJets",nptbins-1,jetptbins);
	h1["eff_pTrel_TaggedBothJets_b"] = new TH1D("eff_pTrel_TaggedBothJets_b","eff_pTrel_TaggedBothJets_b",nptbins-1,jetptbins);
	h1["eff_pTrel_TaggedBothJets_cl"] = new TH1D("eff_pTrel_TaggedBothJets_cl","eff_pTrel_TaggedBothJets_cl",nptbins-1,jetptbins);

	//TH1D *halljets = (TH1D*) h1["eff_pTrel"]->Clone("halljets");
	TH1D *halljets = h2["npT"]->ProjectionX("halljets", -1 , -1,"e");
	TH1D *halljets_b = h2["b_npT"]->ProjectionX("halljets_b", -1 , -1,"e");
	TH1D *halljets_cl = h2["cl_npT"]->ProjectionX("halljets_cl", -1 , -1,"e");

	h1["eff_pTrel"]->Divide( h2["npT"]->ProjectionX("halljets_ptrel", 9 , -1,"e") , halljets ,1.,1.,"B");
	h1["eff_pTrel_b"]->Divide( h2["b_npT"]->ProjectionX("b_halljets_ptrel", 9 , -1,"e") , halljets_b ,1.,1.,"B");
	h1["eff_pTrel_cl"]->Divide( h2["cl_npT"]->ProjectionX("cl_halljets_ptrel", 9 , -1,"e") , halljets_cl ,1.,1.,"B");
	
	h1["eff_TaggedJet"]->Divide( h2["ncmbpT"]->ProjectionX("halljets_tagged", -1 , -1,"e") , halljets ,1.,1.,"B");
	h1["eff_TaggedJet_b"]->Divide( h2["b_ncmbpT"]->ProjectionX("b_halljets_tagged", -1 , -1,"e") , halljets_b ,1.,1.,"B");
	h1["eff_TaggedJet_cl"]->Divide( h2["cl_ncmbpT"]->ProjectionX("cl_halljets_tagged", -1 , -1,"e") , halljets_cl ,1.,1.,"B");

	h1["eff_pTrel_TaggedJet"]->Divide( h2["ncmbpT"]->ProjectionX("halljets_ptreltagged", 9 , -1,"e") , halljets ,1.,1.,"B");
	h1["eff_pTrel_TaggedJet_b"]->Divide( h2["b_ncmbpT"]->ProjectionX("b_halljets_ptreltagged", 9 , -1,"e") , halljets_b ,1.,1.,"B");
	h1["eff_pTrel_TaggedJet_cl"]->Divide( h2["cl_ncmbpT"]->ProjectionX("cl_halljets_ptrelltagged", 9 , -1,"e") , halljets_cl ,1.,1.,"B");


	TH1D *halloppjets = h2["ppT"]->ProjectionX("halloppjets", -1 , -1,"e");
	TH1D *halloppjets_b = h2["b_ppT"]->ProjectionX("halloppjets_b", -1 , -1,"e");
	TH1D *halloppjets_cl = h2["cl_ppT"]->ProjectionX("halloppjets_cl", -1 , -1,"e");

	h1["eff_TaggedBothJets"]->Divide( h2["pcmbpT"]->ProjectionX("halloppjets_tagged", -1 , -1,"e") , halloppjets ,1.,1.,"B");
	h1["eff_TaggedBothJets_b"]->Divide( h2["b_pcmbpT"]->ProjectionX("b_halloppjets_tagged", -1 , -1,"e") , halloppjets_b ,1.,1.,"B");
	h1["eff_TaggedBothJets_cl"]->Divide( h2["cl_pcmbpT"]->ProjectionX("cl_halloppjets_tagged", -1 , -1,"e") , halloppjets_cl ,1.,1.,"B");

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

	cv_map["kappa_b"] = new TCanvas("kappa_b","kappa_b",700,700);
	h1["kappa_b"]->SetMarkerStyle(8);
	h1["kappa_b"]->SetMarkerSize(1.5);
	h1["kappa_b"]->SetLineColor(2);
	h1["kappa_b"]->SetMarkerColor(2);
	h1["kappa_b"]->Draw("PE1");
	h1["kappa_b"]->Fit("pol1","0");
	TF1 *f1_kb = h1["kappa_b"]->GetFunction("pol1");
	f1_kb->SetLineColor(2);
	f1_kb->Draw("same");
	h1["eff_TaggedJet_b"]->Draw("PE1 same");
	h1["eff_pTrel_b"]->Draw("PE1 same");
	
    //______________________________________________________

	cv_map["kappa_cl"] = new TCanvas("kappa_cl","kappa_cl",700,700);
	h1["kappa_cl"]->SetMarkerStyle(8);
	h1["kappa_cl"]->SetMarkerSize(1.5);
	h1["kappa_cl"]->SetLineColor(2);
	h1["kappa_cl"]->SetMarkerColor(2);
	h1["kappa_cl"]->Draw("PE1");
	h1["kappa_cl"]->Fit("pol1","0");
	TF1 *f1_kcl = h1["kappa_cl"]->GetFunction("pol1");
	f1_kcl->SetLineColor(3);
	f1_kcl->Draw("same");
	h1["eff_TaggedJet_cl"]->Draw("PE1 same");
	h1["eff_pTrel_cl"]->Draw("PE1 same");
    //______________________________________________________
	
	cv_map["alpha"] = new TCanvas("alpha","alpha",700,700);
	h1["alpha"]->SetMarkerStyle(8);
	h1["alpha"]->SetMarkerSize(1.5);
	h1["alpha"]->Draw("PE1");
	h1["alpha"]->Fit("pol1","0");
	TF1 *f1_alpha = h1["alpha"]->GetFunction("pol1");
	f1_alpha->SetLineColor(3);
	f1_alpha->Draw("same");
	h1["eff_TaggedBothJets_cl"]->Draw("PE1 same");
	h1["eff_TaggedJet_cl"]->Draw("PE1 same");
    //______________________________________________________
	
	cv_map["beta"] = new TCanvas("beta","beta",700,700);
	h1["beta"]->SetMarkerStyle(8);
	h1["beta"]->SetMarkerSize(1.5);
	h1["beta"]->Draw("PE1");
	h1["beta"]->Fit("pol1","0");
	TF1 *f1_beta = h1["beta"]->GetFunction("pol1");
	f1_beta->SetLineColor(3);
	f1_beta->Draw("same");
	h1["eff_TaggedBothJets_b"]->Draw("PE1 same");
	h1["eff_TaggedJet_b"]->Draw("PE1 same");
		
}




		



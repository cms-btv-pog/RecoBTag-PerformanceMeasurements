
// CMS includes
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "PhysicsTools/FWLite/interface/EventContainer.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
//#include "RecoBTag/PerformanceMeasurements/bin/PMHistograms.h"

// ROOT includes
#include "TAxis.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "TGraphErrors.h"
#if !defined(__CINT__) && !defined(__MAKECINT__)
#endif

// c++ standards
#include <iostream>

using namespace std;

int main (int argc, char* argv[]) 
{

  ////////////////////////////////////////////////
  // // Command Line Options // //
  ////////////////////////////////////////////////
  
  optutl::CommandLineParser parser ("Plots Jet Pt");
  parser.stringValue ("outputFile") = "PM_results"; // .root added automatically
  // Parse the command line arguments
  parser.parseArguments (argc, argv);

  //////////////////////////////////
  // // Create Event Container // //
  //////////////////////////////////

  fwlite::EventContainer eventCont (parser);

  // use this histogram manager for the moment
  // we will change it to use FY or SK managers
  //PMHistograms pmhistos( eventCont );
  //pmhistos.Add();
    
  std::map<std::string, int>  quark_color;
  quark_color[""] = 1;
  quark_color["b"] = 2;
  quark_color["c"] = 3;
  quark_color["uds"] = 4;
  quark_color["g"] = 6;

  TString           ftagger;
  TString           flevel;
  TString           fAwaytagger;
  TString           fAwaylevel;
  std::map< TString, float > fTrackCountingMap;
  TAxis             fJetPtAxis;
  TAxis             fJetEtaAxis;
  TAxis             fCorrPtAxis;
  TAxis             fCorrEtaAxis;

  ftagger = "TrackCounting";
  flevel  = "Loose";
  fAwaytagger = "TrackCounting";
  fAwaylevel = "Loose";
  
  fTrackCountingMap["Loose"]  = 2.0; // use TC2:high eff.
  fTrackCountingMap["Medium"] = 4.2; // use TC2:high eff.
  fTrackCountingMap["Tight"]  = 4.1;
  
  const int nptarray = 4;
  const int netaarray = 10;
  const int ncorrptarray = 5;
  const int ncorretaarray = 5;
  Double_t jetptbins[nptarray] = {30., 70., 120.,230.};
  Double_t jetetabins[netaarray] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5};
  Double_t corrptbins[ncorrptarray] = {20.,40.,60.,80.,230.};
  Double_t corretabins[ncorrptarray] = {0.,0.5,1.,1.5,2.5};
  
  int nptbins = fJetPtAxis.GetNbins();;
  //const Double_t *jetptbins = (fJetPtAxis.GetXbins())->GetArray();
  int netabins = fJetEtaAxis.GetNbins();
  //const Double_t *jetetabins = (fJetEtaAxis.GetXbins())->GetArray();
	// int ncorrptbins = fCorrPtAxis.GetNbins();
	// const Double_t *corrptbins = (fCorrPtAxis.GetXbins())->GetArray();
	// int ncorretabins = fCorrEtaAxis.GetNbins();
	// const Double_t *corretabins = (fCorrEtaAxis.GetXbins())->GetArray();
	
	eventCont.add( new TH2F("npT","MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("ppT","MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("ncmbpT","opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("pcmbpT","opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );

	eventCont.add( new TH2F("qpT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("qcmbpT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	
	
	eventCont.add( new TH2F("nEta","MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("pEta","MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("ncmbEta","opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("pcmbEta","opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );

	eventCont.add( new TH2F("qEta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("qcmbEta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.) );
	
	eventCont.add( new TH2F("b_npT","b MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("b_ppT","b MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("b_ncmbpT","b opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("b_pcmbpT","b opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );

	eventCont.add( new TH2F("b_qpT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("b_qcmbpT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	
	eventCont.add( new TH2F("b_nEta","b MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("b_pEta","b MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("b_ncmbEta","b opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("b_pcmbEta","b opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );

	eventCont.add( new TH2F("b_qEta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("b_qcmbEta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.) );
	
	
	eventCont.add( new TH2F("cl_npT","cl MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("cl_ppT","cl MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("cl_ncmbpT","cl opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("cl_pcmbpT","cl opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );

	eventCont.add( new TH2F("cl_qpT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("cl_qcmbpT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	
	eventCont.add( new TH2F("c_npT","c MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("c_ncmbpT","c opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("g_npT","g MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("g_ncmbpT","g opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("uds_npT","uds MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );
	eventCont.add( new TH2F("uds_ncmbpT","uds opp tag MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.) );

	eventCont.add( new TH2F("cl_nEta","cl MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("cl_pEta","cl MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("cl_ncmbEta","cl opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("cl_pcmbEta","cl opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.) );

	eventCont.add( new TH2F("cl_qEta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.) );
	eventCont.add( new TH2F("cl_qcmbEta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.) );
	
	eventCont.add( new TH2F("ptVsEta","pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );
	eventCont.add( new TH2F("ptVsEta_b","pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );
	eventCont.add( new TH2F("taggedjet_ptVsEta","taggedjet pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );
	eventCont.add( new TH2F("taggedjet_ptVsEta_b","taggedjet pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );

	eventCont.add( new TH1D("alpha","alpha",nptbins,jetptbins) );
	eventCont.add( new TH1D("beta","beta",nptbins,jetptbins) );
	eventCont.add( new TH1D("kappa_cl","kappa_cl",nptbins,jetptbins) );
	eventCont.add( new TH1D("kappa_b","kappa_b",nptbins,jetptbins) );
	eventCont.add( new TH1D("delta","delta",nptbins,jetptbins) );
	eventCont.add( new TH1D("gamma","gamma",nptbins,jetptbins) );

	eventCont.add( new TH1D("alpha_eta","alpha_eta",netabins,jetetabins) );
	eventCont.add( new TH1D("beta_eta","beta_eta",netabins,jetetabins) );
	eventCont.add( new TH1D("kappa_eta_cl","kappa_eta_cl",netabins,jetetabins) );
	eventCont.add( new TH1D("kappa_eta_b","kappa_eta_b",netabins,jetetabins) );
	eventCont.add( new TH1D("delta_eta","delta_eta",netabins,jetetabins) );
	eventCont.add( new TH1D("gamma_eta","gamma_eta",netabins,jetetabins) );

	eventCont.add( new TH1D("jet_deltaR","#Delta R",60,0.,0.55) );
	eventCont.add( new TH1D("jet_deltaR_b","#Delta R",60,0.,0.55) );
	eventCont.add( new TH1D("jet_deltaR_c","#Delta R",60,0.,0.55) );
	eventCont.add( new TH1D("jet_deltaR_udsg","#Delta R",60,0.,0.55) );
	// change colors
	eventCont.add( new TH1D("jet_pTrel","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
	eventCont.add( new TH1D("jet_pTrel_b","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
	eventCont.add( new TH1D("jet_pTrel_c","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
	eventCont.add( new TH1D("jet_pTrel_udsg","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
	
	eventCont.add( new TH1F( "jet_pt", "jet pt", 30, 0, 150) );
	eventCont.add( new TH1F( "muon_pt", "muon pt", 300, 0, 50) );
	eventCont.add( new TH1F( "ptRel", "ptRel", 100, 0, 10) );
	
  // //////////////// //
  // // Event Loop // //
  // //////////////// //

  for (eventCont.toBegin(); ! eventCont.atEnd(); ++eventCont) 
    {

      // load object collections
      fwlite::Handle< vector< pat::Jet > > jetHandle;
      jetHandle.getByLabel (eventCont, "selectedLayer1Jets");
      assert ( jetHandle.isValid() );

      fwlite::Handle< vector< pat::Muon > > muonHandle;
      muonHandle.getByLabel (eventCont, "selectedLayer1Muons");
      assert ( muonHandle.isValid() );

	  TLorentzVector p4Jet;
	  TLorentzVector p4MuJet;
	  TLorentzVector p4Muon;
	  TLorentzVector p4AwayJet;
	  // Loop over jets 
	  for (vector< pat::Jet >::const_iterator jetIter = jetHandle->begin();
	       jetIter != jetHandle->end(); ++jetIter)
	  {
	      // select a good jet
	      if ( jetIter->pt() <= 30.|| std::abs( jetIter->eta() ) >= 2.4 ) continue;
		  
	      eventCont.hist("jet_pt")->Fill (jetIter->pt()); // just for testing
        // get MC flavor of jet
		  
		  int JetFlavor = jetIter->partonFlavour();
		  p4Jet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );
		  int hasLepton = 0;
		  int tmptotmuon = 0;
		  // loop over muons
		  ////////////////////////////////
		  double mu_highest_pt = 0;
		  double ptrel = 0;
		  double ptreltmp = 0;

		  for (vector< pat::Muon >::const_iterator muonIter = muonHandle->begin();
			   muonIter != muonHandle->end(); ++muonIter) 
		  {
			  if (muonIter->isGlobalMuon() == false) continue;
			  int nhit =  (*(muonIter->innerTrack())).numberOfValidHits();

			  // muon cuts
			  double normChi2 = (*(muonIter->combinedMuon())).normalizedChi2();
			  if ( (nhit <= 7 ) || (muonIter->pt()<= 5.0) || (normChi2 >= 5.0 ) ) continue;

			  eventCont.hist("muon_pt")->Fill (muonIter->pt());

			  // find a muon in a jet
			  //check deltar
			  double deltaR  = ROOT::Math::VectorUtil::DeltaR(jetIter->p4().Vect(), muonIter->p4().Vect() );
			  TVector3 tmpvecOrg(jetIter->p4().Vect().X(),jetIter->p4().Vect().Y(),  jetIter->p4().Vect().Z());
			  TVector3 tmpvec;
			  tmpvec = tmpvecOrg;
			  TVector3 leptonvec(muonIter->px(), muonIter->py(),muonIter->pz());
			  tmpvec += leptonvec;
			  ptreltmp = leptonvec.Perp(tmpvec);
			  // muon in jet cuts
			  if ( (deltaR >= 0.4) || (ptreltmp <= -1.0) ) continue;
			  hasLepton = 1;
			  if ( muonIter->pt() > mu_highest_pt )
			    {
			      mu_highest_pt = muonIter->pt();
			      p4Muon.SetPtEtaPhiE(muonIter->pt(),
						  muonIter->eta(),
						  muonIter->phi(),
						  muonIter->energy());
			      // recalculate pTrel
			      tmpvec = tmpvecOrg;
			      leptonvec.SetXYZ(muonIter->px(), muonIter->py(),muonIter->pz());
			      tmpvec += leptonvec;
			      ptrel = leptonvec.Perp(tmpvec);  // maximum
			    }
			  eventCont.hist("jet_pTrel")->Fill (ptrel);
			  eventCont.hist("jet_deltaR")->Fill (deltaR);
			  if ( JetFlavor == 5 ) {
				  eventCont.hist("jet_deltaR_b")->Fill (deltaR);
				  eventCont.hist("jet_pTrel_b")->Fill (ptrel);
			  }
			  if ( JetFlavor == 4 ) {
				  eventCont.hist("jet_deltaR_c")->Fill (deltaR);
				  eventCont.hist("jet_pTrel_c")->Fill (ptrel);
			  }
			  if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor==21 ) {
				  eventCont.hist("jet_deltaR_udsg")->Fill (deltaR);
				  eventCont.hist("jet_pTrel_udsg")->Fill (ptrel);
			  }

			  
		  }//muons
		  if ( hasLepton == 1 )
		    {
		      p4MuJet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );
		    }

		  
		  eventCont.hist("npT")->Fill( p4Jet.Pt() , ptrel );
		  eventCont.hist("nEta")->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		  if ( JetFlavor == 5 ) {
			  eventCont.hist("b_npT")->Fill( p4Jet.Pt() , ptrel );
			  eventCont.hist("b_nEta")->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		  }
		  if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
			  eventCont.hist("c_npT")->Fill( p4Jet.Pt() , ptrel );
			  eventCont.hist("c_nEta")->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
		  }
		  

		  // find away jet
		  ////////////////////////////
		  bool AwayTaggedJet = false;
		  bool AwayMuonJet = false;
		  TLorentzVector p4AwayMuon;
		  for (vector< pat::Jet >::const_iterator awayjetIter = jetHandle->begin();
		       awayjetIter != jetHandle->end(); ++awayjetIter)
		    {
		      if ( hasLepton == 0 ) continue;

		      TLorentzVector p4AwayJet;
		      p4AwayJet.SetPtEtaPhiE(awayjetIter->pt(), awayjetIter->eta(), awayjetIter->phi(), awayjetIter->energy() );
		      // Jet quality cuts
		      if ( (awayjetIter->pt())  <= 30 || std::abs( awayjetIter->eta() ) >= 2.4 ) continue;

		      // skip muon in jet
		      if ( p4AwayJet == p4MuJet ) continue;

		      // now we have an away jet
		      // find an away tagged jet
		      //std::cout << " find an away tagged jet" << std::endl;
			
			if ( !AwayTaggedJet )
			  {
			    double jetBDiscr_track_count_high_eff  = awayjetIter -> bDiscriminator( "trackCountingHighEffBJetTags" );
			    //std::cout <<" bDiscriminator trackCountingHighEffBJetTags = " <<jetBDiscr_track_count_high_eff << std::endl;
			    if (jetBDiscr_track_count_high_eff > 1.66) //loose operating point
			      {
				AwayTaggedJet = true;
			      }

			  }
			// find an away muon in jet
			if ( !AwayMuonJet )
			  {
			    mu_highest_pt = 0;
			    for (vector< pat::Muon >::const_iterator muonIter = muonHandle->begin();
				 muonIter != muonHandle->end(); ++muonIter) 
			      {
				if (muonIter->isGlobalMuon() == false) continue;
				//int nhit = muonIter->numberOfValidHits();
				int nhit =  (*(muonIter->innerTrack())).numberOfValidHits();
				// muon cuts
				double normChi2 = (*(muonIter->combinedMuon())).normalizedChi2();
				if ( (nhit <= 7 ) || (muonIter->pt()<= 5.0) || (normChi2 >= 5.0 ) ) continue;

				// find a muon in a jet
				//check deltar
				double deltaR  = ROOT::Math::VectorUtil::DeltaR(jetIter->p4().Vect(), muonIter->p4().Vect() );
				TVector3 tmpvecOrg(awayjetIter->p4().Vect().X(),awayjetIter->p4().Vect().Y(),  awayjetIter->p4().Vect().Z());
				TVector3 tmpvec;
				tmpvec = tmpvecOrg;
				TVector3 leptonvec(muonIter->px(), muonIter->py(),muonIter->pz());
				tmpvec += leptonvec;
				double awayptrel = leptonvec.Perp(tmpvec);
				// muon in jet cuts
				if ( (deltaR >= 0.4) || (ptreltmp <= -1.0) ) continue;
				// now we have a good muon in a jet
				AwayMuonJet = true;
				// pick the leading muon inside the jet
				if ( muonIter->pt() > mu_highest_pt )
				  {
				    mu_highest_pt = muonIter->pt();
				    p4AwayMuon.SetPtEtaPhiE(muonIter->pt(),
							    muonIter->eta(),
							    muonIter->phi(),
							    muonIter->energy());
				    // recalculate pTrel
				    tmpvec = tmpvecOrg;
				    leptonvec.SetXYZ(muonIter->px(), muonIter->py(),muonIter->pz());
				    tmpvec += leptonvec;
				    awayptrel = leptonvec.Perp(tmpvec);
				  }
			      }	
			  }	
		    } // close away jet loop


	  }//jets

    }//events

  return 0;

}

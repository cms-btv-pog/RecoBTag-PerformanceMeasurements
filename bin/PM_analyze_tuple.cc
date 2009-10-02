
// CMS includes
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "RecoBTag/PerformanceMeasurements/interface/TH1Store.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "RecoBTag/PerformanceMeasurements/interface/PMHistograms.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagHistograms.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

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

  //___________________
  // FWLite setup
  AutoLibraryLoader::enable();

  fwlite::ChainEvent events( parser.stringVector ("inputFiles") );
  
  //fwlite::EventBase *eventBase = & events;

  // for now
  std::string Tagger =  "trackCountingHighPurBJetTags";
  std::string awayTagger = "trackCountingHighPurBJetTags";
  double btag_cut_Tagger = 1.9;//loose
  double btag_cut_awayTagger = 3.99; //medium
  //__________________
  // Histogram manager

  TH1Store *hstore = new TH1Store();
  PMHistograms histos( hstore );
  histos.Add();

  // initialize histograms
  /*
  BTagHistograms *EffHistos     = new BTagHistograms();
  BTagHistograms *PtrelHistos   = new BTagHistograms();
  BTagHistograms *MujetHistos   = new BTagHistograms();
  BTagHistograms *AwayjetHistos = new BTagHistograms();
  BTagHistograms *TaggedMujetHistos   = new BTagHistograms();
  BTagHistograms *TaggedAwayjetHistos = new BTagHistograms();
  BTagHistograms *MujetHistos_mc   = new BTagHistograms();
  BTagHistograms *AwayjetHistos_mc = new BTagHistograms();
  BTagHistograms *TaggedMujetHistos_mc   = new BTagHistograms();
  BTagHistograms *TaggedAwayjetHistos_mc = new BTagHistograms();
  
  EffHistos->Init("efficiencies");
  PtrelHistos->Init("ptrel");
  MujetHistos->Init("n");
  AwayjetHistos->Init("p");
  MujetHistos_mc->Init("n","b");
  AwayjetHistos_mc->Init("p","b");
  MujetHistos_mc->Init("n","cl");
  AwayjetHistos_mc->Init("p","cl");
  MujetHistos_mc->Init("n","c");
  AwayjetHistos_mc->Init("p","c");
  MujetHistos_mc->Init("n","l");
  AwayjetHistos_mc->Init("p","l");

  std::string aalias = "TCHEM";
  EffHistos->Init("efficiencies",aalias);
  PtrelHistos->Init("ptrel",aalias);
			
  TaggedMujetHistos->Init("ntag",aalias);
  TaggedAwayjetHistos->Init("ptag",aalias);
  TaggedMujetHistos_mc->Init("ntag","b",aalias);
  TaggedAwayjetHistos_mc->Init("ptag","b",aalias);
  TaggedMujetHistos_mc->Init("ntag","cl",aalias);
  TaggedAwayjetHistos_mc->Init("ptag","cl",aalias);
  TaggedMujetHistos_mc->Init("ntag","c",aalias);
  TaggedAwayjetHistos_mc->Init("ptag","c",aalias);
  TaggedMujetHistos_mc->Init("ntag","l",aalias);
  TaggedAwayjetHistos_mc->Init("ptag","l",aalias);
  */
  // //////////////// //
  // // Event Loop // //
  // //////////////// //
  int nentry = 0;
  
  for (events.toBegin(); ! events.atEnd(); ++events) 
  {

	  std::cout << "Entry # " << nentry << std::endl;
	  nentry++;
      // load object collections
      fwlite::Handle< vector< pat::Jet > > jetHandle;
      jetHandle.getByLabel ( events, "selectedLayer1Jets");
      assert ( jetHandle.isValid() );

      fwlite::Handle< vector< pat::Muon > > muonHandle;
      muonHandle.getByLabel ( events, "selectedLayer1Muons");
      assert ( muonHandle.isValid() ); // we should always have muons because of the pre-selections

	  TLorentzVector p4Jet;
	  TLorentzVector p4MuJet;
	  TLorentzVector p4Muon;
	  TLorentzVector p4AwayJet;
	 
	  // Loop over jets 
	  for (vector< pat::Jet >::const_iterator jetIter = jetHandle->begin();
	       jetIter != jetHandle->end(); ++jetIter)
	  {
		  //std::cout << " jet pt " << jetIter->pt() << std::endl;
		   bool TaggedJet = false;
	      // select a good jet
	      if ( jetIter->pt() <= 30.|| std::abs( jetIter->eta() ) >= 2.4 ) continue;

		  //std::cout << " jet pt= " << jetIter->pt() << std::endl;
		  
	      hstore->hist("jet_pt")->Fill (jetIter->pt()); // just for testing
        // get MC flavor of jet
		  
		  int JetFlavor = abs( jetIter->partonFlavour() );
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
		    //std::cout << " muon pt " << muonIter->pt() << std::endl;
			  
			  if (muonIter->isGlobalMuon() == false) continue;
			  //int nhit =  (*(muonIter->innerTrack())).numberOfValidHits();

			  // muon cuts
			  double normChi2 = (*(muonIter->combinedMuon())).normalizedChi2();
			  //if ( (nhit <= 7 ) || (muonIter->pt()<= 5.0) || (normChi2 >= 5.0 ) ) continue;
			  if ( (muonIter->pt()<= 5.0) || (normChi2 >= 5.0 ) ) continue;

			  hstore->hist("muon_pt")->Fill (muonIter->pt());

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
			  hstore->hist("jet_pTrel")->Fill (ptrel);
			  hstore->hist("jet_deltaR")->Fill (deltaR);
			  if ( JetFlavor == 5 ) {
				  hstore->hist("jet_deltaR_b")->Fill (deltaR);
				  hstore->hist("jet_pTrel_b")->Fill (ptrel);
			  }
			  if ( JetFlavor == 4 ) {
				  hstore->hist("jet_deltaR_c")->Fill (deltaR);
				  hstore->hist("jet_pTrel_c")->Fill (ptrel);
			  }
			  if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor==21 ) {
				  hstore->hist("jet_deltaR_udsg")->Fill (deltaR);
				  hstore->hist("jet_pTrel_udsg")->Fill (ptrel);
			  }

			  
		  }//muons
		  if ( hasLepton == 1 )
		  {
		      p4MuJet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );

			  double btag   = jetIter -> bDiscriminator( Tagger );
			  //std::cout << "btag = " << btag << std::endl;
			  
			  if (btag > btag_cut_Tagger ) TaggedJet = true;
			  //std::cout << "flavor = " << JetFlavor << std::endl;
			  histos.FillHistos("n", p4MuJet, ptrel, JetFlavor, TaggedJet );
			  
			  //hstore->hist("npT")->Fill( p4Jet.Pt() , ptrel );
			  //MujetHistos->Fill2d("n_pT",p4MuJet.Pt(),ptrel);
			  //MujetHistos->Fill2d("n_eta",TMath::Abs(p4MuJet.Eta()),ptrel);

			  //hstore->hist("nEta")->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			  //if ( JetFlavor == 5 ) {
			  //  hstore->hist("b_npT")->Fill( p4Jet.Pt() , ptrel );
			  //  hstore->hist("b_nEta")->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			  //}
			  //if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor==21 ) {
			  //  hstore->hist("cl_npT")->Fill( p4Jet.Pt() , ptrel );
			  //  hstore->hist("cl_nEta")->Fill(TMath::Abs( p4Jet.Eta() ), ptrel );
			  //}
		  

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
					  btag  = awayjetIter -> bDiscriminator( awayTagger );
					  //std::cout <<" bDiscriminator trackCountingHighEffBJetTags = " <<jetBDiscr_track_count_high_eff << std::endl;
					  if ( btag > btag_cut_awayTagger) //loose operating point
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
						  //int nhit =  (*(muonIter->innerTrack())).numberOfValidHits();
						  // muon cuts
						  double normChi2 = (*(muonIter->combinedMuon())).normalizedChi2();
						  //if ( (nhit <= 7 ) || (muonIter->pt()<= 5.0) || (normChi2 >= 5.0 ) ) continue;
						  if ( (muonIter->pt()<= 5.0) || (normChi2 >= 5.0 ) ) continue;

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
				  } // away muon jet
			  } // close away jet loop

			  if (AwayTaggedJet)
				  histos.FillHistos("p", p4MuJet, ptrel, JetFlavor, TaggedJet );
			  
		  }// muon in jet

	  }//jets

	  //std::cout << " event done. " << std::endl;
	  
    }//events

  std::cout << " done. " << std::endl;
  
  std::string m_outputName  = parser.stringValue  ("outputFile");
  if (optutl::CommandLineParser::kStringVector == 
       parser.hasOption("inputFiles"))
   {
      hstore->write (m_outputName, 
                         parser.argVec(),
                         parser.stringVector ("inputFiles"));
   } else {
      hstore->write (m_outputName, 
                         parser.argVec());
  }
  
   
  return 0;

}

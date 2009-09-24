
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

// ROOT includes
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
		  
//		  int JetFlavor = abs(getMatchedParton(*jetIter).getFlavour());
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
			  eventCont.hist("ptRel")->Fill (ptrel);
			  
		  }//muons
		  if ( hasLepton == 1 )
		    {
		      p4MuJet.SetPtEtaPhiE(jetIter->pt(), jetIter->eta(), jetIter->phi(), jetIter->energy() );
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

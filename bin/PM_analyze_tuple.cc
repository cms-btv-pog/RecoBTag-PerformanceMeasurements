
// CMS includes
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "PhysicsTools/FWLite/interface/EventContainer.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 

// ROOT includes
#include "TFile.h"

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

      // loop over muons

      for (vector< pat::Muon >::const_iterator muonIter = muonHandle->begin();
           muonIter != muonHandle->end(); ++muonIter) 
	{         
	  // select a good muon
	  if ( muonIter->pt() <= 5 ) continue;
	  eventCont.hist("muon_pt")->Fill (muonIter->pt());

	  // find a muon in a jet

	  for (vector< pat::Jet >::const_iterator jetIter = jetHandle->begin();
	       jetIter != jetHandle->end(); ++jetIter)
	    {
	      // select a good jet
	      eventCont.hist("jet_pt")->Fill (jetIter->pt()); // just for testing

	      // check delta R

	    }//jets
  
	}//muons

    }//events

  return 0;

}

// -*- C++ -*-
//
// Package:    EventCounter
// Class:      EventCounter
//
/**\class EventCounter EventCounter.cc RecoBTag/PerformanceMeasurements/plugins/EventCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dinko Ferencek
//         Created:  Sat Jul 13 19:13:45 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
//
// class declaration
//

class EventCounter : public edm::EDAnalyzer {
   public:
      explicit EventCounter(const edm::ParameterSet&);
      ~EventCounter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;

      TH1D *hEventCount;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventCounter::EventCounter(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   hEventCount = fs->make<TH1D>("hEventCount","Even Count", 1, -0.5, 0.5);
}


EventCounter::~EventCounter()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   hEventCount->Fill(0.);
}


// ------------ method called once each job just before starting event loop  ------------
void
EventCounter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EventCounter::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
EventCounter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
EventCounter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
EventCounter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
EventCounter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);

// -*- C++ -*-
//
// Package:    TTbarSelectionFilter
// Class:      TTbarSelectionFilter
// 
/**\class TTbarSelectionFilter TTbarSelectionFilter.cc bTag/TTbarSelectionFilter/src/TTbarSelectionFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Jeremy,B25/117,6262,
//         Created:  Tue Dec  4 15:46:22 CET 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class TTbarSelectionFilter : public edm::EDFilter {
   public:
      explicit TTbarSelectionFilter(const edm::ParameterSet&);
      ~TTbarSelectionFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
     edm::InputTag channel_;
     
     
     
   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      bool select_ee_ ;
      bool select_mumu_ ;
      bool select_emu_ ;
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
TTbarSelectionFilter::TTbarSelectionFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   channel_      = iConfig.getParameter<edm::InputTag> ("channel");
   select_ee_            = iConfig.getParameter<bool > ("select_ee");
   select_mumu_            = iConfig.getParameter<bool > ("select_mumu");
   select_emu_            = iConfig.getParameter<bool > ("select_emu");
}


TTbarSelectionFilter::~TTbarSelectionFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TTbarSelectionFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   Handle<int> pIn;
   iEvent.getByLabel(channel_, pIn);
   
   bool isTTbar = false;

   //std::cout << " valeur de la selection ttbar " << *pIn << std::endl;
   //if(   *pIn>= 0 ) isTTbar = true;
   if (*pIn==0 && select_ee_)   isTTbar = true;
   if (*pIn==1 && select_mumu_) isTTbar = true;
   if (*pIn==2 && select_emu_)  isTTbar = true;


   return isTTbar;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTbarSelectionFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTbarSelectionFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
TTbarSelectionFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TTbarSelectionFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TTbarSelectionFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TTbarSelectionFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTbarSelectionFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TTbarSelectionFilter);

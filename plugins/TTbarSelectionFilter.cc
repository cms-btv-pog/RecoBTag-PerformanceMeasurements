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
// $Id: TTbarSelectionFilter.cc,v 1.1 2013/03/05 11:13:13 ccollard Exp $
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
      edm::EDGetTokenT<int> ttbartop_;   
     
   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
     std::vector<int> selectChannels_;
     bool selectAll_;
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
TTbarSelectionFilter::TTbarSelectionFilter(const edm::ParameterSet& iConfig):
  ttbartop_(consumes<int>(edm::InputTag("ttbarselectionproducer:topChannel")))
{
   //now do what ever initialization is needed
   selectChannels_ = iConfig.getParameter<std::vector<int> > ("selectChannels");
   selectAll_     = iConfig.getParameter<bool > ("selectAll");
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
   iEvent.getByToken(ttbartop_, pIn);   
   std::vector<int>::iterator it = find (selectChannels_.begin(), selectChannels_.end(), *pIn);
   return ( it!=selectChannels_.end() || selectAll_);
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

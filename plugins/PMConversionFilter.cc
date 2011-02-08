// -*- C++ -*-
//
// Package:    PMConversionFilter
// Class:      PMConversionFilter
// 
/**\class PMConversionFilter PMConversionFilter.cc Filter/PMConversionFilter/src/PMConversionFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jason Keller
//         Created:  Wed Oct 20 15:27:14 CDT 2010
// $Id: PMConversionFilter.cc,v 1.1 2010/10/25 14:26:50 kellerjd Exp $
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

class PMConversionFilter : public edm::EDFilter {
   public:
      explicit PMConversionFilter(const edm::ParameterSet&);
      ~PMConversionFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      edm::InputTag electronTag_;
      edm::InputTag conversionTag_;
      bool filter_;
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
PMConversionFilter::PMConversionFilter(const edm::ParameterSet& iConfig) :
  electronTag_(iConfig.getParameter<edm::InputTag>("Electrons")),
  conversionTag_(iConfig.getParameter<edm::InputTag>("Conversions")),
  filter_(iConfig.getParameter<bool>("Filter"))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Electron> >();
}


PMConversionFilter::~PMConversionFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PMConversionFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronTag_, electrons);

   Handle<reco::ConversionCollection> conversions;
   iEvent.getByLabel(conversionTag_, conversions);

   std::auto_ptr<std::vector<pat::Electron> > outputCol(new std::vector<pat::Electron>());

   for(std::vector<pat::Electron>::const_iterator electron = electrons->begin(); 
       electron != electrons->end(); ++electron)
   {
     if(electron->shFracInnerHits() < 0.5)
     {
       outputCol->push_back(*electron);
       continue;
     }
     const reco::TrackBaseRef ctfTrack(electron->closestCtfTrackRef());

     bool isGood = true;
     for(reco::ConversionCollection::const_iterator conversion = conversions->begin();
         conversion != conversions->end(); ++conversion)
     {
       const std::vector<reco::TrackBaseRef>& convTracks = conversion->tracks();
       for(unsigned int i = 0; i != convTracks.size(); ++i)
       {
         if(ctfTrack == convTracks[i])
         {
           isGood = false;
           break;
         }
       }
       if(!isGood)
         break;
     }
     if(isGood)
       outputCol->push_back(*electron);
   }  

   iEvent.put(outputCol); 

   return !filter_ || !outputCol->empty();
}

// ------------ method called once each job just before starting event loop  ------------
void 
PMConversionFilter::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void 
PMConversionFilter::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(PMConversionFilter);

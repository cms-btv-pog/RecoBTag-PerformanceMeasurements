// -*- C++ -*-
//
// Package:    TtProcessIDProducer
// Class:      TtProcessIDProducer
// 
/**\class TtProcessIDProducer TtProcessIDProducer.cc RecoBTag/TtProcessIDProducer/src/TtProcessIDProducer.cc

 Description: Prepare process id for writing into edm:Event (originally designed for CSA07 soups)

 Implementation:
     Wraps csa07::csa07ProcessId(const edm::Event &)
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Tue Apr  8 19:00:04 CDT 2008
// $Id$
//
//



#include "RecoBTag/PerformanceMeasurements/interface/TtProcessIDProducer.h"
#include "CSA07EffAnalyser/CSA07EffAnalyser/interface/CSA07ProcessId.h"





TtProcessIDProducer::TtProcessIDProducer(const edm::ParameterSet& iConfig)
{
  //produces<int>("csa07ProcessId");
  produces<int>();
}


TtProcessIDProducer::~TtProcessIDProducer()
{
 
}


void
TtProcessIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::auto_ptr<int> pOut( new int(1) );
   (*pOut) = csa07::csa07ProcessId(iEvent);
   iEvent.put(pOut);

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
}


void 
TtProcessIDProducer::beginJob(const edm::EventSetup&)
{
}


void 
TtProcessIDProducer::endJob() {
}

DEFINE_FWK_MODULE(TtProcessIDProducer);

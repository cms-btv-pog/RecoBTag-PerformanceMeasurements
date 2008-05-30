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

#ifndef PerformanceMeasurementsTtProcessIDProducer
#define PerformanceMeasurementsTtProcessIDProducer


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"





class TtProcessIDProducer : public edm::EDProducer {
   public:
      explicit TtProcessIDProducer(const edm::ParameterSet&);
      ~TtProcessIDProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
};

#endif

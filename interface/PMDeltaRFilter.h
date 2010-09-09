//
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      PMDeltaRFilter
//
/**\class PerformanceMeasurements/PMDeltaRFilter

 Description:

	 Author: Francisco Yumiceva, Fermilab
*/
//
// $Id: PMDeltaRFilter.h,v 1.2 2009/10/03 20:00:33 yumiceva Exp $
//
//


#ifndef PMDeltaRFilter_H
#define PMDeltaRFilter_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"

class PMDeltaRFilter : public edm::EDFilter
{

public:
    explicit PMDeltaRFilter(const edm::ParameterSet &iConfig);
    virtual ~PMDeltaRFilter() {}
    virtual bool filter(edm::Event& iEvent , const edm::EventSetup & iSetup );

private:
    edm::InputTag jets_;
    edm::InputTag muons_;
    double MaxDeltaR_;

};

#endif

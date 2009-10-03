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
// $Id: PMDeltaRFilter.h,v 1.1 2009/09/22 02:45:37 yumiceva Exp $
//
//


#ifndef PMDeltaRFilter_H
#define PMDeltaRFilter_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

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

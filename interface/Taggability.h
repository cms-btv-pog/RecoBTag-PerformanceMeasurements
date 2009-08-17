//
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      Taggability
// 
/**\class PerformanceMeasurements/Taggability

 Description:

	 Author: Francisco Yumiceva, Fermilab
*/
//
// $Id: Taggability.h,v 1.0 2009/07/13 15:13:36 yumiceva Exp $
//
//


#ifndef Taggability_H
#define Taggability_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class Taggability : public edm::EDFilter {

  public:
	explicit Taggability(const edm::ParameterSet &);
	virtual ~Taggability();
	virtual bool filter(edm::Event& , const edm::EventSetup & );

  private:

	edm::InputTag JetCollection_;
	bool useJetCorr_;
	std::string jetCorrLabel_;
	double MinJetPt_;
	double MaxJetEta_;
	int MinNtrksInJet_;
	double MinTrkPtInJet_;
	int MinNjets_;
	edm::InputTag PVCollection_;
	int MinNPV_;
	std::string bTagTrackEventIPTagInfos_;
	
	

};

#endif

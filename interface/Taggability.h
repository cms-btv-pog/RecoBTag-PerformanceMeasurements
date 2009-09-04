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
// $Id: Taggability.h,v 1.2 2009/08/24 20:37:42 yumiceva Exp $
//
//


#ifndef Taggability_H
#define Taggability_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH2F.h"

class Taggability : public edm::EDFilter
{

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
    bool writeHistos_;

    TH2F *h2_in;
    TH2F *h2_out;

};

#endif

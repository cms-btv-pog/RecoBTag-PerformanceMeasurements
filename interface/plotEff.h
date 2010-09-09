#ifndef plotEff_H
#define plotEff_H

/** \class plotEff
 *
 *  Author: Francisco Yumiceva
 */

// FWK
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Formats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

// Root
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "RecoBTag/PerformanceMeasurements/test/S8Tools/S8bPerformance.h"
#include "RecoBTag/PerformanceMeasurements/interface/PFTools.h"

/*
class WorkingPoint{
 public:
	WorkingPoint() {}
  WorkingPoint(edm::InputTag t,  std::string n, double min, double max, std::map< std::string, double> list) : intag_(t), alias_(n), min_(min), max_(max), wpmap_(list) {}

	std::map<std::string, double > list() const {return wpmap_;}
	edm::InputTag inputTag() const {return intag_;}
	std::string alias () const {return alias_;}
	double Minimum() const { return min_;}
	double Maximum() const { return max_;}
	//void print () const ;
 private:
  edm::InputTag intag_;
  std::string alias_;
  double min_;
  double max_;
  std::map< std::string, double > wpmap_;
};
*/


// class declaration


class plotEff : public edm::EDAnalyzer
{


public:
    explicit plotEff(const edm::ParameterSet&);
    ~plotEff();

    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    int TaggedJet(reco::CaloJet calojet, edm::Handle<reco::JetTagCollection > jetTags );
    reco::JetFlavour getMatchedParton(const reco::CaloJet &jet);

private:

    std::string outputFile_;                   // output file
    std::vector<WorkingPoint> wp;
    std::map< std::string, WorkingPoint > wp_map;
    std::string CaloJetCollectionTags_;
    edm::InputTag flavourSourcef;
    edm::Handle<reco::JetFlavourMatchingCollection> theJetPartonMapf;
    std::string bTagTrackEventIPTagInfos_;

    double MinJetPt_;
    double MaxJetEta_;
    int MinNtrksInJet_;
    double MinTrkPtInJet_;
    int MinNjets_;
    bool debug_;
    bool useJetCorr_;
    std::string jetCorrLabel_;
    bool OPbyMistagRate_;
    TFile*  rootFile_;

    std::map< std::string , TH1* > h_map;
    std::map< std::string , S8bPerformance > TaggerPerformances_;
    std::map< std::string , std::string > TaggerAlias_;

};

#endif



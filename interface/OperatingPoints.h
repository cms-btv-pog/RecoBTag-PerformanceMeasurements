#ifndef OperatingPoints_H
#define OperatingPoints_H

/** \class OperatingPoints
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

#include "RecoBTag/PerformanceMeasurements/test/S8Tools/S8bPerformance.h"
#include "RecoBTag/PerformanceMeasurements/interface/PFTools.h"


// class declaration


class OperatingPoints : public edm::EDAnalyzer
{
public:

    explicit OperatingPoints(const edm::ParameterSet&);

    ~OperatingPoints();

    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    int TaggedJet(reco::CaloJet calojet, edm::Handle<reco::JetTagCollection > jetTags );
    reco::JetFlavour getMatchedParton(const reco::CaloJet &jet);

private:

    std::string outputFile_;                   // output file
    std::vector<WorkingPoint> wp;
    std::map<std::string, WorkingPoint> wp_map;
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

    std::map< std::string , S8bPerformance > TaggerPerformances_;
    std::map< std::string , std::string > TaggerAlias_;
};

#endif



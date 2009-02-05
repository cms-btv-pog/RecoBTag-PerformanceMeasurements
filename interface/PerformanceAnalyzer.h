#ifndef PerformanceAnalyzer_h
#define PerformanceAnalyzer_h





/** \class edm::EDAnalyzer PerformanceAnalyzer
 *
 * Analyzer to select jets together with a muon on it.
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: PerformanceAnalyzer.h,v 1.25 2008/12/19 10:17:04 tboccali Exp $
 *
 */

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// vertex stuff
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"


// Root
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

// event class
#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagHistograms.h"

#include "RecoBTag/PerformanceMeasurements/test/S8Tools/S8bPerformance.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"


class WorkingPoint{
 public:
  WorkingPoint(edm::InputTag t,  std::string n,   double c) : it(t), name_(n), cut_(c) {}
  float cut() const {return cut_;}
  edm::InputTag inputTag() const {return it;}
  std::string name () const {return name_;}
  void print () const ;
 private:
  edm::InputTag it;
  std::string name_;
  double cut_;
};




struct ltstr
{
    bool operator()(const edm::RefToBase<reco::Jet> s1, edm::RefToBase<reco::Jet> s2) const
    {
        if (s1.id() != s2.id()) return s1.id()<s2.id();
        return s1.key()< s2.key();
    }
};


// class declaration


class PerformanceAnalyzer : public edm::EDAnalyzer
{


public:
    explicit PerformanceAnalyzer(const edm::ParameterSet&);
    ~PerformanceAnalyzer();

    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    reco::GenJet GetGenJet(reco::CaloJet calojet, reco::GenJetCollection genJetColl);
    SimTrack GetGenTrk(reco::Track atrack, const edm::SimTrackContainer *simTrkColl, const edm::SimVertexContainer *simVtcs);
    int GetMotherId(const edm::SimVertexContainer *simVtxColl, const edm::SimTrackContainer *simTrkColl, SimTrack muonMC);

    int TaggedJet(reco::CaloJet calojet, edm::Handle<reco::JetTagCollection > jetTags );
    int TaggedJet(reco::CaloJet const &, edm::Handle<std::vector<reco::TrackIPTagInfo> > const &);
    
    std::map< std::string, bool > GetBTaggingMap(reco::CaloJet jet, const edm::Event&, double ptrel=0.);
    void FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
                    int JetFlavor, std::map<std::string, bool> aMap, double weight);
    void FillEff(TLorentzVector p4MuJet, int JetFlavor, std::map<std::string, bool> aMap, double weight);
    void FillPtrel(double ptrel, int JetFlavor, std::map<std::string, bool> aMap, double weight);
    void FillPerformance(reco::CaloJet jet, int JetFlavor, const edm::Event&);
    reco::JetFlavour getMatchedParton(const reco::CaloJet &jet);

private:

    // ----------member data ---------------------------
    std::string outputFile_;                   // output file
    std::string recoTrackList_; // collection of tracks
    std::string recoVtxList_;   // collection of vertices
    //    std::vector< std::string > bTaggerList_;    // list of b-tagggers
    //    std::vector< double > bTagCutList_;
    //    std::vector< std::string > moduleLabel_;
    std::vector<WorkingPoint> wp;
    //std::string JetTrackAssociatorTags_;
    std::string MuonCollectionTags_;
    std::string CaloJetCollectionTags_;
    //std::string CorrCaloJetCollectionTags_;
    std::string GenJetCollectionTags_;
    std::string SimTrkCollectionTags_;
    std::string analyzer_;
    std::string fAwayJetTagger;
    std::string flavourMatchOptionf;
    edm::InputTag flavourSourcef;
    bool fWeightHistograms;
    bool fStoreWeightsInNtuple;
    bool fStorePtHat;
    std::map<edm::RefToBase<reco::Jet>, unsigned int, ltstr> flavoursMapf;
    //  edm::Handle<reco::CandMatchMap> theJetPartonMapf;
    edm::Handle<reco::JetFlavourMatchingCollection> theJetPartonMapf;

    //JetFlavourIdentifier jetFlavourIdentifier2_;
    bool StoreTrackProba_;
    double MinJetPt_;
    double MaxJetEta_;
    double MinDeltaR_;
    double MinPtRel_;
    double MinMuonPt_;
    double MaxMuonEta_;
    double MaxMuonChi2_;
    int MinMuonNHits_;

    TFile*  rootFile_;
    TDirectory *topdir;
    bool fverbose;

    TH1I *histcounterf;

    BTagHistograms *EffHistos;
    BTagHistograms *PtrelHistos;
    BTagHistograms *MujetHistos;
    BTagHistograms *AwayjetHistos;
    BTagHistograms *TaggedMujetHistos;
    BTagHistograms *TaggedAwayjetHistos;
    BTagHistograms *MujetHistos_mc;
    BTagHistograms *AwayjetHistos_mc;
    BTagHistograms *TaggedMujetHistos_mc;
    BTagHistograms *TaggedAwayjetHistos_mc;


    S8bPerformance fperformanceTC2trk;
    S8bPerformance fperformanceTC3trk;
    S8bPerformance fperformanceTP;
	S8bPerformance fperformanceJBP;
    S8bPerformance fperformanceSMT;
	S8bPerformance fperformanceSSV;
	S8bPerformance fperformanceCSV;
	
    S8bPerformance fperformanceMTC2trk;
    S8bPerformance fperformanceMTC3trk;
    bool fWritePerformancePlots;

    //edm::InputTag simG4_;
    //double simUnit_;

    //std::map<std::string, TH1*> h;

    TTree *ftree;
    BTagEvent *fS8evt; // system8 container
    int fnselectors;

    int feventcounter;

    std::map< std::string, float > fOPMap;

    std::string bTagTrackEventIPTagInfos_;

    //  edm::ParameterSet trackHistConfig_ ;
    int fbadeventscounter;

    //
    // jet corrections
    //
    std::string jetCorrLabel_;
    bool useJetCorr_;
};


#endif

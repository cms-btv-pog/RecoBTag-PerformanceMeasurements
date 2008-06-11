#ifndef PerformanceAnalyzer_h
#define PerformanceAnalyzer_h





/** \class edm::EDAnalyzer PerformanceAnalyzer
 *  
 * Analyzer to select jets together with a muon on it.
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: PerformanceAnalyzer.h,v 1.17 2008/04/21 22:47:00 yumiceva Exp $
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

//generator level + CLHEP
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
//#include "CLHEP/HepMC/GenEvent.h"
//#include "CLHEP/HepMC/GenVertex.h"
//#include "CLHEP/HepMC/GenParticle.h"

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
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "RecoBTag/MCTools/interface/JetFlavour.h"
#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"

#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackHistory/interface/TrackCategories.h"
#include "SimTracker/TrackHistory/interface/TrackOrigin.h"

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



struct ltstr{
  bool operator()(const edm::RefToBase<reco::Jet> s1, edm::RefToBase<reco::Jet> s2) const
  {
    if (s1.id() != s2.id()) return s1.id()<s2.id();
    return s1.key()< s2.key();
  }
};




TrackAssociatorBase * associatorByHits;

// class declaration

using BTagMCTools::JetFlavour;

class PerformanceAnalyzer : public edm::EDAnalyzer {


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
  std::map< std::string, bool > GetBTaggingMap(reco::CaloJet jet,std::vector<edm::Handle<reco::JetTagCollection > > jetTags_testManyByType, double ptrel=0.);
  void FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
				  int JetFlavor, std::map<std::string, bool> aMap, double weight);
  void FillEff(TLorentzVector p4MuJet, int JetFlavor, std::map<std::string, bool> aMap, double weight);
  void FillPtrel(double ptrel, int JetFlavor, std::map<std::string, bool> aMap, double weight);
  void FillPerformance(reco::CaloJet jet, int JetFlavor, std::vector<edm::Handle<reco::JetTagCollection > > jetTags_testManyByType);
  JetFlavour getMatchedParton(const reco::CaloJet &jet);

private:

  TrackCategories classifier_;
  TrackOrigin tracer_;
  // ----------member data ---------------------------
  std::string outputFile_;                   // output file
  std::string recoTrackList_; // collection of tracks
  std::string recoVtxList_;   // collection of vertices
  std::vector< std::string > bTaggerList_;    // list of b-tagggers
  std::vector< double > bTagCutList_; 
  std::vector< std::string > moduleLabel_;
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

  JetFlavourIdentifier jetFlavourIdentifier_;
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

  TrackAssociatorBase *associatorByChi2;
  TrackAssociatorByHits *associatorByHits;

  std::map< std::string, float > fOPMap;
  
  bool hasLongLived(const TrackHistory &, int) const;

  bool hasPhotonConversion(const TrackHistory &) const;

  std::vector<bool> getTrackCategories( 
    reco::TrackRef, 
    const reco::RecoToSimCollection &, 
    bool bestMatchByMaxValue
  );
  
  bool bTagTrackEventFlag_;
  //  edm::ParameterSet trackHistConfig_ ;
  int fbadeventscounter;  

};


#endif

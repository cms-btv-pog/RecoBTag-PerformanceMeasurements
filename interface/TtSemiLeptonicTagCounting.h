// -*- C++ -*-
//
// Package:    BTagPerformanceMeasurementFromTop
// Class:      TtSemiLeptonicTagCounting
// 
/**\class TtSemiLeptonicTagCounting TtSemiLeptonicTagCounting.cc TopQuarkAnalysis/BTagPerformanceMeasurementFromTop/src/TtSemiLeptonicTagCounting.cc

 Description: Tag counting method for b-,c- and light tagging efficiency measurement with ttbar semileptonic events

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Fri Jun 29 14:53:10 CDT 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoBTag/PerformanceMeasurements/interface/RooGKCounter.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
//#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/METReco/interface/CaloMET.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "AnalysisDataFormats/TopObjects/interface/TopLepton.h"
#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMET.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/JetTagFwd.h"

#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfoFwd.h"
#include "DataFormats/BTauReco/interface/JetTracksAssociation.h"

#include "DataFormats/Common/interface/RefToBase.h"

#include <fstream>

//
// class decleration
//

class TtSemiLeptonicTagCounting : public edm::EDAnalyzer
{
public:
  explicit TtSemiLeptonicTagCounting(const edm::ParameterSet&);
  ~TtSemiLeptonicTagCounting();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  /// histogram list
  //TList histograms_;
  /// algorithm 
  //TopBTagSelectorWorker algo_;
  
  int findJetMatch( const reco::Jet & theJet, const std::vector<reco::GenJet> & jets, std::vector<int> & jetIndex, double dRCut );
  
  RooGKCounter eventCounter;
  RooGKCounter selectedEvents;

  map< int, map< int, map< int, int > > > Fijk;
  map< int, map< int, map< int, map< int, int > > > > Fxijk;

  map< double, map< int, int> > Nn; // number of events with n tagged jets

  double dRCut;
  int nCaloJetsLow;
  int nLeptonLow;

  int _interFreq;
  string _dataType;
  bool _lookForGenJetMatch;
  string _jetSource, _electronSource, _muonSource, _METSource, _jetTagSource;
  string _genJetSource;
  string _outputFileName;

  double dLow, dHigh, dStep;

  RooGKCounter nOfPassedLJets;
  RooGKCounter nOfPassedBJets;
  RooGKCounter nOfPassedCJets;
  RooGKCounter nOfPassedUnknownJets;
  map<double,RooGKCounter> nOfTaggedPassedLJets;
  map<double,RooGKCounter> nOfTaggedPassedBJets;
  map<double,RooGKCounter> nOfTaggedPassedCJets;
  map<double,RooGKCounter> nOfTaggedPassedUnknownJets;

  RooGKCounter interCounter; // counts intermediate results

  ofstream _outputFile; 
  ofstream _outputFileTable; 
  ofstream _outputFileMC; 
  ofstream _outputFileMC2; 
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

// -*- C++ -*-
//
// Package:    BTagPerformanceMeasurementFromTop
// Class:      TtTagConsistency
// 
/**\class TtTagConsistency TtTagConsistency.cc TopQuarkAnalysis/BTagPerformanceMeasurementFromTop/src/TtTagConsistency.cc

 Description: Tag counting method for b-,c- and light tagging efficiency measurement with ttbar semileptonic events

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Fri Jun 29 14:53:10 CDT 2007
// $Id: TtTagConsistency.h,v 1.1.2.1 2008/04/04 23:54:48 kukartse Exp $
//
//

#ifndef PerformanceMeasurementsTtTagConsistency
#define PerformanceMeasurementsTtTagConsistency

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

class TtTagConsistency : public edm::EDAnalyzer
{
public:
  explicit TtTagConsistency(const edm::ParameterSet&);
  ~TtTagConsistency();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  int findJetMatch( const reco::Jet & theJet, const std::vector<reco::GenJet> & jets, std::vector<int> & jetIndex, double dRCut );
  
  RooGKCounter eventCounter;
  RooGKCounter selectedEvents;

  map< int, map< int, map< int, double > > > Fijk;
  map< int, map< int, map< int, map< int, double > > > > Fxijk;

  map< double, map< int, double> > Nn; // number of events with n tagged jets
  map< double, map< int, double> > Nn_2; // number of events with n tagged jets, second tagger
  map< double, map< int, double> > Nn_3; // number of events with n tagged jets, third tagger

  double dRCut;
  int nCaloJetsLow;
  int nLeptonLow;

  int _interFreq;
  string _dataType;
  bool _lookForGenJetMatch;
  string _jetSource, _electronSource, _muonSource, _METSource;
  string _jetTagSource, _jetTagSource_2, _jetTagSource_3;
  string _genJetSource;
  string _outputFileName;

  double dLow, dHigh, dStep;
  double dLow_2, dHigh_2, dStep_2;
  double dLow_3, dHigh_3, dStep_3;

  RooGKCounter nOfPassedLJets;
  RooGKCounter nOfPassedBJets;
  RooGKCounter nOfPassedCJets;
  RooGKCounter nOfPassedUnknownJets;
  map<double,RooGKCounter> nOfTaggedPassedLJets;
  map<double,RooGKCounter> nOfTaggedPassedBJets;
  map<double,RooGKCounter> nOfTaggedPassedCJets;
  map<double,RooGKCounter> nOfTaggedPassedUnknownJets;
  map<double,RooGKCounter> nOfTaggedPassedLJets_2;
  map<double,RooGKCounter> nOfTaggedPassedBJets_2;
  map<double,RooGKCounter> nOfTaggedPassedCJets_2;
  map<double,RooGKCounter> nOfTaggedPassedUnknownJets_2;
  map<double,RooGKCounter> nOfTaggedPassedLJets_3;
  map<double,RooGKCounter> nOfTaggedPassedBJets_3;
  map<double,RooGKCounter> nOfTaggedPassedCJets_3;
  map<double,RooGKCounter> nOfTaggedPassedUnknownJets_3;

  RooGKCounter interCounter; // counts intermediate results

  ofstream _outputFile; 

  ofstream _outputFileTable; 
  ofstream _outputFileTable_2;
  ofstream _outputFileTable_3;

  ofstream _outputFileMC; 
  ofstream _outputFileMC2; 
  
};

#endif

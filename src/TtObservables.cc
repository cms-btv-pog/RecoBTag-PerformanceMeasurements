// -*- C++ -*-
//
// Package:    BTagPerformanceMeasurementFromTop
// Class:      TtObservables
// 
/**\class TtObservables TtObservables.cc TopQuarkAnalysis/BTagPerformanceMeasurementFromTop/src/TtObservables.cc

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

#include "RecoBTag/PerformanceMeasurements/interface/TtObservables.h"

TtObservables::TtObservables(const edm::ParameterSet& iConfig)
//: histograms_(), 
//algo_( 0, histograms_ )
{
  
  //now do what ever initialization is needed

  //algo_ . begin( iConfig );
    
  dRCut        = iConfig.getParameter<double>("dRGenJetMatch");
  nCaloJetsLow = iConfig.getParameter<int>("nCaloJetsLow");
  nLeptonLow   = iConfig.getParameter<int>("nLeptonLow");

  dLow         = iConfig.getParameter<double>("discriminatorLow");
  dHigh        = iConfig.getParameter<double>("discriminatorHigh");
  dStep        = iConfig.getParameter<double>("discriminatorStep");

  _interFreq   = iConfig.getParameter<int>("intermediateResultFrequency");
  _dataType    = iConfig.getParameter<string>("dataType");
  _lookForGenJetMatch    = iConfig.getParameter<bool>("lookForGenJetMatch");

  _jetSource    = iConfig.getParameter<string>("jetSource");
  _electronSource    = iConfig.getParameter<string>("electronSource");
  _muonSource    = iConfig.getParameter<string>("muonSource");
  _METSource    = iConfig.getParameter<string>("METSource");
  _jetTagSource    = iConfig.getParameter<string>("jetTagSource");

  _genJetSource    = iConfig.getParameter<string>("genJetSource");
  _outputFileName    = iConfig.getParameter<string>("outputFileName");

  //Histograms
  string _rootFileName = _outputFileName + ".root";
  file_ = new TFile( _rootFileName . c_str(),"RECREATE");

  pt_jet = new TH1D("pt_jet","Jet p_{T}",250,0.,250.);
  pt_jet->Sumw2();
  pt_el = new TH1D("pt_el","Electron p_{T}",250,0.,250.);
  pt_el->Sumw2();
  pt_mu = new TH1D("pt_mu","Muon p_{T}",250,0.,250.);
  pt_mu->Sumw2();
  pt_MET = new TH1D("pt_MET","MET p_{T}",250,0.,250.);
  pt_MET->Sumw2();

  et_jet = new TH1D("et_jet","Jet E_{T}",250,0.,250.);
  et_jet->Sumw2();
  et_el = new TH1D("et_el","Electron E_{T}",250,0.,250.);
  et_el->Sumw2();
  et_mu = new TH1D("et_mu","Muon E_{T}",250,0.,250.);
  et_mu->Sumw2();
  et_MET = new TH1D("et_MET","MET E_{T}",250,0.,250.);
  et_MET->Sumw2();

  eta_jet = new TH1D("eta_jet","Jet #eta",200,-10.,10.);
  eta_jet->Sumw2();
  eta_el = new TH1D("eta_el","Electron #eta",200,-10.,10.);
  eta_el->Sumw2();
  eta_mu = new TH1D("eta_mu","Muon #eta",200,-10.,10.);
  eta_mu->Sumw2();
  eta_MET = new TH1D("eta_MET","MET #eta",200,-10.,10.);
  eta_MET->Sumw2();

  n_jet = new TH1I("n_jet","Number of jets",100,0.,100.);
  n_el = new TH1I("n_el","Number of electrons",100,0.,100.);
  n_mu = new TH1I("n_mu","Number of muons",100,0.,100.);

}


TtObservables::~TtObservables()
{
 
  file_->Write();
  file_->Close();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TtObservables::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  
  using namespace reco;
  using namespace std;
  
  Handle< vector< TopJet > > caloJets;
  Handle< vector< TopElectron > > electrons;
  Handle< vector< TopMuon > > muons;
  Handle< vector< TopMET > > METs;

  iEvent . getByLabel( _jetSource, caloJets );
  iEvent . getByLabel( _electronSource, electrons );
  iEvent . getByLabel( _muonSource, muons );
  iEvent . getByLabel( _METSource, METs );
  
  eventCounter . count();

  // loop over jets
  //
  vector<TopJet>::const_iterator cjet;
  RooGKCounter nOfGoodJets;
  for ( cjet = caloJets -> begin(); cjet != caloJets -> end(); cjet++ )  // loop over calo jets
    {
      
      nOfGoodJets . count();

      pt_jet -> Fill( cjet -> pt() );
      et_jet -> Fill( cjet -> et() );
      eta_jet -> Fill( cjet -> eta() );

      
    }
  
  
  // loop over electrons
  vector<TopElectron>::const_iterator el;
  RooGKCounter nOfGoodElectrons( "" );
  for ( el = electrons -> begin(); el != electrons -> end(); el++)
    {
	  
      nOfGoodElectrons . count();
      pt_el -> Fill( el -> pt() );
      et_el -> Fill( el -> et() );
      eta_el -> Fill( el -> eta() );
      
    }
  
  
  // loop over muons
  vector<TopMuon>::const_iterator mu;
  RooGKCounter nOfGoodMuons( "" );
  for ( mu = muons -> begin(); mu != muons -> end(); mu++)
    {
      
      nOfGoodMuons . count();
      pt_mu -> Fill( mu -> pt() );
      et_mu -> Fill( mu -> et() );
      eta_mu -> Fill( mu -> eta() );

    }

  // loop over METs
  vector<TopMET>::const_iterator met;
  RooGKCounter nOfGoodMETs( "" );
  for ( met = METs -> begin(); met != METs -> end(); met++)
    {

      nOfGoodMETs . count();
      pt_MET -> Fill( met -> pt() );
      et_MET -> Fill( met -> et() );
      eta_MET -> Fill( met -> eta() );

    }

  n_jet -> Fill( nOfGoodJets . getCount() );
  n_el -> Fill( nOfGoodElectrons . getCount() );
  n_mu -> Fill( nOfGoodMuons . getCount() );

}


// ------------ method called once each job just before starting event loop  ------------
void 
TtObservables::beginJob(const edm::EventSetup&)
{

  cout << "TtObservables started..." << endl;

  eventCounter . setMessage( "=== Input event counter: " );
  eventCounter . setDivider( 100 );
  eventCounter . setPrintCount( true );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TtObservables::endJob() {

}

//define this as a plug-in
//DEFINE_FWK_MODULE(TtObservables);

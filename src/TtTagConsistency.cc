// -*- C++ -*-
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      TtTagConsistency
// 
/**\class TtTagConsistency TtTagConsistency.cc RecoBTag/PerformanceMeasurements/src/TtTagConsistency.cc

 Description: Tag counting method for b-,c- and light tagging efficiency measurement with ttbar semileptonic events

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Fri Jun 29 14:53:10 CDT 2007
// $Id: TtTagConsistency.cc,v 1.1.2.4 2008/06/02 15:14:38 kukartse Exp $
//
//

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistency.h"

TtTagConsistency::TtTagConsistency(const edm::ParameterSet& iConfig)
{
  nCaloJets_min         = 4;
  nLepton_min           = 1;
  jet_pt_min            = 25.0     ;
  jet_eta_max           = 2.4      ;
  muon_pt_min           = 20.0     ;
  muon_eta_max          = 2.4      ;
  muon_trackIso_max     = 300000.0 ;
  muon_caloIso_max      = 500000.0 ;
  electron_pt_min       = 20.0     ;
  electron_eta_max      = 2.4      ;
  electron_trackIso_max = 300000.0 ;
  electron_caloIso_max  = 600000.0 ;
  met_et_min            = 25.0     ;
  d_min         = 0.0;
  d_max         = 10.0;
  d_step        = 0.1;
  d_min_2         = 0.0;
  d_max_2         = 10.0;
  d_step_2        = 0.1;
  d_min_3         = 0.0;
  d_max_3         = 1.0;
  d_step_3        = 0.01;


  nCaloJets_min         = iConfig.getParameter<int>("nCaloJets_min");
  nLepton_min           = iConfig.getParameter<int>("nLepton_min");

  jet_pt_min            = iConfig.getParameter<double>("jet_pt_min");
  jet_eta_max           = iConfig.getParameter<double>("jet_eta_max");
  muon_pt_min           = iConfig.getParameter<double>("muon_pt_min");
  muon_eta_max          = iConfig.getParameter<double>("muon_eta_max");
  muon_trackIso_max     = iConfig.getParameter<double>("muon_trackIso_max");
  muon_caloIso_max      = iConfig.getParameter<double>("muon_caloIso_max");
  electron_pt_min       = iConfig.getParameter<double>("electron_pt_min");
  electron_eta_max      = iConfig.getParameter<double>("electron_eta_max");
  electron_trackIso_max = iConfig.getParameter<double>("electron_trackIso_max");
  electron_caloIso_max  = iConfig.getParameter<double>("electron_caloIso_max");
  met_et_min            = iConfig.getParameter<double>("met_et_min");

  d_min         = iConfig.getParameter<double>("discriminator_min");
  d_max        = iConfig.getParameter<double>("discriminator_max");
  d_step        = iConfig.getParameter<double>("discriminator_step");
  d_min_2         = iConfig.getParameter<double>("discriminator_min_2");
  d_max_2        = iConfig.getParameter<double>("discriminator_max_2");
  d_step_2        = iConfig.getParameter<double>("discriminator_step_2");
  d_min_3         = iConfig.getParameter<double>("discriminator_min_3");
  d_max_3        = iConfig.getParameter<double>("discriminator_max_3");
  d_step_3        = iConfig.getParameter<double>("discriminator_step_3");
  
  _interFreq   = iConfig.getParameter<int>("intermediateResultFrequency");
  _dataType    = iConfig.getParameter<string>("dataType");

  _jetSource    = iConfig.getParameter<string>("jetSource");
  _electronSource    = iConfig.getParameter<string>("electronSource");
  _muonSource    = iConfig.getParameter<string>("muonSource");
  _METSource    = iConfig.getParameter<string>("METSource");
  _jetTagSource    = iConfig.getParameter<string>("jetTagSource");
  _jetTagSource_2  = iConfig.getParameter<string>("jetTagSource_2");
  _jetTagSource_3  = iConfig.getParameter<string>("jetTagSource_3");

  _outputFileName    = iConfig.getParameter<string>("outputFileName");
}



TtTagConsistency::~TtTagConsistency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TtTagConsistency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  Handle< double> weightHandle;
  double weight = 1.0;

  try{
    iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
    weight = * weightHandle;
  }
  catch( edm::Exception ) {
    cout << "Warning: csa07 soup weight is not available... continuing with weight=1.0" << endl;
  }

  eventCounter . count();
  eventCounter . incrementDouble( weight );

  // loop over jets
  //
  vector<TopJet>::const_iterator cjet;
  RooGKCounter nOfGoodJets;
  RooGKCounter nOfGoodBJets;
  RooGKCounter nOfGoodCJets;
  RooGKCounter nOfGoodLJets;
  RooGKCounter nOfGoodGJets;
  RooGKCounter nOfGoodUnknownJets;
  RooGKCounter nOfGoodOtherJets; // jet that passes selection and is neither quark nor gluon nor id=0
  map<double,RooGKCounter> nTaggedJets;
  map<double,RooGKCounter> nTaggedJets_2;
  map<double,RooGKCounter> nTaggedJets_3;
  map<double,RooGKCounter> nOfTaggedLJets;
  map<double,RooGKCounter> nOfTaggedBJets;
  map<double,RooGKCounter> nOfTaggedCJets;
  map<double,RooGKCounter> nOfTaggedUnknownJets;
  map<double,RooGKCounter> nOfTaggedLJets_2;
  map<double,RooGKCounter> nOfTaggedBJets_2;
  map<double,RooGKCounter> nOfTaggedCJets_2;
  map<double,RooGKCounter> nOfTaggedUnknownJets_2;
  map<double,RooGKCounter> nOfTaggedLJets_3;
  map<double,RooGKCounter> nOfTaggedBJets_3;
  map<double,RooGKCounter> nOfTaggedCJets_3;
  map<double,RooGKCounter> nOfTaggedUnknownJets_3;

  vector<int> goodCaloJetIndex; // indices of the jets that pass quality cuts
  vector<int> goodCaloJetFlavorID; // flavor of the good jets

  // loop over calo jets
  for ( cjet = caloJets -> begin(); cjet != caloJets -> end(); cjet++ ){
    if ( 
	//	  ( fabs( (*cjet) . eta() ) < 2.4 ) &&
	//	  ( (*cjet) . pt() > 25.0 ) &&
	( fabs( (*cjet) . eta() ) < jet_eta_max ) &&
	( (*cjet) . pt() > jet_pt_min ) 
	){
      
      nOfGoodJets . count();
      
      map<double,bool> isTagged;
      map<double,bool> isTagged_2;
      map<double,bool> isTagged_3;
      for(double _d = d_min; _d <= d_max; _d += 0.5) isTagged[_d] = false;
      for(double _d = d_min_2; _d <= d_max_2; _d += 0.5) isTagged_2[_d] = false;
      for(double _d = d_min_3; _d <= d_max_3; _d += 0.5) isTagged_3[_d] = false;
      
      //cjet -> dumpBTagLabels();
      
      double jetBDiscr  = cjet -> getBDiscriminator( _jetTagSource );
      double jetBDiscr2 = cjet -> getBDiscriminator( _jetTagSource_2 );
      double jetBDiscr3 = cjet -> getBDiscriminator( _jetTagSource_3 );
      
      // loop over discriminator values
      for(double _d = d_min; _d < d_max+0.0001; _d += d_step){
	if ( jetBDiscr > _d ){
	  nTaggedJets[_d] . count(); // check the discriminator
	  isTagged[_d] = true;
	}
      }
      for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2){
	if ( jetBDiscr2 > _d ){
	  nTaggedJets_2[_d] . count(); // check the discriminator
	  isTagged_2[_d] = true;
	}
      }
      for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3){
	if ( jetBDiscr3 > _d ){
	  nTaggedJets_3[_d] . count(); // check the discriminator
	  isTagged_3[_d] = true;
	}
      }
      
      
      // ---- jet flavor (MC only!!!)
      if ( _dataType == "MC" ){
	int jetFlavor = cjet -> getPartonFlavour();
	
	goodCaloJetFlavorID . push_back( jetFlavor );
	
	if ( jetFlavor == 5 ){
	  //nOfGoodJets . count();
	  nOfGoodBJets . count();
	  // loop over discriminator values
	  for(double _d = d_min; _d < d_max+0.0001; _d += d_step){
	    if ( isTagged[_d] ) nOfTaggedBJets[_d] . count();
	  }
	  for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2){
	    if ( isTagged_2[_d] ) nOfTaggedBJets_2[_d] . count();
	  }
	  for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3){
	    if ( isTagged_3[_d] ) nOfTaggedBJets_3[_d] . count();
	  }
	}
	else if ( jetFlavor == 4 ){
	  //nOfGoodJets . count();
	  nOfGoodCJets . count();
	  for(double _d = d_min; _d < d_max+0.0001; _d += d_step){
	    if ( isTagged[_d] ) nOfTaggedCJets[_d] . count();
	  }
	  for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2){
	    if ( isTagged_2[_d] ) nOfTaggedCJets_2[_d] . count();
	  }
	  for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3){
	    if ( isTagged_3[_d] ) nOfTaggedCJets_3[_d] . count();
	  }
	}
	else if ( jetFlavor == 1 || jetFlavor == 2 || jetFlavor == 3 || jetFlavor == 21 ){
	  //nOfGoodJets . count();
	  nOfGoodLJets . count();
	  for(double _d = d_min; _d < d_max+0.0001; _d += d_step){
	    if ( isTagged[_d] ) nOfTaggedLJets[_d] . count();
	  }
	  for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2){
	    if ( isTagged_2[_d] ) nOfTaggedLJets_2[_d] . count();
	  }
	  for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3){
	    if ( isTagged_3[_d] ) nOfTaggedLJets_3[_d] . count();
	  }
	}
	else if ( jetFlavor == 0 ){
	  //nOfGoodJets . count();
	  nOfGoodUnknownJets . count();
	  for(double _d = d_min; _d < d_max+0.0001; _d += d_step){
	    if ( isTagged[_d] ) nOfTaggedUnknownJets[_d] . count();
	  }
	  for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2){
	    if ( isTagged_2[_d] ) nOfTaggedUnknownJets_2[_d] . count();
	  }
	  for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3){
	    if ( isTagged_3[_d] ) nOfTaggedUnknownJets_3[_d] . count();
	  }
	}
	else nOfGoodOtherJets . count();
      }
    }
  }
  
  // loop over electrons
  vector<TopElectron>::const_iterator el;
  RooGKCounter nOfGoodElectrons( "" );
  for ( el = electrons -> begin(); el != electrons -> end(); el++){
    if (
	//	  (*el) . pt() > 20.0  &&
	//	  fabs( (*el) . eta() ) < 2.4  
	(*el) . pt() > electron_pt_min  &&
	fabs( (*el) . eta() ) < electron_eta_max &&
	(*el).getTrackIso() < electron_trackIso_max &&
	(*el).getCaloIso() < electron_caloIso_max
	){
      nOfGoodElectrons . count();
    }
  }
  
  
  // loop over muons
  vector<TopMuon>::const_iterator mu;
  RooGKCounter nOfGoodMuons( "" );
  for ( mu = muons -> begin(); mu != muons -> end(); mu++){
    if (
	//	  (*mu) . pt() > 20.0 &&
	//	  fabs( (*mu) . eta() ) < 2.4 
	(*mu) . pt() > muon_pt_min &&
	fabs( (*mu) . eta() ) < muon_eta_max &&
	(*mu).getTrackIso() < muon_trackIso_max &&
	(*mu).getCaloIso() < muon_caloIso_max
	){
      nOfGoodMuons . count();
    }
  }
  
  // loop over METs
  vector<TopMET>::const_iterator met;
  RooGKCounter nOfGoodMETs( "" );
  for ( met = METs -> begin(); met != METs -> end(); met++){
    if ( 
	//	  (*met) . et() > 25.0
	(*met) . et() > met_et_min
	){
      nOfGoodMETs . count();
    }
  }
  
  
  // ====== SELECTION =================
  //
  //
  if ( 
      (int)nOfGoodJets . getCount() >= nCaloJets_min &&
      ( (int)nOfGoodElectrons . getCount() + (int)nOfGoodMuons . getCount() >= nLepton_min ) &&
      nOfGoodMETs . getCount() > 0 )
    {
      
      selectedEvents . count();
      selectedEvents . incrementDouble( weight );

      nOfPassedBJets . incrementDouble( (double)nOfGoodBJets.getCount()*weight );
      nOfPassedCJets . incrementDouble( (double)nOfGoodCJets.getCount()*weight );
      nOfPassedLJets . incrementDouble( (double)nOfGoodLJets.getCount()*weight );
      nOfPassedUnknownJets . incrementDouble( (double)nOfGoodUnknownJets.getCount()*weight );

      Fijk[ (int)nOfGoodBJets.getCount() ][ (int)nOfGoodCJets.getCount() ][ (int)nOfGoodLJets.getCount() ]+=weight;
      Fxijk[ (int)nOfGoodUnknownJets.getCount() ][ (int)nOfGoodBJets.getCount() ][ (int)nOfGoodCJets.getCount() ][ (int)nOfGoodLJets.getCount() ]+=weight;

      for(double _d = d_min; _d < d_max+0.0001; _d += d_step) // loop over discriminator values
	{
	  nOfTaggedPassedBJets[_d] . incrementDouble( (double)nOfTaggedBJets[_d] . getCount()*weight );
	  nOfTaggedPassedCJets[_d] . incrementDouble( (double)nOfTaggedCJets[_d] . getCount()*weight );
	  nOfTaggedPassedLJets[_d] . incrementDouble( (double)nOfTaggedLJets[_d] . getCount()*weight );
	  nOfTaggedPassedUnknownJets[_d] . incrementDouble( (double)nOfTaggedUnknownJets[_d] . getCount()*weight );

	  Nn[_d][ nTaggedJets[_d] . getCount() ]+=weight;
	}
      for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2) // loop over discriminator values
	{
	  nOfTaggedPassedBJets_2[_d] . incrementDouble( (double)nOfTaggedBJets_2[_d] . getCount()*weight );
	  nOfTaggedPassedCJets_2[_d] . incrementDouble( (double)nOfTaggedCJets_2[_d] . getCount()*weight );
	  nOfTaggedPassedLJets_2[_d] . incrementDouble( (double)nOfTaggedLJets_2[_d] . getCount()*weight );
	  nOfTaggedPassedUnknownJets_2[_d] . incrementDouble( (double)nOfTaggedUnknownJets_2[_d] . getCount()*weight );
	  
	  Nn_2[_d][ nTaggedJets_2[_d] . getCount() ]+=weight;
	}
      for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3) // loop over discriminator values
	{
	  nOfTaggedPassedBJets_3[_d] . incrementDouble( (double)nOfTaggedBJets_3[_d] . getCount()*weight );
	  nOfTaggedPassedCJets_3[_d] . incrementDouble( (double)nOfTaggedCJets_3[_d] . getCount()*weight );
	  nOfTaggedPassedLJets_3[_d] . incrementDouble( (double)nOfTaggedLJets_3[_d] . getCount()*weight );
	  nOfTaggedPassedUnknownJets_3[_d] . incrementDouble( (double)nOfTaggedUnknownJets_3[_d] . getCount()*weight );
	  
	  Nn_3[_d][ nTaggedJets_3[_d] . getCount() ]+=weight;
	}
    }  
}


// ------------ method called once each job just before starting event loop  ------------
void 
TtTagConsistency::beginJob(const edm::EventSetup&)
{

  eventCounter . setMessage( "=== Preselected event counter: " );
  eventCounter . setDivider( 100 );
  eventCounter . setPrintCount( true );

  cout << "Opening output files..." << endl;
  string _logFileName = _outputFileName + ".log";

  string _tableFileName = _outputFileName + "_" + _jetTagSource + ".tab";
  string _tableFileName_2 = _outputFileName + "_" + _jetTagSource_2 + ".tab";
  string _tableFileName_3 = _outputFileName + "_" + _jetTagSource_3 + ".tab";

  string _mcFileName = _outputFileName + ".mc";
  string _mcFileName2 = _outputFileName + ".mc2";
  _outputFile . open( _logFileName . c_str() );

  _outputFileTable . open( _tableFileName . c_str() );
 _outputFileTable_2 . open( _tableFileName_2 . c_str() );
  _outputFileTable_3 . open( _tableFileName_3 . c_str() );

  // FIXME: obsolete .mc
  //_outputFileMC . open( _mcFileName . c_str() );
  _outputFileMC2 . open( _mcFileName2 . c_str() );

  for(double _d = d_min; _d < d_max+0.0001; _d += d_step){ // loop over discriminator values
    for (int i=0; i<=4; i++){
      Nn[_d][i] = 0.0;
    }
  }
  for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2){ // loop over discriminator values
    for (int i=0; i<=4; i++){
      Nn_2[_d][i] = 0.0;
    }
  }
  for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3){ // loop over discriminator values
    for (int i=0; i<=4; i++){
      Nn_3[_d][i] = 0.0;
    }
  }
}

int TtTagConsistency::findJetMatch( const reco::Jet & theJet, const std::vector<reco::GenJet> & jets, std::vector<int> & jetIndex, double dRCut ) {

  int currentIndex = -1; // index of the current match, -1 if none

  double currentdR = 1.0e10;

  double theJetEta = theJet . eta();
  double theJetPhi = theJet . phi();

  vector<int>::const_iterator index;
  
  for ( index = jetIndex . begin(); index != jetIndex . end(); index++ )
    {
      
      double dR = sqrt( ( theJetEta - jets[(*index)] . eta() ) * ( theJetEta - jets[(*index)] . eta() )
			+ ( theJetPhi - jets[(*index)] . phi() ) * ( theJetPhi - jets[(*index)] . phi() ));

      if ( dR < dRCut && dR < currentdR )
	{
	  currentIndex = (*index);
	  currentdR = dR;
	}
    }

  //cout << "=============> " << endl;
  //cout << "matched gen jet index: " << currentIndex << endl;
  //cout << "pt = " << jets[currentIndex] . pt() << endl;

  return currentIndex;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TtTagConsistency::endJob() {

  _outputFile << "============= TtTagConsistency summary ================" << endl;
  _outputFile << "Processed events:  " << eventCounter . getCount() << endl;
  _outputFile << "Processed events (weighed):  " << eventCounter . getCountDouble() << endl;
  _outputFile << "Selected events:  " << selectedEvents . getCount() << endl;
  _outputFile << "Selected events (weighed):  " << selectedEvents . getCountDouble() << endl;

  map< int, map< int, map< int, map< int, double > > > >::const_iterator it_x;  
  map< int, map< int, map< int, double > > >::const_iterator it_i;  
  map< int, map< int, double > >::const_iterator it_j;  
  map< int, double >::const_iterator it_k;  

  //===> table
  char buf[10];
  _outputFileTable << "     discr       N_0       N_1       N_2       N_3       N_4       N_b N_btagged       N_c N_ctagged       N_l N_ltagged       N_x N_xtagged" << endl;
  for(double _d = d_min; _d < d_max+0.0001; _d += d_step) // loop over discriminator values
    {
      sprintf( buf, "%10.2f", _d);
      _outputFileTable << buf;
      sprintf( buf, "%10.2f", Nn[_d][0]);
      _outputFileTable << buf;
      sprintf( buf, "%10.2f", Nn[_d][1]);
      _outputFileTable << buf;
      sprintf( buf, "%10.2f", Nn[_d][2]);
      _outputFileTable << buf;
      sprintf( buf, "%10.2f", Nn[_d][3]);
      _outputFileTable << buf;
      sprintf( buf, "%10.2f", Nn[_d][4]);
      _outputFileTable << buf;

      if ( _dataType == "MC" )
	{
	  sprintf( buf, "%10.2f", nOfPassedBJets . getCountDouble() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedBJets[_d] . getCountDouble() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10.2f", nOfPassedCJets . getCountDouble() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedCJets[_d] . getCountDouble() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10.2f", nOfPassedLJets . getCountDouble() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedLJets[_d] . getCountDouble() );
	  _outputFileTable << buf;

	  sprintf( buf, "%10.2f", nOfPassedUnknownJets . getCountDouble() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedUnknownJets[_d] . getCountDouble() );
	  _outputFileTable << buf << endl;

	}
      else
	{
	  for ( int i = 0; i < 8; i++ )
	    {
	      sprintf( buf, "%10.2f", 0.0 );
	      _outputFileTable << buf;
	    }
	  _outputFileTable << endl;
	}
    }

 // table for the second tagger
  _outputFileTable_2 << "     discr       N_0       N_1       N_2       N_3       N_4       N_b N_btagged       N_c N_ctagged       N_l N_ltagged       N_x N_xtagged" << endl;
  for(double _d = d_min_2; _d < d_max_2+0.0001; _d += d_step_2) // loop over discriminator values
    {
      sprintf( buf, "%10.2f", _d);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10.2f", Nn_2[_d][0]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10.2f", Nn_2[_d][1]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10.2f", Nn_2[_d][2]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10.2f", Nn_2[_d][3]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10.2f", Nn_2[_d][4]);
      _outputFileTable_2 << buf;

      if ( _dataType == "MC" )
	{
	  sprintf( buf, "%10.2f", nOfPassedBJets . getCountDouble() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedBJets_2[_d] . getCountDouble() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10.2f", nOfPassedCJets . getCountDouble() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedCJets_2[_d] . getCountDouble() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10.2f", nOfPassedLJets . getCountDouble() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedLJets_2[_d] . getCountDouble() );
	  _outputFileTable_2 << buf;

	  sprintf( buf, "%10.2f", nOfPassedUnknownJets . getCountDouble() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedUnknownJets_2[_d] . getCountDouble() );
	  _outputFileTable_2 << buf << endl;

	}
      else
	{
	  for ( int i = 0; i < 8; i++ )
	    {
	      sprintf( buf, "%10.2f", 0.0 );
	      _outputFileTable_2 << buf;
	    }
	  _outputFileTable_2 << endl;
	}
    }

  // table for the third tagger
  _outputFileTable_3 << "     discr       N_0       N_1       N_2       N_3       N_4       N_b N_btagged       N_c N_ctagged       N_l N_ltagged       N_x N_xtagged" << endl;
  for(double _d = d_min_3; _d < d_max_3+0.0001; _d += d_step_3) // loop over discriminator values
    {
      sprintf( buf, "%10.2f", _d);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10.2f", Nn_3[_d][0]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10.2f", Nn_3[_d][1]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10.2f", Nn_3[_d][2]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10.2f", Nn_3[_d][3]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10.2f", Nn_3[_d][4]);
      _outputFileTable_3 << buf;

      if ( _dataType == "MC" )
	{
	  sprintf( buf, "%10.2f", nOfPassedBJets . getCountDouble() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedBJets_3[_d] . getCountDouble() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10.2f", nOfPassedCJets . getCountDouble() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedCJets_3[_d] . getCountDouble() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10.2f", nOfPassedLJets . getCountDouble() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedLJets_3[_d] . getCountDouble() );
	  _outputFileTable_3 << buf;

	  sprintf( buf, "%10.2f", nOfPassedUnknownJets . getCountDouble() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10.2f", nOfTaggedPassedUnknownJets_3[_d] . getCountDouble() );
	  _outputFileTable_3 << buf << endl;

	}
      else
	{
	  for ( int i = 0; i < 8; i++ )
	    {
	      sprintf( buf, "%10.2f", 0.0 );
	      _outputFileTable_3 << buf;
	    }
	  _outputFileTable_3 << endl;
	}
    }

  
  if ( _dataType == "MC" )
    {
      /* FIXME: obsolete
      for ( it_i = Fijk . begin(); it_i != Fijk . end(); it_i++ )
	{
	  map< int, map< int, double > > map2 = it_i -> second;        
	  for ( it_j = map2 . begin(); it_j != map2 . end(); it_j++ )
	    {
	      map< int, double > map1 = it_j -> second;        
	      for ( it_k = map1 . begin(); it_k != map1 . end(); it_k++ )
		{
		  int ii = it_i -> first;
		  int jj = it_j -> first;
		  int kk = it_k -> first;
		  _outputFileMC << "F_" << ii << "_" << jj << "_" << kk << " = " << Fijk[ii][jj][kk] << endl;
		}
	    }
	}
      */

      // dataset flavor structure with undefined flavor jets counted: F_undefined_b_c_light
      for ( it_x = Fxijk . begin(); it_x != Fxijk . end(); it_x++ )
	{
	  map< int, map< int, map< int, double > > > map3 = it_x -> second;        
	  for ( it_i = map3 . begin(); it_i != map3 . end(); it_i++ )
	    {
	      map< int, map< int, double > > map2 = it_i -> second;        
	      for ( it_j = map2 . begin(); it_j != map2 . end(); it_j++ )
		{
		  map< int, double > map1 = it_j -> second;        
		  for ( it_k = map1 . begin(); it_k != map1 . end(); it_k++ )
		    {
		      int xx = it_x -> first;
		      int ii = it_i -> first;
		      int jj = it_j -> first;
		      int kk = it_k -> first;
		      _outputFileMC2 << "F_" << xx << "_" << ii << "_" << jj << "_" << kk << " = " << Fxijk[xx][ii][jj][kk] << endl;
		    }
		}
	    }
	}

    }
  

  _outputFile.close(); 

  _outputFileTable.close(); 
  _outputFileTable_2 . close();
  _outputFileTable_3 . close();


  //_outputFileMC.close(); 
  _outputFileMC2.close(); 

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TtTagConsistency);

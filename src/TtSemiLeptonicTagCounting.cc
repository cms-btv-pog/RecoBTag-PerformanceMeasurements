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
// $Id: TtSemiLeptonicTagCounting.cc,v 1.1.2.2 2008/03/21 16:59:37 kukartse Exp $
//
//

#include "RecoBTag/PerformanceMeasurements/interface/TtSemiLeptonicTagCounting.h"

TtSemiLeptonicTagCounting::TtSemiLeptonicTagCounting(const edm::ParameterSet& iConfig)
{
  
  dRCut        = iConfig.getParameter<double>("dRGenJetMatch");
  nCaloJetsLow = iConfig.getParameter<int>("nCaloJetsLow");
  nLeptonLow   = iConfig.getParameter<int>("nLeptonLow");

  dLow         = iConfig.getParameter<double>("discriminatorLow");
  dHigh        = iConfig.getParameter<double>("discriminatorHigh");
  dStep        = iConfig.getParameter<double>("discriminatorStep");

  dLow_2         = iConfig.getParameter<double>("discriminatorLow_2");
  dHigh_2        = iConfig.getParameter<double>("discriminatorHigh_2");
  dStep_2        = iConfig.getParameter<double>("discriminatorStep_2");
  
  dLow_3         = iConfig.getParameter<double>("discriminatorLow_3");
  dHigh_3        = iConfig.getParameter<double>("discriminatorHigh_3");
  dStep_3        = iConfig.getParameter<double>("discriminatorStep_3");
  
  _interFreq   = iConfig.getParameter<int>("intermediateResultFrequency");
  _dataType    = iConfig.getParameter<string>("dataType");
  _lookForGenJetMatch    = iConfig.getParameter<bool>("lookForGenJetMatch");

  _jetSource    = iConfig.getParameter<string>("jetSource");
  _electronSource    = iConfig.getParameter<string>("electronSource");
  _muonSource    = iConfig.getParameter<string>("muonSource");
  _METSource    = iConfig.getParameter<string>("METSource");
  _jetTagSource    = iConfig.getParameter<string>("jetTagSource");
  _jetTagSource_2  = iConfig.getParameter<string>("jetTagSource_2");
  _jetTagSource_3  = iConfig.getParameter<string>("jetTagSource_3");

  //_genJetSource    = iConfig.getParameter<string>("genJetSource");
  _outputFileName    = iConfig.getParameter<string>("outputFileName");

}


TtSemiLeptonicTagCounting::~TtSemiLeptonicTagCounting()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TtSemiLeptonicTagCounting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  /* 
  //Handle< vector< GenJet > > genJets; //obsolete since CSA07
  Handle< vector<genParticle > > genParticleCollection;
  vector<int> goodGenJetIndex; // indices of the jets that pass quality cuts
  if ( _dataType == "MC" )
    {
      
      
      // loop over gen jets
      //iEvent . getByLabel( _genJetSource, genJets );
      iEvent . getByLabel( _genJetSource, genParticleCollection );
      vector<reco::GenJet>::const_iterator gjet;
      RooGKCounter nOfGoodGenJets;
      for ( gjet = genJets -> begin(); gjet != genJets -> end(); gjet++ )  // loop over calo jets
	{
	  if ( ( fabs( (*gjet) . eta() ) < 2.5 ) &&
	       ( (*gjet) . pt() > 20.0 ) )
	    {
	      nOfGoodGenJets . count();
	      
	      goodGenJetIndex . push_back( gjet - genJets -> begin() ); // get the iterator index
	    }
	  
	}
    }
  */

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

  for ( cjet = caloJets -> begin(); cjet != caloJets -> end(); cjet++ )  // loop over calo jets
    {
      
      bool _foundGenJetMatch = true;

      /*
      if ( _dataType == "MC" && _lookForGenJetMatch )
	{
	  _foundGenJetMatch = ( findJetMatch( (*cjet), *genJets, goodGenJetIndex, dRCut ) > -1 );
	}
      */

      if ( 
	  ( fabs( (*cjet) . eta() ) < 2.4 ) &&
	  ( (*cjet) . pt() > 25.0 ) &&
	  _foundGenJetMatch )
	{
	  
	  nOfGoodJets . count();

	  map<double,bool> isTagged;
	  map<double,bool> isTagged_2;
	  map<double,bool> isTagged_3;
	  for(double _d = dLow; _d <= dHigh; _d += 0.5) isTagged[_d] = false;
	  for(double _d = dLow_2; _d <= dHigh_2; _d += 0.5) isTagged_2[_d] = false;
	  for(double _d = dLow_3; _d <= dHigh_3; _d += 0.5) isTagged_3[_d] = false;

	  JetTagRef _theJetTagRef = cjet -> getBJetTagRef( _jetTagSource );
	  JetTagRef _theJetTagRef_2 = cjet -> getBJetTagRef( _jetTagSource_2 );
	  JetTagRef _theJetTagRef_3 = cjet -> getBJetTagRef( _jetTagSource_3 );

	  // check if the jet is tagged
	  if ( _theJetTagRef . isNonnull() )
	    {
	      const JetTag * theJetTag = &*( cjet -> getBJetTagRef( _jetTagSource ) );

	      for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
		{
		  if ( theJetTag -> discriminator() > _d )
		    {
		      nTaggedJets[_d] . count(); // check the discriminator
		      
		      isTagged[_d] = true;
		    }
		}
	    }

	  // check if the jet is tagged (second tagger)
	  if ( _theJetTagRef_2 . isNonnull() )
	    {
	      const JetTag * theJetTag_2 = &*( cjet -> getBJetTagRef( _jetTagSource_2 ) );

	      //DEBUG different taggers' info
	      //cout << "DEBUG: b-tag discriminator 2 = " << theJetTag_2->discriminator() << endl << endl;

	      //
	      for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
		{
		  if ( theJetTag_2 -> discriminator() > _d )
		    {
		      nTaggedJets_2[_d] . count(); // check the discriminator
		      
		      isTagged_2[_d] = true;
		    }
		}
	    }

	  // check if the jet is tagged (third tagger)
	  if ( _theJetTagRef_3 . isNonnull() )
	    {
	      const JetTag * theJetTag_3 = &*( cjet -> getBJetTagRef( _jetTagSource_3 ) );

	      //DEBUG different taggers' info
	      //cout << "DEBUG: b-tag discriminator 2 = " << theJetTag_3->discriminator() << endl << endl;

	      //
	      for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
		{
		  if ( theJetTag_3 -> discriminator() > _d )
		    {
		      nTaggedJets_3[_d] . count(); // check the discriminator
		      
		      isTagged_3[_d] = true;
		    }
		}
	    }


	  // ---- jet flavor (MC only!!!)
	  if ( _dataType == "MC" ) 
	    {
	      int jetFlavor = cjet -> getPartonFlavour();
	      
	      goodCaloJetFlavorID . push_back( jetFlavor );
	      
	      if ( jetFlavor == 5 )
		{
		  //nOfGoodJets . count();
		  nOfGoodBJets . count();
		  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
		    {
		      if ( isTagged[_d] ) nOfTaggedBJets[_d] . count();
		    }
		  for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
		    {
		      if ( isTagged_2[_d] ) nOfTaggedBJets_2[_d] . count();
		    }
		  for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
		    {
		      if ( isTagged_3[_d] ) nOfTaggedBJets_3[_d] . count();
		    }
		}
	      else if ( jetFlavor == 4 )
		{
		  //nOfGoodJets . count();
		  nOfGoodCJets . count();
		  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
		    {
		      if ( isTagged[_d] ) nOfTaggedCJets[_d] . count();
		    }
		  for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
		    {
		      if ( isTagged_2[_d] ) nOfTaggedCJets_2[_d] . count();
		    }
		  for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
		    {
		      if ( isTagged_3[_d] ) nOfTaggedCJets_3[_d] . count();
		    }
		  
		}
	      else if ( jetFlavor == 1 || jetFlavor == 2 || jetFlavor == 3 || jetFlavor == 21 )
		{
		  //nOfGoodJets . count();
		  nOfGoodLJets . count();
		  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
		    {
		      if ( isTagged[_d] ) nOfTaggedLJets[_d] . count();
		    }
		  for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
		    {
		      if ( isTagged_2[_d] ) nOfTaggedLJets_2[_d] . count();
		    }
		  for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
		    {
		      if ( isTagged_3[_d] ) nOfTaggedLJets_3[_d] . count();
		    }
		}
	      else if ( jetFlavor == 0 )
		{
		  //nOfGoodJets . count();
		  nOfGoodUnknownJets . count();
		  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
		    {
		      if ( isTagged[_d] ) nOfTaggedUnknownJets[_d] . count();
		    }
		  for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
		    {
		      if ( isTagged_2[_d] ) nOfTaggedUnknownJets_2[_d] . count();
		    }
		  for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
		    {
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
  for ( el = electrons -> begin(); el != electrons -> end(); el++)
    {
      if ( ( (*el) . pt() > 20.0 ) &&
           ( fabs( (*el) . eta() ) < 2.4 ) )
        {
	  
	  nOfGoodElectrons . count();

	}
    }
  
  
  // loop over muons
  vector<TopMuon>::const_iterator mu;
  RooGKCounter nOfGoodMuons( "" );
  for ( mu = muons -> begin(); mu != muons -> end(); mu++)
    {
      if ( ( (*mu) . pt() > 20.0 ) &&
           ( fabs( (*mu) . eta() ) < 2.4 ) )
        {
	  
	  nOfGoodMuons . count();
	}
    }

  // loop over METs
  vector<TopMET>::const_iterator met;
  RooGKCounter nOfGoodMETs( "" );
  for ( met = METs -> begin(); met != METs -> end(); met++)
    {
      if ( ( (*met) . et() > 25.0 ) )
        {
	  
          nOfGoodMETs . count();
	  
        }
    }
  

  // ====== SELECTION =================
  //
  //
  if ( 
      nOfGoodJets . getCount() >= nCaloJetsLow &&
      ( nOfGoodElectrons . getCount() + nOfGoodMuons . getCount() >= nLeptonLow ) &&
      nOfGoodMETs . getCount() > 0 )
    {
      
      selectedEvents . count();

      nOfPassedBJets . increment( nOfGoodBJets . getCount() );
      nOfPassedCJets . increment( nOfGoodCJets . getCount() );
      nOfPassedLJets . increment( nOfGoodLJets . getCount() );
      nOfPassedUnknownJets . increment( nOfGoodUnknownJets . getCount() );

      Fijk[ (int)nOfGoodBJets.getCount() ][ (int)nOfGoodCJets.getCount() ][ (int)nOfGoodLJets.getCount() ]++;
      Fxijk[ (int)nOfGoodUnknownJets.getCount() ][ (int)nOfGoodBJets.getCount() ][ (int)nOfGoodCJets.getCount() ][ (int)nOfGoodLJets.getCount() ]++;

      for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
	{
	  nOfTaggedPassedBJets[_d] . increment( nOfTaggedBJets[_d] . getCount() );
	  nOfTaggedPassedCJets[_d] . increment( nOfTaggedCJets[_d] . getCount() );
	  nOfTaggedPassedLJets[_d] . increment( nOfTaggedLJets[_d] . getCount() );
	  nOfTaggedPassedUnknownJets[_d] . increment( nOfTaggedUnknownJets[_d] . getCount() );

	  Nn[_d][ nTaggedJets[_d] . getCount() ]++;
	}
      for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
	{
	  nOfTaggedPassedBJets_2[_d] . increment( nOfTaggedBJets_2[_d] . getCount() );
	  nOfTaggedPassedCJets_2[_d] . increment( nOfTaggedCJets_2[_d] . getCount() );
	  nOfTaggedPassedLJets_2[_d] . increment( nOfTaggedLJets_2[_d] . getCount() );
	  nOfTaggedPassedUnknownJets_2[_d] . increment( nOfTaggedUnknownJets_2[_d] . getCount() );
	  
	  Nn_2[_d][ nTaggedJets_2[_d] . getCount() ]++;
	}
      for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
	{
	  nOfTaggedPassedBJets_3[_d] . increment( nOfTaggedBJets_3[_d] . getCount() );
	  nOfTaggedPassedCJets_3[_d] . increment( nOfTaggedCJets_3[_d] . getCount() );
	  nOfTaggedPassedLJets_3[_d] . increment( nOfTaggedLJets_3[_d] . getCount() );
	  nOfTaggedPassedUnknownJets_3[_d] . increment( nOfTaggedUnknownJets_3[_d] . getCount() );
	  
	  Nn_3[_d][ nTaggedJets_3[_d] . getCount() ]++;
	}
    }  

  // This is probably obsolete. It used to be utilized to get some intermediate results
  // without interrupting the job.
  /* commented out on October 16, 2007
  if ( fmod( (double)(eventCounter . getCount()), (double)_interFreq ) < (1.0/(double)_interFreq))
    {
      interCounter . count();

      _outputFile << " ---> Intermediate: after " << interCounter . getCount() * _interFreq << " events" << endl;
      _outputFile << "Intermediate: Processed events:  " << eventCounter . getCount() << endl;
      _outputFile << "Intermediate: Selected events:  " << selectedEvents . getCount() << endl;
      _outputFile << " ---> Intermediate:Number of events with n tagged jets: " << endl;
      map< int, int >::const_iterator it_n;
      for( it_n = Nn[dLow] . begin(); it_n != Nn[dLow] . end(); it_n++ )
	{
	  _outputFile << " Intermediate: N_" << it_n -> first << " = " << it_n -> second << endl;
	}
      
      if ( _dataType == "MC" )
	{
	  _outputFile << " ---> Intermediate: Number of good light jets: " << nOfPassedLJets . getCount() << endl;
	  _outputFile << " ---> Intermediate: Number of tagged light jets: " << nOfTaggedPassedLJets[dLow] . getCount() << endl;
	  _outputFile << " ---> Intermediate: Number of good b jets: " << nOfPassedBJets . getCount() << endl;
	  _outputFile << " ---> Intermediate: Number of tagged b jets: " << nOfTaggedPassedBJets[dLow] . getCount() << endl;
	  _outputFile << " ---> Intermediate: Number of good c jets: " << nOfPassedCJets . getCount() << endl;
	  _outputFile << " ---> Intermediate: Number of tagged c jets: " << nOfTaggedPassedCJets[dLow] . getCount() << endl << endl;
	}
    }  
  */
}


// ------------ method called once each job just before starting event loop  ------------
void 
TtSemiLeptonicTagCounting::beginJob(const edm::EventSetup&)
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

  _outputFileMC . open( _mcFileName . c_str() );
  _outputFileMC2 . open( _mcFileName2 . c_str() );

  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
    {
      for (int i=0; i<=4; i++) Nn[_d][i] = 0;
    }

}

int TtSemiLeptonicTagCounting::findJetMatch( const reco::Jet & theJet, const std::vector<reco::GenJet> & jets, std::vector<int> & jetIndex, double dRCut ) {

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
TtSemiLeptonicTagCounting::endJob() {

  _outputFile << "============= TtSemiLeptonicTagCounting summary ================" << endl;
  _outputFile << "Processed events:  " << eventCounter . getCount() << endl;
  _outputFile << "Selected events:  " << selectedEvents . getCount() << endl;

  map< int, map< int, map< int, map< int, int > > > >::const_iterator it_x;  
  map< int, map< int, map< int, int > > >::const_iterator it_i;  
  map< int, map< int, int > >::const_iterator it_j;  
  map< int, int >::const_iterator it_k;  
  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
    {

      _outputFile << "###=== Discriminator > " << _d << " ===###---------------------------------" << endl;
      
      _outputFile << " ---> Number of events with n tagged jets: " << endl;
      map< int, int >::const_iterator it_n;
      for( it_n = Nn[_d] . begin(); it_n != Nn[_d] . end(); it_n++ )
	{
	  _outputFile << " N_" << it_n -> first << " = " << it_n -> second << endl;
	}
      
      if ( _dataType == "MC" )
	{
	  _outputFile << " ---> Number of good light jets: " << nOfPassedLJets . getCount() << endl;
	  _outputFile << " ---> Number of tagged light jets: " << nOfTaggedPassedLJets[_d] . getCount() << endl;
	  _outputFile << " ---> Number of good b jets: " << nOfPassedBJets . getCount() << endl;
	  _outputFile << " ---> Number of tagged b jets: " << nOfTaggedPassedBJets[_d] . getCount() << endl;
	  _outputFile << " ---> Number of good c jets: " << nOfPassedCJets . getCount() << endl;
	  _outputFile << " ---> Number of tagged c jets: " << nOfTaggedPassedCJets[_d] . getCount() << endl;
	  _outputFile << " ---> Number of good unknown jets: " << nOfPassedUnknownJets . getCount() << endl;
	  _outputFile << " ---> Number of tagged unknown jets: " << nOfTaggedPassedUnknownJets[_d] . getCount() << endl;
	}
    }

  //===> table
  char buf[10];
  _outputFileTable << "     discr       N_0       N_1       N_2       N_3       N_4       N_b N_btagged       N_c N_ctagged       N_l N_ltagged       N_x N_xtagged" << endl;
  for(double _d = dLow; _d < dHigh+0.0001; _d += dStep) // loop over discriminator values
    {
      sprintf( buf, "%10.2f", _d);
      _outputFileTable << buf;
      sprintf( buf, "%10d", Nn[_d][0]);
      _outputFileTable << buf;
      sprintf( buf, "%10d", Nn[_d][1]);
      _outputFileTable << buf;
      sprintf( buf, "%10d", Nn[_d][2]);
      _outputFileTable << buf;
      sprintf( buf, "%10d", Nn[_d][3]);
      _outputFileTable << buf;
      sprintf( buf, "%10d", Nn[_d][4]);
      _outputFileTable << buf;

      if ( _dataType == "MC" )
	{
	  sprintf( buf, "%10d", nOfPassedBJets . getCount() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedBJets[_d] . getCount() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10d", nOfPassedCJets . getCount() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedCJets[_d] . getCount() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10d", nOfPassedLJets . getCount() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedLJets[_d] . getCount() );
	  _outputFileTable << buf;

	  sprintf( buf, "%10d", nOfPassedUnknownJets . getCount() );
	  _outputFileTable << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedUnknownJets[_d] . getCount() );
	  _outputFileTable << buf << endl;

	}
      else
	{
	  for ( int i = 0; i < 8; i++ )
	    {
	      sprintf( buf, "%10d", 0 );
	      _outputFileTable << buf;
	    }
	  _outputFileTable << endl;
	}
    }

 // table for the second tagger
  _outputFileTable_2 << "     discr       N_0       N_1       N_2       N_3       N_4       N_b N_btagged       N_c N_ctagged       N_l N_ltagged       N_x N_xtagged" << endl;
  for(double _d = dLow_2; _d < dHigh_2+0.0001; _d += dStep_2) // loop over discriminator values
    {
      sprintf( buf, "%10.2f", _d);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10d", Nn_2[_d][0]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10d", Nn_2[_d][1]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10d", Nn_2[_d][2]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10d", Nn_2[_d][3]);
      _outputFileTable_2 << buf;
      sprintf( buf, "%10d", Nn_2[_d][4]);
      _outputFileTable_2 << buf;

      if ( _dataType == "MC" )
	{
	  sprintf( buf, "%10d", nOfPassedBJets . getCount() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedBJets_2[_d] . getCount() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10d", nOfPassedCJets . getCount() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedCJets_2[_d] . getCount() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10d", nOfPassedLJets . getCount() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedLJets_2[_d] . getCount() );
	  _outputFileTable_2 << buf;

	  sprintf( buf, "%10d", nOfPassedUnknownJets . getCount() );
	  _outputFileTable_2 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedUnknownJets_2[_d] . getCount() );
	  _outputFileTable_2 << buf << endl;

	}
      else
	{
	  for ( int i = 0; i < 8; i++ )
	    {
	      sprintf( buf, "%10d", 0 );
	      _outputFileTable_2 << buf;
	    }
	  _outputFileTable_2 << endl;
	}
    }

  // table for the third tagger
  _outputFileTable_3 << "     discr       N_0       N_1       N_2       N_3       N_4       N_b N_btagged       N_c N_ctagged       N_l N_ltagged       N_x N_xtagged" << endl;
  for(double _d = dLow_3; _d < dHigh_3+0.0001; _d += dStep_3) // loop over discriminator values
    {
      sprintf( buf, "%10.2f", _d);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10d", Nn_3[_d][0]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10d", Nn_3[_d][1]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10d", Nn_3[_d][2]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10d", Nn_3[_d][3]);
      _outputFileTable_3 << buf;
      sprintf( buf, "%10d", Nn_3[_d][4]);
      _outputFileTable_3 << buf;

      if ( _dataType == "MC" )
	{
	  sprintf( buf, "%10d", nOfPassedBJets . getCount() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedBJets_3[_d] . getCount() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10d", nOfPassedCJets . getCount() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedCJets_3[_d] . getCount() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10d", nOfPassedLJets . getCount() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedLJets_3[_d] . getCount() );
	  _outputFileTable_3 << buf;

	  sprintf( buf, "%10d", nOfPassedUnknownJets . getCount() );
	  _outputFileTable_3 << buf;
	  sprintf( buf, "%10d", nOfTaggedPassedUnknownJets_3[_d] . getCount() );
	  _outputFileTable_3 << buf << endl;

	}
      else
	{
	  for ( int i = 0; i < 8; i++ )
	    {
	      sprintf( buf, "%10d", 0 );
	      _outputFileTable_3 << buf;
	    }
	  _outputFileTable_3 << endl;
	}
    }

  
  if ( _dataType == "MC" )
    {
      for ( it_i = Fijk . begin(); it_i != Fijk . end(); it_i++ )
	{
	  map< int, map< int, int > > map2 = it_i -> second;        
	  for ( it_j = map2 . begin(); it_j != map2 . end(); it_j++ )
	    {
	      map< int, int > map1 = it_j -> second;        
	      for ( it_k = map1 . begin(); it_k != map1 . end(); it_k++ )
		{
		  int ii = it_i -> first;
		  int jj = it_j -> first;
		  int kk = it_k -> first;
		  _outputFileMC << "F_" << ii << "_" << jj << "_" << kk << " = " << Fijk[ii][jj][kk] << endl;
		}
	    }
	}

      // dataset flavor structure with undefined flavor jets counted: F_undefined_b_c_light
      for ( it_x = Fxijk . begin(); it_x != Fxijk . end(); it_x++ )
	{
	  map< int, map< int, map< int, int > > > map3 = it_x -> second;        
	  for ( it_i = map3 . begin(); it_i != map3 . end(); it_i++ )
	    {
	      map< int, map< int, int > > map2 = it_i -> second;        
	      for ( it_j = map2 . begin(); it_j != map2 . end(); it_j++ )
		{
		  map< int, int > map1 = it_j -> second;        
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


  _outputFileMC.close(); 
  _outputFileMC2.close(); 

}

//define this as a plug-in
//DEFINE_FWK_MODULE(TtSemiLeptonicTagCounting);

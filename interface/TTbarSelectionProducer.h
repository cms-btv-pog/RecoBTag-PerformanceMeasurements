


// -*- C++ -*-
//
// Package:    TTbarSelectionProducer
// Class:      TTbarSelectionProducer
// 
/**\class TTbarSelectionProducer TTbarSelectionProducer.cc bTag/TTbarSelectionProducer/src/TTbarSelectionProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Jeremy,B25/117,6262,
//         Created:  Mon Nov 26 12:32:34 CET 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//--------------------PAT includes
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "TLorentzVector.h"


#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//----- histo service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"


//
// class declaration
//

class TTbarSelectionProducer : public edm::EDProducer {
   public:
      explicit TTbarSelectionProducer(const edm::ParameterSet&);
      ~TTbarSelectionProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      
      bool isData_;
      
      //Configuration for electrons
      
      edm::InputTag electronColl_;
      
      double electron_cut_pt_ ; 
      double electron_cut_eta_;
      double electron_cut_iso_;
   
      //Configuration for muons
   
   
      edm::InputTag muonColl_ ;
      double muon_cut_pt_ ;
      double muon_cut_eta_ ;
      double muon_cut_iso_;
   
      //Configuration for jets 
   
      edm::InputTag jetColl_;
      double jet_cut_pt_ ;
      double jet_cut_eta_ ;
   
      //Configuration for met 
      edm::InputTag metColl_;
      double met_cut_ ;
   
      // Extract info for BeamSpot
      //bool doBeamSpot_  ;  
      //edm::InputTag beamSpotProducer_;
      
      
      //Configuration for tracks 
      
      edm::InputTag trackColl_;
      // ----------member data ---------------------------
      
      void GetLeptonPair(std::vector<TLorentzVector> electrons, std::vector<TLorentzVector> muons, 
                          std::vector<int> electronCharges, std::vector<int> muonCharges, 
			  int &idxLept1, int &idxLept2, int &thechannel);

      // ----- histo -------
      edm::Service<TFileService> fs;
      TH1F* hcheck_cutflow        ;
      TH1F* hcheck_m_ee           ;
      TH1F* hcheck_m_emu          ;
      TH1F* hcheck_m_mumu         ;
      TH1F* hcheck_met_ee         ;
      TH1F* hcheck_met_emu        ;
      TH1F* hcheck_met_mumu       ;
			  
};

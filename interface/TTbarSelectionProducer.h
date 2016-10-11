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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//--------------------PAT includes
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

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
      virtual void beginRun(const edm::Run & iRun, edm::EventSetup const & iSetup);
      virtual void endRun(const edm::Run & iRun, edm::EventSetup const & iSetup);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      //triggers
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      std::string triggerBitsProc_;

      HLTConfigProvider hltConfig;
      std::vector<int> trigChannels_;
      std::vector<std::string> trigNamesToSel_;
      bool doTrigSel_;
  
      //gen particles
      edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticleCollectionName_;
      edm::EDGetTokenT<GenEventInfoProduct> generatorevt_;     
      edm::EDGetTokenT<LHEEventProduct> generatorlhe_;

      //MET filters
      std::vector<edm::EDGetTokenT<edm::TriggerResults> >metFilters_;
      edm::EDGetTokenT<bool> RecoHBHENoiseFilter_;
      std::vector<std::string> metFiltersToApply_;

      //vertices, beam spot
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::BeamSpot> bsToken_;

      //Configuration for electrons      
      edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > electronIdMapToken_;
      double electron_cut_pt_ ; 
      double electron_cut_eta_;
      double electron_cut_iso_;
   
      //Configuration for muons
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      double muon_cut_pt_ ;
      double muon_cut_eta_ ;
      double muon_cut_iso_;
   
      //Configuration for jets 
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      double jet_cut_pt_ ;
      double jet_cut_eta_ ;
   
      //Configuration for met 
      edm::EDGetTokenT<pat::METCollection> metToken_;
      double met_cut_ ;
   
      // ----- histo -------
      std::map<std::string, TH1F *> histos_;
  
      //verbose level
      int verbose_;

      //return the channel selected
      int AssignChannel(std::vector<pat::Electron> &selElectrons,
			std::vector<pat::Muon> &selMuons,
			int trigWord);
};

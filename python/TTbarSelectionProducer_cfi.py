import FWCore.ParameterSet.Config as cms

ttbarselectionproducer = cms.EDProducer("TTbarSelectionProducer",
   isData           = cms.bool(True),
   electronColl     = cms.InputTag("selectedPatElectronsPF2PAT"),
   electron_cut_pt  = cms.double(20),
   electron_cut_eta = cms.double(2.5),
   electron_cut_iso = cms.double(0.15),
   muonColl         = cms.InputTag("selectedPatMuonsPF2PAT"),
   muon_cut_pt      = cms.double(20),
   muon_cut_eta     = cms.double(2.4),
   muon_cut_iso     = cms.double(0.15),
   jetColl          = cms.InputTag("selectedPatJetsPF2PAT"),
   jet_cut_pt       = cms.double(30),
   jet_cut_eta      = cms.double(2.5),
   metColl          = cms.InputTag("patMETsPF2PAT"),
   met_cut          = cms.double(40),
   trackColl        = cms.InputTag("generalTracks"),
)

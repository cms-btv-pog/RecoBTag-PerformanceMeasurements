import FWCore.ParameterSet.Config as cms

ttbarselectionfilter = cms.EDFilter('TTbarSelectionFilter',
   channel     = cms.InputTag("ttbarselectionproducer"),
   select_ee   = cms.bool(True),
   select_mumu = cms.bool(True),
   select_emu  = cms.bool(True)
)

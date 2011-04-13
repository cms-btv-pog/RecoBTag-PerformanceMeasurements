import FWCore.ParameterSet.Config as cms

HLTFilter = cms.EDFilter(
    'HLTFilter',

    tag = cms.string("TriggerResults::REDIGI"),
    hlt = cms.string("HLT_Mu9")
)

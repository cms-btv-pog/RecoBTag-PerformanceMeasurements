import FWCore.ParameterSet.Config as cms

PMConversionFilter = cms.EDFilter("PMConversionFilter",
  Electrons = cms.InputTag("patElectrons"),
  Conversions = cms.InputTag("trackerOnlyConversions"),
  Filter = cms.bool(False)
)

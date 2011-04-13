import FWCore.ParameterSet.Config as cms

PVFilter = cms.EDFilter(
    'PVFilter',

    pvTag = cms.string("offlinePrimaryVertices"),
    isDataInput = cms.bool(False)
)

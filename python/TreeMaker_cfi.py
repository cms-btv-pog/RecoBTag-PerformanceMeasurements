import FWCore.ParameterSet.Config as cms

TreeMaker = cms.EDAnalyzer(
    'TreeMaker',

    primaryVertices = cms.string("offlinePrimaryVertices"),
    muons = cms.string("selectedPatMuonsForPtRel"),
    jets = cms.string("selectedPatJetsAK5PF")
)

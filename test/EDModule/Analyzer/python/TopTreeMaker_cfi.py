import FWCore.ParameterSet.Config as cms

TopTreeMaker = cms.EDAnalyzer(
    'TopTreeMaker',

    mets = cms.string("patMETs"),
    muons = cms.string("selectedPatMuons"),
    jets = cms.string("selectedPatJets"),
    electrons = cms.string("selectedPatElectrons"),
    beamSpots = cms.string("offlineBeamSpot")
)

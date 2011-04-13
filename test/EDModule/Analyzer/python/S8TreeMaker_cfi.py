import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

S8TreeMaker = cms.EDAnalyzer(
    'S8TreeMaker',

    primaryVertices = cms.string("offlinePrimaryVertices"),
    muons = cms.string("selectedPatMuonsForPtRel"),
    electrons = cms.string("selectedPatElectronsForS8"),
    jets = cms.string("selectedPatJetsAK5PF"),
    triggers = cms.string("TriggerResults::REDIGI"),

    jetSelector = pfJetIDSelector.clone(),

    isPythia = cms.bool(False),
    saveTriggers = cms.bool(False)
)

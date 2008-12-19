import FWCore.ParameterSet.Config as cms

jetcuts = cms.PSet(
    MaxEta = cms.double(2.5),
    MinDeltaR = cms.double(0.4),
    MinPt = cms.double(20.0),
    MinPtRel = cms.double(-1.0) ## no pTrel cut

)



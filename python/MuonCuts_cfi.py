import FWCore.ParameterSet.Config as cms

muoncuts = cms.PSet(
    MinNHits = cms.int32(7),
    MinMuonPt = cms.double(6.0),
    MaxMuonChi2 = cms.double(5.0),
    MaxMuonEta = cms.double(2.5)
)



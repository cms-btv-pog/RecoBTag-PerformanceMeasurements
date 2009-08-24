import FWCore.ParameterSet.Config as cms

Taggability = cms.EDFilter("Taggability",
                           JetCollection = cms.InputTag('iterativeCone5CaloJets'),
                           ApplyJetCorrections = cms.bool(True),
                           jetCorrectionsLabel = cms.string('L2L3JetCorrectorIC5Calo'),
                           MinPt = cms.double(30.),
                           MaxEta = cms.double(2.4),
                           bTagTrackEventIPtagInfos = cms.string('impactParameterTagInfos'),
                           MinNtracks = cms.int32(1),
                           MinTrkPt = cms.double(1),
                           MinNjets = cms.int32(1),
                           PrimaryVertexCollection = cms.InputTag('offlinePrimaryVertices'),
                           MinNPrimaryVertices = cms.int32(0),
			   WriteHistograms = cms.bool(True)
                           )



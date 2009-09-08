import FWCore.ParameterSet.Config as cms

from RecoBTag.PerformanceMeasurements.OperatingPoints import *

OperatingPoints = cms.EDAnalyzer("OperatingPoints",
				 EstimateByMistagRate,
				 OperatingPointsbyMistagRate = cms.bool (True),
				 ApplyJetCorrections = cms.bool (True),
				 jetCorrectionsLabel = cms.string("L2L3JetCorrectorIC5Calo"),
				 jetcuts = cms.PSet(
	MaxEta = cms.double(2.5),
	MinPt = cms.double(30.0),
	MinNtracks = cms.int32(1),
	MinTrkPt = cms.double(1.),
	MinNjets = cms.int32(1)
	),
				 flavourSource = cms.InputTag('IC5byValAlgo'),
				 Jets = cms.string('iterativeCone5CaloJets'),
				 bTagTrackEventIPtagInfos = cms.string('impactParameterTagInfos'),
				 outputFile = cms.untracked.string('OperatingPoints.root'),
				 debug = cms.untracked.bool (False)
				 
				 )


                    

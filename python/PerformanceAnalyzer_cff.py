
import FWCore.ParameterSet.Config as cms

from RecoBTag.PerformanceMeasurements.OperatingPoints import *

Performance = cms.EDAnalyzer("PerformanceAnalyzer",

                             OperatingPoints31X,
                             #
                             # use jet corrections
                             #
                             useJetCorrections = cms.bool (True),
                             jetCorrectionsLabel = cms.string("L2L3JetCorrectorIC5Calo"),                     
    muoncuts = cms.PSet(
        MinNHits = cms.int32(7),
        MinMuonPt = cms.double(6.0),
        MaxMuonChi2 = cms.double(5.0),
        MaxMuonEta = cms.double(2.5)
    ),
    jetIdParameters = cms.PSet(
        vetoFlavour = cms.vstring(),
        rejectBCSplitting = cms.bool(False),
        physicsDefinition = cms.bool(False),
        coneSizeToAssociate = cms.double(0.3),
        fillLeptons = cms.bool(False),
        fillHeavyHadrons = cms.bool(False),
        fillPartons = cms.bool(True),
        mcSource = cms.string('source')
    ),
    Muons = cms.string('muons'),
    StoreTrackProba = cms.bool(False),
    #PSet jetIdParameters2 = {
    #       string mcSource = "source"
    #       bool fillPartons = true
    #       bool fillHeavyHadrons = false
    #       bool fillLeptons =  false
    #       double coneSizeToAssociate = 0.3
    #       bool physicsDefinition = true
    #       bool rejectBCSplitting = true
    #       vstring vetoFlavour = {  }
    #}
    #
    # prepare pset for TrackHistory
    #
    jetcuts = cms.PSet(
        MaxEta = cms.double(2.5),
        MinDeltaR = cms.double(0.4),
        MinPt = cms.double(20.0),
        MinPtRel = cms.double(-1.0)
    ),
    bTagTrackEventIPtagInfos = cms.string(''),
    # bTagTrackEventIPtagInfos = cms.string('impactParameterTagInfos'),
    flavourMatchOption = cms.string('genParticles'),
    WeightHistograms = cms.bool(False),
    TrackCollection = cms.untracked.string('generalTracks'),
    GenJets = cms.string('iterativeCone5GenJets'),
    Jets = cms.string('iterativeCone5CaloJets'),
    bTagTrackEvent = cms.bool(False),
    #InputTag simG4 = g4SimHits
    outputFile = cms.untracked.string('results.root'),
    StorePtHat = cms.bool(False),
    SimTracks = cms.string('g4SimHits'),
    AwayJetTagger = cms.string('TCL'),
    flavourSource = cms.InputTag("IC5byValAlgo"),
    StoreWeightsInNtuple = cms.bool(False),
    PrimaryVertexCollection = cms.untracked.string('offlinePrimaryVerticesFromCTFTracks'),
    WritePerformancePlots = cms.bool(True),
    debug = cms.untracked.bool(False)
)




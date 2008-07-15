# The following comments couldn't be translated into the new config version:

#TCL
#TCM
#TCT
#JBPL
#JBPM
#JBPT
#SLT
#SVM
#SVT
#CSVL 
#CSVM 
#CSVT
import FWCore.ParameterSet.Config as cms

Performance = cms.EDAnalyzer("PerformanceAnalyzer",
    # definition of the Operating Points (L,M,T)
    # cuts estimated either by thomas on 21X, or using old francisco's ones
    # sorted as TCL,TCM,TCT,JPL,JPM,JPT,JBPL, JBPM, JBPT, SLT, SVM, SVT, CSVL, CSVM, CSVT
    bTagCutList = cms.untracked.vdouble(2.0, 4.6, 4.7, 0.26, 0.5, 
        0.76, 1.2, 2.3, 3.2, 0.8, 
        2.0, 3.6, 0.0, 37.0, 0.84, 
        0.96),
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
    associationModule = cms.string('TrackAssociatorByHits'),
    Muons = cms.string('muons'),
    bestMatchByMaxValue = cms.bool(True),
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
    recoTrackModule = cms.string('generalTracks'),
    jetcuts = cms.PSet(
        MaxEta = cms.double(2.5),
        MinDeltaR = cms.double(0.4),
        MinPt = cms.double(20.0),
        MinPtRel = cms.double(-1.0)
    ),
    trackingParticleModule = cms.string('mergedtruth'),
    flavourMatchOption = cms.string('hepMC'),
    WeightHistograms = cms.bool(False),
    TrackCollection = cms.untracked.string('ctfWithMaterialTracks'),
    GenJets = cms.string('iterativeCone5GenJets'),
    Jets = cms.string('iterativeCone5CaloJets'),
    bTagTrackEvent = cms.bool(False),
    #InputTag simG4 = g4SimHits
    outputFile = cms.untracked.string('results.root'),
    StorePtHat = cms.bool(False),
    SimTracks = cms.string('g4SimHits'),
    AwayJetTagger = cms.string('TCL'),
    flavourSource = cms.InputTag("IC5byValAlgo"),
    trackingParticleProduct = cms.string('MergedTrackTruth'),
    StoreWeightsInNtuple = cms.bool(False),
    # b-tagedd jet collection
    bTaggerList = cms.untracked.vstring('trackCountingHighEffBJetTags', 
        'trackCountingHighPurBJetTags', 
        'simpleSecondaryVertexBJetTags', 
        'combinedSecondaryVertexBJetTags', 
        'jetProbabilityBJetTags', 
        'jetBProbabilityBJetTags'),
    PrimaryVertexCollection = cms.untracked.string('offlinePrimaryVerticesFromCTFTracks'),
    WritePerformancePlots = cms.bool(True)
)




import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *


mistag = cms.EDAnalyzer("MistagAnalyzer",
    trackClassifier, 
    useTrackHistory          = cms.bool(True),
    outputFile               = cms.untracked.string('results_50_80.root'),
    
    ntrackMin                = cms.int32(0),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    tagCut                   = cms.double(5.3),
    MinPt                    = cms.double(30.0),
    vetoPos                  = cms.double(4.0),
    isData                   = cms.bool(True),
    Jets                     = cms.string('iterativeCone5CaloJets'),
    
    flavourMatchOption       = cms.string('genParticle'),
    flavourSource            = cms.InputTag("IC5byValAlgo"),
    
    jetIdParameters         = cms.PSet(
        vetoFlavour         = cms.vstring(),
        rejectBCSplitting   = cms.bool(False),
        physicsDefinition   = cms.bool(False),
        coneSizeToAssociate = cms.double(0.3),
        fillLeptons         = cms.bool(False),
        fillHeavyHadrons    = cms.bool(False),
        fillPartons         = cms.bool(True),
        mcSource            = cms.string('source')
    ),
    
    trackCNegHPModuleName   = cms.string('negativeTrackCountingHighPur'),
    trackCNegHEModuleName   = cms.string('negativeTrackCountingHighEffJetTags'),
    trackCHPModuleName      = cms.string('trackCountingHighPurBJetTags'),
    
    jetPModuleName          = cms.string('jetProbabilityBJetTags'),
    jetPPosModuleName       = cms.string('positiveOnlyJetProbabilityJetTags'),
    jetPNegModuleName       = cms.string('negativeOnlyJetProbabilityJetTags'),
    trackCHEModuleName      = cms.string('trackCountingHighEffBJetTags'),
    
    combinedSvtxModuleName    = cms.string(''),
    combinedSvtxNegModuleName = cms.string(''),
    svtxModuleName            = cms.string('simpleSecondaryVertexBJetTags'),
    svtxNegModuleName         = cms.string('simpleSecondaryVertexNegativeBJetTags'),
    softMuonModuleName        = cms.string('positiveSoftMuonBJetTags'),
    softMuonNegModuleName     = cms.string('negativeSoftMuonBJetTags'),
    softMuonTagInfoName       = cms.string('softMuonTagInfos')
)

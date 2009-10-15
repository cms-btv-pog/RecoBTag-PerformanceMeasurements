import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *


mistag = cms.EDAnalyzer("MistagAnalyzer",
    trackClassifier, 
    useTrackHistory          = cms.bool(False),
    #jetCorrector             = cms.string('L2L3JetCorrectorIcone5'),
    #jetCorrector             = cms.string('L2L3JetCorrectorIC5Calo'),
    jetCorrector             = cms.string('L2L3JetCorrectorAK5Calo'),
    outputFile               = cms.untracked.string('Ntuple.root'),
    longLivedDecayLenght     = cms.untracked.double(1e-14),
    
    ntrackMin                = cms.int32(0),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    tagCut                   = cms.double(5.3),
    MinPt                    = cms.double(30.0),
    vetoPos                  = cms.double(4.0),
    isData                   = cms.bool(True),
    #Jets                     = cms.string('iterativeCone5CaloJets'),
    Jets                     = cms.string('antikt5CaloJets'),
    
    #vertexClusteringDistance = cms.untracked.double(0.0001),
    #trackProducer            = cms.untracked.InputTag("generalTracks"),
    #badD0Pull                = cms.untracked.double(3.0),
    #trackAssociator          = cms.untracked.string('TrackAssociatorByHits'),
    flavourMatchOption       = cms.string('genParticle'),
    #flavourSource            = cms.InputTag("IC5byValAlgo"),
    flavourSource            = cms.InputTag("AK5byValAlgo"),
    
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
    
    combinedSvtxModuleName    = cms.string('combinedSecondaryVertexBJetTags'),
    combinedSvtxNegModuleName = cms.string('combinedSecondaryVertexNegativeBJetTags'),
    svtxModuleName            = cms.string('simpleSecondaryVertexBJetTags'),
    svtxNegModuleName         = cms.string('simpleSecondaryVertexNegativeBJetTags'),
    #softMuonModuleName        = cms.string('positiveSoftMuonBJetTags'),
    #softMuonNegModuleName     = cms.string('negativeSoftMuonBJetTags'),
    softMuonModuleName        = cms.string('positiveSoftLeptonByPtBJetTags'),
    softMuonNegModuleName     = cms.string('negativeSoftLeptonByPtBJetTags'),
    softMuonTagInfoName       = cms.string('softMuonTagInfos')
)

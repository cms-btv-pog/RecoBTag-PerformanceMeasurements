import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *


#list of input parameters will be cleanned. Some of them are not used anymore


btagana = cms.EDAnalyzer("BTagAnalyzer",
    trackClassifier, 
    useTrackHistory          = cms.bool(False),
    produceJetProbaTree      = cms.bool(True),
    producePtRelTemplate     = cms.bool(True),
    jetCorrector             = cms.string('ak5PFJetsL2L3'), # obsolete, need to be removed
    longLivedDecayLenght     = cms.untracked.double(1e-14),
    primaryVertexColl        = cms.string('offlinePrimaryVertices'),
    ntrackMin                = cms.int32(0),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    tagCut                   = cms.double(5.3),
    MinPt                    = cms.double(20.0),
    vetoPos                  = cms.double(4.0),
    isData                   = cms.bool(True),
    use_selected_tracks      = cms.bool(True),
    Jets                     = cms.string('ak5PFJets'),
    genJetCollection         = cms.string('ak5GenJets'),
    muonCollectionName       = cms.InputTag('muons'),
    flavourMatchOption       = cms.string('genParticle'),
    flavourSource            = cms.InputTag("AK5byValAlgo"),
    triggerTable             = cms.InputTag("TriggerResults"),
    svComputer               = cms.InputTag( "combinedSecondaryVertex" ),
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
    
    
    #list of netagive tagger for mistag rate
    trackCNegHPModuleName   = cms.string('negativeTrackCountingHighPur'),
    trackCNegHEModuleName   = cms.string('negativeTrackCountingHighEffJetTags'),
    trackCHPModuleName      = cms.string('trackCountingHighPurBJetTags'),
    
    jetBModuleName          = cms.string('jetBProbabilityBJetTags'),
    jetBNegModuleName       = cms.string('negativeOnlyJetBProbabilityJetTags'),
    jetBPosModuleName       = cms.string('positiveOnlyJetBProbabilityJetTags'),
    jetPModuleName          = cms.string('jetProbabilityBJetTags'),
    jetPPosModuleName       = cms.string('positiveOnlyJetProbabilityJetTags'),
    jetPNegModuleName       = cms.string('negativeOnlyJetProbabilityJetTags'),
    trackCHEModuleName      = cms.string('trackCountingHighEffBJetTags'),
    
    combinedSvtxModuleName    = cms.string('combinedSecondaryVertexBJetTags'),
    combinedSvtxNegModuleName = cms.string('combinedSecondaryVertexNegativeBJetTags'),
    combinedSvtxPosModuleName = cms.string('combinedSecondaryVertexPositiveBJetTags'),
        
    combinedSvtxRetrainedModuleName    = cms.string('combinedSecondaryVertexRetrainedBJetTags'),
    combinedSvtxNegRetrainedModuleName = cms.string('combinedSecondaryVertexNegativeRetrainedBJetTags'),
    combinedSvtxPosRetrainedModuleName = cms.string('combinedSecondaryVertexPositiveRetrainedBJetTags'),
			
    simpleIVFModuleNameHighPur      = cms.string('simpleInclusiveSecondaryVertexHighPurBJetTags'),
    simpleIVFModuleNameHighEff      = cms.string('simpleInclusiveSecondaryVertexHighEffBJetTags'),
    doubleIVFModuleNameHighEff      = cms.string('doubleSecondaryVertexHighEffBJetTags'),      
    combinedIVFModuleName	    = cms.string('combinedInclusiveSecondaryVertexBJetTags'),  
    combinedIVFPosModuleName	    = cms.string('combinedInclusiveSecondaryVertexPositiveBJetTags'),

    svtxModuleNameHighPur     = cms.string('simpleSecondaryVertexHighPurBJetTags'),
    svtxNegModuleNameHighPur  = cms.string('simpleSecondaryVertexNegativeHighPurBJetTags'),
    svtxModuleNameHighEff     = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    svtxNegModuleNameHighEff  = cms.string('simpleSecondaryVertexNegativeHighEffBJetTags'),
    softMuonModuleName        = cms.string('positiveSoftLeptonByPtBJetTags'),
    softMuonNegModuleName     = cms.string('negativeSoftLeptonByPtBJetTags'),
    softMuonTagInfoName       = cms.string('softMuonTagInfos'),
)

import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *


#list of input parameters will be cleanned. Some of them are not used anymore


btagana = cms.EDAnalyzer("BTagAnalyzer",
    trackClassifier, 
    useTrackHistory          = cms.bool(False),
    produceJetProbaTree      = cms.bool(True),
    producePtRelTemplate     = cms.bool(True),
    longLivedDecayLenght     = cms.untracked.double(1e-14),
    primaryVertexColl        = cms.string('offlinePrimaryVertices'),
    ntrackMin                = cms.int32(0),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    tagCut                   = cms.double(5.3),
    MinPt                    = cms.double(20.0),
    vetoPos                  = cms.double(4.0),
    use_selected_tracks      = cms.bool(True),
    Jets                     = cms.string('ak5PFJets'),
    genJetCollection         = cms.string('ak5GenJets'),
    muonCollectionName       = cms.InputTag('muons'),
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
        
    combinedSvtxJPModuleName    = cms.string('combinedSecondaryVertexJPBJetTags'),
    combinedSvtxNegJPModuleName = cms.string('combinedSecondaryVertexNegativeJPBJetTags'),
    combinedSvtxPosJPModuleName = cms.string('combinedSecondaryVertexPositiveJPBJetTags'),
        
    combinedSvtxRetrainedModuleName    = cms.string('combinedSecondaryVertexRetrainedBJetTags'),
    combinedSvtxRetrainedNegModuleName = cms.string('combinedSecondaryVertexRetrainedNegativeBJetTags'),
    combinedSvtxRetrainedPosModuleName = cms.string('combinedSecondaryVertexRetrainedPositiveBJetTags'),
        
    combinedCSVJPModuleName    = cms.string('combinedCSVJPBJetTags'),
    combinedCSVJPNegModuleName = cms.string('negativeCombinedCSVJPBJetTags'),
    combinedCSVJPPosModuleName = cms.string('positiveCombinedCSVJPBJetTags'),
        
    combinedCSVSLModuleName    = cms.string('combinedCSVSLBJetTags'),
    combinedCSVSLNegModuleName = cms.string('negativeCombinedCSVSLBJetTags'),
    combinedCSVSLPosModuleName = cms.string('positiveCombinedCSVSLBJetTags'),
        
    combinedCSVJPSLModuleName    = cms.string('combinedCSVJPSLBJetTags'),
    combinedCSVJPSLNegModuleName = cms.string('negativeCombinedCSVJPSLBJetTags'),
    combinedCSVJPSLPosModuleName = cms.string('positiveCombinedCSVJPSLBJetTags'),
        			
    simpleIVFModuleNameHighPur      = cms.string('simpleInclusiveSecondaryVertexHighPurBJetTags'),
    simpleIVFModuleNameHighEff      = cms.string('simpleInclusiveSecondaryVertexHighEffBJetTags'),
    doubleIVFModuleNameHighEff      = cms.string('doubleSecondaryVertexHighEffBJetTags'),      
    combinedIVFModuleName	    = cms.string('combinedInclusiveSecondaryVertexBJetTags'),  
    combinedIVFPosModuleName	    = cms.string('combinedInclusiveSecondaryVertexPositiveBJetTags'),

    svtxModuleNameHighPur     = cms.string('simpleSecondaryVertexHighPurBJetTags'),
    svtxNegModuleNameHighPur  = cms.string('simpleSecondaryVertexNegativeHighPurBJetTags'),
    svtxModuleNameHighEff     = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    svtxNegModuleNameHighEff  = cms.string('simpleSecondaryVertexNegativeHighEffBJetTags'),
                         
    #softMuonModuleName        = cms.string('positiveSoftLeptonByPtBJetTags'),
    #softMuonNegModuleName     = cms.string('negativeSoftLeptonByPtBJetTags'),
    #softMuonTagInfoName       = cms.string('softMuonTagInfos'),
                         
    softPFMuonModuleName        = cms.string('softPFMuonRetrainedBJetsTags'),
    softPFMuonPosModuleName     = cms.string('positiveSoftPFMuonRetrainedBJetsTags'),
    softPFMuonNegModuleName     = cms.string('negativeSoftPFMuonRetrainedBJetsTags'),
                         
    softPFElectronModuleName    = cms.string('softPFElectronRetrainedBJetsTags'),
    softPFElectronPosModuleName = cms.string('positiveSoftPFElectronRetrainedBJetsTags'),
    softPFElectronNegModuleName = cms.string('negativeSoftPFElectronRetrainedBJetsTags'),
                                                 
    softPFMuonTagInfoName       = cms.InputTag('softPFMuonsTagInfos'),
    softPFElectronTagInfoName   = cms.InputTag('softPFElectronsTagInfos'),
                         
    use_ttbar_filter      = cms.bool(False),
    channel       = cms.InputTag("ttbarselectionproducer"),
)

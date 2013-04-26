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
    muonCollectionName       = cms.InputTag('muons'),
    triggerTable             = cms.InputTag("TriggerResults"),
    svComputer               = cms.InputTag( "combinedSecondaryVertex" ),

    # list of taggers
    trackCHEBJetTags      = cms.string('trackCountingHighEffBJetTags'),
    trackCNegHEBJetTags   = cms.string('negativeTrackCountingHighEffJetTags'),

    trackCHPBJetTags      = cms.string('trackCountingHighPurBJetTags'),
    trackCNegHPBJetTags   = cms.string('negativeTrackCountingHighPurJetTags'),

    jetBPBJetTags          = cms.string('jetBProbabilityBJetTags'),
    jetBPNegBJetTags       = cms.string('negativeOnlyJetBProbabilityJetTags'),
    jetBPPosBJetTags       = cms.string('positiveOnlyJetBProbabilityJetTags'),

    jetPBJetTags          = cms.string('jetProbabilityBJetTags'),
    jetPNegBJetTags       = cms.string('negativeOnlyJetProbabilityJetTags'),
    jetPPosBJetTags       = cms.string('positiveOnlyJetProbabilityJetTags'),
    
    combinedSVBJetTags    = cms.string('combinedSecondaryVertexBJetTags'),
    combinedSVNegBJetTags = cms.string('combinedSecondaryVertexNegativeBJetTags'),
    combinedSVPosBJetTags = cms.string('combinedSecondaryVertexPositiveBJetTags'),
        
    combinedSVRetrainedBJetTags    = cms.string('combinedSecondaryVertexRetrainedBJetTags'),
    combinedSVRetrainedNegBJetTags = cms.string('combinedSecondaryVertexRetrainedNegativeBJetTags'),
    combinedSVRetrainedPosBJetTags = cms.string('combinedSecondaryVertexRetrainedPositiveBJetTags'),
        
    combinedCSVJPBJetTags    = cms.string('combinedCSVJPBJetTags'),
    combinedCSVJPNegBJetTags = cms.string('negativeCombinedCSVJPBJetTags'),
    combinedCSVJPPosBJetTags = cms.string('positiveCombinedCSVJPBJetTags'),
        
    combinedCSVSLBJetTags    = cms.string('combinedCSVSLBJetTags'),
    combinedCSVSLNegBJetTags = cms.string('negativeCombinedCSVSLBJetTags'),
    combinedCSVSLPosBJetTags = cms.string('positiveCombinedCSVSLBJetTags'),
        
    combinedCSVJPSLBJetTags    = cms.string('combinedCSVJPSLBJetTags'),
    combinedCSVJPSLNegBJetTags = cms.string('negativeCombinedCSVJPSLBJetTags'),
    combinedCSVJPSLPosBJetTags = cms.string('positiveCombinedCSVJPSLBJetTags'),
        			
    simpleIVFSVHighPurBJetTags      = cms.string('simpleInclusiveSecondaryVertexHighPurBJetTags'),
    simpleIVFSVHighEffBJetTags      = cms.string('simpleInclusiveSecondaryVertexHighEffBJetTags'),
    doubleIVFSVHighEffBJetTags      = cms.string('doubleSecondaryVertexHighEffBJetTags'),
    combinedIVFSVBJetTags	    = cms.string('combinedInclusiveSecondaryVertexBJetTags'),
    combinedIVFSVPosBJetTags	    = cms.string('combinedInclusiveSecondaryVertexPositiveBJetTags'),

    simpleSVHighPurBJetTags     = cms.string('simpleSecondaryVertexHighPurBJetTags'),
    simpleSVNegHighPurBJetTags  = cms.string('simpleSecondaryVertexNegativeHighPurBJetTags'),
    simpleSVHighEffBJetTags     = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    simpleSVNegHighEffBJetTags  = cms.string('simpleSecondaryVertexNegativeHighEffBJetTags'),
                         
    #softMuonBJetTags        = cms.string('positiveSoftLeptonByPtBJetTags'),
    #softMuonNegBJetTags     = cms.string('negativeSoftLeptonByPtBJetTags'),
    #softMuonTagInfoName       = cms.string('softMuonTagInfos'),
                         
    softPFMuonBJetTags        = cms.string('softPFMuonRetrainedBJetTags'),
    softPFMuonNegBJetTags     = cms.string('negativeSoftPFMuonRetrainedBJetTags'),
    softPFMuonPosBJetTags     = cms.string('positiveSoftPFMuonRetrainedBJetTags'),
    
    softPFElectronBJetTags    = cms.string('softPFElectronRetrainedBJetTags'),
    softPFElectronNegBJetTags = cms.string('negativeSoftPFElectronRetrainedBJetTags'),
    softPFElectronPosBJetTags = cms.string('positiveSoftPFElectronRetrainedBJetTags'),
                                                 
    softPFMuonTagInfos       = cms.string('softPFMuons'),     # need to omit the 'TagInfos' part from the label
    softPFElectronTagInfos   = cms.string('softPFElectrons'), # need to omit the 'TagInfos' part from the label
                         
    use_ttbar_filter      = cms.bool(False),
    channel       = cms.InputTag("ttbarselectionproducer")
)

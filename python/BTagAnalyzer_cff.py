import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *


btagana = cms.EDAnalyzer("BTagAnalyzer",
    trackClassifier,
    use_selected_tracks      = cms.bool(True),
    useTrackHistory          = cms.bool(False),
    produceJetProbaTree      = cms.bool(True),
    producePtRelTemplate     = cms.bool(True),
    runSubJets               = cms.bool(False),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    MinPt                    = cms.double(20.0),
    primaryVertexColl        = cms.InputTag('offlinePrimaryVertices'),
    Jets                     = cms.InputTag('ak5PFJets'),
    FatJets                  = cms.InputTag('selectedPatJets'),
    PrunedFatJets            = cms.InputTag('selectedPatJetsCA8PrunedPFPacked'),
    muonCollectionName       = cms.InputTag('muons'),
    patMuonCollectionName    = cms.InputTag('selectedPatMuons'),
    triggerTable             = cms.InputTag('TriggerResults'),
    svComputer               = cms.InputTag('combinedSecondaryVertex'),

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

    simpleIVFSVHighPurBJetTags = cms.string('simpleInclusiveSecondaryVertexHighPurBJetTags'),
    simpleIVFSVHighEffBJetTags = cms.string('simpleInclusiveSecondaryVertexHighEffBJetTags'),
    doubleIVFSVHighEffBJetTags = cms.string('doubleSecondaryVertexHighEffBJetTags'),
    combinedIVFSVBJetTags      = cms.string('combinedInclusiveSecondaryVertexBJetTags'),
    combinedIVFSVPosBJetTags   = cms.string('combinedInclusiveSecondaryVertexPositiveBJetTags'),

    simpleSVHighPurBJetTags     = cms.string('simpleSecondaryVertexHighPurBJetTags'),
    simpleSVNegHighPurBJetTags  = cms.string('simpleSecondaryVertexNegativeHighPurBJetTags'),
    simpleSVHighEffBJetTags     = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    simpleSVNegHighEffBJetTags  = cms.string('simpleSecondaryVertexNegativeHighEffBJetTags'),

    softPFMuonBJetTags        = cms.string('softPFMuonRetrainedBJetTags'),
    softPFMuonNegBJetTags     = cms.string('negativeSoftPFMuonRetrainedBJetTags'),
    softPFMuonPosBJetTags     = cms.string('positiveSoftPFMuonRetrainedBJetTags'),

    softPFElectronBJetTags    = cms.string('softPFElectronRetrainedBJetTags'),
    softPFElectronNegBJetTags = cms.string('negativeSoftPFElectronRetrainedBJetTags'),
    softPFElectronPosBJetTags = cms.string('positiveSoftPFElectronRetrainedBJetTags'),

    softPFMuonTagInfos       = cms.string('softPFMuons'),     # need to omit the 'TagInfos' part from the label
    softPFElectronTagInfos   = cms.string('softPFElectrons'), # need to omit the 'TagInfos' part from the label

    use_ttbar_filter      = cms.bool(False),
    channel       = cms.InputTag("ttbarselectionproducer"),
    TriggerPathNames = cms.vstring(
        "HLT_Jet15U*",
        "HLT_Jet30_v*",
        "HLT_PFJet40_v*",
        "HLT_Jet30U*",
        "HLT_Jet60_v*",
        "HLT_Jet50U*",
        "HLT_Jet80_v*",
        "HLT_PFJet80_v*",
        "HLT_Jet70U*",
        "HLT_Jet110_v*",
        "HLT_Jet100U*",
        "HLT_Jet150_v*",
        "HLT_PFJet140_v*",
        "HLT_Jet140U*",
        "HLT_Jet190_v*",
        "HLT_PFJet200_v*",
        "HLT_Jet240_v*",
        "HLT_PFJet260_v*",
        "HLT_Jet300_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_DiJetAve15U*",
        "HLT_DiJetAve30_v*",
        "HLT_DiPFJetAve40_v*",
        "HLT_DiJetAve30U*",
        "HLT_DiJetAve60_v*",
        "HLT_DiPFJetAve80_v*",
        "HLT_DiJetAve50U*",
        "HLT_DiJetAve80_v*",
        "HLT_DiPFJetAve140_v*",
        "HLT_BTagMu_Jet10U*",
        "HLT_BTagMu_Jet20U*",
        "HLT_BTagMu_DiJet20U*",
        "HLT_BTagMu_DiJet20U_Mu5*",
        "HLT_BTagMu_DiJet20_Mu5*",
        "HLT_BTagMu_DiJet20_L1FastJet_Mu5_v*",
        "HLT_BTagMu_DiJet30U",
        "HLT_BTagMu_DiJet30U_v*",
        "HLT_BTagMu_DiJet30U_Mu5*",
        "HLT_BTagMu_DiJet60_Mu7*",
        "HLT_BTagMu_DiJet40_Mu5*",
        "HLT_BTagMu_DiJet20_L1FastJet_Mu5*",
        "HLT_BTagMu_DiJet80_Mu9*",
        "HLT_BTagMu_DiJet70_Mu5*",
        "HLT_BTagMu_DiJet70_L1FastJet_Mu5*",
        "HLT_BTagMu_DiJet100_Mu9_v*",
        "HLT_BTagMu_DiJet110_Mu5*",
        "HLT_BTagMu_DiJet110_L1FastJet_Mu5*",
        "HLT_BTagMu_Jet300_L1FastJet_Mu5*",
        "HLT_BTagMu_Jet300_Mu5*"
    ),
    TTbarTriggerPathNames = cms.vstring(
        # trigger for ttbar: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Triggers
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",  # 2-electron case
        "HLT_Mu17_Mu8_v*",                                                                          # 2-muon case1
        "HLT_Mu17_TkMu8_v*",                                                                        # 2-muon case2
        "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",                                      # muon + electron case1
        "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                                       # muon + electron case2
    ),
    PFJet80TriggerPathNames = cms.vstring("HLT_PFJet80_v*")
)

import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *
bTagAnalyzerCommon = cms.PSet(
    trackClassifier,
    runOnData                = cms.bool(False),
    runFatJets               = cms.bool(False),
    runSubJets               = cms.bool(False),
    allowJetSkipping         = cms.bool(True),
    useSelectedTracks        = cms.bool(True),
    useTrackHistory          = cms.bool(False),
    produceJetTrackTruthTree = cms.bool(False),
    produceAllTrackTree      = cms.bool(False),
    producePtRelTemplate     = cms.bool(False),
    runEventInfo             = cms.bool(True),
    runJetVariables          = cms.bool(True),
    runQuarkVariables        = cms.bool(False),
    runHadronVariables       = cms.bool(False),
    runGenVariables          = cms.bool(True),
    runPatMuons              = cms.bool(False),
    runTagVariables          = cms.bool(False),
    runTagVariablesSubJets   = cms.bool(False),
    runCSVTagVariables       = cms.bool(False),
    runCSVTagTrackVariables  = cms.bool(False),
    runPrunedDeepFlavourTaggers = cms.bool(False),
    runDeepFlavourTagVariables = cms.bool(False),
    runDeepDoubleXTagVariables = cms.bool(False),
    runDeepBoostedJetTagVariables = cms.bool(False),
    runCSVTagVariablesSubJets = cms.bool(False),
    runPFElectronVariables   = cms.bool(False),
    runPFMuonVariables       = cms.bool(False),
    runCTagVariables         = cms.bool(False),
    fillPU                   = cms.bool(True),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    MinPt                    = cms.double(20.0),
    src                      = cms.InputTag('generator'),
    primaryVertexColl        = cms.InputTag('offlinePrimaryVertices'),
    tracksColl               = cms.InputTag('generalTracks'),
    BranchNamePrefix         = cms.string(''),
    Jets                     = cms.InputTag('selectedPatJets'),
    SubJets                  = cms.VInputTag(),
    SubJetLabels             = cms.vstring(),
    muonCollectionName       = cms.InputTag('muons'),
    patMuonCollectionName    = cms.InputTag('selectedPatMuons'),
    triggerTable             = cms.InputTag('TriggerResults'),
    genParticles             = cms.InputTag('genParticles'),
    prunedGenParticles       = cms.InputTag('prunedGenParticlesBoost'),
    candidates               = cms.InputTag("particleFlow"),
    maxDeltaR                = cms.double(0.4),
    explicitJTA              = cms.bool(False),
    subJetMaxDeltaR          = cms.double(0.4),
    subJetExplicitJTA        = cms.bool(True),
    use_ttbar_filter         = cms.bool(False),
    ttbarproducer            = cms.InputTag("ttbarselectionproducer"),
    clusterTPMap	     = cms.InputTag("tpClusterProducer"),
    rho                      = cms.InputTag("fixedGridRhoFastjetAll"),
    beta                     = cms.double(1.0),
    R0                       = cms.double(0.8),
    maxSVDeltaRToJet         = cms.double(0.7),
    trackPairV0Filter        = cms.PSet(k0sMassWindow = cms.double(0.03)),
    distJetAxis              = cms.double(0.07),
    decayLength              = cms.double(5.0),
    deltaR                   = cms.double(0.3),
    distJetAxisSubJets       = cms.double(0.07),
    decayLengthSubJets       = cms.double(5.0),
    deltaRSubJets            = cms.double(0.3),
    TriggerPathNames = cms.vstring(
        # based on https://cmswbm.cern.ch/cmsdb/servlet/HLTSummary?RUN=297050&NAME=/cdaq/physics/Run2017/2e34/v1.1.1/HLT/V2
        # PF Jets
        "HLT_PFJet40_v*",
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*",
        "HLT_PFJet500_v*",
        "HLT_PFJet550_v*",
        # AK8 PF Jets
        "HLT_AK8PFJet40_v*",
        "HLT_AK8PFJet60_v*",
        "HLT_AK8PFJet80_v*",
        "HLT_AK8PFJet140_v*",
        "HLT_AK8PFJet200_v*",
        "HLT_AK8PFJet260_v*",
        "HLT_AK8PFJet320_v*",
        "HLT_AK8PFJet400_v*",
        "HLT_AK8PFJet450_v*",
        "HLT_AK8PFJet500_v*",
        "HLT_AK8PFJet550_v*",
        # PF HT
        "HLT_PFHT180_v*",
        "HLT_PFHT250_v*",
        "HLT_PFHT370_v*",
        "HLT_PFHT430_v*",
        "HLT_PFHT510_v*",
        "HLT_PFHT590_v*",
        "HLT_PFHT680_v*",
        "HLT_PFHT780_v*",
        "HLT_PFHT890_v*",
        "HLT_PFHT1050_v*",
        # BTagMu
        "HLT_BTagMu_AK4DiJet20_Mu5_v*",
        "HLT_BTagMu_AK4DiJet40_Mu5_v*",
        "HLT_BTagMu_AK4DiJet70_Mu5_v*",
        "HLT_BTagMu_AK4DiJet110_Mu5_v*",
        "HLT_BTagMu_AK4DiJet170_Mu5_v*",
        "HLT_BTagMu_AK4Jet300_Mu5_v*",
        "HLT_BTagMu_AK8DiJet170_Mu5_v*",
        "HLT_BTagMu_AK8Jet300_Mu5_v*",
        # BTagMu Triggers with fix from Xavier
        "HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v*", #triggerIdx=40
        "HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v*", #triggerIdx=41
        "HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v*", #triggerIdx=42
        "HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v*", #triggerIdx=43
        "HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v*", #triggerIdx=44
        "HLT_BTagMu_AK4Jet300_Mu5_noalgo_v*", #triggerIdx+45
        "HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v*",
        "HLT_BTagMu_AK8Jet300_Mu5_noalgo_v*",
        "HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v*"
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

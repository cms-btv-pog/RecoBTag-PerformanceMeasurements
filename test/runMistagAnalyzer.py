##*****************************************************
##*****************************************************
##******** for DATA
##*****************************************************
##*****************************************************
import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *

process = cms.Process("Demo")


process.source = cms.Source(
    "PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(   
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0160/B249DC90-2798-DF11-96D3-001A92811728.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/FACE6E6F-B297-DF11-8C03-002618FDA21D.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/FABF5E30-E897-DF11-8249-002618943836.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/FA906A3C-AA97-DF11-B010-002618943946.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/F650EB46-D697-DF11-A374-00261894396E.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/ECE61B70-B497-DF11-953A-002618943949.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/EA5216EA-AD97-DF11-B7D2-002618943860.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E86D8349-A597-DF11-92FE-002618943946.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E868BC6C-E497-DF11-B07C-002618943967.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E6D0EE42-A897-DF11-AB1E-002618943946.root',
        '/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E2F27541-AA97-DF11-AFD3-00304867BFB0.root', 
    
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_38X_V12::All" # global tag rereco Jun14th
# process.GlobalTag.globaltag = "GR_R_36X_V12B::All" # global tag rereco Jul16th
# process.GlobalTag.globaltag = "GR10_P_V7::All" # global tag prompt reco


process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")

process.load("SimTracker.TrackHistory.TrackHistory_cff")


# for Impact Parameter based taggers
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighPur_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighEffJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetProbabilityBJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetBProbabilityBJetTags_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("JetMETCorrections.Configuration.MCJetCorrections152_cff")


# for Secondary Vertex taggers
process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.secondaryVertexNegativeTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")


# for Soft Muon tagger
process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonBJetTags_cfi")

process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonByPtBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonByPtBJetTags_cfi")


# process.load("SimTracker.TrackHistory.TrackHistory_cff")

# process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")


process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")  

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

process.ak5CaloL2Relative.useCondDB = False
process.ak5CaloL3Absolute.useCondDB = False
process.ak5CaloResidual.useCondDB = False

process.ak5PFL2Relative.useCondDB = False
process.ak5PFL3Absolute.useCondDB = False
process.ak5PFResidual.useCondDB = False

process.ak5JPTL2Relative.useCondDB = False
process.ak5JPTL3Absolute.useCondDB = False
process.ak5JPTResidual.useCondDB = False


process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.PerformanceMeasurements.MistagAnalyzer_cff")


#process.mistag.useTrackHistory = cms.bool(False)
#mistag.jetCorrector    = cms.string('L2L3JetCorrectorIC5Calo')
#process.mcAlgoJetFlavour = cms.Sequence(
#     process.myPartons *
#     process.IC5byRef *
#     process.IC5byValAlgo
# )


# process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree_OldCalib.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree_NEW.root") )
 
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('testtt.root')
)


#-------------------------------------
#Filter for PFJets
process.PFJetsFilter = cms.EDFilter("PFJetSelector",
 src = cms.InputTag("ak5PFJets"),
 cut = cms.string("pt > 10.0 && abs(eta) < 2.4 && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && nConstituents > 1 && chargedHadronEnergyFraction > 0.0 && chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 0.99"),                  
 filter = cms.bool(True)
)
#---------------------------------------
process.load("bTag.CommissioningCommonSetup.caloJetIDFilter_cff")


#Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
                               applyfilter = cms.untracked.bool(True),
                               debugOn = cms.untracked.bool(False),
                               numtrack = cms.untracked.uint32(10),
                               thresh = cms.untracked.double(0.25)
                               )


# Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
# 		      maxAbsZ = cms.double(15), 
 		      maxAbsZ = cms.double(24), 
 		      maxd0 = cms.double(2)	
                     )

#Noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


### from Cristina: calibration JetProb
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_Data14"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_Data14"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)


process.mistag.isData = True
process.mistag.useTrackHistory = False
process.mistag.produceJetProbaTree = False

# process.mistag.Jets = 'ak5CaloJets'
# process.mistag.Jets = 'ak5PFJets'
process.mistag.Jets = 'PFJetsFilter'
# process.mistag.jetCorrector = cms.string('ak5CaloL2L3')
process.mistag.jetCorrector = cms.string('ak5PFL2L3')
# process.ak5JetTracksAssociatorAtVertex.jets = "ak5CaloJets"
# process.ak5JetTracksAssociatorAtVertex.jets = "ak5PFJets"
process.ak5JetTracksAssociatorAtVertex.jets = "PFJetsFilter"
# process.softMuonTagInfos.jets = "ak5CaloJets"
# process.softMuonTagInfos.jets = "ak5PFJets"
process.softMuonTagInfos.jets = "PFJetsFilter"


process.p = cms.Path(
#$$        process.myPartons*process.AK5Flavour
#$$ usefull ?
        process.ak5PFJetsL2L3
#$$
        *process.PFJetsFilter
        *process.HBHENoiseFilter
        *process.noscraping
        *process.primaryVertexFilter
	*process.ak5JetTracksAssociatorAtVertex
	*process.btagging
	*process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags
	*process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
	*process.secondaryVertexTagInfos*process.simpleSecondaryVertexHighEffBJetTags*process.simpleSecondaryVertexHighPurBJetTags
	*process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeHighEffBJetTags*process.simpleSecondaryVertexNegativeHighPurBJetTags
        *process.combinedSecondaryVertexNegativeBJetTags
	#*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
	*process.negativeSoftLeptonByPtBJetTags*process.positiveSoftLeptonByPtBJetTags	
	*process.jetBProbabilityBJetTags*process.negativeOnlyJetBProbabilityJetTags
	*process.mistag
	)


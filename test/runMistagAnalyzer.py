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
        '/store/data/Run2010A/JetMET/RECO/v4/000/143/193/E2686140-2CAB-DF11-AAB9-000423D94908.root',
        '/store/data/Run2010A/JetMET/RECO/v4/000/143/193/C0D6B836-33AB-DF11-BC87-0030487CD906.root',
        '/store/data/Run2010A/JetMET/RECO/v4/000/143/193/0E856940-2CAB-DF11-A2A9-0030487CD6DA.root',
        '/store/data/Run2010A/JetMET/RECO/v4/000/143/192/12F92A25-31AB-DF11-8655-003048F118AA.root',
        '/store/data/Run2010A/JetMET/RECO/v4/000/143/191/E626B382-2BAB-DF11-8C8D-0030487CD178.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_36X_V12::All" # global tag rereco Jun14th
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
 cut = cms.string("pt > 10.0 && abs(eta) < 2.5 && neutralHadronEnergyFraction < 1.0 && neutralEmEnergyFraction < 1.0 && nConstituents > 1 && chargedHadronEnergyFraction > 0.0 && chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 1.0"),
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


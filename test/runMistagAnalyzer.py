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
       # 'file:FCCB5194-257C-E011-83C3-0024E876804B.root'
         '/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0005/FCCB5194-257C-E011-83C3-0024E876804B.root',
         '/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0005/F89858D5-417C-E011-8F20-00266CF256CC.root',
         '/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0005/F4F7BC71-137C-E011-A5D1-00151796C07C.root',
         '/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0005/E4C29D0F-397C-E011-A9C7-0024E87683F8.root',
         '/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0005/E415D3B0-FB7B-E011-9785-0024E876884D.root'

#         '/store/data/Run2010A/JetMET/AOD/Apr21ReReco-v1/0006/08D5B595-0671-E011-A463-001A928116DE.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "FT_R_42_V10A::All" # global tag rerereco 2010
# process.GlobalTag.globaltag = "GR_R_42_V14::All" # global tag rereco 2011
# process.GlobalTag.globaltag = "GR_P_V20::All" # global tag promptreco 2011
process.GlobalTag.globaltag = "GR_R_42_V12::All" # global tag promptreco 2011, with FastJet
#process.GlobalTag.globaltag = "DESIGN42_V11::All" # MC perfect alignment , with FastJet
#process.GlobalTag.globaltag = "START42_V12::All"  # MC startup alignment , with FastJet
#process.GlobalTag.globaltag = "MC_42_V12::All"    # MC default Alignment , with FastJet



##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True









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

# for Secondary Vertex taggers
process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.secondaryVertexNegativeTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveES_cfi")

# for Soft Muon tagger
process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonBJetTags_cfi")

process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonByPtBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonByPtBJetTags_cfi")

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")  

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
#$$
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True
#$$

process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.PerformanceMeasurements.MistagAnalyzer_cff")


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
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")

# Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
 		      maxAbsZ = cms.double(24), 
 		      maxd0 = cms.double(2)	
                     )

### from Cristina JP calibration for cmsRun only : 
# from CondCore.DBCommon.CondDBCommon_cfi import *
# process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
# process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_Data_2010_41X.db")
# # process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_Data_2011_41X.db")
# process.es_prefer_trackProbabilityFakeCond = cms.ESPrefer("PoolDBESSource","trackProbabilityFakeCond")

### from Cristina JP calibration for crab only: 
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_2010Data_v1_offline"),
#        tag = cms.string("TrackProbabilityCalibration_2D_2011Data_v1_offline"),
       connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_2010Data_v1_offline"),
#        tag = cms.string("TrackProbabilityCalibration_3D_2011Data_v1_offline"),
       connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
)

# process.GlobalTag.toGet = cms.VPSet(
#   cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
#        tag = cms.string("TrackProbabilityCalibration_2D_2010_v1_mc"),
#        connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
#   cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
#        tag = cms.string("TrackProbabilityCalibration_3D_2010_v1_mc"),
#        connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
# )


#---------------------------------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )
 
process.mistag.isData              = True
process.mistag.useTrackHistory     = False
process.mistag.produceJetProbaTree = True
process.mistag.triggerTable = 'TriggerResults::HLT'

# process.mistag.Jets = 'ak5PFJets'
process.mistag.Jets = 'PFJetsFilter'
#process.mistag.jetCorrector = cms.string('ak5PFL2L3')
process.mistag.jetCorrector = cms.string('ak5PFL1FastL2L3')
#process.mistag.jetCorrector = cms.string('ak5PFL1FastL2L3Residual')

# process.ak5JetTracksAssociatorAtVertex.jets = "ak5PFJets"
process.ak5JetTracksAssociatorAtVertex.jets = "PFJetsFilter"
# process.softMuonTagInfos.jets = "ak5PFJets"
process.softMuonTagInfos.jets = "PFJetsFilter"


process.p = cms.Path(
#$$
#$$         process.kt6PFJets
#$$        *process.PFJetsFilter
        process.kt6PFJets  
	*process.ak5PFJets  
        *process.PFJetsFilter
        *process.noscraping
        *process.offlinePrimaryVertices 
        *process.primaryVertexFilter
	*process.ak5JetTracksAssociatorAtVertex
	*process.btagging
	*process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags
	*process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
	*process.secondaryVertexTagInfos*process.simpleSecondaryVertexHighEffBJetTags*process.simpleSecondaryVertexHighPurBJetTags
	*process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeHighEffBJetTags*process.simpleSecondaryVertexNegativeHighPurBJetTags
        *process.combinedSecondaryVertexNegativeBJetTags*process.combinedSecondaryVertexPositiveBJetTags
	#*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
	*process.negativeSoftLeptonByPtBJetTags*process.positiveSoftLeptonByPtBJetTags	
	*process.jetBProbabilityBJetTags*process.negativeOnlyJetBProbabilityJetTags
	*process.mistag
	)


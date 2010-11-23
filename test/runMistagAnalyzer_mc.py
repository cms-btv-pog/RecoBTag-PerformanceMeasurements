##*****************************************************
##*****************************************************
##******** for MC
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
        '/store/mc/Fall10/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/START38_V12-v1/0000/FE60B60E-E7DB-DF11-B82F-0030486790C0.root',
        '/store/mc/Fall10/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/START38_V12-v1/0000/FCE36495-FEDB-DF11-8132-0018F3D096DC.root',
        '/store/mc/Fall10/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/START38_V12-v1/0000/FC820950-E7DB-DF11-89C6-00304867D838.root',
        '/store/mc/Fall10/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/START38_V12-v1/0000/FC445F23-F7DB-DF11-926B-001A92810AC0.root',
        '/store/mc/Fall10/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/START38_V12-v1/0000/FC2FC680-E8DB-DF11-9588-0030486790B0.root'

#         '/store/mc/Fall10/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/AODSIM/START38_V12-v1/0006/F2990AEF-EBCA-DF11-B191-485B39800C0C.root',
#         '/store/mc/Fall10/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/AODSIM/START38_V12-v1/0006/E6B32F6B-E7CA-DF11-8429-00261834B5D2.root',
#         '/store/mc/Fall10/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/AODSIM/START38_V12-v1/0006/E27A0FBD-FCCA-DF11-A998-E0CB4E19F987.root',
#         '/store/mc/Fall10/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/AODSIM/START38_V12-v1/0006/D077A8B4-EACA-DF11-AFF8-485B39800C03.root',
#         '/store/mc/Fall10/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/AODSIM/START38_V12-v1/0006/C6ADB927-08CB-DF11-B34A-E0CB4E29C4D4.root'
  )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START38_V13::All" # use a valid global tag here!


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


#Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
 		      maxAbsZ = cms.double(24), 
 		      maxd0 = cms.double(2)	
                     )

#Noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


### from Cristina: calibration JetProb
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_MCQpt30"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_MCQpt30"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)


#---------------------------------------
# process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )
 
process.mistag.isData = False
process.mistag.useTrackHistory = True
process.mistag.produceJetProbaTree = False
#process.mistag.triggerTable = 'TriggerResults::HLT'
#process.mistag.triggerTable = 'TriggerResults::REDIGI38X'
process.mistag.triggerTable = 'TriggerResults::REDIGI38XTP'
#---------------------------------------


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
        process.myPartons*process.AK5Flavour
#$$ usefull ?
        *process.ak5PFJetsL2L3
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
	*process.jetProbabilityBJetTags
	*process.jetBProbabilityBJetTags*process.negativeOnlyJetBProbabilityJetTags
	*process.mistag
	)


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
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/0882B650-69C9-DF11-AE2C-E0CB4E19F990.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/10F7AAC0-BFC8-DF11-A198-E0CB4E4408C6.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/1C007AFC-C2C8-DF11-82E3-E0CB4E553640.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/22504C49-2FC9-DF11-9ABB-0022198F5AEB.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/30A41B2E-32C9-DF11-AE60-90E6BA0D09B4.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/38B48C11-25C9-DF11-BEA3-90E6BA442F3B.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/44C81239-2AC9-DF11-8810-485B39800C0F.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/64DAE9B4-9BC8-DF11-A1E1-001A4BA8AF26.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/64E906EB-1FC9-DF11-B679-00261834B5A7.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/741E9400-51C9-DF11-AF8E-001A4BA950A2.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/7A4987F8-C6C8-DF11-85F2-E0CB4E29C513.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/7AA1AB57-68C9-DF11-8485-E0CB4EA0A937.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/7E745D9A-24C9-DF11-9A69-90E6BA19A1F9.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/8A675863-20C9-DF11-B79F-90E6BA442F1F.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/8AA8CBCB-1CC9-DF11-A492-E0CB4EA0A91E.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/94AC1FE6-AAC8-DF11-9DD2-003048678948.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/A86C9646-A5C8-DF11-98BD-90E6BA442F31.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/AA6CBC07-20C9-DF11-A584-00261834B51D.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/BC6AB206-61C9-DF11-A2EF-001EC9D7F68B.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/BE41F198-20C9-DF11-B1DA-E0CB4E55367B.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/C0A6641C-AAC8-DF11-A844-E0CB4E1A1167.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/D256CB9D-2CC9-DF11-897D-485B39800C16.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/E675BFC6-1DC9-DF11-8B00-0030487CDA66.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/EC0B0ABD-A5C8-DF11-B3A1-E0CB4E19F9A2.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/EC1BBE87-20C9-DF11-9E59-E0CB4E19F972.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/EC961E68-69C9-DF11-8431-E0CB4E4408CB.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/EEACD7FA-1EC9-DF11-85C4-E0CB4E19F9AC.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/F20D663E-BBC8-DF11-8BBC-485B39800C2B.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/F2CEB55F-51C9-DF11-819B-001A4BA900B8.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/F6C04156-C0C8-DF11-91B1-485B39800C30.root",
"rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Fall10/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/FE5A8C92-BFC8-DF11-8CD6-E0CB4E5536F9.root"
  )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START38_V12::All" # use a valid global tag here!
# process.GlobalTag.globaltag = "START36_V10::All" # use a valid global tag here!


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


# process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree_EXT.root") )
 
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


#Filter for good primary vertex
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
       tag = cms.string("TrackProbabilityCalibration_2D_MCQpt30"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_MCQpt30"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)


process.mistag.isData = False
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
        process.myPartons*process.AK5Flavour
#$$ usefull ?
        *process.ak5PFJetsL2L3
#$$
        *process.PFJetsFilter
        #*process.noscraping
        #*process.HBHENoiseFilter
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


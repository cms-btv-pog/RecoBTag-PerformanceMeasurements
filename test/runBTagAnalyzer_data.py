##*****************************************************
##*****************************************************
##******** for DATA
##*****************************************************
##*****************************************************
import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *

process = cms.Process("BTagAna")

process.source = cms.Source(
    "PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(   
# data from /Jet/Run2012A-PromptReco-v1/AOD
#        'file:FE2E21F0-0583-E111-8DA9-001D09F241F0.root'
# data from /BTag/Run2012A-PromptReco-v1/AOD
       #'file:000A540B-76E7-E111-A557-003048C69406.root'
       #'file:FC834104-FF2E-E211-AF15-003048F1BF66.root'
       'file:5AC9D874-D6E8-E111-A768-00261894387D.root'
  )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_53_V19::All"
# process.GlobalTag.globaltag = "GR_R_53_V2B::All"




##############################################
#### Add new retrained CSV b tagging algo ####
##############################################
#process.load("CondCore.DBCommon.CondDBSetup_cfi")

#process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
#   process.CondDBSetup,
#   timetype = cms.string('runnumber'),
#   toGet = cms.VPSet(cms.PSet(
#      record = cms.string('BTauGenericMVAJetTagComputerRcd'),
#                tag = cms.string('MVAComputerContainer_JetTags_v3_offline')
#   )),
#   connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
#   BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
#)

#process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord")







##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#$$
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True
#$$

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")

process.load("SimTracker.TrackHistory.TrackHistory_cff")

############################################################################
# for Impact Parameter based taggers
############################################################################
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighPur_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighEffJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetProbabilityBJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetBProbabilityBJetTags_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

############################################################################
# for Secondary Vertex taggers
############################################################################
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


############################################################################
# for Inclusive Secondary Vertexing
############################################################################
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
process.load("RecoBTag.SecondaryVertex.bToCharmDecayVertexMerger_cfi")

process.combinedInclusiveSecondaryVertexPositiveBJetTags = process.combinedInclusiveSecondaryVertexBJetTags.clone()
process.combinedInclusiveSecondaryVertexPositiveBJetTags.jetTagComputer = cms.string('combinedSecondaryVertexPositive')


############################################################################
# for the Retrained CSV
############################################################################
#process.combinedSecondaryVertexRetrained = process.combinedSecondaryVertex.clone()
#process.combinedSecondaryVertexRetrained.calibrationRecords = cms.vstring(
#      'CombinedSVRetrain2RecoVertex', 
#      'CombinedSVRetrain2PseudoVertex', 
#      'CombinedSVRetrain2NoVertex'
#      )
#process.combinedSecondaryVertexRetrainedBJetTags = process.combinedSecondaryVertexBJetTags.clone()
#process.combinedSecondaryVertexRetrainedBJetTags.jetTagComputer = cms.string('combinedSecondaryVertexRetrained')

#process.combinedSecondaryVertexNegativeRetrained = process.combinedSecondaryVertexNegative.clone()
#process.combinedSecondaryVertexNegativeRetrained.calibrationRecords = cms.vstring(
#      'CombinedSVRetrain2RecoVertex', 
#      'CombinedSVRetrain2PseudoVertex', 
#      'CombinedSVRetrain2NoVertex'
#      )
#process.combinedSecondaryVertexNegativeRetrainedBJetTags = process.combinedSecondaryVertexNegativeBJetTags.clone()
#process.combinedSecondaryVertexNegativeRetrainedBJetTags.jetTagComputer = cms.string('combinedSecondaryVertexNegativeRetrained')

#process.combinedSecondaryVertexPositiveRetrained = process.combinedSecondaryVertexPositive.clone()
#process.combinedSecondaryVertexPositiveRetrained.calibrationRecords = cms.vstring(
#      'CombinedSVRetrain2RecoVertex', 
#      'CombinedSVRetrain2PseudoVertex', 
#      'CombinedSVRetrain2NoVertex'
#      )
#process.combinedSecondaryVertexPositiveRetrainedBJetTags = process.combinedSecondaryVertexPositiveBJetTags.clone()
#process.combinedSecondaryVertexPositiveRetrainedBJetTags.jetTagComputer = cms.string('combinedSecondaryVertexPositiveRetrained')





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

process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff")


process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('testtt.root')
)


#-------------------------------------
#Produce PFJets from PF2PAT

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

# process.load("CommonTools.ParticleFlow.Sources.source_ZtoMus_DBS_cfi")
runOnMC = False

from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PF2PAT"  # to have only PF2PAT
jetAlgo="AK5"

usePFnoPU = True # before any top projection


# JEC set
jecSetBase = jetAlgo
jecSetPF = jecSetBase + 'PF'
if usePFnoPU:
  jecSetPF += 'chs'


usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix, 
jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']) ) 


applyPostfix(process,"pfPileUp",postfix).checkClosestZVertex = cms.bool(False) 

from PhysicsTools.PatAlgos.tools.metTools import *

process.pfPileUpPF2PAT.Enable = True
process.pfPileUpPF2PAT.Vertices = cms.InputTag('goodOfflinePrimaryVertices')

process.load('RecoJets.Configuration.RecoJets_cff')
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets

process.kt6PFJets               = kt4PFJets.clone()
process.kt6PFJets.rParam        = 0.6     
#process.kt6PFJets.src           = cms.InputTag('pfNoElectron'+postfix)
process.kt6PFJets.Rho_EtaMax    = cms.double( 4.4)
process.kt6PFJets.doRhoFastjet  = True
process.kt6PFJets.doAreaFastjet = True
#process.kt6PFJets.voronoiRfact  = 0.9

#process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJets", "rho")
getattr(process,"patJetCorrFactors"+postfix).rho = cms.InputTag("kt6PFJets", "rho")

#-------------------------------------
#Redo the primary vertex

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )
    
#-------------------------------------
#JetID
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi')
process.load('PhysicsTools.SelectorUtils.jetIDSelector_cfi')
process.jetIDSelector.version = cms.string('PURE09')


#-------------------------------------
#Filter for PFJets
process.PATJetsFilter = cms.EDFilter("PATJetSelector",    
                                    src = cms.InputTag("selectedPatJetsPF2PAT"),
                                    cut = cms.string("pt > 10.0 && abs(eta) < 2.5 && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && nConstituents > 1 && chargedHadronEnergyFraction > 0.0 && chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 0.99"),
                                    #filter = cms.bool(True)
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
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
 		      maxAbsZ = cms.double(24), 
 		      maxd0 = cms.double(2)	
                     )


### JP calibration for cmsRun only : 
# from CondCore.DBCommon.CondDBCommon_cfi import *
# process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
# process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_Data_2011_41X.db")
# process.es_prefer_trackProbabilityFakeCond = cms.ESPrefer("PoolDBESSource","trackProbabilityFakeCond")

### JP calibration for crab only: 
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)


#---------------------------------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree_data.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )
 
process.btagana.isData              =True 
process.btagana.use_selected_tracks = False
process.btagana.useTrackHistory     = False
process.btagana.produceJetProbaTree = True
process.btagana.producePtRelTemplate = False
process.btagana.triggerTable = 'TriggerResults::HLT' # Data and MC
#---------------------------------------


process.AK5byRef.jets = "selectedPatJetsPF2PAT"
process.btagana.jetCorrector = cms.string('ak5PFL1FastL2L3Residual')

#process.btagana.Jets = 'selectedPatJetsPF2PAT'
#process.ak5JetTracksAssociatorAtVertex.jets = "selectedPatJetsPF2PAT"
#process.softMuonTagInfos.jets = "selectedPatJetsPF2PAT"


process.btagana.Jets = 'PATJetsFilter'
process.ak5JetTracksAssociatorAtVertex.jets = "PATJetsFilter"
process.softMuonTagInfos.jets = "PATJetsFilter"




#---------------------------------------
# trigger selection !
# import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
# process.JetHLTFilter = hlt.triggerResultsFilter.clone(
#     triggerConditions = cms.vstring(
#	"HLT_PFJet80_v*"
#     ),
#     hltResults = cms.InputTag("TriggerResults","","HLT"),
#     l1tResults = cms.InputTag( "" ),
#     throw = cms.bool( False ) #set to false to deal with missing triggers while running over different trigger menus
#     )
#---------------------------------------






process.p = cms.Path(
#$$
#$$        process.JetHLTFilter
#$$        *process.kt6PFJets  
#$$
        #process.kt6PFJets
        process.noscraping
        *process.primaryVertexFilter
        *process.offlinePrimaryVertices 
	*process.goodOfflinePrimaryVertices
	*getattr(process,"patPF2PATSequence"+postfix)
        *process.PATJetsFilter
        *process.inclusiveVertexing*process.inclusiveMergedVerticesFiltered*process.bToCharmDecayVertexMerged
#$$	*process.myPartons*process.AK5Flavour
#$$
	*process.ak5JetTracksAssociatorAtVertex
	*process.btagging
	*process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags
	*process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
	*process.secondaryVertexTagInfos*process.simpleSecondaryVertexHighEffBJetTags*process.simpleSecondaryVertexHighPurBJetTags
	*process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeHighEffBJetTags*process.simpleSecondaryVertexNegativeHighPurBJetTags
        *process.combinedSecondaryVertexNegativeBJetTags*process.combinedSecondaryVertexPositiveBJetTags
	#*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
	*process.negativeSoftLeptonByPtBJetTags*process.positiveSoftLeptonByPtBJetTags	
	*process.negativeOnlyJetBProbabilityJetTags*process.positiveOnlyJetBProbabilityJetTags
        #*process.combinedSecondaryVertexRetrainedBJetTags*process.combinedSecondaryVertexNegativeRetrainedBJetTags*process.combinedSecondaryVertexPositiveRetrainedBJetTags
        *process.inclusiveSecondaryVertexFinderTagInfosFiltered*process.simpleInclusiveSecondaryVertexHighEffBJetTags*process.simpleInclusiveSecondaryVertexHighPurBJetTags
        *process.doubleSecondaryVertexHighEffBJetTags
        *process.inclusiveSecondaryVertexFinderTagInfos*process.combinedInclusiveSecondaryVertexBJetTags*process.combinedInclusiveSecondaryVertexPositiveBJetTags
	*process.btagana
	)
	
	## Output Module Configuration (expects a path 'p')





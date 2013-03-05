##*****************************************************
##*****************************************************
##******** for DATA & MC
##*****************************************************
##*****************************************************
import sys

# Import configurations

import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#Passing_Command_Line_Arguments_T
# setup 'standard'  options
options = VarParsing.VarParsing ('analysis')
## setup any defaults you want
options.register('FlagData', 0,
                  VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"flag 1 for Data, 0 for MC")
options.register('GlobalTagData', 'FT_53_V6_AN2',
                  VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"global tag for Data")
options.maxEvents = -1



# Hack to be able to pass multiple arguments from multicrab.cfg 
# from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/zbb_louvain/test/patTuple_llbb_cmssw444_cfg_dataMC.py?revision=1.24&view=markup
print sys.argv
if len(sys.argv) > 0:
    last = sys.argv.pop()
    sys.argv.extend(last.split(":"))
    print sys.argv

# get and parse the command line arguments
options.parseArguments()

# Use the options
runOnMC = True # True if MC and False if Data 
if options.FlagData :
  runOnMC = False

print "run On MC? ", runOnMC




from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *

process = cms.Process("BTagAna")

if runOnMC:
  process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
       'file:/opt/sbg/data/data3/cms/ccollard/files/test/FE917359-9CF5-E111-A93F-848F69FD4DCC.root'  #ttbar
    )
  )
else:
  process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/opt/sbg/data/data3/cms/ccollard/files/test/FEEE5F6A-26DA-E111-B08C-00266CFAE8D0.root' ## double-e
    )
  )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if runOnMC:
  process.GlobalTag.globaltag = "START53_V7F::All"  ## MC from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagValidationSamples#MC_samples_at_8_TeV
else:
  process.GlobalTag.globaltag = options.GlobalTagData+'::All' ## Data from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagValidationSamples#Data_samples

print "GlobalTag : ", process.GlobalTag.globaltag





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

##process.outpath = cms.EndPath(process.out)


#-------------------------------------
#Produce PFJets from PF2PAT

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PF2PAT"  # to have only PF2PAT
jetAlgo="AK5"

usePFnoPU = True # before any top projection


# JEC set
jecSetBase = jetAlgo
jecSetPF = jecSetBase + 'PF'
if usePFnoPU:
  jecSetPF += 'chs'


if runOnMC: 
 usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix, 
 jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']) ) 
else:
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

####################################
#  HBHE noise filter
####################################
#Noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


#Filter for good primary vertex
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
 		      maxAbsZ = cms.double(24), 
 		      maxd0 = cms.double(2)	
                     )


### from Cristina JP calibration for cmsRun only : 
# from CondCore.DBCommon.CondDBCommon_cfi import *
# process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
# process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_Data_2011_41X.db")
# process.es_prefer_trackProbabilityFakeCond = cms.ESPrefer("PoolDBESSource","trackProbabilityFakeCond")

### from JP calibration for crab only: 
if runOnMC:
 process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
 )
else:
 process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
 )



#---------------------------------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )
 
if runOnMC:
  process.btagana.isData              = False  ## MC
else:
  process.btagana.isData              = True ## Data
process.btagana.use_selected_tracks = False
process.btagana.useTrackHistory     = False
process.btagana.produceJetProbaTree = True
process.btagana.producePtRelTemplate = False
process.btagana.triggerTable = 'TriggerResults::HLT' # Data and MC
process.btagana.use_ttbar_filter     = True
#---------------------------------------


process.AK5byRef.jets = "selectedPatJetsPF2PAT"
if runOnMC:
  process.btagana.jetCorrector = cms.string('ak5PFL1FastL2L3') ## MC
else:
  process.btagana.jetCorrector = cms.string('ak5PFL1FastL2L3Residual')  ## Data

#process.btagana.Jets = 'selectedPatJetsPF2PAT'
#process.ak5JetTracksAssociatorAtVertex.jets = "selectedPatJetsPF2PAT"
#process.softMuonTagInfos.jets = "selectedPatJetsPF2PAT"


process.btagana.Jets = 'PATJetsFilter'
process.ak5JetTracksAssociatorAtVertex.jets = "PATJetsFilter"
process.softMuonTagInfos.jets = "PATJetsFilter"



#---------------------------------------


process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionFilter_cfi")
process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionProducer_cfi")

process.ttbarselectionproducer.isData              = False
process.ttbarselectionfilter.select_ee             = True
process.ttbarselectionfilter.select_mumu           = True
process.ttbarselectionfilter.select_emu            = True

####################################
#  top projections in PF2PAT:
####################################
getattr(process,"pfNoPileUp"  +postfix).enable = True
getattr(process,"pfNoMuon"    +postfix).enable = False
getattr(process,"pfNoElectron"+postfix).enable = False
getattr(process,"pfNoTau"     +postfix).enable = False
getattr(process,"pfNoJet"     +postfix).enable = False


# change the cone size of electron isolation to 0.3 as default.
applyPostfix(process,"pfIsolatedMuons",postfix).isolationValueMapsCharged = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03PF2PAT' ) )
applyPostfix(process,"pfIsolatedMuons",postfix).isolationValueMapsNeutral = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03PF2PAT' ), cms.InputTag( 'muPFIsoValueGamma03PF2PAT' ) )
applyPostfix(process,"pfIsolatedMuons",postfix).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03PF2PAT' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfNeutralHadrons = cms.InputTag( 'muPFIsoValueNeutral03PF2PAT' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfPhotons = cms.InputTag( 'muPFIsoValueGamma03PF2PAT' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfChargedHadrons = cms.InputTag( 'muPFIsoValueCharged03PF2PAT' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03PF2PAT' )
applyPostfix(process,"pfIsolatedMuons",postfix).combinedIsolationCut = cms.double(9999.)
applyPostfix(process,"pfIsolatedMuons",postfix).isolationCut = cms.double(9999.)

# change the cone size of electron isolation to 0.3 as default.
process.pfIsolatedElectronsPF2PAT.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"))
process.pfIsolatedElectronsPF2PAT.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT")
process.pfIsolatedElectronsPF2PAT.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPF2PAT"), cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT"))

process.pfElectronsPF2PAT.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"))
process.pfElectronsPF2PAT.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT" )
process.pfElectronsPF2PAT.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag( "elPFIsoValueNeutral03PFIdPF2PAT"), cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT"))

process.patElectronsPF2PAT.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPF2PAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPF2PAT"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT")
        )

applyPostfix(process,"pfIsolatedElectrons",postfix).combinedIsolationCut = cms.double(9999.)
applyPostfix(process,"pfIsolatedElectrons",postfix).isolationCut = cms.double(9999.)


process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
#Electron ID
process.patElectronsPF2PAT.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPF2PAT.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePF2PAT.replace( process.patElectronsPF2PAT, process.eidMVASequence * process.patElectronsPF2PAT )

#Convesion Rejection
# this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
process.patConversionsPF2PAT = cms.EDProducer("PATConversionProducer",
                                             electronSource = cms.InputTag("selectedPatElectronsPF2PAT")
                                             )
process.patPF2PATSequencePF2PAT += process.patConversionsPF2PAT


###-------
if runOnMC:
  process.GenHistory = cms.Sequence(
        process.myPartons*process.AK5Flavour
    )
else:
  process.GenHistory = cms.Sequence()
###-------



# process.jetProbabilityMixed = cms.ESProducer("JetProbabilityESProducer",
#     impactParameterType = cms.int32(0), ## 0 = 3D, 1 = 2D
#     deltaR = cms.double(0.3),
#     maximumDistanceToJetAxis = cms.double(0.07),
#     trackIpSign = cms.int32(0), ## 0 = use both, 1 = positive only, -1 = negative only
#     minimumProbability = cms.double(0.005),
#     maximumDecayLength = cms.double(5.0),
#     trackQualityClass = cms.string("any")
# )
# 
# #-------------------------------------
# #Jet Probability
# process.jetBProbabilityMixed = cms.ESProducer("JetBProbabilityESProducer",
#     impactParameterType = cms.int32(0), ## 0 = 3D, 1 = 2D
#     deltaR = cms.double(-1.0), ## use cut from JTA
#     maximumDistanceToJetAxis = cms.double(0.07),
#     trackIpSign = cms.int32(0), ## 0 = use both, 1 = positive only, -1 = negative only
#     minimumProbability = cms.double(0.005),
#     numberOfBTracks = cms.uint32(4),
#     maximumDecayLength = cms.double(5.0),
#     trackQualityClass = cms.string("any")
# )
# 
# process.jetProbabilityBJetTags.jetTagComputer  = 'jetProbabilityMixed'
# process.jetBProbabilityBJetTags.jetTagComputer = 'jetBProbabilityMixed'


#---------------------------------------
# trigger selection !
# process.JetHLTFilter = hlt.triggerResultsFilter.clone(
#     triggerConditions = cms.vstring(
# # 	 " HLT_Jet80_v6"
# 	 " HLT_Jet110_v6"
#     ),
#     hltResults = cms.InputTag("TriggerResults","","HLT"),
#     l1tResults = cms.InputTag( "" ),
#     throw = cms.bool( False) #set to false to deal with missing triggers while running over different trigger menus
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
        *process.HBHENoiseFilter
	*getattr(process,"patPF2PATSequence"+postfix)
        *process.PATJetsFilter
        *process.inclusiveVertexing*process.inclusiveMergedVerticesFiltered*process.bToCharmDecayVertexMerged
        *process.GenHistory
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
        *process.ttbarselectionproducer
        *process.ttbarselectionfilter
	*process.btagana
	)
	
	## Output Module Configuration (expects a path 'p')



# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)



import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import copy

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('outFilename', 'JetTree',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('reportEvery', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('usePFchs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
)
options.register('mcGlobalTag', 'PH2_1K_FB_V6', #2017: DES17_62_V8; 2019: DES19_62_V8, PH1_1K_FB_V2; 20123: DES23_62_V1
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'GR_R_62_V1',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)
options.register('runSubJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)
options.register('processStdAK5Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Process standard AK5 jets"
)
options.register('producePtRelTemplate', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Produce PtRel template"
)
options.register('fatJetPtMin', 150.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum pT for fat jets (default is 150 GeV)"
)
options.register('useTTbarFilter', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use TTbar filter"
)
options.register('usePVSorting', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PV sorting"
)

options.register('useBetterPVSorting', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use better PV sorting"
)

options.register('useBetterPVSortingPF', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use better PV sorting PF"
)

options.register('changeMinNumberOfHits', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Change minimum number of tracker hits"
)
options.register('minNumberOfHits', 8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Minimum number of tracker hits"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if not options.usePFchs:
    jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK5 = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK5[1].append('L2L3Residual')
    jetCorrectionsAK7[1].append('L2L3Residual')

## b-tag infos
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','secondaryVertexNegativeTagInfos','softMuonTagInfos'
     #,'softPFMuonsTagInfos','softPFElectronsTagInfos'
     ,'inclusiveSecondaryVertexFinderTagInfos'
     #,'inclusiveSecondaryVertexFinderFilteredTagInfos'
]
## b-tag discriminators
bTagDiscriminators = ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags'
     #,'negativeOnlyJetBProbabilityJetTags','negativeOnlyJetProbabilityJetTags'
     ,'negativeTrackCountingHighEffJetTags','negativeTrackCountingHighPurJetTags'
     #,'positiveOnlyJetBProbabilityJetTags','positiveOnlyJetProbabilityJetTags'
    ,'simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','simpleSecondaryVertexNegativeHighEffBJetTags'
    ,'simpleSecondaryVertexNegativeHighPurBJetTags','combinedSecondaryVertexBJetTags'
    ,'combinedSecondaryVertexV2BJetTags'
    #,'combinedSecondaryVertexIVFV2Phase1MidPUBJetTags'
    #,'combinedSecondaryVertexPositiveBJetTags',
    #,'combinedSecondaryVertexV1BJetTags','combinedSecondaryVertexV1PositiveBJetTags'
    #,'combinedSecondaryVertexNegativeBJetTags'
    #,'combinedSecondaryVertexV1NegativeBJetTags'
    #,'softPFMuonBJetTags','positiveSoftPFMuonBJetTags','negativeSoftPFMuonBJetTags'
    #,'softPFElectronBJetTags','positiveSoftPFElectronBJetTags','negativeSoftPFElectronBJetTags','simpleInclusiveSecondaryVertexHighEffBJetTags'
    #,'simpleInclusiveSecondaryVertexHighPurBJetTags','doubleSecondaryVertexHighEffBJetTags','combinedInclusiveSecondaryVertexBJetTags'
    #,'combinedInclusiveSecondaryVertexPositiveBJetTags','combinedSecondaryVertexSoftPFLeptonV1BJetTags','positiveCombinedSecondaryVertexSoftPFLeptonV1BJetTags'
    #,'negativeCombinedSecondaryVertexSoftPFLeptonV1BJetTags'
]
bTagDiscriminatorsSubJets = copy.deepcopy(bTagDiscriminators)
#bTagDiscriminatorsSubJets.remove('doubleSecondaryVertexHighEffBJetTags')

process = cms.Process("BTagAna")

## PV collection
pvCollection = 'goodOfflinePrimaryVertices'
if options.useBetterPVSorting and options.usePVSorting:
    print "Conficting PV Sorting configuration, choosing the Better PV Sorting"
    options.usePVSorting = False
if options.useBetterPVSortingPF and options.usePVSorting:
    print "Conficting PV Sorting configuration, choosing the Better PV Sorting"
    options.usePVSorting = False
if options.useBetterPVSorting and options.useBetterPVSortingPF:
    print "Both Better PV sorting and Better PV sorting PF-based required. Choosing the PF-based sorting"
    options.useBetterPVSorting = False
if options.usePVSorting:
    pvCollection = 'sortedGoodOfflinePrimaryVertices'

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10

## Input files
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        # /RelValTTbar_14TeV/CMSSW_6_2_0_SLHC7-DES17_62_V8_UPG2017-v2/GEN-SIM-RECO
        #'/store/relval/CMSSW_6_2_0_SLHC7/RelValTTbar_14TeV/GEN-SIM-RECO/DES17_62_V8_UPG2017-v2/00000/04D56E72-4390-E311-A322-02163E008D7C.root'

        # /PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEM2019Upg14DR-final_phase1_PU50bx25_DES19_62_V8-v1/AODSIM
        ##'/store/mc/GEM2019Upg14DR/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/AODSIM/final_phase1_PU50bx25_DES19_62_V8-v1/00000/009500C4-A821-E411-8E47-02163E00E7E0.root'
        # /PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEM2019Upg14DR-final_phase1_age1k_PU140bx25_PH1_1K_FB_V2-v1/AODSIM
        #'/store/mc/GEM2019Upg14DR/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/AODSIM/final_phase1_age1k_PU140bx25_PH1_1K_FB_V2-v1/00000/00746114-E31F-E411-B4BC-002618FDA211.root'
        '/store/mc/GEM2019Upg14DR/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/AODSIM/final_phase1_age1k_PU140bx25_PH1_1K_FB_V2-v1/00000/007B34CC-331F-E411-BD67-002618943918.root'
        # /RelValTTbar_14TeV/CMSSW_6_2_0_SLHC11-DES23_62_V1_UPG2023Muon-v1/GEN-SIM-RECO
        #'/store/relval/CMSSW_6_2_0_SLHC11/RelValTTbar_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/2AEE7860-67C6-E311-AEB4-0025905964B4.root'
    )
)

if options.runOnData:
    process.source.fileNames = [
        # /JetHT/CMSSW_6_2_1-GR_R_62_V1_dvmc_RelVal_jetHT2012Cdvmc-v1/RECO
        '/store/relval/CMSSW_6_2_1/JetHT/RECO/GR_R_62_V1_dvmc_RelVal_jetHT2012Cdvmc-v1/00000/0030D71C-1E38-E311-BA7A-0025905964BE.root'
    ]

if options.runOnData :
    if options.runSubJets :
        options.outFilename += '_data_subjets.root'
    else :
        options.outFilename += '_data.root'
else :
    if options.runSubJets :
        options.outFilename += '_mc_subjets.root'
    else :
        options.outFilename += '_mc.root'

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outFilename)
)

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    allowUnscheduled = cms.untracked.bool(True)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag + '::All'

###FormultipleCSVv2
##############################################
# Get calibrations for the CSVV2 tagger: Phase I, II and retrained 72X
##############################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
process.CondDBSetup,
timetype = cms.string('runnumber'),
toGet = cms.VPSet(cms.PSet(
record = cms.string('BTauGenericMVAJetTagComputerRcd'),
                tag = cms.string('MVAJetTags_620SLHCX')
)),
connect = cms.string('sqlite_fip:RecoBTag/PerformanceMeasurements/data/MVAJetTags_620SLHCX_Phase1And2Upgrade_v10.db'),
BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
)
process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord")
###FormultipleCSVv2:end

process.load("Configuration.StandardSequences.MagneticField_cff")
#----------------------------------------------------------------------
# 2017 geometry:
#process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2017_cff')
# 2019 geometry:
process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2019_cff')
# 2023 geometry:
#process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaper_cff')
#----------------------------------------------------------------------
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackHistory.TrackHistory_cff")
process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")


#-------------------------------------
## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFilename),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

#-------------------------------------
## PAT Configuration
postfix = "PFlow"
jetAlgo="AK5"

from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
          jetCorrections=jetCorrectionsAK4, pvCollection=cms.InputTag(pvCollection))

## Top projections in PF2PAT
getattr(process,"pfPileUpJME"+postfix).checkClosestZVertex = False
getattr(process,"pfNoPileUpJME"+postfix).enable = options.usePFchs
if options.useTTbarFilter:
    getattr(process,"pfNoMuonJME"+postfix).enable = False
    getattr(process,"pfNoElectronJME"+postfix).enable = False
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoJet"+postfix).enable = False
else:
    getattr(process,"pfNoMuonJME"+postfix).enable = False
    getattr(process,"pfNoElectronJME"+postfix).enable = False
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoJet"+postfix).enable = False

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag('pfJets'+postfix),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag('ak5GenJetsNoNu'+postfix),
    postfix = postfix
)

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
#-------------------------------------

#-------------------------------------
# CA8 jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNu = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"+postfix)
)

from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJets = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(options.fatJetPtMin)
)

## CA8 pruned jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ca8GenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"+postfix),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ca8PFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetPtMin)
)

if options.runSubJets:
    ## PATify CA8 jets
    addJetCollection(
        process,
        labelName = 'CA8',
        jetSource = cms.InputTag('ca8PFJets'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK7,
        postfix = postfix
    )
    process.patJetGenJetMatchPatJetsCA8PFlow.matched = cms.InputTag('ca8GenJetsNoNu')
    process.patJetPartonAssociationPatJetsCA8PFlow.partons = cms.InputTag("patJetPartonsPFlow")

    addJetCollection(
        process,
        labelName = 'CA8Pruned',
        jetSource = cms.InputTag('ca8PFJetsPruned'),
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK7,
        postfix = postfix
    )
    process.patJetGenJetMatchPatJetsCA8PrunedPFlow.matched = cms.InputTag('ca8GenJetsNoNu')
    process.patJetPartonAssociationPatJetsCA8PrunedPFlow.partons = cms.InputTag("patJetPartonsPFlow")

    addJetCollection(
        process,
        labelName = 'CA8PrunedSubJets',
        jetSource = cms.InputTag('ca8PFJetsPruned','SubJets'),
        btagInfos=bTagInfos,
        btagDiscriminators=bTagDiscriminatorsSubJets,
        jetCorrections=jetCorrectionsAK5,
        postfix = postfix
    )
    process.patJetGenJetMatchPatJetsCA8PrunedSubJetsPFlow.matched = cms.InputTag('ca8GenJetsNoNuPruned','SubJets')
    process.patJetPartonAssociationPatJetsCA8PrunedSubJetsPFlow.partons = cms.InputTag("patJetPartonsPFlow")

## Establish references between PAT fat jets and PAT subjets using the BoostedJetMerger
process.selectedPatJetsCA8PrunedPFlowPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8PrunedPFlow"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8PrunedSubJetsPFlow")
)
#-------------------------------------

#-------------------------------------
if options.runOnData and options.runSubJets:
    ## Remove MC matching when running over data
    removeMCMatching( process, ['All'] )

## Add TagInfos to PAT jets
patJets = ['patJets'+postfix]
if options.runSubJets:
    patJets += ['patJetsCA8PFlow','patJetsCA8PrunedSubJetsPFlow']

for m in patJets:
    if hasattr(process,m):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
#-------------------------------------

#-------------------------------------
## Produce a collection of good primary vertices
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')

## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

#-------------------------------------

#-------------------------------------
## For MC PU jet ID
## Select GenJets for MC PU jet ID
process.ak5GenJetsForMCPUJetID = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('ak5GenJetsNoNu'+postfix),
    cut = cms.string('pt > 8.0'),
    filter = cms.bool(False) # in case no GenJets pass the selection, do not filter events, just produce an empty GenJet collection
)
## Adjust dR cut and GenJet collection name
getattr(process,'patJetGenJetMatch'+postfix).maxDeltaR = cms.double(0.25)
getattr(process,'patJetGenJetMatch'+postfix).matched = cms.InputTag("ak5GenJetsForMCPUJetID")
#-------------------------------------

#-------------------------------------
## Load modules for primary vertex sorting
process.load("RecoVertex.PrimaryVertexSorter.sortedOfflinePrimaryVertices_cff")

process.sortedGoodOfflinePrimaryVertices = process.sortedOfflinePrimaryVertices.clone(
    src = cms.InputTag('goodOfflinePrimaryVertices')
)

#-------------------------------------

#-------------------------------------
if options.useTTbarFilter:
    process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionFilter_cfi")
    process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionProducer_cfi")

    process.ttbarselectionproducer.isData       = options.runOnData
    process.ttbarselectionproducer.electronColl = cms.InputTag('selectedPatElectrons'+postfix)
    process.ttbarselectionproducer.muonColl     = cms.InputTag('selectedPatMuons'+postfix)
    process.ttbarselectionproducer.jetColl      = cms.InputTag('selectedPatJets'+postfix)
    process.ttbarselectionproducer.metColl      = cms.InputTag('patMETs'+postfix)
    process.ttbarselectionfilter.select_ee   = True
    process.ttbarselectionfilter.select_mumu = True
    process.ttbarselectionfilter.select_emu  = True
    process.ttbarselectionfilter.Keep_all_events  = False

    ## Change the cone size of muon isolation to 0.3
    getattr(process,"pfIsolatedMuons"+postfix).isolationValueMapsCharged = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03'+postfix ) )
    getattr(process,"pfIsolatedMuons"+postfix).isolationValueMapsNeutral = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03'+postfix ), cms.InputTag( 'muPFIsoValueGamma03'+postfix ) )
    getattr(process,"pfIsolatedMuons"+postfix).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03'+postfix )
    getattr(process,"pfIsolatedMuons"+postfix).combinedIsolationCut = cms.double(9999.)
    getattr(process,"pfIsolatedMuons"+postfix).isolationCut = cms.double(9999.)

    getattr(process,"patMuons"+postfix).isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag('muPFIsoValueNeutral03'+postfix),
        pfPhotons = cms.InputTag('muPFIsoValueGamma03'+postfix),
        pfChargedHadrons = cms.InputTag('muPFIsoValueCharged03'+postfix),
        pfChargedAll = cms.InputTag('muPFIsoValueChargedAll03'+postfix),
        pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU03'+postfix)
    )

    ## Change the cone size of electron isolation to 0.3
    getattr(process,'pfElectrons'+postfix).isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFId'+postfix))
    getattr(process,'pfElectrons'+postfix).deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFId'+postfix)
    getattr(process,'pfElectrons'+postfix).isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFId'+postfix), cms.InputTag('elPFIsoValueGamma03PFId'+postfix))

    getattr(process,'pfIsolatedElectrons'+postfix).isolationValueMapsCharged = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFId'+postfix))
    getattr(process,'pfIsolatedElectrons'+postfix).deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFId'+postfix)
    getattr(process,'pfIsolatedElectrons'+postfix).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFId'+postfix), cms.InputTag('elPFIsoValueGamma03PFId'+postfix))
    getattr(process,'pfIsolatedElectrons'+postfix).combinedIsolationCut = cms.double(9999.)
    getattr(process,'pfIsolatedElectrons'+postfix).isolationCut = cms.double(9999.)

    ## Electron ID
    process.load("EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi")
    process.eidMVASequence = cms.Sequence( process.mvaTrigV0 + process.mvaNonTrigV0 )

    getattr(process,'patElectrons'+postfix).electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    getattr(process,'patElectrons'+postfix).electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    getattr(process,'patElectrons'+postfix).isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag('elPFIsoValueCharged03PFId'+postfix),
        pfChargedAll = cms.InputTag('elPFIsoValueChargedAll03PFId'+postfix),
        pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFId'+postfix),
        pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral03PFId'+postfix),
        pfPhotons = cms.InputTag('elPFIsoValueGamma03PFId'+postfix)
    )

    ## Conversion rejection
    ## This should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
    #setattr(process,'patConversions'+postfix) = cms.EDProducer("PATConversionProducer",
        #electronSource = cms.InputTag('selectedPatElectrons'+postfix)
    #)
#-------------------------------------

#-------------------------------------
## Change the minimum number of tracker hits used in the track selection
if options.changeMinNumberOfHits:
    for m in process.producerNames().split(' '):
        if m.startswith('impactParameterTagInfos'):
            print "Changing 'minimumNumberOfHits' for " + m + " to " + str(options.minNumberOfHits)
            getattr(process, m).minimumNumberOfHits = cms.int32(options.minNumberOfHits)
#-------------------------------------

#-------------------------------------
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag(pvCollection))

## Need to update impactParameterTagInfosForPVSorting to remove circular module dependency
getattr(process,'impactParameterTagInfosForPVSorting').primaryVertex = cms.InputTag("goodOfflinePrimaryVertices")
#-------------------------------------

#-------------------------------------
if not options.runOnData:
    ## JP calibrations:
    process.GlobalTag.toGet = cms.VPSet(
      cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
           tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
      cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
           tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )
#---------------------------------------

###FormultipleCSVv2
#Forcing some selection configurations used during the training
process.inclusiveSecondaryVertexFinderTagInfosPFlowSelection = process.inclusiveSecondaryVertexFinderTagInfosPFlow.clone()
#Uncomment to change the selection on the qualityClass
#process.inclusiveSecondaryVertexFinderTagInfosPFlowSelection.trackSelection.qualityClass = cms.string('any')
## Add additional b-tag discriminators
getattr(process,'patJets'+postfix).discriminatorSources += cms.VInputTag(
cms.InputTag("combinedSecondaryVertexIVFV2Phase1MidPUBJetTags"),
cms.InputTag("combinedSecondaryVertexIVFV2NewBJetTags"),
cms.InputTag("combinedSecondaryVertexIVFV2Phase1HighPUBJetTags"),
cms.InputTag("combinedSecondaryVertexIVFV2Phase2HighPUBJetTags"),
cms.InputTag("combinedSecondaryVertexIVFV2Phase2HighPUMaxEta4BJetTags")
)
#CSVIVFV2: MLP-based
process.combinedSecondaryVertexIVFV2=process.combinedSecondaryVertexV2.clone(
calibrationRecords = cms.vstring(
'CombinedSVIVFV2RecoVertex_53x',
'CombinedSVIVFV2PseudoVertex_53x',
'CombinedSVIVFV2NoVertex_53x'
)
)
getattr(process,'combinedSecondaryVertexV2BJetTags'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexIVFV2')
#CSVIVFV2: 50Pu
process.combinedSecondaryVertexIVFV2Phase1MidPU=process.combinedSecondaryVertexV2.clone(
calibrationRecords = cms.vstring(
'CombinedSVIVFV2Phase1_50PU_RecoVertex',
'CombinedSVIVFV2Phase1_50PU_PseudoVertex',
'CombinedSVIVFV2Phase1_50PU_NoVertex'
)
)
process.combinedSecondaryVertexIVFV2Phase1MidPUBJetTags = process.combinedSecondaryVertexV2BJetTagsPFlow.clone(
jetTagComputer = cms.string('combinedSecondaryVertexIVFV2Phase1MidPU'),
tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPFlow"), cms.InputTag("inclusiveSecondaryVertexFinderTagInfosPFlowSelection"))
)
###CSVIVFV2: 140PU +aging
process.combinedSecondaryVertexIVFV2Phase1HighPU=process.combinedSecondaryVertexV2.clone(
calibrationRecords = cms.vstring(
'CombinedSVIVFV2Phase1_140PU_RecoVertex',
'CombinedSVIVFV2Phase1_140PU_PseudoVertex',
'CombinedSVIVFV2Phase1_140PU_NoVertex'
)
)
process.combinedSecondaryVertexIVFV2Phase1HighPUBJetTags = process.combinedSecondaryVertexV2BJetTagsPFlow.clone(
jetTagComputer = cms.string('combinedSecondaryVertexIVFV2Phase1HighPU'),
tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPFlow"), cms.InputTag("inclusiveSecondaryVertexFinderTagInfosPFlowSelection"))
)
###CSVIVFV2: PHASE2 140PU +aging
process.combinedSecondaryVertexIVFV2Phase2HighPU=process.combinedSecondaryVertexV2.clone(
calibrationRecords = cms.vstring(
'CombinedSVIVFV2Phase2JEC_140PU_RecoVertex',
'CombinedSVIVFV2Phase2JEC_140PU_PseudoVertex',
'CombinedSVIVFV2Phase2JEC_140PU_NoVertex'
)
)
process.combinedSecondaryVertexIVFV2Phase2HighPUBJetTags = process.combinedSecondaryVertexV2BJetTagsPFlow.clone(
jetTagComputer = cms.string('combinedSecondaryVertexIVFV2Phase2HighPU'),
tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPFlow"), cms.InputTag("inclusiveSecondaryVertexFinderTagInfosPFlowSelection"))
)
###CSVIVFV2: PHASE2 140PU +aging MaxEta4  
process.combinedSecondaryVertexIVFV2Phase2HighPUMaxEta4=process.combinedSecondaryVertexV2.clone(
calibrationRecords = cms.vstring(                                                                                                                                                 
'CombinedSVIVFV2Phase2JEC_140PU_MaxEta4_RecoVertex',
'CombinedSVIVFV2Phase2JEC_140PU_MaxEta4_PseudoVertex',
'CombinedSVIVFV2Phase2JEC_140PU_MaxEta4_NoVertex'
)                                                                                                                                                                                 
)                                                                                                                                                                           
process.combinedSecondaryVertexIVFV2Phase2HighPUMaxEta4BJetTags = process.combinedSecondaryVertexV2BJetTagsPFlow.clone(                                                          
jetTagComputer = cms.string('combinedSecondaryVertexIVFV2Phase2HighPUMaxEta4'),                                                                                     
tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPFlow"), cms.InputTag("inclusiveSecondaryVertexFinderTagInfosPFlowSelection"))                                      
)
###CSVIVFV2: new csv from Shimaa
process.combinedSecondaryVertexIVFV2New=process.combinedSecondaryVertexV2.clone(
calibrationRecords = cms.vstring(
'CombinedSVV2MVA_RecoVertex',
'CombinedSVV2MVA_PseudoVertex',
'CombinedSVV2MVA_NoVertex'
)
)
process.combinedSecondaryVertexIVFV2NewBJetTags = process.combinedSecondaryVertexV2BJetTagsPFlow.clone(
jetTagComputer = cms.string('combinedSecondaryVertexIVFV2New'),
tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPFlow"), cms.InputTag("inclusiveSecondaryVertexFinderTagInfosPFlowSelection"))
)
###FormultipleCSVv2:end

#-------------------------------------
process.load("RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff")
# The following combinations should be considered:
# For b-tagging performance measurements:
#   process.btagana.useSelectedTracks    = True
#   process.btagana.useTrackHistory      = False (or True for Mistag systematics with GEN-SIM-RECODEBUG samples)
#   process.btagana.produceJetTrackTree  = False
#   process.btagana.produceAllTrackTree  = False
#   process.btagana.producePtRelTemplate = False (or True for PtRel fit studies)
# or data/MC validation of jets, tracks and SVs:
#   process.btagana.useSelectedTracks    = False (or True for JP calibration)
#   process.btagana.useTrackHistory      = False
#   process.btagana.produceJetTrackTree  = True
#   process.btagana.produceAllTrackTree  = False
#   process.btagana.producePtRelTemplate = False
# or general tracks, PV and jet performance studies:
#   process.btagana.useSelectedTracks    = True
#   process.btagana.useTrackHistory      = False
#   process.btagana.produceJetTrackTree  = False
#   process.btagana.produceAllTrackTree  = True
#   process.btagana.producePtRelTemplate = False
#------------------
process.btagana.useSelectedTracks     = True  ## False if you want to run on all tracks : for commissioning studies
process.btagana.useTrackHistory       = False ## Can only be used with GEN-SIM-RECODEBUG files
process.btagana.produceJetTrackTree   = False ## True if you want to keep info for tracks associated to jets : for commissioning studies
process.btagana.produceAllTrackTree   = False ## True if you want to keep info for tracks associated to jets : for commissioning studies
process.btagana.producePtRelTemplate  = options.producePtRelTemplate  ## True for performance studies
#------------------
process.btagana.primaryVertexColl     = cms.InputTag(pvCollection)
process.btagana.Jets                  = cms.InputTag('selectedPatJets'+postfix)
process.btagana.patMuonCollectionName = cms.InputTag('selectedPatMuons')
process.btagana.use_ttbar_filter      = cms.bool(options.useTTbarFilter)
process.btagana.MaxEta                = cms.double(4.2) ## for extended forward pixel coverage
process.btagana.MinPt                 = cms.double(20.0)
#process.btagana.triggerTable         = cms.InputTag('TriggerResults::HLT') # Data and MC
#Configure here which sv taginfos to store in the SV_xxx branches
process.btagana.svTagInfos = cms.string('inclusiveSecondaryVertexFinder')
process.btagana.svComputer = cms.InputTag('combinedSecondaryVertexIVFV2')

if options.runSubJets:
    process.btaganaSubJets = process.btagana.clone(
        produceJetTrackTree = cms.bool(True),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJetsCA8PrunedSubJetsPFlow'),
        FatJets             = cms.InputTag('selectedPatJetsCA8PFlow'),
        PrunedFatJets       = cms.InputTag('selectedPatJetsCA8PrunedPFlowPacked'),
        runSubJets          = options.runSubJets,
        use_ttbar_filter    = cms.bool(False)
    )
#---------------------------------------

#---------------------------------------
## Trigger selection !
#import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
#process.JetHLTFilter = hlt.triggerResultsFilter.clone(
#    triggerConditions = cms.vstring(
#        "HLT_PFJet80_v*"
#    ),
#    hltResults = cms.InputTag("TriggerResults","","HLT"),
#    l1tResults = cms.InputTag( "" ),
#    throw = cms.bool( False ) #set to false to deal with missing triggers while running over different trigger menus
#)
#---------------------------------------

#---------------------------------------
## Optional MET filters:
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
#process.load("RecoMET.METFilters.metFilters_cff")
#process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
#---------------------------------------

#---------------------------------------
## Filter for HCAL laser events in prompt 2012A+B+C, snippet for "Datasets from the 2013 rereco and Multijet parked":
## https://twiki.cern.ch/twiki/bin/view/CMS/PdmVKnowFeatures#HCAL_laser_events_in_prompt_2012
#process.load("EventFilter.HcalRawToDigi.hcallaserFilterFromTriggerResult_cff")  # ---------->  Does not seem to work with CMSSW_6_2_X rereco !!!
#---------------------------------------

#---------------------------------------
## Event counter
from MyAnalysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone()
process.selectedEvents = eventCounter.clone()
#---------------------------------------

#---------------------------------------
## Switching from AK5 to AK4 jets
print 'Switching from AK5 to AK4 jets'
print '******************************'

## Change cone sizes (note that only the cose sizes are changed, but not the collection names)
process.ak5GenJetsNoNuPFlow.rParam = cms.double(0.4)
process.pfJetsPFlow.rParam = cms.double(0.4)
process.jetTracksAssociatorAtVertexPFlow.coneSize = cms.double(0.4)

## Select JEC version
#JEC from https://hypernews.cern.ch/HyperNews/CMS/get/upgrade-tp-studies/154.html
#jec='PhaseII_Shashlik140PU_V1'
#JEC from https://hypernews.cern.ch/HyperNews/CMS/get/upgrade-tp-studies/157.html
#Note fix with respect to the one circulated in the HN above, which had a bug
#jec='PhaseII_Shashlik140PU_V2'
#JEC from https://hypernews.cern.ch/HyperNews/CMS/get/upgrade-tp-studies/160.html
#jec='PhaseI_140PU_V2'
#JEC from https://hypernews.cern.ch/HyperNews/CMS/get/upgrade-tp-studies/157/1.html
jec='PhaseI_50PU_V2'

jechera=jec

if jec == 'PhaseII_Shashlik140PU_V1':
    jechera='PhaseII_Shashlik140PU'

print 'Using jec', jec ,' with hera name ', jechera

## Get AK4PFchs JECs from a sqlite file
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
    ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
    cms.PSet(
         record = cms.string('JetCorrectionsRecord'),
         tag    = cms.string('JetCorrectorParametersCollection_' + jechera + '_AK4PFchs'),
         label  = cms.untracked.string('AK4PFchs')
    ),
    ## here you add as many jet types as you need
    ## note that the tag name is specific for the particular sqlite file
    ),
    connect = cms.string('sqlite_fip:RecoBTag/PerformanceMeasurements/data/' + jec + '.db')
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
#---------------------------------------

if options.useBetterPVSorting:
    from CommonTools.RecoAlgos.sortedPrimaryVertices_cfi import *
    process.betterOfflinePrimaryVertices=sortedPrimaryVertices.clone(jets = "ak5CaloJets")
    process.vertSeq = cms.Sequence(
        process.betterOfflinePrimaryVertices
        )
    process.primaryVertexFilter.vertexCollection= cms.InputTag('betterOfflinePrimaryVertices')
    process.goodOfflinePrimaryVertices.src = cms.InputTag('betterOfflinePrimaryVertices')

if options.useBetterPVSortingPF:
    from CommonTools.RecoAlgos.sortedPFPrimaryVertices_cfi import *
    process.betterOfflinePrimaryVertices=sortedPFPrimaryVertices.clone(jets = "ak5PFJets")
    process.vertSeq = cms.Sequence(
        process.betterOfflinePrimaryVertices
        )
    process.primaryVertexFilter.vertexCollection= cms.InputTag('betterOfflinePrimaryVertices')
    process.goodOfflinePrimaryVertices.src = cms.InputTag('betterOfflinePrimaryVertices')

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.JetHLTFilter*
    process.noscraping
    *process.primaryVertexFilter
    * process.goodOfflinePrimaryVertices
    #* process.HBHENoiseFilter
    #* process.CSCTightHaloFilter
    #* process.EcalDeadCellTriggerPrimitiveFilter
    #* process.eeBadScFilter
    #* process.ecalLaserCorrFilter
    #* process.trackingFailureFilter
    #* process.trkPOGFilters
)
#if options.runOnData:
    #process.filtSeq *= process.hcalfilter # ---------->  Does not seem to work with CMSSW_6_2_X rereco !!!


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK5Jets:
    process.analyzerSeq += process.btagana
if options.runSubJets:
    process.analyzerSeq += process.btaganaSubJets
if options.processStdAK5Jets and options.useTTbarFilter:
    process.analyzerSeq.replace( process.btagana, process.ttbarselectionproducer * process.ttbarselectionfilter * process.btagana )
#---------------------------------------

if options.useBetterPVSorting or options.useBetterPVSortingPF:
    process.p = cms.Path(
        process.allEvents
        *process.vertSeq
        * process.filtSeq
        * process.selectedEvents
        * process.analyzerSeq
        )
else:
    process.p = cms.Path(
        process.allEvents
        * process.filtSeq
        * process.selectedEvents
        * process.analyzerSeq
        )

# Delete predefined output module (needed for running with CRAB)
del process.out

open('pydump.py','w').write(process.dumpPython())

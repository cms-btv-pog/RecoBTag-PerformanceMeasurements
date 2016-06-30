
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
options.register('reportEvery', 10,
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
options.register('usePuppi', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use Puppi"
)
options.register('usePuppiForBTagging', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use Puppi candidates for b tagging"
)
options.register('mcGlobalTag', '80X_mcRun2_asymptotic_v4',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', '80X_dataRun2_v4', 
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)
options.register('runFatJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run fat jets"
)
options.register('runSubJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)
options.register('processStdAK4Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Process standard AK4 jets"
)
options.register('producePtRelTemplate', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Produce PtRel template"
)
options.register('fatJetRawPtMin', 150.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum raw pT for fat jets (default is 150 GeV)"
)
options.register('fatJetPtMin', 200.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum pT for fat jets (default is 200 GeV)"
)
options.register('fatJetAbsEtaMax', 2.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum |eta| for fat jets (default is 2.5)"
)
options.register('useTTbarFilter', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use TTbar filter"
)
options.register('miniAOD', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Running on miniAOD"
)
options.register('fastSim', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Running using FastSim"
)
options.register('useExplicitJTA', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use explicit jet-track association"
)
options.register('jetAlgo', 'AntiKt',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Jet clustering algorithms (default is AntiKt)"
)
options.register('fatJetRadius', 0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Distance parameter R for fat jet clustering (default is 0.8)"
)
options.register('useLegacyTaggers', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use legacy taggers"
)
options.register('useSoftDrop', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SoftDrop jets"
)
options.register('usePruned', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use pruned jets"
)
#Generally leave to False unless you know what you are doing
options.register('runIVF', False, 
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run IVF, currently leave to False!"
)
### Master switch for boosted b tag commissioning: overrider several other switches
options.register('doBoostedCommissioning', False, 
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Make NTuples with branches for boosted b tag commissioning: overrider several other switches"
)
## Do Ctag
options.register('doCTag', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Make NTuples with branches for CTag"
)
options.register('usePrivateJEC',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'Use JECs from private SQLite files')
options.register('jecDBFile',
    'Spring16_25nsV3', 
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'SQLite filename for JECs')
options.register('isReHLT',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    '80X reHLT samples')

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', -1)

options.parseArguments()

## Use either PFchs or Puppi
if options.usePFchs and options.usePuppi:
    print "WARNING: Both usePFchs and usePuppi set to True. Giving priority to Puppi."
    options.usePFchs = False

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Running using FastSim samples: %s"%('True' if options.fastSim else 'False')
print "Running on MiniAOD: %s"%('True' if options.miniAOD else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')
print "Using Puppi: %s"%('True' if options.usePuppi else 'False')
print "Using Puppi for b tagging: %s"%('True' if (options.usePuppi and options.usePuppiForBTagging) else 'False')

## Subjets only stored when also running over fat jets
if options.runSubJets and not options.runFatJets:
    print "WARNING: You are attempting to store subjet information without running over fat jets. Please enable running over fat jets in order to store the subjet information."
    options.runSubJets = False

if options.doBoostedCommissioning:
    print "**********NTuples will be made for boosted b tag commissioning. The following switches will be reset:**********"
    options.processStdAK4Jets=False
    print "Option processStdAK4Jets will be set to '",options.processStdAK4Jets,"'"
    options.runFatJets=True  
    options.runSubJets = True
    print "Option runFatJets will be set to '",options.runFatJets,"'"
    print "Option runSubJets  will be set to '",options.runSubJets,"'"
    print "********************"
if options.doCTag:
    print "**********You are making NTuple for CTag*************" 

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if not options.usePFchs:
    jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK8 = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK4[1].append('L2L3Residual')
    jetCorrectionsAK8[1].append('L2L3Residual')

trigresults='TriggerResults::HLT'
if options.runOnData: options.isReHLT=False
if options.isReHLT: trigresults = trigresults+'2'

## b-tag infos
bTagInfosLegacy = [
    'impactParameterTagInfos'
   ,'secondaryVertexTagInfos'
   ,'inclusiveSecondaryVertexFinderTagInfos'
   ,'secondaryVertexNegativeTagInfos'
   ,'inclusiveSecondaryVertexFinderNegativeTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
]
bTagInfos = [
    'pfImpactParameterTagInfos'
   ,'pfSecondaryVertexTagInfos'
   ,'pfInclusiveSecondaryVertexFinderTagInfos'
   ,'pfSecondaryVertexNegativeTagInfos'
   ,'pfInclusiveSecondaryVertexFinderNegativeTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
   ,'pfInclusiveSecondaryVertexFinderCvsLTagInfos'
   ,'pfInclusiveSecondaryVertexFinderNegativeCvsLTagInfos'
]
## b-tag discriminators
bTagDiscriminatorsLegacy = [
    'jetBProbabilityBJetTags'
   ,'jetProbabilityBJetTags'
   ,'positiveOnlyJetBProbabilityBJetTags'
   ,'positiveOnlyJetProbabilityBJetTags'
   ,'negativeOnlyJetBProbabilityBJetTags'
   ,'negativeOnlyJetProbabilityBJetTags'
   ,'trackCountingHighPurBJetTags'
   ,'trackCountingHighEffBJetTags'
   ,'negativeTrackCountingHighEffBJetTags'
   ,'negativeTrackCountingHighPurBJetTags'
   ,'simpleSecondaryVertexHighEffBJetTags'
   ,'simpleSecondaryVertexHighPurBJetTags'
   ,'negativeSimpleSecondaryVertexHighEffBJetTags'
   ,'negativeSimpleSecondaryVertexHighPurBJetTags'
   ,'combinedSecondaryVertexV2BJetTags'
   ,'positiveCombinedSecondaryVertexV2BJetTags'
   ,'negativeCombinedSecondaryVertexV2BJetTags'
   ,'combinedInclusiveSecondaryVertexV2BJetTags'
   ,'positiveCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'negativeCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'softPFMuonBJetTags'
   ,'positiveSoftPFMuonBJetTags'
   ,'negativeSoftPFMuonBJetTags'
   ,'softPFElectronBJetTags'
   ,'positiveSoftPFElectronBJetTags'
   ,'negativeSoftPFElectronBJetTags'
   ,'combinedMVAv2BJetTags'
   ,'negativeCombinedMVAv2BJetTags'
   ,'positiveCombinedMVAv2BJetTags'
]
bTagDiscriminators = [
    'pfJetBProbabilityBJetTags'
   ,'pfJetProbabilityBJetTags'
   ,'pfPositiveOnlyJetBProbabilityBJetTags'
   ,'pfPositiveOnlyJetProbabilityBJetTags'
   ,'pfNegativeOnlyJetBProbabilityBJetTags'
   ,'pfNegativeOnlyJetProbabilityBJetTags'
   ,'pfTrackCountingHighPurBJetTags'
   ,'pfTrackCountingHighEffBJetTags'
   ,'pfNegativeTrackCountingHighPurBJetTags'
   ,'pfNegativeTrackCountingHighEffBJetTags'
   ,'pfSimpleSecondaryVertexHighEffBJetTags'
   ,'pfSimpleSecondaryVertexHighPurBJetTags'
   ,'pfNegativeSimpleSecondaryVertexHighEffBJetTags'
   ,'pfNegativeSimpleSecondaryVertexHighPurBJetTags'
   ,'pfCombinedSecondaryVertexV2BJetTags'
   ,'pfPositiveCombinedSecondaryVertexV2BJetTags'
   ,'pfNegativeCombinedSecondaryVertexV2BJetTags'
   ,'pfCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'softPFMuonBJetTags'
   ,'positiveSoftPFMuonBJetTags'
   ,'negativeSoftPFMuonBJetTags'
   ,'softPFElectronBJetTags'
   ,'positiveSoftPFElectronBJetTags'
   ,'negativeSoftPFElectronBJetTags'
   ,'pfCombinedMVAV2BJetTags'
   ,'pfNegativeCombinedMVAV2BJetTags'
   ,'pfPositiveCombinedMVAV2BJetTags'
   ,'pfCombinedCvsBJetTags'
   ,'pfNegativeCombinedCvsBJetTags'
   ,'pfPositiveCombinedCvsBJetTags'
   ,'pfCombinedCvsLJetTags'
   ,'pfNegativeCombinedCvsLJetTags'
   ,'pfPositiveCombinedCvsLJetTags'
]

## Legacy taggers not supported with MiniAOD
if options.miniAOD and options.useLegacyTaggers:
    print "WARNING: Legacy taggers not supported with MiniAOD"
    options.useLegacyTaggers = False

## If using legacy taggers
if options.useLegacyTaggers:
    bTagInfos = bTagInfosLegacy
    bTagDiscriminators = bTagDiscriminatorsLegacy

## Clustering algorithm label
algoLabel = 'CA'
if options.jetAlgo == 'AntiKt':
    algoLabel = 'AK'

## Postfix
postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
jetSource = 'pfJetsPFBRECO'+postfix
genJetCollection = 'ak4GenJetsNoNu'
pfCandidates = 'particleFlow'
pvSource = 'offlinePrimaryVertices'
svSource = 'inclusiveCandidateSecondaryVertices'
muSource = 'muons'
elSource = 'gedGsfElectrons'
patMuons = 'selectedPatMuons'
trackSource = 'generalTracks'
## If running on MiniAOD
if options.miniAOD:
    genParticles = 'prunedGenParticles'
    jetSource = 'ak4PFJets'
    genJetCollection = 'ak4GenJetsNoNu'
    pfCandidates = 'packedPFCandidates'
    pvSource = 'offlineSlimmedPrimaryVertices'
    svSource = 'slimmedSecondaryVertices'
    trackSource = 'unpackedTracksAndVertices'
    muSource = 'slimmedMuons'
    elSource = 'slimmedElectrons'
    patMuons = 'slimmedMuons'

process = cms.Process("BTagAna")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

if options.miniAOD:
    process.source.fileNames = [
        #'/store/relval/CMSSW_8_0_0/RelValTTbar_13/MINIAODSIM/PU25ns_80X_mcRun2_asymptotic_v4-v1/10000/A65CD249-BFDA-E511-813A-0025905A6066.root'
        '/store/mc/RunIISpring16MiniAODv2/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/02A0E7CE-1B35-E611-8612-0CC47A7FC4C8.root'
    ]
    if options.runOnData:
        process.source.fileNames = [
            '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/275/125/00000/86B9123D-7F36-E611-9D13-02163E0139C7.root'
        ]
    if options.fastSim:
        process.source.fileNames = [
            '/store/relval/CMSSW_8_0_0/RelValTTbar_13/MINIAODSIM/PU25ns_80X_mcRun2_asymptotic_v4_FastSim-v2/10000/8E75D08A-3FDE-E511-8374-0CC47A4C8F26.root'
        ]
else:
    process.source.fileNames = [
        '/store/relval/CMSSW_8_0_0/RelValProdTTbar_13/AODSIM/80X_mcRun2_asymptotic_v4-v1/10000/DE81ABBF-1DDA-E511-8AF8-0026189438B5.root'
    ]
    if options.runOnData:
        process.source.fileNames = [
            '/store/data/Run2016B/SingleMuon/AOD/PromptReco-v2/000/275/125/00000/DA2EC189-7E36-E611-8C63-02163E01343B.root'
        ]
    if options.fastSim:
        process.source.fileNames = [
            '/store/relval/CMSSW_8_0_0/RelValTTbar_13/GEN-SIM-DIGI-RECO/PU25ns_80X_mcRun2_asymptotic_v4_FastSim-v2/10000/0400D094-63DD-E511-8B51-0CC47A4C8ED8.root'
        ]

## Define the output file name
if options.runOnData :
    options.outFilename += '_data'
else :
    options.outFilename += '_mc'

if options.runFatJets :
    options.outFilename += '_FatJets'

if options.runSubJets :
    options.outFilename += '_Subjets'

if options.fastSim :
    options.outFilename += '_FastSim'

if options.doBoostedCommissioning:
  options.outFilename += '_BoostedCommissioning' 

options.outFilename += '.root'

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

#Set GT by hand:
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = globalTag
#Choose automatically:
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_' + ('data' if options.runOnData else 'mc'))

#Loading calibrations from db file, example of code for any future use
#process.load("CondCore.DBCommon.CondDBSetup_cfi")
#process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
#    process.CondDBSetup,
#    timetype = cms.string('runnumber'),
#    toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string('BTauGenericMVAJetTagComputerRcd'),
#            tag = cms.string('MVAJetTags')
#        )
#    ),
#    connect = cms.string('sqlite_fip:RecoBTag/PerformanceMeasurements/data/MVAJetTags.db'),
#    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
#)
#process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord")

if options.usePrivateJEC:
    
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    dbfile=''
    if options.runOnData: dbfile=options.jecDBFile+'_DATA'
    else: dbfile=options.jecDBFile+'_MC'
    print "\nUsing private SQLite file", dbfile, "\n"
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
		    connect = cms.string( "sqlite_file:"+dbfile+'.db'),
		    toGet =  cms.VPSet(
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+dbfile+"_AK4PF"),
				    label= cms.untracked.string("AK4PF")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+dbfile+"_AK4PFchs"),
				    label= cms.untracked.string("AK4PFchs")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+dbfile+"_AK4PFPuppi"),
				    label= cms.untracked.string("AK4PFPuppi")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+dbfile+"_AK8PF"),
				    label= cms.untracked.string("AK8PF")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+dbfile+"_AK8PFchs"),
				    label= cms.untracked.string("AK8PFchs")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+dbfile+"_AK8PFPuppi"),
				    label= cms.untracked.string("AK8PFPuppi")
				    ),
			    )
		    )

    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

### to activate the new JP calibration: using the data base
trkProbaCalibTag = "JPcalib_MC800_v1" ###DM"JPcalib_MC76X_25ns_v1"
if options.runOnData:
  trkProbaCalibTag = "JPcalib_Data80X_2016B_v1" ###DM"JPcalib_Data76X_2015D_v1"
# process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
      tag = cms.string(trkProbaCalibTag),
      connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
    )
)

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimTracker.TrackHistory.TrackHistory_cff")
process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
#process.load("RecoBTag.Configuration.RecoBTag_cff")

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
if not options.miniAOD:
    ## PAT Configuration
    jetAlgo="AK4"

    from PhysicsTools.PatAlgos.tools.pfTools import *
    usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
              jetCorrections=jetCorrectionsAK4, pvCollection=cms.InputTag(pvSource))

    ## Top projections in PF2PAT
    getattr(process,"pfPileUpJME"+postfix).checkClosestZVertex = False
    getattr(process,"pfNoPileUpJME"+postfix).enable = options.usePFchs
    getattr(process,"pfNoMuonJMEPFBRECO"+postfix).enable = False
    getattr(process,"pfNoElectronJMEPFBRECO"+postfix).enable = False

    if options.usePuppi:
        process.load('CommonTools.PileupAlgos.Puppi_cff')
        from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
        process.ak4PFJets = ak4PFJets.clone(src = 'puppi', doAreaFastjet = True)
        jetSource = 'ak4PFJets'
        if options.usePuppiForBTagging: pfCandidates = 'puppi'
else:
    ## GenJets
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
    process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

    ## PFchs selection
    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

    ## Reco jets
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    if options.usePFchs:
        process.ak4PFJets = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
    elif options.usePuppi:
        process.load('CommonTools.PileupAlgos.Puppi_cff')
        process.puppi.candName           = cms.InputTag(pfCandidates)
        process.puppi.vertexName         = cms.InputTag(pvSource)
        process.puppi.useExistingWeights = cms.bool(True)
        process.puppi.clonePackedCands   = cms.bool(True)
        process.ak4PFJets = ak4PFJets.clone(src = 'puppi', doAreaFastjet = True)
        if options.usePuppiForBTagging: pfCandidates = 'puppi'
    else:
        process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag(genJetCollection),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = options.useExplicitJTA,
    postfix = postfix
)

#-------------------------------------

#-------------------------------------
## Fat jets (Gen and Reco)
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.genJetsNoNu = ak4GenJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.fatJetRadius),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix))
)
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.PFJetsCHS = ak4PFJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.fatJetRadius),
    src = getattr(process,jetSource).src,
    srcPVs = getattr(process,jetSource).srcPVs,
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(options.fatJetRawPtMin)
)
## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPruned = ak4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.fatJetRadius),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix)),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.PFJetsCHSPruned = ak4PFJetsPruned.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.fatJetRadius),
    src = getattr(process,jetSource).src,
    srcPVs = getattr(process,jetSource).srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetRawPtMin)
)
## SoftDrop fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
process.genJetsNoNuSoftDrop = ak4GenJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.fatJetRadius),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix)),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0 = cms.double(options.fatJetRadius),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsSoftDrop_cfi import ak4PFJetsSoftDrop
process.PFJetsCHSSoftDrop = ak4PFJetsSoftDrop.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.fatJetRadius),
    R0 = cms.double(options.fatJetRadius),
    src = getattr(process,jetSource).src,
    srcPVs = getattr(process,jetSource).srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetRawPtMin)
)

if options.runFatJets:
    ## PATify the above jets
    addJetCollection(
        process,
        labelName='PFCHS',
        jetSource=cms.InputTag('PFJetsCHS'),
        algo=algoLabel,           # needed for jet flavor clustering
        rParam=options.fatJetRadius, # needed for jet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = (bTagDiscriminators + ([] if options.useLegacyTaggers else ['pfBoostedDoubleSecondaryVertex' + ('CA15' if algoLabel=='CA' else 'AK8') + 'BJetTags'])),
        jetCorrections = jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = options.useExplicitJTA,
        runIVF = options.runIVF,
        postfix = postfix
    )
    getattr(process,'selectedPatJetsPFCHS'+postfix).cut = cms.string("pt > %f && abs(eta) < %f"%(float(options.fatJetPtMin), float(options.fatJetAbsEtaMax)))
    addJetCollection(
        process,
        labelName='SoftDropPFCHS',
        jetSource=cms.InputTag('PFJetsCHSSoftDrop'),
        algo=algoLabel,
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False, # jet flavor disabled
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName='SoftDropSubjetsPFCHS',
        jetSource=cms.InputTag('PFJetsCHSSoftDrop','SubJets'),
        algo=algoLabel,           # needed for subjet flavor clustering
        rParam=options.fatJetRadius, # needed for subjet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag('genJetsNoNuSoftDrop','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets = cms.InputTag('PFJetsCHS'),              # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('PFJetsCHSSoftDrop'), # needed for subjet flavor clustering
        runIVF = options.runIVF,
        postfix = postfix
    )

    ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
    process.selectedPatJetsSoftDropPFCHSPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsSoftDropPFCHS"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsSoftDropSubjetsPFCHS"+postfix)
    )

    addJetCollection(
        process,
        labelName='PrunedPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned'),
        algo=algoLabel,
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False, # jet flavor disabled
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName='PrunedSubjetsPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned','SubJets'),
        algo=algoLabel,           # needed for subjet flavor clustering
        rParam=options.fatJetRadius, # needed for subjet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag('genJetsNoNuPruned','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets = cms.InputTag('PFJetsCHS'),              # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('PFJetsCHSPruned'), # needed for subjet flavor clustering
        runIVF = options.runIVF,
        postfix = postfix
    )

    ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
    process.selectedPatJetsPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsPrunedPFCHS"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPFCHS"+postfix)
    )

    ## Pack fat jets with subjets
    process.packedPatJetsPFCHS = cms.EDProducer("JetSubstructurePacker",
            jetSrc = cms.InputTag('selectedPatJetsPFCHS'+postfix),
            distMax = cms.double(options.fatJetRadius),
            algoTags = cms.VInputTag(),
            algoLabels = cms.vstring(),
            fixDaughters = cms.bool(False)
    )
    if options.useSoftDrop:
        process.packedPatJetsPFCHS.algoTags.append( cms.InputTag('selectedPatJetsSoftDropPFCHSPacked') )
        process.packedPatJetsPFCHS.algoLabels.append( 'SoftDrop' )
    if options.usePruned:
        process.packedPatJetsPFCHS.algoTags.append( cms.InputTag('selectedPatJetsPrunedPFCHSPacked') )
        process.packedPatJetsPFCHS.algoLabels.append( 'Pruned' )

    #-------------------------------------
    ## N-subjettiness
    from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

    process.Njettiness = Njettiness.clone(
        src = cms.InputTag("PFJetsCHS"),
        R0  = cms.double(options.fatJetRadius)
    )

    getattr(process,'patJetsPFCHS'+postfix).userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']


    #-------------------------------------
    ## Grooming ValueMaps
    process.SoftDrop = cms.EDProducer("RecoJetDeltaRValueMapProducer",
        src = cms.InputTag("PFJetsCHS"),
        matched = cms.InputTag("selectedPatJetsSoftDropPFCHSPacked"),
        distMax = cms.double(options.fatJetRadius),
        values = cms.vstring('mass','pt','eta','phi','jecFactor(0)'),
        valueLabels = cms.vstring('Mass','Pt','Eta','Phi','jecFactor0'),
        lazyParser = cms.bool(True)
    )
    process.Pruned = cms.EDProducer("RecoJetDeltaRValueMapProducer",
        src = cms.InputTag("PFJetsCHS"),
        matched = cms.InputTag("selectedPatJetsPrunedPFCHSPacked"),
        distMax = cms.double(options.fatJetRadius),
        values = cms.vstring('mass','pt','eta','phi','jecFactor(0)'),
        valueLabels = cms.vstring('Mass','Pt','Eta','Phi','jecFactor0'),
        lazyParser = cms.bool(True)
    )

    getattr(process,'patJetsPFCHS'+postfix).userData.userFloats.src += ['SoftDrop:Mass','SoftDrop:Pt','SoftDrop:Eta','SoftDrop:Phi','SoftDrop:jecFactor0',
                                                                        'Pruned:Mass'  ,'Pruned:Pt'  ,'Pruned:Eta'  ,'Pruned:Phi'  ,'Pruned:jecFactor0']

#-------------------------------------

#Switch to hadron flavour priority for all jet collections
if not options.runOnData:
    process.patJetFlavourAssociation.hadronFlavourHasPriority = cms.bool(True)
    process.patJetFlavourAssociationPFlow.hadronFlavourHasPriority = cms.bool(True)
    if options.runFatJets:
        process.patJetFlavourAssociationPFCHSPFlow.hadronFlavourHasPriority = cms.bool(True)
    if options.runSubJets:
        process.patJetFlavourAssociationPrunedSubjetsPFCHSPFlow.hadronFlavourHasPriority = cms.bool(True)
        process.patJetFlavourAssociationSoftDropSubjetsPFCHSPFlow.hadronFlavourHasPriority = cms.bool(True)

#-------------------------------------
if options.runOnData:
    # Remove MC matching when running over data
    from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
    removeMCMatching( process, ['Photons', 'Electrons','Muons', 'Taus', 'Jets', 'METs', 'PFElectrons','PFMuons', 'PFTaus'] )

#-------------------------------------
## Add GenParticlePruner for boosted b-tagging studies
if not options.runOnData:
    process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
                                                     src = cms.InputTag(genParticles),
                                                     select = cms.vstring("drop  *  ", #by default
                                                                          "keep ( status = 3 || (status>=21 && status<=29) )", #keep hard process particles
                                                                          "keep abs(pdgId) = 13 || abs(pdgId) = 15" #keep muons and taus
                                                                          )
                                                     )

#-------------------------------------

## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag(pvSource),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)
#-------------------------------------

#-------------------------------------
if options.useTTbarFilter:
    process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionFilter_cfi")
    process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionProducer_cfi")

    #electron id
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
        
    if options.miniAOD:
        process.ttbarselectionproducer.electronColl = cms.InputTag('slimmedElectrons')
        process.ttbarselectionproducer.muonColl     = cms.InputTag('slimmedMuons')
        process.ttbarselectionproducer.jetColl      = cms.InputTag('selectedPatJets'+postfix)
        process.ttbarselectionproducer.metColl      = cms.InputTag('slimmedMETs')
        switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    else:
        process.ttbarselectionproducer.electronColl = cms.InputTag('selectedPatElectrons'+postfix)
        process.ttbarselectionproducer.muonColl     = cms.InputTag('selectedPatMuons'+postfix)
        process.ttbarselectionproducer.jetColl      = cms.InputTag('selectedPatJets'+postfix)
        process.ttbarselectionproducer.metColl      = cms.InputTag('patMETs'+postfix)
        switchOnVIDElectronIdProducer(process, DataFormat.AOD)

    # Set up electron ID (VID framework)
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    switchOnVIDElectronIdProducer(process, dataFormat=DataFormat.MiniAOD)
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



    #process.ttbarselectionproducer.isData       = options.runOnData
    #process.ttbarselectionproducer.electronColl = cms.InputTag('selectedPatElectrons'+postfix)
    #process.ttbarselectionproducer.muonColl     = cms.InputTag('selectedPatMuons'+postfix)
    #process.ttbarselectionproducer.jetColl      = cms.InputTag('selectedPatJets'+postfix)
    #process.ttbarselectionproducer.metColl      = cms.InputTag('patMETs'+postfix)
    #process.ttbarselectionfilter.select_ee   = True
    #process.ttbarselectionfilter.select_mumu = True
    #process.ttbarselectionfilter.select_emu  = True
    #process.ttbarselectionfilter.Keep_all_events  = False

    ## Change the cone size of muon isolation to 0.3
    #getattr(process,"pfIsolatedMuons"+postfix).isolationValueMapsCharged = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03'+postfix ) )
    #getattr(process,"pfIsolatedMuons"+postfix).isolationValueMapsNeutral = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03'+postfix ), cms.InputTag( 'muPFIsoValueGamma03'+postfix ) )
    #getattr(process,"pfIsolatedMuons"+postfix).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03'+postfix )
    #getattr(process,"pfIsolatedMuons"+postfix).combinedIsolationCut = cms.double(9999.)
    #getattr(process,"pfIsolatedMuons"+postfix).isolationCut = cms.double(9999.)

    #getattr(process,"patMuons"+postfix).isolationValues = cms.PSet(
    #    pfNeutralHadrons = cms.InputTag('muPFIsoValueNeutral03'+postfix),
    #    pfPhotons = cms.InputTag('muPFIsoValueGamma03'+postfix),
    #    pfChargedHadrons = cms.InputTag('muPFIsoValueCharged03'+postfix),
    #    pfChargedAll = cms.InputTag('muPFIsoValueChargedAll03'+postfix),
    #    pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU03'+postfix)
    #)

    ## Change the cone size of electron isolation to 0.3
    #getattr(process,'pfElectrons'+postfix).isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFId'+postfix))
    #getattr(process,'pfElectrons'+postfix).deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFId'+postfix)
    #getattr(process,'pfElectrons'+postfix).isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFId'+postfix), cms.InputTag('elPFIsoValueGamma03PFId'+postfix))

    #getattr(process,'pfIsolatedElectrons'+postfix).isolationValueMapsCharged = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFId'+postfix))
    #getattr(process,'pfIsolatedElectrons'+postfix).deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFId'+postfix)
    #getattr(process,'pfIsolatedElectrons'+postfix).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFId'+postfix), cms.InputTag('elPFIsoValueGamma03PFId'+postfix))
    #getattr(process,'pfIsolatedElectrons'+postfix).combinedIsolationCut = cms.double(9999.)
    #getattr(process,'pfIsolatedElectrons'+postfix).isolationCut = cms.double(9999.)

    ## Electron ID
    #process.load("EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi")
    #process.eidMVASequence = cms.Sequence( process.mvaTrigV0 + process.mvaNonTrigV0 )

    #getattr(process,'patElectrons'+postfix).electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    #getattr(process,'patElectrons'+postfix).electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    #getattr(process,'patElectrons'+postfix).isolationValues = cms.PSet(
    #    pfChargedHadrons = cms.InputTag('elPFIsoValueCharged03PFId'+postfix),
    #    pfChargedAll = cms.InputTag('elPFIsoValueChargedAll03PFId'+postfix),
    #    pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFId'+postfix),
    #    pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral03PFId'+postfix),
    #    pfPhotons = cms.InputTag('elPFIsoValueGamma03PFId'+postfix)
    #)

    ## Conversion rejection
    ## This should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
    #setattr(process,'patConversions'+postfix) = cms.EDProducer("PATConversionProducer",
        #electronSource = cms.InputTag('selectedPatElectrons'+postfix)
    #)
#-------------------------------------

#-------------------------------------
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag(pvSource))

#-------------------------------------
## Add TagInfos to PAT jets
for m in ['patJets'+postfix, 'patJetsPFCHS'+postfix, 'patJetsSoftDropSubjetsPFCHS'+postfix, 'patJetsPrunedSubjetsPFCHS'+postfix]:
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )

#-------------------------------------
## Adapt fat jet b tagging
if options.runFatJets:
    getattr(process,'softPFElectronsTagInfosPFCHS'+postfix).DeltaRElectronJet = cms.double(options.fatJetRadius) # default is 0.4
    if options.useLegacyTaggers:
        # Set the cone size for the jet-track association to the jet radius
        getattr(process,'jetTracksAssociatorAtVertexPFCHS'+postfix).coneSize = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'secondaryVertexTagInfosPFCHS'+postfix).trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius)   # default is 0.3
        getattr(process,'secondaryVertexTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
        getattr(process,'inclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'inclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).extSVDeltaRToJet = cms.double(options.fatJetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
        process.jetProbabilityComputerFat = process.jetProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.3
        getattr(process,'jetProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('jetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
        process.jetBProbabilityComputerFat = process.jetBProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.4
        getattr(process,'jetBProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('jetBProbabilityComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
        process.combinedSecondaryVertexV2ComputerFat = process.combinedSecondaryVertexV2Computer.clone()
        process.combinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        process.combinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        getattr(process,'combinedInclusiveSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexV2ComputerFat')
    else:
        # Set the cone size for the jet-track association to the jet radius
        getattr(process,'pfImpactParameterTagInfosPFCHS'+postfix).maxDeltaR = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'pfSecondaryVertexTagInfosPFCHS'+postfix).trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius)   # default is 0.3
        getattr(process,'pfSecondaryVertexTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).extSVDeltaRToJet = cms.double(options.fatJetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
        process.candidateJetProbabilityComputerFat = process.candidateJetProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.3
        getattr(process,'pfJetProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateJetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
        process.candidateJetBProbabilityComputerFat = process.candidateJetBProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.4
        getattr(process,'pfJetBProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateJetBProbabilityComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
        process.candidateCombinedSecondaryVertexV2ComputerFat = process.candidateCombinedSecondaryVertexV2Computer.clone()
        process.candidateCombinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        process.candidateCombinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        getattr(process,'pfCombinedInclusiveSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateCombinedSecondaryVertexV2ComputerFat')

if options.miniAOD:
  process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#-------------------------------------
from RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff import *
process.btagana = bTagAnalyzer.clone()
if options.useLegacyTaggers:
    process.btagana = bTagAnalyzerLegacy.clone()
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
process.btagana.tracksColl            = cms.InputTag(trackSource) 
process.btagana.useSelectedTracks     = True  ## False if you want to run on all tracks : for commissioning studies
process.btagana.useTrackHistory       = False ## Can only be used with GEN-SIM-RECODEBUG files
process.btagana.fillsvTagInfo         = False ## True if you want to store information relative to the svTagInfos, set to False if produceJetTrackTree is set to False
process.btagana.produceJetTrackTree   = False ## True if you want to keep info for tracks associated to jets : for commissioning studies
process.btagana.produceAllTrackTree   = False ## True if you want to keep info for all tracks : for commissioning studies
process.btagana.producePtRelTemplate  = options.producePtRelTemplate  ## True for performance studies
#------------------
process.btagana.storeTagVariables     = False  ## True if you want to keep TagInfo TaggingVariables
process.btagana.storeCSVTagVariables  = True   ## True if you want to keep CSV TaggingVariables
process.btagana.primaryVertexColl     = cms.InputTag(pvSource)
process.btagana.Jets                  = cms.InputTag('selectedPatJets'+postfix)
process.btagana.muonCollectionName    = cms.InputTag(muSource)
process.btagana.patMuonCollectionName = cms.InputTag(patMuons)
process.btagana.use_ttbar_filter      = cms.bool(options.useTTbarFilter)
#process.btagana.triggerTable          = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.triggerTable          = cms.InputTag(trigresults) # Data and MC
process.btagana.genParticles          = cms.InputTag(genParticles)
process.btagana.candidates            = cms.InputTag(pfCandidates)

if options.doCTag:
    process.btagana.storeCTagVariables = True
    process.btagana.storeEventInfo = True
    process.btagana.doCTag = options.doCTag

## fillsvTagInfo set to False independently from the choices above, if produceJetTrackTree is set to False
if not process.btagana.produceJetTrackTree:
    process.btagana.fillsvTagInfo = False

if options.runFatJets:
    process.btaganaFatJets = process.btagana.clone(
        storeEventInfo      = cms.bool(not options.processStdAK4Jets),
        fillQuarks = cms.bool(True),
        allowJetSkipping    = cms.bool(False),
        storeTagVariables   = cms.bool(False),
        storeCSVTagVariables = cms.bool(True),
        storeTagVariablesSubJets = cms.bool(False),
        storeCSVTagVariablesSubJets = cms.bool(False),
        useSelectedTracks   = cms.bool(True),
        maxDeltaR           = cms.double(options.fatJetRadius),
        R0                  = cms.double(options.fatJetRadius),
        maxSVDeltaRToJet    = cms.double(options.fatJetRadius-(0.1+(options.fatJetRadius-0.8)*(0.4/0.7))), # linear interpolation from 0.7 at R=0.8 to 1.0 at R=1.5
        weightFile          = cms.FileInPath('RecoBTag/PerformanceMeasurements/data/BoostedDoubleSV_' + ('CA15' if algoLabel=='CA' else 'AK8') + '_BDT_v3.weights.xml.gz'),
        doubleSVBJetTags    = cms.string('pfBoostedDoubleSecondaryVertex' + ('CA15' if algoLabel=='CA' else 'AK8') + 'BJetTags'),
        distJetAxis         = cms.double(9999.),
        decayLength         = cms.double(9999.),
        deltaR              = cms.double(0.8),
        BranchNamePrefix    = cms.string('FatJetInfo'),
        Jets                = cms.InputTag('packedPatJetsPFCHS'),
        SubJets             = cms.VInputTag(),
        SubJetLabels        = cms.vstring(),
        runFatJets          = cms.bool(True),
        runSubJets          = options.runSubJets,
        svComputer          = cms.string('combinedSecondaryVertexV2ComputerFat' if options.useLegacyTaggers else 'candidateCombinedSecondaryVertexV2ComputerFat'),
        use_ttbar_filter    = cms.bool(False)
    )
    if options.useSoftDrop:
        process.btaganaFatJets.SubJets.append( cms.InputTag('selectedPatJetsSoftDropPFCHSPacked','SubJets') )
        process.btaganaFatJets.SubJetLabels.append( 'SoftDrop' )
    if options.usePruned:
        process.btaganaFatJets.SubJets.append( cms.InputTag('selectedPatJetsPrunedPFCHSPacked','SubJets') )
        process.btaganaFatJets.SubJetLabels.append( 'Pruned' )

if options.doBoostedCommissioning:
    process.btaganaFatJets.produceJetTrackTree  = True 
    process.btaganaFatJets.fillsvTagInfo = True  
    process.btaganaFatJets.storeCSVTagVariables = True  
    process.btaganaFatJets.storeCSVTagVariablesSubJets = True 
    print "**********NTuples will be made for boosted b tag commissioning. The following switches will be reset:**********"
    print "produceJetTrackTree set to '",process.btaganaFatJets.produceJetTrackTree,"'" 
    print "fillsvTagInfo set to '",process.btaganaFatJets.fillsvTagInfo,"'" 
    print "For fat jets: storeCSVTagVariables set to '",process.btaganaFatJets.storeCSVTagVariables,"'"
    print "For subjets:  storeCSVTagVariablesSubJet set to '",process.btaganaFatJets.storeCSVTagVariablesSubJets,"'"
    print "********************"

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
## Event counter
from RecoBTag.PerformanceMeasurements.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone()
process.selectedEvents = eventCounter.clone()
#---------------------------------------

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.JetHLTFilter*
    #process.noscraping
    process.primaryVertexFilter
)


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK4Jets:
    process.analyzerSeq += process.btagana
if options.runFatJets:
    process.analyzerSeq += process.btaganaFatJets
if options.processStdAK4Jets and options.useTTbarFilter:
    process.analyzerSeq.replace( process.btagana, process.ttbarselectionproducer * process.ttbarselectionfilter * process.btagana )
#---------------------------------------

process.p = cms.Path(
    process.allEvents
    * process.filtSeq
    * process.selectedEvents
    * process.analyzerSeq
)

# Delete predefined output module (needed for running with CRAB)
del process.out

open('pydump.py','w').write(process.dumpPython())

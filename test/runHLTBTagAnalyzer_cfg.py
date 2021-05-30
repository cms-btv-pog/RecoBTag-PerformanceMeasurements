from FWCore.ParameterSet.VarParsing import VarParsing

###############################
####### Parameters ############
###############################

options = VarParsing ('analysis')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('outFilename', 'JetTreeHLT',
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
options.register('dumpPython', None,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'Dump python config, pass SaveName.py'
)
options.register('mcGlobalTag', 'FIXME',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag, no default value provided"
)
options.register('dataGlobalTag', 'FIXME',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag, no default value provided"
)
options.register('runEventInfo', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run Event Info"
)
options.register('processStdAK4Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Process standard AK4 jets"
)
options.register('useTrackHistory', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "uses track history, for GEN-SIM-RECODEBUG samples only"
)
options.register('produceJetTrackTree', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "True if you want to run info for tracks associated to jets : for commissioning studies"
)
options.register('produceAllTrackTree', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Produce all track tree"
)

options.register('fillPU', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Fill PU"
)
### Options for upgrade studies
# Change hits requirements
options.register('changeMinNumberOfHits', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Change minimum number of tracker hits"
)
options.register('minNumberOfHits', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Minimum number of tracker hits"
)
# Change eta for extended forward pixel coverage
options.register('maxJetEta', 4.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum jet |eta| (default is 4.5)"
)
options.register('minJetPt', 20.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum jet pt (default is 20)"
)
options.register('usePrivateJEC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'Use JECs from private SQLite files')
options.register('jecDBFileMC', 'FIXME',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'SQLite filename for JECs, no default value provided')
options.register('jecDBFileData', 'FIXME',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'SQLite filename for JECs, no default value provided')
options.register('isReHLT', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    '80X reHLT samples')
options.register('JPCalibration', 'FIXME',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'JP Calibration pyload to use')
options.register('runJetVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Jet Variables')
options.register('runCaloJetVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Jet Variables')
options.register('runPuppiJetVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Jet Variables')
options.register('runTagVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Tag Variables')
options.register('runQuarkVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run c/b quark Variables')
options.register('runHadronVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Hadron Variables')
options.register('runGenVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Gen Variables')
options.register('runCSVTagVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run CSV TaggingVariables')
options.register('runCSVTagTrackVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run CSV Tagging Track Variables')
options.register('runDeepFlavourTagVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run DeepFlavour TaggingVariables')
options.register('runPFElectronVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run PF Electron Variables')
options.register('runPFMuonVariables', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run PF Muon Variables')
# options.register('runPatMuons', False,
#     VarParsing.multiplicity.singleton,
#     VarParsing.varType.bool,
#     'True if you want to run Pat Muon Variables')
options.register('defaults', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'baseline default settings to be used')
options.register('eras', [],
    VarParsing.multiplicity.list,
    VarParsing.varType.string,
    'era modifiers to be used to be used')
options.register('groups', [],
    VarParsing.multiplicity.list,
    VarParsing.varType.string,
    'variable groups to be stored')
options.register(
    'skipEvents', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip N events"
)

options.register('reco', 'HLT_GRun',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'keyword to define HLT reconstruction')

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', -1)

options.parseArguments()
if options.defaults:
	from importlib import import_module
	try:
		defaults = import_module('RecoBTag.PerformanceMeasurements.defaults_HLT.%s' % options.defaults)
	except ImportError:
		raise ValueError('The default settings named %s.py are not present in PerformanceMeasurements/python/defaults_HLT/' % options.defaults)
	if not hasattr(defaults, 'common') or not isinstance(defaults.common, dict):
		raise RuntimeError('the default file %s.py does not contain a dictionary named common' % options.defaults)
	items = defaults.common.items()
	if hasattr(defaults, 'data') and options.runOnData:
		if not isinstance(defaults.data, dict):
			raise RuntimeError('the default file %s.py contains an object called "data" which is not a dictionary' % options.defaults)
		items.extend(defaults.data.items())
	if hasattr(defaults, 'mc') and not options.runOnData:
		if not isinstance(defaults.mc, dict):
			raise RuntimeError('the default file %s.py contains an object called "mc" which is not a dictionary' % options.defaults)
		items.extend(defaults.mc.items())
	for key, value in items:
		if key not in options._beenSet:
			raise ValueError('The key set by the defaults: %s does not exist among the cfg options!' % key)
		elif not options._beenSet[key]:
			if key == 'inputFiles' and options.inputFiles: continue #skip input files that for some reason are never considered set
			print 'setting default option for', key
			setattr(options, key, value)

from RecoBTag.PerformanceMeasurements.HLTBTagAnalyzer_cff import *
btagana_tmp = HLTBTagAnalyzer.clone()
print('Storing the variables from the following groups:')
options_to_change = set() #store which swtiches we need on
for requiredGroup in options.groups:
  print(requiredGroup)
  found=False
  for existingGroup in btagana_tmp.groups:
    if(requiredGroup==existingGroup.group):
      existingGroup.store=True
      for var in existingGroup.variables:
        if "CaloJet." in var or "PuppiJet." in var:
            var = var.split(".")[1]
        options_to_change.update([i for i in variableDict[var].runOptions])
      found=True
      break
  if(not found):
    print('WARNING: The group ' + requiredGroup + ' was not found')

#change values accordingly
for switch in options_to_change:
  if switch in ["runTagVariablesSubJets","runCSVTagVariablesSubJets"]:
      continue
  elif switch not in options._beenSet:
    raise ValueError('The option set by the variables: %s does not exist among the cfg options!' % switch)
  elif not options._beenSet[switch]:
    print 'Turning on %s, as some stored variables demands it' % switch
    setattr(options, switch, True)

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

trigresults='TriggerResults::HLT'
if options.runOnData: options.isReHLT=False
if options.isReHLT: trigresults = trigresults+'2'


#pfjets = "hltAK4PFJets" #original ak4PFJetsCHS
#calojets = "hltAK4CaloJets" #original ak4CaloJets
#puppijets = "hltAK4PFPuppiJets"
#PFDeepCSVTags = "hltDeepCombinedSecondaryVertexBPFPatJetTags" # original: pfDeepCSVJetTags
PFDeepFlavourTags = "hltPFDeepFlavourJetTags" # original: pfDeepFlavourJetTagsSlimmedDeepFlavour
rho = "hltFixedGridRhoFastjetAll" #original fixedGridRhoFastjetAll
#hltVertices = "hltVerticesPFFilter" #original offlinePrimaryVertices
#hltVerticesSlimmed = "hltVerticesPFFilter" #original offlineSlimmedPrimaryVertices
#siPixelClusters = "hltSiPixelClusters" #original siPixelClusters
#ecalRecHit = "hltEcalRecHit" #original ecalRecHit
#hbhereco = "hltHbhereco" #original hbhereco
#hfreco = "hltHfreco" #original hfreco
#horeco = "hltHoreco" #original horeco
#rpcRecHits = "hltRpcRecHits" #original rpcRecHits
#tracks = "hltMergedTracks" #original generalTracks
#payload = "AK4PFHLT" #original AK4PFchs




## Postfix
# postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
# jetSource = 'pfJetsPFBRECO'+postfix
# patJetSource = 'selectedPatJets'+postfix
# patJetSource = 'hltSlimmedJets'
patJetSource = 'hltPatJets'
patCaloJetSource = 'hltPatJetsCalo'
patPuppiJetSource = 'hltPatJetsPuppi'
genJetCollection = 'ak4GenJetsNoNu'
# pfCandidates = 'particleFlow'
pfCandidates = 'hltParticleFlow'
# pvSource = 'offlinePrimaryVertices'
pvSource = "hltVerticesPFFilter" 
# svSource = 'inclusiveCandidateSecondaryVertices'
#svSource = 'hltDeepInclusiveMergedVerticesPF'
# muSource = 'muons'
muSource = 'hltMuons'
# elSource = 'gedGsfElectrons'
elSource = 'hltEgammaGsfElectrons'
# patMuons = 'selectedPatMuons'
# trackSource = 'generalTracks'
trackSource = "hltMergedTracks"


if options.reco == 'HLT_BTagROI':
        patJetSource = 'hltPatJetsROI'
        pfCandidates = 'hltParticleFlowForBTag'
        pvSource = "hltVerticesPFFilterForBTag"
        trackSource = "hltMergedTracksForBTag"
        PFDeepFlavourTags = "hltPFDeepFlavourJetTagsROI"
        rho = "hltFixedGridRhoFastjetAllForBTag" #original fixedGridRhoFastjetAll
        patPuppiJetSource = 'hltPatJetsPuppiROI'

def customisePFForPixelTracks(process):

    process.hltPFMuonMerging.TrackProducers = cms.VInputTag("hltIterL3MuonTracks", "hltPixelTracksClean")
    process.hltPFMuonMerging.selectedTrackQuals = cms.VInputTag("hltIterL3MuonTracks", "hltPixelTracksClean")

    return process


###
### HLT configuration
###
if options.reco == 'HLT_GRun':
        #from JMETriggerAnalysis.Common.configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
        from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process

elif options.reco == 'HLT_BTagROI':

        #from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump_johnda import cms, process
        from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process

        from RecoBTag.PerformanceMeasurements.customise_hlt import *
        process = addPaths_PFJetsForBtag(process)


elif options.reco == 'HLT_Run3TRK':
        
    # (a) Run-3 tracking: standard
        from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
        #from JMETriggerAnalysis.Common.configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
        from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
        process = customizeHLTRun3Tracking(process)

elif options.reco == 'HLT_Run3TRKWithPU':
        # (b) Run-3 tracking: all pixel vertices
        from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
        #from JMETriggerAnalysis.Common.configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
        from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3TrackingAllPixelVertices
        process = customizeHLTRun3TrackingAllPixelVertices(process)

elif options.reco == 'HLT_PixelTracks':
        from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
        from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3TrackingAllPixelVertices
        process = customizeHLTRun3TrackingAllPixelVertices(process)

        process = customisePFForPixelTracks(process)

        trackSource = "hltPixelTracksClean" 
else:
        raise RuntimeError('keyword "reco = '+opts.reco+'" not recognised')





# remove cms.OutputModule objects from HLT config-dump
for _modname in process.outputModules_():
    _mod = getattr(process, _modname)
    if type(_mod) == cms.OutputModule:
       process.__delattr__(_modname)
       # if options.verbosity > 0:
       #    print '> removed cms.OutputModule:', _modname

# remove cms.EndPath objects from HLT config-dump
for _modname in process.endpaths_():
    _mod = getattr(process, _modname)
    if type(_mod) == cms.EndPath:
       process.__delattr__(_modname)
       # if options.verbosity > 0:
       #    print '> removed cms.EndPath:', _modname

# remove selected cms.Path objects from HLT config-dump
print '-'*108
print '{:<99} | {:<4} |'.format('cms.Path', 'keep')
print '-'*108
for _modname in sorted(process.paths_()):
    _keepPath = _modname.startswith('MC_') and ('Jets' in _modname or 'MET' in _modname or 'DeepCSV' in _modname or 'DeepFlavour' in _modname or 'AK8Calo' in _modname)
    # _keepPath = _modname.startswith('MC_') and ('Jets' in _modname or 'DeepCSV' in _modname or 'DeepFlavour' in _modname or 'AK8Calo' in _modname)
    # _keepPath = _modname.startswith('MC_')
#    _keepPath |= _modname.startswith('MC_ReducedIterativeTracking')
    if _keepPath:
      print '{:<99} | {:<4} |'.format(_modname, '+')
      continue
    _mod = getattr(process, _modname)
    if type(_mod) == cms.Path:
      process.__delattr__(_modname)
    #   print '{:<99} | {:<4} |'.format(_modname, '')
print '-'*108

# remove FastTimerService
if hasattr(process, 'FastTimerService'):
  del process.FastTimerService

# remove MessageLogger
if hasattr(process, 'MessageLogger'):
  del process.MessageLogger

###
### customizations
###
from JMETriggerAnalysis.Common.customise_hlt import *
process = addPaths_MC_JMECalo(process)
process = addPaths_MC_JMEPF(process)
process = addPaths_MC_JMEPFCluster(process)
process = addPaths_MC_JMEPFPuppi(process)
from RecoBTag.PerformanceMeasurements.PATLikeConfig import customizePFPatLikeJets
process = customizePFPatLikeJets(process)

if options.reco == 'HLT_BTagROI':

        from RecoBTag.PerformanceMeasurements.customise_hlt import *
        process = addPaths_MC_JMEPFPuppiROI(process)

        from RecoBTag.PerformanceMeasurements.ROIPATLikeConfig import customizePFPatLikeJetsROI
        process = customizePFPatLikeJetsROI(process)




# if not options.eras:
# 	process = cms.Process("BTagAna")
# else:
# 	from Configuration.StandardSequences.Eras import eras
# 	eras_to_use = []
# 	for era in options.eras:
# 		if hasattr(eras, era):
# 			eras_to_use.append(getattr(eras, era))
# 		else:
# 			raise ValueError('The requested era (%s) is not available' % era)
# 	process = cms.Process("BTagAna", *eras_to_use)

## ES modules for PF-Hadron Calibrations
import os
from CondCore.CondDB.CondDB_cfi import CondDB as _CondDB

process.pfhcESSource = cms.ESSource('PoolDBESSource',
  _CondDB.clone(connect = 'sqlite_file:'+os.environ['CMSSW_BASE']+'/src/JMETriggerAnalysis/NTuplizers/data/PFHC_Run3Winter20_HLT_v01.db'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('PFCalibrationRcd'),
      tag = cms.string('PFCalibration_HLT_mcRun3_2021'),
      label = cms.untracked.string('HLT'),
    ),
  ),
)

process.pfhcESPrefer = cms.ESPrefer('PoolDBESSource', 'pfhcESSource')
#process.hltParticleFlow.calibrationsLabel = '' # standard label for Offline-PFHC in GT
## ES modules for HLT JECs
process.jescESSource = cms.ESSource('PoolDBESSource',
  _CondDB.clone(connect = 'sqlite_file:'+os.environ['CMSSW_BASE']+'/src/JMETriggerAnalysis/NTuplizers/data/JESC_Run3Winter20_V1_MC.db'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4CaloHLT'),
      label = cms.untracked.string('AK4CaloHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFClusterHLT'),
      label = cms.untracked.string('AK4PFClusterHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFHLT'),
      label = cms.untracked.string('AK4PFHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFHLT'),
      label = cms.untracked.string('AK4PFchsHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFPuppiHLT'),
      label = cms.untracked.string('AK4PFPuppiHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4CaloHLT'),#!!
      label = cms.untracked.string('AK8CaloHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFClusterHLT'),#!!
      label = cms.untracked.string('AK8PFClusterHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFHLT'),#!!
      label = cms.untracked.string('AK8PFHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFHLT'),#!!
      label = cms.untracked.string('AK8PFchsHLT'),
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag = cms.string('JetCorrectorParametersCollection_Run3Winter20_V1_MC_AK4PFPuppiHLT'),#!!
      label = cms.untracked.string('AK8PFPuppiHLT'),
    ),
  ),
)
process.jescESPrefer = cms.ESPrefer('PoolDBESSource', 'jescESSource')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
# process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
# process.MessageLogger.cerr.default.limit = 10

process.MessageLogger.suppressWarning = cms.untracked.vstring(
        'hltPatJetFlavourAssociationCalo'
)
process.MessageLogger.suppressError = cms.untracked.vstring(
        'hltPatJetFlavourAssociationCalo'
)

## Input files
# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring()
# )
#
# process.source.fileNames = [
#     '/store/mc/PhaseIFall16DR/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/AODSIM/PhaseIFall16PUFlat20to50_81X_upgrade2017_realistic_v26-v1/50000/0039E945-35E3-E611-AF8D-001E675A6C2A.root'
# ]
# if options.runOnData:
#     process.source.fileNames = [
#         '/store/data/Run2016B/SingleMuon/AOD/PromptReco-v2/000/275/125/00000/DA2EC189-7E36-E611-8C63-02163E01343B.root'
#     ]
# if options.inputFiles:
#     process.source.fileNames = options.inputFiles

## Define the output file name
if options.runOnData :
    options.outFilename += '_data'
else :
    options.outFilename += '_mc'

options.outFilename += '.root'

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outFilename)
)

## Events to process
process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    # allowUnscheduled = cms.untracked.bool(True)
)

#Set GT by hand:
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag.globaltag = globalTag
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

### to activate the new JP calibration: using the data base
# trkProbaCalibTag = options.JPCalibration
# process.GlobalTag.toGet = cms.VPSet(
#     cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
#       tag = cms.string(trkProbaCalibTag),
#       connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
#     )
# )

# process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
# process.load("Configuration.Geometry.GeometryRecoDB_cff")
# process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#-------------------------------------
## Output Module Configuration (expects a path 'p')
# from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFilename),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    # outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

#-------------------------------------


# if options.runOnData:
#     # Remove MC matching when running over data
#     from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
#     removeMCMatching( process, ['Photons', 'Electrons','Muons', 'Taus', 'Jets', 'METs', 'PFElectrons','PFMuons', 'PFTaus'] )

#-------------------------------------
## Add GenParticlePruner for boosted b-tagging studies
if not options.runOnData:
    process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
        src = cms.InputTag(genParticles),
        select = cms.vstring(
            "drop  *  ", #by default
            "keep ( status = 3 || (status>=21 && status<=29) ) && pt > 0", #keep hard process particles with non-zero pT
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
# process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#     vertexCollection = cms.InputTag(pvSource),
#     minimumNDOF = cms.uint32(4) ,
#     maxAbsZ = cms.double(24),
#     maxd0 = cms.double(2)
# )
#-------------------------------------

#-------------------------------------
## Change the minimum number of tracker hits used in the track selection
if options.changeMinNumberOfHits:
    for m in process.producerNames().split(' '):
        if m.startswith('pfImpactParameterTagInfos'):
            print "Changing 'minimumNumberOfHits' for " + m + " to " + str(options.minNumberOfHits)
            getattr(process, m).minimumNumberOfHits = cms.int32(options.minNumberOfHits)

#-------------------------------------
process.btagana = HLTBTagAnalyzer.clone()
#------------------
#Handle groups
for requiredGroup in process.btagana.groups:
   for storedGroup in btagana_tmp.groups:
     if (requiredGroup.group == storedGroup.group):
       requiredGroup.store = storedGroup.store

process.btagana.MaxEta                = options.maxJetEta ## for extended forward pixel coverage
process.btagana.MinPt                 = options.minJetPt
process.btagana.tracksColl            = cms.InputTag(trackSource)
process.btagana.useTrackHistory       = options.useTrackHistory ## Can only be used with GEN-SIM-RECODEBUG files
process.btagana.produceJetTrackTruthTree = options.useTrackHistory ## can only be used with GEN-SIM-RECODEBUG files and when useTrackHistory is True
process.btagana.produceAllTrackTree   = options.produceAllTrackTree ## True if you want to run info for all tracks : for commissioning studies
#------------------
process.btagana.runTagVariables     = options.runTagVariables  ## True if you want to run TagInfo TaggingVariables
process.btagana.runCSVTagVariables  = options.runCSVTagVariables   ## True if you want to run CSV TaggingVariables
process.btagana.runCSVTagTrackVariables  = options.runCSVTagTrackVariables   ## True if you want to run CSV Tagging Track Variables
process.btagana.runDeepFlavourTagVariables = options.runDeepFlavourTagVariables
process.btagana.primaryVertexColl     = cms.InputTag(pvSource)
process.btagana.Jets                  = cms.InputTag(patJetSource)
process.btagana.CaloJets              = cms.InputTag(patCaloJetSource)
process.btagana.PuppiJets             = cms.InputTag(patPuppiJetSource)
process.btagana.muonCollectionName    = cms.InputTag(muSource)
process.btagana.electronCollectionName= cms.InputTag(elSource)
# process.btagana.patMuonCollectionName = cms.InputTag(patMuons)
process.btagana.rho                   = cms.InputTag(rho)
process.btagana.deepFlavourJetTags    = PFDeepFlavourTags
# process.btagana.triggerTable          = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.triggerTable          = cms.InputTag(trigresults) # Data and MC
process.btagana.genParticles          = cms.InputTag(genParticles)
process.btagana.candidates            = cms.InputTag(pfCandidates)
process.btagana.runJetVariables       = options.runJetVariables
process.btagana.runCaloJetVariables   = options.runCaloJetVariables
process.btagana.runPuppiJetVariables   = options.runPuppiJetVariables
process.btagana.runQuarkVariables     = options.runQuarkVariables
process.btagana.runHadronVariables    = options.runHadronVariables
process.btagana.runGenVariables       = options.runGenVariables
process.btagana.runPFElectronVariables = options.runPFElectronVariables
process.btagana.runPFMuonVariables = options.runPFMuonVariables
# process.btagana.runPatMuons = options.runPatMuons
process.btagana.runEventInfo = options.runEventInfo
process.btagana.runOnData = options.runOnData

if options.runOnData:
  process.btagana.runHadronVariables  = False
  process.btagana.runQuarkVariables   = False
  process.btagana.runGenVariables     = False

if not process.btagana.useTrackHistory  or not options.produceJetTrackTree:
    process.btagana.produceJetTrackTruthTree = False

if process.btagana.useTrackHistory:
    process.load('SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi')
    process.load('SimTracker.TrackerHitAssociation.tpClusterProducer_cfi')

if process.btagana.produceJetTrackTruthTree:
    process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
    process.load("SimTracker.TrackHistory.TrackHistory_cff")
    process.load("SimTracker.TrackHistory.TrackClassifier_cff")
    process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
    process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

#---------------------------------------
process.options.numberOfThreads = cms.untracked.uint32(1)
process.options.numberOfStreams = cms.untracked.uint32(1)
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
# process.filtSeq = cms.Sequence(
#     #process.JetHLTFilter*
#     #process.noscraping
#     process.primaryVertexFilter
# )


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
# if options.processStdAK4Jets:
process.analyzerSeq += process.btagana
#---------------------------------------

process.testOutput = cms.OutputModule("PoolOutputModule",
    # SelectEvents = cms.untracked.PSet(
    #     SelectEvents = cms.vstring(
    #         'HLT_EcalCalibration_v4',
    #         'HLT_HcalCalibration_v5'
    #     )
    # ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RAW'),
        filterName = cms.untracked.string('')
    ),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('testOutput.root'),
    outputCommands = cms.untracked.vstring(
        # 'drop *_hlt*_*_*',
        # 'keep *_hltEcalCalibrationRaw_*_*',
        # 'keep *_hltHcalCalibrationRaw_*_*',
        # 'keep edmTriggerResults_*_*_*',
        # 'keep triggerTriggerEvent_*_*_*'
        'drop *_*_*_*',
        # 'keep hltImpactParameterPatTagInfos_*_*_*',
        'keep *_hltImpactParameterPatTagInfos_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetCaloPatTagInfos_*_*',
        'keep *_hltInclusiveSecondaryVertexFinderPatTagInfos_*_*',
        'keep *_hltDeepBLifetimePFPatTagInfos_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetPatTagInfos_*_*',
        'keep *_hltDeepSecondaryVertexPFPatTagInfos_*_*',
    )
)
# process.myOutput = cms.EndPath(process.testOutput)

process.p = cms.Path(
    process.allEvents
    # * process.filtSeq
    * process.selectedEvents
    * process.analyzerSeq
    # * process.testOutput
)

# Delete predefined output module (needed for running with CRAB)
del process.out
# dump content of cms.Process to python file
if options.dumpPython is not None:
   open(options.dumpPython, 'w').write(process.dumpPython())

print ''
print 'option: output =', options.outFilename
# print 'option: reco =', options.reco
print 'option: dumpPython =', options.dumpPython
print ''
# print 'process.GlobalTag =', process.GlobalTag.dumpPython()
print 'process.GlobalTag =', process.GlobalTag.globaltag
print 'process.source =', process.source.dumpPython()
print 'process.maxEvents =', process.maxEvents.input
# print 'process.options =', process.options.dumpPython()
print '-------------------------------'

from FWCore.ParameterSet.VarParsing import VarParsing
import fnmatch
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
options.register('globalTag', 'FIXME',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag, no default value provided"
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
options.register('minJetPt', 25.0,
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
    "skip N events")

options.register('runTiming', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'run timing instead of rates')
options.register('numThreads', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'number of threads')
options.register('numStreams', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'number of streams')


options.register(
    'reco', 'HLT_GRun',
  VarParsing.multiplicity.singleton,
  VarParsing.varType.string,
  'keyword to define HLT reconstruction'
)

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
			print ('setting default option for', key)
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
    print ('Turning on %s, as some stored variables demands it' % switch)
    setattr(options, switch, True)

## Global tag
globalTag = options.globalTag
# if options.runOnData:
#     globalTag = options.dataGlobalTag

trigresults='TriggerResults::HLT'
if options.runOnData: options.isReHLT=False
if options.isReHLT: trigresults = trigresults+'2'

#pfjets = "hltAK4PFJets" #original ak4PFJetsCHS
#calojets = "hltAK4CaloJets" #original ak4CaloJets
#puppijets = "hltAK4PFPuppiJets"
#PFDeepCSVTags = "hltDeepCombinedSecondaryVertexBPFPatJetTags" # original: pfDeepCSVJetTags
# PFDeepFlavourTags = "hltPFDeepFlavourJetTags" # original: pfDeepFlavourJetTagsSlimmedDeepFlavour
# PFDeepFlavourTagInfos = 'hltPFDeepFlavour'
PFDeepFlavourTags = "hltPFDeepFlavourPatJetTags" # original: pfDeepFlavourJetTagsSlimmedDeepFlavour
PFDeepFlavourTagInfos = 'hltPFDeepFlavourPat'
PFDeepCSVTags = "hltDeepCombinedSecondaryVertexBPFPatJetTags"
IPTagInfos = 'hltDeepBLifetimePFPat'
SVTagInfos = 'hltDeepSecondaryVertexPFPat'

PuppiDeepCSVTags = 'hltDeepCombinedSecondaryVertexBPFPuppiPatJetTags'
PuppiDeepFlavourTags = 'hltPFPuppiDeepFlavourJetTags'
PuppiDeepFlavourTagInfos = 'hltPFPuppiDeepFlavour'
PuppiIPTagInfos = 'hltDeepBLifetimePFPuppiPat'
SVPuppiTagInfos = 'hltDeepSecondaryVertexPFPuppiPat'

rho = "hltFixedGridRhoFastjetAll" #original fixedGridRhoFastjetAll
hltVertices = "hltVerticesPFFilter" #original offlinePrimaryVertices
hltVerticesSlimmed = "hltVerticesPFFilter" #original offlineSlimmedPrimaryVertices
siPixelClusters = "hltSiPixelClusters" #original siPixelClusters
ecalRecHit = "hltEcalRecHit" #original ecalRecHit
hbhereco = "hltHbhereco" #original hbhereco
hfreco = "hltHfreco" #original hfreco
horeco = "hltHoreco" #original horeco
rpcRecHits = "hltRpcRecHits" #original rpcRecHits
tracks = "hltPFMuonMerging" #original generalTracks

genParticles = 'genParticles'
patJetSource = 'hltPatJets'
patCaloJetSource = 'hltPatJetsCalo'
patPuppiJetSource = 'hltPatJetsPuppi'
genJetCollection = 'ak4GenJetsNoNu'
pfCandidates = 'hltParticleFlow'
pvSource = hltVertices
svSource = 'hltDeepInclusiveMergedVerticesPF'
muSource = 'hltMuons'
elSource = 'hltEgammaGsfElectrons'
trackSource = tracks


###
### HLT configuration
###
if options.reco == 'HLT_GRun':
    # from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_12_0_0_GRun_V3_configDump_GT import cms, process

elif options.reco == 'HLT_Run3TRK':
    # (a) Run-3 tracking: standard
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    process = customizeHLTRun3Tracking(process)
elif options.reco == 'HLT_Run3TRKMod':
    # (a) Run-3 tracking: standard
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    process = customizeHLTRun3Tracking(process)
elif options.reco == 'HLT_Run3TRKMod2':
    # (a) Run-3 tracking: standard
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    process = customizeHLTRun3Tracking(process)
elif options.reco == 'HLT_Run3TRKWithPU':
    # (b) Run-3 tracking: all pixel vertices
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3TrackingAllPixelVertices
    process = customizeHLTRun3TrackingAllPixelVertices(process)
elif options.reco == 'HLT_Run3TRKPixelOnly':
    # (c) Run-3 tracking: pixel only tracks
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    from RecoBTag.PerformanceMeasurements.customise_TRK import *
    process = customizeHLTRun3Tracking(process)
    process = customisePFForPixelTracks(process)
elif options.reco == 'HLT_Run3TRKPixelOnlyCleaned':
    # (d) Run-3 tracking: pixel only tracks and trimmed with PVs
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    from RecoBTag.PerformanceMeasurements.customise_TRK import *
    process = customizeHLTRun3Tracking(process)
    process = customisePFForPixelTracksCleaned(process, "hltPixelTracksCleanForBTag")
elif options.reco == 'HLT_Run3TRKPixelOnlyCleaned2':
    # (d) Run-3 tracking: pixel only tracks and trimmed with PVs
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    from RecoBTag.PerformanceMeasurements.customise_TRK import *
    process = customizeHLTRun3Tracking(process)
    process = customisePFForPixelTracksCleaned(process, "hltPixelTracksCleanForBTag", vertex="hltTrimmedPixelVertices", nVertices = 2)
elif options.reco == 'HLT_Run3TRKPixelOnlyCleaned3':
    # (d) Run-3 tracking: pixel only tracks and trimmed with PVs
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    from RecoBTag.PerformanceMeasurements.customise_TRK import *
    process = customizeHLTRun3Tracking(process)
    process = customisePFForPixelTracksCleaned(process, "hltPixelTracksCleanForBTag", vertex="hltPixelVertices", nVertices = 4)
elif options.reco == 'HLT_Run3TRKPixelOnlyCleaned4':
    # (d) Run-3 tracking: pixel only tracks and trimmed with PVs
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    from RecoBTag.PerformanceMeasurements.customise_TRK import *
    process = customizeHLTRun3Tracking(process)
    process = customisePFForPixelTracksCleaned(process, "hltPixelTracksCleanForBTag", vertex="hltPixelVertices", nVertices = 2)
elif options.reco == 'HLT_BTagROI':
    # (e) Run-3 tracking: ROI PF approach
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from RecoBTag.PerformanceMeasurements.customise_hlt import *
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    process = customizeHLTRun3Tracking(process)
    process = addPaths_PFJetsForBtag(process)
    pvSource                 = "hltVerticesPFFilterForBTag"
    pfCandidates             = 'hltParticleFlowForBTag'
    patJetSource             = 'hltPatJetsROI'
    trackSource              = "hltMergedTracksForBTag"
    PFDeepFlavourTags        = "hltPFDeepFlavourROIJetTags"
    PFDeepFlavourTagInfos    = 'hltPFDeepFlavourROI'
    rho                      = "hltFixedGridRhoFastjetAllForBTag"
    patPuppiJetSource        = 'hltPatJetsPuppiROI'
    PFDeepCSVTags            = "hltDeepCombinedSecondaryVertexBPFPatROIJetTags"
    PuppiDeepCSVTags         = 'hltDeepCombinedSecondaryVertexBPFPuppiPatROIJetTags'
    PuppiDeepFlavourTags     = 'hltPFPuppiDeepFlavourROIJetTags'
    PuppiDeepFlavourTagInfos = 'hltPFPuppiDeepFlavourROI'
    PuppiIPTagInfos          = 'hltDeepBLifetimePFPuppiPatROI'
    IPTagInfos               = 'hltDeepBLifetimePFPatROI'
    SVPuppiTagInfos          = 'hltDeepSecondaryVertexPFPuppiPatROI'
    SVTagInfos               = 'hltDeepSecondaryVertexPFPatROI'
elif options.reco == 'HLT_BTagROIPixelTracks':
    # (e) Run-3 tracking: ROI PF approach
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process
    from RecoBTag.PerformanceMeasurements.customise_hlt import *
    from HLTrigger.Configuration.customizeHLTRun3Tracking import customizeHLTRun3Tracking
    process = customizeHLTRun3Tracking(process)
    process = addPatatracksForROI(process)
    process = addPaths_PFJetsForBtag(process)
    pvSource                 = "hltVerticesPFFilterForBTag"
    pfCandidates             = 'hltParticleFlowForBTag'
    patJetSource             = 'hltPatJetsROI'
    trackSource              = "hltMergedTracksForBTag"
    PFDeepFlavourTags        = "hltPFDeepFlavourROIJetTags"
    PFDeepFlavourTagInfos    = 'hltPFDeepFlavourROI'
    rho                      = "hltFixedGridRhoFastjetAllForBTag"
    patPuppiJetSource        = 'hltPatJetsPuppiROI'
    PFDeepCSVTags            = "hltDeepCombinedSecondaryVertexBPFPatROIJetTags"
    PuppiDeepCSVTags         = 'hltDeepCombinedSecondaryVertexBPFPuppiPatROIJetTags'
    PuppiDeepFlavourTags     = 'hltPFPuppiDeepFlavourROIJetTags'
    PuppiDeepFlavourTagInfos = 'hltPFPuppiDeepFlavourROI'
    PuppiIPTagInfos          = 'hltDeepBLifetimePFPuppiPatROI'
    IPTagInfos               = 'hltDeepBLifetimePFPatROI'
    SVPuppiTagInfos          = 'hltDeepSecondaryVertexPFPuppiPatROI'
    SVTagInfos               = 'hltDeepSecondaryVertexPFPatROI'
else:
  raise RuntimeError('keyword "reco = '+options.reco+'" not recognised')



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


# list of patterns to determine paths to keep
keepPaths = [
  'MC_*Jets*',
  'MC_*MET*',
  'MC_*AK8Calo*',
  'MC_*DeepCSV*',
  'MC_*DeepJet*',

  # 'HLT_PFJet*_v*',
  # 'HLT_AK4PFJet*_v*',
  # 'HLT_AK8PFJet*_v*',
  # 'HLT_PFHT*_v*',
  # 'HLT_PFMET*_PFMHT*_v*',
]

# list of paths that are kept
listOfPaths = []
print ("Keep paths:")
print ('-'*108)
# remove selected cms.Path objects from HLT config-dump
for _modname in sorted(process.paths_()):
    _keepPath = False
    for _tmpPatt in keepPaths:
        _keepPath = fnmatch.fnmatch(_modname, _tmpPatt)
        if _keepPath: break
    if _keepPath:
        print ('{:<99} | {:<4} |'.format(_modname, '+'))
        listOfPaths.append(_modname)
        continue
    _mod = getattr(process, _modname)
    if type(_mod) == cms.Path:
        process.__delattr__(_modname)
        # if options.verbosity > 0:
        #     print '{:<99} | {:<4} |'.format(_modname, '')
print ('-'*108)

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
if options.runCaloJetVariables:
    process = addPaths_MC_JMECalo(process)
process = addPaths_MC_JMEPF(process)
# process = addPaths_MC_JMEPFCluster(process)
if options.runPuppiJetVariables:
    process = addPaths_MC_JMEPFPuppi(process)
from RecoBTag.PerformanceMeasurements.customise_TRK import addDeepJet
process = addDeepJet(process, doPF = True, doPuppi = options.runPuppiJetVariables)
from RecoBTag.PerformanceMeasurements.PATLikeConfig import customizePFPatLikeJets
process = customizePFPatLikeJets(process, runPF=True, runCalo=options.runCaloJetVariables, runPuppi=options.runPuppiJetVariables)


if options.reco == 'HLT_Run3TRKMod':
    process = customizeVertices(process)

if options.reco == 'HLT_Run3TRKMod2':
    process = customizeVertices2(process)

if "HLT_Run3TRKPixelOnly" in options.reco:
    process = customizeMinHitsAndPt(process)

if options.reco == 'HLT_BTagROI':
    from RecoBTag.PerformanceMeasurements.customise_hlt import *
    process = addPaths_MC_JMEPFPuppiROI(process)

    from RecoBTag.PerformanceMeasurements.ROIPATLikeConfig import customizePFPatLikeJetsROI
    process = customizePFPatLikeJetsROI(process)



## ES modules for PF-Hadron Calibrations
import os
# from CondCore.DBCommon.CondDBSetup_cfi import *
from CondCore.CondDB.CondDB_cfi import CondDB as _CondDB

process.pfhcESSource = cms.ESSource('PoolDBESSource',
  _CondDB.clone(connect = 'sqlite_fip:JMETriggerAnalysis/NTuplizers/data/PFHC_Run3Winter20_HLT_v01.db'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('PFCalibrationRcd'),
      tag = cms.string('PFCalibration_HLT_mcRun3_2021'),
      label = cms.untracked.string('HLT'),
    ),
  ),
)

process.pfhcESPrefer = cms.ESPrefer('PoolDBESSource', 'pfhcESSource')
## ES modules for HLT JECs
process.jescESSource = cms.ESSource('PoolDBESSource',
  _CondDB.clone(connect = 'sqlite_fip:JMETriggerAnalysis/NTuplizers/data/JESC_Run3Winter20_V1_MC.db'),
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
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
#process.MessageLogger.cerr.default.limit = 10
#process.MessageLogger.cerr.default.limit = 1

process.MessageLogger.suppressWarning = cms.untracked.vstring(
        'hltPatJetFlavourAssociationCalo'
)
process.MessageLogger.suppressError = cms.untracked.vstring(
        'hltPatJetFlavourAssociationCalo'
)
process.MessageLogger.cerr.threshold = "DEBUG"
#process.MessageLogger.debugModules = ["hltBTagPFDeepCSV4p06SingleROI","hltDeepCombinedSecondaryVertexBJetTagsPFROI","hltDeepCombinedSecondaryVertexBPFPatROIJetTags","hltPatJetsROI"]


if options.inputFiles:
    process.source.fileNames = options.inputFiles
process.source.secondaryFileNames = cms.untracked.vstring()

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

if options.runTiming:
    # multi-threading settings
    process.options.numberOfThreads = cms.untracked.uint32(options.numThreads if (options.numThreads > 1) else 1)
    process.options.numberOfStreams = cms.untracked.uint32(options.numStreams if (options.numStreams > 1) else 1)
    process.options.sizeOfStackForThreadsInKB = cms.untracked.uint32(10240)
else:
    # multi-threading settings
    process.options.numberOfThreads = cms.untracked.uint32(options.numThreads if (options.numThreads > 1) else 1)
    process.options.numberOfStreams = cms.untracked.uint32(options.numStreams if (options.numStreams > 1) else 1)

#Set GT by hand:
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# from Configuration.AlCa.GlobalTag import GlobalTag
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
            print ("Changing 'minimumNumberOfHits' for " + m + " to " + str(options.minNumberOfHits))
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


process.btagana.deepCSVBJetTags = PFDeepCSVTags
process.btagana.deepCSVBPuppiJetTags = PuppiDeepCSVTags

process.btagana.deepFlavourJetTags    = PFDeepFlavourTags
process.btagana.deepFlavourTagInfos   = PFDeepFlavourTagInfos

process.btagana.deepFlavourPuppiJetTags    = PuppiDeepFlavourTags
process.btagana.deepFlavourPuppiTagInfos = PuppiDeepFlavourTagInfos

process.btagana.ipPuppiTagInfos = PuppiIPTagInfos
process.btagana.ipTagInfos = IPTagInfos

process.btagana.svPuppiTagInfos = SVPuppiTagInfos
process.btagana.svTagInfos = SVTagInfos

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
#--------


# from JMETriggerAnalysis.Common.TrackHistogrammer_cfi import TrackHistogrammer
# process.TrackHistograms_hltPixelTracks = TrackHistogrammer.clone(src = 'hltPixelTracks')
# process.TrackHistograms_hltTracks = TrackHistogrammer.clone(src = 'hltPFMuonMerging')
# process.TrackHistograms_hltMergedTracks = TrackHistogrammer.clone(src = 'hltMergedTracks')
# # process.TrackHistograms_hltGeneralTracks = TrackHistogrammer.clone(src = 'generalTracks')

# process.trkMonitoringSeq = cms.Sequence(
#    process.TrackHistograms_hltPixelTracks
#  + process.TrackHistograms_hltTracks
#  + process.TrackHistograms_hltMergedTracks
#  # + process.TrackHistograms_hltGeneralTracks
# )

# process.trkMonitoringEndPath = cms.EndPath(process.trkMonitoringSeq)
# process.schedule.extend([process.trkMonitoringEndPath])

# from JMETriggerAnalysis.Common.VertexHistogrammer_cfi import VertexHistogrammer
# process.VertexHistograms_hltDeepInclusiveVertexFinderPF = VertexHistogrammer.clone(src = 'hltDeepInclusiveVertexFinderPF')
# process.VertexHistograms_hltDeepSecondaryVertexPFPatTagInfos = VertexHistogrammer.clone(src = 'hltDeepSecondaryVertexPFPatTagInfos')
# process.VertexHistograms_hltDeepInclusiveMergedVerticesPF = VertexHistogrammer.clone(src = 'hltDeepInclusiveMergedVerticesPF')
# process.VertexHistograms_hltDeepInclusiveSecondaryVerticesPF = VertexHistogrammer.clone(src = 'hltDeepInclusiveSecondaryVerticesPF')
#
# process.vtxMonitoringSeq = cms.Sequence(
#    process.VertexHistograms_hltDeepInclusiveVertexFinderPF
#  + process.VertexHistograms_hltDeepSecondaryVertexPFPatTagInfos
#  + process.VertexHistograms_hltDeepInclusiveMergedVerticesPF
#  + process.VertexHistograms_hltDeepInclusiveSecondaryVerticesPF
# )

if options.runTiming:
    #timing test
    from HLTrigger.Timer.FastTimer import customise_timer_service_print,customise_timer_service,customise_timer_service_singlejob
    # process = customise_timer_service_print(process)
    # process = customise_timer_service(process)
        # remove any instance of the FastTimerService
    if 'FastTimerService' in process.__dict__:
        del process.FastTimerService

    # instrument the menu with the FastTimerService
    process.load("HLTrigger.Timer.FastTimerService_cfi")

    # print a text summary at the end of the job
    process.FastTimerService.printEventSummary        = False
    process.FastTimerService.printRunSummary          = False
    process.FastTimerService.printJobSummary          = True
    # enable DQM plots
    process.FastTimerService.enableDQM                = True
    # enable per-path DQM plots (starting with CMSSW 9.2.3-patch2)
    process.FastTimerService.enableDQMbyPath          = True
    # enable per-module DQM plots
    process.FastTimerService.enableDQMbyModule        = True
    # enable DQM plots vs lumisection
    process.FastTimerService.enableDQMbyLumiSection   = True
    process.FastTimerService.dqmLumiSectionsRange     = 2500    # lumisections (23.31 s)
    # set the time resolution of the DQM plots
    process.FastTimerService.dqmTimeRange             = 1000.   # ms
    process.FastTimerService.dqmTimeResolution        =    5.   # ms
    process.FastTimerService.dqmPathTimeRange         =  100.   # ms
    process.FastTimerService.dqmPathTimeResolution    =    0.5  # ms
    process.FastTimerService.dqmModuleTimeRange       =   40.   # ms
    process.FastTimerService.dqmModuleTimeResolution  =    0.2  # ms
    # set the base DQM folder for the plots
    process.FastTimerService.dqmPath                  = "HLT/TimerService"
    process.FastTimerService.enableDQMbyProcesses     = False
    # enable text dump
    if not hasattr(process,'MessageLogger'):
        process.load('FWCore.MessageService.MessageLogger_cfi')
    process.MessageLogger.categories.append('FastReport')
    process.MessageLogger.cerr.FastReport = cms.untracked.PSet( limit = cms.untracked.int32( 10000000 ) )
    # save the DQM plots in the DQMIO format
    process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
        fileName = cms.untracked.string("DQM.root")
    )
    process.FastTimerOutput = cms.EndPath(process.dqmOutput)
    # process.schedule.append(process.FastTimerOutput)

    # process = customise_timer_service_singlejob(process)
    process.FastTimerService.dqmTimeRange            = 20000.
    process.FastTimerService.dqmTimeResolution       =    10.
    process.FastTimerService.dqmPathTimeRange        = 10000.
    process.FastTimerService.dqmPathTimeResolution   =     5.
    process.FastTimerService.dqmModuleTimeRange      =  1000.
    process.FastTimerService.dqmModuleTimeResolution =     1.
    # process.dqmOutput.fileName = cms.untracked.string(options.output)

process.p = cms.Path(
    process.allEvents
    # * process.filtSeq
    * process.selectedEvents
    * process.analyzerSeq
    # * process.trkMonitoringSeq
    # * process.vtxMonitoringSeq
)

# if options.runTiming:
#     process.p *= process.FastTimerOutput

# Delete predefined output module (needed for running with CRAB)
del process.out
# dump content of cms.Process to python file
if options.dumpPython is not None:
   open(options.dumpPython, 'w').write(process.dumpPython())

print ('')
print ('option: output =', options.outFilename)
print ('option: reco =', options.reco)
print ('option: dumpPython =', options.dumpPython)
print ('')
# print 'process.GlobalTag =', process.GlobalTag.dumpPython()
print ('process.GlobalTag =', process.GlobalTag.globaltag)
print ('process.source =', process.source.dumpPython())
print ('process.maxEvents =', process.maxEvents.input)
# print 'process.options =', process.options.dumpPython()
print ('-------------------------------')

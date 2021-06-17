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
options.register('globalTag', '112X_mcRun3_2021_realistic_v14',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag, no default value provided"
)
options.register('runCaloJetVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Jet Variables'
)
options.register('runPuppiJetVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Jet Variables'
)
options.register('runTiming', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'run timing instead of rates'
)
options.register('numThreads', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'number of threads'
)
options.register('numStreams', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'number of streams'
)
options.register(
    'reco', 'HLT_GRun',
  VarParsing.multiplicity.singleton,
  VarParsing.varType.string,
  'keyword to define HLT reconstruction'
)
options.setDefault('maxEvents', -1)
options.parseArguments()
globalTag = options.globalTag


###
### HLT configuration
###
if options.reco == 'HLT_GRun':
    from RecoBTag.PerformanceMeasurements.Configs.HLT_dev_CMSSW_11_2_0_GRun_V19_configDump import cms, process

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
# _keepPath = _modname.startswith('MC_') and ('Jets' in _modname or 'MET' in _modname or 'DeepCSV' in _modname or 'DeepFlavour' in _modname or 'Tracking' in _modname or 'AK8Calo' in _modname)
keepPaths = [
  # 'MC_*Jets*',
  # 'MC_*PFJets*',
  # 'MC_*CaloJets*',
  # 'MC_*MET*',
  # 'MC_*AK8Calo*',
  # 'MC_*DeepCSV*',
  # 'MC_*CaloBTagDeepCSV*',
  'MC_*PFBTagDeepCSV*',
  # 'MC_*DeepJet*',
  # 'HLT_PFJet*_v*',
  # 'HLT_AK4PFJet*_v*',
  # 'HLT_AK8PFJet*_v*',
  # 'HLT_PFHT*_v*',
  # 'HLT_PFMET*_PFMHT*_v*',

  'HLT_*DeepCSV*_v*',
]
# list of paths that are kept
listOfPaths = []
# remove selected cms.Path objects from HLT config-dump
print "Keep paths:"
print '-'*108
for _modname in sorted(process.paths_()):
    _keepPath = False
    for _tmpPatt in keepPaths:
        _keepPath = fnmatch.fnmatch(_modname, _tmpPatt)
        if _keepPath: break
    if _keepPath:
        print '{:<99} | {:<4} |'.format(_modname, '+')
        listOfPaths.append(_modname)
        continue
    _mod = getattr(process, _modname)
    if type(_mod) == cms.Path:
        process.__delattr__(_modname)
        # if options.verbosity > 0:
        #     print '{:<99} | {:<4} |'.format(_modname, '')
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
# process = addPaths_MC_JMECalo(process)
# process = addPaths_MC_JMEPF(process)
# process = addPaths_MC_JMEPFCluster(process)
if options.runPuppiJetVariables:
    process = addPaths_MC_JMEPFPuppi(process)
from RecoBTag.PerformanceMeasurements.customise_TRK import addDeepJet
# process = addDeepJet(process, doPF = True, doPuppi = options.runPuppiJetVariables)
# from RecoBTag.PerformanceMeasurements.PATLikeConfig import customizePFPatLikeJets
# process = customizePFPatLikeJets(process, runPF=True, runCalo=options.runCaloJetVariables, runPuppi=options.runPuppiJetVariables)


if options.reco == 'HLT_Run3TRKMod':
    process = customizeVertices(process)

if options.reco == 'HLT_Run3TRKMod2':
    process = customizeVertices2(process)

if "HLT_Run3TRKPixelOnly" in options.reco:
    process = customizeMinHitsAndPt(process)

if options.reco == 'HLT_BTagROI':
    if options.runPuppiJetVariables:
        from RecoBTag.PerformanceMeasurements.customise_hlt import *
        process = addPaths_MC_JMEPFPuppiROI(process)

    # from RecoBTag.PerformanceMeasurements.ROIPATLikeConfig import customizePFPatLikeJetsROI
    # process = customizePFPatLikeJetsROI(process)



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


if options.inputFiles:
    process.source.fileNames = options.inputFiles
process.source.secondaryFileNames = cms.untracked.vstring()

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outFilename)
)

## Events to process
# process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
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



# process.schedule_().extend([
#       process.MC_PFBTagDeepCSV_v10,
#       process.MC_PFBTagDeepJet,
# ])















# del process.out
# dump content of cms.Process to python file
if options.dumpPython is not None:
   open(options.dumpPython, 'w').write(process.dumpPython())

print ''
print 'option: output =', options.outFilename
print 'option: reco =', options.reco
print 'option: dumpPython =', options.dumpPython
print ''
# print 'process.GlobalTag =', process.GlobalTag.dumpPython()
print 'process.GlobalTag =', process.GlobalTag.globaltag
print 'process.source =', process.source.dumpPython()
print 'process.maxEvents =', process.maxEvents.input
# print 'process.options =', process.options.dumpPython()
print '-------------------------------'

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import copy
from pdb import set_trace

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
options.register('usePuppiForFatJets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use Puppi for fat jets"
)
options.register('usePuppiForBTagging', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use Puppi candidates for b tagging"
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
options.register('runJetClustering', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Cluster jets from scratch instead of using those already present in the event"
)
options.register('runFatJetClustering', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Cluster fat jets from scratch instead of using those already present in the event"
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
options.register('runEventInfo', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run Event Info"
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
# options.register('fatJetAbsEtaMax', 2.5,
options.register('fatJetAbsEtaMax', 4.5,
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
options.register('remakeAllDiscr', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Remake all b-tag discriminator, including those already stored in MiniAOD"
)
options.register('remakeDoubleB', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Remake double-b-tag discriminator"
)
options.register('fastSim', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Running using FastSim"
)
options.register('useSelectedTracks', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "If you want to run on all tracks: False for commissioning studies"
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
options.register('usePruned', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use pruned jets"
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
options.register('useNegativeDeepFlavourTags', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Include negative deep flavour jet taggers"
)

## Generally leave to False unless you know what you are doing
options.register('runIVF', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run IVF, currently leave to False!"
)
## Master switch for boosted b tag commissioning: overrider several other switches
options.register('doBoostedCommissioning', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Make NTuples with branches for boosted b tag commissioning: overrider several other switches"
)
## Do Ctag
options.register('runCTagVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Make NTuples with branches for CTag"
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
# options.register('maxJetEta', 2.5,
options.register('maxJetEta', 4.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum jet |eta| (default is 2.5)"
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
options.register('runCSVTagVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run CSV TaggingVariables')
options.register('runCSVTagVariablesSubJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run CSV TaggingVariables for SubJets')
options.register('runTagVariablesSubJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run TagVariables for SubJets')
options.register('runCSVTagTrackVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run CSV Tagging Track Variables')
options.register('runDeepFlavourTagVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run DeepFlavour TaggingVariables')
options.register('runPFElectronVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run PF Electron Variables')
options.register('runPFMuonVariables', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run PF Muon Variables')
options.register('runPatMuons', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'True if you want to run Pat Muon Variables')
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
options.register('numThreads', 1,
              VarParsing.multiplicity.singleton,
              VarParsing.varType.int,
              'number of threads')

options.register('numStreams', 1,
              VarParsing.multiplicity.singleton,
              VarParsing.varType.int,
              'number of streams')

options.register('logs', False,
              VarParsing.multiplicity.singleton,
              VarParsing.varType.bool,
              'create log files configured via MessageLogger')

options.register('dumpPython', None,
              VarParsing.multiplicity.singleton,
              VarParsing.varType.string,
              'Path to python file with content of cms.Process')




## 'maxEvents' is already registered by the Framework, changing default value

options.setDefault('maxEvents', -1)
#options.setDefault('maxEvents', 100)


options.parseArguments()
if options.defaults:
	from importlib import import_module
	try:
		defaults = import_module('RecoBTag.PerformanceMeasurements.defaults.%s' % options.defaults)
	except ImportError:
		raise ValueError('The default settings named %s.py are not present in PerformanceMeasurements/python/defaults/' % options.defaults)
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

from RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff import *
btagana_tmp = bTagAnalyzer.clone()
print('Storing the variables from the following groups:')
options_to_change = set() #store which swtiches we need on
for requiredGroup in options.groups:
  print(requiredGroup)
  found=False
  for existingGroup in btagana_tmp.groups:
    if(requiredGroup==existingGroup.group):
      existingGroup.store=True
      for var in existingGroup.variables:
        if "FatJetInfo." in var:
          options_to_change.update({"runFatJets"})
          var = var.split(".")[1]
        if "SubJetInfo." in var:
          options_to_change.update({"runSubJets"})
          var = var.split(".")[1]
        options_to_change.update([i for i in variableDict[var].runOptions])
      found=True
      break
  if(not found):
    print('WARNING: The group ' + requiredGroup + ' was not found')

#change values accordingly
for switch in options_to_change:
  if switch not in options._beenSet:
    raise ValueError('The option set by the variables: %s does not exist among the cfg options!' % switch)
  elif not options._beenSet[switch]:
    print 'Turning on %s, as some stored variables demands it' % switch
    setattr(options, switch, True)


## Use either PFchs or Puppi
if options.usePFchs and options.usePuppi:
    print "WARNING: Both usePFchs and usePuppi set to True. Giving priority to Puppi."
    options.usePFchs = False

## Resolve potential conflicts in Puppi usage
if options.usePuppi and not options.usePuppiForFatJets:
    print "WARNING: usePuppi set to True while usePuppiForFatJets set to False. Puppi will be used for all jet types."
    options.usePuppiForFatJets = True

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Running using FastSim samples: %s"%('True' if options.fastSim else 'False')
print "Running on MiniAOD: %s"%('True' if options.miniAOD else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')
print "Using Puppi: %s"%('True' if options.usePuppi else 'False')
print "Using Puppi for fat jets: %s"%('True' if options.usePuppiForFatJets else 'False')
print "Using Puppi for b tagging: %s"%('True' if (options.usePuppi and options.usePuppiForBTagging) else 'False')

## Subjets only stored when also running over fat jets
if options.runSubJets and not options.runFatJets:
    print "WARNING: You are attempting to store subjet information without running over fat jets. Please enable running over fat jets in order to store the subjet information."
    options.runSubJets = False

if not options.miniAOD and options.runDeepFlavourTagVariables: #FIXME
    print "WARNING: switching off DeepFlavour, as it is not supported in AOD"
    options.runDeepFlavourTagVariables = False

if options.doBoostedCommissioning:
    print "**********NTuples will be made for boosted b tag commissioning. The following switches will be reset:**********"
    options.processStdAK4Jets=False
    print "Option processStdAK4Jets will be set to '",options.processStdAK4Jets,"'"
    options.runFatJets=True
    options.runSubJets = True
    print "Option runFatJets will be set to '",options.runFatJets,"'"
    print "Option runSubJets  will be set to '",options.runSubJets,"'"
    print "********************"
if options.runCTagVariables:
    print "**********You are making NTuple for CTag*************"

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK8 = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
if options.usePFchs:
    jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
if options.usePuppi:
    jetCorrectionsAK4 = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], 'None')
if options.usePuppiForFatJets:
    jetCorrectionsAK8 = ('AK8PFPuppi', ['L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK4[1].append('L2L3Residual')
    jetCorrectionsAK8[1].append('L2L3Residual')

jetCorrectionsSubJets = copy.deepcopy(jetCorrectionsAK4)
if not options.usePuppi and options.usePuppiForFatJets:
    jetCorrectionsSubJets = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], 'None')

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
   ,'pfDeepFlavourTagInfos'
]
bTagInfos_noDeepFlavour = bTagInfos[:-1]
## b-tag discriminators
bTagDiscriminatorsLegacy = set([
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
])
bTagDiscriminators = set([
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
    # DeepCSV
  , 'pfDeepCSVJetTags:probudsg'
  , 'pfDeepCSVJetTags:probb'
  , 'pfDeepCSVJetTags:probc'
  , 'pfDeepCSVJetTags:probbb'
  , 'pfNegativeDeepCSVJetTags:probudsg'
  , 'pfNegativeDeepCSVJetTags:probb'
  , 'pfNegativeDeepCSVJetTags:probc'
  , 'pfNegativeDeepCSVJetTags:probbb'
  , 'pfPositiveDeepCSVJetTags:probudsg'
  , 'pfPositiveDeepCSVJetTags:probb'
  , 'pfPositiveDeepCSVJetTags:probc'
  , 'pfPositiveDeepCSVJetTags:probbb'
    # DeepFlavour
  , 'pfDeepFlavourJetTags:probb'
  , 'pfDeepFlavourJetTags:probbb'
  , 'pfDeepFlavourJetTags:problepb'
  , 'pfDeepFlavourJetTags:probc'
  , 'pfDeepFlavourJetTags:probuds'
  , 'pfDeepFlavourJetTags:probg'
  , 'pfNegativeDeepFlavourJetTags:probb'
  , 'pfNegativeDeepFlavourJetTags:probbb'
  , 'pfNegativeDeepFlavourJetTags:problepb'
  , 'pfNegativeDeepFlavourJetTags:probc'
  , 'pfNegativeDeepFlavourJetTags:probuds'
  , 'pfNegativeDeepFlavourJetTags:probg'
])

## Legacy taggers not supported with MiniAOD
if options.miniAOD and options.useLegacyTaggers:
    print "WARNING: Legacy taggers not supported with MiniAOD"
    options.useLegacyTaggers = False

## If using legacy taggers
if options.useLegacyTaggers:
    bTagInfos = bTagInfosLegacy
    bTagDiscriminators = bTagDiscriminatorsLegacy

## If not including negative deep flavour jet taggers
if not options.useNegativeDeepFlavourTags:
  bTagDiscriminators = {i for i in bTagDiscriminators if 'NegativeDeepFlavourJetTags' not in i}

## Clustering algorithm label
algoLabel = 'CA'
if options.jetAlgo == 'AntiKt':
    algoLabel = 'AK'

## Figure out if jet clustering is needed
if not options.runFatJetClustering:
    options.runFatJetClustering = options.runJetClustering
if not options.miniAOD and options.usePuppi and not options.runJetClustering:
    print "WARNING: You requested Puppi jets which are not stored in AOD. Enabling jet clustering."
    options.runJetClustering = True

if options.miniAOD and not (options.usePFchs or options.usePuppi) and not options.runJetClustering:
    print "WARNING: You requested non-PU-subtracted jets which are not stored in MiniAOD. Enabling jet clustering."
    options.runJetClustering = True

print "Jet clustering: %s"%('True' if options.runJetClustering else 'False')

## Figure out if fat jet clustering is needed
if options.runFatJets and (options.jetAlgo != 'AntiKt' or options.fatJetRadius != 0.8) and not options.runFatJetClustering:
    print "WARNING: You requested fat jets with an algorithm or size not stored in any of the data-tiers. Enabling fat jet clustering."
    options.runFatJetClustering = True

if options.runFatJets and not (options.usePFchs or options.usePuppi or options.usePuppiForFatJets) and not options.runFatJetClustering:
    print "WARNING: You requested non-PU-subtracted fat jets which are not stored in any of the data-tiers. Enabling fat jet clustering."
    options.runFatJetClustering = True

if options.miniAOD and options.runFatJets and options.usePFchs and not options.usePuppiForFatJets and not options.runFatJetClustering:
    print "WARNING: You requested CHS fat jets which are not stored in MiniAOD. Enabling fat jet clustering."
    options.runFatJetClustering = True

if not options.miniAOD and options.runFatJets and (options.usePuppi or options.usePuppiForFatJets) and not options.runFatJetClustering:
    print "WARNING: You requested Puppi fat jets which are not stored in AOD. Enabling fat jet clustering."
    options.runFatJetClustering = True

if options.miniAOD and options.runSubJets and options.usePruned and not options.runFatJetClustering:
    print "WARNING: You requested pruned subjets which are not stored in MiniAOD. Enabling fat jet clustering."
    options.runFatJetClustering = True

if not options.miniAOD and options.runSubJets and options.usePruned and not options.runFatJetClustering:
    print "WARNING: You requested pruned subjets which are not stored in AOD. Will run pruned fat jet clustering."

print "Fat jet clustering: %s"%('True' if options.runFatJetClustering else 'False')

## For fat jets we want to re-run all taggers in order to use the setup adapted to the larger fat jet cone size
bTagInfosFat = copy.deepcopy(bTagInfos_noDeepFlavour)
bTagInfosFat += ([] if options.useLegacyTaggers else ['pfImpactParameter' + ('CA15' if algoLabel=='CA' else 'AK8') + 'TagInfos'])
bTagInfosFat += ([] if options.useLegacyTaggers else ['pfInclusiveSecondaryVertexFinder' + ('CA15' if algoLabel=='CA' else 'AK8') + 'TagInfos'])
bTagInfosFat += ([] if options.useLegacyTaggers else ['pfBoostedDoubleSV' + ('CA15' if algoLabel=='CA' else 'AK8') + 'TagInfos'])

bTagDiscriminators_no_deepFlavour = {i for i in bTagDiscriminators if 'DeepFlavourJetTags' not in i}
bTagDiscriminatorsFat = copy.deepcopy(bTagDiscriminators_no_deepFlavour)
## Add DeepDoubleB tagger to fat jets
bTagDiscriminatorsFat.update(set(['pfDeepDoubleBJetTags:probH']))

if options.runJetClustering:
    options.remakeAllDiscr = True
if options.runFatJetClustering:
    options.remakeDoubleB = True
if options.remakeDoubleB:
    bTagDiscriminatorsFat.update(set([]) if options.useLegacyTaggers else set(['pfBoostedDoubleSecondaryVertex' + ('CA15' if algoLabel=='CA' else 'AK8') + 'BJetTags']))

## Full list of bTagDiscriminators for SoftDrop subjets
bTagDiscriminatorsSubJets  = copy.deepcopy(bTagDiscriminators_no_deepFlavour)
bTagDiscriminatorsSoftDrop = copy.deepcopy(bTagDiscriminators_no_deepFlavour)

## If using MiniAOD and not reclustering jets, only run taggers not already stored (with the exception of JP taggers and DeepCSV)
if options.miniAOD and not options.runJetClustering and not options.remakeAllDiscr:
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets as patJetsDefault
    storedDiscriminators = set([x.value() for x in patJetsDefault.discriminatorSources])
    print "INFO: Removing b-tag discriminators already stored in MiniAOD (with the exception of JP taggers)"
    jptaggers = {i for i in bTagDiscriminators if 'ProbabilityBJetTags' in i or i.startswith('pfDeepCSV')}
    bTagDiscriminators = (bTagDiscriminators - storedDiscriminators) | jptaggers
if options.miniAOD and not options.runFatJetClustering and not options.remakeAllDiscr:
    ## SoftDrop subjets in MiniAOD have only CSVv2AVR and CSVv2IVF discriminators stored
    bTagDiscriminatorsSoftDrop -= {'pfCombinedSecondaryVertexV2BJetTags', 'pfCombinedInclusiveSecondaryVertexV2BJetTags'}

## Postfix
postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
jetSource = 'pfJetsPFBRECO'+postfix
patJetSource = 'selectedPatJets'+postfix
genJetCollection = 'ak4GenJetsNoNu'
pfCandidates = 'particleFlow'
pvSource = 'offlinePrimaryVertices'
svSource = 'inclusiveCandidateSecondaryVertices'
muSource = 'muons'
elSource = 'gedGsfElectrons'
patMuons = 'selectedPatMuons'
trackSource = 'generalTracks'
fatJetSource = 'fatPFJetsCHS'
fatJetSourceSoftDrop = 'fatPFJetsSoftDrop'
fatJetSourcePruned = 'fatPFJetsPruned'
fatGenJetCollection = 'genFatJetsNoNu'
fatGenJetCollectionSoftDrop = 'genFatJetsNoNuSoftDrop'
fatGenJetCollectionPruned = 'genFatJetsNoNuPruned'
patFatJetSource = 'packedPatJetsFatPF'
subJetSourceSoftDrop = 'fatPFJetsSoftDrop:SubJets'
patSubJetSourceSoftDrop = 'selectedPatJetsSoftDropFatPFPacked:SubJets'
if not options.runJetClustering:
    jetSource = ('ak4PFJetsCHS' if options.usePFchs else 'ak4PFJets')
if not options.runFatJetClustering:
    fatJetSource = 'ak8PFJetsCHS'
    fatJetSourceSoftDrop = 'ak8PFJetsCHSSoftDrop'
    fatGenJetCollection = 'ak8GenJetsNoNu'
## If running on MiniAOD
if options.miniAOD:
    genParticles = 'prunedGenParticles'
    jetSource = 'ak4Jets'
    genJetCollection = 'slimmedGenJets'
    pfCandidates = 'packedPFCandidates'
    pvSource = 'offlineSlimmedPrimaryVertices'
    svSource = 'slimmedSecondaryVertices'
    trackSource = 'unpackedTracksAndVertices'
    muSource = 'slimmedMuons'
    elSource = 'slimmedElectrons'
    patMuons = 'slimmedMuons'
    if not options.runJetClustering:
        jetSource = ('slimmedJetsPuppi' if options.usePuppi else 'slimmedJets')
        patJetSource = 'selectedUpdatedPatJets'+postfix
    if not options.runFatJetClustering:
        fatJetSource = 'slimmedJetsAK8'
        patFatJetSource = 'selectedUpdatedPatJetsFatPF'+postfix
        subJetSourceSoftDrop = 'slimmedJetsAK8PFPuppiSoftDropPacked:SubJets'
        patSubJetSourceSoftDrop = 'selectedUpdatedPatJetsSoftDropSubjetsPF'+postfix

if not options.eras:
	process = cms.Process("BTagAna")
else:
	from Configuration.StandardSequences.Eras import eras
	eras_to_use = []
	for era in options.eras:
		if hasattr(eras, era):
			eras_to_use.append(getattr(eras, era))
		else:
			raise ValueError('The requested era (%s) is not available' % era)
	process = cms.Process("BTagAna", *eras_to_use)


## MessageLogger
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# # If you run over many samples and you save the log, remember to reduce
# # the size of the output by prescaling the report of the event number
# process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
# process.MessageLogger.cerr.default.limit = 10

# MessageLogger
if options.logs:
   process.MessageLogger = cms.Service('MessageLogger',
     destinations = cms.untracked.vstring(
       'cerr',
       'logError',
       'logInfo',
       'logDebug',
     ),
     # scram b USER_CXXFLAGS="-DEDM_ML_DEBUG"
     debugModules = cms.untracked.vstring(
       'PixelVerticesSelector',
       'TracksClosestToFirstVerticesSelector',
       'JMETriggerNTuple',
     ),
     categories = cms.untracked.vstring(
       'FwkReport',
     ),
     cerr = cms.untracked.PSet(
       threshold = cms.untracked.string('WARNING'),
       FwkReport = cms.untracked.PSet(
         reportEvery = cms.untracked.int32(1),
       ),
     ),
     logError = cms.untracked.PSet(
       threshold = cms.untracked.string('ERROR'),
       extension = cms.untracked.string('.txt'),
       FwkReport = cms.untracked.PSet(
         reportEvery = cms.untracked.int32(1),
       ),
     ),
     logInfo = cms.untracked.PSet(
       threshold = cms.untracked.string('INFO'),
       extension = cms.untracked.string('.txt'),
       FwkReport = cms.untracked.PSet(
         reportEvery = cms.untracked.int32(1),
       ),
     ),
     logDebug = cms.untracked.PSet(
       threshold = cms.untracked.string('DEBUG'),
       extension = cms.untracked.string('.txt'),
       FwkReport = cms.untracked.PSet(
         reportEvery = cms.untracked.int32(1),
       ),
     ),
   )



## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

if options.miniAOD:
    process.source.fileNames = [
        #/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#        '/store/relval/CMSSW_10_4_0_mtd3/RelValTTbar_Tauola_14TeV/MINIAODSIM/PU25ns_103X_upgrade2023_realistic_v2_2023D35PU200_1-v2/20000/31C1C942-EC8D-1245-B773-2293F5CC87DB.root'
        '/store/relval/CMSSW_10_4_0_mtd3/RelValTTbar_Tauola_14TeV/MINIAODSIM/PU25ns_103X_upgrade2023_realistic_v2_2023D35PU200_5-v2/20000/A8AD8A25-CDC1-2E4B-A404-E9DB14ECF16A.root'

    ]
    if options.runOnData:
        process.source.fileNames = [
            #/JetHT/Run2017A-PromptReco-v2/MINIAOD
            '/store/data/Run2017A/JetHT/MINIAOD/PromptReco-v2/000/296/168/00000/3E20EA58-0F4D-E711-851C-02163E0139CE.root',
        ]
    if options.fastSim:
        process.source.fileNames = [
            '/store/relval/CMSSW_8_0_0/RelValTTbar_13/MINIAODSIM/PU25ns_80X_mcRun2_asymptotic_v4_FastSim-v2/10000/8E75D08A-3FDE-E511-8374-0CC47A4C8F26.root'
        ]
else:
    process.source.fileNames = [
        '/store/mc/PhaseIFall16DR/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/AODSIM/PhaseIFall16PUFlat20to50_81X_upgrade2017_realistic_v26-v1/50000/0039E945-35E3-E611-AF8D-001E675A6C2A.root'
    ]
    if options.runOnData:
        process.source.fileNames = [
            '/store/data/Run2016B/SingleMuon/AOD/PromptReco-v2/000/275/125/00000/DA2EC189-7E36-E611-8C63-02163E01343B.root'
        ]
    if options.fastSim:
        process.source.fileNames = [
            '/store/relval/CMSSW_8_0_0/RelValTTbar_13/GEN-SIM-DIGI-RECO/PU25ns_80X_mcRun2_asymptotic_v4_FastSim-v2/10000/0400D094-63DD-E511-8B51-0CC47A4C8ED8.root'
        ]
if options.inputFiles:
    process.source.fileNames = options.inputFiles

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


# multi-threading settings
process.options.numberOfThreads = cms.untracked.uint32(options.numThreads if (options.numThreads > 1) else 1)
process.options.numberOfStreams = cms.untracked.uint32(options.numStreams if (options.numStreams > 1) else 1)
if hasattr(process, 'DQMStore'):
   process.DQMStore.enableMultiThread = (process.options.numberOfThreads > 1)

## Events to process
process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    allowUnscheduled = cms.untracked.bool(True)
)

#Set GT by hand:
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#~ process.GlobalTag.globaltag = globalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
# process.GlobalTag.globaltag = '110X_mcRun4_realistic_v3' # needed for Pu200 RelVal CMSSW_11_0_0

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
    if options.runOnData: dbfile=options.jecDBFileData
    else: dbfile=options.jecDBFileMC
    print "\nUsing private SQLite file", dbfile, "\n"
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
		    connect = cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/data/"+dbfile+'.db'),
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
# trkProbaCalibTag = options.JPCalibration
# process.GlobalTag.toGet = cms.VPSet(
#     cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
#       tag = cms.string(trkProbaCalibTag),
#       connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
#     )
# )

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")


process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')

# else: process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

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

    if options.usePuppi or options.usePuppiForFatJets:
        process.load('CommonTools.PileupAlgos.Puppi_cff')
        if options.usePuppi:
            from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
            _pfJets = ak4PFJets.clone(src = cms.InputTag('puppi'), doAreaFastjet = True, srcPVs = cms.InputTag(pvSource))
            setattr(process,'pfJetsPFBRECO'+postfix,_pfJets)
        if options.usePuppiForBTagging: pfCandidates = 'puppi'
else:
    ## GenJetsNoNu selection
    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))

    ## PFchs selection
    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

    ## Reco jets
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    if options.usePFchs:
        process.ak4Jets = ak4PFJets.clone(src = cms.InputTag('pfCHS'), doAreaFastjet = True, srcPVs = cms.InputTag(pvSource))
    elif options.usePuppi:
        process.ak4Jets = ak4PFJets.clone(src = cms.InputTag('puppi'), doAreaFastjet = True, srcPVs = cms.InputTag(pvSource))
    else:
        process.ak4Jets = ak4PFJets.clone(src = cms.InputTag('packedPFCandidates'), doAreaFastjet = True, srcPVs = cms.InputTag(pvSource))

    if options.usePuppi or options.usePuppiForFatJets:
        process.load('CommonTools.PileupAlgos.Puppi_cff')
        process.puppi.candName           = cms.InputTag(pfCandidates)
        process.puppi.vertexName         = cms.InputTag(pvSource)
        process.puppi.useExistingWeights = cms.bool(True)
        process.puppi.clonePackedCands   = cms.bool(True)
        if options.usePuppiForBTagging: pfCandidates = 'puppi'

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
## Updated the default jet collection
if options.miniAOD and not options.runJetClustering:
    updateJetCollection(
        process,
        jetSource = cms.InputTag(jetSource),
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = list(bTagDiscriminators),
        explicitJTA = options.useExplicitJTA,
        postfix = postfix
    )
## Switch the default jet collection (done in order to use the above-specified b-tag infos and discriminators)
else:
    #switch off deep flavour on AOD for the moment
    switchJetCollection(
        process,
        jetSource = cms.InputTag(jetSource),
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = list(bTagInfos_noDeepFlavour), #list(bTagInfos),
        btagDiscriminators = list(bTagDiscriminators_no_deepFlavour), #bTagDiscriminators),
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag(genJetCollection),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = options.useExplicitJTA,
        postfix = postfix
    )

#-------------------------------------

#-------------------------------------
if options.runFatJets:
    if options.miniAOD and not options.runFatJetClustering:
        updateJetCollection(
            process,
            labelName='FatPF',
            jetSource=cms.InputTag(fatJetSource),
            jetCorrections = jetCorrectionsAK8,
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos = bTagInfosFat,
            btagDiscriminators = list(bTagDiscriminatorsFat),
            explicitJTA = options.useExplicitJTA,
            runIVF = options.runIVF,
            postfix = postfix
        )
        getattr(process,'selectedUpdatedPatJetsFatPF'+postfix).cut = cms.string("pt > %f && abs(eta) < %f"%(float(options.fatJetPtMin), float(options.fatJetAbsEtaMax)))
        updateJetCollection(
            process,
            labelName='SoftDropSubjetsPF',
            jetSource=cms.InputTag(subJetSourceSoftDrop),
            jetCorrections = jetCorrectionsSubJets,
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos = bTagInfos_noDeepFlavour,
            btagDiscriminators = list(bTagDiscriminatorsSoftDrop),
            explicitJTA = True,          # needed for subjet b tagging
            svClustering = False,        # needed for subjet b tagging (IMPORTANT: Needs to be set to False to disable ghost-association which does not work with slimmed jets)
            fatJets = cms.InputTag(fatJetSource), # needed for subjet b tagging
            rParam=options.fatJetRadius, # needed for subjet b tagging
            algo=algoLabel,              # has to be defined but is not used since svClustering=False
            runIVF = options.runIVF,
            postfix = postfix
        )
    else:
        _src = getattr(process,'ak4Jets').src if options.miniAOD else getattr(process,'pfJetsPFBRECO'+postfix).src
        if options.usePuppiForFatJets and not options.usePuppi:
            _src = cms.InputTag('puppi')
        _srcPVs = getattr(process,'ak4Jets').srcPVs if options.miniAOD else getattr(process,'pfJetsPFBRECO'+postfix).srcPVs
        ## Fat jets (Gen and Reco)
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        process.genFatJetsNoNu = ak4GenJets.clone(
            jetAlgorithm = cms.string(options.jetAlgo),
            rParam = cms.double(options.fatJetRadius),
            src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix))
        )
        from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
        process.fatPFJetsCHS = ak4PFJets.clone(
            jetAlgorithm = cms.string(options.jetAlgo),
            rParam = cms.double(options.fatJetRadius),
            src = _src,
            srcPVs = _srcPVs,
            doAreaFastjet = cms.bool(True),
            jetPtMin = cms.double(options.fatJetRawPtMin)
        )
        ## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
        from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
        process.genFatJetsNoNuPruned = ak4GenJets.clone(
            SubJetParameters,
            jetAlgorithm = cms.string(options.jetAlgo),
            rParam = cms.double(options.fatJetRadius),
            src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix)),
            usePruning = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
        )

        from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJetsCHSPruned
        process.fatPFJetsPruned = ak8PFJetsCHSPruned.clone(
            jetAlgorithm = cms.string(options.jetAlgo),
            rParam = cms.double(options.fatJetRadius),
            src = _src,
            srcPVs = _srcPVs,
            doAreaFastjet = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets"),
            jetPtMin = cms.double(options.fatJetRawPtMin)
        )
        ## SoftDrop fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
        process.genFatJetsNoNuSoftDrop = ak4GenJets.clone(
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
        from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJetsCHSSoftDrop
        process.fatPFJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone(
            jetAlgorithm = cms.string(options.jetAlgo),
            rParam = cms.double(options.fatJetRadius),
            R0 = cms.double(options.fatJetRadius),
            src = _src,
            srcPVs = _srcPVs,
            doAreaFastjet = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets"),
            jetPtMin = cms.double(options.fatJetRawPtMin)
        )

        ## PATify the above jets
        addJetCollection(
            process,
            labelName='FatPF',
            jetSource=cms.InputTag(fatJetSource),
            algo=algoLabel,              # needed for jet flavor clustering
            rParam=options.fatJetRadius, # needed for jet flavor clustering
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos = bTagInfosFat,
            btagDiscriminators = list(bTagDiscriminatorsFat),
            jetCorrections = jetCorrectionsAK8,
            genJetCollection = cms.InputTag(fatGenJetCollection),
            genParticles = cms.InputTag(genParticles),
            explicitJTA = options.useExplicitJTA,
            runIVF = options.runIVF,
            postfix = postfix
        )
        getattr(process,'selectedPatJetsFatPF'+postfix).cut = cms.string("pt > %f && abs(eta) < %f"%(float(options.fatJetPtMin), float(options.fatJetAbsEtaMax)))
        addJetCollection(
            process,
            labelName='SoftDropFatPF',
            jetSource=cms.InputTag(fatJetSourceSoftDrop),
            algo=algoLabel,
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos=['None'],
            btagDiscriminators=['None'],
            jetCorrections=jetCorrectionsAK8,
            genJetCollection = cms.InputTag(fatGenJetCollection),
            genParticles = cms.InputTag(genParticles),
            getJetMCFlavour = False, # jet flavor disabled
            postfix = postfix
        )
        addJetCollection(
            process,
            labelName='SoftDropSubjetsPF',
            jetSource=cms.InputTag(fatJetSourceSoftDrop,'SubJets'),
            algo=algoLabel,              # needed for subjet flavor clustering
            rParam=options.fatJetRadius, # needed for subjet flavor clustering
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos = bTagInfos_noDeepFlavour,
            btagDiscriminators = list(bTagDiscriminatorsSubJets),
            jetCorrections = jetCorrectionsSubJets,
            genJetCollection = cms.InputTag(fatGenJetCollectionSoftDrop,'SubJets'),
            genParticles = cms.InputTag(genParticles),
            explicitJTA = True,  # needed for subjet b tagging
            svClustering = True, # needed for subjet b tagging
            fatJets = cms.InputTag(fatJetSource),                # needed for subjet flavor clustering
            groomedFatJets = cms.InputTag(fatJetSourceSoftDrop), # needed for subjet flavor clustering
            runIVF = options.runIVF,
            postfix = postfix
        )

        ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
        process.selectedPatJetsSoftDropFatPFPacked = cms.EDProducer("BoostedJetMerger",
            jetSrc=cms.InputTag("selectedPatJetsSoftDropFatPF"+postfix),
            subjetSrc=cms.InputTag("selectedPatJetsSoftDropSubjetsPF"+postfix)
        )

        addJetCollection(
            process,
            labelName='PrunedFatPF',
            jetSource=cms.InputTag(fatJetSourcePruned),
            algo=algoLabel,
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos=['None'],
            btagDiscriminators=['None'],
            jetCorrections=jetCorrectionsAK8,
            genJetCollection = cms.InputTag(fatGenJetCollection),
            genParticles = cms.InputTag(genParticles),
            getJetMCFlavour = False, # jet flavor disabled
            postfix = postfix
        )
        addJetCollection(
            process,
            labelName='PrunedSubjetsPF',
            jetSource=cms.InputTag(fatJetSourcePruned,'SubJets'),
            algo=algoLabel,              # needed for subjet flavor clustering
            rParam=options.fatJetRadius, # needed for subjet flavor clustering
            pfCandidates = cms.InputTag(pfCandidates),
            pvSource = cms.InputTag(pvSource),
            svSource = cms.InputTag(svSource),
            muSource = cms.InputTag(muSource),
            elSource = cms.InputTag(elSource),
            btagInfos = bTagInfos_noDeepFlavour,
            btagDiscriminators = list(bTagDiscriminatorsSubJets),
            jetCorrections = jetCorrectionsSubJets,
            genJetCollection = cms.InputTag(fatGenJetCollectionPruned,'SubJets'),
            genParticles = cms.InputTag(genParticles),
            explicitJTA = True,  # needed for subjet b tagging
            svClustering = True, # needed for subjet b tagging
            fatJets = cms.InputTag(fatJetSource),              # needed for subjet flavor clustering
            groomedFatJets = cms.InputTag(fatJetSourcePruned), # needed for subjet flavor clustering
            runIVF = options.runIVF,
            postfix = postfix
        )

        ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
        process.selectedPatJetsPrunedFatPFPacked = cms.EDProducer("BoostedJetMerger",
            jetSrc=cms.InputTag("selectedPatJetsPrunedFatPF"+postfix),
            subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPF"+postfix)
        )

        ## Pack fat jets with subjets
        process.packedPatJetsFatPF = cms.EDProducer("JetSubstructurePacker",
                jetSrc = cms.InputTag('selectedPatJetsFatPF'+postfix),
                distMax = cms.double(options.fatJetRadius),
                algoTags = cms.VInputTag(),
                algoLabels = cms.vstring(),
                fixDaughters = cms.bool(False)
        )
        if options.useSoftDrop:
            process.packedPatJetsFatPF.algoTags.append( cms.InputTag('selectedPatJetsSoftDropFatPFPacked') )
            process.packedPatJetsFatPF.algoLabels.append( 'SoftDropPuppi' )
        if options.usePruned:
            process.packedPatJetsFatPF.algoTags.append( cms.InputTag('selectedPatJetsPrunedFatPFPacked') )
            process.packedPatJetsFatPF.algoLabels.append( 'Pruned' )

        #-------------------------------------
        ## N-subjettiness
        from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

        process.Njettiness = Njettiness.clone(
            src = cms.InputTag(fatJetSource),
            R0  = cms.double(options.fatJetRadius)
        )

        getattr(process,'patJetsFatPF'+postfix).userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']


        #-------------------------------------
        ## Grooming ValueMaps
        process.SoftDrop = cms.EDProducer("RecoJetDeltaRValueMapProducer",
            src = cms.InputTag(fatJetSource),
            matched = cms.InputTag("selectedPatJetsSoftDropFatPFPacked"),
            distMax = cms.double(options.fatJetRadius),
            values = cms.vstring('mass','pt','eta','phi','jecFactor(0)'),
            valueLabels = cms.vstring('Mass','Pt','Eta','Phi','jecFactor0'),
            lazyParser = cms.bool(True)
        )
        process.Pruned = cms.EDProducer("RecoJetDeltaRValueMapProducer",
            src = cms.InputTag(fatJetSource),
            matched = cms.InputTag("selectedPatJetsPrunedFatPFPacked"),
            distMax = cms.double(options.fatJetRadius),
            values = cms.vstring('mass','pt','eta','phi','jecFactor(0)'),
            valueLabels = cms.vstring('Mass','Pt','Eta','Phi','jecFactor0'),
            lazyParser = cms.bool(True)
        )

        getattr(process,'patJetsFatPF'+postfix).userData.userFloats.src += ['SoftDrop:Mass','SoftDrop:Pt','SoftDrop:Eta','SoftDrop:Phi','SoftDrop:jecFactor0',
                                                                               'Pruned:Mass'  ,'Pruned:Pt'  ,'Pruned:Eta'  ,'Pruned:Phi'  ,'Pruned:jecFactor0']


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

    if options.isReHLT and not options.runOnData:
        process.ttbarselectionproducer.triggerColl =  cms.InputTag("TriggerResults","","HLT2")

    #electron id
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

    if options.miniAOD:
        process.ttbarselectionproducer.electronColl = cms.InputTag('slimmedElectrons')
        process.ttbarselectionproducer.muonColl     = cms.InputTag('slimmedMuons')
        process.ttbarselectionproducer.jetColl      = cms.InputTag(patJetSource)
        process.ttbarselectionproducer.metColl      = cms.InputTag('slimmedMETs')
        switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    else:
        process.ttbarselectionproducer.electronColl = cms.InputTag('selectedPatElectrons'+postfix)
        process.ttbarselectionproducer.muonColl     = cms.InputTag('selectedPatMuons'+postfix)
        process.ttbarselectionproducer.jetColl      = cms.InputTag(patJetSource)
        process.ttbarselectionproducer.metColl      = cms.InputTag('patMETs'+postfix)
        switchOnVIDElectronIdProducer(process, DataFormat.AOD)

    # Set up electron ID (VID framework)
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    switchOnVIDElectronIdProducer(process, dataFormat=DataFormat.MiniAOD)
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



    #process.ttbarselectionproducer.isData       = options.runOnData
    #process.ttbarselectionproducer.electronColl = cms.InputTag('selectedPatElectrons'+postfix)
    #process.ttbarselectionproducer.muonColl     = cms.InputTag('selectedPatMuons'+postfix)
    #process.ttbarselectionproducer.jetColl      = cms.InputTag(patJetSource)
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
if options.miniAOD:
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#-------------------------------------
## Change the minimum number of tracker hits used in the track selection
if options.changeMinNumberOfHits:
    for m in process.producerNames().split(' '):
        if m.startswith('pfImpactParameterTagInfos'):
            print "Changing 'minimumNumberOfHits' for " + m + " to " + str(options.minNumberOfHits)
            getattr(process, m).minimumNumberOfHits = cms.int32(options.minNumberOfHits)

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag(pvSource))

#-------------------------------------
## Add TagInfos to PAT jets
for i in ['patJets', 'patJetsFatPF', 'patJetsSoftDropSubjetsPF', 'patJetsPrunedSubjetsPF',
          'updatedPatJetsTransientCorrected', 'updatedPatJetsTransientCorrectedFatPF', 'updatedPatJetsTransientCorrectedSoftDropSubjetsPF']:
    m = i + postfix
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )

#-------------------------------------
## Adapt fat jet b tagging
if options.runFatJets:
    getattr(process,'softPFElectronsTagInfosFatPF'+postfix).DeltaRElectronJet = cms.double(options.fatJetRadius) # default is 0.4
    if options.useLegacyTaggers:
        # Set the cone size for the jet-track association to the jet radius
        getattr(process,'jetTracksAssociatorAtVertexFatPF'+postfix).coneSize = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'secondaryVertexTagInfosFatPF'+postfix).trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius)   # default is 0.3
        getattr(process,'secondaryVertexTagInfosFatPF'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
        getattr(process,'inclusiveSecondaryVertexFinderTagInfosFatPF'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'inclusiveSecondaryVertexFinderTagInfosFatPF'+postfix).extSVDeltaRToJet = cms.double(options.fatJetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
        process.jetProbabilityComputerFat = process.jetProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.3
        getattr(process,'jetProbabilityBJetTagsFatPF'+postfix).jetTagComputer = cms.string('jetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
        process.jetBProbabilityComputerFat = process.jetBProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.4
        getattr(process,'jetBProbabilityBJetTagsFatPF'+postfix).jetTagComputer = cms.string('jetBProbabilityComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
        process.combinedSecondaryVertexV2ComputerFat = process.combinedSecondaryVertexV2Computer.clone()
        process.combinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        process.combinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        getattr(process,'combinedInclusiveSecondaryVertexV2BJetTagsFatPF'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexV2ComputerFat')
    else:
        # Set the cone size for the jet-track association to the jet radius
        getattr(process,'pfImpactParameterTagInfosFatPF'+postfix).maxDeltaR = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'pfSecondaryVertexTagInfosFatPF'+postfix).trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius)   # default is 0.3
        getattr(process,'pfSecondaryVertexTagInfosFatPF'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosFatPF'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.fatJetRadius) # default is 0.4
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosFatPF'+postfix).extSVDeltaRToJet = cms.double(options.fatJetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
        process.candidateJetProbabilityComputerFat = process.candidateJetProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.3
        getattr(process,'pfJetProbabilityBJetTagsFatPF'+postfix).jetTagComputer = cms.string('candidateJetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
        process.candidateJetBProbabilityComputerFat = process.candidateJetBProbabilityComputer.clone( deltaR = cms.double(options.fatJetRadius) ) # default is 0.4
        getattr(process,'pfJetBProbabilityBJetTagsFatPF'+postfix).jetTagComputer = cms.string('candidateJetBProbabilityComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
        process.candidateCombinedSecondaryVertexV2ComputerFat = process.candidateCombinedSecondaryVertexV2Computer.clone()
        process.candidateCombinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        process.candidateCombinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.fatJetRadius) # default is 0.3
        if hasattr(process,'pfCombinedInclusiveSecondaryVertexV2BJetTagsFatPF'+postfix):
            getattr(process,'pfCombinedInclusiveSecondaryVertexV2BJetTagsFatPF'+postfix).jetTagComputer = cms.string('candidateCombinedSecondaryVertexV2ComputerFat')

#-------------------------------------
process.btagana = bTagAnalyzer.clone()
if options.useLegacyTaggers:
    process.btagana = bTagAnalyzerLegacy.clone()
# The following combinations should be considered:
# For b-tagging performance measurements:
#   process.btagana.useSelectedTracks    = True
#   process.btagana.useTrackHistory      = False (or True for Mistag systematics with GEN-SIM-RECODEBUG samples)
#   options.produceJetTrackTree  = False
#   process.btagana.produceAllTrackTree  = False
#   process.btagana.producePtRelTemplate = False (or True for PtRel fit studies)
# or data/MC validation of jets, tracks and SVs:
#   process.btagana.useSelectedTracks    = False (or True for JP calibration)
#   process.btagana.useTrackHistory      = False
#   options.produceJetTrackTree  = True
#   process.btagana.produceAllTrackTree  = False
#   process.btagana.producePtRelTemplate = False
# or general tracks, PV and jet performance studies:
#   process.btagana.useSelectedTracks    = True
#   process.btagana.useTrackHistory      = False
#   options.produceJetTrackTree  = False
#   process.btagana.produceAllTrackTree  = True
#   process.btagana.producePtRelTemplate = False
#------------------

#Handle groups
for requiredGroup in process.btagana.groups:
   for storedGroup in btagana_tmp.groups:
     if (requiredGroup.group == storedGroup.group):
       requiredGroup.store = storedGroup.store

process.btagana.MaxEta                = options.maxJetEta ## for extended forward pixel coverage
process.btagana.MinPt                 = options.minJetPt
process.btagana.tracksColl            = cms.InputTag(trackSource)
process.btagana.useSelectedTracks     = options.useSelectedTracks ## False if you want to run on all tracks : for commissioning studies
process.btagana.useTrackHistory       = options.useTrackHistory ## Can only be used with GEN-SIM-RECODEBUG files
process.btagana.produceJetTrackTruthTree = options.useTrackHistory ## can only be used with GEN-SIM-RECODEBUG files and when useTrackHistory is True
process.btagana.produceAllTrackTree   = options.produceAllTrackTree ## True if you want to run info for all tracks : for commissioning studies
process.btagana.producePtRelTemplate  = options.producePtRelTemplate  ## True for performance studies
#------------------
process.btagana.runTagVariables     = options.runTagVariables  ## True if you want to run TagInfo TaggingVariables
process.btagana.runCSVTagVariables  = options.runCSVTagVariables   ## True if you want to run CSV TaggingVariables
process.btagana.runCSVTagTrackVariables  = options.runCSVTagTrackVariables   ## True if you want to run CSV Tagging Track Variables
process.btagana.runDeepFlavourTagVariables = options.runDeepFlavourTagVariables
process.btagana.primaryVertexColl     = cms.InputTag(pvSource)
process.btagana.Jets                  = cms.InputTag(patJetSource)
process.btagana.muonCollectionName    = cms.InputTag(muSource)
process.btagana.patMuonCollectionName = cms.InputTag(patMuons)
process.btagana.use_ttbar_filter      = cms.bool(options.useTTbarFilter)
#process.btagana.triggerTable          = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.triggerTable          = cms.InputTag(trigresults) # Data and MC
process.btagana.genParticles          = cms.InputTag(genParticles)
process.btagana.candidates            = cms.InputTag(pfCandidates)
process.btagana.runJetVariables     = options.runJetVariables
process.btagana.runQuarkVariables   = options.runQuarkVariables
process.btagana.runHadronVariables  = options.runHadronVariables
process.btagana.runGenVariables     = options.runGenVariables
process.btagana.runPFElectronVariables = options.runPFElectronVariables
process.btagana.runPFMuonVariables = options.runPFMuonVariables
process.btagana.runPatMuons = options.runPatMuons
process.btagana.runCTagVariables = options.runCTagVariables
process.btagana.runEventInfo = options.runEventInfo

process.btagana.runOnData = options.runOnData

if options.runOnData:
  process.btagana.runHadronVariables  = False
  process.btagana.runQuarkVariables   = False
  process.btagana.runGenVariables     = False

if options.runCTagVariables:
    process.btagana.runEventInfo = True

if not process.btagana.useTrackHistory  or not options.produceJetTrackTree:
    process.btagana.produceJetTrackTruthTree = False

if process.btagana.useTrackHistory:
    process.load('SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi')
    process.load('SimTracker.TrackerHitAssociation.tpClusterProducer_cfi')

if options.runFatJets:
    process.btaganaFatJets = process.btagana.clone(
        runEventInfo      = cms.bool(not options.processStdAK4Jets),
        allowJetSkipping    = cms.bool(False),
        runTagVariables     = options.runTagVariables,
        runDeepFlavourTagVariables = cms.bool(False),
        deepFlavourJetTags = cms.string(''),
        deepFlavourNegJetTags = cms.string(''),
        runTagVariablesSubJets = options.runTagVariablesSubJets,
        useSelectedTracks   = cms.bool(True),
        maxDeltaR           = cms.double(options.fatJetRadius),
        R0                  = cms.double(options.fatJetRadius),
        maxSVDeltaRToJet    = cms.double(options.fatJetRadius-(0.1+(options.fatJetRadius-0.8)*(0.4/0.7))), # linear interpolation from 0.7 at R=0.8 to 1.0 at R=1.5
        doubleSVBJetTags    = cms.string('pfBoostedDoubleSecondaryVertex' + ('CA15' if algoLabel=='CA' else 'AK8') + 'BJetTags'),
        distJetAxis         = cms.double(9999.),
        decayLength         = cms.double(9999.),
        deltaR              = cms.double(0.8),
        BranchNamePrefix    = cms.string('FatJetInfo'),
        Jets                = cms.InputTag(patFatJetSource),
        SubJets             = cms.VInputTag(),
        SubJetLabels        = cms.vstring(),
        runFatJets          = cms.bool(True),
        runSubJets          = options.runSubJets,
        svComputer          = cms.string('combinedSecondaryVertexV2ComputerFat' if options.useLegacyTaggers else 'candidateCombinedSecondaryVertexV2ComputerFat'),
        bdsvTagInfos        = cms.string('pfBoostedDoubleSV' + ('CA15' if algoLabel=='CA' else 'AK8')),
        use_ttbar_filter    = cms.bool(False)
    )
    if options.useSoftDrop:
        process.btaganaFatJets.SubJets.append( cms.InputTag(patSubJetSourceSoftDrop) )
        process.btaganaFatJets.SubJetLabels.append( 'SoftDropPuppi' )
    if options.usePruned:
        process.btaganaFatJets.SubJets.append( cms.InputTag('selectedPatJetsPrunedFatPFPacked:SubJets') )
        process.btaganaFatJets.SubJetLabels.append( 'Pruned' )

if options.doBoostedCommissioning:
    process.btaganaFatJets.runHadronVariables = True
    process.btaganaFatJets.runQuarkVariables = True
    process.btaganaFatJets.runPFMuonVariables = True
    process.btaganaFatJets.runCSVTagVariables = True
    process.btaganaFatJets.runCSVTagTrackVariables = True
    process.btaganaFatJets.runCSVTagVariablesSubJets = True
    print "**********NTuples will be made for boosted b tag commissioning. The following switches will be reset:**********"
    print "runHadronVariables set to '",process.btaganaFatJets.runHadronVariables,"'"
    print "runQuarkVariables set to '",process.btaganaFatJets.runQuarkVariables,"'"
    print "runPFMuonVariables set to '",process.btaganaFatJets.runPFMuonVariables,"'"
    print "For fat jets: runCSVTagVariables set to '",process.btaganaFatJets.runCSVTagVariables,"'"
    print "For fat jets: runCSVTagTrackVariables set to '",process.btaganaFatJets.runCSVTagTrackVariables,"'"
    print "For subjets:  runCSVTagVariablesSubJets set to '",process.btaganaFatJets.runCSVTagVariablesSubJets,"'"
    print "********************"

if process.btagana.produceJetTrackTruthTree:
    process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
    process.load("SimTracker.TrackHistory.TrackHistory_cff")
    process.load("SimTracker.TrackHistory.TrackClassifier_cff")
    process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
    process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

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

#Trick to make it work in 9_1_X
process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)

process.p = cms.Path(
    process.allEvents
    * process.filtSeq
    * process.selectedEvents
    * process.analyzerSeq,
    process.tsk
)

# Delete predefined output module (needed for running with CRAB)
del process.out
# dump content of cms.Process to python file
if options.dumpPython is not None:
    open('pydump.py','w').write(process.dumpPython())
# print-outs
print '--- runBTagAnalyzer_cfg.py ---\n'
print 'process.maxEvents.input =', process.maxEvents.input
print 'process.source.skipEvents =', process.source.skipEvents
print 'process.source.fileNames =', process.source.fileNames
print 'numThreads =', options.numThreads
print 'numStreams =', options.numStreams
print 'logs =', options.logs
print 'wantSummary =', options.wantSummary
print 'dumpPython =', options.dumpPython
print '\n-------------------------------'

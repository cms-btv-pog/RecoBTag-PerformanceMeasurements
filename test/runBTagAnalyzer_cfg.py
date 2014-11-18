
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
options.register('mcGlobalTag', 'PHYS14_25_V1',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'GR_R_70_V2',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
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
options.register('producePtRelTemplate', True,
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
options.register('useTopProjections', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use top projections"
)
options.register('miniAOD', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Running on miniAOD"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Running on miniAOD: %s"%('True' if options.miniAOD else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if not options.usePFchs:
    jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK4[1].append('L2L3Residual')
    jetCorrectionsAK7[1].append('L2L3Residual')

## b-tag infos
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','secondaryVertexNegativeTagInfos','softMuonTagInfos'
    #,'softPFMuonsTagInfos','softPFElectronsTagInfos','inclusiveSecondaryVertexFinderTagInfos'
    #,'inclusiveSecondaryVertexFinderFilteredTagInfos'
]
## b-tag discriminators
bTagDiscriminators = ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags',
    'negativeOnlyJetBProbabilityJetTags','negativeOnlyJetProbabilityJetTags','negativeTrackCountingHighEffJetTags',
    'negativeTrackCountingHighPurJetTags','positiveOnlyJetBProbabilityJetTags','positiveOnlyJetProbabilityJetTags',
    'simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','simpleSecondaryVertexNegativeHighEffBJetTags',
    'simpleSecondaryVertexNegativeHighPurBJetTags','combinedSecondaryVertexBJetTags','combinedSecondaryVertexPositiveBJetTags',
    'combinedSecondaryVertexNegativeBJetTags',
    #'softPFMuonBJetTags','positiveSoftPFMuonBJetTags','negativeSoftPFMuonBJetTags',
    #'softPFElectronBJetTags','positiveSoftPFElectronBJetTags','negativeSoftPFElectronBJetTags',
    #'simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags','doubleSecondaryVertexHighEffBJetTags',
    #'combinedInclusiveSecondaryVertexBJetTags','combinedInclusiveSecondaryVertexPositiveBJetTags',
    'combinedInclusiveSecondaryVertexV2BJetTags'
    #'combinedSecondaryVertexSoftPFLeptonV1BJetTags','positiveCombinedSecondaryVertexSoftPFLeptonV1BJetTags',
    #'negativeCombinedSecondaryVertexSoftPFLeptonV1BJetTags'
]
bTagDiscriminatorsSubJets = copy.deepcopy(bTagDiscriminators)
if 'doubleSecondaryVertexHighEffBJetTags' in bTagDiscriminators:
    bTagDiscriminatorsSubJets.remove('doubleSecondaryVertexHighEffBJetTags')

## Postfix
postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
jetSource = 'pfJets'+postfix
genJetCollection = 'ak4GenJetsNoNu'+postfix
trackSource = 'generalTracks'
pvSource = 'offlinePrimaryVertices'
svSource = cms.InputTag('inclusiveSecondaryVertices')
muons = 'muons'
selectedPatMuons = 'selectedPatMuons'
## If running on miniAOD
if options.miniAOD:
    genParticles = 'prunedGenParticles'
    jetSource = 'ak4PFJets'
    genJetCollection = 'ak4GenJetsNoNu'
    trackSource = 'unpackedTracksAndVertices'
    pvSource = 'unpackedTracksAndVertices'
    svSource = cms.InputTag('unpackedTracksAndVertices','secondary')
    muons = 'slimmedMuons'
    selectedPatMuons = muons

process = cms.Process("BTagAna")

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
        # /QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM
        #'/store/mc/Spring14dr/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00A2284F-D5D0-E311-BE96-002590A3A3D2.root'
        # /TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1/AODSIM
        #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/00120F7A-84F5-E311-9FBE-002618943910.root'
        # /RelValTTbar_13/CMSSW_7_2_0-PU50ns_PHYS14_25_V1_Phys14-v2/GEN-SIM-RECO
        '/store/relval/CMSSW_7_2_0/RelValTTbar_13/GEN-SIM-RECO/PU50ns_PHYS14_25_V1_Phys14-v2/00000/0E8D62AF-9059-E411-9066-0025905B85D0.root'
    )
)

if options.miniAOD:
    process.source.fileNames = [
        # /RelValTTbar_13/CMSSW_7_2_0-PU50ns_PHYS14_25_V1_Phys14-v2/MINIAODSIM
        '/store/relval/CMSSW_7_2_0/RelValTTbar_13/MINIAODSIM/PU50ns_PHYS14_25_V1_Phys14-v2/00000/206F2CD1-9C59-E411-B789-00261894388A.root'
    ]

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

##############################################
# Add GenParticlePruner for Boosted b-Tagging Studies
##############################################
process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag(genParticles),
    select = cms.vstring(
    "drop  *  ", #by default
    "keep ( status = 3 || (status>=21 && status<=29) )", #keep hard process particles
    "keep abs(pdgId) = 13 || abs(pdgId) = 15" #keep muons and taus
    )
)

########################################################
# Get calibrations for the CSVV1 and CSVSLV1 taggers
########################################################
#process.load("CondCore.DBCommon.CondDBSetup_cfi")
#process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
   #process.CondDBSetup,
   #timetype = cms.string('runnumber'),
   #toGet = cms.VPSet(cms.PSet(
      #record = cms.string('BTauGenericMVAJetTagComputerRcd'),
                #tag = cms.string('MVAComputerContainer_Retrained53X_JetTags_v2')
   #)),
   #connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
   #BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
#)
#process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord")

##############################################
# Get calibrations for the CSVV2 tagger
##############################################
process.load('CondCore.DBCommon.CondDBSetup_cfi')
process.BTauMVAJetTagComputerRecord = cms.ESSource('PoolDBESSource',
    process.CondDBSetup,
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('BTauGenericMVAJetTagComputerRcd'),
        tag = cms.string('MVAComputerContainer_53X_JetTags_v2')
    )),
    connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
)
process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer('PoolDBESSource','BTauMVAJetTagComputerRecord')


process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
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
if not options.miniAOD:
    ## PAT Configuration
    jetAlgo="AK4"

    from PhysicsTools.PatAlgos.tools.pfTools import *
    usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
	      jetCorrections=jetCorrectionsAK4, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

    ## Top projections in PF2PAT
    getattr(process,"pfPileUpJME"+postfix).checkClosestZVertex = False
    getattr(process,"pfNoPileUpJME"+postfix).enable = options.usePFchs
    if options.useTTbarFilter:
	getattr(process,"pfNoMuonJME"+postfix).enable = False
	getattr(process,"pfNoElectronJME"+postfix).enable = False
	getattr(process,"pfNoTau"+postfix).enable = False
	getattr(process,"pfNoJet"+postfix).enable = False
    else:
	getattr(process,"pfNoMuonJME"+postfix).enable = options.useTopProjections
	getattr(process,"pfNoElectronJME"+postfix).enable = options.useTopProjections
	getattr(process,"pfNoTau"+postfix).enable = False
	getattr(process,"pfNoJet"+postfix).enable = False
else:
    ## Recreate tracks and PVs for b tagging
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    ## Select isolated collections
    process.selectedMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
       (pfIsolationR04().sumChargedHadronPt+
	max(0.,pfIsolationR04().sumNeutralHadronEt+
	pfIsolationR04().sumPhotonEt-
	0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
	(isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>20. &&
	gsfTrack.isAvailable() &&
	gsfTrack.trackerExpectedHitsInner.numberOfLostHits<2 &&
	(pfIsolationVariables().sumChargedHadronPt+
	max(0.,pfIsolationVariables().sumNeutralHadronEt+
	pfIsolationVariables().sumPhotonEt-
	0.5*pfIsolationVariables().sumPUPt))/pt < 0.15'''))

    ## Do projections
    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
    process.pfNoMuonCHS =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfCHS"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectronsCHS = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuonCHS"), veto =  cms.InputTag("selectedElectrons"))

    process.pfNoMuon =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("packedPFCandidates"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectrons = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuon"), veto =  cms.InputTag("selectedElectrons"))

    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
    process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

    if options.useTTbarFilter:
        if options.usePFchs:
            process.ak4PFJets = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
        else:
            process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)
    else:
        if options.usePFchs:
            if options.useTopProjections:
                process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectronsCHS', doAreaFastjet = True)
            else:
                process.ak4PFJets = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
        else:
            if options.useTopProjections:
                process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectrons', doAreaFastjet = True)
            else:
                process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    trackSource = cms.InputTag(trackSource),
    pvSource = cms.InputTag(pvSource),
    svSource = svSource,
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag(genJetCollection),
    postfix = postfix
)
getattr(process,'patJetPartons'+postfix).particles = cms.InputTag(genParticles)
getattr(process,'patJetPartonMatch'+postfix).matched = cms.InputTag(genParticles)
getattr(process,'softMuonTagInfos'+postfix).leptons = cms.InputTag(muons)
process.inclusiveVertexFinder.tracks = cms.InputTag(trackSource)
process.trackVertexArbitrator.tracks = cms.InputTag(trackSource)

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
#-------------------------------------

#-------------------------------------
## CA8 jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNu = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix))
)

from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJets = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJets"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJets"+postfix).srcPVs),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(options.fatJetPtMin)
)

## CA8 pruned jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ca8GenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(0.8),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix)),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.ca8PFJetsPruned = ak4PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJets"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJets"+postfix).srcPVs),
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
        trackSource = cms.InputTag(trackSource),
        pvSource = cms.InputTag(pvSource),
        svSource = svSource,
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK7,
        genJetCollection = cms.InputTag('ca8GenJetsNoNu'),
        algo = 'CA',
        rParam = 0.8,
        postfix = postfix
    )
    getattr(process,'patJetPartonMatchCA8'+postfix).matched = cms.InputTag(genParticles)
    getattr(process,'softMuonTagInfosCA8'+postfix).leptons = cms.InputTag(muons)
    addJetCollection(
        process,
        labelName = 'CA8Pruned',
        jetSource = cms.InputTag('ca8PFJetsPruned'),
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK7,
        genJetCollection = cms.InputTag('ca8GenJetsNoNu'),
        getJetMCFlavour = False,
        postfix = postfix
    )
    getattr(process,'patJetPartonMatchCA8Pruned'+postfix).matched = cms.InputTag(genParticles)
    addJetCollection(
        process,
        labelName = 'CA8PrunedSubJets',
        jetSource = cms.InputTag('ca8PFJetsPruned','SubJets'),
        trackSource = cms.InputTag(trackSource),
        pvSource = cms.InputTag(pvSource),
        svSource = svSource,
        btagInfos=bTagInfos,
        btagDiscriminators=bTagDiscriminatorsSubJets,
        jetCorrections=jetCorrectionsAK4,
        genJetCollection = cms.InputTag('ca8GenJetsNoNuPruned','SubJets'),
        algo = 'CA',
        rParam = 0.8,
        postfix = postfix
    )
    getattr(process,'patJetPartonMatchCA8PrunedSubJets'+postfix).matched = cms.InputTag(genParticles)
    getattr(process,'softMuonTagInfosCA8PrunedSubJets'+postfix).leptons = cms.InputTag(muons)

    ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
    process.selectedPatJetsCA8PrunedPFlowPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsCA8Pruned"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsCA8PrunedSubJets"+postfix)
    )

    ## New jet flavor still requires some cfg-level adjustments for subjets until it is better integrated into PAT
    ## Adjust the jet flavor for pruned subjets
    setattr(process,'patJetFlavourAssociationCA8PrunedSubJets'+postfix, getattr(process,'patJetFlavourAssociationCA8'+postfix).clone(
        groomedJets = cms.InputTag('ca8PFJetsPruned'),
        subjets = cms.InputTag('ca8PFJetsPruned','SubJets')
    ))
    getattr(process,'patJetsCA8PrunedSubJets'+postfix).JetFlavourInfoSource = cms.InputTag('patJetFlavourAssociationCA8PrunedSubJets'+postfix,'SubJets')
#-------------------------------------

#-------------------------------------
## N-subjettiness
if options.runSubJets:
    from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

    process.Njettiness = Njettiness.clone(
        src = cms.InputTag("ca8PFJets"),
        cone = cms.double(0.8)
    )
    getattr(process,'patJetsCA8'+postfix).userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']
#-------------------------------------

#-------------------------------------
if options.runOnData and options.runSubJets:
    ## Remove MC matching when running over data
    removeMCMatching( process, ['All'] )

## Add TagInfos to PAT jets
patJets = ['patJets'+postfix]
if options.runSubJets:
    patJets += ['patJetsCA8'+postfix,'patJetsCA8PrunedSubJets'+postfix]

for m in patJets:
    if hasattr(process,m):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
        print "Switching 'addJetFlavourInfo' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addJetFlavourInfo', cms.bool(True) )
#-------------------------------------

#-------------------------------------
## Produce a collection of good primary vertices
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')
process.goodOfflinePrimaryVertices.src = cms.InputTag(pvSource)

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
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))
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
process.btagana.produceAllTrackTree   = False ## True if you want to keep info for all tracks : for commissioning studies
process.btagana.producePtRelTemplate  = options.producePtRelTemplate  ## True for performance studies
#------------------
process.btagana.storeTagVariables     = False ## True if you want to keep TagInfo TaggingVariables
process.btagana.storeCSVTagVariables  = True  ## True if you want to keep CSV TaggingVariables
process.btagana.primaryVertexColl     = cms.InputTag('goodOfflinePrimaryVertices')
process.btagana.Jets                  = cms.InputTag('selectedPatJets'+postfix)
process.btagana.tracksColl            = cms.InputTag(trackSource)
process.btagana.muonCollectionName    = cms.InputTag(muons)
process.btagana.patMuonCollectionName = cms.InputTag(selectedPatMuons)
process.btagana.use_ttbar_filter      = cms.bool(options.useTTbarFilter)
process.btagana.triggerTable          = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.genParticles          = cms.InputTag(genParticles)

if options.runSubJets:
    process.btaganaSubJets = process.btagana.clone(
        storeEventInfo      = cms.bool(False),
        produceJetTrackTree = cms.bool(True),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJetsCA8PrunedSubJets'+postfix),
        FatJets             = cms.InputTag('selectedPatJetsCA8'+postfix),
        GroomedFatJets      = cms.InputTag('selectedPatJetsCA8PrunedPFlowPacked'),
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
## Event counter
from MyAnalysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone()
process.selectedEvents = eventCounter.clone()
#---------------------------------------

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.JetHLTFilter*
    #process.noscraping
    process.primaryVertexFilter
    * process.goodOfflinePrimaryVertices
    #* process.HBHENoiseFilter
    #* process.CSCTightHaloFilter
    #* process.EcalDeadCellTriggerPrimitiveFilter
    #* process.eeBadScFilter
    #* process.ecalLaserCorrFilter
    #* process.trackingFailureFilter
    #* process.trkPOGFilters
)


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK4Jets:
    process.analyzerSeq += process.btagana
if options.runSubJets:
    process.analyzerSeq += process.btaganaSubJets
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

#open('pydump.py','w').write(process.dumpPython())

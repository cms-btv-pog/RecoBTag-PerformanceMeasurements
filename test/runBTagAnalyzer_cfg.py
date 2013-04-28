
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('outFilename', 'JetTree.root',
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
options.register('runSubJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = 'START53_V20'
if options.runOnData:
    globalTag = 'FT_53_V21_AN3'

## Jet energy corrections
jetCorrectionsAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
jetCorrectionsAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if not options.usePFchs:
    jetCorrectionsAK5 = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    jetCorrectionsAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if options.runOnData:
    jetCorrectionsAK5[1].append('L2L3Residual')
    jetCorrectionsAK7[1].append('L2L3Residual')

## b-tag discriminators
bTagDiscriminators = ['jetBProbabilityBJetTags', 'jetProbabilityBJetTags', 'trackCountingHighPurBJetTags', 'trackCountingHighEffBJetTags',
                      'simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags',
                      'simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags','doubleSecondaryVertexHighEffBJetTags',
                      'combinedInclusiveSecondaryVertexBJetTags','negativeOnlyJetBProbabilityJetTags',
                      'negativeOnlyJetProbabilityJetTags','negativeTrackCountingHighEffJetTags','negativeTrackCountingHighPurJetTags',
                      'simpleSecondaryVertexNegativeHighEffBJetTags','simpleSecondaryVertexNegativeHighPurBJetTags','combinedSecondaryVertexNegativeBJetTags',
                      'positiveOnlyJetBProbabilityJetTags','positiveOnlyJetProbabilityJetTags','combinedSecondaryVertexPositiveBJetTags','combinedInclusiveSecondaryVertexPositiveBJetTags']


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
        # /QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/GEN-SIM-RECODEBUG
        #'/store/mc/Summer12_DR53X/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/GEN-SIM-RECODEBUG/PU_S10_START53_V7A-v2/0000/041C8D66-05F4-E111-B16E-003048D43656.root'
        # /QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/GEN-SIM-RECODEBUG
        #'/store/mc/Summer12_DR53X/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/GEN-SIM-RECODEBUG/PU_S10_START53_V7A-v1/0000/0A4671D2-02F4-E111-9CD8-003048C69310.root'
        # /TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM
        '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v1/00000/FE7C71D8-DB25-E211-A93B-0025901D4C74.root'
    )
)

if options.runOnData:
    process.source.fileNames = [
        # /Jet/Run2012A-22Jan2013-v1/AOD
        #'/store/data/Run2012A/Jet/AOD/22Jan2013-v1/20000/30B21345-4172-E211-9EF3-00304867BEC0.root'
        # /BTag/Run2012A-22Jan2013-v1/AOD
        '/store/data/Run2012A/BTag/AOD/22Jan2013-v1/30000/E6959DA6-7081-E211-ABD3-002590596498.root'
    ]

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outFilename)
)

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag + '::All'


##############################################
# Get calibrations for the CSV taggers
##############################################
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
   process.CondDBSetup,
   timetype = cms.string('runnumber'),
   toGet = cms.VPSet(cms.PSet(
      record = cms.string('BTauGenericMVAJetTagComputerRcd'),
                tag = cms.string('MVAComputerContainer_Retrained53X_JetTags_v2')
   )),
   connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
   BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
)

process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord")


process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackHistory.TrackHistory_cff")
process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")


############################################################################
# For the Retrained CSV
############################################################################
process.combinedSecondaryVertexRetrained = process.combinedSecondaryVertex.clone(
  calibrationRecords = cms.vstring(
    'CombinedSVRetrainRecoVertex',
    'CombinedSVRetrainPseudoVertex',
    'CombinedSVRetrainNoVertex'
  )
)
process.combinedSecondaryVertexRetrainedBJetTags = process.combinedSecondaryVertexBJetTags.clone(
  jetTagComputer = cms.string('combinedSecondaryVertexRetrained'),
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPFlow"), cms.InputTag("secondaryVertexTagInfosAODPFlow"))
)

process.combinedSecondaryVertexRetrainedNegative = process.combinedSecondaryVertexNegative.clone(
  calibrationRecords = cms.vstring(
    'CombinedSVRetrainRecoVertex',
    'CombinedSVRetrainPseudoVertex',
    'CombinedSVRetrainNoVertex'
  )
)
process.combinedSecondaryVertexRetrainedNegativeBJetTags = process.combinedSecondaryVertexNegativeBJetTags.clone(
  jetTagComputer = cms.string('combinedSecondaryVertexRetrainedNegative'),
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPFlow"), cms.InputTag("secondaryVertexNegativeTagInfosAODPFlow"))
)

process.combinedSecondaryVertexRetrainedPositive = process.combinedSecondaryVertexPositive.clone(
  calibrationRecords = cms.vstring(
    'CombinedSVRetrainRecoVertex',
    'CombinedSVRetrainPseudoVertex',
    'CombinedSVRetrainNoVertex'
    )
  )
process.combinedSecondaryVertexRetrainedPositiveBJetTags = process.combinedSecondaryVertexPositiveBJetTags.clone(
  jetTagComputer = cms.string('combinedSecondaryVertexRetrainedPositive'),
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPFlow"), cms.InputTag("secondaryVertexTagInfosAODPFlow"))
)

############################################################################
# For the Soft Lepton taggers
############################################################################
process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonBJetTags_cfi")

process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonByPtBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonByPtBJetTags_cfi")

# for new softElectron tagger
process.softPFElectronRetrained = process.softElectron.clone()
process.softPFElectronRetrainedBJetTags = process.softPFElectronBJetTags.clone(
  jetTagComputer = 'softPFElectronRetrained',
  tagInfos = cms.VInputTag(cms.InputTag("softPFElectronsTagInfosAODPFlow"))
)

process.load("RecoBTag.SoftLepton.negativeSoftElectronES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftPFElectronBJetTags_cfi")
process.negativeSoftPFElectronRetrained = process.negativeSoftElectron.clone()
process.negativeSoftPFElectronRetrainedBJetTags = process.negativeSoftPFElectronBJetTags.clone(
  jetTagComputer = 'negativeSoftPFElectronRetrained',
  tagInfos = cms.VInputTag(cms.InputTag("softPFElectronsTagInfosAODPFlow"))
)

process.load("RecoBTag.SoftLepton.positiveSoftElectronES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftPFElectronBJetTags_cfi")
process.positiveSoftPFElectronRetrained = process.positiveSoftElectron.clone()
process.positiveSoftPFElectronRetrainedBJetTags = process.positiveSoftPFElectronBJetTags.clone(
  jetTagComputer = 'positiveSoftPFElectronRetrained',
  tagInfos = cms.VInputTag(cms.InputTag("softPFElectronsTagInfosAODPFlow"))
)

# for new softMuon tagger
process.softPFMuonRetrained = process.softMuon.clone()
process.softPFMuonRetrainedBJetTags = process.softPFMuonBJetTags.clone(
  jetTagComputer = 'softPFMuonRetrained',
  tagInfos = cms.VInputTag(cms.InputTag("softPFMuonsTagInfosAODPFlow"))
)

process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftPFMuonBJetTags_cfi")
process.negativeSoftPFMuonRetrained = process.negativeSoftMuon.clone()
process.negativeSoftPFMuonRetrainedBJetTags = process.negativeSoftPFMuonBJetTags.clone(
  jetTagComputer = 'negativeSoftPFMuonRetrained',
  tagInfos = cms.VInputTag(cms.InputTag("softPFMuonsTagInfosAODPFlow"))
)

process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftPFMuonBJetTags_cfi")
process.positiveSoftPFMuonRetrained = process.positiveSoftMuon.clone()
process.positiveSoftPFMuonRetrainedBJetTags = process.positiveSoftPFMuonBJetTags.clone(
  jetTagComputer = 'positiveSoftPFMuonRetrained',
  tagInfos = cms.VInputTag(cms.InputTag("softPFMuonsTagInfosAODPFlow"))
)

############################################################################
# For the combined CSV + JP tagger
############################################################################
process.combinedCSVJP = process.combinedMVA.clone(
  calibrationRecord = 'CombinedCSVJP',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer = cms.string('jetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer = cms.string('softPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('softPFElectronRetrained')
    )
  )
)
process.combinedCSVJPBJetTags = process.combinedMVABJetTags.clone(
  jetTagComputer = 'combinedCSVJP',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

process.load("RecoBTau.JetTagComputer.negativeCombinedMVAES_cfi")
process.negativeCombinedCSVJP = process.negativeCombinedMVA.clone(
  calibrationRecord = 'CombinedCSVJP',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeOnlyJetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrainedNegative')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeSoftPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeSoftPFElectronRetrained')
    )
  )
)
process.load("RecoBTau.JetTagComputer.negativeCombinedMVABJetTags_cfi")
process.negativeCombinedCSVJPBJetTags = process.negativeCombinedMVABJetTags.clone(
  jetTagComputer = 'negativeCombinedCSVJP',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

process.load("RecoBTau.JetTagComputer.positiveCombinedMVAES_cfi")
process.positiveCombinedCSVJP = process.positiveCombinedMVA.clone(
  calibrationRecord = 'CombinedCSVJP',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveOnlyJetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrainedPositive')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveSoftPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveSoftPFElectronRetrained')
    )
  )
)
process.load("RecoBTau.JetTagComputer.positiveCombinedMVABJetTags_cfi")
process.positiveCombinedCSVJPBJetTags = process.positiveCombinedMVABJetTags.clone(
  jetTagComputer = 'positiveCombinedCSVJP',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)


############################################################################
# For the combined CSV + JP + SL tagger
############################################################################
process.combinedCSVJPSL = process.combinedMVA.clone(
  calibrationRecord = 'CombinedCSVJPSL',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer = cms.string('jetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer = cms.string('softPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('softPFElectronRetrained')
    )
  )
)
process.combinedCSVJPSLBJetTags = process.combinedMVABJetTags.clone(
  jetTagComputer = 'combinedCSVJPSL',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

process.negativeCombinedCSVJPSL = process.negativeCombinedMVA.clone(
  calibrationRecord = 'CombinedCSVJPSL',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeOnlyJetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrainedNegative')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeSoftPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeSoftPFElectronRetrained')
    )
  )
)
process.negativeCombinedCSVJPSLBJetTags = process.negativeCombinedMVABJetTags.clone(
  jetTagComputer = 'negativeCombinedCSVJPSL',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

process.positiveCombinedCSVJPSL = process.positiveCombinedMVA.clone(
  calibrationRecord = 'CombinedCSVJPSL',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveOnlyJetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrainedPositive')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveSoftPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveSoftPFElectronRetrained')
    )
  )
)
process.positiveCombinedCSVJPSLBJetTags = process.positiveCombinedMVABJetTags.clone(
  jetTagComputer = 'positiveCombinedCSVJPSL',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

############################################################################
# For the combined CSV + SL tagger
############################################################################
process.combinedCSVSL = process.combinedMVA.clone(
  calibrationRecord = 'CombinedCSVSL',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer = cms.string('jetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer = cms.string('softPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('softPFElectronRetrained')
    )
  )
)
process.combinedCSVSLBJetTags = process.combinedMVABJetTags.clone(
  jetTagComputer = 'combinedCSVSL',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

process.negativeCombinedCSVSL = process.negativeCombinedMVA.clone(
  calibrationRecord = 'CombinedCSVSL',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeOnlyJetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrainedNegative')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeSoftPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('negativeSoftPFElectronRetrained')
    )
  )
)
process.negativeCombinedCSVSLBJetTags = process.negativeCombinedMVABJetTags.clone(
  jetTagComputer = 'negativeCombinedCSVSL',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)

process.positiveCombinedCSVSL = process.positiveCombinedMVA.clone(
  calibrationRecord = 'CombinedCSVSL',
  jetTagComputers = cms.VPSet(
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveOnlyJetProbability')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('combinedSecondaryVertexRetrainedPositive')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveSoftPFMuonRetrained')
    ),
    cms.PSet(
      discriminator = cms.bool(True),
      variables = cms.bool(False),
      jetTagComputer =
      cms.string('positiveSoftPFElectronRetrained')
    )
  )
)
process.positiveCombinedCSVSLBJetTags = process.positiveCombinedMVABJetTags.clone(
  jetTagComputer = 'positiveCombinedCSVSL',
  tagInfos = cms.VInputTag(
    cms.InputTag("impactParameterTagInfosAODPFlow"),
    cms.InputTag("secondaryVertexTagInfosAODPFlow"),
    cms.InputTag("softPFMuonsTagInfosAODPFlow"),
    cms.InputTag("softPFElectronsTagInfosAODPFlow")
  )
)
#-------------------------------------

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
## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

postfix = "PFlow"
jetAlgo="AK5"

from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
          jetCorrections=jetCorrectionsAK5, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfNoPileUp"+postfix).enable = options.usePFchs
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above b-tag discriminators)
switchJetCollection(process,
    cms.InputTag('pfNoTau'+postfix),
    doJTA        = True,
    doBTagging   = True,
    #btagInfo=bTagInfos,
    btagdiscriminators = bTagDiscriminators,
    jetCorrLabel = jetCorrectionsAK5,
    doType1MET   = False,
    genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
    doJetID      = False,
    postfix      = postfix
)

## Add additional b-taggers to the PAT sequence
getattr(process,'patPF2PATSequence'+postfix).replace( getattr(process,'combinedInclusiveSecondaryVertexPositiveBJetTagsAOD'+postfix),
    getattr(process,'combinedInclusiveSecondaryVertexPositiveBJetTagsAOD'+postfix)
    + process.softPFMuonRetrainedBJetTags + process.negativeSoftPFMuonRetrainedBJetTags + process.positiveSoftPFMuonRetrainedBJetTags
    + process.softPFElectronRetrainedBJetTags + process.negativeSoftPFElectronRetrainedBJetTags + process.positiveSoftPFElectronRetrainedBJetTags
    + process.combinedSecondaryVertexRetrainedBJetTags + process.combinedSecondaryVertexRetrainedNegativeBJetTags + process.combinedSecondaryVertexRetrainedPositiveBJetTags
    + process.combinedCSVJPBJetTags + process.negativeCombinedCSVJPBJetTags + process.positiveCombinedCSVJPBJetTags
    + process.combinedCSVSLBJetTags + process.negativeCombinedCSVSLBJetTags + process.positiveCombinedCSVSLBJetTags
    + process.combinedCSVJPSLBJetTags + process.negativeCombinedCSVJPSLBJetTags + process.positiveCombinedCSVJPSLBJetTags
)

newDiscriminatorSources = cms.VInputTag(
    cms.InputTag("softPFMuonRetrainedBJetTags"),
    cms.InputTag("negativeSoftPFMuonRetrainedBJetTags"),
    cms.InputTag("positiveSoftPFMuonRetrainedBJetTags"),
    cms.InputTag("softPFElectronRetrainedBJetTags"),
    cms.InputTag("negativeSoftPFElectronRetrainedBJetTags"),
    cms.InputTag("positiveSoftPFElectronRetrainedBJetTags"),
    cms.InputTag("combinedSecondaryVertexRetrainedBJetTags"),
    cms.InputTag("combinedSecondaryVertexRetrainedNegativeBJetTags"),
    cms.InputTag("combinedSecondaryVertexRetrainedPositiveBJetTags"),
    cms.InputTag("combinedCSVJPBJetTags"),
    cms.InputTag("negativeCombinedCSVJPBJetTags"),
    cms.InputTag("positiveCombinedCSVJPBJetTags"),
    cms.InputTag("combinedCSVSLBJetTags"),
    cms.InputTag("negativeCombinedCSVSLBJetTags"),
    cms.InputTag("positiveCombinedCSVSLBJetTags"),
    cms.InputTag("combinedCSVJPSLBJetTags"),
    cms.InputTag("negativeCombinedCSVJPSLBJetTags"),
    cms.InputTag("positiveCombinedCSVJPSLBJetTags")
)

## Add additional b-tag discriminators to the default jet collection
getattr(process,'patJets'+postfix).discriminatorSources = getattr(process,'patJets'+postfix).discriminatorSources + newDiscriminatorSources
#-------------------------------------

#-------------------------------------
# CA8 jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNu = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJets = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True)
)

## CA8 pruned jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ca8GenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
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
    jetCollInstanceName=cms.string("SubJets")
)

if options.runSubJets:
    ## PATify CA8 jets
    addJetCollection(
        process,
        cms.InputTag('ca8PFJets'),
        'CA8','PF',
        doJTA=True,
        doBTagging=True,
        #btagInfo=bTagInfos,
        btagdiscriminators=bTagDiscriminators,
        jetCorrLabel=jetCorrectionsAK7,
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False,
        genJetCollection=cms.InputTag("ca8GenJetsNoNu")
    )
    addJetCollection(
        process,
        cms.InputTag('ca8PFJetsPruned'),
        'CA8Pruned','PF',
        doJTA=False,
        doBTagging=False,
        #btagInfo=bTagInfos,
        btagdiscriminators=bTagDiscriminators,
        jetCorrLabel=jetCorrectionsAK7,
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False,
        genJetCollection=cms.InputTag("ca8GenJetsNoNu")
    )
    addJetCollection(
        process,
        cms.InputTag('ca8PFJetsPruned','SubJets'),
        'CA8PrunedSubJets', 'PF',
        doJTA=True,
        doBTagging=True,
        #btagInfo=bTagInfos,
        btagdiscriminators=bTagDiscriminators,
        jetCorrLabel=jetCorrectionsAK5,
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False,
        genJetCollection=cms.InputTag('ca8GenJetsNoNuPruned','SubJets')
    )

## Establish references between PAT fat jets and PAT subjets using the BoostedJetMerger
process.selectedPatJetsCA8PrunedPFPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8PrunedPF"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8PrunedSubJetsPF")
)
#-------------------------------------

#-------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *
## Remove objects not used from the PAT sequences to speed up processing
removeSpecificPATObjects(process,names=['Photons', 'Electrons', 'Muons', 'Taus'],postfix=postfix)
if options.runSubJets:
    removeAllPATObjectsBut(process, ['Jets'])

if options.runOnData and options.runSubJets:
    ## Remove MC matching when running over data
    removeMCMatching( process, ['All'] )

## Add TagInfos to PAT jets
patJets = ['patJets'+postfix]
if options.runSubJets:
    patJets.append('patJetsCA8PrunedSubJetsPF')

for m in patJets:
    if hasattr(process,m):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
#-------------------------------------

#-------------------------------------
## Produce a collection of good primary vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
)

## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)
#-------------------------------------

#-------------------------------------
from RecoBTag.PerformanceMeasurements.patTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix=postfix, sequence='patPF2PATSequence')
if options.runSubJets:
    adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='patDefaultSequence')
#-------------------------------------

#-------------------------------------
if not options.runOnData:
    ## JP calibration for cmsRun only :
    #from CondCore.DBCommon.CondDBCommon_cfi import *
    #process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
    #process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_Data_2011_41X.db")
    #process.es_prefer_trackProbabilityFakeCond = cms.ESPrefer("PoolDBESSource","trackProbabilityFakeCond")

    ## JP calibration for crab only:
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

process.btagana.use_selected_tracks  = True   ## False if you want to run on all tracks
process.btagana.useTrackHistory      = False
process.btagana.produceJetProbaTree  = False  ## True if you want to keep track and SV info!
process.btagana.producePtRelTemplate = True  ## True for performance studies
if options.runOnData:
    process.btagana.use_selected_tracks  = False   ## False if you want to run on all tracks
    process.btagana.useTrackHistory      = False
    process.btagana.produceJetProbaTree  = True  ## True if you want to keep track and SV info!
    process.btagana.producePtRelTemplate = False  ## True for performance studies
process.btagana.primaryVertexColl = cms.InputTag('goodOfflinePrimaryVertices')
process.btagana.Jets = cms.InputTag('selectedPatJets'+postfix)
process.btagana.triggerTable = cms.InputTag('TriggerResults::HLT') # Data and MC

process.btaganaSubJets = process.btagana.clone( Jets = cms.InputTag('selectedPatJetsCA8PrunedSubJetsPF') )
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
## Filter for HCAL laser events in prompt 2012A+B+C, snippet for "Datasets from the 2013 rereco and Multijet parked":
## https://twiki.cern.ch/twiki/bin/view/CMS/PdmVKnowFeatures#HCAL_laser_events_in_prompt_2012
process.load("EventFilter.HcalRawToDigi.hcallaserFilterFromTriggerResult_cff")
#---------------------------------------

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.JetHLTFilter*
    process.noscraping
    *process.primaryVertexFilter
)
if options.runOnData:
    process.filtSeq = cms.Sequence( process.filtSeq * process.hcalfilter )


## Define jet sequences
process.genJetSeq = cms.Sequence(
    process.ca8GenJetsNoNu
    + process.ca8GenJetsNoNuPruned
)
process.jetSeq = cms.Sequence(
    process.ca8PFJets
    + process.ca8PFJetsPruned
)
if not options.runOnData:
    process.jetSeq = cms.Sequence( process.genJetSeq + process.jetSeq )


## Define combined PF2PAT + subjet sequence
process.combPF2PATSubJetSeq = cms.Sequence( getattr(process,"patPF2PATSequence"+postfix) )
if options.runSubJets:
    process.combPF2PATSubJetSeq = cms.Sequence(
        getattr(process,"patPF2PATSequence"+postfix)
        * process.jetSeq
        * getattr(process,"patDefaultSequence")
        * process.selectedPatJetsCA8PrunedPFPacked
    )

## Define analyzer sequence
process.analyzerSeq = cms.Sequence( process.btagana )
if options.runSubJets:
    process.analyzerSeq = cms.Sequence( process.btagana + process.btaganaSubJets )
#---------------------------------------

process.p = cms.Path(
    process.filtSeq
    * process.goodOfflinePrimaryVertices
    * process.combPF2PATSubJetSeq
    * process.analyzerSeq
)

# Delete predefined output module (needed for running with CRAB)
del process.out

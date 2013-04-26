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
    input = cms.untracked.int32(100)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "FT_53_V21_AN3::All"
#process.GlobalTag.globaltag = "GR_R_53_V19::All"
# process.GlobalTag.globaltag = "GR_R_53_V2B::All"




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


## Output Module Configuration
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('testtt.root')
)

#-------------------------------------
#Produce PFJets from PF2PAT

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

runOnMC = False
usePFnoPU=True
jetCorrectionsAK5=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
jetCorrectionsAK7=('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK5"

usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix,
          jetCorrections=jetCorrectionsAK5, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfNoPileUp"+postfix).enable = usePFnoPU
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

## b-tag discriminators
bTagDiscriminators = ['jetBProbabilityBJetTags', 'jetProbabilityBJetTags', 'trackCountingHighPurBJetTags', 'trackCountingHighEffBJetTags',
                      'simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags',
                      'simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags','doubleSecondaryVertexHighEffBJetTags',
                      'combinedInclusiveSecondaryVertexBJetTags','negativeOnlyJetBProbabilityJetTags',
                      'negativeOnlyJetProbabilityJetTags','negativeTrackCountingHighEffJetTags','negativeTrackCountingHighPurJetTags',
                      'simpleSecondaryVertexNegativeHighEffBJetTags','simpleSecondaryVertexNegativeHighPurBJetTags','combinedSecondaryVertexNegativeBJetTags',
                      'positiveOnlyJetBProbabilityJetTags','positiveOnlyJetProbabilityJetTags','combinedSecondaryVertexPositiveBJetTags','combinedInclusiveSecondaryVertexPositiveBJetTags']

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above b-tag discriminators)
switchJetCollection(process,
    cms.InputTag("pfNoTauPFlow"),
    doJTA        = True,
    doBTagging   = True,
    btagdiscriminators = bTagDiscriminators,
    jetCorrLabel = jetCorrectionsAK5,
    doType1MET   = False,
    genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
    doJetID      = False,
    postfix      = postfix
)

## Add additional b-taggers to the PAT sequence
getattr(process,'patPF2PATSequencePFlow').replace( getattr(process,'combinedInclusiveSecondaryVertexPositiveBJetTagsAODPFlow'),
  getattr(process,'combinedInclusiveSecondaryVertexPositiveBJetTagsAODPFlow')
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
getattr(process,'patJetsPFlow').discriminatorSources = getattr(process,'patJetsPFlow').discriminatorSources + newDiscriminatorSources


## Add TagInfos to PAT jets
for m in ['patJets'+postfix]:
    if hasattr(process,m):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )


from PhysicsTools.PatAlgos.tools.coreTools import *
## Remove taus from the PAT sequence to speed up processing
removeSpecificPATObjects(process,names=['Taus'],postfix=postfix)

#-------------------------------------
#Redo the primary vertex

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

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


from RecoBTag.PerformanceMeasurements.patTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix=postfix, sequence='patPF2PATSequence')


### JP calibration for cmsRun only : 
# from CondCore.DBCommon.CondDBCommon_cfi import *
# process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
# process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_Data_2011_41X.db")
# process.es_prefer_trackProbabilityFakeCond = cms.ESPrefer("PoolDBESSource","trackProbabilityFakeCond")

### JP calibration for crab only: 
## process.GlobalTag.toGet = cms.VPSet(
##   cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
##        tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
##        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
##   cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
##        tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
##        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
## )


#---------------------------------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree_data.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )

#
process.load("RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff")

#$$
process.btagana.use_selected_tracks = False   ## False if you want to run on all tracks
process.btagana.useTrackHistory     = False
process.btagana.produceJetProbaTree = True  ## True if you want to keep track and SV info!
process.btagana.producePtRelTemplate = False  ## True for performance studies
process.btagana.primaryVertexColl = 'goodOfflinePrimaryVertices'
process.btagana.Jets = 'selectedPatJets'+postfix
process.btagana.triggerTable = 'TriggerResults::HLT' # Data and MC
#---------------------------------------

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

# Filter for HCAL laser events in prompt 2012A+B+C, snippet for "Datasets from the 2013 rereco and Multijet parked":
# https://twiki.cern.ch/twiki/bin/view/CMS/PdmVKnowFeatures#HCAL_laser_events_in_prompt_2012
process.load("EventFilter.HcalRawToDigi.hcallaserFilterFromTriggerResult_cff")


process.p = cms.Path(
#$$
#$$        process.JetHLTFilter
#$$
        process.noscraping
        *process.hcalfilter # filter for HCAL laser events
        *process.primaryVertexFilter
	*process.goodOfflinePrimaryVertices
	*getattr(process,"patPF2PATSequence"+postfix)
	*process.btagana
	)
	
	## Output Module Configuration (expects a path 'p')





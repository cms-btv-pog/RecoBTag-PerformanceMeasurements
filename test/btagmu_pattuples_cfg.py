# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('useData',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run this on real data")

options.register ('applyHLT',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Apply HLT Filter")

options.register ('hltProcess',
                  'HLT_BTagMu*',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "HLT process name to use.")

options.register ('use41x',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Use the 41x options")

options.register ('writeS8Tree',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run S8 tree maker")

options.parseArguments()

inputJetCorrLabelnoPF2PAT = ('AK5PF',['L1FastJet', 'L2Relative', 'L3Absolute'])

if not options.useData :
    
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
        
    process.source.fileNames = [
        '/store/relval/CMSSW_4_2_2/RelValTTbar/GEN-SIM-RECO/START42_V11-v1/0005/50AC4DBF-746D-E011-8CF9-00248C55CC62.root'
        ]
    
else :
    
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
        
    process.source.fileNames = [
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/02C65CCC-8F7E-E011-AAEF-00304867BFA8.root',
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/02A210E3-937E-E011-8651-0026189438F3.root',
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/029270DD-E67C-E011-933B-002618943838.root',
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/0286BD59-577C-E011-A09E-0026189438A7.root',
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/006930CD-A17E-E011-B660-002618943833.root',
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/002C519F-967E-E011-B128-0026189438BF.root',
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/001D4EEF-C67E-E011-B46D-00304867904E.root'
        
        ]

print options

print 'Running jet corrections: '
print inputJetCorrLabel

import sys


###############################
####### Global Setup ##########
###############################


# 4.2.x configuration
if not options.use41x :
    if options.useData :
        process.GlobalTag.globaltag = cms.string( 'GR_R_42_V12::All' )
    else :
        process.GlobalTag.globaltag = cms.string( 'START42_V12::All' )
else:
    # 4.1.x configuration
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    process.jec = cms.ESSource("PoolDBESSource",
                               DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
                               timetype = cms.string('runnumber'),
                               toGet = cms.VPSet(
        cms.PSet(
        record = cms.string('JetCorrectionsRecord'),
        tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5PFchs'),
        label  = cms.untracked.string('AK5PFchs')
        )
        ),
                               ## here you add as many jet types as you need (AK5Calo, AK5JPT, AK7PF, AK7Calo, KT4PF, KT4Calo, KT6PF, KT6Calo)
                               connect = cms.string('sqlite_file:Jec10V3.db')
                               )
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
    

    
# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )
# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
# Modify defaults setting to avoid an over-efficiency in the presence of OFT PU
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)

###############################
####### HLT Filter ############
process.HLTfilter = cms.EDFilter("HLTHighLevel",
                                 TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
                                 HLTPaths           = cms.vstring("HLT_BTagMu*"),
                                 #    HLTPaths           = cms.vstring("HLT_Jet*"),
                                 eventSetupPathsKey = cms.string(''),
                                 andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
                                 throw              = cms.bool(True)
                                 )

# switch on PAT trigger
#from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
#switchOnTrigger( process, hltProcess=options.hltProcess )


###############################
####### DAF PV's     ##########
###############################

pvSrc = 'offlinePrimaryVertices'
if options.use41x :
    pvSrc = 'offlinePrimaryVerticesDAF'
    process.offlinePrimaryVerticesDAF = cms.EDProducer("PrimaryVertexProducer",
                                                       verbose = cms.untracked.bool(False),
                                                       algorithm = cms.string('AdaptiveVertexFitter'),
                                                       TrackLabel = cms.InputTag("generalTracks"),
                                                       useBeamConstraint = cms.bool(False),
                                                       beamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                                       minNdof  = cms.double(0.0),
                                                       PVSelParameters = cms.PSet(
        maxDistanceToBeam = cms.double(1.0)
        ),
                                                       TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        maxNormalizedChi2 = cms.double(20.0),
        minPixelLayersWithHits=cms.int32(2),
        minSiliconLayersWithHits = cms.int32(5),
        maxD0Significance = cms.double(5.0),
        minPt = cms.double(0.0),
        trackQuality = cms.string("any")
        ),
                                                       
                                                       TkClusParameters = cms.PSet(
        algorithm   = cms.string("DA"),
        TkDAClusParameters = cms.PSet(
        coolingFactor = cms.double(0.6),  #  moderate annealing speed
        Tmin = cms.double(4.),            #  end of annealing
        vertexSize = cms.double(0.01),    #  ~ resolution / sqrt(Tmin)
        d0CutOff = cms.double(3.),        # downweight high IP tracks
        dzCutOff = cms.double(4.)         # outlier rejection after freeze-out (T<Tmin)
        )
        )
)
    
    
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag(pvSrc),
                                           minimumNDOF = cms.uint32(7) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( maxZ = cms.double(24.0),
                                     minNdof = cms.double(7.0)
                                     ),
    src=cms.InputTag(pvSrc)
    )


###############################
########## Gen Setup ##########
###############################

# prune gen particles
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
                                            src = cms.InputTag("genParticles"),
                                            select = cms.vstring(
                                                "drop  *",
                                                "keep status = 3", #keeps all particles from the hard matrix element
                                                "+keep (abs(pdgId) = 11 | abs(pdgId) = 13) & status = 1" #keeps all stable muons and electrons and their (direct) mothers.
                                                )
                                            )


###############################
########## PF Setup ###########
###############################

# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process,
          runPF2PAT=True,
          jetAlgo='AK5',
          runOnMC=not options.useData,
          postfix=postfix,
          jetCorrections=inputJetCorrLabel)

# Do not remove muons inside jets
getattr(process,"pfNoMuon"+postfix).enable = False
#getattr(process,"pfNoElectron"+postfix).enable = False

process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.Vertices = 'goodOfflinePrimaryVertices'
process.pfPileUpPFlow.checkClosestZVertex = False

process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet = False
process.patJetCorrFactorsPFlow.payload = inputJetCorrLabel[0]
process.patJetCorrFactorsPFlow.levels = inputJetCorrLabel[1]
process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow", "rho")

if not options.use41x : #and not options.forceCheckClosestZVertex :
    process.pfPileUpPFlow.checkClosestZVertex = False



###############################
###### Electron ID ############
###############################

# NOTE: ADDING THE ELECTRON IDs FROM CiC ----- USED WITH 42X 
    

###############################
###### Bare KT 0.6 jets #######
###############################

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
    rParam = cms.double(0.6),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    voronoiRfact = cms.double(0.9)
    )

process.kt6PFJets = kt4PFJets.clone(
    rParam = cms.double(0.6),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    voronoiRfact = cms.double(0.9)
    )


for ipostfix in [postfix] :
    for module in (
        getattr(process,"kt6PFJets"),
        getattr(process,"kt6PFJets" + ipostfix)
        ) :
        getattr(process,"patPF2PATSequence"+ipostfix).replace( getattr(process,"pfNoElectron"+ipostfix), getattr(process,"pfNoElectron"+ipostfix)*module )


# Use the good primary vertices everywhere. 
for imod in [process.patMuonsPFlow,
             process.patElectronsPFlow,
             process.patMuons,
             process.patElectrons] :
    imod.pvSrc = "goodOfflinePrimaryVertices"
    imod.embedTrack = True


addJetCollection(process,
                 cms.InputTag('ak5PFJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
                 'AK5', 'PF',
                 doJTA=True,            # Run Jet-Track association & JetCharge
                 doBTagging=True,       # Run b-tagging
                 jetCorrLabel=inputJetCorrLabelnoPF2PAT,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID = False
                 )
process.patJetCorrFactorsAK5PF.rho = cms.InputTag("kt6PFJetsPFlow", "rho")

###############################
### TagInfo and Matching Setup#
###############################

for jetcoll in (process.patJets,
                process.patJetsPFlow,
                process.patJetsAK5PF
                ) :
    if options.useData == False :
        jetcoll.embedGenJetMatch = True
        jetcoll.getJetMCFlavour = True
        jetcoll.addGenPartonMatch = True
    # Add CATopTag info... piggy-backing on b-tag functionality
    jetcoll.addBTagInfo = True
    # Add the calo towers and PFCandidates.
    # I'm being a little tricksy here, because I only
    # actually keep the products if the "writeFat" switch
    # is on. However, this allows for overlap checking
    # with the Refs so satisfies most use cases without
    # having to add to the object size
    #jetcoll.embedCaloTowers = True
    #jetcoll.embedPFCandidates = True

process.patJetsPFlow.addBTagInfo = True
 
###############################
#### Selections Setup #########
###############################

# AK5 Jets
process.selectedPatJetsPFlow.cut = cms.string("pt > 20 & abs(eta) < 2.4")
process.patJetsPFlow.addTagInfos = True
process.selectedPatJetsAK5PF.cut = cms.string("pt > 20 & abs(eta) < 2.4")

#process.patJetsPFlow.tagInfoSources = cms.VInputTag(
#    cms.InputTag("secondaryVertexTagInfosAODPFlow")
#    )
#process.patJetsPFlow.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
#                                                      "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0")
#process.patJetsPFlow.userData.userFunctionLabels = cms.vstring('secvtxMass')

# electrons
process.selectedPatElectrons.cut = cms.string('pt > 5.0 & abs(eta) < 2.4')
process.patElectrons.embedTrack = cms.bool(True)
process.selectedPatElectronsPFlow.cut = cms.string('pt > 5.0 & abs(eta) < 2.4')
process.patElectronsPFlow.embedTrack = cms.bool(True)
# muons
process.selectedPatMuons.cut = cms.string('pt > 5.0 & abs(eta) < 2.4 & isGlobalMuon() & globalTrack().hitPattern().numberOfValidMuonHits() > 0 & numberOfMatches() > 1 & innerTrack().numberOfValidHits()> 10 & innerTrack().hitPattern().numberOfValidPixelHits()>1 & innerTrack().trackerExpectedHitsOuter().numberOfHits() <3 & innerTrack().normalizedChi2() < 10 & globalTrack().normalizedChi2() < 10')
process.patMuons.embedTrack = cms.bool(True)
process.selectedPatMuonsPFlow.cut = cms.string("pt > 5.0 & abs(eta) < 2.4 & isGlobalMuon() & globalTrack().hitPattern().numberOfValidMuonHits() > 0 & numberOfMatches() > 1 & innerTrack().numberOfValidHits()> 10 & innerTrack().hitPattern().numberOfValidPixelHits()>1 & innerTrack().trackerExpectedHitsOuter().numberOfHits() <3 & innerTrack().normalizedChi2() < 10 & globalTrack().normalizedChi2() < 10")
process.patMuonsPFlow.embedTrack = cms.bool(True)
# taus
process.selectedPatTausPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 3")
process.selectedPatTaus.cut = cms.string("pt > 10.0 & abs(eta) < 3")
process.patTausPFlow.isoDeposits = cms.PSet()
process.patTaus.isoDeposits = cms.PSet()
# photons
process.patPhotonsPFlow.isoDeposits = cms.PSet()
process.patPhotons.isoDeposits = cms.PSet()


# Apply jet ID to all of the jets upstream. We aren't going to screw around
# with this, most likely. So, we don't really to waste time with it
# at the analysis level. 
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJetsPFlow")
                                        )

process.goodPatJetsAK5PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                   filterParams = pfJetIDSelector.clone(),
                                   src = cms.InputTag("selectedPatJetsAK5PF")
                                   )

# let it run

process.patseq = cms.Sequence(
    process.HLTfilter*
    process.scrapingVeto*
    process.HBHENoiseFilter*
    process.goodOfflinePrimaryVertices*
    process.primaryVertexFilter*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.patDefaultSequence*
    process.goodPatJetsPFlow*
    process.goodPatJetsAK5PF*
    process.prunedGenParticles
    )

##########################
## check HLT filter option
if not options.applyHLT:
    process.patseq.remove( process.HLTfilter )
    
#######################
## System8 Tree
#######################
if options.writeS8Tree:
    # load S8 tree maker module
    process.load("EDModule.Analyzer.S8TreeMaker_cfi")
    # configure tree maker
    process.S8TreeMaker.triggers = "TriggerResults::HLT"
    process.S8TreeMaker.primaryVertices = "goodOfflinePrimaryVertices"
    process.S8TreeMaker.jets = "selectedPatJetsAK5PF"
    process.S8TreeMaker.muons = "selectedPatMuons"
    process.S8TreeMaker.electrons = "selectedPatElectrons"
    # configure output
    process.TFileService = cms.Service(
            "TFileService",
            fileName = cms.string("s8_tree.root")
    )
    
    process.patseq.replace( process.goodPatJetsAK5PF,
                            process.goodPatJetsAK5PF *
                            process.S8TreeMaker )
        
#process.patseq.replace( process.goodOfflinePrimaryVertices,
#                        process.goodOfflinePrimaryVertices *
#                        process.eidCiCSequence )

if options.useData == True :
    process.patseq.remove( process.prunedGenParticles )
    removeMCMatching( process, ['All'] )
    
process.p0 = cms.Path(
    process.patseq
    )

if options.useData == True :
    process.p0.remove( process.patJetPartonAssociation )
    process.p0.remove( process.patJetPartonAssociationAK5PF )
    process.p0.remove( process.patJetPartonMatchAK5PF )
    process.p0.remove( process.patJetGenJetMatchAK5PF )
    process.p0.remove( process.patJetPartons )
    process.p0.remove( process.patJetFlavourAssociation )
    process.p0.remove( process.patJetFlavourAssociationAK5PF)
    
    
process.out.SelectEvents.SelectEvents = cms.vstring('p0')

# rename output file
if options.useData :
    process.out.fileName = cms.untracked.string('btagmu_data.root')
else :
    process.out.fileName = cms.untracked.string('btagmu_mc.root')


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)


# process all the events
process.maxEvents.input = 100
process.options.wantSummary = True
process.out.dropMetaData = cms.untracked.string("DROPPED")


process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")


process.out.outputCommands = [
    'drop *_cleanPat*_*_*',
    'keep *_selectedPat*_*_*',
    'keep *_goodPat*_*_*',
    'drop patJets_selectedPat*_*_*',
    'drop *_selectedPatJets_*_*',    
    'keep *_patMETs*_*_*',
    'keep *_goodOfflinePrimaryVertices*_*_*',    
    'drop patPFParticles_*_*_*',
    'drop patTaus_*_*_*',
    'keep *_cleanPatPhotonsTriggerMatch*_*_*',
    'keep *_cleanPatElectronsTriggerMatch*_*_*',
    'keep *_cleanPatMuonsTriggerMatch*_*_*',
    'keep *_cleanPatTausTriggerMatch*_*_*',
    'keep *_cleanPatJetsTriggerMatch*_*_*',
    'keep *_patMETsTriggerMatch*_*_*',
    'keep double_*PFlow*_*_PAT',
    'keep *_TriggerResults_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep *_prunedGenParticles_*_*',
    'keep recoBaseTagInfosOwned_selectedPatJets*_*_*'
    #'keep recoTracks_generalTracks_*_*'
    ]

if options.useData :
    process.out.outputCommands += ['drop *_MEtoEDMConverter_*_*',
                                   'keep LumiSummary_lumiProducer_*_*'
                                   ]
else :
    process.out.outputCommands += ['keep GenRunInfoProduct_generator_*_*',
                                   'keep GenEventInfoProduct_generator_*_*',
                                   'keep PileupSummaryInfos_*_*_*'
                                   ]

#open('junk.py','w').write(process.dumpPython())

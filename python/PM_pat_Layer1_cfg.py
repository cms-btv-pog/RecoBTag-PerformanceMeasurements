#
#  
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *


#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(-1)  ),
    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.MessageLogger.cerr.FwkReport.reportEvery=100

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0027/A0C39391-EB8C-DE11-8EF8-00144F0D84D8.root',
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0027/8A4E6C25-348B-DE11-A857-0030485C6782.root',
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/FAB061EC-078D-DE11-9473-001E0B470AC2.root',
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/F899160A-088D-DE11-A3AC-001CC4A6FB3A.root'
    ]
process.maxEvents.input = 1000

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
#process.GlobalTag.globaltag = 'MC31X_V5::All'

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#-- Tuning of Monte Carlo matching --------------------------------------------
# Also match with leptons of opposite charge
#process.electronMatch.checkCharge = False
process.muonMatch.checkCharge = False


#-- Extra Jet/MET collections -------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import *

addJetCollection(process,cms.InputTag('antikt5CaloJets'),
                 'AK5',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','Calo'),
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("antikt5GenJets")
                 )

#addJetCollection(process,cms.InputTag('sisCone5CaloJets'),
#                 'SC5',
#                 doJTA        = True,
#                 doBTagging   = True,
#                 jetCorrLabel = ('SC5','Calo'),
#                 doType1MET   = True,
#                 doL1Cleaning = True,
#                 doL1Counters = True,
#                 genJetCollection=cms.InputTag("sisCone5GenJets")
#                 )

#-- Tune contents of jet collections  ------- CHECK Do we need this?
#for jetName in ( '', 'SC5' ):
#    module = getattr(process,'allLayer1Jets'+jetName)
#    module.addTagInfos = False    # Remove tag infos
#    module.addJetID    = True     # Add JetID variables

#### PRE-SELECTION

#from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
process.countLayer1Muons.minNumber = cms.uint32(1)

#from BTag algorithm note
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
process.selectedLayer1Jets.cut = 'pt > 30. & abs(eta) < 2.4'

#from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
process.countLayer1Jets.minNumber = cms.uint32(2)

process.selectedLayer1Muons.cut      = cms.string('pt > 5. & abs(eta) < 2.4')
process.selectedLayer1Muons.cut      = cms.string('pt > 15. & abs(eta) < 2.4')

#from BTag algorithm note
#process.patAODJetTracksAssociator.cut  = "quality('highPurity') && pt > 1 && normalizedChi2 < 5"

#0.3 from BTag note. CHECK DO WE REALLY NEED THIS?
process.jetPartonMatch.maxDeltaR  = cms.double(0.3) 
process.jetPartonMatch.maxDPtRel  = cms.double(99999999999)

process.jetGenJetMatch.maxDeltaR  = cms.double(0.5)
process.jetGenJetMatch.maxDPtRel  = cms.double(99999999999)


#-- Execution path --------------------------------------- check do we need this?
# Re-name IC5 collections for uniformity
#process.cleanLayer1JetsIC5 = process.cleanLayer1Jets
#process.layer1METsIC5      = process.layer1METs

# Modify subsequent modules CHECK do we need this?
#process.cleanLayer1Hemispheres.patJets = process.cleanLayer1JetsIC5.label()
#process.countLayer1Jets.src            = process.cleanLayer1JetsIC5.label()

#-- Output module configuration -----------------------------------------------
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out.fileName = 'PM_pattuple.root'
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = [ 'drop *' ]

# Explicit list of collections to keep (basis is default PAT event content)
process.out.outputCommands.extend( [ # PAT Objects
                                     'keep *_cleanLayer1Muons_*_*',
                                     'keep *_cleanLayer1Jets*_*_*',       # All Jets
                                     # Generator information
                                     'keep GenEventInfoProduct_generator_*_*',
                                     # Generator particles/jets/MET
                                     'keep recoGenParticles_genParticles_*_*',
                                     'keep recoGenJets_iterativeCone5GenJets_*_*',
                                     'keep recoGenJets_antikt5GenJets_*_*',
                                     # Trigger information
                                     'keep edmTriggerResults_TriggerResults_*_HLT',
                                     'keep *_hltTriggerSummaryAOD_*_*',
                                     'keep L1GlobalTriggerObjectMapRecord_*_*_*',
                                     # Others
                                     'keep *_offlinePrimaryVertices_*_*',
                                     'keep *_offlineBeamSpot_*_*',
#                                     'keep *_towerMaker_*_*',                 #
                                     'keep recoTracks_generalTracks_*_*',
                                     'keep HcalNoiseSummary_*_*_*'
                                     ] )


#-- Trigger matching ----------------------------------------------------------
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerSequence.remove( process.patTriggerMatcher )
process.patTriggerEvent.patTriggerMatches  = ()


#-- Execution path ------------------------------------------------------------

#for jets in ( 'SC5' ):
#    process.allLayer1Summary.candidates.append(cms.InputTag('allLayer1Jets'+jets))
#    process.selectedLayer1Summary.candidates.append(cms.InputTag('selectedLayer1Jets'+jets))
#    process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1Jets'+jets))

# Full path
#process.p = cms.Path( process.patDefaultSequence*process.patTrigger*process.patTriggerEvent )
process.p = cms.Path( process.patDefaultSequence )


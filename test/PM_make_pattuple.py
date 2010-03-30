#
#  
#

import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

 
#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.FwkReport.reportEvery=100
process.MessageLogger.cerr = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(-1)  ),
    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    FwkReport = cms.untracked.PSet(
	reportEvery = cms.untracked.int32(100)
    )

)
#process.MessageLogger.cerr.FwkReport.reportEvery=100

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0027/A0C39391-EB8C-DE11-8EF8-00144F0D84D8.root',
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0027/8A4E6C25-348B-DE11-A857-0030485C6782.root',
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/FAB061EC-078D-DE11-9473-001E0B470AC2.root',
    '/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/F899160A-088D-DE11-A3AC-001CC4A6FB3A.root'
    ]
process.maxEvents.input = 10

#-- Calibration tag -----------------------------------------------------------
#process.GlobalTag.globaltag = 'MC31X_V5::All'

#-- load pat sequence ----
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_cff")
process.load("RecoBTag.PerformanceMeasurements.PM_pat_Layer1_cfg")
from PhysicsTools.PatAlgos.tools.jetTools import *

#####JEC for the 31X sample
process.jetCorrFactors.corrSample = 'Summer09_7TeV'

#####JEC To use initially while waiting the derived from MC truth from the CMSSW_3_5_5 sample
#process.jetCorrFactors.corrSample = 'Summer09_7TeV_ReReco332'

#switchJECSet(process, newName='Summer09_7TeV', oldName='Summer09')
#switchJECSet_(process.jetCorrFactorsAK5,newName='Summer09_7TeV', oldName='Summer09')

# run the 3.3.x software on Summer 09 MC from 3.1.x:
#   - change the name from "ak" (3.3.x) to "antikt) (3.1.x)
#   - run jet ID (not run in 3.1.x)
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
 
run33xOn31xMC( process,
               jetSrc = cms.InputTag("antikt5CaloJets"),
               jetIdTag = "antikt5"
               )

####### A O D ############
# If we are running in AOD, use the following switch
#
#restrictInputToAOD(process)
process.allLayer1Jets.addJetID = False
#process.allLayer1JetsAK5.addJetID = False

# Full path
#process.p = cms.Path( process.patDefaultSequence*process.patTrigger*process.patTriggerEvent )

# clean up jet filters in path
#process.PM_tuple.remove(process.cleanLayer1JetsAK5)
#process.PM_tuple.remove(process.countLayer1JetsAK5)
process.PM_tuple.remove(process.cleanLayer1Jets)
process.PM_tuple.remove(process.countLayer1Jets)
process.PM_tuple.remove(process.cleanLayer1Hemispheres)

# disable for the moment deltaR filters
process.p = cms.Path( process.PM_tuple*process.cleanLayer1Jets*process.countLayer1Jets) #*process.DeltaRCut )
#process.p2 = cms.Path( process.PM_tuple*process.cleanLayer1JetsAK5*process.countLayer1JetsAK5) #*process.DeltaRCutAK5 ) # for AK5 jets


#-- Output module configuration -----------------------------------------------
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out.fileName = 'PM_pattuple.root'
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = [ 'drop *' ]

#process.out.SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('p1','p2')
#)

# Explicit list of collections to keep (basis is default PAT event content)
process.out.outputCommands.extend( [ # PAT Objects
                                     'keep *_selectedLayer1Muons_*_*',
                                     'keep *_selectedLayer1Jets*_*_*',       # All Jets
                                     # Generator information
                                     'keep GenEventInfoProduct_generator_*_*',
                                     'keep GenRunInfoProduct_generator_*_*',
                                     # Generator particles/jets/MET
                                     'keep recoGenParticles_genParticles_*_*',
                                     #'keep recoGenJets_iterativeCone5GenJets_*_*',
                                     'keep recoGenJets_antikt5GenJets_*_*',
                                     # Trigger information
				     'keep edmTriggerResults_TriggerResults_*_HLT*',
                                     #'keep *_hltTriggerSummaryAOD_*_*',
                                     #'keep L1GlobalTriggerObjectMapRecord_*_*_*',
                                     # Others
                                     'keep *_offlinePrimaryVertices_*_*',
                                     'keep *_offlineBeamSpot_*_*',
#                                    'keep *_towerMaker_*_*',                 #
                                     'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*',
                                     'keep recoTracks_generalTracks_*_*',
				     #'keep *_*JetTrackAssociatorAtVertex*_*_*'
                                     'keep HcalNoiseSummary_*_*_*'
                                     ] )



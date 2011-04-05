#
# see also this example  
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatExamples/test/patLayer1_fromRECO_7TeV_firstdata_cfg.py?revision=1.1.2.2&view=markup&pathrev=V00-02-16
#

# 
#The good collisions runs and "good" lumi section ranges can be found from http://pccmsdqm04.cern.ch/runregistry/
#instruction here https://twiki.cern.ch/twiki/bin/viewauth/CMS/DQMRunRegistry
#

import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *
 
#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(-1)  ),
    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    FwkReport = cms.untracked.PSet(
	reportEvery = cms.untracked.int32(1000)
    )

)


#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
    'file:crab_0_101120_162509/res/test_1_1_vA7.root'
    ]

#-- EndImput Source -----------------------------------------------------------

process.maxEvents.input = 5000

#-- Select good events -----------------------------------------------------------
process.load("RecoBTag.PerformanceMeasurements.getEvent_cff")


#-- Calibration tag -----------------------------------------------------------
process.GlobalTag.globaltag = cms.string('START38_V12::All')
#-- Calibration tag -----------------------------------------------------------
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
             tag = cms.string("TrackProbabilityCalibration_2D_MCQpt30_startup_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("TrackProbabilityCalibration_3D_MCQpt30_startup_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
    )

#-------------------- TO RUN MC -----------------

#from PhysicsTools.PatAlgos.tools.coreTools import *

#removeMCMatching(process,['All'])

#restrictInputToAOD(process)
#not sure is still needed
#process.allPatJets.addJetID = False
#process.allPatJetsAK5.addJetID = False

##------------------------ JEC -----------------------
#from PhysicsTools.PatAlgos.tools.jetTools import *
##switchJECSet( process, "Summer09_7TeV_ReReco332")
#switchJECSet( process, "Spring10")

#----------- clean up jet filters in path -------------------------------------
#need to check if we need this clean up jet filters in path

#process.PM_tuple.remove(process.cleanPatJets)
#process.PM_tuple.remove(process.countPatJets)
#process.PM_tuple.remove(process.cleanPatHemispheres)

#----------- load pat sequence -------------------------------------
# load as last so every setting will be applyed to all jets collection added

from RecoBTag.PerformanceMeasurements.PM_pat_Layer1_cfg import *
process.patTrigger.processName = 'REDIGI38XPU'
process.patTriggerEvent.processName = 'REDIGI38XPU'

from PhysicsTools.PatAlgos.tools.jetTools import *

#process.HLTfilter = cms.EDFilter("HLTHighLevel",
##    TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
#    TriggerResultsTag  = cms.InputTag("TriggerResults","","REDIGI38XPU"),
#    HLTPaths           = cms.vstring("HLT_BTagMu*"),
##    HLTPaths           = cms.vstring("HLT_Jet*"),
#    eventSetupPathsKey = cms.string(''),
#    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
#    throw              = cms.bool(True)
#)

process.countPatJetsAK5Calo = process.countPatJets.clone()
process.countPatJetsAK5Calo.src = cms.InputTag("selectedPatJets")
process.countPatJetsAK5Calo.minNumber = 2
process.countPatJetsAK5Track = process.countPatJets.clone()
process.countPatJetsAK5Track.src = cms.InputTag("selectedPatJetsAK5Track")
process.countPatJetsAK5Track.minNumber = 2
process.countPatJetsAK5PF = process.countPatJets.clone()
process.countPatJetsAK5PF.src = cms.InputTag("selectedPatJetsAK5PF")
process.countPatJetsAK5PF.minNumber = 2

#-------------------- TO RUN ON MC -----------------
process.pCalo = cms.Path(
#    process.HLTfilter*
    process.getEventMC*
    process.btagging*
    process.PM_tuple*
    process.countPatJetsAK5Calo)
process.pTrack = cms.Path(
#    process.HLTfilter*
    process.getEventMC*
    process.btagging*
    process.PM_tuple*
    process.countPatJetsAK5Track)
process.pPF = cms.Path(
#    process.HLTfilter*
    process.getEventMC*
    process.btagging*
    process.PM_tuple*
    process.countPatJetsAK5PF)

#-------- OUTPUT MODULE configuration -----------------------------------------------

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out.fileName = 'PM_pattuple_MC.root'
#process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files)
#process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = [ 'drop *' ]

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('pCalo','pTrack','pPF')
)

# Explicit list of collections to keep (basis is default PAT event content)
process.out.outputCommands.extend( [ # PAT Objects
                                     'keep *_selectedPatMuons*_*_*',
#                                     'keep *_selectedPatElectrons*_*_*',
                                     'keep *_selectedPatJets*_*_*',       # All Jets
                                     # Generator information
                                     'keep GenEventInfoProduct_generator_*_*',
                                     'keep GenRunInfoProduct_generator_*_*',
                                     # Generator particles/jets/MET
                                     'keep recoGenParticles_genParticles_*_*',
                                     'keep recoGenJets_ak5GenJets_*_*',
                                     # MET information
                                     'keep *_patMETs*_*_*',            # All METs
                                     'keep recoGenMETs_*_*_*',
                                     # Luminosity information
                                     'keep edmMergeableCounter_eventCountProducer_*_*',
                                     # Trigger information
				     'keep edmTriggerResults_TriggerResults_*_*',
				     'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
                                     #'keep *_hltTriggerSummaryAOD_*_*',
                                     #'keep L1GlobalTriggerObjectMapRecord_*_*_*',
                                     # Others
                                     'keep *_offlinePrimaryVertices_*_*',
                                     'keep *_offlineBeamSpot_*_*',
#                                    'keep *_towerMaker_*_*',
                                     'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', # waiting the trackJets development for PAT
                                     'keep recoTracks_generalTracks_*_*',
				     #'keep *_*JetTrackAssociatorAtVertex*_*_*'
                                     'keep HcalNoiseSummary_*_*_*'
                                     ] )


process.schedule = cms.Schedule(process.pCalo,process.pTrack,process.pPF
                                ,process.outpath
                                )

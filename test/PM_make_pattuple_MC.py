#
# see also this example  
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatExamples/test/patLayer1_fromRECO_7TeV_firstdata_cfg.py?revision=1.1.2.2&view=markup&pathrev=V00-02-16
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

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [

#'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0007/FCDD5B7A-233C-DF11-AB79-002618943926.root',
#'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0007/FC5FDFDD-233C-DF11-9FC5-002618943964.root',
#'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0007/FA77A575-233C-DF11-9D9A-00248C55CC62.root',
#'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0007/FA7205B6-233C-DF11-9BC2-002618943932.root',
#'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0007/F850EA9F-233C-DF11-8D05-00261894383C.root'
     
'file:/pcussitrk002/data1/dalfonso/MC132440/FCDD5B7A-233C-DF11-AB79-002618943926.root',
'file:/pcussitrk002/data1/dalfonso/MC132440/FC5FDFDD-233C-DF11-9FC5-002618943964.root',
'file:/pcussitrk002/data1/dalfonso/MC132440/FA77A575-233C-DF11-9D9A-00248C55CC62.root'

    ]
process.maxEvents.input = -1

#-- Select good events -----------------------------------------------------------
process.load("RecoBTag.PerformanceMeasurements.getEvent_cff")

#-- Calibration tag -----------------------------------------------------------
process.GlobalTag.globaltag = cms.string('START3X_V25B::All')

#-------------------- TO RUN DATA -----------------

#from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process,['All'])

#restrictInputToAOD(process)
#not sure is still needed
#process.allPatJets.addJetID = False
#process.allPatJetsAK5.addJetID = False

#------------------------ JEC -----------------------

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Summer09_7TeV_ReReco332")


#----------- clean up jet filters in path -------------------------------------
#need to check if we need this clean up jet filters in path

#process.PM_tuple.remove(process.cleanPatJets)
#process.PM_tuple.remove(process.countPatJets)
#process.PM_tuple.remove(process.cleanPatHemispheres)


#----------- load pat sequence -------------------------------------
# load as last so every setting will be applyed to all jets collection added

process.load("RecoBTag.PerformanceMeasurements.PM_pat_Layer1_cfg")
from PhysicsTools.PatAlgos.tools.jetTools import *

#-------------------- TO RUN ON DATA -----------------

process.p = cms.Path( process.getEventMC*process.PM_tuple)


#-------- OUTPUT MODULE configuration -----------------------------------------------

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out.fileName = 'PM_pattuple_MC.root'
#process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files)
#process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = [ 'drop *' ]

#process.out.SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('p1','p2')
#)

# Explicit list of collections to keep (basis is default PAT event content)
process.out.outputCommands.extend( [ # PAT Objects
                                     'keep *_selectedPatMuons_*_*',
                                     'keep *_selectedPatJets*_*_*',       # All Jets
                                     # Generator information
                                     'keep GenEventInfoProduct_generator_*_*',
                                     'keep GenRunInfoProduct_generator_*_*',
                                     # Generator particles/jets/MET
                                     'keep recoGenParticles_genParticles_*_*',
                                     'keep recoGenJets_ak5GenJets_*_*',
                                     # Trigger information
				     'keep edmTriggerResults_TriggerResults_*_HLT*',
                                     #'keep *_hltTriggerSummaryAOD_*_*',
                                     #'keep L1GlobalTriggerObjectMapRecord_*_*_*',
                                     # Others
                                     'keep *_offlinePrimaryVertices_*_*',
                                     'keep *_offlineBeamSpot_*_*',
#                                    'keep *_towerMaker_*_*',                 #
                                     'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', # waiting the trackJets development for PAT
                                     'keep recoTracks_generalTracks_*_*',
				     #'keep *_*JetTrackAssociatorAtVertex*_*_*'
                                     'keep HcalNoiseSummary_*_*_*'
                                     ] )



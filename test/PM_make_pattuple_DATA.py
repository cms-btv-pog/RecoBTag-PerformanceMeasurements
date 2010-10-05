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

 '/store/data/Run2010A/Mu/RECO/Jun14thReReco_v1/0000/FE7AAFE8-A878-DF11-A8CF-0017A4770824.root',
 '/store/data/Run2010A/Mu/RECO/Jun14thReReco_v1/0000/FE86F550-8479-DF11-B056-00237DF345D6.root',
 '/store/data/Run2010A/Mu/RECO/Jun14thReReco_v1/0000/FEFAF52B-AE78-DF11-878C-0017A477103C.root'

    ]

#-- EndImput Source -----------------------------------------------------------

process.maxEvents.input = -1

#-- Select good events -----------------------------------------------------------
process.load("RecoBTag.PerformanceMeasurements.getEvent_cff")


#-- Calibration tag -----------------------------------------------------------
process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')

#-------------------- TO RUN DATA -----------------

from PhysicsTools.PatAlgos.tools.coreTools import *

removeMCMatching(process,['All'])

#restrictInputToAOD(process)
#not sure is still needed
#process.allPatJets.addJetID = False
#process.allPatJetsAK5.addJetID = False

#------------------------ JEC -----------------------
from PhysicsTools.PatAlgos.tools.jetTools import *
##switchJECSet( process, "Summer09_7TeV_ReReco332")
switchJECSet( process, "Spring10")

#----------- clean up jet filters in path -------------------------------------
#need to check if we need this clean up jet filters in path

#process.PM_tuple.remove(process.cleanPatJets)
#process.PM_tuple.remove(process.countPatJets)
#process.PM_tuple.remove(process.cleanPatHemispheres)

#----------- load pat sequence -------------------------------------
# load as last so every setting will be applyed to all jets collection added

from RecoBTag.PerformanceMeasurements.PM_pat_Layer1_cfg import *

from PhysicsTools.PatAlgos.tools.jetTools import *

#-------------------- TO RUN ON DATA -----------------
process.p = cms.Path(
    process.getEventDATA*
    process.PM_tuple)

#-------- OUTPUT MODULE configuration -----------------------------------------------

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out.fileName = 'PM_pattuple_DATA.root'
#process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files)
#process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = [ 'drop *' ]

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
)

# Explicit list of collections to keep (basis is default PAT event content)
process.out.outputCommands.extend( [ # PAT Objects
                                     'keep *_selectedPatMuons_*_*',
	                             'keep *_selectedPatMuonsForPtRel_*_*',
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



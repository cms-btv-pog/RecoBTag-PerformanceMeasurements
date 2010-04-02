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

#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F4C92A98-163C-DF11-9788-0030487C7392.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F427D642-173C-DF11-A909-0030487C60AE.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/E27821C3-0C3C-DF11-9BD9-0030487CD718.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/D87D5469-2E3C-DF11-A470-000423D99896.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/B647CAD9-0E3C-DF11-886F-0030487CD716.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/A860D55E-193C-DF11-BE29-0030487C60AE.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/9884BC11-0C3C-DF11-8F9C-000423D986C4.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/92684831-233C-DF11-ABA0-0030487CD16E.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/90269E76-0D3C-DF11-A1A0-0030487CD840.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/8CAE3014-133C-DF11-A05D-000423D174FE.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/8C51BAC6-1A3C-DF11-A0EE-000423D94A04.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/8C042B04-2D3C-DF11-939F-0030487CD178.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/80471A6B-0E3C-DF11-8DCD-0030487C6A66.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/762824C3-0C3C-DF11-A4FD-0030487CD6D2.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/6A3533F5-103C-DF11-B3AA-00304879BAB2.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/4C8979D2-073C-DF11-B97B-000423D6AF24.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/26C8DED9-0E3C-DF11-9D83-0030487CD7B4.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/181C44F7-093C-DF11-A9CB-001D09F24FEC.root',
#'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/0AA7C390-0F3C-DF11-BD65-000423D998BA.root'

'file:/pcussitrk002/data1/dalfonso/132440/F4C92A98-163C-DF11-9788-0030487C7392.root',
'file:/pcussitrk002/data1/dalfonso/132440/F427D642-173C-DF11-A909-0030487C60AE.root',
'file:/pcussitrk002/data1/dalfonso/132440/E27821C3-0C3C-DF11-9BD9-0030487CD718.root',
'file:/pcussitrk002/data1/dalfonso/132440/D87D5469-2E3C-DF11-A470-000423D99896.root',
'file:/pcussitrk002/data1/dalfonso/132440/B647CAD9-0E3C-DF11-886F-0030487CD716.root',
'file:/pcussitrk002/data1/dalfonso/132440/A860D55E-193C-DF11-BE29-0030487C60AE.root',
'file:/pcussitrk002/data1/dalfonso/132440/9884BC11-0C3C-DF11-8F9C-000423D986C4.root',
'file:/pcussitrk002/data1/dalfonso/132440/92684831-233C-DF11-ABA0-0030487CD16E.root',
'file:/pcussitrk002/data1/dalfonso/132440/90269E76-0D3C-DF11-A1A0-0030487CD840.root',
'file:/pcussitrk002/data1/dalfonso/132440/8CAE3014-133C-DF11-A05D-000423D174FE.root',
'file:/pcussitrk002/data1/dalfonso/132440/8C51BAC6-1A3C-DF11-A0EE-000423D94A04.root',
'file:/pcussitrk002/data1/dalfonso/132440/8C042B04-2D3C-DF11-939F-0030487CD178.root',
'file:/pcussitrk002/data1/dalfonso/132440/80471A6B-0E3C-DF11-8DCD-0030487C6A66.root',
'file:/pcussitrk002/data1/dalfonso/132440/762824C3-0C3C-DF11-A4FD-0030487CD6D2.root',
'file:/pcussitrk002/data1/dalfonso/132440/6A3533F5-103C-DF11-B3AA-00304879BAB2.root',
'file:/pcussitrk002/data1/dalfonso/132440/4C8979D2-073C-DF11-B97B-000423D6AF24.root',
'file:/pcussitrk002/data1/dalfonso/132440/26C8DED9-0E3C-DF11-9D83-0030487CD7B4.root',
'file:/pcussitrk002/data1/dalfonso/132440/181C44F7-093C-DF11-A9CB-001D09F24FEC.root',
'file:/pcussitrk002/data1/dalfonso/132440/0AA7C390-0F3C-DF11-BD65-000423D998BA.root'

    ]
process.maxEvents.input = -1

#-- Select good events -----------------------------------------------------------
process.load("RecoBTag.PerformanceMeasurements.getEvent_cff")

#-- Calibration tag -----------------------------------------------------------
process.GlobalTag.globaltag = cms.string('GR10_P_V2::All')

#-------------------- TO RUN DATA -----------------

from PhysicsTools.PatAlgos.tools.coreTools import *

removeMCMatching(process,['All'])

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

process.p = cms.Path( process.getEventDATA*process.PM_tuple)


#-------- OUTPUT MODULE configuration -----------------------------------------------

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out.fileName = 'PM_pattuple_Run132440.root'
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
                                     # Luminosity information
                                     'keep edmMergeableCounter_eventCountProducer_*_*',
                                     # Trigger information
				     'keep edmTriggerResults_TriggerResults_*_HLT*',
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



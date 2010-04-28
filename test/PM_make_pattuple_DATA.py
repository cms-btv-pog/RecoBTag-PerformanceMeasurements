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

########  run=132440 #######
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0139/3ADE63D6-923E-DF11-B92A-001A92971BD8.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0133/6E74370E-6D3E-DF11-890B-0018F3D09684.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/C851858F-613E-DF11-A403-00261894387E.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/C0938B90-613E-DF11-8F45-00261894386D.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/9A9D3396-613E-DF11-8643-00261894398B.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/02CA6895-613E-DF11-A4F0-003048678B92.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/E024013C-5D3E-DF11-8F34-0018F3D09682.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/DE499234-603E-DF11-94A6-003048678B06.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/BAFFF7F1-5A3E-DF11-9F13-00248C0BE018.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/A8F6AE35-5D3E-DF11-8601-002354EF3BE0.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/8C884021-5B3E-DF11-9D5E-003048678ADA.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/7AA74636-603E-DF11-A5DE-003048678F8A.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/62DE6421-5B3E-DF11-B53B-003048678FA6.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/58C598F2-5A3E-DF11-9B89-00304867BEE4.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/40EA1340-5D3E-DF11-B49E-0026189438D7.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/FAED0412-5A3E-DF11-B5B7-002618943870.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/EC8E3A0D-5A3E-DF11-A4CC-002618943900.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/D2DB4C0A-593E-DF11-97AE-00304867902C.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/C82F420B-5A3E-DF11-B9FD-002618943948.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/BED9F431-593E-DF11-925F-003048D15E24.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/98FEF80A-5A3E-DF11-B1AF-00248C0BE018.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/7CF5670C-5A3E-DF11-8C52-003048678F1C.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/5AAA5A6D-553E-DF11-AEE5-002354EF3BE6.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/50BFA301-593E-DF11-9B2E-00304867918E.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/50178414-593E-DF11-851A-0030486790B0.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/14C3AF58-563E-DF11-B19B-00261894396F.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/1092AD71-973E-DF11-8F4E-0018F3D096D8.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0129/58274A60-523E-DF11-A311-0030486792F0.root',
'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0129/322D2491-533E-DF11-9711-003048679000.root'

    ]

#-- EndImput Source -----------------------------------------------------------

process.maxEvents.input = -1

#-- Select good events -----------------------------------------------------------
process.load("RecoBTag.PerformanceMeasurements.getEvent_cff")

#-- Calibration tag -----------------------------------------------------------
#process.GlobalTag.globaltag = cms.string('GR10_P_V4::All')  # for the /MinimumBias/Commissioning10-GOODCOLL-v8/RAW-RECO (run 132514-today)
#process.GlobalTag.globaltag = cms.string('GR_R_35X_V6::All') # for the /MinimumBias/Commissioning10-Apr1Skim_GOODCOLL-v1/RAW-RECO (run 132440-132513)
process.GlobalTag.globaltag = cms.string('GR_R_35X_V7A::All') # for the /MinimumBias/Commissioning10-Apr20Skim_GOODCOLL-v1/RAW-RECO

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

from RecoBTag.PerformanceMeasurements.PM_pat_Layer1_cfg import *

from PhysicsTools.PatAlgos.tools.jetTools import *

#-------------------- TO RUN ON DATA -----------------
process.p = cms.Path( process.getEventDATA*process.PM_tuple)

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



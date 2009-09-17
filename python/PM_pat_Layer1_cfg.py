#
#  
#

# load PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *

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

#-- Tune contents of jet collections  ------- CHECK Do we need this?
#for jetName in ( '', 'SC5' ):
#    module = getattr(process,'allLayer1Jets'+jetName)
#    module.addTagInfos = False    # Remove tag infos
#    module.addJetID    = True     # Add JetID variables

#### PRE-SELECTION

#from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *

process.selectedLayer1Muons.cut      = cms.string('pt > 5. & abs(eta) < 2.4')
process.selectedLayer1Jets.cut = 'pt > 30. & abs(eta) < 2.4'
process.countLayer1Jets.minNumber = cms.uint32(2)
process.countLayer1Muons.minNumber = cms.uint32(1)

#from BTag algorithm note
#process.patAODJetTracksAssociator.cut  = "quality('highPurity') && pt > 1 && normalizedChi2 < 5"

#0.3 from BTag note. CHECK DO WE REALLY NEED THIS?
process.jetPartonMatch.maxDeltaR  = cms.double(0.3) 
process.jetPartonMatch.maxDPtRel  = cms.double(99999999999)

process.jetGenJetMatch.maxDeltaR  = cms.double(0.5)
process.jetGenJetMatch.maxDPtRel  = cms.double(99999999999)


#-- Trigger matching ----------------------------------------------------------
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerSequence.remove( process.patTriggerMatcher )
process.patTriggerEvent.patTriggerMatches  = ()


#for jets in ( 'SC5' ):
#    process.allLayer1Summary.candidates.append(cms.InputTag('allLayer1Jets'+jets))
#    process.selectedLayer1Summary.candidates.append(cms.InputTag('selectedLayer1Jets'+jets))
#    process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1Jets'+jets))

# Sequence
#process.p = cms.Path( process.patDefaultSequence*process.patTrigger*process.patTriggerEvent )

PM_tuple = cms.Sequence( process.patDefaultSequence )


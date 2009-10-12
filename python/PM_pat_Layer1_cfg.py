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

#-- Tune contents of jet collections  ------- CHECK Do we need this?
#for jetName in ( '', 'AK5' ):
#    module = getattr(process,'allLayer1Jets'+jetName)
#    module.addJetID    = True     # Add JetID variables
#    module.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfos"), cms.InputTag("softMuonTagInfos"), cms.InputTag("impactParameterTagInfos")) # remove softelectronTagInfos
#    module.embedGenJetMatch = False # Only keep reference, since we anyway keep the genJet collections
#    module.discriminatorSources = cms.VInputTag(
#	cms.InputTag("jetProbabilityBJetTags"), 
#	cms.InputTag("simpleSecondaryVertexBJetTags"), 
#        cms.InputTag("softMuonBJetTags"), 
#	cms.InputTag("softMuonByPtBJetTags"), 
#	cms.InputTag("softMuonByIP3dBJetTags"), 
#        cms.InputTag("trackCountingHighEffBJetTags"), 
#	cms.InputTag("trackCountingHighPurBJetTags")
#	)
# something is wrong with the previous loop, that screws up the module name and sequence...
# let's change them by hand
#process.allLayer1Jets.addJetID = True
process.allLayer1Jets.embedGenJetMatch = False
#process.allLayer1Jets.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfos"), cms.InputTag("softMuonTagInfos"), cms.InputTag("impactParameterTagInfos")) # remove softelectronTagInfos
#process.allLayer1Jets.discriminatorSources = cms.VInputTag(
#    cms.InputTag("jetProbabilityBJetTags"), 
#    cms.InputTag("simpleSecondaryVertexBJetTags"), 
#    cms.InputTag("softMuonBJetTags"), 
#    cms.InputTag("softMuonByPtBJetTags"), 
#    cms.InputTag("softMuonByIP3dBJetTags"), 
#    cms.InputTag("trackCountingHighEffBJetTags"), 
#    cms.InputTag("trackCountingHighPurBJetTags")
#    )
#process.allLayer1JetsAK5.addJetID = True
## process.allLayer1JetsAK5.embedGenJetMatch = False
#process.allLayer1JetsAK5.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfos"), cms.InputTag("softMuonTagInfos"), cms.InputTag("impactParameterTagInfos")) # remove softelectronTagInfos
#process.allLayer1JetsAK5.discriminatorSources = cms.VInputTag(
#    cms.InputTag("jetProbabilityBJetTags"), 
#    cms.InputTag("simpleSecondaryVertexBJetTags"), 
#    cms.InputTag("softMuonBJetTags"), 
#    cms.InputTag("softMuonByPtBJetTags"), 
#    cms.InputTag("softMuonByIP3dBJetTags"), 
#    cms.InputTag("trackCountingHighEffBJetTags"), 
#    cms.InputTag("trackCountingHighPurBJetTags")
#    )


#### PRE-SELECTION

#from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *

process.selectedLayer1Muons.cut      = cms.string('pt > 5. & abs(eta) < 2.4 & isGlobalMuon() & innerTrack().numberOfValidHits()> 7')
process.selectedLayer1Jets.cut = 'pt > 30. & abs(eta) < 2.4'
## #&#process.selectedLayer1JetsAK5.cut = 'pt > 30. & abs(eta) < 2.4'
process.countLayer1Jets.minNumber = cms.uint32(2)
## process.countLayer1JetsAK5.minNumber = cms.uint32(2)
process.countLayer1Muons.minNumber = cms.uint32(1)

#from BTag algorithm note
#process.patAODJetTracksAssociator.cut  = "quality('highPurity') && pt > 1 && normalizedChi2 < 5"

#0.3 from BTag note. CHECK DO WE REALLY NEED THIS?
process.jetPartonMatch.maxDeltaR  = cms.double(0.3) 
process.jetPartonMatch.maxDPtRel  = cms.double(99999999999)
process.jetGenJetMatch.maxDeltaR  = cms.double(0.5)
process.jetGenJetMatch.maxDPtRel  = cms.double(99999999999)
## process.jetPartonMatchAK5.maxDeltaR  = cms.double(0.3) 
## process.jetPartonMatchAK5.maxDPtRel  = cms.double(99999999999)
## process.jetGenJetMatchAK5.maxDeltaR  = cms.double(0.5)
## process.jetGenJetMatchAK5.maxDPtRel  = cms.double(99999999999)


# deltaR cut filters
DeltaRCut = cms.EDFilter("PMDeltaRFilter",
     Jets = cms.InputTag("selectedLayer1Jets"),
     Muons = cms.InputTag("selectedLayer1Muons"),
     MaxDeltaR = cms.double(0.4)
)

DeltaRCutAK5 = cms.EDFilter("PMDeltaRFilter",
     Jets = cms.InputTag("selectedLayer1JetsAK5"),
     Muons = cms.InputTag("selectedLayer1Muons"),
     MaxDeltaR = cms.double(0.4)
)
 

#-- Trigger matching ----------------------------------------------------------
#from PhysicsTools.PatAlgos.tools.trigTools import *
#switchOnTrigger( process )
#process.patTriggerSequence.remove( process.patTriggerMatcher )
#process.patTriggerEvent.patTriggerMatches  = ()


#for jets in ( 'SC5' ):
#    process.allLayer1Summary.candidates.append(cms.InputTag('allLayer1Jets'+jets))
#    process.selectedLayer1Summary.candidates.append(cms.InputTag('selectedLayer1Jets'+jets))
#    process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1Jets'+jets))

#
# add new jet collections at the end (in order to copy all modifications
#   to the default jet collection)
#
addJetCollection(process,
		 cms.InputTag('antikt5CaloJets'),
                 'AK5',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','Calo'),
                 doType1MET   = True,
		 doL1Cleaning = True,
		 doL1Counters = True, # NOTE, you need to create two paths and treat carefully the counters
                 genJetCollection=cms.InputTag("antikt5GenJets")
                 )



# Sequence
#process.p = cms.Path( process.patDefaultSequence*process.patTrigger*process.patTriggerEvent )


PM_tuple = cms.Sequence( process.patDefaultSequence )


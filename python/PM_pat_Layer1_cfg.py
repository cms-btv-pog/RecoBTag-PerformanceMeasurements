#
#  
#

# load PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#-- Extra Jet/MET collections -------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.metTools import *

#-----------------------------------------------------  Tune content ----------------------------------------------------------

process.patJets.embedGenJetMatch = False # Only keep reference, since we anyway keep the genJet collections

#from BTag algorithm note
#process.patAODJetTracksAssociator.cut  = "quality('highPurity') && pt > 1 && normalizedChi2 < 5"

#maxDR=0.3 from BTag note. 
process.patJetPartonMatch.maxDeltaR  = cms.double(0.3) 
process.patJetPartonMatch.maxDPtRel  = cms.double(99999999999)
process.patJetGenJetMatch.maxDeltaR  = cms.double(0.5)
process.patJetGenJetMatch.maxDPtRel  = cms.double(99999999999)

# Also match with leptons of opposite charge
#process.electronMatch.checkCharge = False
#process.muonMatch.checkCharge = False

# remove the complex tagger at the beginnig to save space: for example remove softelectronTagInfos and combinedSV

process.patJets.tagInfoSources = cms.VInputTag(
cms.InputTag("secondaryVertexTagInfos"), 
cms.InputTag("softMuonTagInfos"), 
cms.InputTag("impactParameterTagInfos")) 

process.patJets.discriminatorSources = cms.VInputTag(
    cms.InputTag("jetProbabilityBJetTags"), 
    cms.InputTag("simpleSecondaryVertexBJetTags"),
    cms.InputTag("softMuonBJetTags"), 
    cms.InputTag("softMuonByPtBJetTags"), 
    cms.InputTag("softMuonByIP3dBJetTags"), 
    cms.InputTag("trackCountingHighEffBJetTags"), 
    cms.InputTag("trackCountingHighPurBJetTags")
    )

#----------------- Trigger matching ----------------------------------------------------------
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

#----------------- Add other jet collection  ----------------------------------------------------------

addJetCollection(process, cms.InputTag('ak5TrackJets'),
                 'AK5', 'Track',
                 doJTA = False,
                 doBTagging   = True,
                 jetCorrLabel = None,
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )

addJetCollection(process,
                 cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = False,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )

#-----------------------------------------------------  Pre-selection ----------------------------------------------------------

# ususally pt(caloJet)> 30 and pt(muon)>5. Now lowered at the beginning
process.selectedPatMuons.cut= cms.string('pt > 3. & abs(eta) < 2.4 & isGlobalMuon() & innerTrack().numberOfValidHits()> 7')

process.selectedPatJets.cut = cms.string('pt > 10. & abs(eta) < 2.4')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 8. & abs(eta) < 2.4')
process.selectedPatJetsAK5Track.cut = cms.string('pt > 8. & abs(eta) < 2.4')

#process.countPatMuons.minNumber = cms.uint32(1)
#process.countPatJets.minNumber = cms.uint32(2) # commented to avoid bias against other jet collections

# deltaR cut filters
#DeltaRCut = cms.EDFilter("PMDeltaRFilter",
#     Jets = cms.InputTag("selectedPatJets"),
#     Muons = cms.InputTag("selectedPatMuons"),
#     MaxDeltaR = cms.double(0.4)
#)

#DeltaRCutAK5 = cms.EDFilter("PMDeltaRFilter",
#     Jets = cms.InputTag("selectedPatJetsAK5"),
#     Muons = cms.InputTag("selectedPatMuons"),
#     MaxDeltaR = cms.double(0.4)
#)
 

# Sequence
#process.p = cms.Path( process.patDefaultSequence*process.patTrigger*process.patTriggerEvent )

PM_tuple = cms.Sequence( process.patDefaultSequence )


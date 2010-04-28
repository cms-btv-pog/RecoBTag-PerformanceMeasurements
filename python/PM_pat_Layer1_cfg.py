#
#
# load PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#-- Extra Jet/MET collections -------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

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

#----------------- Trigger matching ----------------------------------------------------------
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

#----------------- Add other jet collection  ----------------------------------------------------------
addJetCollection(process, 
                 cms.InputTag('ak5TrackJets'),
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


#-----------------------------------------------------  slim objects ----------------------------------------------------------

#-----------------------------------------------------  B-tagging ----------------------------------------------------------

process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex3TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")

# remove the complex tagger at the beginnig to save space: for example remove softelectronTagInfos and combinedSV, jetBProb

theJetNames = ['','AK5PF','AK5Track']

for jetName in theJetNames:

###### needed temporarily for the Secondary vertex backporting

   module = setattr( process, 'simpleSecondaryVertexHighPurBJetTags'+jetName, process.simpleSecondaryVertexHighPurBJetTags.clone(tagInfos = cms.VInputTag("secondaryVertexTagInfos"+jetName)) )

###### to slim the tagInfo
   module = getattr(process,'patJets'+jetName)

   module.tagInfoSources = cms.VInputTag(
    cms.InputTag("impactParameterTagInfos"+jetName),
    cms.InputTag("secondaryVertexTagInfos"+jetName),
    cms.InputTag("softMuonTagInfos"+jetName), 
   ) 

   module.discriminatorSources = cms.VInputTag(
#    cms.InputTag("jetProbabilityBJetTags"+jetName), 
    cms.InputTag("simpleSecondaryVertexBJetTags"+jetName),
    cms.InputTag("simpleSecondaryVertexHighPurBJetTags"+jetName),
    cms.InputTag("softMuonBJetTags"+jetName), 
    cms.InputTag("softMuonByPtBJetTags"+jetName), 
    cms.InputTag("softMuonByIP3dBJetTags"+jetName), 
    cms.InputTag("trackCountingHighEffBJetTags"+jetName), 
    cms.InputTag("trackCountingHighPurBJetTags"+jetName)
   )

#-----------------------------------------------------  Pre-selection ----------------------------------------------------------

# ususally pt(caloJet)> 30 and pt(muon)>5. Now lowered at the beginning
process.selectedPatMuons.cut= cms.string('pt > 3. & abs(eta) < 2.4 & isGlobalMuon() & innerTrack().numberOfValidHits()> 7')

process.selectedPatJets.cut = cms.string('pt > 15. & abs(eta) < 2.4 & emEnergyFraction() >0.01 & jetID().n90Hits>1 & jetID().fHPD < 0.98  ')
#process.selectedPatJetsAK5PF.cut = cms.string('pt > 15. & abs(eta) < 2.4 & neutralHadronEnergyFraction() < 1.0 & neutralEmEnergyFraction() < 1.0 & nConstituents() > 1 & chargedHadronEnergyFraction() > 0.0 & chargedMultiplicity() > 0.0 & chargedEmEnergyFraction() < 1.0')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 15. & abs(eta) < 2.4')
process.selectedPatJetsAK5Track.cut = cms.string('pt > 10. & abs(eta) < 2.4')

process.countPatMuons.minNumber = cms.uint32(1)
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

process.PM_tuple = cms.Sequence( process.simpleSecondaryVertexHighPurBJetTags*process.patDefaultSequence )

process.PM_tuple.replace(process.simpleSecondaryVertexBJetTagsAK5PF, process.simpleSecondaryVertexBJetTagsAK5PF*process.simpleSecondaryVertexHighPurBJetTagsAK5PF)
process.PM_tuple.replace(process.simpleSecondaryVertexBJetTagsAK5Track, process.simpleSecondaryVertexBJetTagsAK5Track*process.simpleSecondaryVertexHighPurBJetTagsAK5Track)


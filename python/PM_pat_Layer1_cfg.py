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

<<<<<<< PM_pat_Layer1_cfg.py
=======
#
# add new jet collections at the end (in order to copy all modifications
#   to the default jet collection)
#
#addJetCollection(process,
#		 cms.InputTag('antikt5CaloJets'),
#                 'AK5',
#                 doJTA        = True,
#                 doBTagging   = True,
#                 jetCorrLabel = ('AK5','Calo'),
#                 doType1MET   = True,
#		 doL1Cleaning = True,
#		 doL1Counters = True, # NOTE, you need to create two paths and treat carefully the counters
#                genJetCollection=cms.InputTag("antikt5GenJets")
#                 )

#
# Track-Jet Sequence
#
#from RecoJets.JetProducers.ak5TrackJets_cfi import ak5TrackJets
#from CommonTools.RecoAlgos.TrackWithVertexRefSelector_cfi import *
#from RecoJets.JetProducers.TracksForJets_cff import *

#ak7TrackJets = ak5TrackJets.clone( rParam = 0.7 )

process.load('RecoJets.Configuration.RecoTrackJets_cff')

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
ak7GenJets   = ak5GenJets.clone( rParam = 0.7 )

recoTrackJets   = cms.Sequence(genParticlesForJets * ak5GenJets )

#recoTrackJets   = cms.Sequence(genParticlesForJets * ak5GenJets +
#                               trackWithVertexRefSelector+
#                               trackRefsForJets+
#                               ak5TrackJets #ak7TrackJets
#                               )
>>>>>>> 1.9

#----------------- Add other jet collection  ----------------------------------------------------------

addJetCollection(process,
                 cms.InputTag('ak5TrackJets'),
<<<<<<< PM_pat_Layer1_cfg.py
                 'AK5', 'Track',
                 doJTA        = False,
=======
                 'AKT5',
                 doJTA        = False,
>>>>>>> 1.9
                 doBTagging   = True,
                 jetCorrLabel = None,
<<<<<<< PM_pat_Layer1_cfg.py
                 doType1MET   = False,
                 doL1Cleaning = False,                 
=======
                 doType1MET   = True,
                 doL1Cleaning = False,
>>>>>>> 1.9
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
<<<<<<< PM_pat_Layer1_cfg.py
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
=======
                 doJetID      = False,
                 jetIdLabel   = "ak5"
>>>>>>> 1.9
                 )

#-----------------------------------------------------  Pre-selection ----------------------------------------------------------

# ususally pt(caloJet)> 30 and pt(muon)>5. Now lowered at the beginning
process.selectedPatMuons.cut      = cms.string('pt > 3. & abs(eta) < 2.4 & isGlobalMuon() & innerTrack().numberOfValidHits()> 7')

process.selectedPatJets.cut = 'pt > 10. & abs(eta) < 2.4'
process.selectedPatJetsAK5PF.cut = 'pt > 8. & abs(eta) < 2.4'
process.selectedPatJetsAK5Track.cut = 'pt > 8. & abs(eta) < 2.4'

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

<<<<<<< PM_pat_Layer1_cfg.py

PM_tuple = cms.Sequence( process.patDefaultSequence )

=======
PM_tuple = cms.Sequence( process.recoTrackJets *
                         process.patDefaultSequence )
>>>>>>> 1.9


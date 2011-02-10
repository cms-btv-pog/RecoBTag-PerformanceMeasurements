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

#maxDR=0.3 from BTag note. 
process.patJetPartonMatch.maxDeltaR  = cms.double(0.3) 
process.patJetPartonMatch.maxDPtRel  = cms.double(99999999999)
process.patJetGenJetMatch.maxDeltaR  = cms.double(0.5)
process.patJetGenJetMatch.maxDPtRel  = cms.double(99999999999)

#----------------- Trigger matching ----------------------------------------------------------
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

#----------------- Add other jet collection  ----------------------------------------------------------
addJetCollection(process,
                 cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF',['L2Relative','L3Absolute', 'L2L3Residual']),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel = "ak5"
                 )

theJetNames = ['','AK5PF']

for jetName in theJetNames:
   module = getattr(process,'patJets'+jetName)

   module.tagInfoSources = cms.VInputTag()

   module.discriminatorSources = cms.VInputTag(
    cms.InputTag("jetProbabilityBJetTags"+jetName), 
    cms.InputTag("jetBProbabilityBJetTags"+jetName), 
    cms.InputTag("simpleSecondaryVertexHighEffBJetTags"+jetName),
    cms.InputTag("simpleSecondaryVertexHighPurBJetTags"+jetName),
    cms.InputTag("softMuonBJetTags"+jetName), 
    cms.InputTag("softMuonByPtBJetTags"+jetName), 
    cms.InputTag("softMuonByIP3dBJetTags"+jetName), 
    cms.InputTag("trackCountingHighEffBJetTags"+jetName), 
    cms.InputTag("trackCountingHighPurBJetTags"+jetName)
   )


#-----------------------------------------------------  Pre-selection ----------------------------------------------------------

process.selectedPatMuons.cut= cms.string('pt > 5. & abs(eta) < 2.4 & isGlobalMuon() & innerTrack().numberOfValidHits()> 7')

process.selectedPatMuonsForPtRel= cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string('pt > 5. & abs(eta) < 2.4 & isGlobalMuon() & globalTrack().hitPattern().numberOfValidMuonHits() > 0 & numberOfMatches() > 1 & innerTrack().numberOfValidHits()> 10 & innerTrack().hitPattern().numberOfValidPixelHits()>1 & innerTrack().trackerExpectedHitsOuter().numberOfHits() <3 & innerTrack().normalizedChi2() < 10 & globalTrack().normalizedChi2() < 10 ')
)

process.load("RecoBTag.PerformanceMeasurements.PMConversionFilter_cfi")

process.patElectrons.electronIDSources = cms.PSet(
  softElectronCands = cms.InputTag("softElectronCands")
)

process.selectedPatElectrons.cut= cms.string('pt > 5. && abs(eta) < 2.4 && trackerDrivenSeed() && gsfTrack().numberOfValidHits()> 7 && electronID("softElectronCands")')

process.selectedPatElectronsForS8= cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("PMConversionFilter"),
    cut = cms.string('pt > 5. && abs(eta) < 2.4 && trackerDrivenSeed() && gsfTrack().numberOfValidHits() > 10 && gsfTrack().hitPattern().numberOfValidPixelHits() > 1 && gsfTrack().normalizedChi2() < 10 && electronID("softElectronCands") ')
)

process.selectedPatJets.cut = cms.string('pt > 20. & abs(eta) < 2.4')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 20. & abs(eta) < 2.4')

process.countPatLeptons.minNumber = cms.uint32(1)
 
# Sequence
process.PM_tuple = cms.Sequence(
  process.softElectronCands *
  process.patDefaultSequence *
  (
    process.selectedPatMuonsForPtRel +
    process.PMConversionFilter * 
    process.selectedPatElectronsForS8
  ) 
)

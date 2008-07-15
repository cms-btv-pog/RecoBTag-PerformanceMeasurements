import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("JetMETCorrections.MCJet.MCJetCorrections152_cff")

#  include "RecoBTag/PerformanceMeasurements/data/JetPartonAssoc.cff"
#  sequence jetPartonAssoc = {caloJetCollectionClone, tagJet}
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("RecoBTag.PerformanceMeasurements.PerformanceAnalyzer_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:18A80181-5326-DD11-9498-001617E30D0A.root')
)

process.prefer("MCJetCorrectorIcone5")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('testtt.root')
)

process.p = cms.Path(process.caloJetMCFlavour*process.Performance)



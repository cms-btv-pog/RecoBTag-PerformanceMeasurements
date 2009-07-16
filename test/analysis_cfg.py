import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrectorIC5Calo") 


process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

#process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("RecoBTag.PerformanceMeasurements.PerformanceAnalyzer_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('test-tt.root')
)

process.Performance.flavourMatchOption = 'genParticle'
process.p = cms.Path(process.caloJetMCFlavour*process.Performance)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( (
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0009/9885EAE5-876D-DE11-A0EE-001A9243D640.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/DCAE7AE1-2A6D-DE11-9A56-001A9254460C.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/D2FE0B42-2F6D-DE11-9FF3-001A9227D383.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/C2F0B1DF-2C6D-DE11-8849-001A9243D62A.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/9CC08AE1-2A6D-DE11-96C9-001E8CCCE140.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0005/E893E3BA-246D-DE11-9C26-001A9227D32D.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0005/E04DD0BC-246D-DE11-A774-001E8CCCE148.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0005/CC6C4DC1-246D-DE11-8DDB-00304858A675.root',
        '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0005/C671B2AB-266D-DE11-9849-001A925444DA.root'
)
)

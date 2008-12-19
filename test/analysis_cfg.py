import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrectorIC5Calo") 


process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("RecoBTag.PerformanceMeasurements.PerformanceAnalyzer_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
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
    
            '/store/relval/CMSSW_2_2_1/RelValTTbar/GEN-SIM-RECO/IDEAL_V9_v2/0002/00E9B0FB-98C4-DD11-AF51-0030487A322E.root',
        '/store/relval/CMSSW_2_2_1/RelValTTbar/GEN-SIM-RECO/IDEAL_V9_v2/0002/28E6A7DA-98C4-DD11-AE6E-001D09F2512C.root',
        '/store/relval/CMSSW_2_2_1/RelValTTbar/GEN-SIM-RECO/IDEAL_V9_v2/0002/3802DDF4-9FC4-DD11-A9AC-001D09F24600.root',
        '/store/relval/CMSSW_2_2_1/RelValTTbar/GEN-SIM-RECO/IDEAL_V9_v2/0002/5CA764F6-98C4-DD11-A05F-001D09F241B4.root',
        '/store/relval/CMSSW_2_2_1/RelValTTbar/GEN-SIM-RECO/IDEAL_V9_v2/0002/A6D75ADB-98C4-DD11-B05E-0019B9F704D6.root',
        '/store/relval/CMSSW_2_2_1/RelValTTbar/GEN-SIM-RECO/IDEAL_V9_v2/0003/72C99195-FDC4-DD11-8F03-001617C3B79A.root'
)
)

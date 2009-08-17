import FWCore.ParameterSet.Config as cms

process = cms.Process("OP")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrectorIC5Calo") 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("RecoBTag.PerformanceMeasurements.Taggability_cff")
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
process.load("RecoBTag.PerformanceMeasurements.OperatingPoints")

process.p = cms.Path(process.Taggability+process.caloJetMCFlavour * process.OperatingPoints)

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/E2020AB9-CB3D-DE11-8181-001D09F2546F.root'
                             )
                             )

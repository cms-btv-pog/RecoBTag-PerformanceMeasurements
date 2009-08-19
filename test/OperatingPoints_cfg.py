import FWCore.ParameterSet.Config as cms

process = cms.Process("OP")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrectorIC5Calo") 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
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
'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/E2020AB9-CB3D-DE11-8181-001D09F2546F.root',
'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/8A0F0D07-CB3D-DE11-AF2A-001D09F25109.root',
'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/88D54828-CD3D-DE11-90AB-001D09F24D67.root'
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/2AA92D15-083E-DE11-8DE1-001D09F24DA8.root',
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/08F0EE44-CA3D-DE11-A761-001D09F241D2.root'
)
)

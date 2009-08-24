import FWCore.ParameterSet.Config as cms

process = cms.Process("OperatingPoints")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrectorIC5Calo") 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("RecoBTag.PerformanceMeasurements.Taggability_cff")
process.TFileService = cms.Service("TFileService", fileName = cms.string("taggability.root") )

# study taggability
import RecoBTag.PerformanceMeasurements.Taggability_cff
process.TagL = RecoBTag.PerformanceMeasurements.Taggability_cff.Taggability.clone( MinNtracks = cms.int32(2) )
process.TagM = process.Taggability.clone( MinNtracks = cms.int32(2),
					  MinNPrimaryVertices = cms.int32(1) )

process.TagLsequence = cms.Sequence( process.TagL )
process.TagMsequence = cms.Sequence( process.TagM )

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
process.load("RecoBTag.PerformanceMeasurements.OperatingPoints")
process.defaultsequence = cms.Sequence( process.Taggability+process.caloJetMCFlavour * process.OperatingPoints )

#process.load("RecoBTag.PerformanceMeasurements.plotEff")

process.p = cms.Path(process.defaultsequence + process.TagLsequence)
#process.p = cms.Path( process.defaultsequence + process.TagLsequence + process.TagMsequence )
#process.p = cms.Path( process.defaultsequence + process.plotEff + process.TagLsequence )

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/E2020AB9-CB3D-DE11-8181-001D09F2546F.root',
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/8A0F0D07-CB3D-DE11-AF2A-001D09F25109.root',
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/88D54828-CD3D-DE11-90AB-001D09F24D67.root'
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/2AA92D15-083E-DE11-8DE1-001D09F24DA8.root',
#'/store/relval/CMSSW_2_2_10/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_V12_v1/0003/08F0EE44-CA3D-DE11-A761-001D09F241D2.root'

#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0206/0042F729-F9E6-DD11-8E02-0019B9E4FB7C.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0195/B8482E44-08E6-DD11-ADB6-00145EDD757D.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/ECEBB757-02E6-DD11-8BA6-001D0967D4EF.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/E6616C88-02E6-DD11-876B-0019B9E71424.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/E65C4D2F-02E6-DD11-B022-001D0967DF53.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/ACC561C0-00E6-DD11-8FDE-001D0968F26A.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/A6536F74-02E6-DD11-AF87-001D0967DF53.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/9E16C86F-02E6-DD11-B655-001D0967D63E.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/9ADB801C-FFE5-DD11-98D7-001D0967D193.root',
#'/store/mc/Summer08/QCDpt30/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0194/98A97BA9-00E6-DD11-8E8D-001D0968C880.root'

#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/F8E9ABB1-FFFE-DD11-8481-001D0967D6E3.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/F8BB0B8B-03FF-DD11-9DE6-001D0967D6AC.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/EC5C7EBC-FCFE-DD11-91C4-001D0967D684.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/A8635A3B-FDFE-DD11-8760-001D096B1007.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/A4D0A136-F7FE-DD11-967B-0019B9E48B8C.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/983D59C5-FAFE-DD11-917A-0019B9E713C0.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/8E58F2B7-FEFE-DD11-8A05-0019B9E7B7B3.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/8803F8B8-FDFE-DD11-8CF8-001D0967D288.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/76C135AA-FFFE-DD11-93D8-001D0967C649.root',
#'/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0004/482D904D-FBFE-DD11-9471-0019B9E4F925.root'

'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0007/0C9E75BB-28F4-DD11-9284-00144F0D6802.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/E260E8AB-D3F2-DD11-B02B-001F2908AF72.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/BE60A4AF-D3F2-DD11-84B7-001CC4A6AEF0.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/8CCA14AD-D3F2-DD11-BC37-001A4BD33252.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/8AEC25AE-D3F2-DD11-84ED-001CC47D7BE0.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/46BF28AE-D3F2-DD11-8A95-001CC47D43D4.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/F2DFA5C0-72F2-DD11-BA14-001F29086E20.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/EA9C2BF8-55F2-DD11-9D45-001F29087E54.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/B2B12F81-3CF2-DD11-BABD-001E0B475590.root',
'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/A66C1281-3CF2-DD11-8F69-001CC4C0A4A4.root'

#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/FEF013B3-4C96-DD11-8F2A-00E0813405EC.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/FEDA32FC-4596-DD11-BB79-0019B9CAB9CD.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/FE53BBC6-4C96-DD11-810C-00E08134051C.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/FE156A11-4B96-DD11-AD42-0019B9CB0064.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/FAB51779-4996-DD11-BD03-00E0813405EC.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/FA9E51A1-4996-DD11-9D86-001D091C6771.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/F6CACA21-4B96-DD11-B27C-0019B9CABE48.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/F6C9502B-4496-DD11-8F94-001C23BEDBA9.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/F6884DCB-4C96-DD11-B2D0-00E081340514.root',
#'/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RECO/IDEAL_V9_v1/0004/F64FBEE1-4996-DD11-B97E-001C23BED6A2.root'

)
)

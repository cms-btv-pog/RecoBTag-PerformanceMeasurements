import FWCore.ParameterSet.Config as cms

process = cms.Process("OperatingPoints")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrectorIC5Calo") 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(15000)
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
process.load("RecoBTag.PerformanceMeasurements.OperatingPointsAnalyzer")
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

#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0007/0C9E75BB-28F4-DD11-9284-00144F0D6802.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/E260E8AB-D3F2-DD11-B02B-001F2908AF72.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/BE60A4AF-D3F2-DD11-84B7-001CC4A6AEF0.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/8CCA14AD-D3F2-DD11-BC37-001A4BD33252.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/8AEC25AE-D3F2-DD11-84ED-001CC47D7BE0.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0005/46BF28AE-D3F2-DD11-8A95-001CC47D43D4.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/F2DFA5C0-72F2-DD11-BA14-001F29086E20.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/EA9C2BF8-55F2-DD11-9D45-001F29087E54.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/B2B12F81-3CF2-DD11-BABD-001E0B475590.root',
#'/store/mc/Fall08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0004/A66C1281-3CF2-DD11-8F69-001CC4C0A4A4.root'


#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0013/EA4EAE34-8682-DE11-95C6-001E682F8B72.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0013/CABA2D0E-8782-DE11-BE69-0030487CA7C9.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0013/B4BE0220-8A82-DE11-85F7-00096BB5C092.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0006/FC2E5FCA-EB7F-DE11-BD07-001F2965C268.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0006/E0CA92DA-1880-DE11-B3C3-001E68A994C6.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0006/AC02F1B0-6280-DE11-9667-0030487CD7E8.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0006/6E041725-1980-DE11-94AB-001E682F8B72.root',
#'/store/mc/Summer09/QCD_Pt30/GEN-SIM-RECO/MC_31X_V3-v1/0006/681478A6-F07F-DE11-BFDE-001E0BC18A56.root'

#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0027/A0C39391-EB8C-DE11-8EF8-00144F0D84D8.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0027/8A4E6C25-348B-DE11-A857-0030485C6782.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/FAB061EC-078D-DE11-9473-001E0B470AC2.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/F899160A-088D-DE11-A3AC-001CC4A6FB3A.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/F4C0A1F2-078D-DE11-A192-001E0B46B840.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/E851FF14-088D-DE11-94EA-001CC47B2696.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/E4565A6B-F187-DE11-9510-003048762676.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/D0EE24C6-7488-DE11-8B0C-00144F0D6BEC.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/C8EC5F34-EF87-DE11-AE74-00221987EB87.root',
#'/store/mc/Summer09/InclusiveMu5_Pt50/GEN-SIM-RECO/MC_31X_V3-v1/0022/C663680C-088D-DE11-A1F2-001CC4A7D098.root'

'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0060/F206CF7F-0E9C-DE11-AEEA-001A9243D654.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0060/BA3335B2-0D9C-DE11-B366-001E6878F8BA.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0060/3A532854-889C-DE11-92BD-00E081326D54.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0060/24D90990-0F9C-DE11-8D34-0030483344E2.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0040/FC3FB05A-5399-DE11-A115-00304867FD43.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0040/B4A4DEF2-859C-DE11-A737-00163691DC2A.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0040/642DA79E-6899-DE11-B24E-003048C91B0E.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0040/3C1DFEFF-859C-DE11-9A16-001B243DEF3F.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0039/FEDBAD90-799C-DE11-9913-001EC9AA9978.root',
'/store/mc/Summer09/InclusiveMu5_Pt30/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0039/FC6AE6AF-809C-DE11-8765-0016368E0820.root'

)
)

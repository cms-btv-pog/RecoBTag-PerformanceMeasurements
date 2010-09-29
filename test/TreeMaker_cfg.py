import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_12_1_v2E.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_2_1_uDL.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_24_1_UtW.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_11_1_4yB.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_14_1_c89.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_1_1_Q98.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_27_1_T3N.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_13_1_8GN.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_21_1_uIB.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_6_1_a35.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_10_1_SmT.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_17_1_qkS.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_25_1_iXy.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_7_1_Hi6.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_20_1_mYR.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_19_1_um5.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_15_1_Joa.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_8_1_Pzs.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_5_1_YME.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_22_1_6Hw.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_23_1_RMM.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_26_1_Gml.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_4_1_Djp.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_18_1_9Px.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_16_1_NMC.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_3_1_kKm.root',
        'dcap:///pnfs/cms/WAX/11/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_1/48aa5964db1e96bee7081b34ef88dd25/PM_pattuple_MC_9_1_Esj.root'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TFileService = cms.Service(
    "TFileService",

    fileName = cms.string("s8_tree.root")
)

process.load("RecoBTag.PerformanceMeasurements.TreeMaker_cfi")

process.p = cms.Path(
    process.TreeMaker
)

import FWCore.ParameterSet.printPaths as pp      
pp.printPaths(process)

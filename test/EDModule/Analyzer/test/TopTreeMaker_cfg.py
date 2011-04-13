import FWCore.ParameterSet.Config as cms

process = cms.Process("TopTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/ttmuj/samvel/QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen/Fall10_Syst8_2/f63bb0aca8c14dfe0b90c53588e620e8/PM_pattuple_MC_9_1_0NX.root'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TFileService = cms.Service(
    "TFileService",

    fileName = cms.string("top_tree.root")
)

process.load("EDModule.Analyzer.TopTreeMaker_cfi")

process.load("EDModule.Filter.HLTFilter_cfi")
process.load("EDModule.Filter.PVFilter_cfi")

process.p = cms.Path(
    process.HLTFilter *
    process.PVFilter *
    process.TopTreeMaker
)

import FWCore.ParameterSet.printPaths as pp      
pp.printPaths(process)

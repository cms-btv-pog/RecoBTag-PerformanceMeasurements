import FWCore.ParameterSet.Config as cms

process = cms.Process("S8TreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

for file in open("input.txt").readlines():
    process.source.fileNames.append(file.strip())

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TFileService = cms.Service(
    "TFileService",

    fileName = cms.string("s8_tree.root")
)

process.load("EDModule.Analyzer.S8TreeMaker_cfi")

process.p = cms.Path(
    process.S8TreeMaker
)

import FWCore.ParameterSet.printPaths as pp      
pp.printPaths(process)

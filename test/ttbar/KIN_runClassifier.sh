#!/bin/bash

INPUTFILE="analysis/MC13TeV_TTJets.root"
root -b -q "${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/ttbar/KIN_trainClassifier.C+(\"Fisher,BDT\",\"${INPUTFILE}\",1)"

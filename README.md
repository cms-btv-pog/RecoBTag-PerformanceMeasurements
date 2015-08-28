# RecoBTag-PerformanceMeasurements

**ATTENTION:** The execution of the command
```
git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger-v1_from-CMSSW_7_4_1
```
in the following recipe fixes a bug in the soft-electron tagger in 74X. This bug is not fixed in the official 74X reconstruction, so one is effectively running an algorithm different from what is in the standard reconstruction!

## Software setup

```
cmsrel CMSSW_7_4_8
cd CMSSW_7_4_8/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw

git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger-v1_from-CMSSW_7_4_1
git cms-merge-topic -u cms-btv-pog:CSVv2InCombinedMVA-v1_from-CMSSW_7_4_5

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_25ns.00 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```

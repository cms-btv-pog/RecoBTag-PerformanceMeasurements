# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_7_6_5
cd CMSSW_7_6_5/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw

git cms-merge-topic cms-btv-pog:FixHistoryBase_from-CMSSW_7_6_3
git cms-merge-topic cms-btv-pog:CTagNegPos_765

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_6_X_v1.04 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py miniAOD=False maxEvents=100 reportEvery=1 wantSummary=True
```


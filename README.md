# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_8_0_12
cd CMSSW_8_0_12/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw

git cms-merge-topic -u cms-btv-pog:FixHistoryBase-v1_from-CMSSW_8_0_12
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV3-WithWeightFiles-v1_from-CMSSW_8_0_8_patch1
git cms-merge-topic -u cms-btv-pog:FixBoostedTauConfig_from-CMSSW_8_0_12

git clone -b 8_0_X_v1.06 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py miniAOD=False maxEvents=100 reportEvery=1 wantSummary=True
```


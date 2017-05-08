# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_0_0
cd CMSSW_9_0_0/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git

git clone -b 9_0_X --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

the ntuplizer can be run and configured through 

```
RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py
```


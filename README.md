# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_0_3
cd CMSSW_9_0_3/src
cmsenv

git clone -b 9_0_X_v1.00 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8
```

The ntuplizer can be run and configured through

```
cmsRun RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```


# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_11_1_2 
cd CMSSW_11_1_2/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag
git cms-merge-topic emilbols:BTV_CMSSW_11_1_X

git clone -b Phase2_11_1_X --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```

```
cmsRun runBTagAnalyzer_cfg.py defaults=PhaseII runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```

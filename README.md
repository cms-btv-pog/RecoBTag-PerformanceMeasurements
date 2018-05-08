# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git

#Add DeepFlavour and its dependencies
git cms-merge-topic capalmer85:btagSFupdatesForTTbar
git cms-merge-topic cms-btv-pog:DeepFlavourUpdates-from-CMSSW_9_4_4

git clone -b 9_4_X_v1.06 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for the Moriond 2018 SF campaign

```
runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```


# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_10_1_2_patch2
cd CMSSW_10_1_2_patch2/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git

#Add DeepFlavour and its dependencies
git cms-merge-topic cms-btv-pog:DeepFlavourNewTraining-from-CMSSW_10_1_2_patch2

git clone -b https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for 2018 commissioning

```
runBTagAnalyzer_cfg.py defaults=Commissioning18 runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```


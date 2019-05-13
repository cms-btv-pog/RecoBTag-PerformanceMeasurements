# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_4_13
cd CMSSW_9_4_13/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag/TensorFlow
git cherry-pick 94ceae257f846998c357fcad408986cc8a039152

git clone -b 9_4_X_v1.16 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```. 

To run it for the 2016 Legacy SF campaign, run:

```
cmsRun runBTagAnalyzer_cfg.py defaults=2016_SF runOnData=(True or False, depending on your needs)
```

To run it for 2017 data, including the new JEC, run:

```
cmsRun runBTagAnalyzer_cfg.py defaults=2017NewJEC runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```
The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```

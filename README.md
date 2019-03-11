# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_10_2_11
cd CMSSW_10_2_11/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag
git cms-addpkg PhysicsTools/PatAlgos
git cms-merge-topic rauser:PrunedTraining_NoPuppi_10_2_11
git clone -b PrunedTraining_NoPuppi https://github.com/emilbols/RecoBTag-Combined RecoBTag/Combined/data

git clone -b 10_2_X_v1.05 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for 2018 prompt reco

```
cmsRun runBTagAnalyzer_cfg.py defaults=Prompt18 runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```
The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```

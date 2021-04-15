# RecoBTag-PerformanceMeasurements

## Software setup

```
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_11_2_0_Patatrack
cd CMSSW_11_2_0_Patatrack/src
cmsenv

export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-merge-topic missirol:devel_1120pa_kineParticleFilter -u
git cms-merge-topic missirol:devel_puppiPUProxy_1120patatrack -u
git cms-merge-topic mmasciov:tracking-allPVs -u
git clone https://github.com/missirol/JMETriggerAnalysis.git -o missirol -b run3

git clone -b cleanup-devel --recursive https://github.com/SWuchterl/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j12

```

REMARK: OUTDATED BELOW!!!!

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it with the latest defaults

```
cmsRun runBTagAnalyzer_cfg.py defaults=2017_UltraLegacy runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```
The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```

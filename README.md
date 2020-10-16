# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_10_6_8_patch1 
cd CMSSW_10_6_8_patch1/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git clone -b 10_6_X_UL2018_PreliminaryJECs --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```.

NOTE1: due to the structure of the preliminary JECs, there are 5 different "defaults" sets; 

> **_MC_**: ```defaults=2018_UltraLegacy```<br/>
> **_Run2018A_**: ```defaults=2018_UltraLegacy_DataRunA```<br/>
> **_Run2018B_**: ```defaults=2017_UltraLegacy_DataRunB```<br/>
> **_Run2018C_**: ```defaults=2017_UltraLegacy_DataRunC```<br/>
> **_Run2018D_**: ```defaults=2017_UltraLegacy_DataRunD```

This means in your crab configuration file, you will have to check which file you are running on, and pick the correct default set accordingly. As an example, one might do something like:

```
...
config.Data.inputDataset = /BTagMu/Run2018A-12Nov2019_UL2018-v1/MINIAOD
if "Run2018A" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2018_UltraLegacy_DataRunA ,...]
elif "Run2018B" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2018_UltraLegacy_DataRunB ,...]
elif "Run2018C" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2018_UltraLegacy_DataRunC ,...]
elif "Run2018D" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2018_UltraLegacy_DataRunD ,...]
else: 
	config.JobType.pyCfgParams = [defaults=2018_UltraLegacy ,...]
...
```

NOTE2: The preliminary JECs are only available for AK4PFCHS jets, and therefore you **_can not have any FatJet observables listed in you varGroup!_**. Otherwise the BTA will automatically assume you are running over FatJets and it will use the JECs included in the global tag, rather than in the local SQLite .db files!



To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```
The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```

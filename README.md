# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_10_6_8_patch1 
cd CMSSW_10_6_8_patch1/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git clone -b 10_6_X_v1.05 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```.

NOTE1: due to the structure of the preliminary JECs, there are 6 different "defaults" sets; 

> **_MC_**: ```defaults=2017_UltraLegacy```<br/>
> **_Run2017B_**: ```defaults=2017_UltraLegacy_DataRunB```<br/>
> **_Run2017C_**: ```defaults=2017_UltraLegacy_DataRunC```<br/>
> **_Run2017D_**: ```defaults=2017_UltraLegacy_DataRunD```<br/>
> **_Run2017E_**: ```defaults=2017_UltraLegacy_DataRunE```<br/>
> **_Run2017F_**: ```defaults=2017_UltraLegacy_DataRunF```

This means in your crab configuration file, you will have to check which file you are running on, and pick the correct default set accordingly. As an example, one might do something like:

```
...
config.Data.inputDataset = /BTagMu/Run2017D-09Aug2019_UL2017-v1/MINIAOD
if "Run2017B" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2017_UltraLegacy_DataRunB ,...]
elif "Run2017C" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2017_UltraLegacy_DataRunC ,...]
elif "Run2017D" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2017_UltraLegacy_DataRunD ,...]
elif "Run2017E" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2017_UltraLegacy_DataRunE ,...]
elif "Run2017F" in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2017_UltraLegacy_DataRunF ,...]
else: 
	config.JobType.pyCfgParams = [defaults=2017_UltraLegacy ,...]
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

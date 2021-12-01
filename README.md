# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_12_1_0_pre2
cd CMSSW_12_1_0_pre2/src
cmsenv

export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily

for tcsh
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily

git cms-init

git clone -b 12_1_X --recursive https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j12

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```.

NOTE1: due to the structure of the preliminary JECs, there are 2 different "defaults" sets; 

> **_preVFP_**: ```defaults=2016_UltraLegacy_APV```<br/>
> **_postVFP_**: ```defaults=2018_UltraLegacy_nonAPV```<br/>

This means in your crab configuration file, you will have to check which file you are running on, and pick the correct default set accordingly. As an example, one might do something like:

# Make Target Reference BTagAnalyzer NTuple
```
For MC samples: 
...
config.Data.inputDataset = /*/*RunIISummer20UL16*/MINIAODSIM
...

For collision data samples:
...
config.Data.inputDataset = /BTagMu*/Run2016*-21Feb2020_UL2016*/MINIAOD
...

...
if "HIPM" not in config.Data.inputDataset and "APV" not in config.Data.inputDataset: 
	config.JobType.pyCfgParams = [defaults=2016_UltraLegacy_nonAPV,...]
else: 
	config.JobType.pyCfgParams = [defaults=2016_UltraLegacy_APV, ...]
...
```

NOTE2: The preliminary JECs are only available for AK4PFCHS jets, and therefore you **_can not have any FatJet observables listed in you varGroup!_**. Otherwise the BTA will automatically assume you are running over FatJets and it will use the JECs included in the global tag, rather than in the local SQLite .db files!

# To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```


## BTagAnalyzer General information

The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```


## How to get the latest HLT configuration
For MC:
```
hltGetConfiguration /dev/CMSSW_12_0_0/GRun/V2 \
--full \
--offline \
--unprescale \
--mc \
--process HLT2 \
--globaltag auto:phase1_2021_realistic \
--max-events 10 \
> tmp.py
```
```
edmConfigDump tmp.py > HLT_dev_CMSSW_12_0_0_GRun_V3_configDump_MC.py
```
For data:
```
hltGetConfiguration /dev/CMSSW_12_0_0/GRun \
--full \
--offline \
--unprescale \
--data \
--process HLT2 \
--globaltag auto:run3_hlt \
--max-events 10 \
--customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2018Input \
> tmp_data.py
```
```
edmConfigDump tmp_data.py > HLT_dev_CMSSW_12_0_0_GRun_V3_configDump_Data.py
```

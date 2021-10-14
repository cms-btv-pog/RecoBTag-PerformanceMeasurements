# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_12_1_0_pre2
cd CMSSW_12_1_0_pre2/src
cmsenv

export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git clone -b 12_1_X --recursive https://github.com/johnalison/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements
git cms-merge-topic johnalison:MakePy3Bind11ParameterSetsIncludingCommandLineArguments
git clone -b BTVDQM git@github.com:patrickbryant/nTupleAnalysis.git
git clone -b 12_1_X git@github.com:johnalison/TriggerStudies.git

scram b -j12

```


# Make Target Reference BTagAnalyzer NTuple
cp RecoBTag/PerformanceMeasurements/python/defaults/Run3.py RecoBTag/PerformanceMeasurements/python/defaults/Run3Reference.py 
[ edit Run3Refernce to point to the right file ]
cd RecoBTag/PerformanceMeasurements/test/
cmsRun runBTagAnalyzer_cfg.py defaults=Run3Reference runOnData=(True or False, depending on your needs) 
cd -

# Make Target Reference BTagAnalyzer NTuple
cp RecoBTag/PerformanceMeasurements/python/defaults/Run3.py RecoBTag/PerformanceMeasurements/python/defaults/Run3Target.py 
[ edit Run3Target to point to the right file ]
cd RecoBTag/PerformanceMeasurements/test/
cmsRun runBTagAnalyzer_cfg.py defaults=Run3Target runOnData=(True or False, depending on your needs) 
cd -

### Process Ntuples To make ROOT file with histograms

```
BTagAnalyzer TriggerStudies/NtupleAna/scripts/BTagAnalyzer_cfg.py \
    --inputRAW  path/to/target/BTagAnalyzer/Ntuple/JetTree.root \
    --inputAOD path/to/reference/BTagAnalyzer/Ntuple/JetTree.root \
    -o $PWD \
    -y Run3 \
    --histogramming 1 \
    --histFile hists_Run3_offlineValidation.root \
    --isMC \
    --doTracks \
    --pfJetName "" \ 
    --nevents -1
```



### Make plots presentation
``` 
git clone -b 12_1_X  git@github.com:johnalison/ROOTHelp.git
source ROOTHelp/setup.sh
source TriggerStudies/plotting/offlineDQMPlots.sh hists_Run3_offlineValidation.root hists_Run3_offlineValidation targetName
```
Outpus a slides in hists_Run3_offlineValidation directory




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

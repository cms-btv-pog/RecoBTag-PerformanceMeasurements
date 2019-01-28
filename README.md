# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_4_12
cd CMSSW_9_4_12/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag/TensorFlow
git cherry-pick 94ceae257f846998c357fcad408986cc8a039152
git cms-addpkg RecoBTag/Combined
wget https://github.com/daseith/RecoBTag-Combined/raw/master/DeepFlavourPhaseIIV01_94X_training.pb -P RecoBTag/Combined/data/DeepFlavourPhaseIIV01/
sed 's/DeepFlavourV03_10X_training\/constant_graph.pb/DeepFlavourPhaseIIV01\/DeepFlavourPhaseIIV01_94X_training.pb/' RecoBTag/TensorFlow/plugins/DeepFlavourTFJetTagsProducer.cc -i
sed 's/\"cpf_input_batchnorm\/keras_learning_phase\"//g' RecoBTag/TensorFlow/plugins/DeepFlavourTFJetTagsProducer.cc -i


git clone -b 9_4_X_PhaseII_v1.04 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for PhaseII studies

```
cmsRun runBTagAnalyzer_cfg.py defaults=PhaseII runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```
The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```

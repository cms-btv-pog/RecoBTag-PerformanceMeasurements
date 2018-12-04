# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_10_1_8
cd CMSSW_10_1_8/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git

# Add DeepDoubleB tagger and updated DeepFlavour and DeepDoubleB training
git cms-merge-topic cms-btv-pog:DeepDoubleB_from-CMSSW_10_1_8
git cms-merge-topic cms-btv-pog:DeepFlavourNewTraining-from-CMSSW_10_1_2_patch2
mkdir -p RecoBTag/Combined/data/DeepDoubleB/V01/
wget https://github.com/cms-data/RecoBTag-Combined/raw/V01-00-13/DeepDoubleB/V01/constant_graph_PtCut.pb -P RecoBTag/Combined/data/DeepDoubleB/V01/
wget https://github.com/cms-data/RecoBTag-Combined/raw/V01-00-13/DeepDoubleB/V01/constant_graph_PtCut_MassSculptPen.pb -P RecoBTag/Combined/data/DeepDoubleB/V01/

# OPTIONAL: Add new DeepFlavour training for PhaseII studies
wget https://github.com/daseith/RecoBTag-Combined/raw/master/DeepFlavourPhaseIIV01_94X_training.pb -P RecoBTag/Combined/data/DeepFlavourPhaseIIV01/
sed 's/DeepFlavourV02_realistic_training\/constant_graph.pb/DeepFlavourPhaseIIV01\/DeepFlavourPhaseIIV01_94X_training.pb/' RecoBTag/TensorFlow/plugins/DeepFlavourTFJetTagsProducer.cc -i
sed 's/\"cpf_input_batchnorm\/keras_learning_phase\"//g' RecoBTag/TensorFlow/plugins/DeepFlavourTFJetTagsProducer.cc -i


git clone -b 10_1_X_v1.03 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for 2018 commissioning

```
cmsRun runBTagAnalyzer_cfg.py defaults=Commissioning18 runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```
The content of the output ntuple is by default empty and has to be configured according to your needs. The ```store*Variables``` options have been removed.
The new variable configuration can be customized in the file ```RecoBTag/PerformanceMeasurements/python/varGroups_cfi.py```.
New variables need also to be added (apart from adding them in the code) in ```RecoBTag/PerformanceMeasurements/python/variables_cfi.py```

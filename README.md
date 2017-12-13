# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_4_1
cd CMSSW_9_4_1/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git

#Add DeepFlavour and its dependencies
git cms-addpkg DataFormats/BTauReco
git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg RecoBTag/Combined
git cms-addpkg RecoBTag/Configuration
git cms-merge-topic pablodecm:DeepFlavour_9_4_1_backport

git clone -b 9_4_X --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for the Moriond 2018 SF campaign

```
runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=(True or False, depending on your needs)
```



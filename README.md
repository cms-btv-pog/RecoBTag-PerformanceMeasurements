# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_8_0_23
cd CMSSW_8_0_23/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-merge-topic 16377
git cms-merge-topic mverzett:DeepFlavour-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
cd -

git clone -b 8_0_X_v2.01 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```

Note: if you run on ICHEP MC you need to set by hand the relaxed BTV track selection, setting the options remakeAllDiscr=True and runIVF=True


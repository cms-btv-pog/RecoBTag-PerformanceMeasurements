# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_8_0_23
cd CMSSW_8_0_23/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw

git cms-merge-topic 16377
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
cd -

git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21

git clone -b 8_0_X_v2.04 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```
Note: the "-u" option in git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21, prevents to download additional 61 packages that have some dependences on the changes introduced in mverzett:DeepFlavour-from-CMSSW_8_0_21. However, if you want to use additional packages in your installation, for some reason, the code may crash. In that case it's safer to remove the -u option and ask: git cms-merge-topic mverzett:DeepFlavour-from-CMSSW_8_0_21

Note: if you run on ICHEP MC you need to set by hand the relaxed BTV track selection, setting the options remakeAllDiscr=True and runIVF=True


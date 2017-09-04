# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_8_0_28_patch1
cd CMSSW_8_0_28_patch1/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw
```

DeeFlavour taggers to be implemented following instructions from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepJet
```
git cms-merge-topic -u mverzett:Experimental_DeepFlavour_80X

cd RecoBTag/DeepFlavour/scripts/

./setup_legacy.sh

cd ../../..

git clone -b 8_0_X --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```

OLD note for the BTagAnalyzer 8_0_X_v2.05: the "-u" option in git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21, prevents to download additional 61 packages that have some dependences on the changes introduced in mverzett:DeepFlavour-from-CMSSW_8_0_21. However, if you want to use additional packages in your installation, for some reason, the code may crash. In that case it's safer to remove the -u option and ask: git cms-merge-topic mverzett:DeepFlavour-from-CMSSW_8_0_21

Note: if you run on ICHEP MC you need to set by hand the relaxed BTV track selection, setting the options remakeAllDiscr=True and runIVF=True


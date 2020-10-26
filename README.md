# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_11_1_0_pre6
cd CMSSW_11_1_0_pre6/src
cmsenv

export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag
git cms-addpkg RecoBTag/TensorFlow
git cms-addpkg RecoBTag/Combined

git clone -b PrunedTraining_NoPuppi https://github.com/emilbols/RecoBTag-Combined RecoBTag/Combined/data
git cms-merge-topic emilbols:BTV_CMSSW_11_1_X

git clone -b PhaseIIOffline --depth 1 https://github.com/johnalison/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

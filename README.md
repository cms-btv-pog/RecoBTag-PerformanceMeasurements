# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_11_0_0
cd CMSSW_11_0_0/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag/TensorFlow
git cms-addpkg RecoBTag/Combined
wget https://raw.githubusercontent.com/cms-data/RecoBTag-Combined/master/DeepCSV_PhaseII.json -P RecoBTag/Combined/data/

git clone -b PhaseIIOffline --depth 1 https://github.com/johnalison/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

git cms-merge-topic missirol:devel_pixvtx_buildfile_1100pre13

scram b -j8

```

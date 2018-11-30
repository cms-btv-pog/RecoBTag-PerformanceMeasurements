# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_4_9
cd CMSSW_9_4_9/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag/TensorFlow
git cherry-pick 94ceae257f846998c357fcad408986cc8a039152

git clone -b 9_4_X_v1.10 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

```

The ntuplizer can be run and configured through ```RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py```, to run it for the Summer 2018 SF campaign

```
cmsRun runBTagAnalyzer_cfg.py defaults=Summer18 runOnData=(True or False, depending on your needs)
```

To run the tests for integrating changes run:

```
cd RecoBTag/PerformanceMeasurements/test/
./run_tests.sh
```

To run DeepAK8:
  * The algorithm is merged in CMSSW_9_4_11_cand1.
  * For previous releases:
~~~
# setup MXNet and DeepAK8 models
scram setup /cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_11_cand1/config/toolbox/slc6_amd64_gcc630/tools/selected/mxnet-predict.xml
# get DeepAK8 PR for 94X
git cms-merge-topic -u hqucms:deep-boosted-jets-rebase-94X
# setup the data files
mkdir -p $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoBTag/Combined/data
cp -r /cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_11_cand1/external/slc6_amd64_gcc630/data/RecoBTag/Combined/data/DeepBoostedJet $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoBTag/Combined/data
# compile
scram b -j 24
~~~

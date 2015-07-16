# RecoBTag-PerformanceMeasurements
ATTENTION: the excution of the command "git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger_from-CMSSW_7_4_5" in the following recipe fixes a bug in the soft-electron tagger in 74X. This bug is not fixed in the official 74X reconstruction, so one is effectively running an algorithm different from what is in the standard reconstruction! 
=========================
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src
cmsenv
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch btv-cmssw
git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger-v1_from-CMSSW_7_4_1
git cms-merge-topic -u cms-btv-pog:CSVv2InCombinedMVA-v1_from-CMSSW_7_4_5

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_v1.09 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py

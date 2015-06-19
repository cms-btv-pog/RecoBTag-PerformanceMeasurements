# RecoBTag-PerformanceMeasurements
ATTENTION: the excution of the command "git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger_from-CMSSW_7_4_5" in the following recipe fixes a bug in the soft-electron tagger in 74X. This bug is not fixed in the official 74X reconstruction, so one is effectively running an algorithm different from what is in the standard reconstruction! 
=========================
cmsrel CMSSW_7_4_5
cd CMSSW_7_4_5/src
cmsenv
git cms-init
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTagger-WithWeightFiles-v2_from-CMSSW_7_4_1
git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger_from-CMSSW_7_4_5

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_v1.06 git@github.com:cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py

# RecoBTag-PerformanceMeasurements
=========================
cmsrel CMSSW_7_4_0_pre7
cd CMSSW_7_4_0_pre7/src
cmsenv

git cms-merge-topic -u cms-btv-pog:PATBTaggingUpdates_from-CMSSW_7_4_0_pre7

git cms-merge-topic -u cms-btv-pog:TrueTrackQualityFix_from-CMSSW_7_4_0_pre7

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_dev_full git@github.com:cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py

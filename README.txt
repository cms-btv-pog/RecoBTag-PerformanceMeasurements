# RecoBTag-PerformanceMeasurements
=========================
cmsrel CMSSW_7_4_1
cd CMSSW_7_4_1/src
cmsenv

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X git@github.com:cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

scram b -j8

cd RecoBTag/PerformanceMeasurements/test/

cmsRun runBTagAnalyzer_cfg.py

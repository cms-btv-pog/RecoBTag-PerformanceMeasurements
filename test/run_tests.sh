#! /bin/bash

set -o errexit
function die { echo 'FAILED': status $1 ; }

echo 'Running on data APV'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_APV runOnData=True maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on MC APV'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_APV runOnData=False maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on data APV -- ctag'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_APV runOnData=True maxEvents=20 runCTagVariables=True runFatJets=True runSubJets=True groups='testfat' &> last.log || die $?

echo 'Running on MC APV -- ctag'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_APV runOnData=False maxEvents=20 runCTagVariables=True runFatJets=True runSubJets=True groups='testfat' &> last.log || die $?

echo 'Running on data APV -- FatJets'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_APV runOnData=True maxEvents=20 runFatJets=True groups='testfat' &> last.log || die $?

echo 'Running on MC APV -- FatJets'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_APV runOnData=False maxEvents=20 runFatJets=True groups='testfat' &> last.log || die $?

echo 'Running on data nonAPV'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_nonAPV runOnData=True maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on MC nonAPV'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_nonAPV runOnData=False maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on data nonAPV -- ctag'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_nonAPV runOnData=True maxEvents=20 runCTagVariables=True runFatJets=True runSubJets=True groups='testfat' &> last.log || die $?

echo 'Running on MC nonAPV -- ctag'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_nonAPV runOnData=False maxEvents=20 runCTagVariables=True runFatJets=True runSubJets=True groups='testfat' &> last.log || die $?

echo 'Running on data nonAPV -- FatJets'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_nonAPV runOnData=True maxEvents=20 runFatJets=True groups='testfat' &> last.log || die $?

echo 'Running on MC nonAPV -- FatJets'
cmsRun runBTagAnalyzer_cfg.py defaults=2016_UltraLegacy_nonAPV runOnData=False maxEvents=20 runFatJets=True groups='testfat' &> last.log || die $?


# echo 'Running on AODSIM'
# cmsRun runBTagAnalyzer_cfg.py defaults=2017_UltraLegacy runOnData=False maxEvents=20 miniAOD=False inputFiles=/store/mc/RunIIAutumn18DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/1010000/FA97679F-774B-7F43-9258-D8A0AE3A7A01.root useSelectedTracks=False produceJetTrackTree=True runCTagVariables=False groups='test' &> last.log || die $?


#echo 'Running on data -- commissioning'
#cmsRun runBTagAnalyzer_cfg.py defaults=2018_UltraLegacy useSelectedTracks=False produceJetTrackTree=True runOnData=True maxEvents=20  groups='test'&> last.log || die $?

#echo 'Running on MC -- commissioning'
#cmsRun runBTagAnalyzer_cfg.py defaults=2018_UltraLegacy useSelectedTracks=False produceJetTrackTree=True runOnData=False maxEvents=20 groups='test' &> last.log || die $?


# echo 'RECODEBUG'
# cmsRun runBTagAnalyzer_cfg.py defaults=2017_UltraLegacy runOnData=False maxEvents=20 miniAOD=False runFatJets=True runSubJets=True useTrackHistory=True produceJetTrackTree=True inputFiles=/store/mc/RunIIFall17DRPremix/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/70000/04F991F0-C2DD-E711-B302-0CC47A78A3F8.root groups='testfat' &> last.log || die $?

#echo 'Running on MINIAODSIM boosted'
#cmsRun runBTagAnalyzer_cfg.py defaults=Moriond19Boosted runOnData=False maxEvents=20 &> boosted.log || die $? 

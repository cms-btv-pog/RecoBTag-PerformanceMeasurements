#! /bin/bash

# set -o errexit
function die { echo 'FAILED': status $1 ;  exit $1; }

echo 'Running on data'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=True maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on MC'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=False maxEvents=20 groups='test' &> last.log || die $?

echo 'Rinning on AODSIM'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=False maxEvents=20 miniAOD=False inputFiles=/store/mc/RunIIFall17DRPremix/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/AODSIM/94X_mc2017_realistic_v10-v1/60000/2AE14B27-6FE7-E711-B564-00266CF94C44.root useSelectedTracks=False produceJetTrackTree=True runCTagVariables=False groups='test' &> last.log || die $?

echo 'Running on data -- FatJets'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=True maxEvents=20 runFatJets=True groups='testfat' &> last.log || die $?

echo 'Running on MC -- FatJets'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=False maxEvents=20 runFatJets=True groups='testfat' &> last.log || die $?

echo 'Running on data -- commissioning'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 useSelectedTracks=False produceJetTrackTree=True runOnData=True maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on MC -- commissioning'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 useSelectedTracks=False produceJetTrackTree=True runOnData=False maxEvents=20 groups='test' &> last.log || die $?

echo 'Running on data -- ctag'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=True maxEvents=20 runCTagVariables=True runFatJets=True runSubJets=True groups='testfat' &> last.log || die $?

echo 'Running on MC -- ctag'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=False maxEvents=20 runCTagVariables=True runFatJets=True runSubJets=True groups='testfat' &> last.log || die $?

echo 'Running on MINIAODSIM boosted'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond19Boosted runOnData=False maxEvents=20 &> boosted.log || die $? 

echo 'RECODEBUG'
cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=False maxEvents=20 miniAOD=False runFatJets=True runSubJets=True useTrackHistory=True produceJetTrackTree=True inputFiles=/store/mc/RunIIFall17DRPremix/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/70000/04F991F0-C2DD-E711-B302-0CC47A78A3F8.root groups='testfat' &> last.log || die $?





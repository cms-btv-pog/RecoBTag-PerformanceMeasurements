#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_GRun \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_GRun_v0
#
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRK \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_Run3TRK_v0
# #
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKWithPU \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_Run3TRKwithPU_v0

# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKMod \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_Run3TRKMod_v0
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKMod2 \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_Run3TRKMod2_v0
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnly \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlySorted_v0
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnlyCleaned \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlyCleaned_v0

# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnlyCleaned2 \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlyCleaned2_v1
#
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnlyCleaned3 \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlyCleaned3_v1
#
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnlyCleaned4 \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlyCleaned4_v1

# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnlyCleaned2 \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlyCleaned2Mod2Pt_v1

# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKPixelOnlyCleaned2 \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_HLT_Run3TRKPixelOnlyCleaned2Mod3OnlyPt08_v1
python ../submit_all.py \
  ../runHLTBTagAnalyzer_cfg.py \
  -f tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_BTagROI \
  -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
  -v crab_projects_Run3Study_HLT_BTagROI_v1

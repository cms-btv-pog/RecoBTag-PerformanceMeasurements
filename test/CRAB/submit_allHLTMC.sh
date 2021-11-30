#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=False reco=HLT_GRun  \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_120X_GRun_v0
# #
# python ../submit_all.py \
#   ../runHLTBTagAnalyzer_cfg.py \
#   -f tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=False reco=HLT_Run3TRK  \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_120X_Run3TRK_v0

python ../submit_all.py \
  ../runHLTBTagAnalyzer_cfg.py \
  -f tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=False reco=HLT_Run3TRKPixelOnlyCleaned2  \
  -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
  -v crab_projects_Run3Study_120X_HLT_Run3TRKPixelOnlyCleaned2_v0

python ../submit_all.py \
  ../runHLTBTagAnalyzer_cfg.py \
  -f tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=False reco=HLT_Run3TRKForBTag  \
  -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
  -v crab_projects_Run3Study_120X_HLT_Run3TRKForBTag_v0

python ../submit_all.py \
  ../runHLTBTagAnalyzer_cfg.py \
  -f tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=False reco=HLT_Run3TRKForBTag_2  \
  -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
  -v crab_projects_Run3Study_120X_HLT_Run3TRKForBTag_2_v1

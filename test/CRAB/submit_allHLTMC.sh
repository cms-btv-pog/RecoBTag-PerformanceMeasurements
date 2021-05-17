#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
python ../submit_all.py \
  ../runHLTBTagAnalyzer_cfg.py \
  -f tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_GRun \
  -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
  -v crab_projects_Run3Study_GRun_v0

# python submit_all.py \
#   runHLTBTagAnalyzer_cfg.py \
#   -f CRAB/tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRK \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_Run3TRK_v0
#
# python submit_all.py \
#   runHLTBTagAnalyzer_cfg.py \
#   -f CRAB/tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=Run3 runOnData=False runCaloJetVariables=False runPuppiJetVariables=True reco=HLT_Run3TRKWithPU \
#   -o /store/user/sewuchte/BTagServiceWork/Run3/Online/ \
#   -v crab_projects_Run3Study_Run3TRKwithPU_v0

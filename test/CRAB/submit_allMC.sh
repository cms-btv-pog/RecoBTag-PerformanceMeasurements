#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
python ../submit_all.py \
  ../runBTagAnalyzer_cfg.py \
  -f tosubmit_miniaod.txt \
  -s T2_DE_DESY \
  -p defaults=Run3 runOnData=False \
  -o /store/user/sewuchte/BTagServiceWork/Run3/Offline/ \
  -v crab_projects_Run3Study_120X_Offline_v0

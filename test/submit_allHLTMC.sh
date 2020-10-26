#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
# python submit_all.py \
#   runBTagAnalyzer_cfg.py \
#   -f CRAB/tosubmit.txt \
#   -s T2_DE_DESY \
#   -p defaults=PhaseII runOnData=False \
#   -o /store/user/sewuchte/BTagServiceWork/PhaseII/Offline/ \
#   -v crab_projects_HLTTDR_v4

python submit_all.py \
  runBTagAnalyzer_cfg.py \
  -f CRAB/tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=PhaseII_puppi runOnData=False \
  -o /store/user/sewuchte/BTagServiceWork/PhaseII/Offline/ \
  -v crab_projects_SeptemberL1_v4

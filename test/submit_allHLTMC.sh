
#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
python submit_all.py \
  runBTagAnalyzer_cfg.py \
  -f CRAB/tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=PhaseII runOnData=False  \
  -o /store/user/sewuchte/BTagServiceWork/PhaseII/Offline/ \
  -v NewGT_v1_TT

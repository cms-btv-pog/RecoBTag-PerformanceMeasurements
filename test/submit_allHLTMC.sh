
#!/bin/bash

echo "!!!! WARNING: Submitting for MC!!!"
python submit_all.py \
  runBTagAnalyzer_cfg.py \
  -f CRAB/tosubmit.txt \
  -s T2_DE_DESY \
  -p defaults=PhaseII_CMSSW_11_0_0 runOnData=False runCMSSW11Sample=True \
  -o /store/user/sewuchte/BTagServiceWork/PhaseII/Offline/ \
  -v crab_NewOfflineSequence_CMSSW11_0_0_TrackingV0_1

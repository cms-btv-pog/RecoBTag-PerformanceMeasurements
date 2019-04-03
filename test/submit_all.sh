#!/bin/bash

###echo "!!!! WARNING: Submitting for Data!!!!"
echo "!!!! WARNING: Submitting for MC!!!"
python submit_all.py \
  runBTagAnalyzer_cfg.py \
  -f CRAB/tosubmit.txt \
  -s T2_CH_CERN \
  -p defaults=Moriond19Boosted runOnData=False \
  -o /store/group/phys_btag/BoostedBTag/BTagNTuples/2018/ \
  -v 10_2_X_9bc62860


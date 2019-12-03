#!/bin/bash
# Discriminant for scale factors:
tagger=DeepCSVBDisc #DeepFlavourBDisc #DeepCSVBDisc
# Fit directory
folders=(UltraLegacy_fit_dir)
#UltraLegacy_fit_dir
#UltraLegacy_fit_dir_TuneCP5down
#UltraLegacy_fit_dir_TuneCP5up
#UltraLegacy_fit_dir_hdampDOWN
#UltraLegacy_fit_dir_hdampUP
#UltraLegacy_fit_dir_mtop171p5
#UltraLegacy_fit_dir_mtop173p5
#UltraLegacy_fit_dir_nonTTXSecDown

# Loop over chosen folders above.
for f in ${folders[@]}; do
  python createSFbSummaryReport.py -i "kin":$f/kindisc_templates/.$tagger\_fits.pck -o $folders\/$tagger\_fits;
done

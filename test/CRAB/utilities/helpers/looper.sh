#!/bin/bash

SAMPLES=(

 crab_*
)

for idx in {1..35}; do

  sleep 1800

  crab3_monitor.py -t ${SAMPLES[*]} -r

done
unset -v idx

unset -v SAMPLES

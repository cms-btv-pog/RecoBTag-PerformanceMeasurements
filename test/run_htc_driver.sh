#!/bin/bash

set -e

if [ $# -lt 1 ]; then
  echo ">> argument missing - specify path to output directory"
  exit 1
fi

IDIR=das_jsondumps
ODIR=$1

DSETS=(
 # QCD_Pt_15to3000_Flat_14TeV_NoPU
 # QCD_Pt_15to3000_Flat_14TeV_PU140
 # QCD_Pt_15to3000_Flat_14TeV_PU200

# TT_14TeV_NoPU
# TT_14TeV_PU140
TT_14TeV_PU200

 # VBF_HToInvisible_M125_14TeV_NoPU
 # VBF_HToInvisible_M125_14TeV_PU140
 # VBF_HToInvisible_M125_14TeV_PU200
)

EXE="htc_driver -c runHLTBTagAnalyzer_PhaseII_cfg.py -n 350 numThreads=1 --cpus 1 --memory 3000 --jobflavour workday"

for dset in "${DSETS[@]}"; do

  if [ ! -f ${IDIR}/${dset}.json ]; then
    echo "input file does not exist: ${IDIR}/${dset}.json"
    continue
  fi

  if [ -d ${ODIR}/${dset} ]; then
    echo "output directory already exists: ${ODIR}/${dset}"
    continue

  elif [ ! -d ${ODIR} ]; then
    mkdir -p ${ODIR}
  fi

  ${EXE} -d ${IDIR}/${dset}.json -o ${ODIR}/${dset} "${@:2}"

done
unset -v dset

unset -v IDIR ODIR DSETS EXE

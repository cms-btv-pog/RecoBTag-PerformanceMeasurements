#!/bin/bash

set -e

ODIR=das_jsondumps

if [ -d ${ODIR} ]; then
  echo "output directory already exists: ${ODIR}"
  exit 1
fi

mkdir -p ${ODIR}

# das_jsondump -v -o ${ODIR}/QCD_Pt_15to3000_Flat_14TeV_NoPU.json  -d /QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/PhaseIITDRSpring19MiniAOD-NoPU_castor_106X_upgrade2023_realistic_v3-v2/MINIAODSIM
# das_jsondump -v -o ${ODIR}/QCD_Pt_15to3000_Flat_14TeV_PU140.json -d /QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/PhaseIITDRSpring19MiniAOD-PU140_castor_106X_upgrade2023_realistic_v3-v2/MINIAODSIM
# das_jsondump -v -o ${ODIR}/QCD_Pt_15to3000_Flat_14TeV_PU200.json -d /QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/PhaseIITDRSpring19MiniAOD-PU200_castor_106X_upgrade2023_realistic_v3-v2/MINIAODSIM

# das_jsondump -v -o ${ODIR}/TT_14TeV_NoPU.json  -d /TT_TuneCP5_14TeV-powheg-pythia8/PhaseIITDRSpring19MiniAOD-NoPU_106X_upgrade2023_realistic_v3-v2/MINIAODSIM
# das_jsondump -v -o ${ODIR}/TT_14TeV_PU140.json -d /TT_TuneCP5_14TeV-powheg-pythia8/PhaseIITDRSpring19MiniAOD-PU140_106X_upgrade2023_realistic_v3-v1/MINIAODSIM
das_jsondump -v -o ${ODIR}/TT_14TeV_PU200.json -d /TT_TuneCP5_14TeV-powheg-pythia8/PhaseIITDRSpring19MiniAOD-PU200_106X_upgrade2023_realistic_v3-v1/MINIAODSIM
# das_jsondump -v -o ${ODIR}/TT_14TeV_PU200.json -d /TT_TuneCP5_14TeV-powheg-pythia8/PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v1/GEN-SIM-DIGI-RAW

# das_jsondump -v -o ${ODIR}/VBF_HToInvisible_M125_14TeV_NoPU.json  -d /VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRSpring19MiniAOD-NoPU_106X_upgrade2023_realistic_v3-v2/MINIAODSIM
# das_jsondump -v -o ${ODIR}/VBF_HToInvisible_M125_14TeV_PU140.json -d /VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRSpring19MiniAOD-PU140_106X_upgrade2023_realistic_v3-v1/MINIAODSIM
# das_jsondump -v -o ${ODIR}/VBF_HToInvisible_M125_14TeV_PU200.json -d /VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRSpring19MiniAOD-PU200_106X_upgrade2023_realistic_v3-v1/MINIAODSIM

unset -v ODIR

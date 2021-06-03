#!/bin/bash

set -e

if [ $# -lt 1 ]; then
  printf "\n%s\n\n" ">> argument missing - specify path to output directory"
  exit 1
fi

NEVT=-1

if [ $# -eq 1 ]; then
  ODIR=${1}
  ODIR_cmsRun=$1
else
  ODIR=${1}
  ODIR_cmsRun=${2}
fi

if [ -d ${ODIR} ]; then
  printf "%s\n" "output directory already exists: ${ODIR}"
  exit 1
fi

declare -A samplesMap

samplesMap["RelValTTbar_14TeV"]="/RelValTTbar_14TeV/CMSSW_11_2_0_pre9-PU_112X_mcRun3_2021_realistic_v11-v1/GEN-SIM-DIGI-RAW"

recoKeys=(
    HLT_Run3TRKPixelOnlyCleaned
)

# options (JobFlavour and AccountingGroup)
opts=""
if [[ ${HOSTNAME} == lxplus* ]]; then
  opts+="--JobFlavour longlunch"
  # if [[ ${USER} == sewuchte ]]; then
    # opts+=" --AccountingGroup group_u_CMS.CAF.PHYS"
  # fi
fi
    for recoKey in "${recoKeys[@]}"; do
        for sampleKey in ${!samplesMap[@]}; do
            sampleName=${samplesMap[${sampleKey}]}
            numEvents=${NEVT}
            python ../runHLTBTagAnalyzer_cfg.py dumpPython=/tmp/${USER}/BTV_Run3_Online_cfg.py defaults=Run3 runOnData=False reco=${recoKey} runPuppiJetVariables=True runCaloJetVariables=False

            htc_driver -c /tmp/${USER}/BTV_Run3_Online_cfg.py --customize-cfg -m ${numEvents} -n 1000 --cpus 1 --memory 2000 --runtime 10800 ${opts} \
            -d ${sampleName} -p 0 -o ${ODIR}/${sampleKey} --cmsRun-output-dir ${ODIR_cmsRun}/${sampleKey}
        done
    done
    unset sampleKey sampleName numEvents
unset opts samplesMap NEVT ODIR ODIR_cmsRun recoKey

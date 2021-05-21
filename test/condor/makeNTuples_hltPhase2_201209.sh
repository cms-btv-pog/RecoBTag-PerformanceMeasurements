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

samplesMap["Phase2HLTTDR_TTbar_14TeV_PU140"]="/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"

recoKeys=(
    HLT_TRKv00_TICL
)
# /eos/home-s/sewuchte/BTV-Phase2/December_TDR/Condor_Prod_v01

# options (JobFlavour and AccountingGroup)
opts=""
if [[ ${HOSTNAME} == lxplus* ]]; then
  opts+="--JobFlavour longlunch"
  # if [[ ${USER} == sewuchte ]]; then
    # opts+=" --AccountingGroup group_u_CMS.CAF.PHYS"
  # fi
fi

# for recoKey in "${recoKeys[@]}"; do
        for sampleKey in ${!samplesMap[@]}; do
            sampleName=${samplesMap[${sampleKey}]}
            numEvents=${NEVT}
            # python ../runHLTBTagAnalyzer_PhaseII_cfg.py dumpPython=/tmp/${USER}/BTV_PhaseII_Offline_cfg.py numThreads=1 reco=${recoKey} globalTag=111X_mcRun4_realistic_T15_v2
            python ../runHLTBTagAnalyzer_cfg.py dumpPython=/tmp/${USER}/BTV_Run3_Online_cfg.py defaults=Run3 runOnData=False reco=${recoKey} runPuppiJetVariables=True runCaloJetVariables=False

            # htc_driver -c /tmp/${USER}/BTV_PhaseII_Offline_cfg.py --customize-cfg -m ${numEvents} -n 250 --cpus 1 --memory 2000 --runtime 10800 ${opts} \
            htc_driver -c /tmp/${USER}/BTV_PhaseII_Offline_cfg.py --customize-cfg -m ${numEvents} -n 1000 --cpus 1 --memory 2000 --runtime 10800 ${opts} \
            -d ${sampleName} -p 0 -o ${ODIR}/${sampleKey} --cmsRun-output-dir ${ODIR_cmsRun}/${sampleKey}
        done
    # done
    unset sampleKey sampleName numEvents
# done
# unset recoKey opts recoKeys samplesMap NEVT ODIR ODIR_cmsRun
unset opts samplesMap NEVT ODIR ODIR_cmsRun

#!/bin/bash

set -e

if [ $# -lt 1 ]; then
  printf "\n%s\n\n" ">> argument missing - specify path to output directory"
  exit 1
fi

NEVT=-1

if [ $# -eq 1 ]; then
  ODIR=${1}
  # ODIR_cmsRun=$1
else
  ODIR=${1}
  # ODIR_cmsRun=${2}
fi

if [ -d ${ODIR} ]; then
  printf "%s\n" "output directory already exists: ${ODIR}"
  exit 1
fi

declare -A samplesMap

# samplesMap["RelValTTbar_14TeV"]="/RelValTTbar_14TeV/CMSSW_11_2_0_pre9-PU_112X_mcRun3_2021_realistic_v11-v1/GEN-SIM-DIGI-RAW"
samplesMap["TTbar_14TeV"]="/TT_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW"

recoKeys=(
    # HLT_GRun
    HLT_Run3TRK
    # HLT_Run3TRKPixelOnly
    # HLT_Run3TRKPixelOnlyCleaned
    # HLT_Run3TRKPixelOnlyCleaned2
    # HLT_Run3TRKPixelOnlyCleaned3
    # HLT_Run3TRKPixelOnlyCleaned4
    # HLT_Run3TRKForBTag
    # HLT_Run3TRKForBTag_2
    # HLT_Run3TRKForBTag_3
)

# options (JobFlavour and AccountingGroup)
opts=""
if [[ ${HOSTNAME} == lxplus* ]]; then
  opts+="--JobFlavour longlunch"
  # opts+="--JobFlavour workday"
  if [[ ${USER} == sewuchte ]]; then
    opts+=" --AccountingGroup group_u_CMS.CAF.PHYS"
  fi
fi

for recoKey in "${recoKeys[@]}"; do
    python3 ../runHLTBTagAnalyzer_cfg.py dumpPython=.tmp_${recoKey}_cfg.py defaults=Run3 runOnData=False reco=${recoKey} runPuppiJetVariables=False runCaloJetVariables=False
    for sampleKey in ${!samplesMap[@]}; do
        sampleName=${samplesMap[${sampleKey}]}
        numEvents=${NEVT}
        # htc_driver -c /tmp/${USER}/BTV_Run3_Online_cfg.py --customize-cfg -m ${numEvents} -n 1000 --cpus 1 --memory 2000 --runtime 10800 ${opts} \
        # htc_driver -c /tmp/${USER}/BTV_Run3_Online_cfg.py --customize-cfg -m ${numEvents} -n 1000 ${opts} \
        bdriver -c .tmp_${recoKey}_cfg.py --customize-cfg -m ${numEvents} -n 1000 ${opts} --cpus 1 --mem 1999 --time 03:00:00 \
        -d ${sampleName} -p 0 -o ${ODIR}/${recoKey}/${sampleKey} \
        --customise-commands \
       '# output [TFileService]' \
       "if hasattr(process, 'TFileService'):" \
       '  process.TFileService.fileName = opts.output'
        # --cmsRun-output-dir ${ODIR_cmsRun}/${recoKey}/${sampleKey}
    done
    unset sampleKey sampleName numEvents
    rm -f .tmp_${recoKey}_cfg.py
done
# unset opts samplesMap NEVT ODIR ODIR_cmsRun recoKey
unset opts samplesMap NEVT ODIR recoKey

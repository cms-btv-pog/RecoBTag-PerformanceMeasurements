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

# MinBias
# samplesMap["Phase2HLTTDR_MinBias_14TeV_PU140"]="/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT"
# samplesMap["Phase2HLTTDR_MinBias_14TeV_PU200"]="/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT"
# samplesMap["Phase2HLTTDR_MinBias_orig_14TeV_PU200"]="/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT"
#
# # QCD Pt-binned
# samplesMap["Phase2HLTTDR_QCD_Pt020to030_14TeV_PU140"]="/QCD_Pt_20to30_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt030to050_14TeV_PU140"]="/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt030to050_14TeV_PU140_ext1"]="/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt050to080_14TeV_PU140"]="/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt050to080_14TeV_PU140_ext1"]="/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt080to120_14TeV_PU140"]="/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt120to170_14TeV_PU140"]="/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt170to300_14TeV_PU140"]="/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt300to470_14TeV_PU140"]="/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt470to600_14TeV_PU140"]="/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt600toInf_14TeV_PU140"]="/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
#
# samplesMap["Phase2HLTTDR_QCD_Pt020to030_14TeV_PU200"]="/QCD_Pt_20to30_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt030to050_14TeV_PU200"]="/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt030to050_14TeV_PU200_ext1"]="/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt050to080_14TeV_PU200"]="/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt050to080_14TeV_PU200_ext1"]="/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v3/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt080to120_14TeV_PU200"]="/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt120to170_14TeV_PU200"]="/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt170to300_14TeV_PU200"]="/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt300to470_14TeV_PU200"]="/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt470to600_14TeV_PU200"]="/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_QCD_Pt600toInf_14TeV_PU200"]="/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
#
# # V+Jets
# samplesMap["Phase2HLTTDR_WJetsToLNu_14TeV_PU140"]="/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_WJetsToLNu_14TeV_PU200"]="/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
#
# samplesMap["Phase2HLTTDR_DYJetsToLL_M010to050_14TeV_PU140"]="/DYToLL_M-10To50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_pilot_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
# samplesMap["Phase2HLTTDR_DYJetsToLL_M010to050_14TeV_PU200"]="/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"
#
# samplesMap["Phase2HLTTDR_DYJetsToLL_M050toInf_14TeV_PU140"]="/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT"
# samplesMap["Phase2HLTTDR_DYJetsToLL_M050toInf_14TeV_PU200"]="/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT"

# for btv signal
# samplesMap["Phase2HLTTDR_HH4b_14TeV_PU200"]="/GluGluToHHTo4B_node_SM_TuneCP5_14TeV-madgraph_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"
samplesMap["Phase2HLTTDR_TTbar_14TeV_NoPU"]="/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/FEVT"
samplesMap["Phase2HLTTDR_TTbar_14TeV_PU200"]="/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/FEVT"
# samplesMap["Phase2HLTTDR_TTbar2L_14TeV_PU200"]="/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"
# samplesMap["Phase2HLTTDR_TTbar1L_14TeV_PU200"]="/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"
#
# samplesMap["Phase2HLTTDR_TTbar2L_14TeV_PU140"]="/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"
# samplesMap["Phase2HLTTDR_TTbar1L_14TeV_PU140"]="/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"
samplesMap["Phase2HLTTDR_TTbar_14TeV_PU140"]="/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"
# samplesMap["Phase2HLTTDR_HH4b_14TeV_PU140"]="/GluGluToHHTo4B_node_SM_TuneCP5_14TeV-madgraph_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"

# samplesMap["Phase2HLTTDR_ZZTo4bQ01j_14TeV_PU200"]="/ZZTo4bQ01j_5f_TuneCP5_amcatNLO_FXFX_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"

# recoKeys=(
    # HLT_TRKv00_TICL
    # HLT_TRKv06p1_TICL
    # HLT_TRKv06p3_TICL
    # HLT_TRKv07p2_TICL
# )
# /eos/home-s/sewuchte/BTV-Phase2/December_TDR/Condor_Prod_v01
# btvRecoKeys=(
    # default
    # cutsV2
# )

# options (JobFlavour and AccountingGroup)
opts=""
if [[ ${HOSTNAME} == lxplus* ]]; then
  opts+="--JobFlavour longlunch"
  if [[ ${USER} == sewuchte ]]; then
    opts+=" --AccountingGroup group_u_CMS.CAF.PHYS"
  fi
fi

# for recoKey in "${recoKeys[@]}"; do
    # for btvRecoKey in "${btvRecoKeys[@]}"; do
        for sampleKey in ${!samplesMap[@]}; do
            sampleName=${samplesMap[${sampleKey}]}
            numEvents=${NEVT}
            # if [[ ${sampleKey} == *MinBias* ]]; then
                # numEvents=2000000
            # fi
            # python ../runHLTBTagAnalyzer_PhaseII_cfg.py dumpPython=/tmp/${USER}/BTV_PhaseII_Offline_cfg.py numThreads=1 reco=${recoKey} BTVreco=${btvRecoKey} globalTag=111X_mcRun4_realistic_T15_v2
            python ../runBTagAnalyzer_cfg.py dumpPython=/tmp/${USER}/BTV_PhaseII_Offline_cfg.py defaults=PhaseII_puppi runOnData=False

            # htc_driver -c /tmp/${USER}/BTV_PhaseII_Offline_cfg.py --customize-cfg -m ${numEvents} -n 250 --cpus 1 --memory 2000 --runtime 10800 ${opts} \
            htc_driver -c /tmp/${USER}/BTV_PhaseII_Offline_cfg.py --customize-cfg -m ${numEvents} -n 1000 --cpus 1 --memory 2000 --runtime 10800 ${opts} \
            -d ${sampleName} -p 0 -o ${ODIR}/${sampleKey} --cmsRun-output-dir ${ODIR_cmsRun}/${sampleKey}
        done
    # done
    unset sampleKey sampleName numEvents
# done
# unset recoKey opts recoKeys samplesMap btvRecoKeys btvRecoKey NEVT ODIR ODIR_cmsRun
unset opts samplesMap NEVT ODIR ODIR_cmsRun

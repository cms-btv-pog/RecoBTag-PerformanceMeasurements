#!/bin/sh

if   [ $(basename ${SHELL}) = "bash" ]; then BASE_DIR=$(dirname $(readlink -f ${BASH_SOURCE}))
elif [ $(basename ${SHELL}) = "zsh"  ]; then BASE_DIR=$(dirname $(readlink -f $0))
else

  printf "\n%s\n\n" " >>> WARNING -- unsupported type of shell (${SHELL}), no action taken)"
fi

if [ ! "$(command -v scram)" ]; then

  printf "\n%s\n\n" " >>> WARNING -- executable \"scram\" not available, loading CMSSW utilities first"

  if [ ! "$(command -v module)" ]; then

    printf "\n%s\n\n" " >>> ERROR -- executable \"module\" (used to load CMSSW utilities at DESY) not available, installation stopped (enable \"scram\" before installing or update the install script)"
    return

  else

    if [ ! -d /afs/desy.de/group/cms/modulefiles/ ]; then

      printf "\n%s\n\n" " >>> ERROR -- directory /afs/desy.de/group/cms/modulefiles/ (used to load CMSSW utilities at DESY via module) does not exist (please enable \"scram\")"
      return
    fi

    module use -a /afs/desy.de/group/cms/modulefiles/
    module load cmssw
  fi
fi

eval `scram runtime -sh`

source /cvmfs/cms.cern.ch/crab3/crab.sh

voms-proxy-init -voms cms

echo ${BASE_DIR}

export PATH=${PATH}:"${BASE_DIR}"/utilities
export PYTHONPATH=${PYTHONPATH}:"${BASE_DIR}"

unset -v BASE_DIR

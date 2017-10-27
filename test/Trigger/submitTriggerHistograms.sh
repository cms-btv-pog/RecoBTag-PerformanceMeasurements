#!/bin/bash

export WORKDIRECTORY=$PWD
export DATASET=$1
export YEAR=$2
export PERIODS=$3

mkdir -p jobs

cd jobs/

bsub -q cmscaf1nd -o $WORKDIRECTORY/jobs $WORKDIRECTORY/runTriggerHistograms.sh $WORKDIRECTORY $DATASET $YEAR $PERIODS 


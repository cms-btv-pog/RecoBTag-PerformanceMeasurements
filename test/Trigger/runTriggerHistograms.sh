#!/bin/bash

echo "Processing " $2 " for periods " $3 

cd $1 

eval `scramv1 runtime -sh`

./mkTriggerHistograms.py --Dataset=$2 --Year=$3 --Period=$4 --Verbose=3


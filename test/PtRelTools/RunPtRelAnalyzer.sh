#!/bin/sh
echo ""
cd $WORKINGDIRECTORY/
echo "WORKINGDIRECTORY" $PWD
echo ""
eval `scramv1 runtime -sh`
source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh
echo "LSB_JOBINDEX" $LSB_JOBINDEX
#if [$?LSB_JOBINDEX]; then
if [ -z ${LSB_JOBINDEX+x} ]; then
    echo "LSB_JOBINDEX not set"
else
    export DATARANGEINDEX=$LSB_JOBINDEX
fi
echo "TEMPLATEVARIABLE " $TEMPLATEVARIABLE
echo "PUWEIGHTING      " $PUWEIGHTING
echo "KINWEIGHTING     " $KINWEIGHTING
echo "SELECTION        " $SELECTION
echo "MACRONAME        " $MACRONAME
echo "DATARANGEINDEX   " $DATARANGEINDEX
echo ""
root -l -b -q 'RunPtRelAnalyzer.C'
rm core.* 
   

#!/bin/csh

if ($#argv == 1) root -l -b -q CompilePtRelAnalyzer.C

setenv WORKINGDIRECTORY $PWD

setenv TEMPLATEVARIABLE 'PtRel'
#setenv TEMPLATEVARIABLE 'System8' 

setenv PUWEIGHTING '_PSRun2016ICHEP'

#setenv KINWEIGHTING ''
setenv KINWEIGHTING '_KinPtBinsCentral'

#setenv SELECTION ''
#setenv SELECTION '_TrgEmul'
setenv SELECTION '_TrgConf'

#setenv MACRONAME 'ComputePileUpWeightsPSQCDMuQCDXJetHT'
#setenv MACRONAME 'ComputePileUpWeights_PSVRun2016BJuly15:JetHTQCDMuQCDX'
#setenv MACRONAME 'ComputePileUpWeightsTestPV'
setenv MACRONAME 'BuildTemplatesAll'
#setenv MACRONAME 'ComputeKinematicWeightsQCDMuJetHTQCDXAll'
#setenv MACRONAME 'CompareDataToMC_anyEta'
#setenv MACRONAME 'ComputePtRelScaleFactorsCSVv2:_Central:'
#setenv MACRONAME 'PlotBTagPerformanceICHEP2016'
#setenv MACRONAME 'AnalyzeSystematics'
#setenv MACRONAME 'StoreScaleFactors'

setenv DATARANGEINDEX '0'
echo ""

cd $WORKINGDIRECTORY/
echo "WORKINGDIRECTORY" $PWD
echo ""

eval `scramv1 runtime -csh`
source /afs/cern.ch/project/eos/installation/cms/etc/setup.csh

echo "TEMPLATEVARIABLE " $TEMPLATEVARIABLE
echo "PUWEIGHTING      " $PUWEIGHTING
echo "KINWEIGHTING     " $KINWEIGHTING
echo "SELECTION        " $SELECTION
echo "MACRONAME        " $MACRONAME
echo "DATARANGEINDEX   " $DATARANGEINDEX
echo ""

mkdir -p Weights/KinematicWeights
mkdir -p Weights/PileUp
mkdir -p Tables
mkdir -p Plots

root -l -b -q 'RunPtRelAnalyzer.C'
 
   

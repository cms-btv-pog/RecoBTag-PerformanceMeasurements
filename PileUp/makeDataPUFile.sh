#!/bin/sh
JSONALL2017=Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_ALL2017.txt           
JSONRunBCDEF=Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunBCDEF.txt

JSONRUNB=Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RUNB.txt      
JSONRunCDE304670=Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunCDE304670.txt
JSONRunE304671PlusF=Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunE304671PlusF.txt

PUJSON=pileup_latest.txt

MINBIAS=69200

#for JSON in $JSONALL2017 $JSONRUNB $JSONRunCDE304670 $JSONRunE304671PlusF $JSONRunBCDEF;
for JSON in JSONRUNB;
    do
        echo "JSON",$JSON
        pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON  --calcMode true --minBiasXsec $MINBIAS --maxPileupBin 100 --numPileupBins 100  pileup_${JSON}.root
    done

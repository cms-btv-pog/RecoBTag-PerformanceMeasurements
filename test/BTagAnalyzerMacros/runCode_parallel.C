#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include "stdio.h"
using namespace std;

// for  MuEnriched QCD pythia at 13TeV, provide : 
//  -->  m.Fill_nevent(    0.,n20_30,n30_50,n50_80,n80_120,n120_170,n170-300,n300-470,n470-600,n600-800,n800-1000,n1000-inf);
//
// for Inclusive QCD pythia at 13TeV, provide : 
//  -->  m.Fill_nevent(    0.,n15_30,n30_50,n50_80,n80_120,n120_170,n170-300,n300-470,n470-600,n600-800,n800-1000,n1000-inf);

void runCode_parallel(TString charForChain, Int_t nFile, TString outFileName, bool isMuEn, Int_t TrigerPt, Int_t jetPtMin, Int_t jetPtMax, TString PUreweightingFile = "", bool officialPU = false){
    //TO DO if you change something in CommPlot : dans root gROOT->ProcessLine(".L /home/fynu/bfrancois/bTag/CMSSW_7_4_5/src/RecoBTag/PerformanceMeasurements/test/BTagAnalyzerMacros/CommPlotProducer.C++");
    cout << "LOADING CommPlotProducer_C.so..." << endl;
    gSystem->Load("/home/fynu/bfrancois/bTag/CMSSW_7_4_5/src/RecoBTag/PerformanceMeasurements/test/BTagAnalyzerMacros/CommPlotProducer_C.so");
    TChain *superTree = new TChain("btagana/ttree");
    if (nFile == 0)
    {
        superTree->Add(charForChain);
    }
    else
    {
        for(Int_t i = 0; i < nFile; i++)
        {
            TString index, fileName = charForChain;
            index.Form("%d",i+1);
            cout << index << endl;
            fileName.ReplaceAll("INDEX",index);
            cout << fileName << endl;
            superTree->Add(fileName);
        }
    }
    //Int_t nEvent = superTree->GetEntries();
    CommPlotProducer m(superTree);
    TString trigName = "jet";
    if(!isMuEn)
    {
        m.SetInfo("pythia",0,13);
        m.Fill_nevent(0., 9800608.,9930948.,9897538.,6986123.,6778942.,6914063.,5874780. , 3928855., 3869318., 3843252., 2999055.);   // Inclusive QCD runIISpring16_MiniAODv2 v3 HIP !!! Correct up to 470 GeV !
        //m.Fill_nevent(0., 9800608., 9930948., 9968391., 6986123., 6339488., 6914063., 5874780., 3928855., 3959746., 3883812., 2999055.);   // Inclusive QCD, runIISpring16_MiniAODv2 v3
        cout << "Inclusive QCD" << endl;
    }
    else
    {   
        m.SetInfo("pythia",1,13);
        m.Fill_nevent(0., 579682., 432168., 313919., 371584., 410021., 414704., 0., 0., 0., 0., 0.); // for newJP star 1
        cout << "MuEnriched QCD" << endl;
        trigName = "btag";
    }
    m.SetXS();     // Assign the correct x-sections to QCD pthat bins, depending on SetInfo(), default = use inclusive pythia x-sections for 8 TeV.
    m.SetSumXS();
    TString name_root = outFileName; 
    m.Loop(trigName, TrigerPt, jetPtMin, jetPtMax, name_root, PUreweightingFile, officialPU);    
}



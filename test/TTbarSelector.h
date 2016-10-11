#ifndef TTbarSelector_h
#define TTbarSelector_h

#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <vector>
#include <map>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

class TTbarSelector
{
        public:

                // Constructor
                TTbarSelector(){};

                // Copy constructor 
                TTbarSelector(const TTbarSelector &){};

                // Destructor 
                ~TTbarSelector(){};
                
                // Apply trigger bits efficiency
                bool passTrigger(bool isData, Int_t ttbar_chan, Int_t ttbar_trigWord);
                // semilepton trigger
                bool passSingleTrigger(bool isData, Int_t ttbar_chan, Int_t ttbar_trigWord);
                
                // Return true if the event pass the selection
                bool passTTbarSelection(bool isData, vector<TLorentzVector> theLeptColl, vector<Int_t> theLeptIds, vector< pair< TLorentzVector, Float_t> > theJetColl, Int_t ttbar_trigWord, Float_t ttbar_w[250], int ttbar_nw, TH1F* wgtcounter, TString syst, bool computeEvtWgtOnly);
                //Reture true if the event pass the semilep selection
		bool passSemiLepTTbarSelection(bool isData, vector<TLorentzVector> theLeptColl, vector<Int_t> theLeptIds, vector< pair< TLorentzVector, Float_t> > theJetColl, Int_t ttbar_trigWord, Float_t ttbar_w[250], int ttbar_nw, TH1F* wgtcounter, TString syst, bool computeEvtWgtOnly); 

                // Apply object selection on event (fill p4, variables... )	
                pair<float,float> getTriggerEfficiency(int channel);
                pair<float,float> getLeptonSelectionEfficiencyScaleFactor(int id,float pt,float eta);
                vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
	
                TLorentzVector  lept1_, lept2_;
                Float_t mll_ = 0;	
                Float_t met_ = 0;	
                vector< pair< TLorentzVector, Int_t> > theSelJetColl;
                Float_t evWgt;

        	//jet uncertainty parameterization
	        string thePathToJECfile_ = string(getenv("CMSSW_BASE")) + "/src/RecoBTag/PerformanceMeasurements/test/ttbar/data/Summer15_25nsV5_DATA_Uncertainty_AK4PFchs.txt";
	        JetCorrectionUncertainty *theJECuncertainty_ = new JetCorrectionUncertainty(thePathToJECfile_);
	
		
	private :
	
		TString thechannel_;
		
		
};

#endif // #ifdef TTbarSelector_cxx

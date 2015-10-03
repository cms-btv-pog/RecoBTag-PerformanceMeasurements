#ifndef _TTbarEventAnalysis_h_
#define _TTbarEventAnalysis_h_

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TSystem.h"
#include "Math/VectorUtil.h"

#include "TMVA/Reader.h"

#include <iostream>
#include <map>
#include <vector>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

struct LJKinematics_t
{
  Float_t dr,dphi,deta,ptrel,mlj,kindisc;
  Float_t lj2ll_deta,lj2ll_dphi;
};

bool sortLJKinematicsByDR (LJKinematics_t i,LJKinematics_t j) { return (i.dr<j.dr); }

class TTbarEventAnalysis
{
 public:
 TTbarEventAnalysis() : 
  tmvaReader_(0),
    useOnlySignOfGenWeight_(false),
    readTTJetsGenWeights_(false),
    applyMETFilters_(false),
    applyTriggerEff_(true),
    applyLepSelEff_(true) 
      {
	//jet uncertainty parameterization
	TString jecUncUrl("${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/ttbar/Summer15_50nsV5_DATA_Uncertainty_AK4PFchs.txt");
	gSystem->ExpandPathName(jecUncUrl);
	jecUnc_ = new JetCorrectionUncertainty(jecUncUrl.Data());
      }
  ~TTbarEventAnalysis(){}
  void setUseOnlySignOfGenWeight(bool useOnlySignOfGenWeight) { useOnlySignOfGenWeight_=useOnlySignOfGenWeight; }
  void setReadTTJetsGenWeights(bool readTTJetsGenWeights)     { readTTJetsGenWeights_=readTTJetsGenWeights; }
  void setApplyMETFilters(bool applyMETFilters)               { applyMETFilters_=applyMETFilters; }
  void setApplyTriggerEff(bool applyTriggerEff)               { applyTriggerEff_=applyTriggerEff; }
  void setApplyLepSelEff(bool applyLepSelEff)                 { applyLepSelEff_=applyLepSelEff; }
  void setTMVAWeightsBaseDir(TString url)                     { weightsDir_=url; gSystem->ExpandPathName(weightsDir_); }
  void addTriggerBit(Int_t bit,Int_t ch)                      { triggerBits_.push_back(std::pair<Int_t,Int_t>(bit,ch)); }
  void addVarForTMVA(TString varName)                         { tmvaVarNames_.push_back(varName); }
  void prepareOutput(TString outFile);
  void processFile(TString inFile,float xsecWgt);
  void finalizeOutput();

 private:
  JetCorrectionUncertainty *jecUnc_;
  std::pair<float,float> getTriggerEfficiency(int id1,float pt1,float eta1,int id2,float pt2,float eta2,int ch);
  std::pair<float,float> getLeptonSelectionEfficiencyScaleFactor(int id,float pt,float eta);
  std::vector<float> getJetEnergyScales(float pt,float eta,float rawsf,float area,float rho);
  std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);

  bool useOnlySignOfGenWeight_, readTTJetsGenWeights_, applyMETFilters_, applyTriggerEff_, applyLepSelEff_;
  TString weightsDir_;
  std::vector<TString> tmvaVarNames_;
  TMVA::Reader *tmvaReader_;
  TFile *outF_;
  Int_t eventInfo_[3],ttbar_chan_;
  Float_t weight_[15];
  Int_t jetFlavour_[2],jetmult_,jetrank_;
  Float_t jetPt_[2],jetEta_[2];
  Float_t close_mlj_[5],close_deta_,close_dphi_,close_ptrel_,close_lj2ll_deta_, close_lj2ll_dphi_;
  Float_t far_mlj_, far_deta_, far_dphi_, far_ptrel_,far_lj2ll_deta_, far_lj2ll_dphi_;
  Float_t  j2ll_deta_,j2ll_dphi_;
  Float_t kinDisc_[5];
  Float_t jp_[2],svhe_[2],csv_[2];
  std::vector<std::pair<Int_t,Int_t> > triggerBits_;
  TTree *kinTree_,*ftmTree_;
  std::map<TString,TH1F *> histos_;
  
};

#endif

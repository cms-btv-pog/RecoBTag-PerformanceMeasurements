#ifndef _TTbarEventAnalysis_h_
#define _TTbarEventAnalysis_h_

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TSystem.h"
#include "Math/VectorUtil.h"
#include "TGraph.h"
#include "TMVA/Reader.h"

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

struct LJKinematics_t
{
  Float_t dr,dphi,deta,ptrel,mlj,kindisc;
  Float_t lj2ll_deta,lj2ll_dphi;
};

bool sortLJKinematicsByDR (LJKinematics_t i,LJKinematics_t j) { return (i.dr<j.dr); }

//prepare to read the tree (for jets only interested in a couple of variables)
class MyEventInfoBranches_t 
{
  public:
    Int_t Run,Evt,LumiBlock,nPV;
    Int_t   ttbar_chan, ttbar_trigWord, ttbar_metfilterWord;
    Int_t   ttbar_nl, ttbar_lid[10], ttbar_lgid[10], ttbar_lch[10];
    Float_t ttbar_lpt[10], ttbar_leta[10], ttbar_lphi[10], ttbar_lm[10];
    Float_t ttbar_metpt,ttbar_metphi;
    Float_t ttbar_rho;
    Int_t   ttbar_nw;
    Float_t nPUtrue;
    Int_t nPU;
    Float_t ttbar_w[1095];
    Int_t nJet;
    Float_t Jet_pt[100],Jet_genpt[100],Jet_area[100],Jet_jes[100],Jet_eta[100],Jet_phi[100],Jet_mass[100];
    Float_t Jet_Svx[100],Jet_CombIVF[100],Jet_Proba[100],Jet_Ip2P[100];
    Float_t Jet_DeepCSVb[100], Jet_DeepCSVc[100], Jet_DeepCSVl[100], Jet_DeepCSVbN[100], Jet_DeepCSVcN[100], Jet_DeepCSVlN[100];
    Float_t Jet_DeepCSVBDisc[100],Jet_DeepCSVBDiscN[100],Jet_DeepCSVCvsLDisc[100],Jet_DeepCSVCvsLDiscN[100],Jet_DeepCSVCvsBDisc[100],Jet_DeepCSVCvsBDiscN[100];
    Float_t Jet_DeepFlavourBDisc[100], Jet_DeepFlavourCvsLDisc[100], Jet_DeepFlavourCvsBDisc[100];
    Float_t Jet_DeepFlavourB[100];
    Int_t Jet_nseltracks[100];
    Int_t Jet_flavour[100];
};
    
class TTbarEventAnalysis
{
 public:
 TTbarEventAnalysis() : 
    tmvaReader_(0),
    readTTJetsGenWeights_(true),   
    puWgtGr_(0),puWgtDownGr_(0),puWgtUpGr_(0)
      {
  //jet uncertainty parameterization
  TString jecUncUrl("${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/ttbar/data/Autumn18_V8_MC_Uncertainty_AK4PF.txt");
  gSystem->ExpandPathName(jecUncUrl);
  jecUnc_ = new JetCorrectionUncertainty(jecUncUrl.Data());

  //pileup weights
  TString stdTarget("${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/ttbar/data/pileupWgts.root");
  gSystem->ExpandPathName(stdTarget);
  SetPUWeightTarget(stdTarget,"");
      }
  ~TTbarEventAnalysis(){}
  void setReadTTJetsGenWeights(bool readTTJetsGenWeights)     { readTTJetsGenWeights_=readTTJetsGenWeights; }
  void setTwoTagCount(bool runTwoTagCount_)                   { runTwoTagAnalysis=runTwoTagCount_; }
  void setTMVAWeightsBaseDir(TString url)                     { weightsDir_=url; gSystem->ExpandPathName(weightsDir_); }
  void addTriggerBit(Int_t bit,Int_t ch)                      { triggerBits_.push_back(std::pair<Int_t,Int_t>(bit,ch)); }
  void addVarForTMVA(TString varName)                         { tmvaVarNames_.push_back(varName); }
  void prepareOutput(TString outFile);
  void processFile(TString inFile,TH1F *xsecWgt,Bool_t isData);
  void finalizeOutput();
  void SetPUWeightTarget(TString targetFile,TString sampleName){
      TFile *fIn=TFile::Open(targetFile);
      if(fIn){
        std::string nom("puwgts_nom");  
        std::string up("puwgts_down");
        std::string down("puwgts_up");
        if(!sampleName.IsNull()){
            nom.append("_");
            nom.append(sampleName);
            up.append("_");
            up.append(sampleName);
            down.append("_");
            down.append(sampleName);
        }else{
            std::cout << "Warning: PUWeight target histogram " << sampleName << " does not exist. Check naming convention of samples matches that of PU histograms." << std::endl;
        }

        puWgtGr_     = (TGraph *)fIn->Get(nom.c_str());
        puWgtDownGr_ = (TGraph *)fIn->Get(down.c_str());
        puWgtUpGr_   = (TGraph *)fIn->Get(up.c_str());
      }else{
        std::cout << "Unable to find data/pileupWgts.root, no PU reweighting will be applied" << std::endl;
      }
      fIn->Close();
  }

 private:
  JetCorrectionUncertainty *jecUnc_;
  std::pair<float,float> getTriggerEfficiency(int id1,float pt1,float eta1,int id2,float pt2,float eta2,int ch);
  std::pair<float,float> getLeptonSelectionEfficiencyScaleFactor(int id,float pt,float eta);
  std::vector<float> getJetEnergyScales(float pt,float eta,float rawsf,float area,float rho);
  std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);

  MyEventInfoBranches_t ev;
  std::vector<Int_t> selJets;
  Float_t evWgt;
  
  std::vector<std::string> wpLabel;
  std::vector<std::string> systName;
  std::map<std::string,double> systWeight;

  void TwoTag(std::string tagName, std::string discriminator, std::pair<int, int>, int ptBin);  
  std::map<std::string,std::vector<float>> btaggingWPs;
  void GetBestJetPair(std::pair<int, int>& myIndices, std::string discriminator="deepCSV");
  std::map<std::string,std::pair<int, int>> bestJetPairs;
  float ReturnVarAtIndex(std::string varName, unsigned int index);
  std::string ReturnPtLabel(int iPT);
  std::vector<unsigned int> lowerPtBinEdges;
  std::vector<std::string> twoTagNames; 
 
  TGraph *puWgtGr_,*puWgtDownGr_,*puWgtUpGr_;
  bool readTTJetsGenWeights_;
  TString weightsDir_;
  std::vector<TString> tmvaVarNames_;
  TMVA::Reader *tmvaReader_;
  TFile *outF_;
  Int_t eventInfo_[3],ttbar_chan_,npv_;
  Float_t weight_[27];
  Int_t jetFlavour_[2],jetmult_,jetrank_;
  Float_t jetPt_[2],jetEta_[2];
  Float_t close_mlj_[5],close_deta_,close_dphi_,close_ptrel_,close_lj2ll_deta_, close_lj2ll_dphi_;
  Float_t far_mlj_, far_deta_, far_dphi_, far_ptrel_,far_lj2ll_deta_, far_lj2ll_dphi_;
  Float_t  j2ll_deta_,j2ll_dphi_;
  Float_t kinDisc_[5];
  Float_t jp_[2],svhe_[2],csv_[2];
  Float_t DeepCSVb_[2], DeepCSVc_[2], DeepCSVl_[2], DeepCSVbb_[2], DeepCSVcc_[2], DeepCSVbN_[2], DeepCSVcN_[2], DeepCSVlN_[2], DeepCSVbbN_[2], DeepCSVccN_[2], DeepCSVbP_[2], DeepCSVcP_[2], DeepCSVlP_[2], DeepCSVbbP_[2], DeepCSVccP_[2];
  Float_t DeepCSVBDisc_[2], DeepCSVBDiscN_[2], DeepCSVBDiscP_[2], DeepCSVCvsLDisc_[2], DeepCSVCvsLDiscN_[2], DeepCSVCvsLDiscP_[2], DeepCSVCvsBDisc_[2], DeepCSVCvsBDiscN_[2], DeepCSVCvsBDiscP_[2];
  Float_t DeepFlavourBDisc_[2], DeepFlavourCvsLDisc_[2], DeepFlavourCvsBDisc_[2];
  Float_t DeepFlavourB_[2], DeepFlavourBB_[2], DeepFlavourLEPB_[2];
  std::vector<std::pair<Int_t,Int_t> > triggerBits_;
  TTree *kinTree_,*ftmTree_;
  std::map<TString,TH1F *> histos_;
  std::map<TString,TH2F *> histos2d_;
  bool noEventsSelected;
  unsigned int runTwoTagAnalysis;
  
};

#endif

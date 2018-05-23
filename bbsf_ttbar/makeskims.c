#include <TH1D.h>
#include <TVector3.h>
#include <iostream>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

class makeskims{

   public:

   void run(const TString intag, const TString infiles, const float Weight, const bool ismc);
   
   private:

   bool isGoodMuon(const int idx);
   bool isGoodAK4(const int idx);
   bool isGoodAK8(const int idx);
   bool passTrigger();
   int nSubjetOverlap(const int id4, const int id8);
   TChain *chain4, *chain8;
   void chainInit(const bool ismc);

   //trigger
   int ttbar_trigWord;
   TBranch * b_ttbar_trigWord;
   //muons
   int nPatMuon;
   float PatMuon_pt[100], PatMuon_eta[100], PatMuon_phi[100];
   int PatMuon_isTightMuon[100];
   float PatMuon_iso[100];
   TBranch * b_nPatMuon;
   TBranch * b_PatMuon_pt;
   TBranch * b_PatMuon_eta;
   TBranch * b_PatMuon_phi;
   TBranch * b_PatMuon_isTightMuon;
   TBranch * b_PatMuon_iso;
   //AK4 jets
   int nJet;
   float Jet_pt[100], Jet_eta[100], Jet_phi[100];
   float Jet_DeepCSVb[100];
   int Jet_tightID[100];
   TBranch * b_nJet;
   TBranch * b_Jet_pt;
   TBranch * b_Jet_eta;
   TBranch * b_Jet_phi;
   TBranch * b_Jet_DeepCSVb;
   TBranch * b_Jet_tightID;
   //AK8 jets
   int nJet8;
   float Jet8_pt[100], Jet8_eta[100], Jet8_phi[100];
   float Jet8_massSoftDrop[100], Jet8_tau2[100], Jet8_tau1[100], Jet8_DoubleSV[100];
   float Jet8_massPruned[100];
   int Jet8_tightID[100];
   TBranch * b_nJet8;
   TBranch * b_Jet8_pt;
   TBranch * b_Jet8_eta;
   TBranch * b_Jet8_phi;
   TBranch * b_Jet8_DoubleSV;
   TBranch * b_Jet8_massSoftDrop;
   TBranch * b_Jet8_tau2;
   TBranch * b_Jet8_tau1;
   TBranch * b_Jet8_tightID;
   //subjets
   int sub_nJet;
   int sub_Jet_FatJetIdx[1000];
   float sub_Jet_pt[1000];
   float sub_Jet_eta[1000];
   float sub_Jet_phi[1000];
   TBranch * b_sub_nJet;
   TBranch * b_sub_Jet_FatJetIdx;
   TBranch * b_sub_Jet_pt;
   TBranch * b_sub_Jet_eta;
   TBranch * b_sub_Jet_phi;
   //truth
   float nPUtrue;
   TBranch * b_nPUtrue;
};

void makeskims::chainInit(const bool ismc=false)
{
   //triggers
   chain4->SetBranchAddress("ttbar_trigWord", &ttbar_trigWord, &b_ttbar_trigWord);
   //muons
   chain4->SetBranchAddress("nPatMuon", &nPatMuon, &b_nPatMuon);
   chain4->SetBranchAddress("PatMuon_pt", PatMuon_pt, &b_PatMuon_pt);
   chain4->SetBranchAddress("PatMuon_eta", PatMuon_eta, &b_PatMuon_eta);
   chain4->SetBranchAddress("PatMuon_phi", PatMuon_phi, &b_PatMuon_phi);
   chain4->SetBranchAddress("PatMuon_isTightMuon", PatMuon_isTightMuon, &b_PatMuon_isTightMuon);
   chain4->SetBranchAddress("PatMuon_iso", PatMuon_iso, &b_PatMuon_iso);
   //AK4 jets
   chain4->SetBranchAddress("nJet", &nJet, &b_nJet);
   chain4->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   chain4->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   chain4->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   chain4->SetBranchAddress("Jet_tightID", Jet_tightID, &b_Jet_tightID);
   chain4->SetBranchAddress("Jet_DeepCSVb", Jet_DeepCSVb, &b_Jet_DeepCSVb);
   //AK8 jets   
   chain8->SetBranchAddress("FatJetInfo.nJet", &nJet8, &b_nJet8);
   chain8->SetBranchAddress("FatJetInfo.Jet_pt", Jet8_pt, &b_Jet8_pt);
   chain8->SetBranchAddress("FatJetInfo.Jet_eta", Jet8_eta, &b_Jet8_eta);
   chain8->SetBranchAddress("FatJetInfo.Jet_phi", Jet8_phi, &b_Jet8_phi);
   chain8->SetBranchAddress("FatJetInfo.Jet_massSoftDrop", Jet8_massSoftDrop, &b_Jet8_massSoftDrop);
   chain8->SetBranchAddress("FatJetInfo.Jet_tau2", Jet8_tau2, &b_Jet8_tau2);
   chain8->SetBranchAddress("FatJetInfo.Jet_tau1", Jet8_tau1, &b_Jet8_tau1);
   chain8->SetBranchAddress("FatJetInfo.Jet_DoubleSV", Jet8_DoubleSV, &b_Jet8_DoubleSV);
   chain8->SetBranchAddress("FatJetInfo.Jet_tightID", Jet8_tightID, &b_Jet8_tightID);
   //subjets
   chain8->SetBranchAddress("SoftDropPuppiSubJetInfo.nJet", &sub_nJet, &b_sub_nJet);
   chain8->SetBranchAddress("SoftDropPuppiSubJetInfo.Jet_FatJetIdx", sub_Jet_FatJetIdx, &b_sub_Jet_FatJetIdx);
   chain8->SetBranchAddress("SoftDropPuppiSubJetInfo.Jet_pt", sub_Jet_pt, &b_sub_Jet_pt);
   chain8->SetBranchAddress("SoftDropPuppiSubJetInfo.Jet_eta", sub_Jet_eta, &b_sub_Jet_eta); 
   chain8->SetBranchAddress("SoftDropPuppiSubJetInfo.Jet_phi", sub_Jet_phi, &b_sub_Jet_phi);
   //puweights
   if (ismc) chain4->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);     
}

int makeskims::nSubjetOverlap(const int id4, const int id8)
{
   TVector3 ak4jet;
   ak4jet.SetPtEtaPhi(Jet_pt[id4], Jet_eta[id4], Jet_phi[id4]); 
   int n = 0;
   for (int i = 0; i < sub_nJet; ++i) {
      if (sub_Jet_FatJetIdx[i] == id8) {
         if (isGoodAK4(i)) {
            TVector3 subjet;
            subjet.SetPtEtaPhi(sub_Jet_pt[i], sub_Jet_eta[i], sub_Jet_phi[i]);
            const float dr = subjet.DeltaR(ak4jet);
            if (dr < 0.4) ++n;
         }
      }
   }
   return n;
}

//https://twiki.cern.ch/CMS/SWGuideMuonIdRun2
bool makeskims::isGoodMuon(const int idx)
{
   if (PatMuon_pt[idx] >= 50. && std::abs(PatMuon_eta[idx]) < 2.1) {
      if (PatMuon_iso[idx] < 0.15) {
         if (PatMuon_isTightMuon[idx]==1) {
            return true;
         }
      }
   }
   return false;
}

//https://twiki.cern.ch/CMS/JetID13TeVRun2017
bool makeskims::isGoodAK4(const int idx)
{
   if (Jet_pt[idx] >= 30. && std::abs(Jet_eta[idx]) < 2.4) {
      if (Jet_tightID[idx]==1) {
         return true;
      }
   }
   return false;
}

//https://twiki.cern.ch/CMS/JetID13TeVRun2017
bool makeskims::isGoodAK8(const int idx)
{
   if (Jet8_pt[idx] >= 250. && std::abs(Jet8_eta[idx]) < 2.4) {
      if (Jet8_tau2[idx]/Jet8_tau1[idx] < 0.6) {
         if (Jet8_massSoftDrop[idx] >= 50. && Jet8_massSoftDrop[idx] < 200.) {
            if (Jet8_tightID[idx]==1) {
               return true;
            }
         }
      }
   }
   return false;
}

bool makeskims::passTrigger()
{
   const int pass = (ttbar_trigWord >> 8) & 1;
   return pass;
}

void makeskims::run(const TString intag, TString infiles, float xsWeight=1., const bool ismc=false)
{
   TH1D * sf[3];
   TH1D * sfclone[3];
   if (ismc) {
      const TString errname[3] = {"down", "central", "up"};
      TFile * f_pu = TFile::Open("pu.root");
      for (int i = 0; i < 3; ++i) {
         TH1D * temphist = (TH1D*)f_pu->Get("sf_" + errname[i] + "_" + intag);
         sf[i] = (TH1D*)temphist->Clone("sf0_" + errname[i] + "_" + intag);
         sfclone[i] = (TH1D*)temphist->Clone("sf0clone_" + errname[i] + "_" + intag);
      }
   }
   //f_pu->Close();
 
   chain4 = new TChain("btagana/ttree");
   chain8 = new TChain("btaganaFatJets/ttree");
   int nfiles = 0;
   for (int i = 0; i < 3; ++i) {
      TString tempfile;
      if (i==0) tempfile = infiles + "0000/*.root";
      if (i==1) tempfile = infiles + "0001/*.root";
      if (i==2) tempfile = infiles + "0002/*.root";
      std::cout << "adding " << tempfile << std::endl;
      nfiles += chain4->Add(tempfile);
      chain8->Add(tempfile);
   }
   
   //nfiles += chain4->Add("/uscms_data/d3/fojensen/bbsfProduction_10052018/CMSSW_9_4_4/src/JetTree_mc_FatJets_Subjets.root");
   //          chain8->Add("/uscms_data/d3/fojensen/bbsfProduction_10052018/CMSSW_9_4_4/src/JetTree_mc_FatJets_Subjets.root");
   
   std::cout << nfiles << " files in the dataset..." << std::endl; 
   std::cout << chain4->GetEntries() << " entries in the chain..." << std::endl;
   chainInit(ismc);

   const TString name = "./skims/"+intag+".root";
   std::cout << "making file " << name << std::endl;
   TFile *f = new TFile(name, "recreate");
   chain4->LoadTree(0);
   TTree *tree4 = chain4->GetTree()->CloneTree(0);
   tree4->SetName("treeAK4");
   chain8->LoadTree(0);
   TTree *tree8 = chain8->GetTree()->CloneTree(0);
   tree8->SetName("treeAK8");

   int idx_muon = -1;   TBranch * b_idx_muon = tree4->Branch("idx_muon", &idx_muon, "idx_muon/I");
   int idx_ak8 = -1;    TBranch * b_idx_ak8 = tree8->Branch("idx_ak8", &idx_ak8, "idx_ak8/I");
   int idx_lepak4 = -1; TBranch * b_idx_lepak4 = tree4->Branch("idx_lepak4", &idx_lepak4, "idx_lepak4/I");
   int idx_hadak4 = -1; TBranch * b_idx_hadak4 = tree4->Branch("idx_hadak4", &idx_hadak4, "idx_hadak4/I");
   // add special mc branches
   TBranch * b_xsWeight;
   float puWeight = -1;      TBranch * b_puWeight;
   float puWeight_up = -1;   TBranch * b_puWeight_up;
   float puWeight_down = -1; TBranch * b_puWeight_down; 
   if (ismc) {
      b_xsWeight = tree4->Branch("xsWeight", &xsWeight, "xsWeight/F");
      b_puWeight = tree4->Branch("puWeight", &puWeight, "puWeight/F");
      b_puWeight_up = tree4->Branch("puWeight_up", &puWeight_up, "puWeight_up/F");
      b_puWeight_down = tree4->Branch("puWeight_down", &puWeight_down, "puWeight_down/F");
   }
  
   const float twothirdspi = (2./3.)*TMath::Pi();

   int nmuonevents = 0;
   int nak8events = 0;
   int ntwojetevents = 0;
   int nisoak8events = 0;
   int npasstriggerevents = 0;

   for (int j = 0; j < chain4->GetEntries(); ++j) {
   
      if (j%100000==0) std::cout << "now beginning entry " << j << std::endl;
      chain4->GetEntry(j);
      chain8->GetEntry(j);

      //single muon
      idx_muon = -1;
      int nmuon = 0;
      for (int k = 0; k < nPatMuon; ++k) {
         if (isGoodMuon(k)) {
            idx_muon = k;
            ++nmuon;
         }
      }
      if (!(nmuon==1)) continue;
      ++nmuonevents;
      TVector3 muon;
      muon.SetPtEtaPhi(PatMuon_pt[idx_muon], PatMuon_eta[idx_muon], PatMuon_phi[idx_muon]);

      //AK8 jet (W) far from the muon
      idx_ak8 = -1;
      for (int k = 0; k < nJet8; ++k) {
         if (isGoodAK8(k)) {
            TVector3 tempvec8;
            tempvec8.SetPtEtaPhi(Jet8_pt[k], Jet8_eta[k], Jet8_phi[k]);
            if (std::abs(muon.DeltaPhi(tempvec8)) >= twothirdspi) {
               idx_ak8 = k;
               break;
            }
         }
      }
      if (idx_ak8==-1) continue;
      ++nak8events;

      //b jets from top decays
      idx_lepak4 = -1;
      idx_hadak4 = -1;
      bool haveLepAK4 = false;
      bool haveHadAK4 = false;
      for (int k = 0; k < nJet; ++k) {
         if (isGoodAK4(k)) {
            TVector2 tempvec4;
            tempvec4.SetMagPhi(Jet_pt[k], Jet_phi[k]);
            if (std::abs(muon.XYvector().DeltaPhi(tempvec4)) < twothirdspi) {
               if (!haveLepAK4) {    
                  if (Jet_DeepCSVb[k] >= 0.1522) {
                     idx_lepak4 = k;
                     haveLepAK4 = true;
                  }
               }
            } else {
               if (!haveHadAK4) {
                  if (nSubjetOverlap(k, idx_ak8)==0) {
                     idx_hadak4 = k;
                     haveHadAK4 = true;
                  }
               }
            }
         }
      }
      if (!haveLepAK4) continue;
      if (!haveHadAK4) continue;
      ++ntwojetevents;

      if (!passTrigger()) continue;
      ++npasstriggerevents;

      //require the AK8 jet (W) to be isolated
      TVector3 Wak8;
      Wak8.SetPtEtaPhi(Jet8_pt[idx_ak8], Jet8_eta[idx_ak8], Jet8_phi[idx_ak8]);
      double minAK8AK4dr = 999.;
      for (int k = 0; k < nJet; ++k) {
         if (isGoodAK4(k)) {
            if (nSubjetOverlap(k, idx_ak8)==0) { 
               TVector3 tempvec4;
               tempvec4.SetPtEtaPhi(Jet_pt[k], Jet_eta[k], Jet_phi[k]);
               const float tempdr = tempvec4.DeltaR(Wak8);
               if (tempdr < minAK8AK4dr) minAK8AK4dr = tempdr;
            }
         }
      }
      if (!(minAK8AK4dr>=0.8)) continue;
      ++nisoak8events;

      if (ismc) {
         puWeight_down = sf[0]->GetBinContent(sfclone[0]->Fill(nPUtrue));
         puWeight      = sf[1]->GetBinContent(sfclone[1]->Fill(nPUtrue));
         puWeight_up   = sf[2]->GetBinContent(sfclone[2]->Fill(nPUtrue));
      }

      tree4->Fill();
      tree8->Fill(); 
   }

   f->Write();
   f->Close();

   std::cout << nmuonevents << std::endl;
   std::cout << nak8events << std::endl;
   std::cout << ntwojetevents << std::endl;
   std::cout << npasstriggerevents << std::endl;
   std::cout << nisoak8events << std::endl;
}

//const float wlep = 3.*0.1086; const float whad = 0.6741;

//makeskims d; d.run("TTToSemiLeptonic", "/eos/uscms/store/user/fojensen/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_BTagAnalyzer_TTToSemiLeptonic/180510_210825/", (831.76*2.*.3258*0.6741)/111381888., true);
//makeskims d; d.run("TTTo2L2Nu", "/eos/uscms/store/user/fojensen/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_BTagAnalyzer_TTTo2L2Nu/180510_210135/", (831.76*.3258*.3258)/69705626., true);
//makeskims d; d.run("TTToHadronic", "/eos/uscms/store/user/fojensen/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_BTagAnalyzer_TTToHadronic/180510_210014/", (831.76*0.6741*0.6741)/129985840., true);

//makeskims d; d.run("SingleMuonB", "/eos/uscms/store/user/fojensen/SingleMuon/crab_BTagAnalyzer_SingleMuonB/180510_211932/");
//makeskims d; d.run("SingleMuonC", "/eos/uscms/store/user/fojensen/SingleMuon/crab_BTagAnalyzer_SingleMuonC/180510_211756/");
//makeskims d; d.run("SingleMuonD", "/eos/uscms/store/user/fojensen/SingleMuon/crab_BTagAnalyzer_SingleMuonD/180510_211529/");
//makeskims d; d.run("SingleMuonE", "/eos/uscms/store/user/fojensen/SingleMuon/crab_BTagAnalyzer_SingleMuonE/180510_211415/");
//makeskims d; d.run("SingleMuonF", "/eos/uscms/store/user/fojensen/SingleMuon/crab_BTagAnalyzer_SingleMuonF/180510_211317/");

//makeskims d; d.run("ST_tW_antitop_5f", "/eos/uscms/store/user/fojensen/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_BTagAnalyzer_ST_tW_antitop_5f/180510_210949/", (35.85*(.6741*.3258+.3258*.6741+.3258*.3258))/5365983., true);
//makeskims d; d.run("ST_tW_top_5f", "/eos/uscms/store/user/fojensen/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_BTagAnalyzer_ST_tW_top_5f/180510_211052/", (35.85*(.6741*.3258+.3258*.6741+.3258*.3258))/4656358., true);
//makeskims d; d.run("W3JetsToLNu", "/eos/uscms/store/user/fojensen/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/crab_BTagAnalyzer_W3JetsToLNu/180510_210600/", 942./19644745., true);
//makeskims d; d.run("W4JetsToLNu", "/eos/uscms/store/user/fojensen/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/crab_BTagAnalyzer_W4JetsToLNu/180510_210453/", 524.2/11285729., true);
//makeskims d; d.run("WWToLNuQQ_NNPDF31", "/eos/uscms/store/user/fojensen/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/crab_BTagAnalyzer_WWToLNuQQ_NNPDF31/180510_210254/", (117.8*0.6741*.3258)/8516920., true);

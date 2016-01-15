#include "TTbarEventAnalysis.h"
#include "TLorentzVector.h"

using namespace std;

//
void TTbarEventAnalysis::prepareOutput(TString outFile)
{

  //prepare output file
  outF_=TFile::Open(outFile,"RECREATE");

  //init KIN tree
  kinTree_=new TTree("kin","kinematics analysis");
  kinTree_->SetDirectory(outF_);
  kinTree_->Branch("EventInfo",    eventInfo_,         "EventInfo[3]/I");
  kinTree_->Branch("ttbar_chan",    &ttbar_chan_,      "ttbar_chan/I");
  kinTree_->Branch("npvn",    &npv_,      "npv/I");
  kinTree_->Branch("flavour",        jetFlavour_,      "flavour/I");
  kinTree_->Branch("jetmult",       &jetmult_,         "jetmult/I");
  kinTree_->Branch("jetpt",          jetPt_,           "jetpt/F");
  kinTree_->Branch("jeteta",         jetEta_,          "jeteta/F");
  kinTree_->Branch("jetrank",       &jetrank_,          "jetrank/I");       
  kinTree_->Branch("close_mlj",      close_mlj_,       "close_mlj[5]/F");
  kinTree_->Branch("close_deta",    &close_deta_,      "close_deta/F");
  kinTree_->Branch("close_dphi",    &close_dphi_,      "close_dphi/F");
  kinTree_->Branch("close_ptrel",   &close_ptrel_,     "close_ptrel/F");
  kinTree_->Branch("close_lj2ll_deta",    &close_lj2ll_deta_,      "close_lj2ll_deta/F");
  kinTree_->Branch("close_lj2ll_dphi",    &close_lj2ll_dphi_,      "close_lj2ll_dphi/F");
  kinTree_->Branch("far_mlj",       &far_mlj_,         "far_mlj/F");
  kinTree_->Branch("far_deta",      &far_deta_,        "far_deta/F");
  kinTree_->Branch("far_dphi",      &far_dphi_,        "far_dphi/F");
  kinTree_->Branch("far_ptrel",     &far_ptrel_,       "far_ptrel/F");
  kinTree_->Branch("far_lj2ll_deta",    &far_lj2ll_deta_,      "far_lj2ll_deta/F");
  kinTree_->Branch("far_lj2ll_dphi",    &far_lj2ll_dphi_,      "far_lj2ll_dphi/F");
  kinTree_->Branch("j2ll_deta",     &j2ll_deta_,       "j2ll_deta/F");
  kinTree_->Branch("j2ll_dphi",     &j2ll_dphi_,       "j2ll_dphi/F");
  kinTree_->Branch("kindisc",        kinDisc_,         "kindisc[5]/F");
  kinTree_->Branch("jp",             &(jp_[0]),              "jp/F");
  kinTree_->Branch("svhe",           &(svhe_[0]),            "svhe/F");
  kinTree_->Branch("csv",            &(csv_[0]),             "csv/F");
  kinTree_->Branch("weight",         weight_,          "weight[15]/F");

  ftmTree_=new TTree("ftm","flavour tag matching");
  ftmTree_->SetDirectory(outF_);
  ftmTree_->Branch("EventInfo",      eventInfo_,  "EventInfo[3]/I");
  ftmTree_->Branch("ttbar_chan",    &ttbar_chan_, "ttbar_chan/I");
  ftmTree_->Branch("jetmult",       &jetmult_,    "jetmult/I");
  ftmTree_->Branch("flavour",        jetFlavour_, "flavour[2]/I");
  ftmTree_->Branch("jetpt",          jetPt_,      "jetpt[2]/F");
  ftmTree_->Branch("jeteta",         jetEta_,     "jeteta[2]/F");
  ftmTree_->Branch("jp",             jp_,         "jp[2]/F");
  ftmTree_->Branch("svhe",           svhe_,       "svhe[2]/F");
  ftmTree_->Branch("csv",            csv_,        "csv[2]/F");
  ftmTree_->Branch("kindisc",        kinDisc_,    "kindisc[2]/F");
  ftmTree_->Branch("weight",         weight_,     "weight[15]/F");

  //prepare histograms
  std::map<TString,TH1F *> baseHistos;
  baseHistos["npvinc" ]  = new TH1F("npvinc", ";N_{PV,good}-N_{HS};Events",              50, 0, 50);
  baseHistos["npv"    ]  = new TH1F("npv",    ";N_{PV,good}-N_{HS};Events",              50, 0, 50);
  baseHistos["rho"    ]  = new TH1F("rho",    ";#rho [GeV];Events",                      50, 0, 30);
  baseHistos["mll"    ]  = new TH1F("mll",    ";Dilepton invariant mass [GeV];Events",   20, 0, 300);
  baseHistos["mllinc" ]  = new TH1F("mllinc", ";Dilepton invariant mass [GeV];Events",   20, 0, 300);
  baseHistos["met"    ]  = new TH1F("met",    ";Missing transverse energy [GeV];Events", 15, 0, 300);
  baseHistos["njets"  ]  = new TH1F("njets",  ";Jet multiplicity;Events;",               6,  2, 8);
  baseHistos["leadjpt"]  = new TH1F("leadjpt",";Leading jet p_{T} [GeV];Events;",        14,30,300);
  baseHistos["leadlpt"]  = new TH1F("leadlpt",";Leading lepton p_{T} [GeV];Events;",     9,20,200);
  baseHistos["trailjpt"] = new TH1F("trailjpt",";Trailing jet p_{T} [GeV];Events;",      14,30,300);
  baseHistos["traillpt"] = new TH1F("traillpt",";Trailing lepton p_{T} [GeV];Events;",   9,20,200);
  baseHistos["leadjeta"]    = new TH1F("leadjeta",    ";Pseudo-rapidity; Jets",              25, 0, 2.5);
  baseHistos["trailjeta"]   = new TH1F("trailjeta",    ";Pseudo-rapidity; Jets",              25, 0, 2.5);
  baseHistos["evsel"]    = new TH1F("evsel",   ";Event selection;Events;",               4,0,4);
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(1,"#geq 2j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(2,"=2j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(3,"=3j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(4,"=4j");
  baseHistos["jp"]=new TH1F("jp",";Jet probability;Jets",50,0,3);
  baseHistos["svhe"]=new TH1F("svhe",";Simple secondary vertex (HE);Jets",50,0,6);
  baseHistos["csv"]=new TH1F("csv",";Combined secondary vertex (IVF);Jets",50,0,1.1);
  baseHistos["tche"]=new TH1F("tche",";Track Counting High Efficiency;Jets",50,-20,50);
  baseHistos["jetseltrk"]=new TH1F("jetseltrk",";Selected track multiplicity;Jets",15,0,15);
  baseHistos["jp_leadkin"]=new TH1F("jp_leadkin",";Jet probability;Jets",50,0,3);
  baseHistos["svhe_leadkin"]=new TH1F("svhe_leadkin",";Simple secondary vertex (HE);Jets",50,0,6);
  baseHistos["csv_leadkin"]=new TH1F("csv_leadkin",";Combined secondary vertex (IVF);Jets",50,0,1.1);
  baseHistos["tche_leadkin"]=new TH1F("tche_leadkin",";Track Counting High Efficiency;Jets",50,-20,50);
  baseHistos["jetseltrk_leadkin"]=new TH1F("jetseltrk_leadkin",";Selected track multiplicity;Jets",15,0,15);
  baseHistos["flavour"] = new TH1F("flavour",     ";Jet flavour;Jets",                   4,  0, 4 );
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(1,"unmatched");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(2,"udsg");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(3,"c");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(4,"b");
  baseHistos["close_mlj"]   = new TH1F("close_mlj",   ";M(lepton,jet) [GeV]; Jets",          50, 0, 250);
  baseHistos["close_deta"]  = new TH1F("close_deta",  ";#Delta#eta(lepton,jet); Jets",       50, 0, 4);
  baseHistos["close_dphi"]  = new TH1F("close_dphi",  ";#Delta#phi(lepton,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["close_ptrel"] = new TH1F("close_ptrel", ";p_{T}^{rel}(lepton,jet) [GeV];Jets", 50, 0, 1);
  baseHistos["close_lj2ll_deta"] = new TH1F("close_lj2ll_deta",  ";#Delta#eta(lj,ll); Jets",       50, 0, 4);
  baseHistos["close_lj2ll_dphi"]  = new TH1F("close_l2jll_dphi",  ";#Delta#phi(lj,ll) [rad]; Jets", 50, 0, 3.15);
  baseHistos["far_mlj"]     = new TH1F("far_mlj",     ";M(lepton,jet) [GeV]; Jets",          50, 0, 250);
  baseHistos["far_deta"]    = new TH1F("far_deta",    ";#Delta#eta(lepton,jet); Jets",       50, 0, 4);
  baseHistos["far_dphi"]    = new TH1F("far_dphi",    ";#Delta#phi(lepton,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["far_ptrel"]   = new TH1F("far_ptrel",   ";p_{T}^{rel}(lepton,jet) [GeV];Jets", 50, 0, 1);
  baseHistos["far_lj2ll_deta"] = new TH1F("far_lj2ll_deta",  ";#Delta#eta(lepton+jet,ll); Jets",       50, 0, 4);
  baseHistos["far_lj2ll_dphi"]  = new TH1F("far_l2jll_dphi",  ";#Delta#phi(lj,ll) [rad]; Jets", 50, 0, 3.15);
  baseHistos["j2ll_deta"]    = new TH1F("j2ll_deta",    ";#Delta#eta(ll,jet); Jets",       50, 0, 4);
  baseHistos["j2ll_dphi"]    = new TH1F("j2ll_dphi",    ";#Delta#phi(lll,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["kindisc"]     = new TH1F("kindisc",     ";Kinematics discriminator;Jets",      100, -1, 1);

  //replicate histos per channel
  TString ch[]={"emu","ll","zll"};
  for(size_t i=0; i<sizeof(ch)/sizeof(TString); i++)
    {
      for(std::map<TString,TH1F *>::iterator it=baseHistos.begin(); it!=baseHistos.end(); it++)
	{
	  TString tag=ch[i]+"_"+it->first;
	  histos_[tag]=(TH1F *)it->second->Clone(tag);
	  histos_[tag]->Sumw2();
	  histos_[tag]->SetDirectory(outF_);
	}
    }

  histos_["puwgtnorm"] = new TH1F("puwgtnorm", ";puwgtnorm;Events",              4, 0, 4);
  histos_["puwgtnorm"]->Sumw2();
  histos_["puwgtnorm"]->SetDirectory(outF_);
}  

//
void TTbarEventAnalysis::processFile(TString inFile,TH1F *xsecWgt,Bool_t isData)
{
  //loop over events
  TFile *inF=TFile::Open(inFile);
  TTree *tree=(TTree *)inF->Get("btagana/ttree");
  Int_t nentries=tree->GetEntriesFast();
  std::cout << "...opening " << inFile << " -> analysing " << nentries << " events -> " << outF_->GetName();
  if(xsecWgt) std::cout << " xsec weight=" << xsecWgt->GetBinContent(1);
  if(isData)  std::cout << " is data";
  std::cout << std::endl;

  //prepare reader
  std::vector<Float_t> tmvaVars( tmvaVarNames_.size(), 0. );
  if(weightsDir_!="")
    {
      tmvaReader_=new TMVA::Reader( "!Color:!Silent" );
      for(size_t ivar=0; ivar<tmvaVarNames_.size(); ivar++)   
	tmvaReader_->AddVariable( tmvaVarNames_[ivar], &tmvaVars[ivar] );
      
      TString jranks[]={"leading",  "others",  "subleading" };
      for(size_t i=0; i<sizeof(jranks)/sizeof(TString); i++)
	tmvaReader_->BookMVA("BDT_"+jranks[i], weightsDir_+"/"+jranks[i]+"/TMVAClassification_BDT.weights.xml");
    }

  //prepare to read the tree (for jets only interested in a couple of variables)
  struct MyEventInfoBranches_t
  {
    Int_t Run,Evt,LumiBlock,nPV;
    Int_t   ttbar_chan, ttbar_trigWord, ttbar_metfilterWord;
    Int_t   ttbar_nl, ttbar_lid[10], ttbar_lgid[10], ttbar_lch[10];
    Float_t ttbar_lpt[10], ttbar_leta[10], ttbar_lphi[10], ttbar_lm[10];
    Float_t ttbar_metpt,ttbar_metphi;
    Float_t ttbar_rho;
    Int_t   ttbar_nw;
    Float_t nPUtrue;
    Float_t ttbar_w[500];
    Int_t nJet;
    Float_t Jet_pt[100],Jet_genpt[100],Jet_area[100],Jet_jes[100],Jet_eta[100],Jet_phi[100],Jet_mass[100];
    Float_t Jet_Svx[100],Jet_CombIVF[100],Jet_Proba[100],Jet_Ip2P[100];
    Int_t Jet_nseltracks[100];
    Int_t Jet_flavour[100];
  };
  MyEventInfoBranches_t ev;
  tree->SetBranchAddress("Run"        , &ev.Run        );
  tree->SetBranchAddress("Evt"        , &ev.Evt        );
  tree->SetBranchAddress("LumiBlock"  , &ev.LumiBlock  );
  tree->SetBranchAddress("nPV"        , &ev.nPV        );
  tree->SetBranchAddress("nPUtrue",     &ev.nPUtrue );
  tree->SetBranchAddress("ttbar_chan" , &ev.ttbar_chan);
  tree->SetBranchAddress("ttbar_metfilterWord", &ev.ttbar_metfilterWord);
  tree->SetBranchAddress("ttbar_trigWord", &ev.ttbar_trigWord);
  tree->SetBranchAddress("ttbar_nl"   ,  &ev.ttbar_nl);
  tree->SetBranchAddress("ttbar_lpt"  ,   ev.ttbar_lpt); 
  tree->SetBranchAddress("ttbar_leta" ,   ev.ttbar_leta);
  tree->SetBranchAddress("ttbar_lphi" ,   ev.ttbar_lphi);
  tree->SetBranchAddress("ttbar_lm"   ,   ev.ttbar_lm);
  tree->SetBranchAddress("ttbar_lid"  ,   ev.ttbar_lid);
  tree->SetBranchAddress("ttbar_lgid" ,   ev.ttbar_lgid);
  tree->SetBranchAddress("ttbar_lch"  ,   ev.ttbar_lch);
  tree->SetBranchAddress("ttbar_metpt",  &ev.ttbar_metpt);
  tree->SetBranchAddress("ttbar_metphi", &ev.ttbar_metphi);
  tree->SetBranchAddress("ttbar_rho",    &ev.ttbar_rho);
  tree->SetBranchAddress("ttbar_nw",     &ev.ttbar_nw);
  tree->SetBranchAddress("ttbar_w",      ev.ttbar_w);
  tree->SetBranchAddress("nJet",            &ev.nJet);
  tree->SetBranchAddress("Jet_pt",          ev.Jet_pt);
  tree->SetBranchAddress("Jet_genpt",       ev.Jet_genpt);
  tree->SetBranchAddress("Jet_area",        ev.Jet_area);
  tree->SetBranchAddress("Jet_jes",         ev.Jet_jes);
  tree->SetBranchAddress("Jet_eta",         ev.Jet_eta);
  tree->SetBranchAddress("Jet_phi",         ev.Jet_phi);
  tree->SetBranchAddress("Jet_mass",        ev.Jet_mass);
  tree->SetBranchAddress("Jet_Svx",         ev.Jet_Svx);
  tree->SetBranchAddress("Jet_CombIVF",     ev.Jet_CombIVF);
  tree->SetBranchAddress("Jet_Proba",       ev.Jet_Proba);
  tree->SetBranchAddress("Jet_Ip2P",        ev.Jet_Ip2P);
  tree->SetBranchAddress("Jet_nseltracks",  ev.Jet_nseltracks);
  tree->SetBranchAddress("Jet_flavour",     ev.Jet_flavour);

  for(Int_t i=0; i<nentries; i++)
    {
      tree->GetEntry(i);
      
      //progress bar
      if(i%100==0) std::cout << "\r[ " << int(100.*i/nentries) << "/100 ] to completion" << std::flush;
      
      //generator level weights
      Float_t genWgt=ev.ttbar_nw==0 ? 1.0 : ev.ttbar_w[0];
      Float_t qcdScaleLo(1.0),qcdScaleHi(1.0),hdampLo(1.0),hdampHi(1.0);
      if(readTTJetsGenWeights_ && ev.ttbar_nw>17)
      {
	qcdScaleLo=ev.ttbar_w[9]*(xsecWgt->GetBinContent(10)/xsecWgt->GetBinContent(1));
	qcdScaleHi=ev.ttbar_w[5]*(xsecWgt->GetBinContent(6)/xsecWgt->GetBinContent(1));
	hdampLo=ev.ttbar_w[ev.ttbar_nw-17]*(xsecWgt->GetBinContent(ev.ttbar_nw-17+1)/xsecWgt->GetBinContent(1));
	hdampHi=ev.ttbar_w[ev.ttbar_nw-9]*(xsecWgt->GetBinContent(ev.ttbar_nw-9+1)/xsecWgt->GetBinContent(1));
      }

      //pileup weights
      Float_t puWgtLo(1.0), puWgtNom(1.0), puWgtHi(1.0);
      if(!isData)
	{
	  if(puWgtGr_)     puWgtNom = puWgtGr_->Eval(ev.nPUtrue);
	  if(puWgtDownGr_) puWgtLo  = puWgtDownGr_->Eval(ev.nPUtrue);
	  if(puWgtUpGr_)   puWgtHi  = puWgtUpGr_->Eval(ev.nPUtrue);
	}
      histos_["puwgtnorm" ]->Fill(0.,1.0);
      histos_["puwgtnorm" ]->Fill(1.,puWgtNom);
      histos_["puwgtnorm" ]->Fill(2.,puWgtLo);
      histos_["puwgtnorm" ]->Fill(3.,puWgtHi);
        
      //
      //CHANNEL ASSIGNMENT 
      //
      if(ev.ttbar_nl<2 || ev.nJet<2) continue;
      ev.ttbar_chan=ev.ttbar_lid[0]*ev.ttbar_lch[0]*ev.ttbar_lid[1]*ev.ttbar_lch[1]; 

      TString ch("");
      if(ev.ttbar_chan==-11*13) ch="emu";
      if(ev.ttbar_chan==-11*11 || ev.ttbar_chan==-13*13) ch="ll";            
      if(ch=="") continue;

      //
      //TRIGGER
      //
      bool hasTrigger( triggerBits_.size()==0  ? true : false);
      for(size_t ibit=0; ibit<triggerBits_.size(); ibit++)
	{
	  if(triggerBits_[ibit].second!=ev.ttbar_chan) continue;
	  hasTrigger |= ((ev.ttbar_trigWord>>triggerBits_[ibit].first) & 1);
	}
      if(!hasTrigger) continue;

      //trigger efficiency weight
      Float_t trigWgtLo(1.0), trigWgtNom(1.0), trigWgtHi(1.0);
      if(!isData)
	{
	  std::pair<float,float> eff=getTriggerEfficiency(ev.ttbar_lid[0],ev.ttbar_lpt[0],ev.ttbar_leta[0],
							  ev.ttbar_lid[1],ev.ttbar_lpt[1],ev.ttbar_leta[1],
							  ev.ttbar_chan);
	  trigWgtLo=eff.first-eff.second;
	  trigWgtNom=eff.first;
	  trigWgtHi=eff.first+eff.second;
	}

      //lepton selection efficiency
      Float_t lepSelEffLo(1.0), lepSelEffNom(1.0), lepSelEffHi(1.0);
      if(!isData)
	{
	  for(size_t il=0; il<2; il++)
	    {
	      std::pair<float,float> lepSF = getLeptonSelectionEfficiencyScaleFactor(ev.ttbar_lid[il],ev.ttbar_lpt[il],ev.ttbar_leta[il]);
	      lepSelEffLo  *= (lepSF.first-lepSF.second);
	      lepSelEffNom *= lepSF.first;
	      lepSelEffHi  *= (lepSF.first+lepSF.second);
	    }
	}

      //dilepton invariant mass
      std::vector<TLorentzVector> lp4;
      for(Int_t il=0; il<ev.ttbar_nl; il++)
	{
	  lp4.push_back( TLorentzVector(0,0,0,0) );
          lp4[il].SetPtEtaPhiM(ev.ttbar_lpt[il],ev.ttbar_leta[il],ev.ttbar_lphi[il],0.);
	}

      TLorentzVector dilepton(lp4[0]+lp4[1]);
      Float_t mll=dilepton.M();
      if(mll<12) continue;

      //nominal event weight
      Float_t evWgt(1.0);
      if(!isData){
	evWgt *= puWgtNom*trigWgtNom*lepSelEffNom*genWgt;
	if(xsecWgt) evWgt *= xsecWgt->GetBinContent(1);
      }
      histos_[ch+"_npvinc"]->Fill(ev.nPV-1,evWgt);
      npv_=ev.nPV;

      //
      //JET/MET SELECTION
      //
      Int_t jetCount[5]={0,0,0,0,0};
      std::vector<Int_t> selJets;
      std::vector<std::vector<Float_t> > selJetsKINDisc;
      std::vector< std::vector<TLorentzVector> > selJetsP4;
      std::vector< std::vector< std::vector<LJKinematics_t> > > selJetsLJKinematics;
      for(Int_t ij=0; ij<ev.nJet; ij++)
	{      
	  //convert to P4
	  TLorentzVector jp4(0,0,0,0);
	  jp4.SetPtEtaPhiM(ev.Jet_pt[ij],ev.Jet_eta[ij],ev.Jet_phi[ij],ev.Jet_mass[ij]);

	  //cross clean wrt to leptons
	  Float_t minDRlj(9999.);
	  for(size_t il=0; il<2; il++) minDRlj = TMath::Min( (Float_t)minDRlj, (Float_t)lp4[il].DeltaR(jp4) );
	  if(minDRlj<0.4) continue;
	  
	  //update jet energy scale/resolution
	  Float_t jrawsf=1./ev.Jet_jes[ij];
	  Float_t jarea=ev.Jet_area[ij];
	  Float_t genjpt=ev.Jet_genpt[ij];
      
	  // update JES+JER for this jet
	  std::vector<float> jesSF(3,1.0);
	  jecUnc_->setJetEta(fabs(jp4.Eta()));
	  jecUnc_->setJetPt(jp4.Pt());
	  float unc = jecUnc_->getUncertainty(true);
	  jesSF[1]=(1.+fabs(unc));
	  jesSF[2]=(1.-fabs(unc));
	  
	  std::vector<float> jerSF= getJetResolutionScales(jesSF[0]*jp4.Pt(), jp4.Eta(), genjpt);

	  TLorentzVector oldjp4(jp4);
	  jp4 = jp4*jesSF[0]*jerSF[0];
	  
	  // apply energy shifts according to systematic variation
	  Bool_t canBeSelected(false);
	  std::vector<TLorentzVector> varjp4;
	  std::vector<Float_t> varkindisc;
	  std::vector< std::vector<LJKinematics_t> > varLJKinematics;
	  for(Int_t iSystVar=0; iSystVar<5; iSystVar++)
	    {
	      varjp4.push_back( jp4 );
	      if(iSystVar==1) varjp4[iSystVar] *= jesSF[1]/jesSF[0];
	      if(iSystVar==2) varjp4[iSystVar] *= jesSF[2]/jesSF[0];
	      if(iSystVar==3) varjp4[iSystVar] *= jerSF[1]/jerSF[0];
	      if(iSystVar==4) varjp4[iSystVar] *= jerSF[2]/jerSF[0];
	          
	      //prepare variables for MVA
	      std::vector< LJKinematics_t > ljkinematics;
	      for(Int_t il=0; il<2; il++)
		{
		  LJKinematics_t iljkin;
		  iljkin.dr         = lp4[il].DeltaR(varjp4[iSystVar]);
		  iljkin.dphi       = fabs(lp4[il].DeltaPhi(varjp4[iSystVar]));
		  iljkin.deta       = fabs(lp4[il].Eta()-varjp4[iSystVar].Eta());
		  iljkin.ptrel      = ROOT::Math::VectorUtil::Perp(lp4[il].Vect(),varjp4[iSystVar].Vect().Unit())/lp4[il].P();
		  TLorentzVector ljP4(lp4[il]+varjp4[iSystVar]);
		  iljkin.mlj        = ljP4.M();
		  iljkin.lj2ll_deta = fabs(ljP4.Eta()-dilepton.Eta());
		  iljkin.lj2ll_dphi = fabs(ljP4.DeltaPhi(dilepton));		  
		  ljkinematics.push_back(iljkin);
		}
	      sort(ljkinematics.begin(),ljkinematics.end(),sortLJKinematicsByDR);
	      varLJKinematics.push_back(ljkinematics);
      
	      //evaluate the MVA
	      Float_t kindisc(0.0);
	      if(tmvaReader_)
		{
		  for(size_t ivar=0; ivar<tmvaVarNames_.size(); ivar++)
		    {
		      if( tmvaVarNames_[ivar].Contains("j2ll_") )
			{
			  if(tmvaVarNames_[ivar]=="j2ll_deta") tmvaVars[ivar]=fabs(varjp4[iSystVar].Eta()-dilepton.Eta());
			  if(tmvaVarNames_[ivar]=="j2ll_phi")  tmvaVars[ivar]=fabs(varjp4[iSystVar].DeltaPhi(dilepton));
			}
		      else
			{
			  int ljidx( tmvaVarNames_[ivar].Contains("close") ? 0 : 1);
			  if(tmvaVarNames_[ivar].Contains("_dr"))    tmvaVars[ivar]=ljkinematics[ljidx].dr;			  
			  if(tmvaVarNames_[ivar].Contains("_dphi"))  
			    {
			      if(tmvaVarNames_[ivar].Contains("lj2ll_")) tmvaVars[ivar]=ljkinematics[ljidx].lj2ll_dphi;
			      else                                       tmvaVars[ivar]=ljkinematics[ljidx].dphi;
			    }
			  if(tmvaVarNames_[ivar].Contains("_deta"))  
			    {
			      if(tmvaVarNames_[ivar].Contains("lj2ll_")) tmvaVars[ivar]=ljkinematics[ljidx].lj2ll_deta;
			      else                                       tmvaVars[ivar]=ljkinematics[ljidx].deta;
			    }			      
			  if(tmvaVarNames_[ivar].Contains("_ptrel")) tmvaVars[ivar]=ljkinematics[ljidx].ptrel;
			  if(tmvaVarNames_[ivar].Contains("_mlj"))   tmvaVars[ivar]=ljkinematics[ljidx].mlj;
			}
		    }

		  TString methodPFix("_others");
		  if(selJets.size()==0) methodPFix="_leading";
		  else if(selJets.size()==1) methodPFix="_subleading";
		  varkindisc.push_back( tmvaReader_->EvaluateMVA("BDT"+methodPFix) );
		}	      
	      
	      //check if can be selected for this variation
	      if(varjp4[iSystVar].Pt()<30 || TMath::Abs(varjp4[iSystVar].Eta())>2.5) continue;
	      canBeSelected=true;
	      jetCount[iSystVar]++;
	    }

	  //add jet if it is selectable
	  if(!canBeSelected) continue;
	  selJets.push_back(ij);
	  selJetsP4.push_back(varjp4);
	  selJetsLJKinematics.push_back( varLJKinematics );
	  if(tmvaReader_) selJetsKINDisc.push_back(varkindisc);
	}
      
      //
      // base selection and n-1 plots
      //
      bool zCand( (ch.Contains("ll") && TMath::Abs(mll-91)<15) ? true : false );
      bool passMet( (ch.Contains("emu") || ev.ttbar_metpt>40) ?  true : false);
      bool passJets(selJets.size()>=2 ? true : false);
      
      if(!passJets) continue;
      if(zCand)  ch="z"+ch;
      histos_[ch+"_mllinc"]->Fill(mll,evWgt);
      histos_[ch+"_met"]->Fill(ev.ttbar_metpt,evWgt);

      if(!passMet) continue;
      histos_[ch+"_evsel"]->Fill(0.,evWgt);
      if(selJets.size()<5)   histos_[ch+"_evsel"]->Fill(selJets.size()-1,evWgt);
      histos_[ch+"_rho"]->Fill(ev.ttbar_rho,evWgt);
      histos_[ch+"_npv"]->Fill(ev.nPV-1,evWgt);       
      histos_[ch+"_mll"]->Fill(mll,evWgt);
      histos_[ch+"_njets"]->Fill(selJets.size(),evWgt);
      histos_[ch+"_leadjpt"]->Fill(selJetsP4[0][0].Pt(),evWgt);
      histos_[ch+"_leadjeta"]->Fill((selJetsP4[0][0].Eta()),evWgt);
      histos_[ch+"_leadlpt"]->Fill(lp4[0].Pt(),evWgt);
      histos_[ch+"_trailjpt"]->Fill(selJetsP4[1][0].Pt(),evWgt);
      histos_[ch+"_trailjeta"]->Fill(fabs(selJetsP4[1][0].Eta()),evWgt);
      histos_[ch+"_traillpt"]->Fill(lp4[1].Pt(),evWgt);

      std::vector<float> leadingkindisc(2,-9999);
      std::vector<int> leadingkindiscIdx(2,-1);
      for(size_t ij=0; ij<selJets.size(); ij++)
	{
	  Int_t jetIdx(selJets[ij]);
	  histos_[ch+"_close_mlj"]->Fill(selJetsLJKinematics[ij][0][0].mlj,evWgt);
	  histos_[ch+"_close_deta"]->Fill(selJetsLJKinematics[ij][0][0].deta,evWgt);
	  histos_[ch+"_close_dphi"]->Fill(selJetsLJKinematics[ij][0][0].dphi,evWgt);
	  histos_[ch+"_close_ptrel"]->Fill(selJetsLJKinematics[ij][0][0].ptrel,evWgt);
	  histos_[ch+"_close_lj2ll_deta"]->Fill(selJetsLJKinematics[ij][0][0].lj2ll_deta,evWgt);
	  histos_[ch+"_close_lj2ll_dphi"]->Fill(selJetsLJKinematics[ij][0][0].lj2ll_dphi,evWgt);
	  histos_[ch+"_far_mlj"]->Fill(selJetsLJKinematics[ij][0][1].mlj,evWgt);
	  histos_[ch+"_far_deta"]->Fill(selJetsLJKinematics[ij][0][1].deta,evWgt);
	  histos_[ch+"_far_dphi"]->Fill(selJetsLJKinematics[ij][0][1].dphi,evWgt);
	  histos_[ch+"_far_ptrel"]->Fill(selJetsLJKinematics[ij][0][1].ptrel,evWgt);
	  histos_[ch+"_far_lj2ll_deta"]->Fill(selJetsLJKinematics[ij][0][1].lj2ll_deta,evWgt);
	  histos_[ch+"_far_lj2ll_dphi"]->Fill(selJetsLJKinematics[ij][0][1].lj2ll_dphi,evWgt);
	  histos_[ch+"_j2ll_deta"]->Fill(fabs(selJetsP4[ij][0].Eta()-dilepton.Eta()),evWgt);
	  histos_[ch+"_j2ll_dphi"]->Fill(fabs(selJetsP4[ij][0].DeltaPhi(dilepton)),evWgt);
	  if(tmvaReader_) histos_[ch+"_kindisc"]->Fill(selJetsKINDisc[ij][0],evWgt);
	  histos_[ch+"_jp"]->Fill(ev.Jet_Proba[jetIdx],evWgt);
	  histos_[ch+"_svhe"]->Fill(ev.Jet_Svx[jetIdx],evWgt);
	  histos_[ch+"_csv"]->Fill(ev.Jet_CombIVF[jetIdx],evWgt);
	  histos_[ch+"_tche"]->Fill(ev.Jet_Ip2P[jetIdx],evWgt);
	  histos_[ch+"_jetseltrk"]->Fill(ev.Jet_nseltracks[jetIdx],evWgt);
	  

	  Int_t flavBin(0),partonFlav(abs(ev.Jet_flavour[jetIdx]));
	  if(partonFlav==21 || (partonFlav>0 && partonFlav<4)) flavBin=1;
	  if(partonFlav==4) flavBin=2;
	  if(partonFlav==5) flavBin=3;
	  histos_[ch+"_flavour"]->Fill(flavBin,evWgt);

	  //rank jets by kinematics discriminator
	  if(tmvaReader_)
	    {
	      if(selJetsKINDisc[ij][0]>leadingkindisc[0])
		{
		  leadingkindisc[1]=leadingkindisc[0];     leadingkindiscIdx[1]=leadingkindiscIdx[0];
		  leadingkindisc[0]=selJetsKINDisc[ij][0]; leadingkindiscIdx[0]=jetIdx;
		}
	      else if(selJetsKINDisc[ij][0]>leadingkindisc[1])
		{
		  leadingkindisc[1]=selJetsKINDisc[ij][0]; leadingkindiscIdx[1]=jetIdx;
		}
	    }
	}
      
      //control b-tagging quantities for the most promising jets in the KIN discriminator
      if(tmvaReader_)
	{
	  for(size_t ij=0; ij<2; ij++)
            {
              size_t jetIdx=leadingkindiscIdx[ij];
	      histos_[ch+"_jp_leadkin"]->Fill(ev.Jet_Proba[jetIdx],evWgt);
	      histos_[ch+"_svhe_leadkin"]->Fill(ev.Jet_Svx[jetIdx],evWgt);
	      histos_[ch+"_csv_leadkin"]->Fill(ev.Jet_CombIVF[jetIdx],evWgt);
	      histos_[ch+"_tche_leadkin"]->Fill(ev.Jet_Ip2P[jetIdx],evWgt);
	      histos_[ch+"_jetseltrk_leadkin"]->Fill(ev.Jet_nseltracks[jetIdx],evWgt);
	    }
	}

      //
      //prepare to store trees
      //
      eventInfo_[0]=ev.Run;
      eventInfo_[1]=ev.Evt;
      eventInfo_[2]=ev.LumiBlock;

      jetmult_=selJets.size();
      ttbar_chan_=ev.ttbar_chan;
      if(zCand) ttbar_chan_+= 230000;

      //weights for systematic uncertainties
      for(Int_t iSystVar=0; iSystVar<5; iSystVar++)
	{
	  Float_t selWeight(jetCount[iSystVar]>=2 ? 1.0 : 0.0);
	  weight_[iSystVar]=evWgt*selWeight;
	}
      weight_[5] = puWgtNom>0 ? evWgt*puWgtLo/puWgtNom : evWgt;
      weight_[6] = puWgtLo>0  ? evWgt*puWgtHi/puWgtNom : evWgt;
      weight_[7] = evWgt*trigWgtLo/trigWgtNom;
      weight_[8] = evWgt*trigWgtHi/trigWgtNom;
      weight_[9] = evWgt*lepSelEffLo/lepSelEffNom;
      weight_[10]= evWgt*lepSelEffHi/lepSelEffNom;
      weight_[11]= evWgt*qcdScaleLo/genWgt;
      weight_[12]= evWgt*qcdScaleHi/genWgt;
      weight_[13]= evWgt*hdampLo/genWgt;
      weight_[14]= evWgt*hdampHi/genWgt;

      //fill trees
      for(size_t ij=0; ij<selJets.size(); ij++)
	{
	  Int_t jetIdx(selJets[ij]);

	  jetrank_ = ij;
	  jetFlavour_[0] = ev.Jet_flavour[jetIdx];
	  jetPt_[0]      = selJetsP4[ij][0].Pt();
	  jetEta_[0]     = selJetsP4[ij][0].Eta();
	  for(size_t iSystVar=0; iSystVar<5; iSystVar++)
	    {
	      close_mlj_[iSystVar] = selJetsLJKinematics[ij][iSystVar][0].mlj;
	      if(tmvaReader_) kinDisc_[iSystVar] = selJetsKINDisc[ij][iSystVar];
	      else            kinDisc_[iSystVar] = -999;
	    }
	  close_deta_ =selJetsLJKinematics[ij][0][0].deta;
	  close_dphi_ =selJetsLJKinematics[ij][0][0].dphi;
	  close_ptrel_=selJetsLJKinematics[ij][0][0].ptrel;
	  close_lj2ll_deta_ = selJetsLJKinematics[ij][0][0].lj2ll_deta;
	  close_lj2ll_dphi_ = selJetsLJKinematics[ij][0][0].lj2ll_dphi;

	  far_mlj_    =selJetsLJKinematics[ij][0][1].mlj;
	  far_deta_   =selJetsLJKinematics[ij][0][1].deta;
	  far_dphi_   =selJetsLJKinematics[ij][0][1].dphi;
	  far_ptrel_  =selJetsLJKinematics[ij][0][1].ptrel;
	  far_lj2ll_deta_ = selJetsLJKinematics[ij][0][1].lj2ll_deta;
	  far_lj2ll_dphi_ = selJetsLJKinematics[ij][0][1].lj2ll_dphi;	  

	  j2ll_deta_  =fabs(selJetsP4[ij][0].Eta()-dilepton.Eta());
	  j2ll_dphi_  =fabs(selJetsP4[ij][0].DeltaPhi( dilepton ));
	  
	  jp_[0]=ev.Jet_Proba[jetIdx];
	  svhe_[0]=ev.Jet_Svx[jetIdx];
	  csv_[0]=ev.Jet_CombIVF[jetIdx];

	  kinTree_->Fill();
	}

      //FtM tree is filled with the two leading jets in the KIN discriminator
      if(tmvaReader_)
	{
	  for(size_t ij=0; ij<2; ij++)
	    {
	      size_t jetIdx=leadingkindiscIdx[ij];
	      jetFlavour_[ij] = ev.Jet_flavour[jetIdx];
	      jetPt_[ij]      = ev.Jet_pt[jetIdx];
	      jetEta_[ij]     = ev.Jet_eta[jetIdx];
	      jp_[ij]         = ev.Jet_Proba[jetIdx];
	      svhe_[ij]       = ev.Jet_Svx[jetIdx];
	      csv_[ij]        = ev.Jet_CombIVF[jetIdx];
	      kinDisc_[ij]    = leadingkindisc[ij];
	      //std::cout << ij << " " <<  jetFlavour_[ij] << " "
	      //<< jetPt_[ij] << " " << csv_[ij]  << " " << kinDisc_[ij] 
	      //		<< std::endl;
	    }
	  ftmTree_->Fill();
	}
    }

  //all done with this file
  inF->Close();
}

//Sources                                                                                                                                                          
// CMS AN 022/2015 v15                                                                                                                                     
// https://indico.cern.ch/event/434078/#preview:1614815  
std::pair<float,float> TTbarEventAnalysis::getTriggerEfficiency(int id1,float pt1,float eta1,int id2,float pt2,float eta2,int ch)
{
  std::pair<float,float>res(1.0,0.0);
  if(ch==-11*13) { res.first=1.0; res.second=0.05; }
  if(ch==-11*11) { res.first=1.0; res.second=0.05; }
  if(ch==-13*13) { res.first=1.0; res.second=0.05; } // this is a single muon trigger
 return res;
}

//Sources
// CMS AN 022/2015 v15
std::pair<float,float> TTbarEventAnalysis::getLeptonSelectionEfficiencyScaleFactor(int id,float pt,float eta)
{
  std::pair<float,float>res(1.0,0.0);
 
  //electrons
  if(abs(id)==11)
    {
      if (fabs(eta)<0.8)
	{
	  if (pt<30)      { res.first=0.927; res.second=0.073; }
	  else if (pt<40) { res.first=0.975; res.second=0.018; }
	  else if (pt<50) { res.first=0.962; res.second=0.036; }
	  else            { res.first=0.955; res.second=0.022; }
	}
      else if (fabs(eta)<1.5)
	{
	  if (pt<30)      { res.first=0.891; res.second=0.074; }
	  else if (pt<40) { res.first=0.965; res.second=0.020; }
	  else if (pt<50) { res.first=0.968; res.second=0.018; }
	  else            { res.first=0.955; res.second=0.018; }
	}
      else
	{
	  if (pt<30)      { res.first=0.956; res.second=0.059; }
	  else if (pt<40) { res.first=0.995; res.second=0.018; }
	  else if (pt<50) { res.first=0.993; res.second=0.019; }
	  else            { res.first=0.985; res.second=0.023; }
	}
    }

  //muons
  if (abs(id)==13)
    {
      if (fabs(eta)<0.9)
	{
	  if (pt<30)      { res.first=1.003; res.second=0.019; }
	  else if (pt<40) { res.first=1.014; res.second=0.015; }
	  else if (pt<50) { res.first=1.001; res.second=0.014; }
	  else            { res.first=0.983; res.second=0.014; }
	}
      else if(fabs(eta)<1.2)
	{
	  if (pt<30)      { res.first=0.993; res.second=0.019; }
	  else if (pt<40) { res.first=0.994; res.second=0.015; }
	  else if (pt<50) { res.first=0.980; res.second=0.014; }
	  else            { res.first=0.987; res.second=0.015; }
	}
      else
	{
	  if (pt<30)      { res.first=1.023; res.second=0.028; }
	  else if (pt<40) { res.first=0.994; res.second=0.014; }
	  else if (pt<50) { res.first=0.996; res.second=0.014; }
	  else            { res.first=0.979; res.second=0.014; }
	}
    }

  return res;
}

//Sources
//  Assuming nominal JER but uncertainties from Run I
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> TTbarEventAnalysis::getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.06);
  if(TMath::Abs(eta)<0.5) 
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2));
    }
  else if(TMath::Abs(eta)<1.1)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2));
    }
  else if(TMath::Abs(eta)<1.7)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2));
    }
  else if(TMath::Abs(eta)<2.3)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    }
  else
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2));
    }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}


//
void TTbarEventAnalysis::finalizeOutput()
{
  //dump results to file
  outF_->cd();
  
  //pileup weighting screws up a bit normalization - fix it a posteriori
  float puwgtSF(histos_["puwgtnorm" ]->GetBinContent(1)/histos_["puwgtnorm" ]->GetBinContent(2));

  for(std::map<TString,TH1F *>::iterator it = histos_.begin(); it != histos_.end(); it++) 
    {
      if(it->first!="puwgtnorm") 
	it->second->Scale(puwgtSF);
      it->second->Write();
    }
  kinTree_->Write();
  if(ftmTree_) ftmTree_->Write();
  outF_->Close();
}

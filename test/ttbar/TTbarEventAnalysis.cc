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
  kinTree_->Branch("flavour",        jetFlavour_,      "flavour/I");
  kinTree_->Branch("jetmult",       &jetmult_,         "jetmult/I");
  kinTree_->Branch("jetpt",          jetPt_,           "jetpt/F");
  kinTree_->Branch("jeteta",         jetEta_,          "jeteta/F");
  kinTree_->Branch("close_mlj",      close_mlj_,       "close_mlj[5]/F");
  kinTree_->Branch("close_deta",    &close_deta_,      "close_deta/F");
  kinTree_->Branch("close_dphi",    &close_dphi_,      "close_dphi/F");
  kinTree_->Branch("close_ptrel",   &close_ptrel_,     "close_ptrel/F");
  kinTree_->Branch("far_mlj",       &far_mlj_,         "far_mlj/F");
  kinTree_->Branch("far_deta",      &far_deta_,        "far_deta/F");
  kinTree_->Branch("far_dphi",      &far_dphi_,        "far_dphi/F");
  kinTree_->Branch("far_ptrel",     &far_ptrel_,       "far_ptrel/F");
  kinTree_->Branch("kindisc",        kinDisc_,         "kindisc[5]/F");
  kinTree_->Branch("jp",             jp_,              "jp/F");
  kinTree_->Branch("svhe",           svhe_,            "svhe/F");
  kinTree_->Branch("csv",            csv_,             "csv/F");
  kinTree_->Branch("weight",         weight_,          "weight[15]/F");

  ftmTree_=new TTree("ftm","flavour tag matching");
  ftmTree_->SetDirectory(outF_);
  ftmTree_->Branch("EventInfo",      eventInfo_,  "EventInfo[3]/I");
  ftmTree_->Branch("ttbar_chan",    &ttbar_chan_, "ttbar_chan/I");
  ftmTree_->Branch("jetmult",       &jetmult_,    "jetmult/I");
  ftmTree_->Branch("flavour",        jetFlavour_, "flavour[2]/I");
  ftmTree_->Branch("jetpt",          jetPt_,      "jetpt[2]/F");
  ftmTree_->Branch("jeteta",         jetEta_,     "jeteta[2]/F");
  kinTree_->Branch("jp",             jp_,         "jp[2]/F");
  kinTree_->Branch("svhe",           svhe_,       "svhe[2]/F");
  kinTree_->Branch("csv",            csv_,        "csv[2]/F");
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
  baseHistos["evsel"]    = new TH1F("evsel",   ";Event selection;Events;",               4,0,4);
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(1,"#geq 2j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(2,"=2j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(3,"=3j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(4,"=4j");
  baseHistos["jp"]=new TH1F("jp",";Jet probability;Jets",50,0,3);
  baseHistos["svhe"]=new TH1F("svhe",";Simple secondary vertex (HE);Jets",50,0,6);
  baseHistos["csv"]=new TH1F("csv",";Combined secondary vertex (IVF);Jets",50,0,1.1);
  baseHistos["flavour"]     = new TH1F("flavour",     ";Jet flavour;Jets",                   4,  0, 4 );
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(1,"unmatched");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(2,"udsg");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(3,"c");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(4,"b");
  baseHistos["jeteta"]      = new TH1F("jeteta",      ";Pseudo-rapidity; Jets",              25, 0, 2.5);
  baseHistos["close_mlj"]   = new TH1F("close_mlj",   ";M(lepton,jet) [GeV]; Jets",          50, 0, 250);
  baseHistos["close_deta"]  = new TH1F("close_deta",  ";#Delta#eta(lepton,jet); Jets",       50, 0, 4);
  baseHistos["close_dphi"]  = new TH1F("close_dphi",  ";#Delta#phi(lepton,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["close_ptrel"] = new TH1F("close_ptrel", ";p_{T}^{rel}(lepton,jet) [GeV];Jets", 50, 0, 1);
  baseHistos["far_mlj"]     = new TH1F("far_mlj",     ";M(lepton,jet) [GeV]; Jets",          50, 0, 250);
  baseHistos["far_deta"]    = new TH1F("far_deta",    ";#Delta#eta(lepton,jet); Jets",       50, 0, 4);
  baseHistos["far_dphi"]    = new TH1F("far_dphi",    ";#Delta#phi(lepton,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["far_ptrel"]   = new TH1F("far_ptrel",   ";p_{T}^{rel}(lepton,jet) [GeV];Jets", 50, 0, 1);
  baseHistos["kindisc"]     = new TH1F("kindisc",     ";Kinematics discriminator;Jets",      50, -1, 1);

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
}  

//
void TTbarEventAnalysis::processFile(TString inFile,float xsecWgt)
{
  //loop over events
  TFile *inF=TFile::Open(inFile);
  TTree *tree=(TTree *)inF->Get("btagana/ttree");
  Int_t nentries=tree->GetEntriesFast();
  std::cout << "...opening " << inFile << " -> analysing " << nentries << " events -> " << outF_->GetName() << std::endl;

  //prepare reader
  std::vector<Float_t> tmvaVars( tmvaVarNames_.size(), 0. );
  if(weightsFile_!="")
    {
      tmvaReader_=new TMVA::Reader( "!Color:!Silent" );
      for(size_t ivar=0; ivar<tmvaVarNames_.size(); ivar++)   
	tmvaReader_->AddVariable( tmvaVarNames_[ivar], &tmvaVars[ivar] );
      tmvaReader_->BookMVA("BDT", weightsFile_);
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
    Float_t ttbar_w[500];
    Int_t nJet;
    Float_t Jet_pt[100],Jet_genpt[100],Jet_area[100],Jet_jes[100],Jet_eta[100],Jet_phi[100],Jet_mass[100];
    Float_t Jet_Svx[100],Jet_CombIVF[100],Jet_Proba[100];
    Int_t Jet_flavour[100];
  };
  MyEventInfoBranches_t ev;
  tree->SetBranchAddress("Run"        , &ev.Run        );
  tree->SetBranchAddress("Evt"        , &ev.Evt        );
  tree->SetBranchAddress("LumiBlock"  , &ev.LumiBlock  );
  tree->SetBranchAddress("nPV"        , &ev.nPV        );
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
  tree->SetBranchAddress("Jet_flavour",     ev.Jet_flavour);

  for(Int_t i=0; i<nentries; i++)
    {
      tree->GetEntry(i);
      
      //progress bar
      if(i%100==0) std::cout << "\r[ " << int(100.*i/nentries) << "/100 ] to completion" << std::flush;

      //generator level weights
      Float_t genWgt=ev.ttbar_nw==0 ? 1.0 : ev.ttbar_w[0];
      if(useOnlySignOfGenWeight_) genWgt=ev.ttbar_w[0]<0 ? -1.0 : 1.0;
      Float_t qcdScaleLo(1.0),qcdScaleHi(1.0),hdampLo(1.0),hdampHi(1.0);
      if(readTTJetsGenWeights_ && ev.ttbar_nw>17)
      {
	qcdScaleLo=ev.ttbar_w[9];
	qcdScaleHi=ev.ttbar_w[5];
	hdampLo=ev.ttbar_w[ev.ttbar_nw-17];
	hdampHi=ev.ttbar_w[ev.ttbar_nw-9];
      }

      //pileup weights
      Float_t puWgtLo(1.0), puWgtNom(1.0), puWgtHi(1.0);

      //
      //MET FILTERS
      //
      if(applyMETFilters_ && ev.ttbar_metfilterWord!=0) continue;
        
      //
      //CHANNEL ASSIGNMENT 
      //
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


      //
      //LEPTONS
      //
      if(ev.ttbar_nl<2) continue;

      //trigger efficiency weight
      Float_t trigWgtLo(1.0), trigWgtNom(1.0), trigWgtHi(1.0);
      if(applyTriggerEff_)
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
      if(applyLepSelEff_)
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
      Float_t mll=(lp4[0]+lp4[1]).M();

      //nominal event weight
      Float_t evWgt(puWgtNom*trigWgtNom*lepSelEffNom*genWgt);
      histos_[ch+"_npvinc"]->Fill(ev.nPV-1,evWgt);

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
	  std::vector<float> jesSF= getJetEnergyScales(jp4.Pt(), jp4.Eta(), jrawsf,jarea,ev.ttbar_rho);
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
		  iljkin.dr=lp4[il].DeltaR(varjp4[iSystVar]);
		  iljkin.dphi=TMath::Abs(lp4[il].DeltaPhi(varjp4[iSystVar]));
		  iljkin.deta=TMath::Abs(lp4[il].Eta()-varjp4[iSystVar].Eta());
		  iljkin.ptrel=ROOT::Math::VectorUtil::Perp(lp4[il].Vect(),varjp4[iSystVar].Vect().Unit())/lp4[il].P();
		  iljkin.mlj=(varjp4[iSystVar]+lp4[il]).M();
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
		      int ljidx( tmvaVarNames_[ivar].Contains("close") ? 0 : 1);
		      if(tmvaVarNames_[ivar].Contains("_dr"))    tmvaVars[ivar]=ljkinematics[ljidx].dr;
		      if(tmvaVarNames_[ivar].Contains("_dphi"))  tmvaVars[ivar]=ljkinematics[ljidx].dphi;
		      if(tmvaVarNames_[ivar].Contains("_deta"))  tmvaVars[ivar]=ljkinematics[ljidx].deta;
		      if(tmvaVarNames_[ivar].Contains("_ptrel")) tmvaVars[ivar]=ljkinematics[ljidx].ptrel;
		      if(tmvaVarNames_[ivar].Contains("_mlj"))   tmvaVars[ivar]=ljkinematics[ljidx].mlj;
		    }
		  varkindisc.push_back( tmvaReader_->EvaluateMVA("BDT") );
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
      histos_[ch+"_mllinc"]->Fill(mll,evWgt);
      if(passMet) histos_[ch+"_mll"]->Fill(mll,evWgt);
      if(zCand) ch="z"+ch;
      histos_[ch+"_met"]->Fill(ev.ttbar_metpt,evWgt);
      if(!passMet) continue;
      histos_[ch+"_evsel"]->Fill(0.,evWgt);
      if(selJets.size()<5)   histos_[ch+"_evsel"]->Fill(selJets.size()-1,evWgt);
      histos_[ch+"_rho"]->Fill(ev.ttbar_rho,evWgt);
      histos_[ch+"_npv"]->Fill(ev.nPV-1,evWgt);       
      histos_[ch+"_njets"]->Fill(selJets.size(),evWgt);
      histos_[ch+"_leadjpt"]->Fill(selJetsP4[0][0].Pt(),evWgt);
      histos_[ch+"_leadlpt"]->Fill(lp4[0].Pt(),evWgt);
      histos_[ch+"_trailjpt"]->Fill(selJetsP4[1][0].Pt(),evWgt);
      histos_[ch+"_traillpt"]->Fill(lp4[1].Pt(),evWgt);
      
      for(size_t ij=0; ij<selJets.size(); ij++)
	{
	  Int_t jetIdx(selJets[ij]);
	  histos_[ch+"_close_mlj"]->Fill(selJetsLJKinematics[ij][0][0].mlj,evWgt);
	  histos_[ch+"_close_deta"]->Fill(selJetsLJKinematics[ij][0][0].deta,evWgt);
	  histos_[ch+"_close_dphi"]->Fill(selJetsLJKinematics[ij][0][0].dphi,evWgt);
	  histos_[ch+"_close_ptrel"]->Fill(selJetsLJKinematics[ij][0][0].ptrel,evWgt);
	  histos_[ch+"_far_mlj"]->Fill(selJetsLJKinematics[ij][0][1].mlj,evWgt);
	  histos_[ch+"_far_deta"]->Fill(selJetsLJKinematics[ij][0][1].deta,evWgt);
	  histos_[ch+"_far_dphi"]->Fill(selJetsLJKinematics[ij][0][1].dphi,evWgt);
	  histos_[ch+"_far_ptrel"]->Fill(selJetsLJKinematics[ij][0][1].ptrel,evWgt);
	  if(tmvaReader_) histos_[ch+"_kindisc"]->Fill(selJetsKINDisc[ij][0],evWgt);
	  histos_[ch+"_jp"]->Fill(ev.Jet_Proba[jetIdx],evWgt);
	  histos_[ch+"_svhe"]->Fill(ev.Jet_Svx[jetIdx],evWgt);
	  histos_[ch+"_csv"]->Fill(ev.Jet_CombIVF[jetIdx],evWgt);

	  Int_t flavBin(0),partonFlav(abs(ev.Jet_flavour[jetIdx]));
	  if(partonFlav==21 || (partonFlav>0 && partonFlav<4)) flavBin=1;
	  if(partonFlav==4) flavBin=2;
	  if(partonFlav==5) flavBin=3;
	  histos_[ch+"_flavour"]->Fill(flavBin,evWgt);
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
      weight_[5] = evWgt*puWgtLo/puWgtNom;
      weight_[6] = evWgt*puWgtHi/puWgtNom;
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
	  jetFlavour_[0] = ev.Jet_flavour[jetIdx];
	  jetPt_[0]      = selJetsP4[ij][0].Pt();
	  jetEta_[0]     = selJetsP4[ij][0].Eta();
	  for(size_t iSystVar=0; iSystVar<5; iSystVar++)
	    {
	      close_mlj_[iSystVar] = selJetsLJKinematics[ij][iSystVar][0].mlj;
	      if(tmvaReader_) kinDisc_[iSystVar]   = selJetsKINDisc[ij][iSystVar];
	      else            kinDisc_[iSystVar]=-999;
	    }
	  close_deta_=selJetsLJKinematics[ij][0][0].deta;
	  close_dphi_=selJetsLJKinematics[ij][0][0].dphi;
	  close_ptrel_=selJetsLJKinematics[ij][0][0].ptrel;
	  far_mlj_=selJetsLJKinematics[ij][0][0].mlj;
	  far_deta_=selJetsLJKinematics[ij][0][1].deta;
	  far_dphi_=selJetsLJKinematics[ij][0][1].dphi;
	  far_ptrel_=selJetsLJKinematics[ij][0][1].ptrel;
	  jp_[0]=ev.Jet_Proba[jetIdx];
	  svhe_[0]=ev.Jet_Svx[jetIdx];
	  csv_[0]=ev.Jet_CombIVF[jetIdx];

	  kinTree_->Fill();
	}

      //      ftmTree_->Fill();
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
  
  if(ch==-11*13) { res.first=0.91; res.second=0.05; }
  if(ch==-11*11) { res.first=0.95; res.second=0.05; }
  if(ch==-13*13) { res.first=0.92; res.second=0.05; }

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

//to be updated
std::vector<float> TTbarEventAnalysis::getJetEnergyScales(float pt,float eta,float rawsf,float area,float rho)
{
  std::vector<float> res(3,1.0);
  res[0]=1.0;
  res[1]=0.96;
  res[2]=1.04;
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
  for(std::map<TString,TH1F *>::iterator it = histos_.begin(); it != histos_.end(); it++) it->second->Write();
  kinTree_->Write();
  ftmTree_->Write();
  outF_->Close();
}

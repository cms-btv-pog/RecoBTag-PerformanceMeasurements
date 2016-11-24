#include "PtRelAnalyzer.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include "TFractionFitter.h"
 
#include "BTagAnalyzerUtilities.C"

PtRelAnalyzer::PtRelAnalyzer(TString TemplateFlag, TString PileUpReweighting, TString KinematicWeighting, TString SelectionFlag) {

  PUWeighting = PileUpReweighting;

  KinWeighting = KinematicWeighting;

  Selection = SelectionFlag;
  if (Selection=="") Selection = "_";
  Selection.ReplaceAll("_", "_" + AwayJetTreatment);

  TemplateVariable = TemplateFlag;
 
  if (TemplateVariable=="PtRel") {

    nBinsForTemp = 100;
    LowerEdgeForTemp = 0.;
    UpperEdgeForTemp = 5.;
    TemplateRebinning = 2;
    TemplateXTitle = "p_{T}^{rel} [GeV]";
    TemplateYTitle = "Jets / 0.1 GeV";
 
  } else if (TemplateVariable=="IP3D") {
        
    nBinsForTemp = 100.;
    LowerEdgeForTemp = -10.;
    UpperEdgeForTemp = 0.;
    TemplateRebinning = 2;
    TemplateXTitle = "log(| IP [cm] |)";
    TemplateYTitle = "Jets / 0.4";

  } else if (TemplateVariable=="IP3DSig") {

    nBinsForTemp = 120.;
    LowerEdgeForTemp = -6.;
    UpperEdgeForTemp = 6.;
    TemplateRebinning = 2;
    TemplateXTitle = "";
    TemplateYTitle = "";

  } else if (TemplateVariable=="SoftMu") {
    
    nBinsForTemp = 140.;
    LowerEdgeForTemp = -0.2;
    UpperEdgeForTemp = 1.2;
    TemplateRebinning = 2;
    TemplateXTitle = "";
    TemplateYTitle = "";

  } else if (TemplateVariable=="System8") {

    nBinsForTemp = 70;
    LowerEdgeForTemp = 0.;
    UpperEdgeForTemp = 7.;
    TemplateRebinning = 2;
    TemplateXTitle = "";
    TemplateYTitle = "";
    
  }
  
  TH1::SetDefaultSumw2();

}

PtRelAnalyzer::~PtRelAnalyzer() {

}

TString PtRelAnalyzer::HistogramName(TString VariableName, TString PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx) {

  TString ThisHistogramName = VariableName + "_" + PtBin + "_" + PtRelEtaBin[EtaBin];

  if (TriggerIdx>=0) ThisHistogramName += TriggerName[TriggerIdx];
  if (SystematicIdx>=0) ThisHistogramName += SystematicName[SystematicIdx];
  if (TaggerIdx>=0) { ThisHistogramName += "_"; ThisHistogramName += TaggerName[TaggerIdx]; }

  if (FlavourIdx>=0) {

    if (FlavourIdx<10) ThisHistogramName += "_Tag";
    else ThisHistogramName += "_Untag";
    
    if (FlavourIdx==0 || FlavourIdx==10)      ThisHistogramName += "";
    else if (FlavourIdx==1 || FlavourIdx==11) ThisHistogramName += "_b";
    else if (FlavourIdx==2 || FlavourIdx==12) ThisHistogramName += "_c";
    else if (FlavourIdx==3 || FlavourIdx==13) ThisHistogramName += "_lg";
    else if (FlavourIdx==4 || FlavourIdx==14) ThisHistogramName += "_trk";
    else std::cout << "PtRelAnalyzer::HistogramName: bad choice for paratemeter FlavourIdx " << FlavourIdx << std::endl;

    if (FlavourIdx==4 || FlavourIdx==14) {
      
      ThisHistogramName.ReplaceAll("_BTagMu", "_JetHT");
      ThisHistogramName.ReplaceAll( "_QCDMu", "_QCD");
      
    }
    
  }
  
  return ThisHistogramName;
  
}

TString PtRelAnalyzer::HistogramName(TString VariableName, int PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx) {

  return HistogramName(VariableName, PtRelPtBin[PtBin], EtaBin, TriggerIdx, SystematicIdx, TaggerIdx, FlavourIdx);

}

void PtRelAnalyzer::BookHistograms() {

  std::cout << "Booking histograms" << std::endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) {
     
      TString ThisHistoName;

      for (int is = 0; is<nSystematics; is++) {
	
	ThisHistoName = HistogramName("Observed_DataType", ptb, nb, -1, is, -1, -1);  
	Observed[ptb][nb][is] = new TH2D(ThisHistoName, ThisHistoName, nTriggers, 0., nTriggers, nTriggers, 0., nTriggers);

      }

      for (int tr = 0; tr<nTriggers; tr++)  {
	
	for (int is = 0; is<nSystematics; is++) {

	  ThisHistoName = HistogramName("muonPt_DataType", ptb, nb, tr, is, -1, -1);	  
	  MuPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	  
	  ThisHistoName = HistogramName("muonDR_DataType", ptb, nb, tr, is, -1, -1);
	  MuDRForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 50, 0., 0.5); 
	  
	  ThisHistoName = HistogramName("jetPt_DataType", ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, int(MaxPtRelPtEdge), 0., MaxPtRelPtEdge); 
	  
	  ThisHistoName = HistogramName("jetEta_DataType", ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.);
	  
	  ThisHistoName = HistogramName("nGoodPV_DataType", ptb, nb, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.); 
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	    	   
	    for (int fl = 0; fl<2; fl++) {
 
	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, fl);
	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 

	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 10 + fl);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	      
	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg,  2 + fl);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	    
	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 12 + fl);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	    	    
	    }
	  
	  }
	  
	}
	
      }

    }
    
}

void PtRelAnalyzer::ResetHistograms(TString DataType) {

  std::cout << "Resetting histograms" << std::endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int etab = 0; etab<nPtRelEtaBins; etab++) 
      for (int is = 0; is<nSystematics; is++) {

	TString ThisHistoName = HistogramName("Observed_" + DataType, ptb, etab, -1, is, -1, -1);
	Observed[ptb][etab][is]->SetNameTitle(ThisHistoName, ThisHistoName);
	Observed[ptb][etab][is]->Reset();
	
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  ThisHistoName = HistogramName("muonPt_" + DataType, ptb, etab, tr, is, -1, -1);
	  MuPtForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  MuPtForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("muonDR_" + DataType, ptb, etab, tr, is, -1, -1);
	  MuDRForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  MuDRForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetPtForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetEtaForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, etab, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  PVMultiplicity[tr][is][ptb][etab]->Reset();
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	  	    
	    for (int fl = 0; fl<2; fl++) {
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, fl);
	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset();
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 10 + fl);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset();
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  2 + fl);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset(); 
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 12 + fl);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset();
	      
	    }
	    
	  }
	  
	}
	
      }
  
}

void PtRelAnalyzer::NormalizeHistogramsToEntries() {

  std::cout << "Normalizing histograms to entries" << std::endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int etab = 0; etab<nPtRelEtaBins; etab++) 
      for (int is = 0; is<nSystematics; is++) {

	float HistoEntries, HistoIntegral;

	HistoEntries = Observed[ptb][etab][is]->GetEntries();
	HistoIntegral = Observed[ptb][etab][is]->Integral();
	Observed[ptb][etab][is]->Scale(HistoEntries/HistoIntegral);
	
	for (int tr = 0; tr<nTriggers; tr++) {

	  HistoEntries = MuPtForWeighting[tr][is][ptb][etab]->GetEntries();
	  HistoIntegral = MuPtForWeighting[tr][is][ptb][etab]->Integral(0,  MuPtForWeighting[tr][is][ptb][etab]->GetNbinsX()+1);

	  MuPtForWeighting[tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	  MuDRForWeighting[tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	  JetPtForWeighting[tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	  JetEtaForWeighting[tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	  PVMultiplicity[tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	  	    
	    for (int fl = 0; fl<2; fl++) {

	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Scale(HistoEntries/HistoIntegral);
	      
	    }
	    
	  }
	  
	}
	
      }
  
}

void PtRelAnalyzer::FillHistograms(TString DataType, int DataRange) {

  ResetHistograms(DataType);
  
  std::cout << "Filling histograms" << std::endl;

  GetKinematicWeights(DataType);
  GetPileUpWeights(DataType);
  GetBTemplateCorrections();
  GetTriggerPrescales(DataType);
  if (!DataType.Contains("QCD")) ReadJetEnergyUncertainty("_DATA", JECVersionData, "_AK4PFchs");
  else ReadJetEnergyUncertainty("_MC", JECVersionMC, "_AK4PFchs");

  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
  GetDataRangeInformation(DataType, DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
  
  int nEventsFromTrees = 0;
   
  for (int tf = FirstTree; tf<=nTrees; tf++) {

    TString FileName = GetFileName(FileDirectoryName, tf);
	        
    TFile *ThisTree = TFile::Open(FileName);
    
    if (!ThisTree) continue;

    TTree *tchain = GetChain(ThisTree, false);
    
    std::cout << "  Getting " << DataType << " sample: " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
    
    Int_t nentries = (Int_t)tchain->GetEntries();
	    
    nEventsFromTrees += nentries;
    
    for (Int_t i = 0; i<nentries; i++) {
      
      tchain->GetEntry(i);
   
      if (!DataType.Contains("QCD")) {
	if (PUWeighting.Contains("ICHEP2016") && Run>276811) continue;
      }

      int iMu; 
      int jMu = GetPFMuonJet(&iMu);
      
      if (jMu>=0 && nPV>0 && Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1] && Jet_looseID[jMu]==1) {
	
	int ptBin = -1;
	for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	  if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	
	int tr = -1;
	
	bool PassTrigger[nTriggers];
	
	for (int trg = 0; trg<nTriggers; trg++) {
	  
	  PassTrigger[trg] = false;
	  
	  if (Jet_pt[jMu]>MinPtJetTrigger[trg] && Jet_pt[jMu]<MaxPtJetTrigger[trg]) {
	    
	    bool FiredTrigger = false;
	    if (DataType=="BTagMu" || !Selection.Contains("TrgEmul")) {
	      if (PassTriggerBit(TriggerName[trg])) FiredTrigger = true;
	      if (FiredTrigger && (Selection.Contains("TrgConf") || Selection.Contains("TrgEmul"))) 
		if (!PassTriggerEmulation(trg, jMu)) FiredTrigger = false;
	    } else {
	      if (PassTriggerEmulation(trg, jMu)) FiredTrigger = true;
	    }
	    
	    if (FiredTrigger) {
	      
	      PassTrigger[trg] = true;
	      tr = trg;
	      
	    }
	    
	  }
	  
	}
	
	if (DataType=="QCDMu") 
	  if (!PassPtHat("MuEnrichedPt5", jMu)) tr = -1;	  
	if (DataType=="QCD") 
	  if (!PassPtHat("Inclusive", jMu)) tr = -1;
				     
	if (tr>=0) {
	  
	  bool JetAwayPass[nSystematics], EventJetAwayPass = false;
	  for (int is = 0; is<nSystematics; is++) {
	    
	    TString ThisAwayTaggerName = "JBPM4";//"TCHPM";
	    //if (Jet_pt[jMu]<50.) ThisAwayTaggerName = "JBPM7"; 
	    if (Jet_pt[jMu]>=140.) ThisAwayTaggerName = "JBPM5"; 
	    for (int at = 0; at<nAwayTaggers; at++) 
	      if (SystematicName[is].Contains(AwayTaggerName[at]))
		ThisAwayTaggerName = AwayTaggerName[at];
	    
	    int aJet = GetAwayJet(ThisAwayTaggerName, jMu, 1.5, true);
	    
	    JetAwayPass[is] = false;
	    
	    if (aJet>=0) {

	      float PtAwayJetCut = PtAwayJet[tr];

	      if (TriggerName[tr]=="_DiJet20" && Jet_pt[jMu]>=50.) PtAwayJetCut = 30.;

	      if (Jet_pt[aJet]>PtAwayJetCut) {

		JetAwayPass[is] = true;
		EventJetAwayPass = true;
		
	      }
	    
	    }

	  }
	  
	  if (EventJetAwayPass) {
		
	    int JetEtaBin = 0;
	    for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) 
	      if (fabs(Jet_eta[jMu])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[jMu])<PtRelEtaEdge[jeta]) JetEtaBin = jeta;
	    
	    int SavePtBin = ptBin;
		    
	    for (int is = 0; is<nSystematics; is++) {
	      
	      ptBin = SavePtBin;
	      
	      float TrackPtCut =  MuonPtCut[ptBin];
	      if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	      if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
	      
	      float TrackDRCut = 999., TrackMinDRCut = 0.;
	      GetDRCuts(SystematicName[is], Jet_pt[jMu], &TrackDRCut, &TrackMinDRCut);
	      
	      double JetMuonDR = DeltaR(Jet_eta[jMu], Jet_phi[jMu], PFMuon_eta[iMu], PFMuon_phi[iMu]);
	      
	      if (PFMuon_pt[iMu]>TrackPtCut && JetMuonDR<TrackDRCut  && JetMuonDR>=TrackMinDRCut && JetAwayPass[is]) {
		
		double ThisJetWeight = 1.;

		if (DataType!="BTagMu") {

		  int ipt = Jet_pt[jMu]; 
		  if (ipt>=PtRelPtEdge[nPtRelPtBins]) ipt = PtRelPtEdge[nPtRelPtBins] - 1;
		  int imu = PFMuon_pt[iMu]; if (imu>=60) imu = 59;
		  int ieta = JetEtaBin;
		  if (KinWeighting.Contains("KinEtaAfterPtBins") || KinWeighting.Contains("KinEtaBins"))
		    ieta = (Jet_eta[jMu]+2.4)/0.1;
		  
		  int EventPileUp = nPUtrue;
		  if (PUWeighting.Contains("PV") || PUWeighting.Contains("PSV")) EventPileUp = nPV;
		  if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
		  
		  ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][imu][ieta]*PileUpWeight[tr][EventPileUp][is];
		  
		} else {

		  ThisJetWeight = TriggerPrescaleWeight(DataType, tr);
		  
		}

		if (DataType!="BTagMu" && SystematicName[is].Contains("_GluonSplitting")) {
		  if (fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_GluonSplittingB")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 5)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==4 && SystematicName[is].Contains("_GluonSplittingC")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 4)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==21 && SystematicName[is].Contains("_GluonSplittingG")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 21)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } 
		}
		
		if (DataType!="BTagMu" && fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_BTemplates")) {
		  int iB = GetBHadron(jMu);
		  if (iB>=0) {
		    float EnergyFraction = BHadron_pT[iB]/Jet_genpt[jMu]; int efbin = EnergyFraction/0.02;
		    if (efbin>=0 && efbin<100) {
		      if (SystematicName[is].Contains("Minus")) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][0];
		      if (SystematicName[is].Contains("Plus") ) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][1];
		    }
		  }
		}
		/*
		if (DataType=="QCDMu" && SystematicName[is].Contains("_MuPtWeight")) {
		  
		  int ThisMuonPtBin = PFMuon_pt[iMu]; if (ThisMuonPtBin>=60) ThisMuonPtBin = 59;
		  ThisJetWeight *= MuPtWeight[ptBin][ThisMuonPtBin];
		  
		  }*/
			
		if (SystematicName[is].Contains("_JEU")) {
		  
		  int td = 0;
		  if (DataType=="QCDMu" || DataType=="QCD") td = 1;
		  float ThisJEU = GetJEU(td, Jet_eta[jMu], Jet_pt[jMu]);
		  if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		  if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		  
		  int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		  float ScaledJetPt = Jet_pt[jMu]*(1. + JEUSign*ThisJEU);
		  
		  ptBin = -1;
		  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		    if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		  
		}
		
		if (ptBin==-1) continue;
		
		for (int nb = 0; nb<nPtRelEtaBins; nb++) 
		  if (nb==0 || JetEtaBin==nb) {
		    
		    double Discriminator = PFMuon_ptrel[iMu];
		    		    
		    if (TemplateVariable=="IP3D")
		      Discriminator = log(fabs(PFMuon_IP[iMu]));
		    
		    if (TemplateVariable=="LinearIP3D") 
		      Discriminator = fabs(PFMuon_IP[iMu]);
		    
		    if (TemplateVariable=="IP3DSig")
		      Discriminator = log(fabs(PFMuon_IP[iMu]));
		    
		    if (TemplateVariable=="JetProb")	      
		      Discriminator = Jet_Proba[jMu];		    
		    
		    if (TemplateVariable=="SoftMu")
		      Discriminator = Jet_SoftMu[jMu];	
		    
		    if (Discriminator==0.) continue;

		    if (DataType=="BTagMu")
		      for (int trg1 = 0; trg1<nTriggers; trg1++) 
			for (int trg2 = 0; trg2<nTriggers; trg2++) 
			  if (PassTrigger[trg1] || PassTrigger[trg2]) 
			    Observed[ptBin][nb][is]->Fill(trg1, trg2, ThisJetWeight); 	     
		 
		    PVMultiplicity[tr][is][ptBin][nb]->Fill(nPV, ThisJetWeight);   
		    MuDRForWeighting[tr][is][ptBin][nb]->Fill(JetMuonDR, ThisJetWeight);    
		    MuPtForWeighting[tr][is][ptBin][nb]->Fill(PFMuon_pt[iMu], ThisJetWeight);
		    JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[jMu], ThisJetWeight);
		    JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[jMu], ThisJetWeight);
		    
		    for (int tg = 0; tg<nTaggers; tg++) 
		      if (IsTaggedJet(jMu, TaggerName[tg])) {
			
			if (DataType=="BTagMu") 
			  PtRelTagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else if (fabs(Jet_flavour[jMu])==5)	  
			  PtRelTagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
		        else if (fabs(Jet_flavour[jMu])==4) 
			  PtRelLightTagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else 	  
			  PtRelLightTagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
				
		      } else {
						
			if (DataType=="BTagMu")
			  PtRelUntagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else if (fabs(Jet_flavour[jMu])==5)
			  PtRelUntagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			
			else if (fabs(Jet_flavour[jMu])==4)
			  PtRelLightUntagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else 	  
			  PtRelLightUntagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			
		      }

		  } // Loop on eta bins

	      } // Muon pt and DR cuts
	      
	    } // Loop on systematics
	    
	  } // Has jet away
	  
	} // Pass trigger
	
      } // Muon-jet selection
      
    } // Loop on events
    
    ThisTree->Close();
	    
  } // Loop on files
  
  std::cout << "  nEventsFromTrees " << nEventsFromTrees << std::endl;
  
  TString Configuration = PUWeighting + KinWeighting; 
  if (!DataType.Contains("QCD")) {
    if (DataType=="BTagMu") Configuration = PUWeighting;
    Configuration.ReplaceAll("_PU", "_");
    if (DataType=="BTagMu") Configuration.ReplaceAll("_PV", "_");
    Configuration.ReplaceAll("_PSV", "_PS");
  }

  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    if ((DataType=="BTagMu" && nBTagMuRanges==1) || (DataType=="JetHT" && nJetRunRanges==1))
      NormalizeHistogramsToEntries();

  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + "_" + DataRangeName + Configuration + Selection + ".root";

  SaveHistograms(OutputFileName, DataType);

}

void PtRelAnalyzer::FillSubjetHistograms(TString DataType, int DataRange) {

  ResetHistograms(DataType);
  
  std::cout << "Filling histograms" << std::endl;

  TString DataTypeForKinReweighting = DataType;
  //if (DataType=="QCD") DataTypeForKinReweighting = "QCDMu";
  GetKinematicWeights(DataTypeForKinReweighting);

  TString DataTypeForPUReweighting = DataType;
  if (DataType=="QCD") DataTypeForPUReweighting += "BTagMu";
  GetPileUpWeights(DataTypeForPUReweighting);
  
  GetBTemplateCorrections();
  if (!DataType.Contains("QCD")) ReadJetEnergyUncertainty("_DATA", JECVersionData, "_AK8PFchs");
  else ReadJetEnergyUncertainty("_MC", JECVersionMC, "_AK8PFchs");
  
  //int ph = DataRange;

  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
  GetDataRangeInformation(DataType, DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
  
  int nEventsFromTrees = 0;
   
  for (int tf = FirstTree; tf<=nTrees; tf++) {
    
    TString FileName = GetFileName(FileDirectoryName, tf);
	        
    TFile *ThisTree = TFile::Open(FileName);
    
    if (!ThisTree) continue;

    TTree *tchain = GetChain(ThisTree, false, "btaganaFatJets");
    TTree *tstdchain = GetChain(ThisTree, false, "btagana", true);
    
    if (DataType=="BTagMu") std::cout << "  Getting " << DataType << "   sample: BTagMu    " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
    if (DataType=="JetHT") std::cout << "  Getting " << DataType << "   sample: JetHT    " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
    if (DataType=="QCDMu") std::cout << "  Getting " << DataType << "    sample: MuPt5 " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
    
    Int_t nentries = (Int_t)tchain->GetEntries();
	    
    nEventsFromTrees += nentries;
    
    for (Int_t i = 0; i<nentries; i++) {
      
      tchain->GetEntry(i);
      tstdchain->GetEntry(i);
      
      if (!DataType.Contains("QCD")) {
	if (PUWeighting.Contains("ICHEP2016") && Run>276811) continue;
      }
      
      int iMu; 
      int jMu = GetPFMuonJet(&iMu, false);
      
      if (jMu>=0 && nPV>0 && Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1]) {

	int ajet = -1;
	int fatjet = GetFatJet(jMu, &ajet);

	if (fatjet<0 || (ajet<0 && !Selection.Contains("StdAway"))) continue;

	int ptBin = -1;
	for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	  if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	
	int tr = -1;
	
	bool PassTrigger[nTriggers];
	
	for (int trg = 0; trg<nTriggers; trg++) {
	  
	  PassTrigger[trg] = false;
	 
	  if (FatJet_pt[fatjet]>MinPtJetTrigger[trg] && FatJet_pt[fatjet]<MaxPtJetTrigger[trg]) {
	  
	    bool FiredTrigger = false;
	    if (DataType=="BTagMu" || DataType=="JetHT" || !Selection.Contains("TrgEmul")) {
	      if (PassTriggerBit(TriggerName[trg])) FiredTrigger = true;
	      if (FiredTrigger && (Selection.Contains("TrgConf") || Selection.Contains("TrgEmul"))) 
		if (!PassTriggerEmulation(trg, jMu)) FiredTrigger = false;
	    } else {
	      if (PassTriggerEmulation(trg, jMu)) FiredTrigger = true;
	    }
	    
	    if (FiredTrigger) {
	      
	      PassTrigger[trg] = true;
	      tr = trg;
	      
	    }
	    
	  }
	  
	}
	
	if (DataType=="QCDMu") 
	  if (!PassPtHat("MuEnrichedPt5", jMu)) tr = -1;	  
	if (DataType=="QCD") 
	  if (!PassPtHat("Inclusive", jMu)) tr = -1;
	
	if (tr>=0) {

	  bool JetAwayPass[nSystematics], EventJetAwayPass = false;
	  for (int is = 0; is<nSystematics; is++) {

	    JetAwayPass[is] = false;
	    
	    TString ThisAwayTaggerName = "JBPL";
	    for (int at = 0; at<nAwayTaggers; at++) 
	      if (SystematicName[is].Contains(AwayTaggerName[at]))
		ThisAwayTaggerName = AwayTaggerName[at];
	    
	    int aJet = -1;
	    if (Selection.Contains("StdAway"))
	      aJet = GetAwayStdJet(ThisAwayTaggerName, jMu, 1.5, true, tr);
	    else
	      aJet = IsTaggedJet(ajet, ThisAwayTaggerName) - 1;
	    
	    if (aJet>=0 || ThisAwayTaggerName=="NONE") {
	      
	      JetAwayPass[is] = true;
	      EventJetAwayPass = true;
	      
	    }
	    
	  }
	  
	  if (EventJetAwayPass) {

	    int JetEtaBin = 0;
	    for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) 
	      if (fabs(Jet_eta[jMu])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[jMu])<PtRelEtaEdge[jeta]) JetEtaBin = jeta;
	    
	    int SavePtBin = ptBin;
	    
	    for (int is = 0; is<nSystematics; is++) {
	      
	      ptBin = SavePtBin;
	      
	      float TrackPtCut =  MuonPtCut[ptBin];
	      if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	      if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
	      
	      float TrackDRCut = 999., TrackMinDRCut = 0.;
	      GetDRCuts(SystematicName[is], Jet_pt[jMu], &TrackDRCut, &TrackMinDRCut);
	      
	      double JetMuonDR = DeltaR(Jet_eta[jMu], Jet_phi[jMu], PFMuon_eta[iMu], PFMuon_phi[iMu]);
	      
	      if (PFMuon_pt[iMu]>TrackPtCut && JetMuonDR<TrackDRCut  && JetMuonDR>=TrackMinDRCut && JetAwayPass[is]) {
		
		double ThisJetWeight = 1.;

		if (DataType=="QCDMu" || DataType=="QCD") {

		  int ipt = Jet_pt[jMu]; 
		  if (ipt>=PtRelPtEdge[nPtRelPtBins]) ipt = PtRelPtEdge[nPtRelPtBins] - 1;
		  int imu = PFMuon_pt[iMu]; if (imu>=60) imu = 59;
		  int ieta = JetEtaBin;
		  if (KinWeighting.Contains("KinEtaAfterPtBins") || KinWeighting.Contains("KinEtaBins"))
		    ieta = (Jet_eta[jMu]+2.4)/0.1;

		  int EventPileUp = nPUtrue;
		  if (PUWeighting.Contains("PV") || PUWeighting.Contains("PSV")) EventPileUp = nPV;
		  if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
		  
		  ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][imu][ieta]*PileUpWeight[tr][EventPileUp][is];
		  		  
		}

		if (DataType!="BTagMu" && DataType!="JetHT" && SystematicName[is].Contains("_GluonSplitting")) {
		  if (fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_GluonSplittingB")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 5)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==4 && SystematicName[is].Contains("_GluonSplittingC")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 4)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==21 && SystematicName[is].Contains("_GluonSplittingG")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 21)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } 
		}
		
		if (DataType!="BTagMu" && DataType!="JetHT" && fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_BTemplates")) {
		  int iB = GetBHadron(jMu);
		  if (iB>=0) {
		    float EnergyFraction = BHadron_pT[iB]/Jet_genpt[jMu]; int efbin = EnergyFraction/0.02;
		    if (efbin>=0 && efbin<100) {
		      if (SystematicName[is].Contains("Minus")) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][0];
		      if (SystematicName[is].Contains("Plus") ) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][1];
		    }
		  }
		}
		/*
		if (DataType=="QCDMu" && SystematicName[is].Contains("_MuPtWeight")) {
		  
		  int ThisMuonPtBin = PFMuon_pt[iMu]; if (ThisMuonPtBin>=60) ThisMuonPtBin = 59;
		  ThisJetWeight *= MuPtWeight[ptBin][ThisMuonPtBin];
		  
		  }*/
			
		if (SystematicName[is].Contains("_JEU")) {
		  
		  int td = 0;
		  if (DataType=="QCDMu" || DataType=="QCD") td = 1;
		  float ThisJEU = GetJEU(td, Jet_eta[jMu], Jet_pt[jMu]);
		  if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		  if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		  
		  int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		  float ScaledJetPt = Jet_pt[jMu]*(1. + JEUSign*ThisJEU);
		  
		  ptBin = -1;
		  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		    if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		  
		}
		
		if (ptBin==-1) continue;
		
		for (int nb = 0; nb<nPtRelEtaBins; nb++) 
		  if (nb==0 || nb==JetEtaBin) {
		    
		    double Discriminator = PFMuon_ptrel[iMu];
		    		    
		    if (TemplateVariable=="IP3D")
		      Discriminator = log(fabs(PFMuon_IP[iMu]));
		    
		    if (TemplateVariable=="LinearIP3D") 
		      Discriminator = fabs(PFMuon_IP[iMu]);
		    
		    if (TemplateVariable=="IP3DSig")
		      Discriminator = log(fabs(PFMuon_IP[iMu]));
		    
		    if (TemplateVariable=="JetProb")	      
		      Discriminator = Jet_Proba[jMu];		    
		    
		    if (TemplateVariable=="SoftMu")
		      Discriminator = Jet_SoftMu[jMu];	
		    
		    if (Discriminator==0.) continue;

		    if (DataType=="BTagMu" || DataType=="JetHT")
		      for (int trg1 = 0; trg1<nTriggers; trg1++) 
			for (int trg2 = 0; trg2<nTriggers; trg2++) 
			  if (PassTrigger[trg1] || PassTrigger[trg2]) 
			    Observed[ptBin][nb][is]->Fill(trg1, trg2, ThisJetWeight); 	     
		 
		    PVMultiplicity[tr][is][ptBin][nb]->Fill(nPV, ThisJetWeight);   
		    MuDRForWeighting[tr][is][ptBin][nb]->Fill(JetMuonDR, ThisJetWeight);    
		    MuPtForWeighting[tr][is][ptBin][nb]->Fill(PFMuon_pt[iMu], ThisJetWeight);
		    JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[jMu], ThisJetWeight);
		    JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[jMu], ThisJetWeight);
		 
		    for (int tg = 0; tg<nTaggers; tg++) 
		      if (IsTaggedJet(jMu, TaggerName[tg])) {
			
			if (DataType=="BTagMu" || DataType=="JetHT") 
			  PtRelTagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else if (fabs(Jet_flavour[jMu])==5)	  
			  PtRelTagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
		        else if (fabs(Jet_flavour[jMu])==4) 
			  PtRelLightTagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else 	  
			  PtRelLightTagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
				
		      } else {
						
			if (DataType=="BTagMu" || DataType=="JetHT")
			  PtRelUntagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else if (fabs(Jet_flavour[jMu])==5)
			  PtRelUntagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			
			else if (fabs(Jet_flavour[jMu])==4)
			  PtRelLightUntagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else 	  
			  PtRelLightUntagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			
		      }

		  } // Loop on eta bins

	      } // Muon pt and DR cuts
	      
	    } // Loop on systematics
	    
	  } // Has jet away
	  
	} // Pass trigger
	
      } // Muon-jet selection
      
    } // Loop on events
    
    ThisTree->Close();
	    
  } // Loop on files
  
  std::cout << "  nEventsFromTrees " << nEventsFromTrees << std::endl;

  TString Configuration = PUWeighting + KinWeighting; 
  if (!DataType.Contains("QCD")) {
    if (DataType=="BTagMu") Configuration = PUWeighting;
    Configuration.ReplaceAll("_PU", "_");
    if (DataType=="BTagMu") Configuration.ReplaceAll("_PV", "_");
    Configuration.ReplaceAll("_PSV", "_PS");
  }
  
  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    if ((DataType=="BTagMu" && nBTagMuRanges==1) || (DataType=="JetHT" && nJetRunRanges==1))
      NormalizeHistogramsToEntries();

  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + "_" + DataRangeName + Configuration + Selection + ".root";
  
  SaveHistograms(OutputFileName, DataType);

}

void PtRelAnalyzer::SaveHistograms(TString OutputFileName, TString DataType) {
  
  std::cout << "Saving histograms" << std::endl;
  
  TFile *OutFile = new TFile(OutputFileName, "recreate");
  
  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int etab = 0; etab<nPtRelEtaBins; etab++) 
      for (int is = 0; is<nSystematics; is++) {

	if (!DataType.Contains("QCD")) {
	  
	  Observed[ptb][etab][is]->Write();
	  
	}

	for (int tr = 0; tr<nTriggers; tr++) {

	  MuDRForWeighting[tr][is][ptb][etab]->Write();
	  MuPtForWeighting[tr][is][ptb][etab]->Write();
	  JetPtForWeighting[tr][is][ptb][etab]->Write();
	  JetEtaForWeighting[tr][is][ptb][etab]->Write();
	  PVMultiplicity[tr][is][ptb][etab]->Write();
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	   
	    if (!DataType.Contains("QCD")) {
	      
	      PtRelTagForWeighting[tg][tr][is][ptb][etab]->Write();
	      PtRelUntagForWeighting[tg][tr][is][ptb][etab]->Write();
	    
	    } else {
	      
	      PtRelTagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Write();
	      PtRelUntagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Write();

	      for (int fl = 0; fl<2; fl++) {
		  
		PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Write(); 
		PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Write();
		  
	      }
	      
	    }
	    
	  }
	  
	}

      }

  OutFile->Close();

  std::cout << "Exiting" << std::endl;
  
}
 
void PtRelAnalyzer::MergeHistograms(TString DataType, int FirstBin, int LastBin) {

  ResetHistograms(DataType);

  std::cout << "Merging histograms: " << DataType << std::endl;

  if (LastBin==100) LastBin = GetNumberOfDataRanges(DataType);
  
  for (int ph = FirstBin; ph<LastBin; ph++) {

    TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
    GetDataRangeInformation(DataType, ph, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
    
    std::cout << "  Merging " << DataRangeName << std::endl;
  
    TString InputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + "_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";

    if (DataType=="BTagMu") {
      InputFileName.ReplaceAll("_PV", "_");
      InputFileName.ReplaceAll("_PU", "_");
      InputFileName.ReplaceAll("_PSV", "_PS");
      InputFileName.ReplaceAll(KinWeighting, "");
    }
    TFile *HistogramFile = TFile::Open(InputFileName); 

    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

      std::cout << "    Merging " << PtRelPtBin[ptb] << std::endl;

      for (int etab = 0; etab<nPtRelEtaBins; etab++) 
	for (int is = 0; is<nSystematics; is++) {
	  
	  TString ThisHistoName;

	  if (DataType=="BTagMu" || DataType=="JetHT") {
	    
	    ThisHistoName = HistogramName("Observed_" + DataType, ptb, etab, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) HistogramFile->Get(ThisHistoName);
	    Observed[ptb][etab][is]->Add(ThisObserved);
	  
	  }

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    ThisHistoName = HistogramName("muonPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisMuPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    MuPtForWeighting[tr][is][ptb][etab]->Add(ThisMuPt);
	    
	    ThisHistoName = HistogramName("muonDR_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisMuDR = (TH1D*) HistogramFile->Get(ThisHistoName);
	    MuDRForWeighting[tr][is][ptb][etab]->Add(ThisMuDR);
	    
	    ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetPtForWeighting[tr][is][ptb][etab]->Add(ThisJetPt);
	    
	    ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetEta = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetEtaForWeighting[tr][is][ptb][etab]->Add(ThisJetEta);
	    
	    ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisPVMult = (TH1D*) HistogramFile->Get(ThisHistoName);
	    PVMultiplicity[tr][is][ptb][etab]->Add(ThisPVMult);
	    
	    for (int tg = 0; tg<nTaggers; tg++) {
	      
	      if (DataType=="BTagMu" || DataType=="JetHT") {
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  0);
		TH1D *ThisTag = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelTagForWeighting[tg][tr][is][ptb][etab]->Add(ThisTag);
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 10);
		TH1D *ThisUntag = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelUntagForWeighting[tg][tr][is][ptb][etab]->Add(ThisUntag);
		
	      } else {
		
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  1);
		TH1D *ThisTagB = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelTagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Add(ThisTagB);
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 11);
		TH1D *ThisUntagB = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelUntagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Add(ThisUntagB);
		
		for (int fl = 0; fl<2; fl++) {
		  
		
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  2 + fl);
		  TH1D *ThisTagL = (TH1D*) HistogramFile->Get(ThisHistoName);
		  PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Add(ThisTagL);
		
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 12 + fl);
		  TH1D *ThisUntagL = (TH1D*) HistogramFile->Get(ThisHistoName); 
		  PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Add(ThisUntagL);
		  
		}
		
	      }
	      
	    }
	    
	  }
	  
	}

    }

  }

  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    NormalizeHistogramsToEntries();
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + PUWeighting + KinWeighting + Selection + ".root";

  if (DataType=="BTagMu") {
    OutputFileName.ReplaceAll("_PV", "_");
    OutputFileName.ReplaceAll("_PU", "_");
    OutputFileName.ReplaceAll("_PSV", "_PS");
    OutputFileName.ReplaceAll(KinWeighting, "");
  }
  
  SaveHistograms(OutputFileName, DataType);
  
}

void PtRelAnalyzer::BookLightHistograms() {
  
  std::cout << "Booking light histograms" << std::endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int tr = 0; tr<nTriggers; tr++) 
	for (int is = 0; is<nSystematics; is++)	{
	  
	  TString ThisHistoName = HistogramName("jetPt_DataType", ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, int(MaxPtRelPtEdge), 0., MaxPtRelPtEdge); 
	  
	  ThisHistoName = HistogramName("jetEta_DataType", ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.);
  
	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 4);
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	    
	    ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 14);
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
 
	    bool DoErrors = false;
  
	    if (DoErrors) {
	      
	      PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Sumw2();
	      PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Sumw2();
	      
	    }
	    
	  } 

        } 
  
}

void PtRelAnalyzer::ResetLightHistograms(TString DataType) {

  std::cout << "Resetting light histograms" << std::endl; 

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
      for (int is = 0; is<nSystematics; is++)	
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  TString ThisHistoName = HistogramName("jetPt_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetPtForWeighting[tr][is][ptb][nb]->Reset();
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetEtaForWeighting[tr][is][ptb][nb]->Reset();
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, nb, tr, is, tg, 4);
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);;
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Reset();
	    
	    ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, nb, tr, is, tg, 14);
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);;
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Reset();
	    
	  }

        } 
  
}

void PtRelAnalyzer::NormalizeLightHistogramsToEntries() {

  std::cout << "Normalizing light histograms to entries" << std::endl; 

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
      for (int is = 0; is<nSystematics; is++)	
	for (int tr = 0; tr<nTriggers; tr++) {

	  float HistoEntries, HistoIntegral;
	  
	  HistoEntries = JetPtForWeighting[tr][is][ptb][nb]->GetEntries();
	  HistoIntegral = JetPtForWeighting[tr][is][ptb][nb]->Integral(0, JetPtForWeighting[tr][is][ptb][nb]->GetNbinsX()+1);

	  JetPtForWeighting[tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	  JetEtaForWeighting[tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	  	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	  
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	    
	  }

        } 
  
}

void PtRelAnalyzer::FillLightHistograms(TString DataType, int DataRange) {

  ResetLightHistograms(DataType);

  std::cout << "Filling light histograms" << std::endl;
  
  GetKinematicWeights(DataType);
  GetPileUpWeights(DataType);
  GetTriggerPrescales(DataType);
  if (!DataType.Contains("QCD")) ReadJetEnergyUncertainty("_DATA", JECVersionData, "_AK4PFchs");
  else ReadJetEnergyUncertainty("_MC", JECVersionMC, "_AK4PFchs");

  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
  GetDataRangeInformation(DataType, DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);

  int nEventsFromTrees = 0;
 
  for (int tf = FirstTree; tf<=nTrees; tf++) {
	
    TString FileName = GetFileName(FileDirectoryName, tf);
    	 
    TFile *ThisTree = TFile::Open(FileName);
	
    if (!ThisTree) continue;

    TString DataRangeName;
    if (DataType=="JetHT")DataRangeName = JetRunRangeName[DataRange];
    if (DataType=="QCD")DataRangeName = MCInclusivePtHatRange[DataRange];

    TTree *tchain = GetChain(ThisTree, true);
    
    std::cout << "  Getting " << DataType << " sample: " << DataRangeName << " " << tf << "/" << nTrees << std::endl;
	
    Int_t nentries = (Int_t)tchain->GetEntries();

    nEventsFromTrees += nentries;

    for (Int_t i = 0; i<nentries; i++) {

      tchain->GetEntry(i);
      
      if (!DataType.Contains("QCD")) {
	if (PUWeighting.Contains("ICHEP2016") && Run>276811) continue;
      }

      for (int ijet = 0; ijet<nJet; ijet++) {
	    
	if (Jet_pt[ijet]>PtRelPtEdge[0] && !HasPFMuon(ijet, false) && !HasTaggedJet(ijet) && nPV>0 && 
	    Jet_pt[ijet]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[ijet])<PtRelEtaEdge[nPtRelEtaBins-1]
	    && Jet_looseID[ijet]==1)  {
	      
	  int ptBin = -1;
	  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
            if (Jet_pt[ijet]>PtRelPtEdge[ptb]) ptBin = ptb;
	   
	  int IdxAwayJet = GetAwayJet("NONE", ijet, 1.5, false);
	      
          if (IdxAwayJet<0) continue;
	      
	  int tr = -1;

	  for (int trg = 0; trg<nTriggers; trg++)
            if (Jet_pt[IdxAwayJet]>PtAwayJet[trg] && Jet_pt[ijet]>MinPtJetTrigger[trg] && Jet_pt[ijet]<MaxPtJetTrigger[trg]) {

	      if (TriggerName[trg]=="_DiJet20" && Jet_pt[ijet]>=50.) 
		if (Jet_pt[IdxAwayJet]<30.) continue; 
		  
       	      bool PassTrigger = PassTriggerBit(JetTriggerName[trg]);
		  
	      if (DataType=="QCD")// && Selection.Contains("TrgEmul"))
		PassTrigger = true;

	      if (!PassEventTriggerEmulation(trg, DataType)) PassTrigger = false;
	      
	      if (trg==1 && MaxPtJetTrigger[trg-1]>MinPtJetTrigger[trg])
		if (Jet_pt[ijet]<MaxPtJetTrigger[0] && i%2!=1) PassTrigger = false;

	      if (PassTrigger && (Selection.Contains("TrgConf") || Selection.Contains("TrgEmul"))) 
		if (!PassTriggerEmulation(trg, ijet)) PassTrigger = false;
		
	      if (PassTrigger) tr = trg;
		  
	   }
		  
	   if (DataType=="QCD") 
             if (!PassPtHat("Inclusive", ijet)) tr = -1;
	      
	   if (tr>=0) {

	     int JetEtaBin = 0;
	     for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) 
	       if (fabs(Jet_eta[ijet])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[ijet])<PtRelEtaEdge[jeta]) JetEtaBin = jeta;
	     
	     int SavePtBin = ptBin;
		
	     for (int is = 0; is<nSystematics; is++) {
		  
	       ptBin = SavePtBin;
		  
               double ThisJetWeight = 1.;
	       int ipt = Jet_pt[ijet]; if (ipt>=PtRelPtEdge[nPtRelPtBins]) ipt = PtRelPtEdge[nPtRelPtBins] - 1;
	       int ieta = JetEtaBin;
	       if (KinWeighting.Contains("KinEtaAfterPtBins") || KinWeighting.Contains("KinEtaBins"))
		 ieta = (Jet_eta[ijet]+2.4)/0.1;
	       int EventPileUp = (DataType.Contains("QCD")) ? nPUtrue : 1;
               if (PUWeighting.Contains("PV") || PUWeighting.Contains("PSV")) EventPileUp = nPV;
	       if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
	       ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][0][ieta]*PileUpWeight[tr][EventPileUp][is];
	       ThisJetWeight *= TriggerPrescaleWeight(DataType, tr);

	       float TrackPtCut = MuonPtCut[ptBin];
	       if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	       if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
	       
	       if (SystematicName[is].Contains("_JEU")) {
		 
		 int td = 0;
		 if (DataType=="QCDMu" || DataType=="QCD") td = 1;
		 float ThisJEU = GetJEU(td, Jet_eta[ijet], Jet_pt[ijet]);
		    if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		    if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		    
		    int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		    float ScaledJetPt = Jet_pt[ijet]*(1. + JEUSign*ThisJEU);
		    
		    ptBin = -1;
		    for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		      if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		    
	       }
	       
	       if (ptBin==-1) continue;
		  
	       float TrackDRCut = 999., TrackMinDRCut = 0.;
	       GetDRCuts(SystematicName[is], Jet_pt[ijet], &TrackDRCut, &TrackMinDRCut);
		  
	       bool FilledJet[nPtRelEtaBins]; 
	       for (int jeta = 0; jeta<nPtRelEtaBins; jeta++) 
		 FilledJet[jeta] = false; 

	       int nLightTracks = 0;//Jet_nLastTrkInc[ijet] - Jet_nFirstTrkInc[ijet];
		  
	       for (int tk = Jet_nFirstTrkInc[ijet]; tk<Jet_nLastTrkInc[ijet]; tk++) {
		    
		 double JetTrkDR = DeltaR(Jet_eta[ijet], Jet_phi[ijet], TrkInc_eta[tk], TrkInc_phi[tk]);
		    
		 if (TrkInc_pt[tk]>TrackPtCut && JetTrkDR<TrackDRCut && JetTrkDR>=TrackMinDRCut) nLightTracks++;
		 //if (fabs(Jet_eta[ijet]-TrkInc_eta[tk])<TrackDRCut && DeltaPhi(Jet_phi[ijet], TrkInc_phi[tk])<TrackDRCut) nLightTracks++;
		    
	       }
	     
	       for (int tk = Jet_nFirstTrkInc[ijet]; tk<Jet_nLastTrkInc[ijet]; tk++) {
		    
		 double JetTrkDR = DeltaR(Jet_eta[ijet], Jet_phi[ijet], TrkInc_eta[tk], TrkInc_phi[tk]);
		    
		 if (TrkInc_pt[tk]>TrackPtCut && JetTrkDR<TrackDRCut && JetTrkDR>=TrackMinDRCut) {
		 //if (TrkInc_pt[tk]>TrackPtCut && fabs(Jet_eta[ijet]-TrkInc_eta[tk])<TrackDRCut && DeltaPhi(Jet_phi[ijet], TrkInc_phi[tk])<TrackDRCut) {
		      
		 double ThisTrackWeight = ThisJetWeight;

		 /*if (DataType=="QCD" && ApplyDeltaRCorrections) {
		 int idrpt = 0, idreta = 0, idrdr = 0;
			float drjetpt = Jet_pt[ijet], drjeteta = fabs(Jet_eta[ijet]);
			if (drjetpt<30) idrpt = 0;
			else if (drjetpt<50) idrpt = 1;
			else if (drjetpt<80) idrpt = 2;
			else if (drjetpt<120) idrpt = 3;
			else if (drjetpt<160) idrpt = 4;
			else if (drjetpt<320) idrpt = 5;
			else if (drjetpt<500) idrpt = 6;
			else if (drjetpt<800) idrpt = 7;
			if (drjeteta<0.4) idreta = 0;
			else if (drjeteta<0.8) idreta = 1;
			else if (drjeteta<1.2) idreta = 2;
			else if (drjeteta<1.8) idreta = 3;
			else if (drjeteta<12.4) idreta = 4;
			idrdr = JetTrkDR/0.02;
			int iDEta = fabs(Jet_eta[ijet]-TrkInc_eta[tk])/0.01;
			//if (iDEta>=0 && iDEta<50)
			//ThisTrackWeight *= DeltaRCorrections[iDEta][0][idrpt][idreta];
			int ThisDPhi = fabs(Jet_phi[ijet]-TrkInc_phi[tk]);
			if (ThisDPhi>=3.14159) ThisDPhi = 2*3.14159 - ThisDPhi;
			int iDPhi = ThisDPhi/0.01;
			//if (iDPhi>=0 && iDPhi<50)
			//ThisTrackWeight *= DeltaRCorrections[iDPhi][0][idrpt][idreta];
			if (idrdr<50)
			  ThisTrackWeight *= DeltaRCorrections[idrdr][0][idrpt][idreta];
		      }*/
		      
		 double Discriminator = TrkInc_ptrel[tk];
		      
		 if (TemplateVariable=="IP2D")
		   Discriminator = log(fabs(TrkInc_IP[tk]));
		      
		 if (TemplateVariable=="IP2DSig")
		   Discriminator = log(fabs(TrkInc_IPsig[tk]));
		      
		 if (TemplateVariable=="IP3D")
	           Discriminator = log(fabs(TrkInc_IP[tk]));
		      
		 if (TemplateVariable=="LinearIP3D")
	           Discriminator = TrkInc_IP[tk];
		      
		 if (TemplateVariable=="IP3DSig")
	           Discriminator = log(fabs(TrkInc_IPsig[tk]));
		      		      
		 for (int nb = 0; nb<nPtRelEtaBins; nb++) 
	           if (nb==0 || JetEtaBin==nb) {

                     if (!FilledJet[nb]) {			  
		       JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[ijet], ThisJetWeight);
                       JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[ijet], ThisJetWeight);
                       FilledJet[nb] = true;
		
                     }

	             PtRelLightUntagForSystematic[nTaggers-1][tr][is][ptBin][nb]->Fill(Discriminator, ThisTrackWeight/nLightTracks);
			  
	           } // Loop on eta bins
		      
		 } // Muon pt and DR cuts
			
	       } // Loop on tracks
		      
	     } // Loop on systematics
		
	   } // Pass trigger
	    
	 } // Light jet selection
	    
       } // Loop on jets
	  
     } // Loop on events
	
     ThisTree->Close();
	
   } // Loop on files
      
   std::cout << "  nEventsFromTrees " << nEventsFromTrees << std::endl;
    
   for (int tr = 0; tr<nTriggers; tr++)
     for (int ptb = 0; ptb<nPtRelPtBins; ptb++) 
       for (int nb = 0; nb<nPtRelEtaBins; nb++)
	 for (int is = 0; is<nSystematics; is++) 
	   for (int tg = 0; tg<nTaggers; tg++) {

             PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Add(PtRelLightUntagForSystematic[nTaggers-1][tr][is][ptb][nb]);

             if (tg<nTaggers-1)
	       PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Add(PtRelLightUntagForSystematic[nTaggers-1][tr][is][ptb][nb]);

	   }

  std::cout << "  nEventsFromTrees " << nEventsFromTrees << std::endl; 
  
  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    if ((DataType=="BTagMu" && nBTagMuRanges==1) || (DataType=="JetHT" && nJetRunRanges==1))
      NormalizeLightHistogramsToEntries();

  TString Configuration = PUWeighting + KinWeighting; 
  if (!DataType.Contains("QCD")) {
    if (DataType=="BTagMu") Configuration = PUWeighting;
    Configuration.ReplaceAll("_PU", "_");
    if (DataType=="BTagMu") Configuration.ReplaceAll("_PV", "_");
    Configuration.ReplaceAll("_PSV", "_PS");
  }

  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_" + DataRangeName + Configuration + Selection + ".root";
  
  SaveLightHistograms(OutputFileName);
  
}

void PtRelAnalyzer::SaveLightHistograms(TString TemplateFileName) {

  std::cout << "Saving light histograms" << std::endl;

  TFile *OutFile = new TFile(TemplateFileName, "recreate");

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int tr = 0; tr<nTriggers; tr++) 
	for (int is = 0; is<nSystematics; is++)	{

	  JetPtForWeighting[tr][is][ptb][nb]->Write();
	  JetEtaForWeighting[tr][is][ptb][nb]->Write();

	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Write();
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Write();
	    
	  }

        }
  
  OutFile->Close();

  std::cout << "Exiting" << std::endl;

}

void PtRelAnalyzer::MergeLightHistograms(TString DataType, int FirstBin, int LastBin) {

  ResetLightHistograms(DataType);

  std::cout << "Merging light histograms: " << DataType << std::endl;

  if (LastBin==100) LastBin = GetNumberOfDataRanges(DataType);

  TString Configuration = PUWeighting + KinWeighting; 
  if (!DataType.Contains("QCD")) {
    if (DataType=="BTagMu") Configuration = PUWeighting;
    Configuration.ReplaceAll("_PU", "_");
    if (DataType=="BTagMu") Configuration.ReplaceAll("_PV", "_");
    Configuration.ReplaceAll("_PSV", "_PS");
  }
  
  for (int ph = FirstBin; ph<LastBin; ph++) {
    
    TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
    GetDataRangeInformation(DataType, ph, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);

    std::cout << "  Merging " << DataRangeName << std::endl;
  
    TString InputFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_" + DataRangeName + Configuration + Selection + ".root";

    TFile *HistogramFile = TFile::Open(InputFileName); 

    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

      std::cout << "    Merging " << PtRelPtBin[ptb] << std::endl;

      for (int etab = 0; etab<nPtRelEtaBins; etab++) 
	for (int is = 0; is<nSystematics; is++) {
	  
	  TString ThisHistoName;

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetPt = (TH1D*) HistogramFile->Get(ThisHistoName); 
	    if (ThisJetPt->GetEntries()==0) continue;
	    JetPtForWeighting[tr][is][ptb][etab]->Add(ThisJetPt);
	    
	    ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetEta = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetEtaForWeighting[tr][is][ptb][etab]->Add(ThisJetEta);
	    
	    for (int tg = 0; tg<nTaggers; tg++) {

	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 4);
              TH1D *ThisTag = (TH1D*) HistogramFile->Get(ThisHistoName);
	      PtRelLightTagForSystematic[tg][tr][is][ptb][etab]->Add(ThisTag);

	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 14);
              TH1D *ThisUntag = (TH1D*) HistogramFile->Get(ThisHistoName);
	      PtRelLightUntagForSystematic[tg][tr][is][ptb][etab]->Add(ThisUntag);
	      
	    }
	    
	  }
	  
	}

    }

  }

  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    NormalizeLightHistogramsToEntries();

  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_" + DataType + Configuration + Selection + ".root";
  
  SaveLightHistograms(OutputFileName);

}

void PtRelAnalyzer::BookSystem8Histograms() {
  
  std::cout << "Booking System8 histograms" << std::endl;

  TString ThisHistoName;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) {
      
      for (int is = 0; is<nSystematics; is++) {
	
	ThisHistoName = HistogramName("Observed_DataType", ptb, nb, -1, is, -1, -1);  
	Observed[ptb][nb][is] = new TH2D(ThisHistoName, ThisHistoName, nTriggers, 0., nTriggers, nTriggers, 0., nTriggers);
	
      }
      
      for (int tr = 0; tr<nTriggers; tr++) {
	
	for (int is = 0; is<nSystematics; is++)	{
	  
	  ThisHistoName = HistogramName("muonPt_DataType", ptb, nb, tr, is, -1, -1);	  
	  MuPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	  
	  ThisHistoName = HistogramName("jetPt_DataType", ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, int(MaxPtRelPtEdge), 0., MaxPtRelPtEdge); 
	  
	  ThisHistoName = HistogramName("jetEta_DataType", ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.);
	  
	  ThisHistoName = HistogramName("nGoodPV_DataType", ptb, nb, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.); 
		  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {

	      if (tg==0)
		ThisHistoName = HistogramName("ntag_pT", ptb, nb, tr, is, -1, -1);
	      else ThisHistoName = HistogramName("ntag_pT", ptb, nb, tr, is, tg-1, -1);
	      if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
	      if (tg==0) ThisHistoName.ReplaceAll("tag", "");
	      
	      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName,         ThisHistoName,         nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName + "_b",  ThisHistoName + "_b",  nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName + "_cl", ThisHistoName + "_cl", nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);

	    } 
	  
        }

      }
      
    }

}

void PtRelAnalyzer::ResetSystem8Histograms(TString DataType) {

  std::cout << "Resetting System8 histograms" << std::endl; 

  TString ThisHistoName;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int is = 0; is<nSystematics; is++)	{

	ThisHistoName = HistogramName("Observed_" + DataType, ptb, nb, -1, is, -1, -1);
	Observed[ptb][nb][is]->SetNameTitle(ThisHistoName, ThisHistoName);
	Observed[ptb][nb][is]->Reset();
	
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  ThisHistoName = HistogramName("muonPt_" + DataType, ptb, nb, tr, is, -1, -1);	
	  MuPtForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  MuPtForWeighting[tr][is][ptb][nb]->Reset();
	  	  
	  ThisHistoName = HistogramName("jetPt_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetPtForWeighting[tr][is][ptb][nb]->Reset();
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetEtaForWeighting[tr][is][ptb][nb]->Reset();
	  
	  ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, nb, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  PVMultiplicity[tr][is][ptb][nb]->Reset();
		  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {

	      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Reset();
	      System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Reset();
	      System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Reset();

	    } 
	  
        } 

      }
  
}

void PtRelAnalyzer::NormalizeSystem8HistogramsToEntries() {

  std::cout << "Normalizing System8 histograms to entries" << std::endl; 

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int is = 0; is<nSystematics; is++)	{

	float HistoEntries, HistoIntegral;

	HistoEntries = Observed[ptb][nb][is]->GetEntries();
	HistoIntegral = Observed[ptb][nb][is]->Integral();
	Observed[ptb][nb][is]->Scale(HistoEntries/HistoIntegral);
	
	for (int tr = 0; tr<nTriggers; tr++) {

	  HistoEntries = MuPtForWeighting[tr][is][ptb][nb]->GetEntries();
	  HistoIntegral = MuPtForWeighting[tr][is][ptb][nb]->Integral(0,  MuPtForWeighting[tr][is][ptb][nb]->GetNbinsX()+1);

	  MuPtForWeighting[tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	  JetPtForWeighting[tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	  JetEtaForWeighting[tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	  PVMultiplicity[tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	  		  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {
	  
	      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	      System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	      System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Scale(HistoEntries/HistoIntegral);
	      
	    } 
	  
        } 
	
      }
  
}

void PtRelAnalyzer::FillSystem8Histograms(TString DataType, int DataRange) {

  ResetSystem8Histograms(DataType);

  std::cout << "Filling System8 histograms" << std::endl;

  GetKinematicWeights(DataType);
  GetPileUpWeights(DataType);
  GetBTemplateCorrections();
  GetTriggerPrescales(DataType);
  if (!DataType.Contains("QCD")) ReadJetEnergyUncertainty("_DATA", JECVersionData, "_AK4PFchs");
  else ReadJetEnergyUncertainty("_MC", JECVersionMC, "_AK4PFchs");
  
  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
  GetDataRangeInformation(DataType, DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);

  int nEventsFromTrees = 0;
  
  for (int tf = FirstTree; tf<=nTrees; tf++) {
    
    TString FileName = GetFileName(FileDirectoryName, tf);
    
    TFile *ThisTree = TFile::Open(FileName);
    
    if (!ThisTree) continue;

    TTree *tchain = GetChain(ThisTree, false);
    
    if (DataType=="BTagMu") std::cout << "  Getting " << DataType << "   sample: BTagMu    " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
    if (DataType=="QCDMu") std::cout << "  Getting " << DataType << "    sample: MuPt5 " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
    
    Int_t nentries = (Int_t)tchain->GetEntries();
	    
    nEventsFromTrees += nentries;
    
    for (Int_t i = 0; i<nentries; i++) {
      
      tchain->GetEntry(i);
      
      if (!DataType.Contains("QCD")) {
	if (PUWeighting.Contains("ICHEP2016") && Run>276811) continue;
      }
      
      int iMu; 
      int jMu = GetPFMuonJet(&iMu);
      
      if (jMu>=0 && nPV>0 && Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1] && Jet_looseID[jMu]==1) {
	
	int ptBin = -1;
	for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	  if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	
	int tr = -1;
	
	bool PassTrigger[nTriggers];
	
	for (int trg = 0; trg<nTriggers; trg++) {
	  
	  PassTrigger[trg] = false;
	  
	  if (Jet_pt[jMu]>MinPtJetTrigger[trg] && Jet_pt[jMu]<MaxPtJetTrigger[trg]) {
	    
	    bool FiredTrigger = false;
	    if (DataType=="BTagMu" || !Selection.Contains("TrgEmul")) {
	      if (PassTriggerBit(TriggerName[trg])) FiredTrigger = true;
	      if (FiredTrigger && (Selection.Contains("TrgConf") || Selection.Contains("TrgEmul"))) 
		if (!PassTriggerEmulation(trg, jMu)) FiredTrigger = false;
	    } else {
	      if (PassTriggerEmulation(trg, jMu)) FiredTrigger = true;
	    }
	    
	    if (FiredTrigger) {
	      
	      PassTrigger[trg] = true;
	      tr = trg;
	      
	    }
	    
	  }
	  
	}
	
	if (DataType=="QCDMu") 
	  if (!PassPtHat("MuEnrichedPt5", jMu)) tr = -1;
	
	if (tr>=0) {

	  bool EventJetAwayPass = false;

	  int aJet = GetAwayJet("NONE", jMu, 0., false);
	  
	  if (aJet>=0) {
	    
	    float PtAwayJetCut = PtAwayJet[tr];
	    if (TriggerName[tr]=="_DiJet20" && Jet_pt[jMu]>=50.) PtAwayJetCut = 30.;

	    if (Jet_pt[aJet]>PtAwayJetCut)
	      EventJetAwayPass = true;
	    
	  }
	  
	  if (EventJetAwayPass) {

	    int JetEtaBin = 0;
	    for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) 
	      if (fabs(Jet_eta[jMu])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[jMu])<PtRelEtaEdge[jeta]) JetEtaBin = jeta;
	    
	    int SavePtBin = ptBin;
	    
	    for (int is = 0; is<nSystematics; is++) {
	      
	      ptBin = SavePtBin;
	      
	      float TrackPtCut =  MuonPtCut[ptBin];
	      if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	      if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
	      
	      float TrackDRCut = 999., TrackMinDRCut = 0.;
	      GetDRCuts(SystematicName[is], Jet_pt[jMu], &TrackDRCut, &TrackMinDRCut);
	      
	      double JetMuonDR = DeltaR(Jet_eta[jMu], Jet_phi[jMu], PFMuon_eta[iMu], PFMuon_phi[iMu]);
	      
	      if (PFMuon_pt[iMu]>TrackPtCut && JetMuonDR<TrackDRCut  && JetMuonDR>=TrackMinDRCut) {
		
		double ThisJetWeight = 1.;

		if (DataType=="QCDMu") {
		  
		  int ipt = Jet_pt[jMu]; 
		  if (ipt>=PtRelPtEdge[nPtRelPtBins]) ipt = PtRelPtEdge[nPtRelPtBins] - 1;
		  int imu = PFMuon_pt[iMu]; if (imu>=60) imu = 59;
		  int ieta = JetEtaBin;
		  if (KinWeighting.Contains("KinEtaAfterPtBins") || KinWeighting.Contains("KinEtaBins"))
		    ieta = (Jet_eta[jMu]+2.4)/0.1;
		  
		  int EventPileUp = nPUtrue;
		  if (PUWeighting.Contains("PV") || PUWeighting.Contains("PSV")) EventPileUp = nPV;
		  if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
		  
		  ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][imu][ieta]*PileUpWeight[tr][EventPileUp][is];
		  
		} else {

		  ThisJetWeight = TriggerPrescaleWeight(DataType, tr);

		}
		
		if (DataType=="QCDMu" && SystematicName[is].Contains("_GluonSplitting")) {
		  if (fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_GluonSplittingB")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 5)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==4 && SystematicName[is].Contains("_GluonSplittingC")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 4)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==21 && SystematicName[is].Contains("_GluonSplittingG")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 21)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } 
		}
		
		if (DataType=="QCDMu" && fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_BTemplates")) {
		  int iB = GetBHadron(jMu);
		  if (iB>=0) {
		    float EnergyFraction = BHadron_pT[iB]/Jet_genpt[jMu]; int efbin = EnergyFraction/0.02;
		    if (efbin>=0 && efbin<100) {
		      if (SystematicName[is].Contains("Minus")) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][0];
		      if (SystematicName[is].Contains("Plus") ) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][1];
		    }
		  }
		}
		/*
		  if (DataType=="QCDMu" && SystematicName[is].Contains("_MuPtWeight")) {
		  
		  int ThisMuonPtBin = PFMuon_pt[iMu]; if (ThisMuonPtBin>=60) ThisMuonPtBin = 59;
		  ThisJetWeight *= MuPtWeight[ptBin][ThisMuonPtBin];
		  
		  }*/
			
		if (SystematicName[is].Contains("_JEU")) {
		  
		  int td = 0;
		  if (DataType=="QCDMu" || DataType=="QCD") td = 1;
		  float ThisJEU = GetJEU(td, Jet_eta[jMu], Jet_pt[jMu]);
		  if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		  if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		  
		  int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		  float ScaledJetPt = Jet_pt[jMu]*(1. + JEUSign*ThisJEU);
		  
		  ptBin = -1;
		  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		    if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		  
		}
		
		if (ptBin==-1) continue;
		
		double Discriminator = PFMuon_ptrel[iMu];
		
		if (Discriminator==0.) continue;
	       
		TString ThisAwayTaggerName = "JBPM4";//"TCHPM";
		//if (Jet_pt[jMu]<50.) ThisAwayTaggerName = "JBPM7"; 
		if (Jet_pt[jMu]>=140.) ThisAwayTaggerName = "JBPM5"; 
		for (int at = 0; at<nAwayTaggers; at++) 
		  if (SystematicName[is].Contains(AwayTaggerName[at]))
		    ThisAwayTaggerName = AwayTaggerName[at];
		
		bool IsPSet = IsTaggedJet(aJet, ThisAwayTaggerName);
		
		for (int nb = 0; nb<nPtRelEtaBins; nb++) 
		  if (nb==0 || JetEtaBin==nb) {
		    
		    if (DataType=="BTagMu")
		      for (int trg1 = 0; trg1<nTriggers; trg1++) 
			for (int trg2 = 0; trg2<nTriggers; trg2++) 
			  if (PassTrigger[trg1] || PassTrigger[trg2]) 
			    Observed[ptBin][nb][is]->Fill(trg1, trg2, ThisJetWeight); 	     
	     
		    PVMultiplicity[tr][is][ptBin][nb]->Fill(nPV, ThisJetWeight);      
		    MuPtForWeighting[tr][is][ptBin][nb]->Fill(PFMuon_pt[iMu], ThisJetWeight);
		    JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[jMu], ThisJetWeight); 
		    JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[jMu], ThisJetWeight);
		    
		    for (int st = 0; st<2; st++)
		      if (st==0 || IsPSet) 
			for (int tg = 0; tg<nTaggers+1; tg++) 
			  if (tg==0 || IsTaggedJet(jMu, TaggerName[tg-1])) {
			    
			    if (DataType=="BTagMu")
			      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			    else {
			      
			      if (fabs(Jet_flavour[jMu])==5) 
				System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			      
			      else if (fabs(Jet_flavour[jMu])!=5) 
				System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			      
			    }
			    
			  }
		    
		  } // Loop on eta bins
		
	      } // Muon pt and DR cuts
	      
	    } // Loop on systematics
	    
	  } // Has jet away
	  
	} // Pass trigger
	
      } // Muon-jet selection
      
    } // Loop on events
    
    ThisTree->Close();
	    
  } // Loop on files
  
  std::cout << "  nEventsFromTrees " << nEventsFromTrees << std::endl;
  
  TString Configuration = PUWeighting + KinWeighting; 
  if (!DataType.Contains("QCD")) {
    if (DataType=="BTagMu") Configuration = PUWeighting;
    Configuration.ReplaceAll("_PU", "_");
    if (DataType=="BTagMu") Configuration.ReplaceAll("_PV", "_");
    Configuration.ReplaceAll("_PSV", "_PS");
  }

  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    if ((DataType=="BTagMu" && nBTagMuRanges==1) || (DataType=="JetHT" && nJetRunRanges==1))
      NormalizeSystem8HistogramsToEntries();
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataRangeName + Configuration + Selection + ".root";
  
  SaveSystem8Histograms(OutputFileName, DataType);
  
}

void PtRelAnalyzer::SaveSystem8Histograms(TString OutputFileName, TString DataType) {

  std::cout << "Saving System8 histograms" << std::endl;
 
  TFile *OutFile = new TFile(OutputFileName, "recreate");

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int is = 0; is<nSystematics; is++)	{

	if (DataType=="BTagMu") 
	  Observed[ptb][nb][is]->Write();

	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  MuPtForWeighting[tr][is][ptb][nb]->Write();
	  JetPtForWeighting[tr][is][ptb][nb]->Write();
	  JetEtaForWeighting[tr][is][ptb][nb]->Write();
	  PVMultiplicity[tr][is][ptb][nb]->Write();
	  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {
	      
	      if (DataType=="BTagMu") System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Write();
	      else {

		System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Write();
		System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Write();
		
	      }
	      
	    } 
	  
        } 
	
      }
  
  OutFile->Close();

  std::cout << "Exiting" << std::endl;
  
}

void PtRelAnalyzer::MergeSystem8Histograms(TString DataType, int FirstBin, int LastBin) {

  ResetSystem8Histograms(DataType);

  std::cout << "Merging System8 histograms: " << DataType << std::endl;
  
  if (LastBin==100) LastBin = GetNumberOfDataRanges(DataType);
  
  for (int ph = FirstBin; ph<LastBin; ph++) {

    TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
    GetDataRangeInformation(DataType, ph, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
    
    std::cout << "  Merging " << DataRangeName << std::endl;
    
    TString InputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";
    if (DataType=="BTagMu") {
      InputFileName.ReplaceAll("_PV", "_");
      InputFileName.ReplaceAll("_PU", "_");
      InputFileName.ReplaceAll("_PSV", "_PS");
      InputFileName.ReplaceAll(KinWeighting, "");
    }
    TFile *HistogramFile = TFile::Open(InputFileName); 

    TString ThisHistoName;

    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

      std::cout << "    Merging " << PtRelPtBin[ptb] << std::endl;

      for (int etab = 0; etab<nPtRelEtaBins; etab++) 
	for (int is = 0; is<nSystematics; is++) {

	  if (DataType=="BTagMu") {
	    
	    ThisHistoName = HistogramName("Observed_" + DataType, ptb, etab, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) HistogramFile->Get(ThisHistoName);
	    Observed[ptb][etab][is]->Add(ThisObserved);
	  
	  }

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    ThisHistoName = HistogramName("muonPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisMuPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    MuPtForWeighting[tr][is][ptb][etab]->Add(ThisMuPt);
	    
	    ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetPtForWeighting[tr][is][ptb][etab]->Add(ThisJetPt);
	    
	    ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetEta = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetEtaForWeighting[tr][is][ptb][etab]->Add(ThisJetEta);
	    
	    ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisPVMult = (TH1D*) HistogramFile->Get(ThisHistoName);
	    PVMultiplicity[tr][is][ptb][etab]->Add(ThisPVMult);
		  
	    for (int st = 0; st<2; st++)
	      for (int tg = 0; tg<nTaggers+1; tg++) {
		
		if (tg==0)
		  ThisHistoName = HistogramName("ntag_pT", ptb, etab, tr, is, -1, -1);
		else ThisHistoName = HistogramName("ntag_pT", ptb, etab, tr, is, tg-1, -1);
		if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
		if (tg==0) ThisHistoName.ReplaceAll("tag", "");
		
		if (DataType=="BTagMu") {

		  TH1D *ThisData = (TH1D*) HistogramFile->Get(ThisHistoName);
		  System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][etab]->Add(ThisData);
		
		} else {
		  
		  TH1D *ThisBJet = (TH1D*) HistogramFile->Get(ThisHistoName + "_b");
		  System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][etab]->Add(ThisBJet);
		  
		  TH1D *ThisLJet = (TH1D*) HistogramFile->Get(ThisHistoName + "_cl");
		  System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][etab]->Add(ThisLJet);
		  
		}
		
	      }
	    
	  }
      
	}
      
    }

  }

  if (PUWeighting.Contains("_PS") && !DataType.Contains("QCD"))
    NormalizeSystem8HistogramsToEntries();
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + PUWeighting + KinWeighting + Selection + ".root";
  if (DataType=="BTagMu") {
    OutputFileName.ReplaceAll("_PV", "_");
    OutputFileName.ReplaceAll("_PU", "_");
    OutputFileName.ReplaceAll("_PSV", "_PS");
    OutputFileName.ReplaceAll(KinWeighting, "");
  }

  SaveSystem8Histograms(OutputFileName, DataType);

}

void PtRelAnalyzer::BuildSystem8Templates() {

  BookSystem8Templates();

  std::cout << "Building System8 templates" << std::endl;

  TString BTagMuConfiguration = PUWeighting;
  BTagMuConfiguration.ReplaceAll("_PU", "_");
  BTagMuConfiguration.ReplaceAll("_PV", "_");
  BTagMuConfiguration.ReplaceAll("_PSV", "_PS");
 
  TString DataFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_BTagMu" + BTagMuConfiguration + Selection + ".root";
  TFile *DataFile = TFile::Open(DataFileName); 
  
  TString QCDFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_QCDMu" + PUWeighting + KinWeighting + Selection + ".root";
  TFile *QCDFile = TFile::Open(QCDFileName); 

  TString ThisHistoName;
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++) {
    
    std::cout << "    Merging " << FitPtBin[fpt] << std::endl;

    float LowFitPtBinEdge = FitPtEdge[fpt];
    float HighFitPtBinEdge = FitPtEdge[fpt+1];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      float MiddlePtBin = (PtRelPtEdge[bpt] + PtRelPtEdge[bpt+1])/2.;
      
      if (MiddlePtBin>LowFitPtBinEdge && MiddlePtBin<HighFitPtBinEdge) {
	    
	for (int nb = 0; nb<nPtRelEtaBins; nb++) { 
	      
	  for (int is = 0; is<nSystematics; is++) {
	    
	    ThisHistoName = HistogramName("Observed_BTagMu", bpt, nb, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) DataFile->Get(ThisHistoName);
	    
	    for (int tr = 0; tr<nTriggers; tr++)
	      if (AllowedTrigger[bpt][tr]) {

		float LuminosityWeight =  1.;
		
		ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetEta = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataPVMult = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetEta = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDPVMult = (TH1D*) QCDFile->Get(ThisHistoName);
		
		//if (dt==1)
		//if (ThisMuPt->Integral(0, 61)>0.)
		//  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		
		int tr1 = -1, tr2 = -1;
		for (int trg = 0; trg<nTriggers; trg++)
		  if (trg!=tr)
		    if (AllowedTrigger[bpt][trg]) {
		      
		      if (tr1==-1) tr1 = trg;
		      else if (tr2==-1) tr2 = trg;
		      else 
			std::cout << "PtRelAnalyzer::BuildSystem8Templates Warning:" 
			     << " too many triggers merged in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		      
		    }
		
		if (tr1==-1) {
		  
		  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		  
		} else if (tr2==-1) {
		  
		  int T1 = tr, T2 = tr1;
		  if (tr1<tr) { T1 = tr1; T2 = tr; } 
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT2 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) - nT1;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else {
		    
		    std::cout << "PtRelAnalyzer::BuildSystem8Templates Warning:" 
			 << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		    
		  }
		  
		} else {
		  
		  int T1 = tr, T2 = tr1, T3 = tr2;
		  
		  if (tr<tr1 && tr<tr2 && tr1>tr2) { T1 = tr;  T2 = tr2; T3 = tr1; }
		  if (tr>tr1 && tr<tr2 && tr1<tr2) { T1 = tr1; T2 = tr;  T3 = tr2; }
		  //if (tr<tr1 && tr>tr2 && tr1<tr2) { }
		  if (tr>tr1 && tr>tr2 && tr1<tr2) { T1 = tr1; T2 = tr2; T3 = tr;  }
		  //if (tr>tr1 && tr<tr2 && tr1>tr2) { }
		  if (tr<tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr;  T3 = tr1; }
		  if (tr>tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr1; T3 = tr;  }
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  float nT2 = (ThisObserved->GetBinContent(T1+1, T2+1) - nT1)*(1-(ThisObserved->GetBinContent(T3+1, T3+1)/ThisObserved->GetBinContent(T2+1, T2+1))*(TriggerLuminosity[T2]/TriggerLuminosity[T3]));
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T3, is, -1, -1);
		  TH1D *ThisDataMuPt3 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT3 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) + ThisDataMuPt3->Integral(0, 61) - nT1 - nT2;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else if (tr==T3) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T3, is, -1, -1);
		    TH1D *ThisQCDMuPt3 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT3/ThisQCDMuPt3->Integral(0, 61);
		    
		  } else {
		    
		    std::cout << "PtRelAnalyzer::BuildSystem8Templates Warning:" 
			   << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		    
		  }
		  
		} 
		
		MuonPt[0][is][fpt][nb]->Add(ThisDataMuPt);
		JetPt[0][is][fpt][nb]->Add(ThisDataJetPt);
		JetEta[0][is][fpt][nb]->Add(ThisDataJetEta);
		PVEvent[0][is][fpt][nb]->Add(ThisDataPVMult);
		
		if (LuminosityWeight<0. || LuminosityWeight>1000000000.) {

		  std::cout << "PtRelAnalyzer::BuildSystem8Templates: LuminosityWeight " << LuminosityWeight
		       << " for " << TriggerName[tr] << " in " << PtRelPtBin[bpt] << std::endl;

		  LuminosityWeight = 0.;

		}

		ThisQCDMuPt->Scale(LuminosityWeight);
		MuonPt[1][is][fpt][nb]->Add(ThisQCDMuPt);

		ThisQCDJetPt->Scale(LuminosityWeight);
		JetPt[1][is][fpt][nb]->Add(ThisQCDJetPt);

		ThisQCDJetEta->Scale(LuminosityWeight);
		JetEta[1][is][fpt][nb]->Add(ThisQCDJetEta);

		ThisQCDPVMult->Scale(LuminosityWeight);
		PVEvent[1][is][fpt][nb]->Add(ThisQCDPVMult);
		
		for (int st = 0; st<2; st++)
		  for (int tg = 0; tg<nTaggers+1; tg++) {
		
		    if (tg==0)
		      ThisHistoName = HistogramName("ntag_pT", bpt, nb, tr, is, -1, -1);
		    else ThisHistoName = HistogramName("ntag_pT", bpt, nb, tr, is, tg-1, -1);
		    if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
		    if (tg==0) ThisHistoName.ReplaceAll("tag", "");
		    
		    TH1D *ThisDataSystem8 = (TH1D*) DataFile->Get(ThisHistoName);
		    System8[st*nTaggers+st+tg][is][fpt][nb][0]->Add(ThisDataSystem8);
		    
		    TH1D *ThisBJetSystem8 = (TH1D*) QCDFile->Get(ThisHistoName + "_b");
		    ThisBJetSystem8->Scale(LuminosityWeight);
		    System8[st*nTaggers+st+tg][is][fpt][nb][1]->Add(ThisBJetSystem8);
		    
		    TH1D *ThisLJetSystem8 = (TH1D*) QCDFile->Get(ThisHistoName + "_cl");
		    ThisLJetSystem8->Scale(LuminosityWeight);
		    System8[st*nTaggers+st+tg][is][fpt][nb][2]->Add(ThisLJetSystem8);
		    
		  }		
		
	      }
	    
	  }
	  
	}
	
      }
      
    }
    
  }

  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + PUWeighting + KinWeighting + Selection + Production + ".root";
  
  SaveSystem8Templates(TemplateFileName);

}
  
void PtRelAnalyzer::BookSystem8Templates() {
  
  std::cout << "Booking System8 templates" << std::endl;

  TString ThisHistoName;
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)       
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<2; dt++) {

	  TString DataType = "BTagMu";
	  if (dt==1) DataType = "QCDMu";
	  
	  ThisHistoName =  ThisHistoName = HistogramName("jetPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, int(MaxPtRelPtEdge), 0., MaxPtRelPtEdge); 
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetEta[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.); 

	  ThisHistoName = HistogramName("muonPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  MuonPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	  
	  ThisHistoName = HistogramName("PV_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  PVEvent[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.);
	  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {

	      if (tg==0)
		ThisHistoName = HistogramName("ntag_pT", FitPtBin[fpt], nb, -1, is, -1, -1);
	      else ThisHistoName = HistogramName("ntag_pT", FitPtBin[fpt], nb, -1, is, tg-1, -1);
	      if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
	      if (tg==0) ThisHistoName.ReplaceAll("tag", "");
	      
	      if (dt==0) System8[st*nTaggers+st+tg][is][fpt][nb][0] = new TH1D(ThisHistoName,         ThisHistoName,         nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      else {
		System8[st*nTaggers+st+tg][is][fpt][nb][1] = new TH1D(ThisHistoName + "_b",  ThisHistoName + "_b",  nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
		System8[st*nTaggers+st+tg][is][fpt][nb][2] = new TH1D(ThisHistoName + "_cl", ThisHistoName + "_cl", nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      }

	    }
	  
	}
  
}

void PtRelAnalyzer::SaveSystem8Templates(TString TemplateFileName) {
  
  std::cout << "Saving System8 templates" << std::endl;

  TFile *OutFile = new TFile(TemplateFileName, "recreate");

  TString ThisHistoName;
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)       
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<2; dt++) {
	  	  
	  JetPt[dt][is][fpt][nb]->Write();
	  JetEta[dt][is][fpt][nb]->Write();
	  MuonPt[dt][is][fpt][nb]->Write();
	  PVEvent[dt][is][fpt][nb]->Write();
	  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {
	      
	      if (dt==0) System8[st*nTaggers+st+tg][is][fpt][nb][0]->Write();
	      else {
		System8[st*nTaggers+st+tg][is][fpt][nb][1]->Write();
		System8[st*nTaggers+st+tg][is][fpt][nb][2]->Write();
	      }

	    }
	  
	}
  
  OutFile->Close();
  
  std::cout << "Exiting" << std::endl;
  
}

int PtRelAnalyzer::GetAwayJet(TString ThisAwayTaggerName, int jMu, float AwayDRCut, bool IsUnique) {
  
  int aJet = -1;

  int nAwayJets = 0;
  float PtLeadingAwayJet = -1;
    
  for (int ijet = 0; ijet<nJet; ijet++) 
    if (Jet_pt[ijet]>PtRelPtEdge[0] && ijet!=jMu && fabs(Jet_eta[ijet])<2.4) {
    
      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Jet_eta[jMu], Jet_phi[jMu])>AwayDRCut) {
	
	bool pass = IsTaggedJet(ijet, ThisAwayTaggerName);
		
	if (pass) {
	
	  if (Jet_pt[ijet]>PtLeadingAwayJet) {
	    
	    aJet = ijet;
	    PtLeadingAwayJet = Jet_pt[ijet];

	  }

	  nAwayJets++;

	}
      
      }

    }

  if (IsUnique && nAwayJets>1) aJet = -1;

  return aJet;
   
}

int PtRelAnalyzer::GetAwayStdJet(TString ThisAwayTaggerName, int jMu, float AwayDRCut, bool IsUnique, int TriggerIdx) {
  
  int aJet = -1;

  int nAwayJets = 0;
  float PtLeadingAwayJet = -1;
    
  for (int ijet = 0; ijet<nStdJet; ijet++) 
    if (StdJet_pt[ijet]>PtAwayJet[TriggerIdx] && fabs(StdJet_eta[ijet])<2.4) {
      
      if (DeltaR(StdJet_eta[ijet], StdJet_phi[ijet], Jet_eta[jMu], Jet_phi[jMu])>AwayDRCut) {
	
	bool pass = IsTaggedJet(ijet, ThisAwayTaggerName, "Standard");
		
	if (pass) {
	
	  if (StdJet_pt[ijet]>PtLeadingAwayJet) {
	    
	    aJet = ijet;
	    PtLeadingAwayJet = StdJet_pt[ijet];

	  }

	  nAwayJets++;

	}
      
      }

    }

  if (IsUnique && nAwayJets>1) aJet = -1;

  return aJet;
   
}

void PtRelAnalyzer::GetDRCuts(TString ThisSystematicName, float MuonJetPt, float *TrackDRCut, float *TrackMinDRCut) {

  *TrackDRCut    = 999.; 
  *TrackMinDRCut = 0.;
  
  if (!ThisSystematicName.Contains("_MuDR")) {

    if (TemplateVariable=="PtRel") {

      if (MuonJetPt<30) *TrackDRCut    = 0.20; 
      else if (MuonJetPt<80) *TrackDRCut    = 0.15; 
      else *TrackDRCut    = 0.12;
    
    } else if (TemplateVariable=="System8") {

      *TrackDRCut    = 0.40;
    
    } else {

      *TrackDRCut    = 999.; 
      *TrackMinDRCut = 0.;
      
    } 

  } else if (ThisSystematicName=="_MuDRMinus") {
    
    if (TemplateVariable=="PtRel") {

      if (MuonJetPt<30) *TrackDRCut    = 0.15; 
      else if (MuonJetPt<80) *TrackDRCut    = 0.12; 
      else *TrackDRCut    = 0.09; 

    } else if (TemplateVariable=="System8") {

      *TrackDRCut    = 0.30;
    
    } else {
    
      if (MuonJetPt<30) *TrackDRCut    = 0.20; 
      else if (MuonJetPt<80) *TrackDRCut    = 0.15; 
      else *TrackDRCut    = 0.12; 
      
    }
    
  } else if (ThisSystematicName=="_MuDRPlus") {
    
    *TrackDRCut    = 999.; 
    *TrackMinDRCut = 0.;
    
  } else {
    
    if (ThisSystematicName=="_MuDR003") *TrackDRCut = 0.03;
    else if (ThisSystematicName=="_MuDR004") *TrackDRCut = 0.04;
    else if (ThisSystematicName=="_MuDR005") *TrackDRCut = 0.05;
    else if (ThisSystematicName=="_MuDR006") *TrackDRCut = 0.06;
    else if (ThisSystematicName=="_MuDR007") *TrackDRCut = 0.07;
    else if (ThisSystematicName=="_MuDR008") *TrackDRCut = 0.08;
    else if (ThisSystematicName=="_MuDR009") *TrackDRCut = 0.09;
    else if (ThisSystematicName=="_MuDR010") *TrackDRCut = 0.10;
    else if (ThisSystematicName=="_MuDR012") *TrackDRCut = 0.12;
    else if (ThisSystematicName=="_MuDR015") *TrackDRCut = 0.15;
    else if (ThisSystematicName=="_MuDR020") *TrackDRCut = 0.20;
    else if (ThisSystematicName=="_MuDR050") *TrackDRCut = 999.;
    else if (ThisSystematicName=="_MuDRp08") *TrackMinDRCut = 0.08;
    else if (ThisSystematicName=="_MuDRp10") *TrackMinDRCut = 0.10; 
    else if (ThisSystematicName=="_MuDRp15") *TrackMinDRCut = 0.15; 
    else std::cout << "Wrong choise for MuDRCut" << std::endl;
    
  }
  
}
 
void PtRelAnalyzer::GetPileUpWeights(TString DataType) { 

  cout << "Getting pileup weights" << endl;

  if (DataType=="BTagMu") {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int pu = 0; pu<nMaxPU; pu++)
	for (int is = 0; is<nSystematics; is++)
	  PileUpWeight[tr][pu][is] = 1.;
      
  } else if (PUWeighting.Contains("_PV") || PUWeighting.Contains("_PU") || PUWeighting.Contains("_PS")) {
    
    for (int tr = 0; tr<nTriggers; tr++) {
      
      TString PUWeightFileName = "./Weights/PileUp/PileUpWeights" + TriggerName[tr] + "_" + DataType + PUWeighting + ".txt";
      if (PUWeighting.Contains("_PV") || PUWeighting.Contains("_PSV")) PUWeightFileName.ReplaceAll("PileUpWeights", "PVMultTriggered");

      ifstream PUWeightsFile; PUWeightsFile.open(PUWeightFileName);
      if (!PUWeightsFile)
	throw std::invalid_argument("PU weights " + PUWeightFileName + " not found!");

      int nMaxPUFile = nMaxPU;
      if (PUWeighting.Contains("_PV")) nMaxPUFile = 60;
      for (int npu = 0; npu<nMaxPUFile; npu++) {

	int pum; float weightpu;
	PUWeightsFile >> pum >> weightpu;

	for (int is = 0; is<nSystematics; is++)  
	  PileUpWeight[tr][pum][is] = weightpu;
	
      }

      if ((PUWeighting.Contains("_PS") || PUWeighting.Contains("_PU")) && !PUWeighting.Contains("_PSV"))
	for (int is = 0; is<nSystematics; is++)  
	  if (SystematicName[is].Contains("_PileUp")) {
	    
	    TString PUSyst = "_p05";
	    if (SystematicName[is].Contains("Down")) PUSyst = "_m05";
	    ifstream PUSystWeightsFile; PUSystWeightsFile.open("./Weights/PileUp/PileUpWeights" + TriggerName[tr] + "_" + DataType + PUWeighting + PUSyst + ".txt");
	    if (!PUSystWeightsFile)
	      throw std::invalid_argument("PU weights for systematics not found!");
	    
	    for (int npu = 0; npu<nMaxPU; npu++) {
	      
	      int pum; float weightpu;
	      PUSystWeightsFile >> pum >> weightpu;
	    
	      PileUpWeight[tr][pum][is] = weightpu;
	      
	    }
	    
	  }
      
    }

  } else {
    
    throw std::invalid_argument("PU reweighting option " + PUWeighting + " not supported!");
    
  }
  
}

void PtRelAnalyzer::GetKinematicWeights(TString DataType) {

  std::cout << "Getting kinematic weights" << std::endl;

  if (DataType=="BTagMu") {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<PtRelPtEdge[nPtRelPtBins]; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 1; nb<nPtRelEtaBins; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb-1] = 1.;
    
  } else if (KinWeighting.Contains("_KinPtBins")) {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<PtRelPtEdge[nPtRelPtBins]; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      TString CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinWeighting + "_Pt_" + PtRelPtBin[bpt] + "_anyEta_" + DataType + PUWeighting + Selection + ".txt";
      
      ifstream PtWeightsFile; PtWeightsFile.open(CorrectionFileName);
      if (!PtWeightsFile)
	throw std::invalid_argument("Kinemtic weights not found!");
      
      int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
      for (int stpt = 0; stpt<nPtSteps; stpt++) {
	
	float PtStep, ThisPtWeight; 
	PtWeightsFile >> PtStep >> ThisPtWeight;	 
	
	for (int tr = 0; tr<nTriggers; tr++)
	  for (int is = 0; is<nSystematics; is++) 
	    for (int imu = 0; imu<60; imu++) 
	      for (int nb = 0; nb<nPtRelEtaBins; nb++) 
		KinematicWeight[tr][is][int(PtStep)][imu][nb] = ThisPtWeight;
	
      }
      
    }
    
  } else if (KinWeighting.Contains("_KinEtaAfterPtBins")) {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<PtRelPtEdge[nPtRelPtBins]; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins-1; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      TString KinPtWeighting = KinWeighting; KinPtWeighting.ReplaceAll("EtaAfter", "");
      TString CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinPtWeighting + "_Pt_" + PtRelPtBin[bpt] + "_anyEta_" + DataType + PUWeighting + Selection + ".txt";

      ifstream PtWeightsFile; PtWeightsFile.open(CorrectionFileName);
      if (!PtWeightsFile)
	throw std::invalid_argument("Kinemtic weights not found!");
      
      int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
      for (int stpt = 0; stpt<nPtSteps; stpt++) {
	
	float PtStep, ThisPtWeight; 
	PtWeightsFile >> PtStep >> ThisPtWeight;	 
	
	for (int tr = 0; tr<nTriggers; tr++)
	  for (int is = 0; is<nSystematics; is++) 
	    for (int imu = 0; imu<60; imu++) 
	      for (int ieta = 0; ieta<48; ieta++)
		KinematicWeight[tr][is][int(PtStep)][imu][ieta] = ThisPtWeight;
	
      }
      
      CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinWeighting + "_Eta_" + PtRelPtBin[bpt] + "_anyEta_" + DataType + PUWeighting + Selection + ".txt";

      ifstream EtaWeightsFile; EtaWeightsFile.open(CorrectionFileName);
      if (!EtaWeightsFile)
	throw std::invalid_argument("Kinemtic weights not found!");
      
      int nEtaSteps = 48;
      for (int steta = 0; steta<nEtaSteps; steta++) {
	
	float EtaStep, ThisEtaWeight; 
	EtaWeightsFile >> EtaStep >> ThisEtaWeight;	 
	
	int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
	for (int stpt = PtRelPtEdge[bpt]; stpt<int(PtRelPtEdge[bpt+1]); stpt++) 
	  for (int tr = 0; tr<nTriggers; tr++)
	    for (int is = 0; is<nSystematics; is++) 
	      for (int imu = 0; imu<60; imu++)	      
		KinematicWeight[tr][is][stpt][imu][steta] *= ThisEtaWeight;
	
      }
      
    }

  } else if (KinWeighting.Contains("_KinMuonPtAndPtBins")) {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<PtRelPtEdge[nPtRelPtBins]; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins-1; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      TString KinPtWeighting = KinWeighting; KinPtWeighting.ReplaceAll("MuonPtAnd", "");
      TString CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinPtWeighting + "_Pt_" + PtRelPtBin[bpt] + "_anyEta_" + DataType + PUWeighting + Selection + ".txt";

      ifstream PtWeightsFile; PtWeightsFile.open(CorrectionFileName);
      if (!PtWeightsFile)
	throw std::invalid_argument("Kinemtic weights not found!");
      
      int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
      for (int stpt = 0; stpt<nPtSteps; stpt++) {
	
	float PtStep, ThisPtWeight; 
	PtWeightsFile >> PtStep >> ThisPtWeight;	 
	
	for (int tr = 0; tr<nTriggers; tr++)
	  for (int is = 0; is<nSystematics; is++) 
	    for (int imu = 0; imu<60; imu++) 
	      for (int ieta = 0; ieta<48; ieta++)
		KinematicWeight[tr][is][int(PtStep)][imu][ieta] = ThisPtWeight;
	
      }

      if (DataType=="QCDMu") {

	TString KinMuonPtWeighting = KinWeighting; KinMuonPtWeighting.ReplaceAll("AndPt", "");
	CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinMuonPtWeighting + "_MuonPt_" + PtRelPtBin[bpt] + "_anyEta_" + DataType + PUWeighting + Selection + ".txt";
	
	ifstream MuonPtWeightsFile; MuonPtWeightsFile.open(CorrectionFileName);
	if (!MuonPtWeightsFile)
	  throw std::invalid_argument("Kinemtic weights not found!");
	
	int nMuonPtSteps = 60;
	for (int stmpt = 0; stmpt<nMuonPtSteps; stmpt++) {
	  
	  float MuonPtStep, ThisMuonPtWeight; 
	  MuonPtWeightsFile >> MuonPtStep >> ThisMuonPtWeight;	 
	  
	  for (int stpt = PtRelPtEdge[bpt]; stpt<int(PtRelPtEdge[bpt+1]); stpt++) 
	    for (int tr = 0; tr<nTriggers; tr++)
	      for (int is = 0; is<nSystematics; is++) 
		for (int ieta = 0; ieta<48; ieta++)
		  KinematicWeight[tr][is][stpt][int(MuonPtStep)][ieta] *= ThisMuonPtWeight;
	  
	}

      }
      
    }

  } else if (KinWeighting.Contains("_KinPtInEtaBins")) {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<PtRelPtEdge[nPtRelPtBins]; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins-1; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;

    TString KinPtWeighting = KinWeighting; KinPtWeighting.ReplaceAll("InEta", "");

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      for (int nb = 0; nb<nPtRelEtaBins; nb++) {
	
        TString CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinPtWeighting + "_Pt_" + PtRelPtBin[bpt] + "_" + PtRelEtaBin[nb] + "_" + DataType + PUWeighting + Selection + ".txt";
	std::cout << CorrectionFileName << std::endl;

	ifstream PtWeightsFile; PtWeightsFile.open(CorrectionFileName);
	
	int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
	for (int stpt = 0; stpt<nPtSteps; stpt++) {
	  
	  float PtStep, ThisPtWeight; 
	  PtWeightsFile >> PtStep >> ThisPtWeight;	 

	  for (int tr = 0; tr<nTriggers; tr++)
	    for (int is = 0; is<nSystematics; is++) 
	      for (int imu = 0; imu<60; imu++) {
		
		KinematicWeight[tr][is][int(PtStep)][imu][nb] = ThisPtWeight;
		
	      }  
	  
	}
	
      }
      
    }

  } else if (KinWeighting.Contains("_KinWeights") || KinWeighting.Contains("_FullKinWeights")) {

  } else {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<PtRelPtEdge[nPtRelPtBins]; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;
    
  }
  
}
 
void PtRelAnalyzer::GetTriggerPrescales(TString DataType) { 

  cout << "Getting trigger prescales" << endl;

  if (PUWeighting.Contains("_PS") && (DataType=="BTagMu" || DataType=="JetHT")) {

    for (int tr = 0; tr<nTriggers; tr++) {

      TString ThisTriggerName = TriggerName[tr];
      if (DataType=="JetHT") ThisTriggerName = JetTriggerName[tr];

      ifstream PrescalesFile; PrescalesFile.open("./Weights/Prescales/Prescales" + ThisTriggerName + ".txt");
      if (!PrescalesFile) throw std::invalid_argument("Prescale weights not found!");

      int nPrescales = 0;

      while (PrescalesFile) {

	int ThisRunNumber, ThisLumiSection, ThisPrescale; TString LumiSectionString;
	PrescalesFile >> ThisRunNumber >> LumiSectionString >> ThisPrescale;
	ThisLumiSection = LumiSectionString.Atoi();
	
	if (ThisLumiSection>0) {

	  TriggerPrescaleRunNumber[tr][nPrescales] = ThisRunNumber;
	  TriggerPrescaleLumiSection[tr][nPrescales] = ThisLumiSection;
	  TriggerPrescale[tr][nPrescales] = ThisPrescale;
	  
	  nPrescales++;

	}

      }

      if (nPrescales>=nMaxPrescales-1) 
	std::cout << "PtRelAnalyzer::GetTriggerPrescales: Error -> need to increase nMaxPrescales " << std::endl;
      else 
	for (int ps = nPrescales; ps<nMaxPrescales; ps++)
	  TriggerPrescaleRunNumber[tr][ps] = -1;

    }

  }

}

float PtRelAnalyzer::TriggerPrescaleWeight(TString DataType, int TriggerIndex) { 

  if (!PUWeighting.Contains("_PS") || (DataType!="BTagMu" && DataType!="JetHT")) return 1.;

  float ThisTriggerPrescaleWeight = -1.;
  
  for (int ps = 0; ps<nMaxPrescales; ps++) {
    if (TriggerPrescaleRunNumber[TriggerIndex][ps]>0) {
      if (Run==TriggerPrescaleRunNumber[TriggerIndex][ps]) {
	if (LumiBlock>=TriggerPrescaleLumiSection[TriggerIndex][ps])
	  ThisTriggerPrescaleWeight = 1.*TriggerPrescale[TriggerIndex][ps];
      } else if (ThisTriggerPrescaleWeight>0.) return ThisTriggerPrescaleWeight;
    }
  }

  if (ThisTriggerPrescaleWeight==-1.)
    std::cout << "PtRelAnalyzer::TriggerPrescaleWeight: prescale not found for run " << Run << " and lumi section " << LumiBlock << ", " << TriggerName[TriggerIndex] << std::endl;
  
  return -1.;

}

void PtRelAnalyzer::GetBTemplateCorrections() {

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

    for (int ib = 0; ib<100; ib++) {

      BTemplateCorrections[ib][ptb][0] = 1.;
      BTemplateCorrections[ib][ptb][1] = 1.;

    }

    ifstream MnusCorrectionsFile; MnusCorrectionsFile.open("./Weights/BTemplatesCorrections/EnergyFraction_" + PtRelPtBin[ptb] + "_m5.txt");

    if (!MnusCorrectionsFile)
      throw std::invalid_argument("B energy weights minus not found!");

    while (MnusCorrectionsFile) {

      float xBin, efcorr;
      MnusCorrectionsFile >> xBin >> efcorr;
      
      if (efcorr>4.) efcorr = 1.;

      int ib = xBin/0.02;
      BTemplateCorrections[ib][ptb][0] = efcorr;

    }

    ifstream PlusCorrectionsFile; PlusCorrectionsFile.open("./Weights/BTemplatesCorrections/EnergyFraction_" + PtRelPtBin[ptb] + "_p5.txt");

    if (!PlusCorrectionsFile)
      throw std::invalid_argument("B energy weights plus not found!");

    while (PlusCorrectionsFile) {

      float xBin, efcorr;
      PlusCorrectionsFile >> xBin >> efcorr;
      
      if (efcorr>4.) efcorr = 1.;

      int ib = xBin/0.02;
      BTemplateCorrections[ib][ptb][1] = efcorr;

    }

  }

}

void PtRelAnalyzer::BookTemplates() {

  std::cout << "Booking templates" << std::endl;

  TString ThisHistoName;

  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)       
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<4; dt++) {

	  TString DataType = "BTagMu";
	  if (dt==1) DataType = "QCDMu";
	  if (dt==2) DataType = "JetHT";
	  if (dt==3) DataType = "QCD";
	  
	  ThisHistoName =  ThisHistoName = HistogramName("jetPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, int(MaxPtRelPtEdge), 0., MaxPtRelPtEdge); 
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetEta[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.); 

	  if (dt<2) {

	    ThisHistoName = HistogramName("muonPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	    MuonPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	    
	    ThisHistoName = HistogramName("muonDR_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	    MuonDR[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 50, 0., 0.5);
	  
	    ThisHistoName = HistogramName("PV_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	    PVEvent[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.);

	    for (int tg = 0; tg<nTaggers; tg++)
	      for (int tp = 0; tp<2; tp++)
		for (int fl = 0; fl<4; fl++) {
		  
		  if (dt==0 && fl!=0 && fl!=3) continue;
  
		  int lf = 0;
		  if (dt==0 && fl==3) lf = 4;
		  if (dt==1 && fl==0) lf = 1;
		  if (dt==1 && fl==1) lf = 2;
		  if (dt==1 && fl==2) lf = 3;
		  if (dt==1 && fl==3) lf = 4;
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, FitPtBin[fpt], nb, -1, is, tg, 10*tp + lf);
		  PtRel[dt*nTaggers+tg][is][fpt][nb][4*tp+fl] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
		  
		}

	  }

	}
  
}

void PtRelAnalyzer::BuildTemplates(bool AddTrackTemplates, TString DataFlag, TString MCFlag) {
  
  if (KinWeighting.Contains("EtaSplit")) {

    BuildTemplatesEtaSplit(AddTrackTemplates);
    return;

  }

  BookTemplates();

  std::cout << "Building templates" << std::endl;
  
  TFile *DataFile, *QCDFile, *JetFile, *IncFile;

  TString BTagMuConfiguration = PUWeighting;
  BTagMuConfiguration.ReplaceAll("_PU", "_");
  BTagMuConfiguration.ReplaceAll("_PV", "_");
  BTagMuConfiguration.ReplaceAll("_PSV", "_PS");

  TString JetHTConfiguration = PUWeighting + KinWeighting;
  JetHTConfiguration.ReplaceAll("_PU", "_");
  JetHTConfiguration.ReplaceAll("_PSV", "_PS");

  TString DataFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataFlag + BTagMuConfiguration + Selection + ".root"; 
  DataFile = TFile::Open(DataFileName); 

  TString QCDFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + MCFlag + PUWeighting + KinWeighting + Selection + ".root"; 
  QCDFile = TFile::Open(QCDFileName); 

  if (AddTrackTemplates) {

    TString JetFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_JetHT" + JetHTConfiguration + Selection + ".root";
    JetFile = TFile::Open(JetFileName); 
    
    TString IncFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_QCD" + PUWeighting + KinWeighting + Selection + ".root";
    IncFile = TFile::Open(IncFileName); 

  }

  TString ThisHistoName;

  for (int fpt = 0; fpt<nFitPtBins; fpt++) {
    
    std::cout << "    Merging " << FitPtBin[fpt] << std::endl;

    float LowFitPtBinEdge = FitPtEdge[fpt];
    float HighFitPtBinEdge = FitPtEdge[fpt+1];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      float MiddlePtBin = (PtRelPtEdge[bpt] + PtRelPtEdge[bpt+1])/2.;
      
      if (MiddlePtBin>LowFitPtBinEdge && MiddlePtBin<HighFitPtBinEdge) {
	
	for (int nb = 0; nb<nPtRelEtaBins; nb++) { 
	  
	  for (int is = 0; is<nSystematics; is++) {
	    
	    ThisHistoName = HistogramName("Observed_" + DataFlag, bpt, nb, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) DataFile->Get(ThisHistoName);
	    
	    for (int tr = 0; tr<nTriggers; tr++)
	      if (AllowedTrigger[bpt][tr]) {

		float LuminosityWeight =  1.;
		
		ThisHistoName = HistogramName("muonPt_" + DataFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuPt = (TH1D*) DataFile->Get(ThisHistoName);

		ThisHistoName = HistogramName("muonDR_" + DataFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuDR = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_" + DataFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_" + DataFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetEta = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_" + DataFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataPVMult = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonPt_" + MCFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_" + MCFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuDR = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_" + MCFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_" + MCFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetEta = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_" + MCFlag, bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDPVMult = (TH1D*) QCDFile->Get(ThisHistoName);
		
		//if (dt==1)
		//if (ThisMuPt->Integral(0, 61)>0.)
		//  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
	
		int tr1 = -1, tr2 = -1;
		for (int trg = 0; trg<nTriggers; trg++)
		  if (trg!=tr)
		    if (AllowedTrigger[bpt][trg]) {
		      
		      if (tr1==-1) tr1 = trg;
		      else if (tr2==-1) tr2 = trg;
		      else 
			std::cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			     << " too many triggers merged in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		      
		    }
		
		if (tr1==-1) {

		  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);

		} else if (tr2==-1) {
		  
		  int T1 = tr, T2 = tr1;
		  if (tr1<tr) { T1 = tr1; T2 = tr; } 
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  ThisHistoName = HistogramName("muonPt_" + DataFlag, bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_" + DataFlag, bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT2 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) - nT1;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_" + MCFlag, bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_" + MCFlag, bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else {
		    
		    std::cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			 << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		    
		  }
		  
		} else {
		  
		  int T1 = tr, T2 = tr1, T3 = tr2;
		  
		  if (tr<tr1 && tr<tr2 && tr1>tr2) { T1 = tr;  T2 = tr2; T3 = tr1; }
		  if (tr>tr1 && tr<tr2 && tr1<tr2) { T1 = tr1; T2 = tr;  T3 = tr2; }
		  //if (tr<tr1 && tr>tr2 && tr1<tr2) { }
		  if (tr>tr1 && tr>tr2 && tr1<tr2) { T1 = tr1; T2 = tr2; T3 = tr;  }
		  //if (tr>tr1 && tr<tr2 && tr1>tr2) { }
		  if (tr<tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr;  T3 = tr1; }
		  if (tr>tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr1; T3 = tr;  }
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  float nT2 = (ThisObserved->GetBinContent(T1+1, T2+1) - nT1)*(1-(ThisObserved->GetBinContent(T3+1, T3+1)/ThisObserved->GetBinContent(T2+1, T2+1))*(TriggerLuminosity[T2]/TriggerLuminosity[T3]));
		  
		  ThisHistoName = HistogramName("muonPt_" + DataFlag, bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_" + DataFlag, bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_" + DataFlag, bpt, nb, T3, is, -1, -1);
		  TH1D *ThisDataMuPt3 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT3 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) + ThisDataMuPt3->Integral(0, 61) - nT1 - nT2;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_" + MCFlag, bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_" + MCFlag, bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else if (tr==T3) {
		    
		    ThisHistoName = HistogramName("muonPt_" + MCFlag, bpt, nb, T3, is, -1, -1);
		    TH1D *ThisQCDMuPt3 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT3/ThisQCDMuPt3->Integral(0, 61);
		    
		  } else {
		    
		    std::cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			   << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		    
		  }
		  
		} 
		
		MuonPt[0][is][fpt][nb]->Add(ThisDataMuPt);
		MuonDR[0][is][fpt][nb]->Add(ThisDataMuDR);
		JetPt[0][is][fpt][nb]->Add(ThisDataJetPt);
		JetEta[0][is][fpt][nb]->Add(ThisDataJetEta);
		PVEvent[0][is][fpt][nb]->Add(ThisDataPVMult);

		if (LuminosityWeight<0. || LuminosityWeight>1000000000.) {

		  std::cout << "PtRelAnalyzer::BuildTemplates: LuminosityWeight " << LuminosityWeight
		       << " for " << TriggerName[tr] << " in " << PtRelPtBin[bpt] << std::endl;

		  LuminosityWeight = 0.;

		}

		ThisQCDMuPt->Scale(LuminosityWeight);
		MuonPt[1][is][fpt][nb]->Add(ThisQCDMuPt);

		ThisQCDMuDR->Scale(LuminosityWeight);
		MuonDR[1][is][fpt][nb]->Add(ThisQCDMuDR);

		ThisQCDJetPt->Scale(LuminosityWeight);
		JetPt[1][is][fpt][nb]->Add(ThisQCDJetPt);

		ThisQCDJetEta->Scale(LuminosityWeight);
		JetEta[1][is][fpt][nb]->Add(ThisQCDJetEta);

		ThisQCDPVMult->Scale(LuminosityWeight);
		PVEvent[1][is][fpt][nb]->Add(ThisQCDPVMult);
		
		if (AddTrackTemplates) {
		 
		  int nScaleBins = ThisQCDJetPt->GetNbinsX();

		  ThisHistoName = HistogramName("jetPt_JetHT", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetPt = (TH1D*) JetFile->Get(ThisHistoName);

                  float JetLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisJetJetPt->Integral(0, nScaleBins);

                  ThisJetJetPt->Scale(JetLuminosityWeight);
		  JetPt[2][is][fpt][nb]->Add(ThisJetJetPt);
		  
		  ThisHistoName = HistogramName("jetEta_JetHT", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetEta = (TH1D*) JetFile->Get(ThisHistoName);

                  ThisJetJetEta->Scale(JetLuminosityWeight);
		  JetEta[2][is][fpt][nb]->Add(ThisJetJetEta);
		  
		  ThisHistoName = HistogramName("jetPt_QCD", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetPt = (TH1D*) IncFile->Get(ThisHistoName);
		  
                  float IncLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisIncJetPt->Integral(0, nScaleBins);

                  ThisIncJetPt->Scale(IncLuminosityWeight);
		  JetPt[3][is][fpt][nb]->Add(ThisIncJetPt);

		  ThisHistoName = HistogramName("jetEta_QCD", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetEta = (TH1D*) IncFile->Get(ThisHistoName);

                  ThisIncJetEta->Scale(IncLuminosityWeight);
		  JetEta[3][is][fpt][nb]->Add(ThisIncJetEta);
		  
		}
		
		for (int tg = 0; tg<nTaggers; tg++) {
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataFlag, bpt, nb, tr, is, tg,  0);
		  TH1D *ThisDataPtRelTag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][0]->Add(ThisDataPtRelTag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataFlag, bpt, nb, tr, is, tg, 10);
		  TH1D *ThisDataPtRelUntag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][4]->Add(ThisDataPtRelUntag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + MCFlag, bpt, nb, tr, is, tg,  1);
		  TH1D *ThisQCDPtRelTagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][0]->Add(ThisQCDPtRelTagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + MCFlag, bpt, nb, tr, is, tg, 11);
		  TH1D *ThisQCDPtRelUntagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][4]->Add(ThisQCDPtRelUntagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + MCFlag, bpt, nb, tr, is, tg,  2);
		  TH1D *ThisQCDPtRelTagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][1]->Add(ThisQCDPtRelTagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + MCFlag, bpt, nb, tr, is, tg, 12);
		  TH1D *ThisQCDPtRelUntagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][5]->Add(ThisQCDPtRelUntagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + MCFlag, bpt, nb, tr, is, tg,  3);
		  TH1D *ThisQCDPtRelTagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][2]->Add(ThisQCDPtRelTagLG);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + MCFlag, bpt, nb, tr, is, tg, 13);
		  TH1D *ThisQCDPtRelUntagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][6]->Add(ThisQCDPtRelUntagLG);
		  
		  if (AddTrackTemplates) {

		    int nScaleBins = ThisQCDPtRelTagLG->GetNbinsX();
		    float TagLuminosityWeight   = ThisQCDPtRelTagLG->Integral(0, nScaleBins);
		    float UntagLuminosityWeight = ThisQCDPtRelUntagLG->Integral(0, nScaleBins);

		    ThisHistoName = HistogramName(TemplateVariable + "_JetHT", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisJetPtRelTagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetTagLuminosityWeight = TagLuminosityWeight/ThisJetPtRelTagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelTagTrk->Scale(JetTagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][3]->Add(ThisJetPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_JetHT", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisJetPtRelUntagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetUntagLuminosityWeight = UntagLuminosityWeight/ThisJetPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelUntagTrk->Scale(JetUntagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][7]->Add(ThisJetPtRelUntagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisIncPtRelTagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncTagLuminosityWeight = TagLuminosityWeight/ThisIncPtRelTagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelTagTrk->Scale(IncTagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][3]->Add(ThisIncPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisIncPtRelUntagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncUntagLuminosityWeight = UntagLuminosityWeight/ThisIncPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelUntagTrk->Scale(IncUntagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][7]->Add(ThisIncPtRelUntagTrk);

		  }
		  
		}

	      }

	    }
		    
	  }
		  
      }
		
    }
    
  }
    
  TString TemplateFlag = "";
  if (DataFlag!="BTagMu") TemplateFlag += DataFlag;
  if (MCFlag!="QCDMu") TemplateFlag += MCFlag;	
  if (AddTrackTemplates) TemplateFlag += "All";	
  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + TemplateFlag + PUWeighting + KinWeighting + Selection + Production + ".root";
  
  SaveTemplates(TemplateFileName, AddTrackTemplates);
  
}

void PtRelAnalyzer::BuildTemplatesEtaSplit(bool AddTrackTemplates) {
  
  BookTemplates();

  std::cout << "Building templates" << std::endl;
  
  TFile *DataFile, *QCDFile, *JetFile, *IncFile;

  TString HistoKinWeighting = KinWeighting;
  HistoKinWeighting.ReplaceAll("EtaSplit", "");

  TString BTagMuConfiguration = PUWeighting;
  BTagMuConfiguration.ReplaceAll("_PU", "_");
  BTagMuConfiguration.ReplaceAll("_PV", "_");
  BTagMuConfiguration.ReplaceAll("_PSV", "_PS");

  TString JetHTConfiguration = PUWeighting + HistoKinWeighting;
  JetHTConfiguration.ReplaceAll("_PU", "_");
  JetHTConfiguration.ReplaceAll("_PSV", "_PS");

  TString DataFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_BTagMu" + BTagMuConfiguration + Selection + ".root";
  DataFile = TFile::Open(DataFileName); 
  
  TString QCDFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_QCDMu" + PUWeighting + HistoKinWeighting + Selection + ".root";
  QCDFile = TFile::Open(QCDFileName); 

  if (AddTrackTemplates) {

    TString JetFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_JetHT" + PUWeighting + HistoKinWeighting + Selection + ".root";
    JetFile = TFile::Open(JetFileName); 
    
    TString IncFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_QCD" + JetHTConfiguration + Selection + ".root";
    IncFile = TFile::Open(IncFileName); 

  }
  
  TString ThisHistoName;

  for (int fpt = 0; fpt<nFitPtBins; fpt++) {
    
    std::cout << "    Merging " << FitPtBin[fpt] << std::endl;

    float LowFitPtBinEdge = FitPtEdge[fpt];
    float HighFitPtBinEdge = FitPtEdge[fpt+1];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      float MiddlePtBin = (PtRelPtEdge[bpt] + PtRelPtEdge[bpt+1])/2.;
      
      if (MiddlePtBin>LowFitPtBinEdge && MiddlePtBin<HighFitPtBinEdge) {
	    
	for (int nb = 1; nb<nPtRelEtaBins; nb++) { 
	      
	  for (int is = 0; is<nSystematics; is++) {

	    ThisHistoName = HistogramName("Observed_BTagMu", bpt, nb, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) DataFile->Get(ThisHistoName);
	    
	    for (int tr = 0; tr<nTriggers; tr++)
	      if (AllowedTrigger[bpt][tr]) {
		
		float LuminosityWeight =  1.;
		
		ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuDR = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetEta = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_BTagMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataPVMult = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuDR = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetEta = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_QCDMu", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDPVMult = (TH1D*) QCDFile->Get(ThisHistoName);
		
		//if (dt==1)
		//if (ThisMuPt->Integral(0, 61)>0.)
		//  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		
		int tr1 = -1, tr2 = -1;
		for (int trg = 0; trg<nTriggers; trg++)
		  if (trg!=tr)
		    if (AllowedTrigger[bpt][trg]) {
		      
		      if (tr1==-1) tr1 = trg;
		      else if (tr2==-1) tr2 = trg;
		      else 
			std::cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			     << " too many triggers merged in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		      
		    }
		
		if (tr1==-1) {
		  
		  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		  
		} else if (tr2==-1) {
		  
		  int T1 = tr, T2 = tr1;
		  if (tr1<tr) { T1 = tr1; T2 = tr; } 
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT2 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) - nT1;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else {
		    
		    std::cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			 << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		    
		  }
		  
		} else {
		  
		  int T1 = tr, T2 = tr1, T3 = tr2;
		  
		  if (tr<tr1 && tr<tr2 && tr1>tr2) { T1 = tr;  T2 = tr2; T3 = tr1; }
		  if (tr>tr1 && tr<tr2 && tr1<tr2) { T1 = tr1; T2 = tr;  T3 = tr2; }
		  //if (tr<tr1 && tr>tr2 && tr1<tr2) { }
		  if (tr>tr1 && tr>tr2 && tr1<tr2) { T1 = tr1; T2 = tr2; T3 = tr;  }
		  //if (tr>tr1 && tr<tr2 && tr1>tr2) { }
		  if (tr<tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr;  T3 = tr1; }
		  if (tr>tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr1; T3 = tr;  }
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  float nT2 = (ThisObserved->GetBinContent(T1+1, T2+1) - nT1)*(1-(ThisObserved->GetBinContent(T3+1, T3+1)/ThisObserved->GetBinContent(T2+1, T2+1))*(TriggerLuminosity[T2]/TriggerLuminosity[T3]));
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_BTagMu", bpt, nb, T3, is, -1, -1);
		  TH1D *ThisDataMuPt3 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT3 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) + ThisDataMuPt3->Integral(0, 61) - nT1 - nT2;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else if (tr==T3) {
		    
		    ThisHistoName = HistogramName("muonPt_QCDMu", bpt, nb, T3, is, -1, -1);
		    TH1D *ThisQCDMuPt3 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT3/ThisQCDMuPt3->Integral(0, 61);
		    
		  } else {
		    
		    std::cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			   << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << std::endl; 
		    
		  }
		  
		} 
		
		MuonPt[0][is][fpt][nb]->Add(ThisDataMuPt);
		MuonDR[0][is][fpt][nb]->Add(ThisDataMuDR);
		JetPt[0][is][fpt][nb]->Add(ThisDataJetPt);
		JetEta[0][is][fpt][nb]->Add(ThisDataJetEta);
		PVEvent[0][is][fpt][nb]->Add(ThisDataPVMult);

		if (LuminosityWeight<0. || LuminosityWeight>1000000000.) {

		  std::cout << "PtRelAnalyzer::BuildTemplates: LuminosityWeight " << LuminosityWeight
		       << " for " << TriggerName[tr] << " in " << PtRelPtBin[bpt] << std::endl;

		  LuminosityWeight = 0.;

		}

		ThisQCDMuPt->Scale(LuminosityWeight);
		MuonPt[1][is][fpt][nb]->Add(ThisQCDMuPt);

		ThisQCDMuDR->Scale(LuminosityWeight);
		MuonDR[1][is][fpt][nb]->Add(ThisQCDMuDR);

		ThisQCDJetPt->Scale(LuminosityWeight);
		JetPt[1][is][fpt][nb]->Add(ThisQCDJetPt);

		ThisQCDJetEta->Scale(LuminosityWeight);
		JetEta[1][is][fpt][nb]->Add(ThisQCDJetEta);

		ThisQCDPVMult->Scale(LuminosityWeight);
		PVEvent[1][is][fpt][nb]->Add(ThisQCDPVMult);
 
		if (AddTrackTemplates) {
		 
		  int nScaleBins = ThisQCDJetPt->GetNbinsX();

		  ThisHistoName = HistogramName("jetPt_JetHT", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetPt = (TH1D*) JetFile->Get(ThisHistoName);

                  float JetLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisJetJetPt->Integral(0, nScaleBins);

                  ThisJetJetPt->Scale(JetLuminosityWeight);
		  JetPt[2][is][fpt][nb]->Add(ThisJetJetPt);
		  
		  ThisHistoName = HistogramName("jetEta_JetHT", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetEta = (TH1D*) JetFile->Get(ThisHistoName);

                  ThisJetJetEta->Scale(JetLuminosityWeight);
		  JetEta[2][is][fpt][nb]->Add(ThisJetJetEta);
		  
		  ThisHistoName = HistogramName("jetPt_QCD", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetPt = (TH1D*) IncFile->Get(ThisHistoName);
		  
                  float IncLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisIncJetPt->Integral(0, nScaleBins);

                  ThisIncJetPt->Scale(IncLuminosityWeight);
		  JetPt[3][is][fpt][nb]->Add(ThisIncJetPt);

		  ThisHistoName = HistogramName("jetEta_QCD", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetEta = (TH1D*) IncFile->Get(ThisHistoName);

                  ThisIncJetEta->Scale(IncLuminosityWeight);
		  JetEta[3][is][fpt][nb]->Add(ThisIncJetEta);
		  
		}

		for (int tg = 0; tg<nTaggers; tg++) {
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_BTagMu", bpt, nb, tr, is, tg,  0);
		  TH1D *ThisDataPtRelTag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][0]->Add(ThisDataPtRelTag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_BTagMu", bpt, nb, tr, is, tg, 10);
		  TH1D *ThisDataPtRelUntag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][4]->Add(ThisDataPtRelUntag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCDMu", bpt, nb, tr, is, tg,  1);
		  TH1D *ThisQCDPtRelTagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][0]->Add(ThisQCDPtRelTagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCDMu", bpt, nb, tr, is, tg, 11);
		  TH1D *ThisQCDPtRelUntagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][4]->Add(ThisQCDPtRelUntagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCDMu", bpt, nb, tr, is, tg,  2);
		  TH1D *ThisQCDPtRelTagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][1]->Add(ThisQCDPtRelTagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCDMu", bpt, nb, tr, is, tg, 12);
		  TH1D *ThisQCDPtRelUntagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][5]->Add(ThisQCDPtRelUntagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCDMu", bpt, nb, tr, is, tg,  3);
		  TH1D *ThisQCDPtRelTagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][2]->Add(ThisQCDPtRelTagLG);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCDMu", bpt, nb, tr, is, tg, 13);
		  TH1D *ThisQCDPtRelUntagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][6]->Add(ThisQCDPtRelUntagLG);

		  if (AddTrackTemplates) {

		    int nScaleBins = ThisQCDPtRelTagLG->GetNbinsX();
		    float TagLuminosityWeight   = ThisQCDPtRelTagLG->Integral(0, nScaleBins);
		    float UntagLuminosityWeight = ThisQCDPtRelUntagLG->Integral(0, nScaleBins);

		    ThisHistoName = HistogramName(TemplateVariable + "_JetHT", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisJetPtRelTagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetTagLuminosityWeight = TagLuminosityWeight/ThisJetPtRelTagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelTagTrk->Scale(JetTagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][3]->Add(ThisJetPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_JetHT", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisJetPtRelUntagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetUntagLuminosityWeight = UntagLuminosityWeight/ThisJetPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelUntagTrk->Scale(JetUntagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][7]->Add(ThisJetPtRelUntagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisIncPtRelTagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncTagLuminosityWeight = TagLuminosityWeight/ThisIncPtRelTagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelTagTrk->Scale(IncTagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][3]->Add(ThisIncPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisIncPtRelUntagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncUntagLuminosityWeight = UntagLuminosityWeight/ThisIncPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelUntagTrk->Scale(IncUntagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][7]->Add(ThisIncPtRelUntagTrk);

		  }

		}

	      }
		      
	  }
	  
	}

      }

    }
	    
    for (int nb = 1; nb<nPtRelEtaBins; nb++) { 
      
      for (int is = 0; is<nSystematics; is++) {
	
	for (int dt = 0; dt<2; dt++) {
	  
	  MuonPt[dt][is][fpt][0]->Add(MuonPt[dt][is][fpt][nb]);
	  MuonDR[dt][is][fpt][0]->Add(MuonDR[dt][is][fpt][nb]);
	  JetPt[dt][is][fpt][0]->Add(JetPt[dt][is][fpt][nb]);
	  JetEta[dt][is][fpt][0]->Add(JetEta[dt][is][fpt][nb]);
	  PVEvent[dt][is][fpt][0]->Add(PVEvent[dt][is][fpt][nb]);
	  
	}
	
	if (AddTrackTemplates) {
	  
	  JetPt[2][is][fpt][0]->Add(JetPt[2][is][fpt][nb]);
	  JetEta[2][is][fpt][0]->Add(JetEta[2][is][fpt][nb]);
	  
	  JetPt[3][is][fpt][0]->Add(JetPt[3][is][fpt][nb]);
	  JetEta[3][is][fpt][0]->Add(JetEta[3][is][fpt][nb]);
	  
	}
	
	for (int tg = 0; tg<nTaggers; tg++) {
	  
	  PtRel[tg][is][fpt][0][0]->Add(PtRel[tg][is][fpt][nb][0]);
	  PtRel[tg][is][fpt][0][4]->Add(PtRel[tg][is][fpt][nb][4]);
	  
	  PtRel[nTaggers+tg][is][fpt][0][0]->Add(PtRel[nTaggers+tg][is][fpt][nb][0]);
	  PtRel[nTaggers+tg][is][fpt][0][4]->Add(PtRel[nTaggers+tg][is][fpt][nb][4]);
	  
	  PtRel[nTaggers+tg][is][fpt][0][1]->Add(PtRel[nTaggers+tg][is][fpt][nb][1]);
	  PtRel[nTaggers+tg][is][fpt][0][5]->Add(PtRel[nTaggers+tg][is][fpt][nb][5]);
	  
	  PtRel[nTaggers+tg][is][fpt][0][2]->Add(PtRel[nTaggers+tg][is][fpt][nb][2]);
	  PtRel[nTaggers+tg][is][fpt][0][6]->Add(PtRel[nTaggers+tg][is][fpt][nb][6]);
  
	  if (AddTrackTemplates) {
	    
	    PtRel[tg][is][fpt][0][3]->Add(PtRel[tg][is][fpt][nb][3]);
	    PtRel[tg][is][fpt][0][7]->Add(PtRel[tg][is][fpt][nb][7]);
	    
	    PtRel[nTaggers+tg][is][fpt][0][3]->Add(PtRel[nTaggers+tg][is][fpt][nb][3]);
	    PtRel[nTaggers+tg][is][fpt][0][7]->Add(PtRel[nTaggers+tg][is][fpt][nb][7]);
	    
	  }
	  
	}
	
      }
      
    }
      
  }

  TString LightTemplates = "";
  if (AddTrackTemplates) LightTemplates = "All";	
  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + KinWeighting + Selection + Production + ".root";
  
  SaveTemplates(TemplateFileName, AddTrackTemplates);
  
}

void PtRelAnalyzer::SaveTemplates(TString TemplateFileName, bool AddTrackTemplates) {

  std::cout << "Saving templates" << std::endl;

  TFile *OutFile = new TFile(TemplateFileName, "recreate");
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<2; dt++) {
	  
	  JetPt[dt][is][fpt][nb]->Write();
	  JetEta[dt][is][fpt][nb]->Write();

	  if (AddTrackTemplates) {

	    JetPt[dt+2][is][fpt][nb]->Write();
	    JetEta[dt+2][is][fpt][nb]->Write();
	    
	  }

	  MuonPt[dt][is][fpt][nb]->Write();
	  MuonDR[dt][is][fpt][nb]->Write();
	  PVEvent[dt][is][fpt][nb]->Write();
	    
	  for (int tg = 0; tg<nTaggers; tg++)
	    for (int tp = 0; tp<2; tp++)
	      for (int fl = 0; fl<4; fl++)
		if ((dt==1 && fl<3) || fl==0 || (fl==3 && AddTrackTemplates))
		  PtRel[dt*nTaggers+tg][is][fpt][nb][4*tp+fl]->Write();
	  
	}
  
  OutFile->Close();

  std::cout << "Exiting" << std::endl;
  
}

void PtRelAnalyzer::ComputeKinematicWeights(TString SystematicFlag, TString DataType, TString LightTemplates) {

  std::cout << "Producing kinematic weights" << std::endl;
  
  if (KinWeighting.Contains("_KinPtBins")) {
	
    TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + Selection + "_BaseProduction.root";
    TFile *TemplateFile = TFile::Open(TemplateFileName); 
    std::cout << TemplateFileName << std::endl;
    bool GoodSystematicFlag = false;

    for (int is = 0; is<nSystematics; is++) 
      if (SystematicFlag==SystematicName[is]) {
	
	for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
	  for (int nb = 0; nb<nPtRelEtaBins; nb++) {

	    std::cout << "  Weights for " << PtRelPtBin[bpt] << " " << PtRelEtaBin[nb] << std::endl;
	    
	    TString ThisRefPtName = HistogramName("jetPt_" + DataType, PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisRefPt = (TH1D*) TemplateFile->Get(ThisRefPtName);
	    	    
	    TString ThisDataPtName = HistogramName("jetPt_BTagMu", PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisDataPt = (TH1D*) TemplateFile->Get(ThisDataPtName);
	    	    
	    float DataEvents = ThisDataPt->Integral();
	    float RefEvents = ThisRefPt->Integral();
	    ThisRefPt->Scale(DataEvents/RefEvents);
	    
	    int PtRebin = 2;
	    if (PtRelPtEdge[bpt]>=120.) PtRebin = 4;
	    if (PtRelPtEdge[bpt]>=160.) PtRebin = 10;
	    if (PtRelPtEdge[bpt]>=320.) PtRebin = 20;

	    ThisDataPt->Rebin(PtRebin);
	    ThisRefPt->Rebin(PtRebin);
	    ThisDataPt->Divide(ThisRefPt);
	    std::cout << "Pt " << DataEvents << " " << RefEvents << std::endl;
	    WriteKinematicWeights(ThisDataPt, "_Pt", bpt, nb, DataType);
	   	    
	    GoodSystematicFlag = true;
	    
	  }
	
      }

  } else if (KinWeighting.Contains("_KinEtaBins")) {
	
    TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + Selection + "_BaseProduction.root";
    TFile *TemplateFile = TFile::Open(TemplateFileName); 
    std::cout << TemplateFileName << std::endl;
    bool GoodSystematicFlag = false;

    for (int is = 0; is<nSystematics; is++) 
      if (SystematicFlag==SystematicName[is]) {
	
	for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
	  for (int nb = 0; nb<nPtRelEtaBins; nb++) {

	    std::cout << "  Weights for " << PtRelPtBin[bpt] << " " << PtRelEtaBin[nb] << std::endl;
	    	    
	    TString ThisRefEtaName = HistogramName("jetEta_" + DataType, PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisRefEta = (TH1D*) TemplateFile->Get(ThisRefEtaName);
	    
	    TString ThisDataEtaName = HistogramName("jetEta_BTagMu", bpt, nb, -1, is, -1, -1);
	    TH1D *ThisDataEta = (TH1D*) TemplateFile->Get(ThisDataEtaName);
	    
	    float DataEvents = ThisDataEta->Integral();
	    float RefEvents = ThisRefEta->Integral();
	    ThisRefEta->Scale(DataEvents/RefEvents);
	    
	    int EtaRebin = 6;
	    ThisDataEta->Rebin(EtaRebin);
	    ThisRefEta->Rebin(EtaRebin);
	    ThisDataEta->Divide(ThisRefEta);
	    std::cout << "Eta " << DataEvents << " " << RefEvents << std::endl;
	    WriteKinematicWeights(ThisDataEta, "_Eta", bpt, nb, DataType);
	   	    
	    GoodSystematicFlag = true;
	    
	  }
	
      }

  } else if (KinWeighting.Contains("_KinMuonPtBins")) {
	
    if (DataType=="QCDMu") {

      TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + Selection + "_BaseProduction.root";
      TFile *TemplateFile = TFile::Open(TemplateFileName); 
      std::cout << TemplateFileName << std::endl;
      bool GoodSystematicFlag = false;
      
      for (int is = 0; is<nSystematics; is++) 
	if (SystematicFlag==SystematicName[is]) {
	  
	  for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
	    for (int nb = 0; nb<nPtRelEtaBins; nb++) {
	      
	      std::cout << "  Weights for " << PtRelPtBin[bpt] << " " << PtRelEtaBin[nb] << std::endl;
	      
	      TString ThisRefMuonPtName = HistogramName("muonPt_" + DataType, PtRelPtBin[bpt], nb, -1, is, -1, -1);
	      TH1D *ThisRefMuonPt = (TH1D*) TemplateFile->Get(ThisRefMuonPtName);
	      
	      TString ThisDataMuonPtName = HistogramName("muonPt_BTagMu", bpt, nb, -1, is, -1, -1);
	      TH1D *ThisDataMuonPt = (TH1D*) TemplateFile->Get(ThisDataMuonPtName);
	      
	      float DataEvents = ThisDataMuonPt->Integral();
	      float RefEvents = ThisRefMuonPt->Integral();
	      ThisRefMuonPt->Scale(DataEvents/RefEvents);
	      
	      int MuonPtRebin = 1;
	      ThisDataMuonPt->Rebin(MuonPtRebin);
	      ThisRefMuonPt->Rebin(MuonPtRebin);
	      ThisDataMuonPt->Divide(ThisRefMuonPt);
	      std::cout << "MuonPt " << DataEvents << " " << RefEvents << std::endl;
	      WriteKinematicWeights(ThisDataMuonPt, "_MuonPt", bpt, nb, DataType);
	      
	      GoodSystematicFlag = true;
	      
	    }
	  
	}

    }

  } else if (KinWeighting.Contains("_KinEtaAfterPtBins")) {
	
    TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + "_KinPtBinsCentral" + Selection + "_BaseProduction.root";
    TFile *TemplateFile = TFile::Open(TemplateFileName); 
    std::cout << TemplateFileName << std::endl;
    bool GoodSystematicFlag = false;

    for (int is = 0; is<nSystematics; is++) 
      if (SystematicFlag==SystematicName[is]) {
	
	for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
	  for (int nb = 0; nb<nPtRelEtaBins; nb++) {

	    std::cout << "  Weights for " << PtRelPtBin[bpt] << " " << PtRelEtaBin[nb] << std::endl;
	    	    
	    TString ThisRefEtaName = HistogramName("jetEta_" + DataType, PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisRefEta = (TH1D*) TemplateFile->Get(ThisRefEtaName);
	    	    
	    TString ThisDataEtaName = HistogramName("jetEta_BTagMu", bpt, nb, -1, is, -1, -1);
	    TH1D *ThisDataEta = (TH1D*) TemplateFile->Get(ThisDataEtaName);
	    	    
	    float DataEvents = ThisDataEta->Integral();
	    float RefEvents = ThisRefEta->Integral();
	    ThisRefEta->Scale(DataEvents/RefEvents);
	    
	    int EtaRebin = 1;
	    ThisDataEta->Rebin(EtaRebin);
	    ThisRefEta->Rebin(EtaRebin);
	    ThisDataEta->Divide(ThisRefEta);
	    std::cout << "Eta " << DataEvents << " " << RefEvents << std::endl;
	    WriteKinematicWeights(ThisDataEta, "_Eta", bpt, nb, DataType);
	   	    
	    GoodSystematicFlag = true;
	    
	  }
	
      }
	
    if (!GoodSystematicFlag) 
      std::cout << "PtRelAnalyzer::ComputeKinematicWeights: SystematicFlag " << SystematicFlag << " not supported" << std::endl;
    
  }
  
}

void PtRelAnalyzer::WriteKinematicWeights(TH1D *&HistoDivide, TString Variable, int PtBin, int EtaBin, TString DataType) {

  std::cout << "Writing kinematic weights" << std::endl;

  ofstream KinWeightsFile; 
  KinWeightsFile.open("./Weights/KinematicWeights/" + TemplateVariable + KinWeighting + Variable + "_" + PtRelPtBin[PtBin] + "_" + PtRelEtaBin[EtaBin] + "_" + DataType + PUWeighting + Selection + ".txt");
  
  int nSteps = 1; float FirstStep = 0., StepSize = 1.;
  if (Variable=="_Pt") { 

    nSteps = PtRelPtEdge[PtBin+1] - PtRelPtEdge[PtBin];
    FirstStep = PtRelPtEdge[PtBin];
    StepSize = 1.;

  } else if (Variable=="_Eta") { 

    nSteps = 48;
    StepSize = 0.1;
    FirstStep = -PtRelEtaEdge[nPtRelEtaBins-1] + StepSize/2.;

  } else if (Variable=="_MuonPt") { 

    nSteps = 60;
    StepSize = 1.;
    FirstStep = 0.;

  } else 
    cout << "PtRelAnalyzer::WriteKinematicWeights: " << Variable << " not supported" << endl;

  for (int st = 0; st<nSteps; st++) {

    float ThisStep = FirstStep + st*StepSize;
    int ThisStepBin = HistoDivide->FindBin(ThisStep);
    
    float ThisWeight = HistoDivide->GetBinContent(ThisStepBin);
    if (ThisWeight<0.) ThisWeight = 1.;
    if (PtBin==nPtRelPtBins-1) ThisWeight = 1.;

    KinWeightsFile << ThisStep << " " << ThisWeight << std::endl; 

  }

  KinWeightsFile.close();

}

bool PtRelAnalyzer::PassTriggerEmulation(int TriggerIdx, int MuonJetIdx) {

  if (TriggerName[TriggerIdx].Contains("_Jet")) return true; 
  else if (TriggerName[TriggerIdx].Contains("_DiJet")) {
    
    int tJet = GetAwayJet("NONE", MuonJetIdx, 0., false);

    if (tJet>=0) 
      if (Jet_pt[tJet]>=PtEmulationTrigger[TriggerIdx]) return true;

  }

  return false;

}

bool PtRelAnalyzer::PassEventTriggerEmulation(int TriggerIdx, TString DataType) {

  if (DataType=="BTagMu" || DataType=="QCDMu") {

    int MuonIdx;
    int MuonJetIdx = GetPFMuonJet(&MuonIdx);
    
    if (MuonJetIdx<0) return false;
    
    if (Jet_pt[MuonJetIdx]<PtEmulationTrigger[TriggerIdx]) return false;
    
    if (TriggerName[TriggerIdx].Contains("_Jet")) return true; 
    else if (TriggerName[TriggerIdx].Contains("_DiJet")) {
      
      int tJet = GetAwayJet("NONE", MuonJetIdx, 0., false);
      
      if (tJet>=0) 
	if (Jet_pt[tJet]>=PtEmulationTrigger[TriggerIdx]) return true;
      
    }

  } else if (DataType=="JetHT" || DataType=="QCD") {
    
    float MaxJetPt = -1.;
    for (int ijet = 0; ijet<nJet; ijet++)
      if (Jet_pt[ijet]>MaxJetPt) MaxJetPt = Jet_pt[ijet];
    
    if (JetTriggerName[TriggerIdx]=="_PFJet40"   && MaxJetPt> 50.) return true;
    if (JetTriggerName[TriggerIdx]=="_PFJet60"   && MaxJetPt> 70.) return true;
    if (JetTriggerName[TriggerIdx]=="_PFJet80"   && MaxJetPt>100.) return true;
    if (JetTriggerName[TriggerIdx]=="_PFJet260"  && MaxJetPt>310.) return true;
    
  }
  
  return false;

}

TH1D *PtRelAnalyzer::GetPtRelTemplate(bool IsMC, TString Flavour, TString Tagger, TString EtaBin, TString PtBin, TString FitOption) {

  TH1D *ThisPtRel; 

  TString SystematicCode = "_Central";

  for (int is = 0; is<nFitSystematics; is++)
   if (FitOption.Contains(FitSystematicName[is]) && !FitSystematicName[is].Contains("_Central"))
      SystematicCode = FitSystematicName[is];

  if (SystematicCode.Contains("_JEUMC")   && !IsMC) SystematicCode = "_Central";
  if (SystematicCode.Contains("_JEUMC")   &&  IsMC) SystematicCode.ReplaceAll("_JEUMC",   "_JEU");
  if (SystematicCode.Contains("_JEUData") && !IsMC) SystematicCode.ReplaceAll("_JEUData", "_JEU");
  if (SystematicCode.Contains("_JEUData") &&  IsMC) SystematicCode = "_Central";
  
  TString HistoName = TemplateVariable + "_BTagMu_" + PtBin + "_" + EtaBin + SystematicCode + "_" + Tagger + "_Tag" + Flavour;
  if (IsMC) HistoName.ReplaceAll("BTagMu", "QCDMu");
  
  if (!FitOption.Contains("_Pretag") || !Tagger.Contains("Anti")) {

    if (Tagger.Contains("Anti")) { HistoName.ReplaceAll("Anti", ""); HistoName.ReplaceAll("_Tag", "_Untag"); }
    
    if (Flavour=="_lg" && FitOption.Contains("_LightTemplatesJetHT")) {
      
      HistoName.ReplaceAll("_QCDMu", "_JetHT"); HistoName.ReplaceAll("_lg", "_trk");     
      
    }

    if (Flavour=="_lg" && Tagger.Contains("Anti") && FitOption.Contains("_IP3DBias")) {

      HistoName.ReplaceAll("_QCDMu", "_JetHT"); HistoName.ReplaceAll("_lg", "_trk");     
      
    }
    
    if (Flavour=="_lg" && FitOption.Contains("_LightTemplatesQCD")) {
      
      HistoName.ReplaceAll("_QCDMu", "_QCD"); HistoName.ReplaceAll("_lg", "_trk");     
      
    }
    
    if (Flavour=="_b" && FitOption.Contains("_bTemplatesBTagMu") && !FitOption.Contains("_bCorr")) {

      HistoName.ReplaceAll("_QCDMu", "_BTagMu"); 
      HistoName.ReplaceAll("_b", ""); 
      HistoName.ReplaceAll("_Untag", "_Tag");
      HistoName.ReplaceAll("CSVv2L", "CSVv2T");
      HistoName.ReplaceAll("CSVv2M", "CSVv2T");

    }
   
    ThisPtRel = (TH1D*) PtRelTemplateFile->Get(HistoName);
    ThisPtRel->SetDirectory(0); 
    ThisPtRel->Rebin(TemplateRebinning);
   
    if (FitOption.Contains("_LightTemplatesRatio") && !FitOption.Contains("_lCorr")) {
       
      if (Flavour=="_lg" || FitOption.Contains("_LightTemplatesRatioB")) {
	
	TString HistoQCDName = HistoName; HistoQCDName.ReplaceAll("_QCDMu", "_QCD");  
        HistoQCDName.ReplaceAll(Flavour, "_trk");
	TH1D *ThisQCD = (TH1D*) PtRelTemplateFile->Get(HistoQCDName);
	ThisQCD->SetDirectory(0); 
	ThisQCD->Rebin(TemplateRebinning);
	
	TString HistoJetName = HistoName; HistoJetName.ReplaceAll("_QCDMu", "_JetHT");  
	HistoJetName.ReplaceAll(Flavour, "_trk");
	TH1D *ThisJet = (TH1D*) PtRelTemplateFile->Get(HistoJetName);
	ThisJet->SetDirectory(0); 
	ThisJet->Rebin(TemplateRebinning);
	if (ThisJet->Integral()==0) std::cout << "??? " << PtBin << std::endl;
	ThisPtRel->Divide(ThisQCD); ThisPtRel->Multiply(ThisJet); 
	
      }
      
    }
    
    if (Flavour=="_b" && FitOption.Contains("_bTemplatesRatio") && !FitOption.Contains("_bCorr")) {
      
      TString HistoQCDName = HistoName;
      HistoQCDName.ReplaceAll("_Untag", "_Tag");
      HistoQCDName.ReplaceAll("CSVv2L", "CSVv2T");
      HistoQCDName.ReplaceAll("CSVv2M", "CSVv2T");
      HistoQCDName.ReplaceAll("JPL", "JPT");
      HistoQCDName.ReplaceAll("JPM", "JPT");
      
      TH1D *ThisQCD = (TH1D*) PtRelTemplateFile->Get(HistoQCDName); 
      ThisQCD->SetDirectory(0); 
      ThisQCD->Rebin(TemplateRebinning);
      //if (HistoName.Contains("JP") && HistoName.Contains("Pt3050")) {
      //int BinTail = ThisQCD->FindBin(2.99);
      //ThisQCD->SetBinContent(BinTail, 4.);
      //}

      TString HistoDataName = HistoName; 
      HistoDataName.ReplaceAll("_QCDMu", "_BTagMu"); 
      HistoDataName.ReplaceAll("_b", ""); 
      HistoDataName.ReplaceAll("_Untag", "_Tag");
      HistoDataName.ReplaceAll("CSVv2L", "CSVv2T");
      HistoDataName.ReplaceAll("CSVv2M", "CSVv2T");
      HistoDataName.ReplaceAll("JPL", "JPT");
      HistoDataName.ReplaceAll("JPM", "JPT");
      TH1D *ThisData = (TH1D*) PtRelTemplateFile->Get(HistoDataName);
      ThisData->SetDirectory(0); 
      ThisData->Rebin(TemplateRebinning);

      if (FitOption.Contains("_bTemplatesRatioCorr")) {
	
	TString HistoQCDNameC = HistoQCDName;
	HistoQCDNameC.ReplaceAll("_b", "_c");
	TH1D *ThisQCD_C = (TH1D*) PtRelTemplateFile->Get(HistoQCDNameC);
	ThisQCD_C->SetDirectory(0); 
	ThisQCD_C->Rebin(TemplateRebinning);

	ThisData->Add(ThisQCD_C, -1);
       
      }

      float NormQCD = ThisQCD->Integral();
      float NormData = ThisData->Integral();
      ThisQCD->Scale(NormData/NormQCD);
      float NormPtRel = ThisPtRel->Integral();
      //ThisData->Divide(ThisQCD);
      //int BinTail = ThisData->FindBin(3.7);
      //for (int tb = BinTail; tb<=ThisData->GetNbinsX(); tb++)
	//if (ThisData->GetBinContent(tb)>500.) 
        //ThisData->SetBinContent(tb, 1.); 
      ThisPtRel->Divide(ThisQCD); 
      ThisPtRel->Multiply(ThisData); 
      //int BinTail = ThisPtRel->FindBin(3.5);
      //float NormPtRelCorr = ThisPtRel->Integral();
      //for (int tb = BinTail; tb<=ThisPtRel->GetNbinsX(); tb++)
      //if (ThisPtRel->GetBinContent(tb)/NormPtRelCorr>0.1) {
      //  float PrevBin = ThisPtRel->GetBinContent(tb-1);
	  //ThisPtRel->SetBinContent(tb, 0.9*PrevBin); 
	  //ThisPtRel->SetBinError(tb, 0.); 
      //}

      float NormPtRelCorr = ThisPtRel->Integral();
      ThisPtRel->Scale(NormPtRel/NormPtRelCorr);

    }
    
    if (Flavour=="_b" && FitOption.Contains("_bTempRatio") && !FitOption.Contains("_bCorr")) {

      TString HistoQCDName = HistoName;
      HistoQCDName.ReplaceAll("_Untag", "_Tag");
      if (HistoName.Contains("CSVv2")) {
	HistoQCDName.ReplaceAll("CSVv2L", "JPT");
	HistoQCDName.ReplaceAll("CSVv2M", "JPT");
	HistoQCDName.ReplaceAll("CSVv2T", "JPT");
      } else if (HistoName.Contains("JP")) {
	HistoQCDName.ReplaceAll("JPL", "CSVv2T");
	HistoQCDName.ReplaceAll("JPM", "CSVv2T");
	HistoQCDName.ReplaceAll("JPT", "CSVv2T");
      }
      TH1D *ThisQCD = (TH1D*) PtRelTemplateFile->Get(HistoQCDName); 
      ThisQCD->SetDirectory(0); 
      ThisQCD->Rebin(TemplateRebinning);
      if (HistoName.Contains("CSVv2") && HistoName.Contains("Pt3050")) {
	int BinTail = ThisQCD->FindBin(2.99);
	ThisQCD->SetBinContent(BinTail, 4.);
      }

      TString HistoDataName = HistoName; 
      HistoDataName.ReplaceAll("_QCDMu", "_BTagMu"); 
      HistoDataName.ReplaceAll("_b", ""); 
      HistoDataName.ReplaceAll("_Untag", "_Tag");
      if (HistoName.Contains("CSVv2")) {
	HistoDataName.ReplaceAll("CSVv2L", "JPT");
	HistoDataName.ReplaceAll("CSVv2M", "JPT");
	HistoDataName.ReplaceAll("CSVv2T", "JPT");
      } else if (HistoName.Contains("JP")) {
	HistoDataName.ReplaceAll("JPL", "CSVv2T");
	HistoDataName.ReplaceAll("JPM", "CSVv2T");
	HistoDataName.ReplaceAll("JPT", "CSVv2T");
      }
      TH1D *ThisData = (TH1D*) PtRelTemplateFile->Get(HistoDataName);
      ThisData->SetDirectory(0); 
      ThisData->Rebin(TemplateRebinning);
      
      if (FitOption.Contains("_bTempRatioCorr")) {
	
	TString HistoQCDNameC = HistoQCDName;
	HistoQCDNameC.ReplaceAll("_b", "_c");
	TH1D *ThisQCD_C = (TH1D*) PtRelTemplateFile->Get(HistoQCDNameC);
	ThisQCD_C->SetDirectory(0); 
	ThisQCD_C->Rebin(TemplateRebinning);

	ThisData->Add(ThisQCD_C, -1);
       
      }

      float NormQCD = ThisQCD->Integral();
      float NormData = ThisData->Integral();
      ThisQCD->Scale(NormData/NormQCD);
      float NormPtRel = ThisPtRel->Integral();
      //ThisData->Divide(ThisQCD);
      //int BinTail = ThisData->FindBin(3.7);
      //for (int tb = BinTail; tb<=ThisData->GetNbinsX(); tb++)
	//if (ThisData->GetBinContent(tb)>500.) 
        //ThisData->SetBinContent(tb, 1.); 
      ThisPtRel->Divide(ThisQCD); 
      ThisPtRel->Multiply(ThisData); 
      //int BinTail = ThisPtRel->FindBin(3.5);
      //float NormPtRelCorr = ThisPtRel->Integral();
      //for (int tb = BinTail; tb<=ThisPtRel->GetNbinsX(); tb++)
      //if (ThisPtRel->GetBinContent(tb)/NormPtRelCorr>0.1) {
      //  float PrevBin = ThisPtRel->GetBinContent(tb-1);
	  //ThisPtRel->SetBinContent(tb, 0.9*PrevBin); 
	  //ThisPtRel->SetBinError(tb, 0.); 
      //}

      float NormPtRelCorr = ThisPtRel->Integral();
      ThisPtRel->Scale(NormPtRel/NormPtRelCorr);
     
    }
    
    if (Flavour=="_b" && FitOption.Contains("_bTempTightRatio") && !FitOption.Contains("_bCorr")) {

      TString HistoQCDName = HistoName;
      HistoQCDName.ReplaceAll("_Untag", "_Tag");
      if (HistoName.Contains("CSVv2")) {
	HistoQCDName.ReplaceAll("CSVv2L", "JPT");
	HistoQCDName.ReplaceAll("CSVv2M", "JPT");
	HistoQCDName.ReplaceAll("CSVv2T", "JPT");
      } else if (HistoName.Contains("JP")) {
	HistoQCDName.ReplaceAll("JPL", "CSVv2T");
	HistoQCDName.ReplaceAll("JPM", "CSVv2T");
	HistoQCDName.ReplaceAll("JPT", "CSVv2T");
      }
      HistoQCDName.ReplaceAll("_Central", "_JBPT");
      TH1D *ThisQCD = (TH1D*) PtRelTemplateFile->Get(HistoQCDName); 
      ThisQCD->SetDirectory(0); 
      ThisQCD->Rebin(TemplateRebinning);
      if (HistoName.Contains("CSVv2") && HistoName.Contains("Pt3050")) {
	int BinTail = ThisQCD->FindBin(2.99);
	ThisQCD->SetBinContent(BinTail, 4.);
      }

      TString HistoDataName = HistoName; 
      HistoDataName.ReplaceAll("_QCDMu", "_BTagMu"); 
      HistoDataName.ReplaceAll("_b", ""); 
      HistoDataName.ReplaceAll("_Untag", "_Tag");
      if (HistoName.Contains("CSVv2")) {
	HistoDataName.ReplaceAll("CSVv2L", "JPT");
	HistoDataName.ReplaceAll("CSVv2M", "JPT");
	HistoDataName.ReplaceAll("CSVv2T", "JPT");
      } else if (HistoName.Contains("JP")) {
	HistoDataName.ReplaceAll("JPL", "CSVv2T");
	HistoDataName.ReplaceAll("JPM", "CSVv2T");
	HistoDataName.ReplaceAll("JPT", "CSVv2T");
      }
      HistoDataName.ReplaceAll("_Central", "_JBPT");
      TH1D *ThisData = (TH1D*) PtRelTemplateFile->Get(HistoDataName);
      ThisData->SetDirectory(0); 
      ThisData->Rebin(TemplateRebinning);
      
      if (FitOption.Contains("_bTempTightRatioCorr")) {
	
	TString HistoQCDNameC = HistoQCDName;
	HistoQCDNameC.ReplaceAll("_b", "_c");
	TH1D *ThisQCD_C = (TH1D*) PtRelTemplateFile->Get(HistoQCDNameC);
	ThisQCD_C->SetDirectory(0); 
	ThisQCD_C->Rebin(TemplateRebinning);

	ThisData->Add(ThisQCD_C, -1);
       
      }

      float NormQCD = ThisQCD->Integral();
      float NormData = ThisData->Integral();
      ThisQCD->Scale(NormData/NormQCD);
      float NormPtRel = ThisPtRel->Integral();
      //ThisData->Divide(ThisQCD);
      //int BinTail = ThisData->FindBin(3.7);
      //for (int tb = BinTail; tb<=ThisData->GetNbinsX(); tb++)
	//if (ThisData->GetBinContent(tb)>500.) 
        //ThisData->SetBinContent(tb, 1.); 
      
      ThisPtRel->Divide(ThisQCD); 
      ThisPtRel->Multiply(ThisData); 
      //int BinTail = ThisPtRel->FindBin(3.5);
      //float NormPtRelCorr = ThisPtRel->Integral();
      //for (int tb = BinTail; tb<=ThisPtRel->GetNbinsX(); tb++)
      //if (ThisPtRel->GetBinContent(tb)/NormPtRelCorr>0.1) {
      //  float PrevBin = ThisPtRel->GetBinContent(tb-1);
	  //ThisPtRel->SetBinContent(tb, 0.9*PrevBin); 
	  //ThisPtRel->SetBinError(tb, 0.); 
      //}

      float NormPtRelCorr = ThisPtRel->Integral();
      ThisPtRel->Scale(NormPtRel/NormPtRelCorr);
     
    }

  } else {

    HistoName.ReplaceAll("Anti", ""); 
    TString HistoName1 = HistoName, HistoName2 = HistoName;
    HistoName2.ReplaceAll("_Tag", "_Untag"); 

    if (Flavour=="_lg" && FitOption.Contains("_LightTemplatesJetHT")) {
      
      HistoName1.ReplaceAll("_QCDMu", "_JetHT"); HistoName1.ReplaceAll("_lg", "_trk");  
      HistoName2.ReplaceAll("_QCDMu", "_JetHT"); HistoName2.ReplaceAll("_lg", "_trk");     
      
    }

    if (Flavour=="_lg" && FitOption.Contains("_LightTemplatesQCD")) {
      
      HistoName1.ReplaceAll("_QCDMu", "_QCD"); HistoName1.ReplaceAll("_lg", "_trk");     
      HistoName2.ReplaceAll("_QCDMu", "_QCD"); HistoName2.ReplaceAll("_lg", "_trk");     
       
    }

    if (FitOption.Contains("_PT2")) {

      TH1D *ThisFirstPtRel; ThisFirstPtRel = (TH1D*) PtRelTemplateFile->Get(HistoName1); 
      ThisFirstPtRel->SetDirectory(0); 
      ThisFirstPtRel->Rebin(TemplateRebinning);
      //ThisPtRel = ReshapeHistogram(ThisFirstPtRel, "JP");
      
    } else  { 
      ThisPtRel = (TH1D*) PtRelTemplateFile->Get(HistoName1); 
      ThisPtRel->SetDirectory(0); 
      ThisPtRel->Rebin(TemplateRebinning);
    }
    TH1D *ThisOtherPtRel; ThisOtherPtRel =  (TH1D*) PtRelTemplateFile->Get(HistoName2);  
    ThisOtherPtRel->SetDirectory(0); 
    ThisOtherPtRel->Rebin(TemplateRebinning);
    ThisPtRel->Add( ThisOtherPtRel );

    if (FitOption.Contains("_LightTemplatesRatio")) {
      
      if (Flavour=="_lg" || FitOption.Contains("_LightTemplatesRatioB")) {
	
	TString HistoQCDName1 = HistoName1; HistoQCDName1.ReplaceAll("_QCDMu", "_QCD");   
	HistoQCDName1.ReplaceAll(Flavour, "_trk");
	TH1D *ThisQCD1 = (TH1D*) PtRelTemplateFile->Get(HistoQCDName1);
	ThisQCD1->SetDirectory(0); 
	ThisQCD1->Rebin(TemplateRebinning);

	TString HistoQCDName2 = HistoName2; HistoQCDName2.ReplaceAll("_QCDMu", "_QCD");    
	HistoQCDName2.ReplaceAll(Flavour, "_trk");
	TH1D *ThisQCD2 = (TH1D*) PtRelTemplateFile->Get(HistoQCDName2);
	ThisQCD2->SetDirectory(0); 
	ThisQCD2->Rebin(TemplateRebinning);

	ThisQCD1->Add( ThisQCD2 );
	
	TString HistoJetName1 = HistoName1; HistoJetName1.ReplaceAll("_QCDMu", "_JetHT");  
	HistoJetName1.ReplaceAll(Flavour, "_trk");
	TH1D *ThisJet1 = (TH1D*) PtRelTemplateFile->Get(HistoJetName1);
	ThisJet1->SetDirectory(0); 
	ThisJet1->Rebin(TemplateRebinning);
	if (ThisJet1->Integral()==0) std::cout << "??? " << PtBin << std::endl;
	
	TString HistoJetName2 = HistoName2; HistoJetName2.ReplaceAll("_QCDMu", "_JetHT");  
	HistoJetName2.ReplaceAll(Flavour, "_trk");
	TH1D *ThisJet2 = (TH1D*) PtRelTemplateFile->Get(HistoJetName2);
	ThisJet2->SetDirectory(0); 
	ThisJet2->Rebin(TemplateRebinning);
	if (ThisJet2->Integral()==0) std::cout << "??? " << PtBin << std::endl;

	ThisJet1->Add( ThisJet2 );

	ThisPtRel->Divide(ThisQCD1); ThisPtRel->Multiply(ThisJet1); 
	
      }

    }

  }

  if (Flavour=="_lg" && FitOption.Contains("_LightFraction")) ThisPtRel->Scale(1.27);

  if (Flavour=="_lg" && FitOption.Contains("_LightCharmRatio")) {

    if (FitOption.Contains("_LightCharmRatio-m30")) ThisPtRel->Scale(0.7);
    if (FitOption.Contains("_LightCharmRatio-p30")) ThisPtRel->Scale(1.3);
    if (FitOption.Contains("_LightCharmRatio-m20")) ThisPtRel->Scale(0.8);
    if (FitOption.Contains("_LightCharmRatio-p20")) ThisPtRel->Scale(1.2);

  }
  
  if (!ThisPtRel) {
    
    std::cout << "Warning: ptrel " << HistoName << " template not found in histogram file!" << std::endl;
    TH1D *ThisTemplate = new TH1D("ThisTemplate", "", nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
    ThisTemplate->Rebin(TemplateRebinning);
    
    return ThisTemplate;
    
  } else {
  
    TH1D *ThisTemplate = (TH1D*)ThisPtRel->Clone();
    delete ThisPtRel;
    return ThisTemplate;
	      
  }

}

struct PtRelFitResult PtRelAnalyzer::PtRelFit(TH1D* DataTemplate, TH1D* bTemplate, TH1D* cTemplate, TH1D* lgTemplate, TString FitOption) {

  struct PtRelFitResult FitContent;
  
  double bIntegral = bTemplate->Integral();
  double cIntegral = cTemplate->Integral();
  double lgIntegral = lgTemplate->Integral();
  
  double bEntries = bTemplate->GetEntries();
  double cEntries = cTemplate->GetEntries();
  double lgEntries = lgTemplate->GetEntries();

  double DataContent = DataTemplate->Integral();

  TObjArray *templates; int nFitComponents;
  if (FitOption.Contains("_TCF")) {

    nFitComponents = 3;

    bTemplate->Scale(bEntries/bIntegral);
    cTemplate->Scale(cEntries/cIntegral);
    lgTemplate->Scale(lgEntries/lgIntegral);

    templates = new TObjArray(nFitComponents);
    templates->Add(bTemplate);
    templates->Add(cTemplate);
    templates->Add(lgTemplate);

  } else {

    nFitComponents = 2;

    lgTemplate->Add(cTemplate);

    bTemplate->Scale(bEntries/bIntegral);
    lgTemplate->Scale((lgEntries+cEntries)/(lgIntegral+cIntegral));

    templates = new TObjArray(nFitComponents);
    templates->Add(bTemplate);
    templates->Add(lgTemplate);

  } 

  TFractionFitter* PtRelFitter = new TFractionFitter(DataTemplate, templates);

  TVirtualFitter::SetErrorDef(0.5);

  PtRelFitter->SetRangeX(1, DataTemplate->FindBin(UpperEdgeForTemp));

  double LowerLimit[3] = {0.001, 0.001, 0.001};
  double UpperLimit[3] = {DataContent, DataContent, DataContent};

  for (int fc = 0; fc<nFitComponents; fc++) 
    PtRelFitter->Constrain(fc, LowerLimit[fc],  UpperLimit[fc]);
   
  double FittedFraction[3] = {-1., -1., -1.};
  double FittedFractionError[3] = {999., 999., 999.};

  //1 sigma : 1. for chi2 and 0.5 for chi2
  int status = PtRelFitter->Fit();

  FitContent.Status = status;
 
  if( (!status)==0) {

    std::cout << "Fit failed for data " << DataTemplate->GetName() << "\n";
      
  } else {
    
    // need to check the effec od this ErrorAnalysis
    FitContent.Chi2N = PtRelFitter->GetChisquare()/PtRelFitter->GetNDF(); 
    std::cout << "Fit done, Chi2 = " << PtRelFitter->GetChisquare()/PtRelFitter->GetNDF() << std::endl;;
    PtRelFitter->ErrorAnalysis(0.5);

    for (int fc = 0; fc<nFitComponents; fc++) 
      PtRelFitter->GetResult(fc, FittedFraction[fc],  FittedFractionError[fc]);
    
    bool LimitIsReached = false;
    for (int fc = 0; fc<nFitComponents; fc++) 
      if (FittedFraction[fc]<=LowerLimit[fc] ||  FittedFraction[fc]>=UpperLimit[fc]) LimitIsReached = true;
    
    if (LimitIsReached) { 
   
       std::cout << "  PtRelAnalyzer::PtRelFit: limit reached !!!!!" << std::endl;
       
    } else {
    
      // What a hell was that???
//       double c1_1 = clg_frac/((clg_frac+b_frac)*(clg_frac+b_frac));
//       double c2_1 = -b_frac/((clg_frac+b_frac)*(clg_frac+b_frac));
            
//       //std::cout << "  clg_frac/((clg_frac+b_frac)*(clg_frac+b_frac)): " << c1_1  << std::endl 
//       //   << "  -b_frac/((clg_frac+b_frac)*(clg_frac+b_frac)): " << c2_1 << std::endl;

//       b_newError = sqrt(c1_1*c1_1*fitter->GetCovarianceMatrixElement(0,0) + 
// 			c2_1*c2_1*fitter->GetCovarianceMatrixElement(1,1) + 
// 			2*c1_1*c2_1*fitter->GetCovarianceMatrixElement(0,1));
      
//       clg_newError = b_newError;

    }

  }
  
  bool FitIsFailed = false; double TotalFraction = 0;
  for (int fc = 0; fc<nFitComponents; fc++) 
    if (FittedFraction[fc]>0.) TotalFraction += FittedFraction[fc]; // Should be 1 anyway ...
    else FitIsFailed = true;

  if (!FitIsFailed && FittedFraction[0]>=0.01) {

      // What a hell was that???
    //bContent.bFrac = b_frac/((clg_frac+b_frac)*(clg_frac+b_frac));
    //bContent.clgFrac = clg_frac/((clg_frac+b_frac)*(clg_frac+b_frac));

    FitContent.bFraction = FittedFraction[0]/TotalFraction;
    FitContent.lightFraction = FittedFraction[1]/TotalFraction;
    if (FitOption.Contains("_TCF")) FitContent.cFraction = FittedFraction[2]/TotalFraction; 
    
    //bContent.bFrac = b_frac; 
    //bContent.clgFrac = clg_frac; 

    // New errors
    //bContent.bError = b_newError;
    //bContent.clgError = clg_newError;

    // Default TFractionFitter Error
    FitContent.bFractionError  = FittedFractionError[0]/TotalFraction;
    FitContent.lightFractionError  = FittedFractionError[1]/TotalFraction; 
    if (FitOption.Contains("_TCF")) FitContent.cFractionError = FittedFractionError[2]/TotalFraction; 

  } else {

    double mcIntegral = bIntegral + cIntegral + lgIntegral;

    FitContent.bFraction = bIntegral/mcIntegral;
    FitContent.lightFraction = (lgIntegral+cIntegral)/mcIntegral;
    if (FitOption.Contains("_TCF")) {
      FitContent.lightFraction = lgIntegral/mcIntegral;
      FitContent.cFraction = cIntegral/mcIntegral;
    }

    FitContent.bFractionError  = 0.3;
    FitContent.lightFractionError = 0.3;
    FitContent.cFractionError  = 0.3;

  }

  if (FitOption.Contains("_TCF")) {

    lgTemplate->Scale(DataContent*FitContent.lightFraction/lgEntries);

  } else {

    FitContent.cFraction = 0.;
    lgTemplate->Scale(DataContent*FitContent.lightFraction/(lgEntries+cEntries));

  } 

  bTemplate->Scale(DataContent*FitContent.bFraction/bEntries);
  cTemplate->Scale(DataContent*FitContent.cFraction/cEntries);
 
  return FitContent;
  
}

double PtRelAnalyzer::BinomialError(double a, double b) {
  double e=a/b;
  double error=sqrt(e*(1-e)/b);
  if(e==1.)error=sqrt(a)/b;
  if(e==0.)error=0;
  return error;
}

double PtRelAnalyzer::EfficiencyStatisticalError(double a, double b) {
  double e_a=sqrt(a);
  double e_b=sqrt(b);
  double error = fabs(1/(a+b)-a/((a+b)*(a+b)))*e_a+(a/((a+b)*(a+b)))*e_b;
  return error;
}

double PtRelAnalyzer::PtRelEfficiencyError(double a, double b, double fa, double fb, double e_fa, double e_fb) {
  double afa   = a*fa;
  double bfb   = b*fb;
  double e_afa = sqrt(a*e_fa*a*e_fa + fa*fa*a);
  double e_bfb = sqrt(b*e_fb*b*e_fb + fb*fb*b);
  double error_afa = (1/(afa+bfb)-afa/((afa+bfb)*(afa+bfb)))*e_afa;
  double error_bfb = (afa/((afa+bfb)*(afa+bfb)))*e_bfb;
  double error = sqrt(error_afa*error_afa + error_bfb*error_bfb);
  //double error = fabs(1/(afa+bfb)-afa/((afa+bfb)*(afa+bfb)))*e_afa+(afa/((afa+bfb)*(afa+bfb)))*e_bfb;
  //double error = e_afa/(b*fb) + (a*fa)/(b*fb*b*fb)*e_bfb;
  return error;
}

double PtRelAnalyzer::RatioPropagationError(double a, double b, double e_a, double e_b) {

  double error_a = e_a/b;
  double error_b = a/(b*b)*e_b;
  double error = sqrt(error_a*error_a + error_b*error_b);

  return error;

}

struct ScaleFactorResult PtRelAnalyzer::PtRelScaleFactor(TString Tagger, TString EtaBin, TString PtBin, TString SystematicFlag, TString FitOption, TString PrintPlot) {

  PtRelFitResult FitContent[2];

  double DataContent[2], MonteCarloBottomContent[2], MonteCarloBottomEntries[2];
  
  for (int t = 0; t<2; t++) {
    
    TString ThisTagger = Tagger;
    if (t==0) ThisTagger.Insert(0, "Anti");
    
    TH1D *DataTemplate = GetPtRelTemplate(false, "",    ThisTagger, EtaBin, PtBin, SystematicFlag + FitOption);
    TH1D *bTemplate    = GetPtRelTemplate(true,  "_b",  ThisTagger, EtaBin, PtBin, SystematicFlag + FitOption);
    TH1D *cTemplate    = GetPtRelTemplate(true,  "_c",  ThisTagger, EtaBin, PtBin, SystematicFlag + FitOption);
    TH1D *lgTemplate   = GetPtRelTemplate(true,  "_lg", ThisTagger, EtaBin, PtBin, SystematicFlag + FitOption);

    DataContent[t] = DataTemplate->Integral();
    MonteCarloBottomContent[t] = bTemplate->Integral();
    MonteCarloBottomEntries[t] = bTemplate->GetEntries();

    FitContent[t] = PtRelFit(DataTemplate, bTemplate, cTemplate, lgTemplate, FitOption);

    TCanvas *CanvasPtRelFit = new TCanvas();// á la Kirill, it was --> "PtRelFit", "", 900, 600); 
    CanvasPtRelFit->Divide(1, 1);
        
    //CanvasPtRelFit->cd(1);

    TPad *PadPtRelFit = (TPad*)CanvasPtRelFit->GetPad(1);

    PadPtRelFit->cd();

    // Run2015B Setting
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    CanvasPtRelFit->Range(0,0,1,1);
    CanvasPtRelFit->SetFillColor(10);
    CanvasPtRelFit->SetBorderMode(0);
    CanvasPtRelFit->SetBorderSize(2);
    CanvasPtRelFit->SetTickx(1);
    CanvasPtRelFit->SetTicky(1);
    CanvasPtRelFit->SetLeftMargin(0.16);
    CanvasPtRelFit->SetRightMargin(0.02);
    CanvasPtRelFit->SetTopMargin(0.05);
    CanvasPtRelFit->SetBottomMargin(0.13);
    CanvasPtRelFit->SetFrameFillColor(0);
    CanvasPtRelFit->SetFrameFillStyle(0);
    CanvasPtRelFit->SetFrameBorderMode(0);

    PadPtRelFit->SetFillColor(0);
    PadPtRelFit->SetBorderMode(0);
    PadPtRelFit->SetBorderSize(2);
    //PadPtRelFit->SetGridy();
    //PadPtRelFit->SetLogx();
    PadPtRelFit->SetTickx(1);
    PadPtRelFit->SetTicky(1);
    PadPtRelFit->SetLeftMargin(0.16);
    PadPtRelFit->SetRightMargin(0.02);
    //PadPtRelFit->SetTopMargin(0.05);
    //PadPtRelFit->SetBottomMargin(0.31);
    PadPtRelFit->SetTopMargin(0.065);
    PadPtRelFit->SetBottomMargin(0.13);
    PadPtRelFit->SetFrameFillStyle(0);
    PadPtRelFit->SetFrameBorderMode(0);
    PadPtRelFit->SetFrameFillStyle(0);
    PadPtRelFit->SetFrameBorderMode(0);
    PadPtRelFit->Draw();
    // End Run2015B Setting 
    
    DataTemplate->GetXaxis()->SetLabelFont(42);
    DataTemplate->GetYaxis()->SetLabelFont(42);
    DataTemplate->GetXaxis()->SetTitleFont(42);
    DataTemplate->GetYaxis()->SetTitleFont(42);

    DataTemplate->SetTitle("");
    DataTemplate->SetXTitle(TemplateXTitle);
    DataTemplate->SetYTitle(TemplateYTitle);

    DataTemplate->GetYaxis()->SetTitleSize(0.06);
    DataTemplate->GetYaxis()->SetLabelSize(0.05);
    DataTemplate->GetXaxis()->SetTitleSize(0.06);
    DataTemplate->GetXaxis()->SetLabelSize(0.05);
    DataTemplate->GetXaxis()->SetTitleOffset(0.95);
    DataTemplate->GetYaxis()->SetTitleOffset(1.25);

    DataTemplate->SetMarkerStyle(20);
    float HistoMaximum = DataTemplate->GetMaximum();
    DataTemplate->SetMaximum(1.4*HistoMaximum);
    DataTemplate->Draw("pe1x0");

    TH1D *tTemplate; tTemplate = (TH1D*) bTemplate->Clone(); 
    tTemplate->Add(cTemplate); 
    tTemplate->Add(lgTemplate);

    //tTemplate->SetFillColor(kBlue);
    tTemplate->SetFillColor(kRed); 
    //tTemplate->SetLineStyle(2); 
    tTemplate->SetLineWidth(2);
    tTemplate->Draw("histsame");

    TH1D *hTemplate; hTemplate = (TH1D*) lgTemplate->Clone(); hTemplate->Add(cTemplate);   
    hTemplate->SetFillColor(kGreen);
    //hTemplate->SetLineStyle(2); 
    hTemplate->SetLineWidth(2);
    if (FitOption.Contains("_TCF")) hTemplate->Draw("histsame");
    
    //bTemplate->SetFillColor(kRed); 
    lgTemplate->SetFillColor(kBlue); 
    lgTemplate->SetLineWidth(2);
    lgTemplate->Draw("histsame");

    DataTemplate->Draw("pe1x0same");

    TString TitleHeader = "#splitline{ThisTagger}{tagged jets}";
    TitleHeader.ReplaceAll("ThisTagger", ThisTagger);
    if (PtBin=="Pt3050") TitleHeader.ReplaceAll("jets", "jets (30<p_{T}<50 GeV)");
    if (PtBin=="Pt5070") TitleHeader.ReplaceAll("jets", "jets (50<p_{T}<70 GeV)");
    if (PtBin=="Pt70100") TitleHeader.ReplaceAll("jets", "jets (70<p_{T}<100 GeV)");
    if (PtBin=="Pt100140") TitleHeader.ReplaceAll("jets", "jets (100<p_{T}<140 GeV)");
    if (ThisTagger.Contains("Anti")) {
      TitleHeader.ReplaceAll("Anti", "");
      TitleHeader.ReplaceAll("tagged", "vetoed");
    }

    TLegend  *leg_fit = new TLegend(0.45, 0.7, 0.7, 0.5);
    leg_fit->SetFillColor(kWhite); leg_fit->SetBorderSize(0.);
    leg_fit->SetTextColor(1); leg_fit->SetTextSize(0.045);
    leg_fit->SetTextFont(62); 
    leg_fit->SetHeader(TitleHeader); //leg_fit->SetMargin(0.2); 
    leg_fit->AddEntry((TObject*)0, " ", "");
    leg_fit->AddEntry(DataTemplate, " Data", "e1px0");
    leg_fit->AddEntry(tTemplate, " b", "f");
    if (!FitOption.Contains("_TCF")) leg_fit->AddEntry(lgTemplate, " udsg + c", "f");
    else {
      leg_fit->AddEntry(hTemplate, " c", "f");
      leg_fit->AddEntry(tTemplate, " udsg", "f");
    }
    leg_fit->Draw();

    // Run2015B Style
    TLatex *tex = new TLatex(0.2,0.885,"CMS"); 
    tex->SetNDC(); 
    tex->SetTextAlign(13);
    tex->SetTextFont(61);
    tex->SetTextSize(0.07475);
    tex->SetLineWidth(2); 
    tex->Draw();           
                                                                            
    TLatex *tex2 = new TLatex(0.2,0.805,"Preliminary"); 
    tex2->SetNDC();
    tex2->SetTextAlign(13);
    tex2->SetTextFont(52);
    tex2->SetTextSize(0.05681);
    tex2->SetLineWidth(2);   
    tex2->Draw();      
    
    TLatex *text1 = new TLatex(0.98,0.95125, CampaignLuminosity); 
    text1->SetNDC();                                              
    text1->SetTextAlign(31);                          
    text1->SetTextFont(42);    
    text1->SetTextSize(0.04875);   
    text1->SetLineWidth(2);    
    text1->Draw();  
    // End Run2015B Style

    TString PlotName = "./Plots/" + TemplateVariable + "Fit_" + ThisTagger + "_" + EtaBin + "_" + PtBin + SystematicFlag + FitOption + PUWeighting + KinWeighting + Selection + Production;
    if (PrintPlot.Contains("png")) CanvasPtRelFit->Print(PlotName + ".png");
    if (PrintPlot.Contains("pdf")) CanvasPtRelFit->Print(PlotName + ".pdf");
    if (PrintPlot.Contains(  "C")) CanvasPtRelFit->Print(PlotName + ".C");
    
  }

  ScaleFactorResult PtRelSF;

  PtRelSF.Chi2NTag = FitContent[1].Chi2N;
  PtRelSF.FracTag = FitContent[1].bFraction;
  PtRelSF.FracTagError = FitContent[1].bFractionError;
  PtRelSF.Chi2NUntag = FitContent[0].Chi2N;
  PtRelSF.FracUntag = FitContent[0].bFraction;
  PtRelSF.FracUntagError = FitContent[0].bFractionError;

  if (!FitOption.Contains("_Pretag")) 
    PtRelSF.EffData = DataContent[1]*FitContent[1].bFraction/(DataContent[1]*FitContent[1].bFraction+DataContent[0]*FitContent[0].bFraction);
  else
    PtRelSF.EffData = DataContent[1]*FitContent[1].bFraction/(DataContent[0]*FitContent[0].bFraction);
  PtRelSF.EffDataError = PtRelEfficiencyError(DataContent[1], DataContent[0], 
					      FitContent[1].bFraction, FitContent[0].bFraction, 
					      FitContent[1].bFractionError, FitContent[0].bFractionError);

  if (!FitOption.Contains("_Pretag")) {
    PtRelSF.EffMC = MonteCarloBottomContent[1]/(MonteCarloBottomContent[1]+MonteCarloBottomContent[0]);
    PtRelSF.EffMCError = EfficiencyStatisticalError(MonteCarloBottomEntries[1], MonteCarloBottomEntries[0]);
  } else {
    PtRelSF.EffMC = MonteCarloBottomContent[1]/MonteCarloBottomContent[0];
    PtRelSF.EffMCError = BinomialError(MonteCarloBottomEntries[1], MonteCarloBottomEntries[0]);
  }

  PtRelSF.SF = PtRelSF.EffData/PtRelSF.EffMC;
  PtRelSF.SFError = RatioPropagationError(PtRelSF.EffData, PtRelSF.EffMC, PtRelSF.EffDataError, PtRelSF.EffMCError);

  PtRelSF.DataEventTag   = DataContent[1];
  PtRelSF.DataEventUntag = DataContent[0];

  return PtRelSF;

}

void PtRelAnalyzer::ComputePtRelScaleFactors(TString Taggers, TString SystematicFlag, TString FitOption, TString PrintPlot, TString PtBin, TString EtaBin, TString TemplateFlag) {

  TString TemplateFileName =  "./Templates/" + TemplateVariable + "_Templates" + TemplateFlag + PUWeighting + KinWeighting + Selection + Production + ".root";
  PtRelTemplateFile = TFile::Open(TemplateFileName); 

  for (int tg = 0; tg<nTaggers; tg++) 
    if (Taggers=="All" || TaggerName[tg].Contains(Taggers)) 
      for (int fpt = 0; fpt<nFitPtBins; fpt++) 
	if (PtBin=="All" || PtBin==FitPtBin[fpt])
	  for (int nb = 0; nb<nPtRelEtaBins; nb++) 
	    if (EtaBin==PtRelEtaBin[nb]) {
	      
	      ScaleFactorResult PtRelSF = PtRelScaleFactor(TaggerName[tg], EtaBin, FitPtBin[fpt], SystematicFlag, FitOption, PrintPlot);
	      
	      TString TableName = "./Tables/" + TemplateVariable + "Fit_" + TaggerName[tg] + "_" + EtaBin + "_" + FitPtBin[fpt] + SystematicFlag + FitOption  + PUWeighting + KinWeighting + Selection + Production + ".txt";
	      ofstream Table; Table.open(TableName);
	      
	      Table << EtaBin << " " << FitPtBin[fpt] << " " << std::endl
		    << "Efficiency MC = " << PtRelSF.EffMC << " +/- " <<  PtRelSF.EffMCError << std::endl
		    << "Chi2N         = " << PtRelSF.Chi2NUntag << " " <<  PtRelSF.Chi2NTag << std::endl
		    << "b fraction    = " << PtRelSF.FracUntag << " +/- " << PtRelSF.FracUntagError << std::endl
		    << "b fract. tag  = " << PtRelSF.FracTag << " +/- " << PtRelSF.FracTagError << std::endl
		    << "Eff. data     = " << PtRelSF.EffData << " +/- " <<  PtRelSF.EffDataError 
		    << " (Raw eff. = " <<  PtRelSF.DataEventTag 
		    << "/(" << PtRelSF.DataEventTag << " + " <<  PtRelSF.DataEventUntag  << ") = " 
		    << float(PtRelSF.DataEventTag)/(PtRelSF.DataEventTag+PtRelSF.DataEventUntag) << ") " << std::endl
		    << "Scale factor  = " << PtRelSF.SF << " +/- " << PtRelSF.SFError << std::endl << std::endl;
	      
	      Table.close(); 
	      
	    }
  
}

ScaleFactorResult PtRelAnalyzer::ScaleFactorResultFromTable(TString TableName) {

  struct ScaleFactorResult ThisScaleFactorResult;
  
  ifstream Table; Table.open(TableName);
  
  TString EtaBin, PtBin;
  Table >> EtaBin >> PtBin;
  
  TString Eff, MC, Eq, PM; 
  double ThisEff, ThisEffError;
  Table >> Eff >> MC >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisScaleFactorResult.EffMC = ThisEff; ThisScaleFactorResult.EffMCError = ThisEffError;

  TString Chi2NString; 
  TString /*double*/ Chi2N_Untag, Chi2N_Tag;
  Table >> Chi2NString >> Eq >> Chi2N_Untag >> Chi2N_Tag; 
    
  //ThisScaleFactorResult.Chi2NTag = Chi2N_Tag; ThisScaleFactorResult.Chi2NUntag = Chi2N_Untag;

  TString bottom, frac;
  Table >> bottom >> frac >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisScaleFactorResult.FracUntag = ThisEff; ThisScaleFactorResult.FracUntagError = ThisEffError;
  
  TString tag;
  Table >> bottom >> frac >> tag >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisScaleFactorResult.FracTag = ThisEff; ThisScaleFactorResult.FracTagError = ThisEffError;
  
  TString Data;
  Table >> Eff >> Data >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisScaleFactorResult.EffData = ThisEff; ThisScaleFactorResult.EffDataError = ThisEffError;
  
  TString OpenPar, Raw1, Raw2, Raw3, Eq2, EffString;
  Table >> OpenPar >> Eff >> Eq >> Raw1 >> Raw2 >> Raw3 >> Eq2 >> EffString;
  
  EffString.ReplaceAll(")", "");
  
  ThisScaleFactorResult.EffDataRaw = EffString.Atof();
  
  TString Scale, Factor;
  Table >> Scale >> Factor >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisScaleFactorResult.SF = ThisEff; ThisScaleFactorResult.SFError = ThisEffError;

  return ThisScaleFactorResult;
  
}

TString PtRelAnalyzer::ScaleFactorResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString KinWeightingFlag, TString SelectionFlag, TString ProductionFlag) {

  return ScaleFactorResultTableName(Tagger, EtaBin, PtBin, Systematic, PUWeightingFlag, KinWeightingFlag + SelectionFlag, ProductionFlag);

}

TString PtRelAnalyzer::ScaleFactorResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString Configuration, TString ProductionFlag) {
 
  if (Configuration.Contains("_PU") || Configuration.Contains("_PV") || Configuration.Contains("_PS"))
    PUWeightingFlag = "";
  
  TString TableName = "./Tables/" + TemplateVariable + "Fit_" + Tagger + "_" + EtaBin + "_" + PtBin + Systematic + PUWeightingFlag + Configuration + ProductionFlag + ".txt";
  
  return TableName;

}

TCanvas *PtRelAnalyzer::BTagPerformanceCanvas(TString Type, float CanvasHeight, TString LogAxes) {

  gStyle->SetTitleBorderSize( 0);
  gStyle->SetTitleFillColor (10);
  gStyle->SetStatFont       (42);
  gStyle->SetTitleFont      (42);

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat("");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

  //float CanvasHeight = 400.;
  if (Type=="Performance") CanvasHeight = 800.;

  TCanvas* ThisCanvas = new TCanvas("BTagPerformance", "BTagPerformance", 1200., CanvasHeight);
  
  if (Type=="Performance") {

    ThisCanvas->Divide(1,2,0,0);

    TPad *EfficiencyPad = (TPad*) ThisCanvas->GetPad(1);
    EfficiencyPad->SetPad(0.01, 0.51, 0.95, 0.95);

    TPad *ScaleFactorsPad = (TPad*) ThisCanvas->GetPad(2);
    ScaleFactorsPad->SetPad(0.01, 0.07, 0.95, 0.51);

    if (LogAxes.Contains("X")) { EfficiencyPad->SetLogx(); ScaleFactorsPad->SetLogx(); }
    if (LogAxes.Contains("Y")) { EfficiencyPad->SetLogy(); ScaleFactorsPad->SetLogy(); }

  } else if (Type=="ScaleFactors") {

    ThisCanvas->Divide(1,1,0,0);

    TPad *ScaleFactorsPad = (TPad*) ThisCanvas->GetPad(1);
    ScaleFactorsPad->SetPad(0.01, 0.07, 0.95, 0.95);

    if (LogAxes.Contains("X")) ScaleFactorsPad->SetLogx(); 
    if (LogAxes.Contains("Y")) ScaleFactorsPad->SetLogy();

  }
   
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.19);
  //gPad->SetGrid();
    
  gPad->SetTickx();
  gPad->SetTicky();

  return ThisCanvas;
  
}

void PtRelAnalyzer::FillBTagPerformanceHistograms(TString Tagger, TString EtaBin, TString Configuration, TString Systematic, int ColorIdx, TString DrawedSystematics) {

  MCEfficiency->Reset();
  DataEfficiency->Reset();
  DataMCSF->Reset();

  float MinEfficiency = 1.;
  float MaxEfficiency = 0.;

  TString ThisPUWeighting = PUWeighting;

  for (int ppt = 2; ppt<MCEfficiency->GetNbinsX(); ppt++) {

    TString TableName = ScaleFactorResultTableName(Tagger, EtaBin, FitPtBin[ppt-1], Systematic, ThisPUWeighting, Configuration, Production);

    ScaleFactorResult ThisScaleFactorResult;

    ThisScaleFactorResult.EffMC = 0.; ThisScaleFactorResult.EffData = 0.; ThisScaleFactorResult.SF = 0.;

    //if (DrawedSystematics!="NONE")

    if (VerboseSystematics)
      std::cout << TableName << " " << ThisScaleFactorResult.SF << std::endl;

    ThisScaleFactorResult = ScaleFactorResultFromTable(TableName);

    if (VerboseSystematics)
      std::cout << TableName << " " << ThisScaleFactorResult.SF << std::endl;

    MCEfficiency->SetBinContent(ppt+1, ThisScaleFactorResult.EffMC);
    MCEfficiency->SetBinError(ppt+1, ThisScaleFactorResult.EffMCError);

    DataEfficiency->SetBinContent(ppt+1, ThisScaleFactorResult.EffData);
    DataEfficiency->SetBinError(ppt+1, ThisScaleFactorResult.EffDataError);

    DataMCSF->SetBinContent(ppt+1, ThisScaleFactorResult.SF);
    DataMCSF->SetBinError(ppt+1, ThisScaleFactorResult.SFError);

    if (ThisScaleFactorResult.EffMC<MinEfficiency) MinEfficiency = ThisScaleFactorResult.EffMC;
    if (ThisScaleFactorResult.EffData<MinEfficiency) MinEfficiency = ThisScaleFactorResult.EffData;

    if (ThisScaleFactorResult.EffMC>MaxEfficiency) MaxEfficiency = ThisScaleFactorResult.EffMC;
    if (ThisScaleFactorResult.EffData>MaxEfficiency) MaxEfficiency = ThisScaleFactorResult.EffData;

  }

  MCEfficiency->GetYaxis()->SetRangeUser(0.05*(int(MinEfficiency/0.05)-2), 0.05*(int(MaxEfficiency/0.05)+3));

  MCEfficiency->SetMarkerSize(1.);
  MCEfficiency->SetMarkerStyle(25);
  MCEfficiency->SetMarkerColor(ColorIdx);
  
  MCEfficiency->GetXaxis()->SetLabelFont(42);
  MCEfficiency->GetYaxis()->SetLabelFont(42);
  MCEfficiency->GetXaxis()->SetTitleFont(42);
  MCEfficiency->GetYaxis()->SetTitleFont(42);
  MCEfficiency->SetTitle("");
  MCEfficiency->SetXTitle("Jet p_{T} [GeV]");
  MCEfficiency->SetYTitle("b-tag Efficiency #epsilon_{b}");
  MCEfficiency->GetYaxis()->SetTitleSize(0.09);
  MCEfficiency->GetYaxis()->SetTitleOffset(0.5);
  MCEfficiency->GetYaxis()->SetLabelSize(0.08);
  
  DataEfficiency->SetMarkerSize(1.);
  DataEfficiency->SetMarkerStyle(20);
  DataEfficiency->SetMarkerColor(ColorIdx);
  
  DataMCSF->SetMarkerSize(1.);
  DataMCSF->SetMarkerStyle(20);
  DataMCSF->SetMarkerColor(ColorIdx);
  
  DataMCSF->GetXaxis()->SetLabelFont(42);
  DataMCSF->GetYaxis()->SetLabelFont(42);
  DataMCSF->GetXaxis()->SetTitleFont(42);
  DataMCSF->GetYaxis()->SetTitleFont(42);
  DataMCSF->SetTitle("");
  DataMCSF->SetXTitle("Jet p_{T} [GeV]");
  DataMCSF->SetYTitle(Tagger + "Data/Sim. SF_{b}");
  //DataMCSF->GetXaxis()->SetTitleOffset(1.3);
  DataMCSF->GetXaxis()->SetTitleOffset(0.8);
  DataMCSF->GetXaxis()->SetTitleSize(0.09);
  DataMCSF->GetXaxis()->SetLabelSize(0.08);
  DataMCSF->GetYaxis()->SetTitle("Data/Sim. SF_{b}");
  DataMCSF->GetYaxis()->SetTitleSize(0.09);
  DataMCSF->GetYaxis()->SetTitleOffset(0.5);
  DataMCSF->GetYaxis()->CenterTitle(true);
  DataMCSF->GetYaxis()->SetLabelSize(0.09); 
  
  DataMCSF->SetMinimum(0.6);
  DataMCSF->SetMaximum(1.4);

  if (DrawedSystematics!="NONE") {

    DataMCSFSystematics->Add(DataMCSF);
    
    ComputeScaleFactorSystematics(Tagger, EtaBin, Configuration, DrawedSystematics);
 
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      if (EtaBin==PtRelEtaBin[nb]) 
	for (int ppt = 1; ppt<MCEfficiency->GetNbinsX(); ppt++)
	  DataMCSFSystematics->SetBinError(ppt+1, TotalScaleFactorSystematic[ppt-1][nb]);

    DataMCSFSystematics->SetMarkerStyle(20);
    DataMCSFSystematics->SetMarkerColor(2);
    DataMCSFSystematics->SetFillStyle(3005);
    DataMCSFSystematics->SetFillColor(ColorIdx+1);
    
  }
  
}

void PtRelAnalyzer::PlotBTagPerformance(TString PlotFlag, TString Tagger, TString EtaList[], TString ConfigurationList[], TString SystematicList[], float MaxJetPt, TString DrawedSystematics, TString Type) {

  TString LogAxes = "";
  if (MaxJetPt>300.) LogAxes = "X";
  TCanvas *ThisCanvas = BTagPerformanceCanvas(Type, 400., LogAxes);
  
  int nPlotPtBins = 0;
  float PlotPtEdge[1000];
  
  PlotPtEdge[0] = 0.;
  
  for (int fpt = 0; fpt<nFitPtBins+1; fpt++) 
    if (FitPtEdge[fpt]<=MaxJetPt) {
      
      nPlotPtBins++;
      PlotPtEdge[nPlotPtBins] = FitPtEdge[fpt];
      
    }
  
  MCEfficiency = new TH1D("MCEfficiency", "", nPlotPtBins, PlotPtEdge);
  DataEfficiency = new TH1D("DataEfficiency", "", nPlotPtBins, PlotPtEdge);
  DataMCSF = new TH1D("DataMCSF", "", nPlotPtBins, PlotPtEdge);
  DataMCSFSystematics = new TH1D("DataMCSFSystematics", "", nPlotPtBins, PlotPtEdge);
  
  TString Option;
  
  if (ConfigurationList[0]=="") ConfigurationList[0] = KinWeighting + Selection;

  int ColorIndex = 1;

  for (int nb = 0; nb<100 && EtaList[nb]!="-"; nb++)   
    for (int cfg = 0; cfg<100 && ConfigurationList[cfg]!="-"; cfg++)
      for (int is = 0; is<100 && SystematicList[is]!="-"; is++) {
	
	FillBTagPerformanceHistograms(Tagger, EtaList[nb], ConfigurationList[cfg], SystematicList[is], ColorIndex, DrawedSystematics);
	
	ThisCanvas->cd(1);

	gPad->SetTickx();
	gPad->SetTicky();
	
	if (Type!="ScaleFactors") {
	  
	  MCEfficiency->DrawCopy("pe1" + Option);
	  DataEfficiency->DrawCopy("pe1same");

	  if (PlotFlag=="Final") {

	    TLegend  *EffLegend = new TLegend(0.6, 0.95, 0.9, 0.75);
	    EffLegend->SetFillColor(kWhite); EffLegend->SetBorderSize(0.);
	    EffLegend->SetTextColor(1); EffLegend->SetTextSize(0.08);
	    EffLegend->SetTextFont(62); 
	    //EffLegend->SetHeader(TitleHeader); //EffLegend->SetMargin(0.2); 
	    //EffLegend->AddEntry((TObject*)0, " ", "");
	    EffLegend->AddEntry(DataEfficiency, " Data", "e1px0");
	    EffLegend->AddEntry(MCEfficiency, " Simulation", "e1px0");
	    EffLegend->Draw();
	  
	  } else {

	    double YLegend = 0.99 - cfg*0.08;
	    TLegend  *ThisLegend = new TLegend(0.2, YLegend, 0.4, YLegend-0.04);
	    ThisLegend->SetFillColor(kWhite); ThisLegend->SetBorderSize(0.);
	    ThisLegend->SetTextColor(cfg+1); ThisLegend->SetTextSize(0.06);
	    ThisLegend->SetTextFont(62);
	    //ThisLegend->AddEntry(DataEfficiency, " Data", "e1px0");
	    //ThisLegend->AddEntry(MCEfficiency, " Simulation", "e1px0");
	    TString LegendText = ConfigurationList[cfg];
	    LegendText.ReplaceAll("AwayTrgEmul", "AwayDiJet20TrgEmul");
	    LegendText.ReplaceAll("_KinPtBinsCentral_LowPtAway", "");
	    LegendText.ReplaceAll("TrgConf", "");
	    LegendText.ReplaceAll("TrgEmul", "");
	    ThisLegend->AddEntry(MCEfficiency, LegendText, "e1px0");
	    ThisLegend->Draw();

	  }

	}
	
	if (Type=="Performance") {

	  ThisCanvas->cd(2);
	
	  gPad->SetTickx();
	  gPad->SetTicky();

	}

	if (Type!="Efficiency") {

	  DataMCSF->DrawCopy("pe1" + Option);
	  
	  if (DrawedSystematics!="NONE" && ColorIndex==1)  {

	    DataMCSFSystematics->DrawCopy("pe2same");
	    DataMCSF->DrawCopy("pe1same");

	  }

	  if (PlotFlag=="Final") {

	    TLegend  *SFLegend = new TLegend(0.6, 0.95, 0.9, 0.75);
	    SFLegend->SetFillColor(kWhite); SFLegend->SetBorderSize(0.);
	    SFLegend->SetTextColor(1); SFLegend->SetTextSize(0.08);
	    SFLegend->SetTextFont(62); 
	    //SFLegend->SetHeader(TitleHeader); //SFLegend->SetMargin(0.2); 
	    //SFLegend->AddEntry((TObject*)0, " ", "");
	    SFLegend->AddEntry(DataMCSF, " Scale factors", "e1px0");
	    SFLegend->AddEntry(DataMCSFSystematics, " Scale factor uncertainty", "f");
	    SFLegend->Draw();
	  
	  }

	}

	Option = "same";
	
	ColorIndex++;

      }

  TString ThisPlotName = "BTag" + Type + "_DepOn" + PlotFlag + "_" + Tagger + "_" + EtaList[0] + PUWeighting + ConfigurationList[0] + Production;
  
  if (EtaList[1]!="-") ThisPlotName.ReplaceAll("_" + EtaList[0], "");
  if (ConfigurationList[1]!="-") ThisPlotName.ReplaceAll(ConfigurationList[0], "");

  ThisCanvas->Print("./Plots/" + ThisPlotName + ".png");
  ThisCanvas->Print("./Plots/" + ThisPlotName + ".pdf");
  
}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString Configuration, TString CentralFlag) {

  TString TableNameCentral = ScaleFactorResultTableName(Tagger, PtRelEtaBin[EtaBin], FitPtBin[PtBin], CentralFlag, PUWeighting, Configuration, Production);
  ScaleFactorResult CentralScaleFactorResult = ScaleFactorResultFromTable(TableNameCentral);
  
  float CentralScaleFactor = CentralScaleFactorResult.SF;
  
  ScaleFactorValue[PtBin][EtaBin] = CentralScaleFactorResult.SF;
  ScaleFactorSystematic[PtBin][EtaBin][0] = CentralScaleFactorResult.SFError;
  
  TotalScaleFactorSystematic[PtBin][EtaBin] = ScaleFactorSystematic[PtBin][EtaBin][0]*ScaleFactorSystematic[PtBin][EtaBin][0];
  
  for (int sfs = 1; sfs<nScaleFactorSystematics; sfs++) {

    ScaleFactorSystematic[PtBin][EtaBin][sfs] = 0.;
    
    TString VariationSign = "down";

    for (int is = 2*sfs-1; is<2*sfs+1; is++) {

      if (ScaleFactorSystematicName[sfs]=="_MuPt" || ScaleFactorSystematicName[sfs]=="_bCorr" || ScaleFactorSystematicName[sfs]=="_lCorr")
	VariationSign = "up";
      
      TString SystematicTableName = TableNameCentral; 
      SystematicTableName.ReplaceAll("_Central", FitSystematicName[is]);
      ScaleFactorResult SystematicScaleFactorResult = ScaleFactorResultFromTable(SystematicTableName);
      
      float SystematicScaleFactor = SystematicScaleFactorResult.SF - CentralScaleFactor;
      float SystematicSFError = SystematicScaleFactorResult.SFError;
      
      if (SystematicSFError<0.1 && fabs(SystematicScaleFactor)<0.1)
	if (fabs(SystematicScaleFactor)>fabs(ScaleFactorSystematic[PtBin][EtaBin][sfs])) {

	  if (VariationSign=="down") ScaleFactorSystematic[PtBin][EtaBin][sfs] = -1.*SystematicScaleFactor;
	  else ScaleFactorSystematic[PtBin][EtaBin][sfs] = SystematicScaleFactor;

	}

      VariationSign = "up";

    }
    
    if (VerboseSystematics)
      std::cout << Tagger << " " << FitPtBin[PtBin] << " " << ScaleFactorSystematicName[sfs] << " " << ScaleFactorSystematic[PtBin][EtaBin][sfs] <<std::endl;

    TotalScaleFactorSystematic[PtBin][EtaBin] += ScaleFactorSystematic[PtBin][EtaBin][sfs]*ScaleFactorSystematic[PtBin][EtaBin][sfs];
    
  }
  
  TotalScaleFactorSystematic[PtBin][EtaBin] = sqrt(TotalScaleFactorSystematic[PtBin][EtaBin]);
  
}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString CentralFlag) {

  ComputeScaleFactorSystematics(Tagger, PtBin, EtaBin, KinWeighting + Selection, CentralFlag);

}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Taggers, TString EtaBin, TString CentralFlag) {
  
  ComputeScaleFactorSystematics(Taggers, EtaBin, KinWeighting + Selection, CentralFlag);
  
}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Taggers, TString EtaBin, TString Configuration, TString CentralFlag) {
  
  for (int tg = 0; tg<nTaggers; tg++)
    if (Taggers=="All" || TaggerName[tg].Contains(Taggers) || Taggers.Contains(TaggerName[tg]))
      for (int fpt = 0; fpt<nFitPtBins; fpt++)
	for (int nb = 0; nb<nPtRelEtaBins; nb++)
	  if (EtaBin==PtRelEtaBin[nb])
	    ComputeScaleFactorSystematics(TaggerName[tg], fpt, nb, Configuration, CentralFlag);
  
}

void PtRelAnalyzer::EventCounter(TString DataType, int DataRange) {

  int nEventsFromTrees = 0;
  
  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
  GetDataRangeInformation(DataType, DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
  
  for (int tf = FirstTree; tf<=nTrees; tf++) {
    
    TString FileName = GetFileName(FileDirectoryName, tf);
    
    TFile *ThisTree = TFile::Open(FileName);
    
    if (!ThisTree) continue;
    
    TTree *tchain = GetChain(ThisTree, false);
    
    Int_t nentries = (Int_t)tchain->GetEntries();
    
    nEventsFromTrees += nentries;
    
  }

  float ExpectedEvents = 0;
  if (DataType=="QCDMu") ExpectedEvents = GeneratedEvents[DataRange];
  if (DataType=="QCD") ExpectedEvents = GeneratedEventsInclusive[DataRange];

  std::cout << "PtRelAnalyzer::EventCounter: Dataset " << DataType << " " << DataRangeName << " has " << nEventsFromTrees << " events (" << 100.*nEventsFromTrees/ExpectedEvents << "%)" << std::endl;

}

void PtRelAnalyzer::ComputePileUpWeights(TString PUMethod, TString DataType, int DataRange, TString PrimaryDataset) {
      
  const int nDataTypes = 4;
  TString DataTypeName[nDataTypes] = {"BTagMu", "QCDMu", "JetHT", "QCD"};
  
  TString PUFlagBTagMu = PUWeighting;
  PUFlagBTagMu.ReplaceAll("_PU", "_");
  PUFlagBTagMu.ReplaceAll("_PV", "_");
  PUFlagBTagMu.ReplaceAll("_PSV", "_PS");

  TString PUFlagJetHT = PUWeighting;
  PUFlagJetHT.ReplaceAll("_PU", "_");
  PUFlagJetHT.ReplaceAll("_PSV", "_PS");

  if (DataType.Contains("Loop")) {
    
    for (int dt = 0; dt<nDataTypes; dt++) 
      if (DataType.Contains(DataTypeName[dt])) {
	
	if (DataTypeName[dt]=="QCD" && !DataType.Contains("QCDX")) continue;

	std::cout << "PtRelAnalyzer::ComputePileUpWeights: Dataset " << DataTypeName[dt] << std::endl;

	TString DataTrees = DataTypeName[dt];

	GetTriggerPrescales(DataTypeName[dt]);
	GetPileUpWeights(DataTypeName[dt]);

	TH1F *PVMultDefault[nTriggers], *PVMultWeighted[nTriggers];
	for (int tr = 0; tr<nTriggers; tr++) {

	  PVMultDefault[tr] = new TH1F("PVMultDefault_" + DataTypeName[dt] + TriggerName[tr], "", 60, 0., 60.);
	  PVMultWeighted[tr] = new TH1F("PVMultWeighted_" + DataTypeName[dt] + TriggerName[tr], "", 60, 0., 60.);
	  
	}
	
	int nDataRanges = GetNumberOfDataRanges(DataTrees);

	for (int ph = 0; ph<nDataRanges; ph++) {
	 
	  if (DataRange!=-1 && DataRange!=ph) continue; 

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    PVMultDefault[tr]->Reset();
	    PVMultWeighted[tr]->Reset();
	    
	  }
	  
	  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
	  GetDataRangeInformation(DataTrees, ph, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
	
	  std::cout << "PtRelAnalyzer::ComputePileUpWeights:   Data Range Name " << DataRangeName << std::endl;
	  
	  for (int tf = FirstTree; tf<=nTrees; tf++) {
	    
	    TString FileName = GetFileName(FileDirectoryName, tf);
	    
	    TFile *ThisTree = TFile::Open(FileName);
	  
	    std::cout << "PtRelAnalyzer::ComputePileUpWeights:     Tree " << tf << "/" << nTrees << std::endl;
  
	    if (!ThisTree) continue;

	    TTree *tchain = GetChain(ThisTree, false);

	    Int_t nentries = (Int_t)tchain->GetEntries();
	        
	    for (Int_t i = 0; i<nentries; i++) {
      
	      tchain->GetEntry(i);
      
	      if (!DataTypeName[dt].Contains("QCD")) {
		if (PUWeighting.Contains("ICHEP2016") && Run>276811) continue;
	      }
	    
	      if (nPV>0) {
		
		int EventPileUp = (DataTypeName[dt].Contains("QCD")) ? nPUtrue : 1;
		if (PUWeighting.Contains("PV") || PUWeighting.Contains("PSV")) EventPileUp = nPV;
		if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
		
		for (int tr = 0; tr<nTriggers; tr++) {

		  float ThisPrescaleWeight = TriggerPrescaleWeight(DataTypeName[dt], tr);
		  
		  TString ThisTriggerName = TriggerName[tr];
		  if (DataTypeName[dt]=="JetHT" || DataTypeName[dt]=="QCD")
		    ThisTriggerName = JetTriggerName[tr];

		  bool PassTrigger = true;
		  if (DataTypeName[dt]!="QCD") 
		    PassTrigger = PassTriggerBit(ThisTriggerName);
		   
		  if (DataTypeName[dt]=="QCD" || DataTypeName[dt]=="JetHT")
		    if (!PassEventTriggerEmulation(tr, DataTypeName[dt])) PassTrigger = false;

		  if (PassTrigger) {
		    
		    PVMultDefault[tr]->Fill(nPV, PtHatWeight);
		    PVMultWeighted[tr]->Fill(nPV, PtHatWeight*PileUpWeight[tr][EventPileUp][0]*ThisPrescaleWeight); 
		    
		  }

		}
		
	      }

	    }
	   
	    ThisTree->Close();
 
	  }
	  
	  std::cout << "PtRelAnalyzer::ComputePileUpWeights:   Saving PV Multiplicity Histograms" << std::endl;
	  
	  TString OutFileName = "./Weights/PileUp/PVMultTriggered_" + DataTypeName[dt] + "_" + DataRangeName + PUWeighting + ".root";
	  if (DataTypeName[dt]=="BTagMu") 
	    OutFileName.ReplaceAll(PUWeighting, PUFlagBTagMu);
	  if (DataTypeName[dt]=="JetHT") 
	    OutFileName.ReplaceAll(PUWeighting, PUFlagJetHT);
	  TFile *OutFile = new TFile(OutFileName, "recreate");
	  
	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    PVMultDefault[tr]->Write();
	    PVMultWeighted[tr]->Write();
	    
	  }
	
	  OutFile->Close();
	  
	}
	
      }
    
  }

  if (DataType.Contains("Merge")) {
    
    for (int dt = 0; dt<nDataTypes; dt++) 
      if (DataType.Contains(DataTypeName[dt])) {
	
	if (DataTypeName[dt]=="QCD" && !DataType.Contains("QCDX")) continue;
	 
	std::cout << "PtRelAnalyzer::ComputePileUpWeights: Merging dataset " << DataTypeName[dt] << std::endl;

	TString DataTrees = DataTypeName[dt];

	TH1F *PVMultDefault[nTriggers], *PVMultWeighted[nTriggers];
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  PVMultDefault[tr] = new TH1F("PVMultDefault_" + DataTypeName[dt] + TriggerName[tr], "", 60, 0., 60.);
	  PVMultWeighted[tr] = new TH1F("PVMultWeighted_" + DataTypeName[dt] + TriggerName[tr], "", 60, 0., 60.);
	  
	}
	
	int nDataRanges = GetNumberOfDataRanges(DataTrees);
	
	int MaxDataRange = nDataRanges;
	if (DataRange>=0) MaxDataRange = DataRange + 1;
	for (int ph = 0; ph<MaxDataRange; ph++) {
	  
	  TString DataRangeName = GetDataRangeName(DataTrees, ph);
	  
	  TString InFileName = "./Weights/PileUp/PVMultTriggered_" + DataTypeName[dt] + "_" + DataRangeName + PUWeighting + ".root";
	  if (DataTypeName[dt]=="BTagMu") 
	    InFileName.ReplaceAll(PUWeighting, PUFlagBTagMu);
	  if (DataTypeName[dt]=="JetHT") 
	    InFileName.ReplaceAll(PUWeighting, PUFlagJetHT);
	  TFile *InFile = TFile::Open(InFileName);
	  
	  for (int tr = 0; tr<nTriggers; tr++) {

	    TH1F *ThisPVMultDefault = (TH1F*) InFile->Get("PVMultDefault_" + DataTypeName[dt] + TriggerName[tr]);	    
	    PVMultDefault[tr]->Add(ThisPVMultDefault);
	    
	    TH1F *ThisPVMultWeighted = (TH1F*) InFile->Get("PVMultWeighted_" + DataTypeName[dt] + TriggerName[tr]);	
	    PVMultWeighted[tr]->Add(ThisPVMultWeighted);
	    
	  }
	  
	  InFile->Close();
	  
	}
	
	TString OutFileName = "./Weights/PileUp/PVMultTriggered_" + DataTypeName[dt] + PUWeighting + ".root";
	if (DataTypeName[dt]=="BTagMu") OutFileName.ReplaceAll(PUWeighting, PUFlagBTagMu);
	if (DataTypeName[dt]=="JetHT") OutFileName.ReplaceAll(PUWeighting, PUFlagJetHT);
	TFile *OutFile = new TFile(OutFileName, "recreate");
	  
	for (int tr = 0; tr<nTriggers; tr++) {
	 
	  if (PUWeighting.Contains("_PS") && !DataTypeName[dt].Contains("QCD")) {
	  
	    float TotalIntegral = PVMultWeighted[tr]->Integral(0, PVMultWeighted[tr]->GetNbinsX()+1);
	    PVMultWeighted[tr]->Scale(PVMultWeighted[tr]->GetEntries()/TotalIntegral);
	    
	  }
 
	  PVMultDefault[tr]->Write();
	  PVMultWeighted[tr]->Write();
	  
	}
	
	OutFile->Close();
	
      }

  }
  
  if (PUMethod!="") {

    if (PUMethod.Contains("TestPV")) {
      
      TFile *DataFile = TFile::Open("./Weights/PileUp/PVMultTriggered_BTagMu" + PUFlagBTagMu + ".root");
      TString QCDFileName = "./Weights/PileUp/PVMultTriggered_" + DataType + PUWeighting + ".root";
      if (DataType=="JetHT") QCDFileName.ReplaceAll(PUWeighting, PUFlagJetHT);
      TFile *QCD_File = TFile::Open(QCDFileName);
      
      TCanvas *ThisCanvas = BTagPerformanceCanvas("PSTest", 800.);
      ThisCanvas->Divide(1, 1);	
      
      for (int tr = 0; tr<nTriggers; tr++) {
	
	TH1F *NumberOfPVs[2][2];

	NumberOfPVs[0][0] = (TH1F*) DataFile->Get("PVMultDefault_BTagMu" + TriggerName[tr]);
	NumberOfPVs[1][0] = (TH1F*) QCD_File->Get("PVMultDefault_" + DataType + TriggerName[tr]);
	NumberOfPVs[0][0]->SetXTitle("Number of PVs");
	NumberOfPVs[1][0]->SetXTitle("Number of PVs");
	
	NumberOfPVs[0][1] = (TH1F*) DataFile->Get("PVMultWeighted_BTagMu" + TriggerName[tr]);
	NumberOfPVs[1][1] = (TH1F*) QCD_File->Get("PVMultWeighted_" + DataType + TriggerName[tr]);
	NumberOfPVs[0][1]->SetXTitle("Number of PVs");
	NumberOfPVs[1][1]->SetXTitle("Number of PVs");
	
	ThisCanvas->cd(1);
	
	float DataEvents = NumberOfPVs[0][0]->Integral();
	float QCD_Events = NumberOfPVs[1][0]->Integral();
	NumberOfPVs[1][0]->Scale(DataEvents/QCD_Events);
	NumberOfPVs[1][1]->Scale(DataEvents/QCD_Events);
	
	NumberOfPVs[0][0]->SetLineColor(2);
	NumberOfPVs[0][0]->SetLineWidth(2);
	NumberOfPVs[0][1]->SetMarkerStyle(20);
	NumberOfPVs[0][1]->SetMarkerColor(2);
	NumberOfPVs[1][0]->SetLineColor(4);
	NumberOfPVs[1][0]->SetLineWidth(2);
	NumberOfPVs[1][1]->SetMarkerStyle(24);
	NumberOfPVs[1][1]->SetMarkerColor(4);
	
	NumberOfPVs[0][0]->Draw("hist");
	NumberOfPVs[1][0]->Draw("histsame");
	if (PUFlagBTagMu.Contains("_PS")) NumberOfPVs[0][1]->Draw("pesame");
	NumberOfPVs[1][1]->Draw("pesame");
	
	TLegend  *leg = new TLegend(0.55, 0.85, 0.75, 0.65);
	leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
	leg->SetTextColor(1); leg->SetTextSize(0.04);
	leg->AddEntry(NumberOfPVs[0][0], " BTagMu", "f");
	leg->AddEntry(NumberOfPVs[1][0], " QCDMu", "f");
	if (PUFlagBTagMu.Contains("_PS")) leg->AddEntry(NumberOfPVs[0][1], " BTagMu weighted", "lep");
	leg->AddEntry(NumberOfPVs[1][1], " QCDMu weighted", "lep");
	leg->SetMargin(0.2); leg->Draw();
	
	ThisCanvas->Print("./Plots/" + PUMethod + "_" + DataType + TriggerName[tr] + PUWeighting + ".png");
	ThisCanvas->Print("./Plots/" + PUMethod + "_" + DataType + TriggerName[tr] + PUWeighting + ".pdf");

      }
      
    } else if (PUMethod.Contains("PV") || PUMethod.Contains("PSV")) {
      
      for (int pd = 0; pd<nDataTypes; pd++)
	if (PrimaryDataset.Contains(DataTypeName[pd]))
	  for (int dt = 0; dt<nDataTypes; dt++)
	    if (pd!=dt && DataType.Contains(DataTypeName[dt])) {
	
	      if (DataTypeName[dt]=="QCD" && !DataType.Contains("QCDX")) continue;

	      if (DataTypeName[dt]!="JetHT" || !PUMethod.Contains("PSV")) {

		TString RefDataFile = PUFlagBTagMu;
		if (PUMethod.Contains("_PSV")) {
		  RefDataFile = PUMethod;
		  RefDataFile.ReplaceAll("_PSV", "_PS");
		}
		TFile *DataFile = TFile::Open("./Weights/PileUp/PVMultTriggered_" + DataTypeName[pd] + RefDataFile + ".root");
		TFile *Den_File = TFile::Open("./Weights/PileUp/PVMultTriggered_" + DataTypeName[dt] + PUWeighting + ".root");
		
		for (int tr = 0; tr<nTriggers; tr++) {
		  
		  TH1F *NumberOfPVs[2][2];
		  
		  TString RefHistogram = "PVMultDefault_";
		  if (PUMethod.Contains("_PSV")) RefHistogram = "PVMultWeighted_";
		  NumberOfPVs[0][0] = (TH1F*) DataFile->Get(RefHistogram + DataTypeName[pd] + TriggerName[tr]);
		  NumberOfPVs[1][0] = (TH1F*) Den_File->Get("PVMultDefault_" + DataTypeName[dt] + TriggerName[tr]);
		  
		  float DataEvents = NumberOfPVs[0][0]->Integral(0, 61);
		  float Den_Events = NumberOfPVs[1][0]->Integral(0, 61);
		  NumberOfPVs[1][0]->Scale(DataEvents/Den_Events);
		  
		  NumberOfPVs[0][0]->Divide(NumberOfPVs[1][0]);
		  
		  ofstream OutFile; OutFile.open("./Weights/PileUp/PVMultTriggered" + TriggerName[tr] + "_" + DataTypeName[dt] + PUMethod + ".txt");
		  
		  for (int mpv = 1; mpv<60; mpv++)
		    OutFile << mpv << " " << NumberOfPVs[0][0]->GetBinContent(mpv+1) << std::endl;  
		  
		}
	      
	      } else {

		for (int tr = 0; tr<nTriggers; tr++) {

		  ofstream OutFile; OutFile.open("./Weights/PileUp/PVMultTriggered" + TriggerName[tr] + "_" + DataTypeName[dt] + PUMethod + ".txt");
		  
		  for (int mpv = 1; mpv<60; mpv++)
		    OutFile << mpv << " " << 1. << std::endl; 
		  
		}

	      }

	    }
      
    } else if (PUMethod.Contains("PS")) {

      TString PUSyst = "";
      if (PUMethod.Contains("m05")) PUSyst = "_m05";
      if (PUMethod.Contains("p05")) PUSyst = "_p05";

      TString PileUpRootFileName = "./Weights/PileUp/PileUp" + PUWeighting + PUSyst + ".root";
      PileUpRootFileName.ReplaceAll("_PS", "_");
      TFile *PileUpDataFile = TFile::Open(PileUpRootFileName);
      TH1F *PileUpData = (TH1F*) PileUpDataFile->Get("pileup");
      
      int PileUpDataBins = PileUpData->GetNbinsX();
      float PileUpDataIntegral = PileUpData->Integral(0, PileUpDataBins+1);
      PileUpData->Scale(1./PileUpDataIntegral);

      for (int dt = 0; dt<nDataTypes; dt++)
	if (DataType.Contains(DataTypeName[dt])) {
	  
	  if (DataTypeName[dt]=="QCD" && !DataType.Contains("QCDX")) continue;

	  for (int tr = 0; tr<nTriggers; tr++) {
		 
	    ofstream OutFile; OutFile.open("./Weights/PileUp/PileUpWeights" + TriggerName[tr] + "_" + DataTypeName[dt] + PUWeighting + PUSyst + ".txt");
	
	    if (DataTypeName[dt]=="JetHT") {
	      for (int npu = 0; npu<nMaxPU; npu++) 
		OutFile << npu << " " << 1. << std::endl;
	    } else {
	      for (int npu = 0; npu<nMaxPU; npu++) {
		if (PileUpScenario[npu]>0.) 
		  OutFile << npu << " " << PileUpData->GetBinContent(npu+1)/PileUpScenario[npu] << std::endl;  
		else 
		  OutFile << npu << " " << 0. << std::endl;
	      }
	    }
	    
	  }
	  
	}

    }
	    
  }
  
}

void PtRelAnalyzer::CompareDataToMC(TString Variable, TString EtaBin, TString PtBin, int Rebin, TString LightTemplates, TString DataType) {

  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + KinWeighting + Selection + Production + ".root";
  TFile *TemplateFile = TFile::Open(TemplateFileName);

  TCanvas *ThisCanvas = BTagPerformanceCanvas("DataVsMC", 800.);
  ThisCanvas->Divide(1, 1); ThisCanvas->cd(1);
 
  const int nVariables = 5;
  TString VariableName[nVariables] = {"jetPt", "jetEta", "PV", "muonPt", "muonDR"};
  TString XTitle[nVariables] = {"muon-jet p_{T} (GeV/c)", "muon-jet #eta", "Number of PVs", "muon p_{T} (GeV/c)", "#DeltaR(#mu, jet)"};

  TString Systematic = "_Central";

  for (int vr = 0; vr<nVariables; vr++)
    if (VariableName[vr]==Variable) {

      if (DataType=="JetHT" && VariableName[vr].Contains("muon")) continue;

      for (int bpt = 0; bpt<nFitPtBins; bpt++)
	if (PtBin=="All" || PtBin==FitPtBin[bpt]) 
	  for (int beta = 0; beta<nPtRelEtaBins; beta++) 
	    if(EtaBin=="All" || EtaBin==PtRelEtaBin[beta]) {

	      TString HistoName = "Variable_" + DataType + "_PtBin_EtaBin" + Systematic;
	      HistoName.ReplaceAll("Variable", VariableName[vr]);
	      HistoName.ReplaceAll("PtBin", FitPtBin[bpt]); 
	      HistoName.ReplaceAll("EtaBin", PtRelEtaBin[beta]);
	      TH1D *DataHisto = (TH1D*)TemplateFile->Get(HistoName);
  
	      HistoName.ReplaceAll("BTagMu", "QCDMu");
	      HistoName.ReplaceAll("JetHT", "QCD");
	      TH1D *QCD_Histo = (TH1D*)TemplateFile->Get(HistoName);
 
	      DataHisto->Rebin(Rebin);
	      QCD_Histo->Rebin(Rebin);

	      float DataIntegral = DataHisto->Integral();
	      float QCD_Integral = QCD_Histo->Integral();

	      QCD_Histo->Scale(DataIntegral/QCD_Integral);

	      int first_bin = 1; bool Stop = false;
	      for (int bn = 1; bn<=DataHisto->GetNbinsX() && !Stop; bn++)
		if (DataHisto->GetBinContent(bn)==0.) first_bin = bn;
		else Stop = true;

	      int last_bin = DataHisto->GetNbinsX(); Stop = false;
	      for (int bn = DataHisto->GetNbinsX(); bn>=1 && !Stop; bn--)
		if (DataHisto->GetBinContent(bn)==0.) last_bin = bn;
		else Stop = true;

	      DataHisto->GetXaxis()->SetRange(first_bin, last_bin);
	      QCD_Histo->GetXaxis()->SetRange(first_bin, last_bin);
  
	      TPad *PadCC1 = (TPad*)ThisCanvas->GetPad(1); 

	      PadCC1->cd();
    
	      QCD_Histo->SetTitle("");
	      QCD_Histo->SetXTitle(XTitle[vr]);
    
	      DataHisto->SetTitle("");
	      DataHisto->SetXTitle(XTitle[vr]);

	      QCD_Histo->SetFillColor(38);
	      
	      DataHisto->SetMarkerStyle(20);
	      DataHisto->SetMarkerColor(46);

	      if (DataHisto->GetMaximum()>QCD_Histo->GetMaximum())
		DataHisto->Draw("pe");
	      else
		QCD_Histo->Draw("hist");
	      
	      QCD_Histo->Draw("histsame");
	      DataHisto->Draw("pesame");
    
	      TLegend  *leg = new TLegend(0.65, 0.35, 0.9, 0.25);
	      leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
	      leg->SetTextColor(1); leg->SetTextSize(0.04);
	      leg->AddEntry(DataHisto, " Data", "lep");
	      leg->AddEntry(QCD_Histo, " MC", "f");
	      leg->SetMargin(0.2); leg->Draw();
 
	      ThisCanvas->Print("./Plots/DataVsMC_" + VariableName[vr] + "_" + FitPtBin[bpt] + "_" + PtRelEtaBin[beta] + Systematic + PUWeighting + KinWeighting + Selection + Production + ".png");
	      
	    }
      
    }
  
}

void PtRelAnalyzer::StudyAwayJetSelection(TString Step, int DataRange) {

  TString AwayJetCut[2] = {"_Medium", "_Tight"};
  TString MuonJetFlavour[3] = {"_b", "_c", "_lg"};
  TString AwayJetFlavour[5] = {"_b", "_c", "_lg", "_GS", "_All"};
  TString ThirdJetDR[2] = {"_0", "_1p5"};

  if (Step=="FillHistograms") {

    TH1F *AwayJetDiscriminator[nPtRelPtBins][2][3][5];
    TH1F *ThirdJetDiscriminator[nPtRelPtBins][2][3][5][2];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
      for (int ajc = 0; ajc<2; ajc++)
	for (int mjf = 0; mjf<3; mjf++)
	  for (int ajf = 0; ajf<5; ajf++) {
	    
	    TString JetName = PtRelPtBin[bpt] + AwayJetCut[ajc] + MuonJetFlavour[mjf] + AwayJetFlavour[ajf];
	    AwayJetDiscriminator[bpt][ajc][mjf][ajf] = new TH1F("AwayJet" + JetName, "", 400, 0., 4.);
				  
	    for (int tjdr= 0; tjdr<2; tjdr++)
	      ThirdJetDiscriminator[bpt][ajc][mjf][ajf][tjdr] = new TH1F("ThirdJet" + JetName + ThirdJetDR[tjdr], "", 400, 0., 4.);
				    
	  }
    
    TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
    GetDataRangeInformation("QCDMu", DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
  
    for (int tf = FirstTree; tf<=nTrees; tf++) {
    
      TString FileName = GetFileName(FileDirectoryName, tf);
      
      TFile *ThisTree = TFile::Open(FileName);
      
      if (!ThisTree) continue;
      
      TTree *tchain = GetChain(ThisTree, false);
      
      std::cout << "  Getting sample: MuPt5 " << DataRangeName << " File " << tf << "/" << nTrees << std::endl;
      
      Int_t nentries = (Int_t)tchain->GetEntries();
      
      for (Int_t i = 0; i<nentries; i++) {
	
	tchain->GetEntry(i);
	
	int iMu; 
	int jMu = GetPFMuonJet(&iMu);
	
	if (jMu>=0 && nPV>0 && Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1]) {
	  
	  int ptBin = -1;
	  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	    if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	  
	  bool PassTrigger = false;
	  
	  for (int trg = 0; trg<nTriggers; trg++) 
	    if (Jet_pt[jMu]>MinPtJetTrigger[trg] && Jet_pt[jMu]<MaxPtJetTrigger[trg])     
	      if (PassTriggerBit(TriggerName[trg]))
		PassTrigger = true;
	  
	  if (PassTrigger) {
	    
	    for (int ajc = 0; ajc<2; ajc++) {
	      
	      TString ThisAwayTaggerName = "JBPM";
	      if (ajc==1) ThisAwayTaggerName = "JBPT";
	      
	      int aJet = GetAwayJet(ThisAwayTaggerName, jMu, 1.5, true);
	      
	      if (aJet>=0) {
		
		int mjf = 2;
		if (fabs(Jet_flavour[jMu])==5) mjf = 0;
		if (fabs(Jet_flavour[jMu])==4) mjf = 1;
		
		int ajf = 2;
		if (fabs(Jet_flavour[aJet])==5) ajf = 0;
		if (fabs(Jet_flavour[aJet])==4) ajf = 1;
		if (IsFromGluonSplittingFromHadron(aJet, 5)) ajf = 3;
		
		AwayJetDiscriminator[ptBin][ajc][mjf][4]->Fill(Jet_Bprob[aJet], PtHatWeight);
		AwayJetDiscriminator[ptBin][ajc][mjf][ajf]->Fill(Jet_Bprob[aJet], PtHatWeight);
		
		int tJet0 = -1, tJet1 = -1;
		float DiscJet0 = -1., DiscJet1 = -1.;
		
		for (int ijet = 0; ijet<nJet; ijet++)
		  if (ijet!=aJet && ijet!=jMu) 
		    if (Jet_pt[ijet]>30. && fabs(Jet_eta[ijet])<2.4) {
		      
		      if (Jet_Bprob[ijet]>DiscJet0) {
			
			tJet0 = ijet;
			DiscJet0 = Jet_Bprob[ijet];
			
		      }
		      
		      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Jet_eta[jMu], Jet_phi[jMu])>1.5)
			if (Jet_Bprob[ijet]>DiscJet1) {
			  
			  tJet1 = ijet;
			  DiscJet1 = Jet_Bprob[ijet];
			  
			}
		      
		    }
		
		ThirdJetDiscriminator[ptBin][ajc][mjf][4][0]->Fill(Jet_Bprob[tJet0], PtHatWeight);
		ThirdJetDiscriminator[ptBin][ajc][mjf][ajf][0]->Fill(Jet_Bprob[tJet0], PtHatWeight);
		
		ThirdJetDiscriminator[ptBin][ajc][mjf][4][1]->Fill(Jet_Bprob[tJet1], PtHatWeight);
		ThirdJetDiscriminator[ptBin][ajc][mjf][ajf][1]->Fill(Jet_Bprob[tJet1], PtHatWeight);
		
	      }
	      
	    } 
	    
	  }
	  
	}
	
      }
      
    }

    std::cout << "Saving histograms" << std::endl;
  
    TString OutputFileName = "./Templates/Histograms/AwayJetSelection_" + DataRangeName + ".root";

    TFile *OutFile = new TFile(OutputFileName, "recreate");

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
      for (int ajc = 0; ajc<2; ajc++)
	for (int mjf = 0; mjf<3; mjf++)
	  for (int ajf = 0; ajf<5; ajf++) {
	    
	    AwayJetDiscriminator[bpt][ajc][mjf][ajf]->Write();
				  
	    for (int tjdr= 0; tjdr<2; tjdr++)
	      ThirdJetDiscriminator[bpt][ajc][mjf][ajf][tjdr]->Write();
				    
	  }
    
    OutFile->Close();

  } else if (Step=="MergeHistograms") {

    std::cout << "Merging histograms: " << std::endl;

    TH1F *AwayJetDiscriminator[nPtRelPtBins][2][3][5];
    TH1F *ThirdJetDiscriminator[nPtRelPtBins][2][3][5][2];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
      for (int ajc = 0; ajc<2; ajc++)
	for (int mjf = 0; mjf<3; mjf++)
	  for (int ajf = 0; ajf<5; ajf++) {
	    
	    TString JetName = PtRelPtBin[bpt] + AwayJetCut[ajc] + MuonJetFlavour[mjf] + AwayJetFlavour[ajf];
	    AwayJetDiscriminator[bpt][ajc][mjf][ajf] = new TH1F("AwayJet" + JetName, "", 400, 0., 4.);
				  
	    for (int tjdr= 0; tjdr<2; tjdr++)
	      ThirdJetDiscriminator[bpt][ajc][mjf][ajf][tjdr] = new TH1F("ThirdJet" + JetName + ThirdJetDR[tjdr], "", 400, 0., 4.);
				    
	  }
  
    for (int ph = 0; ph<DataRange; ph++) {
      
      TString InputFileName = "./Templates/Histograms/AwayJetSelection_" + MonteCarloPtHatRange[ph] + ".root";
      TFile *HistogramFile = TFile::Open(InputFileName);

      for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
	for (int ajc = 0; ajc<2; ajc++)
	  for (int mjf = 0; mjf<3; mjf++)
	    for (int ajf = 0; ajf<5; ajf++) {

	      TString JetName = PtRelPtBin[bpt] + AwayJetCut[ajc] + MuonJetFlavour[mjf] + AwayJetFlavour[ajf];
	      TH1F *ThisAwayJet = (TH1F*) HistogramFile->Get("AwayJet" + JetName);
	      AwayJetDiscriminator[bpt][ajc][mjf][ajf]->Add(ThisAwayJet);
				  
	      for (int tjdr= 0; tjdr<2; tjdr++) {

		TH1F *ThisThirdJet = (TH1F*) HistogramFile->Get("ThirdJet" + JetName + ThirdJetDR[tjdr]);
		ThirdJetDiscriminator[bpt][ajc][mjf][ajf][tjdr]->Add(ThisThirdJet);
	  
	      }

	    }

    }

    std::cout << "Saving histograms" << std::endl;
  
    TString OutputFileName = "./Templates/Histograms/AwayJetSelection.root";

    TFile *OutFile = new TFile(OutputFileName, "recreate");

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
      for (int ajc = 0; ajc<2; ajc++)
	for (int mjf = 0; mjf<3; mjf++)
	  for (int ajf = 0; ajf<5; ajf++) {
	    
	    AwayJetDiscriminator[bpt][ajc][mjf][ajf]->Write();
				  
	    for (int tjdr= 0; tjdr<2; tjdr++)
	      ThirdJetDiscriminator[bpt][ajc][mjf][ajf][tjdr]->Write();
				    
	  }
    
    OutFile->Close();

  }

}

void PtRelAnalyzer::CompareHistograms(TCanvas *PlotCanvas, TString Variable, TH1D *FirstHisto, TH1D *SecondHisto, TString PlotName, int Rebin) {

  const int nVariables = 6;
  TString VariableName[nVariables] = {"jetPt", "jetEta", "PV", "muonPt", "muonDR", "PtRel"};
  TString XTitle[nVariables] = {"muon-jet p_{T} [GeV]", "muon-jet #eta", "Number of PVs", "muon p_{T} [GeV/c]", "#DeltaR(#mu, jet)", "p_{T}^{rel} [GeV]"};
  
  bool GoodVariable = false;

  for (int vr = 0; vr<nVariables; vr++)
    if (VariableName[vr]==Variable) {
      
      GoodVariable = true;
      
      FirstHisto->Rebin(Rebin);
      SecondHisto->Rebin(Rebin);
      
      float FirstIntegral = FirstHisto->Integral();
      float SecondIntegral = SecondHisto->Integral();

      SecondHisto->Scale(FirstIntegral/SecondIntegral);

      int first_bin = 1; bool Stop = false;
      for (int bn = 1; bn<=FirstHisto->GetNbinsX() && !Stop; bn++)
	if (FirstHisto->GetBinContent(bn)==0.) first_bin = bn;
	else Stop = true;
      
      int last_bin = FirstHisto->GetNbinsX(); Stop = false;
      for (int bn = SecondHisto->GetNbinsX(); bn>=1 && !Stop; bn--)
	if (FirstHisto->GetBinContent(bn)==0.) last_bin = bn;
	else Stop = true;
      
      FirstHisto->GetXaxis()->SetRange(first_bin, last_bin);
      SecondHisto->GetXaxis()->SetRange(first_bin, last_bin);
  
      TPad *PadCC1 = (TPad*)PlotCanvas->GetPad(1); 
      
      PadCC1->cd();
      
      SecondHisto->SetTitle("");
      SecondHisto->SetXTitle(XTitle[vr]);
      
      FirstHisto->SetTitle("");
      FirstHisto->SetXTitle(XTitle[vr]);
      
      SecondHisto->SetFillColor(38);
      
      FirstHisto->SetMarkerStyle(20);
      FirstHisto->SetMarkerColor(46);
      
      if (FirstHisto->GetMaximum()>SecondHisto->GetMaximum())
	FirstHisto->Draw("pe");
      else
	SecondHisto->Draw("hist");
      
      SecondHisto->Draw("histsame");
      FirstHisto->Draw("pesame");

      TLegend  *leg = new TLegend(0.65, 0.35, 0.9, 0.25);
      leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
      leg->SetTextColor(1); leg->SetTextSize(0.04);
      leg->AddEntry(FirstHisto, FirstHisto->GetName(), "lep");
      leg->AddEntry(SecondHisto, SecondHisto->GetName(), "f");
      leg->SetMargin(0.2); leg->Draw(); 

      PlotCanvas->Print(PlotName); 
      
    }

  if (!GoodVariable) 
    std::cout << "PtRelAnalyzer::CompareHistograms: variable " << Variable << " no supported" << std::endl; 

}

void PtRelAnalyzer::CompareTemplates(TString Variable, TString FirstFileName, TString SecondFileName, TString FirstName, TString SecondName) {

  TFile *FirstFile = TFile::Open("./Templates/" + FirstFileName + ".root");
  TFile *SecondFile = TFile::Open("./Templates/" + SecondFileName + ".root");

  TCanvas *PlotCanvas = BTagPerformanceCanvas("CompareTemplates", 800.);
  PlotCanvas->Divide(1, 1); PlotCanvas->cd(1);

  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<1/*nPtRelEtaBins*/; nb++)
      for (int is = 0; is<nSystematics; is++)
	if (SystematicName[is]=="_Central") 
	  for (int dt = 0; dt<2; dt++) {

	    TString DataType = "BTagMu";
	    if (dt==1) DataType = "QCDMu";
	    
	    if (Variable!=TemplateVariable) {

	      TString ThisHistoName = HistogramName(Variable + "_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	      std::cout << DataType << std::endl;
	      TH1D *FirstHisto = (TH1D*) FirstFile->Get(ThisHistoName); FirstHisto->SetName(FirstName);
	      TH1D *SecondHisto = (TH1D*) SecondFile->Get(ThisHistoName);  SecondHisto->SetName(SecondName);

	      CompareHistograms(PlotCanvas, Variable, FirstHisto, SecondHisto, "./Plots/TemplateComparison_" + DataType + "_" + Variable + "_" + FirstName + "_" + SecondName + "_" + FitPtBin[fpt] + "_" + PtRelEtaBin[nb] + ".png");

	    } else {
	  
	      for (int tg = 0; tg<nTaggers; tg++)
		for (int tp = 0; tp<2; tp++)
		  for (int fl = 0; fl<3; fl++) {
		
		    if (dt==0 && fl!=0) continue;
		    
		    int lf = 0;
		    if (dt==1 && fl==0) lf = 1;
		    if (dt==1 && fl==1) lf = 2;
		    if (dt==1 && fl==2) lf = 3;
		    
		    TString ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, FitPtBin[fpt], nb, -1, is, tg, 10*tp + lf);
		    TH1D *FirstHisto = (TH1D*) FirstFile->Get(ThisHistoName); FirstHisto->SetName(FirstName);
		    TH1D *SecondHisto = (TH1D*) SecondFile->Get(ThisHistoName);  SecondHisto->SetName(SecondName);
		    
		    TString PlotFlag = "-" + TaggerName[tg] + "Tag";
		    if (ThisHistoName.Contains("Untag")) PlotFlag.ReplaceAll("Tag", "Untag");
		    if (ThisHistoName.Contains("_b")) PlotFlag += "_b";
		    if (ThisHistoName.Contains("_c")) PlotFlag += "_c";
		    if (ThisHistoName.Contains("_lg")) PlotFlag += "_lg";
		    

		    CompareHistograms(PlotCanvas, Variable, FirstHisto, SecondHisto, "./Plots/TemplateComparison_" + DataType + "_" + Variable + PlotFlag + "_" + FirstName + "_" + SecondName + "_" + FitPtBin[fpt] + "_" + PtRelEtaBin[nb] + ".png", 4);
		    
		  }
	      
	    }	    
	    
	  }
  
}

void PtRelAnalyzer::ComputeMethodsOverlap() {

  int const nMethods = 3;
  TString MethodName[nMethods] = {"LT", "System8", "PtRel"};

  int NumberOfSelectedEvents[nPtRelPtBins][2][2][2];

  for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
    for (int pm0 = 0; pm0<2; pm0++)
      for (int pm1 = 0; pm1<2; pm1++)
	for (int pm2 = 0; pm2<2; pm2++) 
	  NumberOfSelectedEvents[bpt][pm0][pm1][pm2] = 0;
  
  int nDataRanges = GetNumberOfDataRanges("BTagMu");
  
  for (int ph = 0; ph<nDataRanges; ph++) {
    
    TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
    GetDataRangeInformation("BTagMu", ph, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
    
    std::cout << "PtRelAnalyzer::ComputeMethodsOverlap:   Data Range Name " << DataRangeName << std::endl;
    
    for (int tf = FirstTree; tf<=nTrees; tf++) {
	    
      TString FileName = GetFileName(FileDirectoryName, tf);
      
      TFile *ThisTree = TFile::Open(FileName);
      
      std::cout << "PtRelAnalyzer::ComputeMethodsOverlap:     Tree " << tf << "/" << nTrees << std::endl;
  
      if (!ThisTree) continue;
      
      TTree *tchain = GetChain(ThisTree, false);

      Int_t nentries = (Int_t)tchain->GetEntries();
      
      for (Int_t i = 0; i<nentries; i++) {
	
	tchain->GetEntry(i);
	  
	if (nPV>0) {
      
	  int iMu; 
	  int jMu = GetPFMuonJet(&iMu);
	  
	  if (jMu>=0 && Jet_pt[jMu]>PtRelPtEdge[0] &&  Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1]) {

	    int PassMethod[3] = {0, 0, 0};
	    
	    int ptBin = -1;
	    for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	      if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	    
	    bool PassTrigger[nMethods] = {false, false, false};
	    
	    for (int trg = 0; trg<nTriggers; trg++) {
	      
	      if (Jet_pt[jMu]>MinPtJetTrigger[trg] && Jet_pt[jMu]<MaxPtJetTrigger[trg]) {
		
		if (PassTriggerBit(TriggerName[trg])) {
		  
		  PassTrigger[1] = true;
		  
		  if (PassTriggerEmulation(trg, jMu)) {
		    
		    PassTrigger[0] = true;
		    PassTrigger[2] = true;
		    
		  }
		  
		}
		
	      }
	      
	    }
	    
	    if (PassTrigger[0] || PassTrigger[1] || PassTrigger[2]) {
	      
	      bool PassAwayJet[3] = {false, false, false};
	      
	      PassAwayJet[0] = true;
	      
	      int aJet = GetAwayJet("JBPM4", jMu, 1.5, true);
	      if (aJet>=0) {
		
		if (Jet_pt[aJet]>30.)  PassAwayJet[2] = true;
		else if (Jet_pt[aJet]>20. && Jet_pt[jMu]<50.)  PassAwayJet[2] = true;
		
	      }
	      
	      int sJet = GetAwayJet("NONE", jMu, 0., false);
	      if (sJet>=0) {
		
		if (Jet_pt[sJet]>30.)  PassAwayJet[1] = true;
		else if (Jet_pt[sJet]>20. && Jet_pt[jMu]<50.)  PassAwayJet[1] = true;
		
	      }
	     
	      for (int mt = 0; mt<nMethods; mt++)
		if (PassTrigger[mt] && PassAwayJet[mt]) 
		  PassMethod[mt] = 1;
 
	    }

	    NumberOfSelectedEvents[ptBin][PassMethod[0]][PassMethod[1]][PassMethod[2]]++;

	  }

	}
	
      }
      
    }
    
  }
  
  std::cout << "PtRelAnalyzer::ComputeMethodsOverlap " << std::endl; 
       
  for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
    
    std::cout << "   Pt bin " << PtRelPtBin[bpt] << std::endl;
    
    float TotalMuonJets = 0., TotalAllMethods = 0.;
    float TotalMethod[nMethods] = {0., 0., 0.};
    float TotalOtherMethod[nMethods] = {0., 0., 0.};

    for (int pm0 = 0; pm0<2; pm0++)
      for (int pm1 = 0; pm1<2; pm1++)
	for (int pm2 = 0; pm2<2; pm2++) {

	  TotalMuonJets += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm0==1) TotalMethod[0] += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm1==1) TotalMethod[1] += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm2==1) TotalMethod[2] += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm0==1 && pm1==1) TotalOtherMethod[2] += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm0==1 && pm2==1) TotalOtherMethod[1] += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm2==1 && pm1==1) TotalOtherMethod[0] += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];
	  if (pm0==1 && pm1==1 && pm2==1) TotalAllMethods += NumberOfSelectedEvents[bpt][pm0][pm1][pm2];

	}
    
    std::cout << "         Events with a muon-jet: " << TotalMuonJets << std::endl;
    for (int mt = 0; mt<nMethods; mt++)
      std::cout << "         Events selected by " << MethodName[mt] << ": " << TotalMethod[mt] << std::endl;
    for (int mt = 0; mt<nMethods; mt++) {
      std::cout << "         Events selected by ";
      int nMethodsCounted = 0;
      for (int mt2 = 0; mt2<nMethods; mt2++)
	if (mt2!=mt) {
	  std::cout << MethodName[mt2];
	  if (nMethodsCounted<nMethods-2) std::cout << " and ";
	  nMethodsCounted++;
	}
      std::cout << ": " << TotalOtherMethod[mt] << std::endl;
    }
    std::cout << "         Events selected by all methods: " << TotalAllMethods << std::endl;
    
  }
  
}

void PtRelAnalyzer::ComputeBTaggingEfficiency(TString DataType, int DataRange) {

  TString DataRangeName, FileDirectoryName; float PtHatWeight; int nTrees, FirstTree;
  GetDataRangeInformation(DataType, DataRange, &DataRangeName, &PtHatWeight, &FileDirectoryName, &nTrees, &FirstTree);
  
  for (int tf = FirstTree; tf<=nTrees; tf++) {

    TString FileName = GetFileName(FileDirectoryName, tf);
	        
    TFile *ThisTree = TFile::Open(FileName);
    
    if (!ThisTree) continue;

    TTree *tchain = GetChain(ThisTree, false);
      
    Int_t nentries = (Int_t)tchain->GetEntries();
    
    for (Int_t i = 0; i<nentries; i++) {
      
      tchain->GetEntry(i);

      if (nPV>0) {

	for (int ijet = 0; ijet<nJet; ijet++) 
	  if (Jet_pt[ijet]>20. && fabs(Jet_eta[ijet])<2.4) {
     
	    int JetFlavour = 4;
	    if (fabs(Jet_flavour[ijet])==5) {
	      JetFlavour = 0;
	      if (IsFromGluonSplittingFromHadron(ijet, 5)) JetFlavour++;
	    } else if (fabs(Jet_flavour[ijet])==4) {
	      JetFlavour = 2;
	      if (IsFromGluonSplittingFromHadron(ijet, 4)) JetFlavour++ ;
	    }
	    
	  }
	
      }
      
    }
    
  }

}



#include "RecoBTag/PerformanceMeasurements/interface/MistagAnalyzer.h"



MistagAnalyzer::MistagAnalyzer(const edm::ParameterSet& iConfig): classifier_(iConfig)  
{
  outputFile_  = iConfig.getUntrackedParameter<std::string>("outputFile");
  rootFile_ = new TFile(outputFile_.c_str(),"RECREATE");
  
  rootFile_->cd();
  
  //now do what ever initialization is needed
  // Flavour identification
  flavourMatchOptionf = iConfig.getParameter<std::string>( "flavourMatchOption" );
  if (flavourMatchOptionf == "fastMC") {
    flavourSourcef = iConfig.getParameter<edm::InputTag>("flavourSource");
  } else if (flavourMatchOptionf == "genParticle") {
    flavourSourcef = iConfig.getParameter<edm::InputTag> ("flavourSource");
  }
  
  minJetPt_  = iConfig.getParameter<double>("MinPt");
  maxJetEta_ = iConfig.getParameter<double>("MaxEta");
  selTagger_ = iConfig.getParameter<int>("selTagger");
  tagCut_    = iConfig.getParameter<double>("tagCut");
  vetoPos_   = iConfig.getParameter<double>("vetoPos");
  ntrackMin_ = iConfig.getParameter<int>("ntrackMin");
  isData_    = iConfig.getParameter<bool>("isData");
  
  
  
  
  useTrackHistory_= iConfig.getParameter<bool>("useTrackHistory");
  
  CaloJetCollectionTags_ = iConfig.getParameter<std::string>("Jets");
  
  jetPModuleName_        = iConfig.getParameter<std::string>("jetPModuleName");
  jetPPosModuleName_     = iConfig.getParameter<std::string>("jetPPosModuleName");
  jetPNegModuleName_     = iConfig.getParameter<std::string>("jetPNegModuleName");
  
  trackCHEModuleName_    = iConfig.getParameter<std::string>("trackCHEModuleName");
  trackCNegHEModuleName_ = iConfig.getParameter<std::string>("trackCNegHEModuleName");
  
  trackCHPModuleName_    = iConfig.getParameter<std::string>("trackCHPModuleName");
  trackCNegHPModuleName_ = iConfig.getParameter<std::string>("trackCNegHPModuleName");
  
  
  combinedSvtxModuleName_     = iConfig.getParameter<std::string>("combinedSvtxModuleName");
  combinedSvtxNegModuleName_  = iConfig.getParameter<std::string>("combinedSvtxNegModuleName");
  
  svtxModuleName_       = iConfig.getParameter<std::string>("svtxModuleName");
  svtxNegModuleName_    = iConfig.getParameter<std::string>("svtxNegModuleName");
  
  softMuonModuleName_       = iConfig.getParameter<std::string>("softMuonModuleName");
  softMuonNegModuleName_    = iConfig.getParameter<std::string>("softMuonNegModuleName");
  softMuonTagInfoName_      = iConfig.getParameter<std::string>("softMuonTagInfoName");
  
  
  
  
  
  
  hData_All_NJets       = new TH1F("hData_All_NJets","nb. of jets",21,-0.5,20.5);
  hData_All_NTracks     = new TH1F("hData_All_NTracks","nb. of tracks",20,0.5,20.5);
  hData_All_JetPt       = new TH1F("hData_All_JetPt","pt(jet)",21,20.,230.);
  hData_All_JetEta      = new TH1F("hData_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_NJets           = new TH1F("hData_NJets","nb. of jets",21,-0.5,20.5);
  hData_NTracks         = new TH1F("hData_NTracks","nb. of tracks",20,0.5,20.5);
  hData_JetPt           = new TH1F("hData_JetPt","pt(jet)",21,20.,230.);
  hData_JetEta          = new TH1F("hData_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_NegTag_NTracks  = new TH1F("hData_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hData_NegTag_JetPt    = new TH1F("hData_NegTag_JetPt","pt(jet)",21,20.,230.);
  hData_NegTag_JetEta   = new TH1F("hData_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_PosTag_NTracks  = new TH1F("hData_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hData_PosTag_JetPt    = new TH1F("hData_PosTag_JetPt","pt(jet)",21,20.,230.);
  hData_PosTag_JetEta   = new TH1F("hData_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_Tagger          = new TH1F("hData_Tagger","Tagger",100,-50.,50.);
  hData_Veto_NTracks        = new TH1F("hData_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hData_Veto_JetPt          = new TH1F("hData_Veto_JetPt","pt(jet)",21,20.,230.);
  hData_Veto_JetEta         = new TH1F("hData_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_Veto_NegTag_NTracks = new TH1F("hData_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hData_Veto_NegTag_JetPt   = new TH1F("hData_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hData_Veto_NegTag_JetEta  = new TH1F("hData_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_Veto_Tagger         = new TH1F("hData_Veto_Tagger","Tagger",100,-50.,50.);
  
  //**********************************
  // All flavours in Monte Carlo
  //**********************************
  hAllFlav_All_NJets       = new TH1F("hAllFlav_All_NJets","nb. of jets",21,-0.5,20.5);
  hAllFlav_All_NTracks     = new TH1F("hAllFlav_All_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_All_JetPt       = new TH1F("hAllFlav_All_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_All_JetEta      = new TH1F("hAllFlav_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_All_Flavour     = new TH1F("hAllFlav_All_Flavour","Flavour",22,-0.5,21.5);
  hAllFlav_NJets           = new TH1F("hAllFlav_NJets","nb. of jets",21,-0.5,20.5);
  hAllFlav_NTracks         = new TH1F("hAllFlav_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_JetPt           = new TH1F("hAllFlav_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_JetEta          = new TH1F("hAllFlav_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Flavour         = new TH1F("hAllFlav_Flavour","Flavour",22,-0.5,21.5);
  hAllFlav_NegTag_NTracks  = new TH1F("hAllFlav_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_NegTag_JetPt    = new TH1F("hAllFlav_NegTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_NegTag_JetEta   = new TH1F("hAllFlav_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_PosTag          = new TH1F("hAllFlav_PosTag","Tag(+)",100,0.,25.);
  hAllFlav_PosTag_NTracks  = new TH1F("hAllFlav_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_PosTag_JetPt    = new TH1F("hAllFlav_PosTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_PosTag_JetEta   = new TH1F("hAllFlav_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Tagger          = new TH1F("hAllFlav_Tagger","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Gam      = new TH1F("hAllFlav_Tagger_Gam","Tagger",100,-50.,50.);
  hAllFlav_Tagger_K0s      = new TH1F("hAllFlav_Tagger_K0s","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Lam      = new TH1F("hAllFlav_Tagger_Lam","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Bwd      = new TH1F("hAllFlav_Tagger_Bwd","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Cwd      = new TH1F("hAllFlav_Tagger_Cwd","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Tau      = new TH1F("hAllFlav_Tagger_Tau","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Int      = new TH1F("hAllFlav_Tagger_Int","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Fak      = new TH1F("hAllFlav_Tagger_Fak","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Bad      = new TH1F("hAllFlav_Tagger_Bad","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Oth      = new TH1F("hAllFlav_Tagger_Oth","Tagger",100,-50.,50.);
  hAllFlav_Veto_NTracks        = new TH1F("hAllFlav_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Veto_JetPt          = new TH1F("hAllFlav_Veto_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Veto_JetEta         = new TH1F("hAllFlav_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Veto_Flavour        = new TH1F("hAllFlav_Veto_Flavour","Flavour",22,-0.5,21.5);
  hAllFlav_Veto_NegTag_NTracks = new TH1F("hAllFlav_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Veto_NegTag_JetPt   = new TH1F("hAllFlav_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Veto_NegTag_JetEta  = new TH1F("hAllFlav_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Veto_Tagger         = new TH1F("hAllFlav_Veto_Tagger","Tagger",100,-50.,50.);
  
  hAllFlav_K0s_NTracks             = new TH1F("hAllFlav_K0s_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_K0s_JetPt               = new TH1F("hAllFlav_K0s_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_K0s_JetEta              = new TH1F("hAllFlav_K0s_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_K0s_PosTag_NTracks      = new TH1F("hAllFlav_K0s_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_K0s_PosTag_JetPt        = new TH1F("hAllFlav_K0s_PosTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_K0s_PosTag_JetEta       = new TH1F("hAllFlav_K0s_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_K0s_Veto_NegTag_NTracks = new TH1F("hAllFlav_K0s_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_K0s_Veto_NegTag_JetPt   = new TH1F("hAllFlav_K0s_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_K0s_Veto_NegTag_JetEta  = new TH1F("hAllFlav_K0s_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Lam_NTracks             = new TH1F("hAllFlav_Lam_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Lam_JetPt               = new TH1F("hAllFlav_Lam_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Lam_JetEta              = new TH1F("hAllFlav_Lam_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Lam_PosTag_NTracks      = new TH1F("hAllFlav_Lam_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Lam_PosTag_JetPt        = new TH1F("hAllFlav_Lam_PosTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Lam_PosTag_JetEta       = new TH1F("hAllFlav_Lam_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Lam_Veto_NegTag_NTracks = new TH1F("hAllFlav_Lam_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Lam_Veto_NegTag_JetPt   = new TH1F("hAllFlav_Lam_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Lam_Veto_NegTag_JetEta  = new TH1F("hAllFlav_Lam_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Gam_NTracks             = new TH1F("hAllFlav_Gam_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Gam_JetPt               = new TH1F("hAllFlav_Gam_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Gam_JetEta              = new TH1F("hAllFlav_Gam_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Gam_PosTag_NTracks      = new TH1F("hAllFlav_Gam_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Gam_PosTag_JetPt        = new TH1F("hAllFlav_Gam_PosTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Gam_PosTag_JetEta       = new TH1F("hAllFlav_Gam_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Gam_Veto_NegTag_NTracks = new TH1F("hAllFlav_Gam_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Gam_Veto_NegTag_JetPt   = new TH1F("hAllFlav_Gam_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Gam_Veto_NegTag_JetEta  = new TH1F("hAllFlav_Gam_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Fak_NTracks             = new TH1F("hAllFlav_Fak_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Fak_JetPt               = new TH1F("hAllFlav_Fak_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Fak_JetEta              = new TH1F("hAllFlav_Fak_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Fak_PosTag_NTracks      = new TH1F("hAllFlav_Fak_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Fak_PosTag_JetPt        = new TH1F("hAllFlav_Fak_PosTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Fak_PosTag_JetEta       = new TH1F("hAllFlav_Fak_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hAllFlav_Fak_Veto_NegTag_NTracks = new TH1F("hAllFlav_Fak_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hAllFlav_Fak_Veto_NegTag_JetPt   = new TH1F("hAllFlav_Fak_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hAllFlav_Fak_Veto_NegTag_JetEta  = new TH1F("hAllFlav_Fak_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  
  //**********************************
  // udsg-jets
  //**********************************
  hLightFlav_All_NTracks    = new TH1F("hLightFlav_All_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_All_JetPt      = new TH1F("hLightFlav_All_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_All_JetEta     = new TH1F("hLightFlav_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_NTracks        = new TH1F("hLightFlav_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_JetPt          = new TH1F("hLightFlav_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_JetEta         = new TH1F("hLightFlav_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_NegTag_NTracks = new TH1F("hLightFlav_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_NegTag_JetPt   = new TH1F("hLightFlav_NegTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_NegTag_JetEta  = new TH1F("hLightFlav_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_PosTag         = new TH1F("hLightFlav_PosTag","Tag(+)",100,0.,25.);
  hLightFlav_PosTag_NTracks = new TH1F("hLightFlav_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_PosTag_JetPt   = new TH1F("hLightFlav_PosTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_PosTag_JetEta  = new TH1F("hLightFlav_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Tagger         = new TH1F("hLightFlav_Tagger","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Gam     = new TH1F("hLightFlav_Tagger_Gam","Tagger",100,-50.,50.);
  hLightFlav_Tagger_K0s     = new TH1F("hLightFlav_Tagger_K0s","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Lam     = new TH1F("hLightFlav_Tagger_Lam","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Bwd     = new TH1F("hLightFlav_Tagger_Bwd","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Cwd     = new TH1F("hLightFlav_Tagger_Cwd","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Tau     = new TH1F("hLightFlav_Tagger_Tau","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Int     = new TH1F("hLightFlav_Tagger_Int","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Fak     = new TH1F("hLightFlav_Tagger_Fak","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Bad     = new TH1F("hLightFlav_Tagger_Bad","Tagger",100,-50.,50.);
  hLightFlav_Tagger_Oth     = new TH1F("hLightFlav_Tagger_Oth","Tagger",100,-50.,50.);
  hLightFlav_Veto_NTracks        = new TH1F("hLightFlav_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Veto_JetPt          = new TH1F("hLightFlav_Veto_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Veto_JetEta         = new TH1F("hLightFlav_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Veto_NegTag_NTracks = new TH1F("hLightFlav_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Veto_NegTag_JetPt   = new TH1F("LighthFlav_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Veto_NegTag_JetEta  = new TH1F("hLightFlav_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Veto_Tagger         = new TH1F("hLightFlav_Veto_Tagger","Tagger",100,-50.,50.);
  
  hLightFlav_K0s_NTracks        = new TH1F("hLightFlav_K0s_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_K0s_JetPt          = new TH1F("hLightFlav_K0s_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_K0s_JetEta         = new TH1F("hLightFlav_K0s_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_K0s_PosTag_NTracks = new TH1F("hLightFlav_K0s_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_K0s_PosTag_JetPt   = new TH1F("hLightFlav_K0s_PosTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_K0s_PosTag_JetEta  = new TH1F("hLightFlav_K0s_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_K0s_Veto_Tagger    = new TH1F("hLightFlav_K0s_Veto_Tagger","Tagger",100,-50.,50.);
  hLightFlav_Lam_NTracks        = new TH1F("hLightFlav_Lam_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Lam_JetPt          = new TH1F("hLightFlav_Lam_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Lam_JetEta         = new TH1F("hLightFlav_Lam_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Lam_PosTag_NTracks = new TH1F("hLightFlav_Lam_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Lam_PosTag_JetPt   = new TH1F("hLightFlav_Lam_PosTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Lam_PosTag_JetEta  = new TH1F("hLightFlav_Lam_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Lam_Veto_Tagger    = new TH1F("hLightFlav_Lam_Veto_Tagger","Tagger",100,-50.,50.);
  hLightFlav_Gam_NTracks        = new TH1F("hLightFlav_Gam_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Gam_JetPt          = new TH1F("hLightFlav_Gam_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Gam_JetEta         = new TH1F("hLightFlav_Gam_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Gam_PosTag_NTracks = new TH1F("hLightFlav_Gam_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Gam_PosTag_JetPt   = new TH1F("hLightFlav_Gam_PosTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Gam_PosTag_JetEta  = new TH1F("hLightFlav_Gam_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Gam_Veto_Tagger    = new TH1F("hLightFlav_Gam_Veto_Tagger","Tagger",100,-50.,50.);
  hLightFlav_Fak_NTracks        = new TH1F("hLightFlav_Fak_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Fak_JetPt          = new TH1F("hLightFlav_Fak_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Fak_JetEta         = new TH1F("hLightFlav_Fak_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Fak_PosTag_NTracks = new TH1F("hLightFlav_Fak_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hLightFlav_Fak_PosTag_JetPt   = new TH1F("hLightFlav_Fak_PosTag_JetPt","pt(jet)",21,20.,230.);
  hLightFlav_Fak_PosTag_JetEta  = new TH1F("hLightFlav_Fak_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hLightFlav_Fak_Veto_Tagger    = new TH1F("hLightFlav_Fak_Veto_Tagger","Tagger",100,-50.,50.);

  //**********************************
  // gluon-jets
  //**********************************
  hGluonFlav_All_NTracks     = new TH1F("hGluonFlav_All_NTracks","nb. of tracks",20,0.5,20.5);
  hGluonFlav_All_JetPt       = new TH1F("hGluonFlav_All_JetPt","pt(jet)",21,20.,230.);
  hGluonFlav_All_JetEta      = new TH1F("hGluonFlav_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hGluonFlav_NTracks         = new TH1F("hGluonFlav_NTracks","nb. of tracks",20,0.5,20.5);
  hGluonFlav_JetPt           = new TH1F("hGluonFlav_JetPt","pt(jet)",21,20.,230.);
  hGluonFlav_JetEta          = new TH1F("hGluonFlav_JetEta","|#eta(jet)|",25,0.,2.5);
  hGluonFlav_PosTag          = new TH1F("hGluonFlav_PosTag","Tag(+)",100,0.,25.);
  hGluonFlav_PosTag_NTracks  = new TH1F("hGluonFlav_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hGluonFlav_PosTag_JetPt    = new TH1F("hGluonFlav_PosTag_JetPt","pt(jet)",21,20.,230.);
  hGluonFlav_PosTag_JetEta   = new TH1F("hGluonFlav_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hGluonFlav_Tagger          = new TH1F("hGluonFlav_Tagger","Tagger",100,-50.,50.);
  hGluonFlav_Veto_NTracks        = new TH1F("hGluonFlav_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hGluonFlav_Veto_JetPt          = new TH1F("hGluonFlav_Veto_JetPt","pt(jet)",21,20.,230.);
  hGluonFlav_Veto_JetEta         = new TH1F("hGluonFlav_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hGluonFlav_Veto_NegTag_NTracks = new TH1F("hGluonFlav_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hGluonFlav_Veto_NegTag_JetPt   = new TH1F("hGluonFlav_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hGluonFlav_Veto_NegTag_JetEta  = new TH1F("hGluonFlav_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hGluonFlav_Veto_Tagger         = new TH1F("hGluonFlav_Veto_Tagger","Tagger",100,-50.,50.);
  
  //**********************************
  // uds-jets
  //**********************************
  hUDSFlav_All_NTracks    = new TH1F("hUDSFlav_All_NTracks","nb. of tracks",20,0.5,20.5);
  hUDSFlav_All_JetPt      = new TH1F("hUDSFlav_All_JetPt","pt(jet)",21,20.,230.);
  hUDSFlav_All_JetEta     = new TH1F("hUDSFlav_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hUDSFlav_NTracks        = new TH1F("hUDSFlav_NTracks","nb. of tracks",20,0.5,20.5);
  hUDSFlav_JetPt          = new TH1F("hUDSFlav_JetPt","pt(jet)",21,20.,230.);
  hUDSFlav_JetEta         = new TH1F("hUDSFlav_JetEta","|#eta(jet)|",25,0.,2.5);
  hUDSFlav_PosTag         = new TH1F("hUDSFlav_PosTag","Tag(+)",100,0.,25.);
  hUDSFlav_PosTag_NTracks = new TH1F("hUDSFlav_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hUDSFlav_PosTag_JetPt   = new TH1F("hUDSFlav_PosTag_JetPt","pt(jet)",21,20.,230.);
  hUDSFlav_PosTag_JetEta  = new TH1F("hUDSFlav_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hUDSFlav_Tagger         = new TH1F("hUDSFlav_Tagger","Tagger",100,-50.,50.);
  hUDSFlav_Veto_NTracks        = new TH1F("hUDSFlav_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hUDSFlav_Veto_JetPt          = new TH1F("hUDSFlav_Veto_JetPt","pt(jet)",21,20.,230.);
  hUDSFlav_Veto_JetEta         = new TH1F("hUDSFlav_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hUDSFlav_Veto_NegTag_NTracks = new TH1F("hUDSFlav_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hUDSFlav_Veto_NegTag_JetPt   = new TH1F("hUDSFlav_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hUDSFlav_Veto_NegTag_JetEta  = new TH1F("hUDSFlav_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hUDSFlav_Veto_Tagger         = new TH1F("hUDSFlav_Veto_Tagger","Tagger",100,-50.,50.);
  
  //**********************************
  // c-jets
  //**********************************
  hCFlav_All_NTracks        = new TH1F("hCFlav_All_NTracks","nb. of tracks",20,0.5,20.5);
  hCFlav_All_JetPt          = new TH1F("hCFlav_All_JetPt","pt(jet)",21,20.,230.);
  hCFlav_All_JetEta         = new TH1F("hCFlav_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hCFlav_NTracks            = new TH1F("hCFlav_NTracks","nb. of tracks",20,0.5,20.5);
  hCFlav_JetPt              = new TH1F("hCFlav_JetPt","pt(jet)",21,20.,230.);
  hCFlav_JetEta             = new TH1F("hCFlav_JetEta","|#eta(jet)|",25,0.,2.5);
  hCFlav_PosTag             = new TH1F("hCFlav_PosTag","Tag(+)",100,0.,25.);
  hCFlav_PosTag_NTracks     = new TH1F("hCFlav_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hCFlav_PosTag_JetPt       = new TH1F("hCFlav_PosTag_JetPt","pt(jet)",21,20.,230.);
  hCFlav_PosTag_JetEta      = new TH1F("hCFlav_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hCFlav_Tagger             = new TH1F("hCFlav_Tagger","Tagger",100,-50.,50.);
  hCFlav_Veto_NTracks        = new TH1F("hCFlav_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hCFlav_Veto_JetPt          = new TH1F("hCFlav_Veto_JetPt","pt(jet)",21,20.,230.);
  hCFlav_Veto_JetEta         = new TH1F("hCFlav_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hCFlav_Veto_NegTag_NTracks = new TH1F("hCFlav_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hCFlav_Veto_NegTag_JetPt   = new TH1F("hCFlav_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hCFlav_Veto_NegTag_JetEta  = new TH1F("hCFlav_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hCFlav_Veto_Tagger         = new TH1F("hCFlav_Veto_Tagger","Tagger",100,-50.,50.);
  
  //**********************************
  // b-jets
  //**********************************
  hBFlav_All_NTracks        = new TH1F("hBFlav_All_NTracks","nb. of tracks",20,0.5,20.5);
  hBFlav_All_JetPt          = new TH1F("hBFlav_All_JetPt","pt(jet)",21,20.,230.);
  hBFlav_All_JetEta         = new TH1F("hBFlav_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hBFlav_NTracks            = new TH1F("hBFlav_NTracks","nb. of tracks",20,0.5,20.5);
  hBFlav_JetPt              = new TH1F("hBFlav_JetPt","pt(jet)",21,20.,230.);
  hBFlav_JetEta             = new TH1F("hBFlav_JetEta","|#eta(jet)|",25,0.,2.5);
  hBFlav_PosTag             = new TH1F("hBFlav_PosTag","Tag(+)",100,0.,25.);
  hBFlav_PosTag_NTracks     = new TH1F("hBFlav_PosTag_NTracks","nb. of tracks",20,0.5,20.5);
  hBFlav_PosTag_JetPt       = new TH1F("hBFlav_PosTag_JetPt","pt(jet)",21,20.,230.);
  hBFlav_PosTag_JetEta      = new TH1F("hBFlav_PosTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hBFlav_Tagger             = new TH1F("hBFlav_Tagger","Tagger",100,-50.,50.);
  hBFlav_Veto_NTracks        = new TH1F("hBFlav_Veto_NTracks","nb. of tracks",20,0.5,20.5);
  hBFlav_Veto_JetPt          = new TH1F("hBFlav_Veto_JetPt","pt(jet)",21,20.,230.);
  hBFlav_Veto_JetEta         = new TH1F("hBFlav_Veto_JetEta","|#eta(jet)|",25,0.,2.5);
  hBFlav_Veto_NegTag_NTracks = new TH1F("hBFlav_Veto_NegTag_NTracks","nb. of tracks",20,0.5,20.5);
  hBFlav_Veto_NegTag_JetPt   = new TH1F("hBFlav_Veto_NegTag_JetPt","pt(jet)",21,20.,230.);
  hBFlav_Veto_NegTag_JetEta  = new TH1F("hBFlav_Veto_NegTag_JetEta","|#eta(jet)|",25,0.,2.5);
  hBFlav_Veto_Tagger         = new TH1F("hBFlav_Veto_Tagger","Tagger",100,-50.,50.);
  
  nTuplesJets =  new TNtuple("Jets","All Jets",			     
  
"Njets:Ijet:Flavour:Ntagtracks:Ptjet:Etajet:Phijet:Pxjet:Pyjet:Pzjet:Ejet:Ip1N:Ip1P:Ip2N:Ip2P:Ip3N:Ip3P:ProbaN:ProbaP:Proba:Svtx:SvtxN:SoftM:SoftMN:CategoryN:CategoryP:CategoryJetN:CategoryJetP:CategorySVxN:CategorySVxP:CategorySoftMuonN:CategorySoftMuonP");   
			    
//"Njets:Ijet:Flavour:Ntagtracks:Ptjet:Etajet:Phijet:Pxjet:Pyjet:Pzjet:Ejet:Ip1N:Ip1P:Ip2N:Ip2P:Ip3N:Ip3P:ProbaN:ProbaP:Proba:Svtx:SvtxN:CombinedSvtx:CombinedSvtxN:SoftM:SoftMN:CategoryN:CategoryP:CategoryJetN:CategoryJetP:CategorySVxN:CategorySVxP");
  
  
}


MistagAnalyzer::~MistagAnalyzer()
{
  
  rootFile_->cd();
  //**********************************
  // Data
  //**********************************
  hData_All_NJets       ->Write();
  hData_All_NTracks     ->Write();
  hData_All_JetPt       ->Write();
  hData_All_JetEta      ->Write();
  hData_NJets           ->Write();
  hData_NTracks         ->Write();
  hData_JetPt           ->Write();
  hData_JetEta          ->Write();
  hData_NegTag_NTracks  ->Write();
  hData_NegTag_JetPt    ->Write();
  hData_NegTag_JetEta   ->Write();
  hData_PosTag_NTracks  ->Write();
  hData_PosTag_JetPt    ->Write();
  hData_PosTag_JetEta   ->Write();
  hData_Tagger          ->Write();
  hData_Veto_NTracks        ->Write();
  hData_Veto_JetPt          ->Write();
  hData_Veto_JetEta         ->Write();
  hData_Veto_NegTag_NTracks ->Write();
  hData_Veto_NegTag_JetPt   ->Write();
  hData_Veto_NegTag_JetEta  ->Write();
  hData_Veto_Tagger         ->Write();
  
  //**********************************
  // All flavours in Monte Carlo
  //**********************************
  hAllFlav_All_NJets       ->Write();
  hAllFlav_All_NTracks     ->Write();
  hAllFlav_All_JetPt       ->Write();
  hAllFlav_All_JetEta      ->Write();
  hAllFlav_All_Flavour     ->Write();
  hAllFlav_NJets           ->Write();
  hAllFlav_NTracks         ->Write();
  hAllFlav_JetPt           ->Write();
  hAllFlav_JetEta          ->Write();
  hAllFlav_Flavour         ->Write();
  hAllFlav_NegTag_NTracks  ->Write();
  hAllFlav_NegTag_JetPt    ->Write();
  hAllFlav_NegTag_JetEta   ->Write();
  hAllFlav_PosTag          ->Write();
  hAllFlav_PosTag_NTracks  ->Write();
  hAllFlav_PosTag_JetPt    ->Write();
  hAllFlav_PosTag_JetEta   ->Write();
  hAllFlav_Tagger          ->Write();
  hAllFlav_Tagger_Gam      ->Write();
  hAllFlav_Tagger_K0s      ->Write();
  hAllFlav_Tagger_Lam      ->Write();
  hAllFlav_Tagger_Bwd      ->Write();
  hAllFlav_Tagger_Cwd      ->Write();
  hAllFlav_Tagger_Tau      ->Write();
  hAllFlav_Tagger_Int      ->Write();
  hAllFlav_Tagger_Fak      ->Write();
  hAllFlav_Tagger_Bad      ->Write();
  hAllFlav_Tagger_Oth      ->Write();
  hAllFlav_Veto_NTracks        ->Write();
  hAllFlav_Veto_JetPt          ->Write();
  hAllFlav_Veto_JetEta         ->Write();
  hAllFlav_Veto_Flavour        ->Write();
  hAllFlav_Veto_NegTag_NTracks ->Write();
  hAllFlav_Veto_NegTag_JetPt   ->Write();
  hAllFlav_Veto_NegTag_JetEta  ->Write();
  hAllFlav_Veto_Tagger         ->Write();
  
  hAllFlav_K0s_NTracks        ->Write();
  hAllFlav_K0s_JetPt          ->Write();
  hAllFlav_K0s_JetEta         ->Write();
  hAllFlav_K0s_PosTag_NTracks ->Write();
  hAllFlav_K0s_PosTag_JetPt   ->Write();
  hAllFlav_K0s_PosTag_JetEta  ->Write();
  hAllFlav_K0s_Veto_NegTag_NTracks ->Write();
  hAllFlav_K0s_Veto_NegTag_JetPt   ->Write();
  hAllFlav_K0s_Veto_NegTag_JetEta  ->Write();
  hAllFlav_Lam_NTracks        ->Write();
  hAllFlav_Lam_JetPt          ->Write();
  hAllFlav_Lam_JetEta         ->Write();
  hAllFlav_Lam_PosTag_NTracks ->Write();
  hAllFlav_Lam_PosTag_JetPt   ->Write();
  hAllFlav_Lam_PosTag_JetEta  ->Write();
  hAllFlav_Lam_Veto_NegTag_NTracks ->Write();
  hAllFlav_Lam_Veto_NegTag_JetPt   ->Write();
  hAllFlav_Lam_Veto_NegTag_JetEta  ->Write();
  hAllFlav_Gam_NTracks        ->Write();
  hAllFlav_Gam_JetPt          ->Write();
  hAllFlav_Gam_JetEta         ->Write();
  hAllFlav_Gam_PosTag_NTracks ->Write();
  hAllFlav_Gam_PosTag_JetPt   ->Write();
  hAllFlav_Gam_PosTag_JetEta  ->Write();
  hAllFlav_Gam_Veto_NegTag_NTracks ->Write();
  hAllFlav_Gam_Veto_NegTag_JetPt   ->Write();
  hAllFlav_Gam_Veto_NegTag_JetEta  ->Write();
  hAllFlav_Fak_NTracks        ->Write();
  hAllFlav_Fak_JetPt          ->Write();
  hAllFlav_Fak_JetEta         ->Write();
  hAllFlav_Fak_PosTag_NTracks ->Write();
  hAllFlav_Fak_PosTag_JetPt   ->Write();
  hAllFlav_Fak_PosTag_JetEta  ->Write();
  hAllFlav_Fak_Veto_NegTag_NTracks ->Write();
  hAllFlav_Fak_Veto_NegTag_JetPt   ->Write();
  hAllFlav_Fak_Veto_NegTag_JetEta  ->Write();
  
  //**********************************
  // udsg-jets
  //**********************************
  hLightFlav_All_NTracks    ->Write();
  hLightFlav_All_JetPt      ->Write();
  hLightFlav_All_JetEta     ->Write();
  hLightFlav_NTracks        ->Write();
  hLightFlav_JetPt          ->Write();
  hLightFlav_JetEta         ->Write();
  hLightFlav_NegTag_NTracks ->Write();
  hLightFlav_NegTag_JetPt   ->Write();
  hLightFlav_NegTag_JetEta  ->Write();
  hLightFlav_PosTag         ->Write();
  hLightFlav_PosTag_NTracks ->Write();
  hLightFlav_PosTag_JetPt   ->Write();
  hLightFlav_PosTag_JetEta  ->Write();
  hLightFlav_Tagger         ->Write();
  hLightFlav_Tagger_Gam     ->Write();
  hLightFlav_Tagger_K0s     ->Write();
  hLightFlav_Tagger_Lam     ->Write();
  hLightFlav_Tagger_Bwd     ->Write();
  hLightFlav_Tagger_Cwd     ->Write();
  hLightFlav_Tagger_Tau     ->Write();
  hLightFlav_Tagger_Int     ->Write();
  hLightFlav_Tagger_Fak     ->Write();
  hLightFlav_Tagger_Bad     ->Write();
  hLightFlav_Tagger_Oth     ->Write();
  hLightFlav_Veto_NTracks        ->Write();
  hLightFlav_Veto_JetPt          ->Write();
  hLightFlav_Veto_JetEta         ->Write();
  hLightFlav_Veto_NegTag_NTracks ->Write();
  hLightFlav_Veto_NegTag_JetPt   ->Write();
  hLightFlav_Veto_NegTag_JetEta  ->Write();
  hLightFlav_Veto_Tagger         ->Write();
  
  hLightFlav_K0s_NTracks        ->Write();
  hLightFlav_K0s_JetPt          ->Write();
  hLightFlav_K0s_JetEta         ->Write();
  hLightFlav_K0s_PosTag_NTracks ->Write();
  hLightFlav_K0s_PosTag_JetPt   ->Write();
  hLightFlav_K0s_PosTag_JetEta  ->Write();
  hLightFlav_K0s_Veto_Tagger    ->Write();
  hLightFlav_Lam_NTracks        ->Write();
  hLightFlav_Lam_JetPt          ->Write();
  hLightFlav_Lam_JetEta         ->Write();
  hLightFlav_Lam_PosTag_NTracks ->Write();
  hLightFlav_Lam_PosTag_JetPt   ->Write();
  hLightFlav_Lam_PosTag_JetEta  ->Write();
  hLightFlav_Lam_Veto_Tagger    ->Write();
  hLightFlav_Gam_NTracks        ->Write();
  hLightFlav_Gam_JetPt          ->Write();
  hLightFlav_Gam_JetEta         ->Write();
  hLightFlav_Gam_PosTag_NTracks ->Write();
  hLightFlav_Gam_PosTag_JetPt   ->Write();
  hLightFlav_Gam_PosTag_JetEta  ->Write();
  hLightFlav_Gam_Veto_Tagger    ->Write();
  hLightFlav_Fak_NTracks        ->Write();
  hLightFlav_Fak_JetPt          ->Write();
  hLightFlav_Fak_JetEta         ->Write();
  hLightFlav_Fak_PosTag_NTracks ->Write();
  hLightFlav_Fak_PosTag_JetPt   ->Write();
  hLightFlav_Fak_PosTag_JetEta  ->Write();
  hLightFlav_Fak_Veto_Tagger    ->Write();
  
  //**********************************
  // gluon-jets
  //**********************************
  hGluonFlav_All_NTracks     ->Write();
  hGluonFlav_All_JetPt       ->Write();
  hGluonFlav_All_JetEta      ->Write();
  hGluonFlav_NTracks         ->Write();
  hGluonFlav_JetPt           ->Write();
  hGluonFlav_JetEta          ->Write();
  hGluonFlav_PosTag          ->Write();
  hGluonFlav_PosTag_NTracks  ->Write();
  hGluonFlav_PosTag_JetPt    ->Write();
  hGluonFlav_PosTag_JetEta   ->Write();
  hGluonFlav_Tagger          ->Write();
  hGluonFlav_Veto_NTracks        ->Write();
  hGluonFlav_Veto_JetPt          ->Write();
  hGluonFlav_Veto_JetEta         ->Write();
  hGluonFlav_Veto_NegTag_NTracks ->Write();
  hGluonFlav_Veto_NegTag_JetPt   ->Write();
  hGluonFlav_Veto_NegTag_JetEta  ->Write();
  hGluonFlav_Veto_Tagger         ->Write();
  
  //**********************************
  // uds-jets
  //**********************************
  hUDSFlav_All_NTracks    ->Write();
  hUDSFlav_All_JetPt      ->Write();
  hUDSFlav_All_JetEta     ->Write();
  hUDSFlav_NTracks        ->Write();
  hUDSFlav_JetPt          ->Write();
  hUDSFlav_JetEta         ->Write();
  hUDSFlav_PosTag         ->Write();
  hUDSFlav_PosTag_NTracks ->Write();
  hUDSFlav_PosTag_JetPt   ->Write();
  hUDSFlav_PosTag_JetEta  ->Write();
  hUDSFlav_Tagger         ->Write();
  hUDSFlav_Veto_NTracks        ->Write();
  hUDSFlav_Veto_JetPt          ->Write();
  hUDSFlav_Veto_JetEta         ->Write();
  hUDSFlav_Veto_NegTag_NTracks ->Write();
  hUDSFlav_Veto_NegTag_JetPt   ->Write();
  hUDSFlav_Veto_NegTag_JetEta  ->Write();
  hUDSFlav_Veto_Tagger         ->Write();
  
  //**********************************
  // c-jets
  //**********************************
  hCFlav_All_NTracks        ->Write();
  hCFlav_All_JetPt          ->Write();
  hCFlav_All_JetEta         ->Write();
  hCFlav_NTracks            ->Write();
  hCFlav_JetPt              ->Write();
  hCFlav_JetEta             ->Write();
  hCFlav_PosTag             ->Write();
  hCFlav_PosTag_NTracks     ->Write();
  hCFlav_PosTag_JetPt       ->Write();
  hCFlav_PosTag_JetEta      ->Write();
  hCFlav_Tagger             ->Write();
  hCFlav_Veto_NTracks        ->Write();
  hCFlav_Veto_JetPt          ->Write();
  hCFlav_Veto_JetEta         ->Write();
  hCFlav_Veto_NegTag_NTracks ->Write();
  hCFlav_Veto_NegTag_JetPt   ->Write();
  hCFlav_Veto_NegTag_JetEta  ->Write();
  hCFlav_Veto_Tagger         ->Write();
  
  //**********************************
  // b-jets
  //**********************************
  hBFlav_All_NTracks        ->Write();
  hBFlav_All_JetPt          ->Write();
  hBFlav_All_JetEta         ->Write();
  hBFlav_NTracks            ->Write();
  hBFlav_JetPt              ->Write();
  hBFlav_JetEta             ->Write();
  hBFlav_PosTag             ->Write();
  hBFlav_PosTag_NTracks     ->Write();
  hBFlav_PosTag_JetPt       ->Write();
  hBFlav_PosTag_JetEta      ->Write();
  hBFlav_Tagger             ->Write();
  hBFlav_Veto_NTracks        ->Write();
  hBFlav_Veto_JetPt          ->Write();
  hBFlav_Veto_JetEta         ->Write();
  hBFlav_Veto_NegTag_NTracks ->Write();
  hBFlav_Veto_NegTag_JetPt   ->Write();
  hBFlav_Veto_NegTag_JetEta  ->Write();
  hBFlav_Veto_Tagger         ->Write();
  
  rootFile_->Write();
  //rootFile_->Close();
  
  
}


static std::vector<std::size_t> sortedIndexes(std::vector<TrackIPTagInfo::TrackIPData > const & values)
{
  std::multimap<float, std::size_t> sortedIdx;
  std::vector<std::size_t> result;
  
  for(size_t i = 0; i < values.size(); ++i)
    sortedIdx.insert( std::pair<float, std::size_t>(values[i].ip3d.significance() , i) );
  
  for(std::multimap<float, std::size_t>::reverse_iterator it = sortedIdx.rbegin(); it != sortedIdx.rend(); ++it)
    result.push_back(it->second);
  
  return result;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MistagAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  // Tag Jets
  classifier_.newEvent(iEvent, iSetup);
  
  const JetCorrector *acorrector = JetCorrector::getJetCorrector("MCJetCorrectorIcone5",iSetup);
  // Calo Jets
  Handle<reco::CaloJetCollection> jetsColl;
  iEvent.getByLabel(CaloJetCollectionTags_, jetsColl);
  
  const reco::CaloJetCollection recoJets =   *(jetsColl.product());
  
  // initialize flavour identifiers
  edm::Handle<JetFlavourMatchingCollection> jetMC;
  
  
  Handle<std::vector<reco::TrackIPTagInfo> > tagInfo;
  iEvent.getByLabel("impactParameterTagInfos", tagInfo);
  
  Handle<std::vector<reco::SecondaryVertexTagInfo> > tagInfoSVx;
  iEvent.getByLabel("secondaryVertexTagInfos", tagInfoSVx);
  
  Handle<std::vector<reco::SecondaryVertexTagInfo> > tagInfoNegSVx;
  iEvent.getByLabel("secondaryVertexNegativeTagInfos", tagInfoNegSVx);
  
  
  
  //------------------------------------------------------
  //Jet Probability tagger
  //------------------------------------------------------
  edm::Handle<reco::JetTagCollection> jetTags_JP;
  iEvent.getByLabel(jetPModuleName_, jetTags_JP);
  
  edm::Handle<reco::JetTagCollection> jetTags_PosJP;
  iEvent.getByLabel(jetPPosModuleName_, jetTags_PosJP);
  
  edm::Handle<reco::JetTagCollection> jetTags_NegJP;
  iEvent.getByLabel(jetPNegModuleName_, jetTags_NegJP);
  
  //------------------------------------------------------
  //TrackCounting taggers
  //------------------------------------------------------
  edm::Handle<reco::JetTagCollection> jetTags_TCHighEff;
  iEvent.getByLabel(trackCHEModuleName_, jetTags_TCHighEff);
  edm::Handle<reco::JetTagCollection> jetTags_NegTCHighEff;
  iEvent.getByLabel(trackCNegHEModuleName_, jetTags_NegTCHighEff);
  
  edm::Handle<reco::JetTagCollection> jetTags_TCHighPur;
  iEvent.getByLabel(trackCHPModuleName_, jetTags_TCHighPur);
  edm::Handle<reco::JetTagCollection> jetTags_NegTCHighPur;
  iEvent.getByLabel(trackCNegHPModuleName_, jetTags_NegTCHighPur);
  
  //------------------------------------------------------
  //Secondary vertex taggers
  //------------------------------------------------------
  //edm::Handle<reco::JetTagCollection> jetTags_CombinedSvtx;
  //iEvent.getByLabel(combinedSvtxModuleName_, jetTags_CombinedSvtx);
  //edm::Handle<reco::JetTagCollection> jetTags_negCombinedSvtx;
  //iEvent.getByLabel(combinedSvtxNegModuleName_, jetTags_negCombinedSvtx);
  
  edm::Handle<reco::JetTagCollection> jetTags_Svtx;
  iEvent.getByLabel(svtxModuleName_, jetTags_Svtx);
  edm::Handle<reco::JetTagCollection> jetTags_negSvtx;
  iEvent.getByLabel(svtxNegModuleName_, jetTags_negSvtx);
  
  
  //------------------------------------------------------
  //Soft muon tagger
  //------------------------------------------------------
  
  edm::Handle<reco::JetTagCollection> jetTags_softM;
  iEvent.getByLabel(softMuonModuleName_, jetTags_softM);
  edm::Handle<reco::JetTagCollection> jetTags_softMneg;
  iEvent.getByLabel(softMuonNegModuleName_, jetTags_softMneg);
  
  edm::Handle<reco::SoftLeptonTagInfoCollection> tagInos_softmuon;
  iEvent.getByLabel(softMuonTagInfoName_, tagInos_softmuon);
  
  
  
  
  
  
  if (flavourMatchOptionf == "fastMC" && isData_) {
    iEvent.getByLabel(flavourSourcef, jetMC);
    for(JetFlavourMatchingCollection::const_iterator iter =
	  jetMC->begin(); iter != jetMC->end(); iter++)
      flavoursMapf.insert(std::pair <const edm::RefToBase<reco::Jet>, unsigned int>((*iter).first, (unsigned int)((*iter).second).getFlavour()));
  } else if (flavourMatchOptionf == "genParticle") {
    iEvent.getByLabel (flavourSourcef, theJetPartonMapf);
  }
  
  
  bool TagPos = false, TagNeg = false, Veto = false;
  float varpos, varneg;
  int ntagtracks = 0;
  int numjet = 0;
  //int allevents = 0;
  float etajet, phijet;
  //float ptjet, ejet;
  
  //float JES = 1.;
  
  
  //$$$$ int njets = recoJets.size();
  //$$$$ int nselectedJet = 0;
  //$$$$ for(unsigned int i = 0; i++;  njets){
  //$$$$   if( (recoJets)[i].pt() > 20 && fabs((recoJets)[i].eta()) < 2.5 ) nselectedJet++;
  //$$$$ }
  
  //*********************************
  // Loop over the jets
  
  int ijet = 0;
  CaloJetCollection::const_iterator jet;
  
  //$$$$
  int njets = 0;
  //$$$$
  
  for( jet = recoJets.begin(); jet != recoJets.end(); ++jet ) {
    
    
    int ith_tagged = -1;      
    //$$$$
    float Flavour  ;    
    int flavour    ; 
    if(isData_!=0){
      //      Flavour  = getMatchedParton(*jet).flavour();    
      //       flavour    = int(TMath::Abs( Flavour ));    
      flavour = abs(getMatchedParton(*jet).getFlavour());
      //$$$$
      //$$$$  int flavour    = getMatchedParton(*jet).flavour();     
      if ( flavour >= 1 && flavour <= 3 ) flavour = 1;  
    }
    
    double JES     =  acorrector->correction(*jet, iEvent, iSetup);
    
    double jetpt =   (*jet).pt()  ;
    double jeteta=   (*jet).eta() ;
    double ptjet = jetpt* JES;
    double ejet  = (*jet).energy() * JES;
    
    
    if ( (jetpt * JES) >= minJetPt_ || std::fabs( jeteta ) <= maxJetEta_ )  njets++;
    
    
    TagPos = false;
    TagNeg = false;
    Veto = false;
    etajet = TMath::Abs( (*jet).eta());
    phijet = (*jet).phi();
    if (phijet < 0.) phijet += 2*TMath::Pi();
    
    //get ntracks
    ith_tagged = this->TaggedJet(*jet,jetTags_JP);
    
    ntagtracks = (*tagInfo)[ith_tagged].probabilities(0).size();
    float Proba  = (*jetTags_JP)[ith_tagged].second;
    
    ith_tagged = this->TaggedJet(*jet,jetTags_PosJP);
    float ProbaP = (*jetTags_PosJP)[ith_tagged].second;
    
    ith_tagged = this->TaggedJet(*jet,jetTags_NegJP);
    float ProbaN = (*jetTags_NegJP)[ith_tagged].second;
    
    
    //ith_tagged              = this->TaggedJet(*jet,jetTags_CombinedSvtx);
    //float   CombinedSvtx    = (*jetTags_CombinedSvtx)[ith_tagged].second;
    //ith_tagged              = this->TaggedJet(*jet,jetTags_negCombinedSvtx);
    //float   CombinedSvtxN   = (*jetTags_negCombinedSvtx)[ith_tagged].second;
    
    ith_tagged            = this->TaggedJet(*jet,jetTags_Svtx);
    float   Svtx          = (*jetTags_Svtx)[ith_tagged].second;
    ith_tagged            = this->TaggedJet(*jet,jetTags_negSvtx);
    float   SvtxN         = (*jetTags_negSvtx)[ith_tagged].second;
    
    
    ith_tagged            = this->TaggedJet(*jet,jetTags_softM);
    float   SoftM         = (*jetTags_Svtx)[ith_tagged].second;
    ith_tagged            = this->TaggedJet(*jet,jetTags_softMneg);
    float   SoftMN        = (*jetTags_negSvtx)[ith_tagged].second;
    
    
    
    
    ith_tagged = this->TaggedJet(*jet,jetTags_TCHighEff);
    std::vector<TrackIPTagInfo::TrackIPData>  ipdata = (*tagInfo)[ith_tagged].impactParameterData();
    TrackRefVector tracks( (*tagInfo)[ith_tagged].selectedTracks() );
    std::vector<std::size_t> indexes( sortedIndexes(ipdata) );
    
    
    TrackCategories::Flags flags1P;
    TrackCategories::Flags flags2P;
    TrackCategories::Flags flags3P;
    TrackCategories::Flags flags1N;
    TrackCategories::Flags flags2N;
    TrackCategories::Flags flags3N;
    int idSize = 0;
    if(useTrackHistory_ && indexes.size() != 0 && isData_!=0){
      idSize = indexes.size();
      flags1P = classifier_.evaluate( tracks[indexes[0]] ).flags();
      if(idSize > 1) flags2P = classifier_.evaluate( tracks[indexes[1]] ).flags();
      if(idSize > 2) flags3P = classifier_.evaluate( tracks[indexes[2]] ).flags();
      flags1N = classifier_.evaluate( tracks[indexes[idSize-1]] ).flags();
      if(idSize > 1) flags2N = classifier_.evaluate( tracks[indexes[idSize-2]] ).flags();
      if(idSize > 2) flags3N = classifier_.evaluate( tracks[indexes[idSize-3]] ).flags();
    }
    
    
    //*****************************************************************
    //get track histories of trakcs in jets (for Jet Proba) 
    //*****************************************************************
    int CategoryJetP = 0;
    int CategoryJetN = 0;
    
    if(useTrackHistory_ && isData_!=0){
      ith_tagged = this->TaggedJet(*jet,jetTags_PosJP);
      TrackRefVector jetProbTracks( (*tagInfo)[ith_tagged].selectedTracks() );  
      
      for(unsigned int i=0; i<jetProbTracks.size(); i++){
	reco::TrackIPTagInfo::TrackIPData ip      = ((*tagInfo)[ith_tagged].impactParameterData())[i];
	if(ip.ip3d.significance()> 0) {
	  TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
	  if      ( theFlag[TrackCategories::Conversion] )   CategoryJetP += int(pow(10., -1 + 1)); 
	  else if ( theFlag[TrackCategories::KsDecay] )      CategoryJetP += int(pow(10., -1 + 2)); 
	  else if ( theFlag[TrackCategories::LambdaDecay] )  CategoryJetP += int(pow(10., -1 + 3)); 
	  else if ( theFlag[TrackCategories::BWeakDecay] )   CategoryJetP += int(pow(10., -1 + 4)); 
	  else if ( theFlag[TrackCategories::CWeakDecay] )   CategoryJetP += int(pow(10., -1 + 5)); 
	  else if ( theFlag[TrackCategories::TauDecay] )     CategoryJetP += int(pow(10., -1 + 6)); 
	  else if ( theFlag[TrackCategories::Interaction] )  CategoryJetP += int(pow(10., -1 + 7)); 
	  else if ( theFlag[TrackCategories::Fake] )         CategoryJetP += int(pow(10., -1 + 8)); 
	  else if ( theFlag[TrackCategories::Bad] )          CategoryJetP += int(pow(10., -1 + 9)); 
	}
	else {
	  TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
	  if      ( theFlag[TrackCategories::Conversion] )   CategoryJetN += int(pow(10., -1 + 1)); 
	  else if ( theFlag[TrackCategories::KsDecay] )      CategoryJetN += int(pow(10., -1 + 2)); 
	  else if ( theFlag[TrackCategories::LambdaDecay] )  CategoryJetN += int(pow(10., -1 + 3)); 
	  else if ( theFlag[TrackCategories::BWeakDecay] )   CategoryJetN += int(pow(10., -1 + 4)); 
	  else if ( theFlag[TrackCategories::CWeakDecay] )   CategoryJetN += int(pow(10., -1 + 5)); 
	  else if ( theFlag[TrackCategories::TauDecay] )     CategoryJetN += int(pow(10., -1 + 6)); 
	  else if ( theFlag[TrackCategories::Interaction] )  CategoryJetN += int(pow(10., -1 + 7)); 
	  else if ( theFlag[TrackCategories::Fake] )         CategoryJetN += int(pow(10., -1 + 8)); 
	  else if ( theFlag[TrackCategories::Bad] )          CategoryJetN += int(pow(10., -1 + 9)); 
	}
      }
    }
    
    
    
    //*****************************************************************
    //get track histories  associated to sec. vertex (for simple SV)
    //*****************************************************************
    int CategorySVxP = 0;
    int CategorySVxN = 0;
    
    
    if(useTrackHistory_ && isData_!=0){
      ith_tagged =    this->TaggedJet(*jet,jetTags_Svtx);
      TrackRefVector  svxPostracks( (*tagInfoSVx)[ith_tagged].vertexTracks(0) );
      for(unsigned int i=0; i<svxPostracks.size(); i++){
	TrackCategories::Flags theFlag = classifier_.evaluate( svxPostracks[i] ).flags();
	if      ( theFlag[TrackCategories::Conversion] )   CategorySVxP += int(pow(10., -1 + 1)); 
	else if ( theFlag[TrackCategories::KsDecay] )      CategorySVxP += int(pow(10., -1 + 2)); 
	else if ( theFlag[TrackCategories::LambdaDecay] )  CategorySVxP += int(pow(10., -1 + 3)); 
	else if ( theFlag[TrackCategories::BWeakDecay] )   CategorySVxP += int(pow(10., -1 + 4)); 
	else if ( theFlag[TrackCategories::CWeakDecay] )   CategorySVxP += int(pow(10., -1 + 5)); 
	else if ( theFlag[TrackCategories::TauDecay] )     CategorySVxP += int(pow(10., -1 + 6)); 
	else if ( theFlag[TrackCategories::Interaction] )  CategorySVxP += int(pow(10., -1 + 7)); 
	else if ( theFlag[TrackCategories::Fake] )         CategorySVxP += int(pow(10., -1 + 8)); 
	else if ( theFlag[TrackCategories::Bad] )          CategorySVxP += int(pow(10., -1 + 9)); 
      }
      
      
      ith_tagged  =    this->TaggedJet(*jet,jetTags_negSvtx);
      TrackRefVector   svxNegtracks( (*tagInfoNegSVx)[ith_tagged].vertexTracks(0) );
      for(unsigned int i=0; i<svxNegtracks.size(); i++){
	TrackCategories::Flags theFlag = classifier_.evaluate( svxNegtracks[i] ).flags();
	if      ( theFlag[TrackCategories::Conversion] )   CategorySVxN += int(pow(10., -1 + 1)); 
	else if ( theFlag[TrackCategories::KsDecay] )      CategorySVxN += int(pow(10., -1 + 2)); 
	else if ( theFlag[TrackCategories::LambdaDecay] )  CategorySVxN += int(pow(10., -1 + 3)); 
	else if ( theFlag[TrackCategories::BWeakDecay] )   CategorySVxN += int(pow(10., -1 + 4)); 
	else if ( theFlag[TrackCategories::CWeakDecay] )   CategorySVxN += int(pow(10., -1 + 5)); 
	else if ( theFlag[TrackCategories::TauDecay] )     CategorySVxN += int(pow(10., -1 + 6)); 
	else if ( theFlag[TrackCategories::Interaction] )  CategorySVxN += int(pow(10., -1 + 7)); 
	else if ( theFlag[TrackCategories::Fake] )         CategorySVxN += int(pow(10., -1 + 8)); 
	else if ( theFlag[TrackCategories::Bad] )          CategorySVxN += int(pow(10., -1 + 9)); 
      }
    }
    
     
    
     //*****************************************************************
    //get track histories of the muon (SoftMuon tagger)
    //*****************************************************************
    int CategorySoftMuonP   = 0;
    int CategorySoftMuonN   = 0;
    if(useTrackHistory_ && isData_!=0){
     ith_tagged = this->TaggedJet(*jet,jetTags_softM);
     
     if(  SoftM > 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ){
      TrackCategories::Flags theFlagP = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
      if      ( theFlagP[TrackCategories::Conversion] )   CategorySoftMuonP += int(pow(10., -1 + 1)); 
      else if ( theFlagP[TrackCategories::KsDecay] )      CategorySoftMuonP += int(pow(10., -1 + 2)); 
      else if ( theFlagP[TrackCategories::LambdaDecay] )  CategorySoftMuonP += int(pow(10., -1 + 3)); 
      else if ( theFlagP[TrackCategories::BWeakDecay] )   CategorySoftMuonP += int(pow(10., -1 + 4)); 
      else if ( theFlagP[TrackCategories::CWeakDecay] )   CategorySoftMuonP += int(pow(10., -1 + 5)); 
      else if ( theFlagP[TrackCategories::TauDecay] )     CategorySoftMuonP += int(pow(10., -1 + 6)); 
      else if ( theFlagP[TrackCategories::Interaction] )  CategorySoftMuonP += int(pow(10., -1 + 7)); 
      else if ( theFlagP[TrackCategories::Fake] )         CategorySoftMuonP += int(pow(10., -1 + 8)); 
      else if ( theFlagP[TrackCategories::Bad] )          CategorySoftMuonP += int(pow(10., -1 + 9)); 
     }
     ith_tagged            = this->TaggedJet(*jet,jetTags_softMneg);
     if(   -SoftMN > 1 && (*tagInos_softmuon)[ith_tagged].leptons()!=0  ){
      TrackCategories::Flags theFlagN = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
      if      ( theFlagN[TrackCategories::Conversion] )   CategorySoftMuonN += int(pow(10., -1 + 1)); 
      else if ( theFlagN[TrackCategories::KsDecay] )      CategorySoftMuonN += int(pow(10., -1 + 2)); 
      else if ( theFlagN[TrackCategories::LambdaDecay] )  CategorySoftMuonN += int(pow(10., -1 + 3)); 
      else if ( theFlagN[TrackCategories::BWeakDecay] )   CategorySoftMuonN += int(pow(10., -1 + 4)); 
      else if ( theFlagN[TrackCategories::CWeakDecay] )   CategorySoftMuonN += int(pow(10., -1 + 5)); 
      else if ( theFlagN[TrackCategories::TauDecay] )     CategorySoftMuonN += int(pow(10., -1 + 6)); 
      else if ( theFlagN[TrackCategories::Interaction] )  CategorySoftMuonN += int(pow(10., -1 + 7)); 
      else if ( theFlagN[TrackCategories::Fake] )         CategorySoftMuonN += int(pow(10., -1 + 8)); 
      else if ( theFlagN[TrackCategories::Bad] )          CategorySoftMuonN += int(pow(10., -1 + 9)); 
     }
    }
    
    
    float Ip1P   = -1000;
    float Ip1N   = 1000;
    int index = 0;
    for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
	 itipdata != ipdata.end(); itipdata++){
      double ip3D = (*itipdata).ip3d.significance() ;
      //negatif tracks
      if( ip3D > Ip1P  ) Ip1P = ip3D;         
      if( ip3D < Ip1N  ) Ip1N = ip3D;
      index++;
    }
    
    
    Ip1N = -Ip1N;
    ith_tagged = this->TaggedJet(*jet,jetTags_TCHighEff);
    float Ip2P   = (*jetTags_TCHighEff)[ith_tagged].second;
    ith_tagged = this->TaggedJet(*jet,jetTags_TCHighPur);
    float Ip3P   = (*jetTags_TCHighPur)[ith_tagged].second;
    ith_tagged = this->TaggedJet(*jet,jetTags_NegTCHighEff);
    float Ip2N   = (*jetTags_NegTCHighEff)[ith_tagged].second;
    ith_tagged = this->TaggedJet(*jet,jetTags_NegTCHighPur);
    float Ip3N   = (*jetTags_NegTCHighPur)[ith_tagged].second;
    varpos = -1000.;
    varneg = -1000.;
    
    //***************************
    // V0s
    bool K0sP = false, K0sN = false;
    bool LamP = false, LamN = false;
    bool GamP = false, GamN = false;
    bool FakP = false, FakN = false;
    
    int cat1P = 0, cat2P = 0, cat3P = 0, catP = 0;
    int cat1N = 0, cat2N = 0, cat3N = 0, catN = 0;
    
    ijet++;
    
    //*****************************************************************
    //define positive and negative tags
    //*****************************************************************
    
    
    if ( selTagger_ == 0 ) {       // jet proba
      if ( ProbaP > 0 ) varpos = 20.*ProbaP;
      if ( ProbaN > 0 ) varneg = 20.*ProbaN;
      if ( ProbaP > tagCut_ ) TagPos = true;
      if ( ProbaN > tagCut_ ) TagNeg = true;
    }
    else if ( selTagger_ == 1 ) {     // jet fisrt track
      if ( Ip1P > 0. ) varpos = Ip1P;
      if ( Ip1N > 0. ) varneg = Ip1N;
    }
    else if ( selTagger_ == 2 ) {     // TC High Eff.
      if ( Ip2P > 0. ) varpos = Ip2P;
      if ( Ip2N > 0. ) varneg = Ip2N;
    }
    else if ( selTagger_ == 3 ) {     // TC High Pure.
      if ( Ip3P > 0. ) varpos = Ip3P;
      if ( Ip3N > 0. ) varneg = Ip3N;
    }
    else if ( selTagger_ == 4 ) {    // SV simple
      if  (Svtx > 0) varpos = Svtx;
      if  (Svtx < 0) varneg = SvtxN;
      if ( Svtx  > tagCut_ ) TagPos = true;
      if (-SvtxN > tagCut_ ) TagNeg = true; 
    }
    /*else if ( selTagger_ == 5 ) {    // SV combined
    //  if  (CombinedSvtx > 0) varpos = CombinedSvtx;
    //  if  (CombinedSvtx < 0) varpos = CombinedSvtxN; 
      if  (CombinedSvtx > 0) varpos = CombinedSvtx;
      if  (CombinedSvtx < 0) varneg = CombinedSvtxN;
      if ( CombinedSvtx  > tagCut_ ) TagPos = true;
      if (-CombinedSvtxN > tagCut_ ) TagNeg = true;   
    }*/
    else if ( selTagger_ == 6 ) {    // SV combined
    //  if  (CombinedSvtx > 0) varpos = CombinedSvtx;
    //  if  (CombinedSvtx < 0) varpos = CombinedSvtxN; 
      if  ( SoftM  > 0) varpos = 5*SoftM;
      if  ( SoftMN < 0) varneg = 5*SoftM;
      if  (  SoftM > tagCut_ ) TagPos = true;
      if  (- SoftMN> tagCut_ ) TagNeg = true;   
    }
    
    
    
    if ( selTagger_ >= 1 && selTagger_ <= 3) {
      if ( varpos > tagCut_ ) TagPos = true;
      if ( varneg > tagCut_ ) TagNeg = true;
    }
    
    // veto on positive tag
    if ( Ip1P < vetoPos_ ) Veto = true;
    
    //*****************************************************************
    //get track histories for 1st, 2nd and 3rd track (TC)
    //*****************************************************************
    
    // Track history
    if (useTrackHistory_ && indexes.size()!=0 && isData_!=0) {
      flags1P[TrackCategories::Conversion] ;
      if ( flags1P[TrackCategories::Conversion] )  cat1P = 1; 
      else if ( flags1P[TrackCategories::KsDecay] )     cat1P = 2; 
      else if ( flags1P[TrackCategories::LambdaDecay] ) cat1P = 3; 
      else if ( flags1P[TrackCategories::BWeakDecay] )  cat1P = 4; 
      else if ( flags1P[TrackCategories::CWeakDecay] )  cat1P = 5; 
      else if ( flags1P[TrackCategories::TauDecay] )    cat1P = 6; 
      else if ( flags1P[TrackCategories::Interaction] ) cat1P = 7; 
      else if ( flags1P[TrackCategories::Fake] )        cat1P = 8; 
      else if ( flags1P[TrackCategories::Bad] )         cat1P = 9; 
      if(idSize > 1){
	if ( flags2P[TrackCategories::Conversion] )  cat2P = 1;
	else if ( flags2P[TrackCategories::KsDecay] )     cat2P = 2;
	else if ( flags2P[TrackCategories::LambdaDecay] ) cat2P = 3;
	else if ( flags2P[TrackCategories::BWeakDecay] )  cat2P = 4;
	else if ( flags2P[TrackCategories::CWeakDecay] )  cat2P = 5;
	else if ( flags2P[TrackCategories::TauDecay] )    cat2P = 6;
	else if ( flags2P[TrackCategories::Interaction] ) cat2P = 7;
	else if ( flags2P[TrackCategories::Fake] )        cat2P = 8;
	else if ( flags2P[TrackCategories::Bad] )         cat2P = 9;
      }
      if(idSize > 2){
	if ( flags3P[TrackCategories::Conversion] )  cat3P = 1;
	else if ( flags3P[TrackCategories::KsDecay] )     cat3P = 2;
	else if ( flags3P[TrackCategories::LambdaDecay] ) cat3P = 3;
	else if ( flags3P[TrackCategories::BWeakDecay] )  cat3P = 4;
	else if ( flags3P[TrackCategories::CWeakDecay] )  cat3P = 5;
	else if ( flags3P[TrackCategories::TauDecay] )    cat3P = 6;
	else if ( flags3P[TrackCategories::Interaction] ) cat3P = 7;
	else if ( flags3P[TrackCategories::Fake] )        cat3P = 8;
	else if ( flags3P[TrackCategories::Bad] )         cat3P = 9;
      }
      if ( flags1N[TrackCategories::Conversion] )  cat1N = 1;
      else if ( flags1N[TrackCategories::KsDecay] )     cat1N = 2;
      else if ( flags1N[TrackCategories::LambdaDecay] ) cat1N = 3;
      else if ( flags1N[TrackCategories::BWeakDecay] )  cat1N = 4;
      else if ( flags1N[TrackCategories::CWeakDecay] )  cat1N = 5;
      else if ( flags1N[TrackCategories::TauDecay] )    cat1N = 6;
      else if ( flags1N[TrackCategories::Interaction] ) cat1N = 7;
      else if ( flags1N[TrackCategories::Fake] )        cat1N = 8;
      else if ( flags1N[TrackCategories::Bad] )         cat1N = 9;
      if(idSize > 1){
	if ( flags2N[TrackCategories::Conversion] )  cat2N = 1;
	else if ( flags2N[TrackCategories::KsDecay] )     cat2N = 2;
	else if ( flags2N[TrackCategories::LambdaDecay] ) cat2N = 3;
	else if ( flags2N[TrackCategories::BWeakDecay] )  cat2N = 4;
	else if ( flags2N[TrackCategories::CWeakDecay] )  cat2N = 5;
	else if ( flags2N[TrackCategories::TauDecay] )    cat2N = 6;
	else if ( flags2N[TrackCategories::Interaction] ) cat2N = 7;
	else if ( flags2N[TrackCategories::Fake] )        cat2N = 8;
	else if ( flags2N[TrackCategories::Bad] )         cat2N = 9;
      }
      if(idSize > 2){
	if ( flags3N[TrackCategories::Conversion] )  cat3N = 1;
	else if ( flags3N[TrackCategories::KsDecay] )     cat3N = 2;
	else if ( flags3N[TrackCategories::LambdaDecay] ) cat3N = 3;
	else if ( flags3N[TrackCategories::BWeakDecay] )  cat3N = 4;
	else if ( flags3N[TrackCategories::CWeakDecay] )  cat3N = 5;
	else if ( flags3N[TrackCategories::TauDecay] )    cat3N = 6;
	else if ( flags3N[TrackCategories::Interaction] ) cat3N = 7;
	else if ( flags3N[TrackCategories::Fake] )        cat3N = 8;
	else if ( flags3N[TrackCategories::Bad] )         cat3N = 9;
      }
      
      
      if ( selTagger_ == 0 || selTagger_ == 3 ) {
	catP = cat3P;
	catN = cat3N;
	if ( cat1P == 1 || cat2P == 1 || cat3P == 1 ) GamP = true;
	if ( cat1N == 1 || cat2N == 1 || cat3N == 1 ) GamN = true;
	if ( cat1P == 2 || cat2P == 2 || cat3P == 2 ) K0sP = true;
	if ( cat1N == 2 || cat2N == 2 || cat3N == 2 ) K0sN = true;
	if ( cat1P == 3 || cat2P == 3 || cat3P == 3 ) LamP = true;
	if ( cat1N == 3 || cat2N == 3 || cat3N == 3 ) LamN = true;
	if ( cat1P == 8 || cat2P == 8 || cat3P == 8 ) FakP = true;
	if ( cat1N == 8 || cat2N == 8 || cat3N == 8 ) FakN = true;
      }
      else if ( selTagger_ == 2 ) {
	catP = cat2P;
	catN = cat2N;
	if ( cat1P == 1 || cat2P == 1 ) GamP = true;
	if ( cat1N == 1 || cat2N == 1 ) GamN = true;
	if ( cat1P == 2 || cat2P == 2 ) K0sP = true;
	if ( cat1N == 2 || cat2N == 2 ) K0sN = true;
	if ( cat1P == 3 || cat2P == 3 ) LamP = true;
	if ( cat1N == 3 || cat2N == 3 ) LamN = true;
	if ( cat1P == 8 || cat2P == 8 ) FakP = true;
	if ( cat1N == 8 || cat2N == 8 ) FakN = true;
      }
      else if ( selTagger_ == 1 ) {
	catP = cat1P;
	catN = cat1N;
	if ( cat1P == 1 ) GamP = true;
	if ( cat1N == 1 ) GamN = true;
	if ( cat1P == 2 ) K0sP = true;
	if ( cat1N == 2 ) K0sN = true;
	if ( cat1P == 3 ) LamP = true;
	if ( cat1N == 3 ) LamN = true;
	if ( cat1P == 8 ) FakP = true;
	if ( cat1N == 8 ) FakN = true;
      }
    }
    
    //*********************************
    // Jet selection
    
    if ( (jetpt * JES) <= minJetPt_ || std::fabs( jeteta ) >= maxJetEta_ ) continue;
    numjet++;
    
    float Njets =      njets;
    float Ijet =       numjet;
    //$$$$  float Flavour =    flavour;
    float Ntagtracks = ntagtracks;
    
    
    float CategoryP = 10*cat1P + 100*cat2P + 1000*cat3P;
    float CategoryN = 10*cat1N + 100*cat2N + 1000*cat3N;
    
    
    
    //*****************************************************************
    //fill the ntuple
    //*****************************************************************
    
    float jet_input[32] =
      {Njets, Ijet, Flavour, Ntagtracks, ptjet, (*jet).eta(), phijet, 
       (*jet).px() * JES, (*jet).py() * JES, (*jet).pz() * JES, ejet,
       //$$$$
       Ip1N, Ip1P, Ip2N, Ip2P, Ip3N, Ip3P, ProbaN, ProbaP, Proba,
       //Svtx, SvtxN, CombinedSvtx, CombinedSvtxN,
       Svtx, SvtxN,
       SoftM, SoftMN,
       CategoryN, CategoryP,
       CategoryJetN, CategoryJetP,
       CategorySVxN, CategorySVxP,
       CategorySoftMuonN, CategorySoftMuonP
      };
    
    
    nTuplesJets->Fill(jet_input);
    
    //*****************************************************************
    //fill the histograms
    //*****************************************************************
    
    
    
    
    hData_All_NTracks->Fill( ntagtracks ) ;
    hData_All_JetPt->Fill( ptjet );
    hData_All_JetEta->Fill( etajet );
    
    hAllFlav_All_NTracks->Fill( ntagtracks ) ;
    hAllFlav_All_JetPt->Fill( ptjet );
    hAllFlav_All_JetEta->Fill( etajet );
    hAllFlav_All_Flavour->Fill( flavour );
    
    if(isData_!=0){
      if (flavour == 1 || flavour == 21  ) {
	hLightFlav_All_NTracks->Fill( ntagtracks );
	hLightFlav_All_JetPt->Fill( ptjet );
	hLightFlav_All_JetEta->Fill( etajet );
      }
      
      if (flavour == 21) {
	hGluonFlav_All_NTracks->Fill( ntagtracks );
	hGluonFlav_All_JetPt->Fill( ptjet );
	hGluonFlav_All_JetEta->Fill( etajet );
      }
      else if (flavour == 1) {
	hUDSFlav_All_NTracks->Fill( ntagtracks );
	hUDSFlav_All_JetPt->Fill( ptjet );
	hUDSFlav_All_JetEta->Fill( etajet );
      }
      else if (flavour == 4) {
	hCFlav_All_NTracks->Fill( ntagtracks );
	hCFlav_All_JetPt->Fill( ptjet );
	hCFlav_All_JetEta->Fill( etajet );
      }
      else if (flavour == 5) {
	hBFlav_All_NTracks->Fill( ntagtracks );
	hBFlav_All_JetPt->Fill( ptjet );
	hBFlav_All_JetEta->Fill( etajet );
      }
    }
    //*********************************
    // Taggability
    //$$
    if ( ntagtracks < ntrackMin_ ) continue;
    //$$
    
    hData_NTracks->Fill( ntagtracks ) ;
    hData_JetPt->Fill( ptjet );
    hData_JetEta->Fill( etajet );
    
    if(isData_!=0){
      hAllFlav_NTracks->Fill( ntagtracks ) ;
      hAllFlav_JetPt->Fill( ptjet );
      hAllFlav_JetEta->Fill( etajet );
      hAllFlav_Flavour->Fill( flavour );
      if ( K0sP || K0sN ) {
	hAllFlav_K0s_NTracks->Fill( ntagtracks ) ;
	hAllFlav_K0s_JetPt->Fill( ptjet );
	hAllFlav_K0s_JetEta->Fill( etajet );
      }
      if ( LamP || LamN ) {
	hAllFlav_Lam_NTracks->Fill( ntagtracks ) ;
	hAllFlav_Lam_JetPt->Fill( ptjet );
	hAllFlav_Lam_JetEta->Fill( etajet );
      }
      if ( GamP || GamN ) {
	hAllFlav_Gam_NTracks->Fill( ntagtracks ) ;
	hAllFlav_Gam_JetPt->Fill( ptjet );
	hAllFlav_Gam_JetEta->Fill( etajet );
      }
      if ( FakP || FakN ) {
	hAllFlav_Fak_NTracks->Fill( ntagtracks ) ;
	hAllFlav_Fak_JetPt->Fill( ptjet );
	hAllFlav_Fak_JetEta->Fill( etajet );
      }
      
      if (flavour == 1 || flavour == 21  ) {
	hLightFlav_NTracks->Fill( ntagtracks );
	hLightFlav_JetPt->Fill( ptjet );
	hLightFlav_JetEta->Fill( etajet );
	if ( K0sP || K0sN ) {
	  hLightFlav_K0s_NTracks->Fill( ntagtracks ) ;
	  hLightFlav_K0s_JetPt->Fill( ptjet );
	  hLightFlav_K0s_JetEta->Fill( etajet );
	}
	if ( LamP || LamN ) {
	  hLightFlav_Lam_NTracks->Fill( ntagtracks ) ;
	  hLightFlav_Lam_JetPt->Fill( ptjet );
	  hLightFlav_Lam_JetEta->Fill( etajet );
	}
	if ( GamP || GamN ) {
	  hLightFlav_Gam_NTracks->Fill( ntagtracks ) ;
	  hLightFlav_Gam_JetPt->Fill( ptjet );
	  hLightFlav_Gam_JetEta->Fill( etajet );
	}
	if ( FakP || FakN ) {
	  hLightFlav_Fak_NTracks->Fill( ntagtracks ) ;
	  hLightFlav_Fak_JetPt->Fill( ptjet );
	  hLightFlav_Fak_JetEta->Fill( etajet );
	}
      }
      
      if (flavour == 21) {
	hGluonFlav_NTracks->Fill( ntagtracks );
	hGluonFlav_JetPt->Fill( ptjet );
	hGluonFlav_JetEta->Fill( etajet );
      }
      else if (flavour == 1) {
	hUDSFlav_NTracks->Fill( ntagtracks );
	hUDSFlav_JetPt->Fill( ptjet );
	hUDSFlav_JetEta->Fill( etajet );
      }
      else if (flavour == 4) {
	hCFlav_NTracks->Fill( ntagtracks );
	hCFlav_JetPt->Fill( ptjet );
	hCFlav_JetEta->Fill( etajet );
      }
      else if (flavour == 5) {
	hBFlav_NTracks->Fill( ntagtracks );
	hBFlav_JetPt->Fill( ptjet );
	hBFlav_JetEta->Fill( etajet );
      }
    }
    //*********************************
    // Tagging
    
    if ( TagNeg ) {
      hData_NegTag_NTracks->Fill( ntagtracks );
      hData_NegTag_JetPt->Fill( ptjet );
      hData_NegTag_JetEta->Fill( etajet );
    }
    if ( TagPos ) {
      hData_PosTag_NTracks->Fill( ntagtracks );
      hData_PosTag_JetPt->Fill( ptjet );
      hData_PosTag_JetEta->Fill( etajet );
    }
    if ( varneg > 0 ) hData_Tagger->Fill(-varneg );
    if ( varpos > 0 ) {
      hData_Tagger->Fill( varpos );
      hData_Veto_Tagger->Fill( varpos );
    }
    if(isData_!=0){
      if ( TagNeg ) {
	hAllFlav_NegTag_NTracks->Fill( ntagtracks );
	hAllFlav_NegTag_JetPt->Fill( ptjet );
	hAllFlav_NegTag_JetEta->Fill( etajet );
      }
      if ( varpos > 0 ) hAllFlav_PosTag->Fill( varpos );
      if ( TagPos ) {
	hAllFlav_PosTag_NTracks->Fill( ntagtracks );
	hAllFlav_PosTag_JetPt->Fill( ptjet );
	hAllFlav_PosTag_JetEta->Fill( etajet );
	if ( K0sP ) {
	  hAllFlav_K0s_PosTag_NTracks->Fill( ntagtracks );
	  hAllFlav_K0s_PosTag_JetPt->Fill( ptjet );
	  hAllFlav_K0s_PosTag_JetEta->Fill( etajet );
	}
	if ( LamP ) {
	  hAllFlav_Lam_PosTag_NTracks->Fill( ntagtracks );
	  hAllFlav_Lam_PosTag_JetPt->Fill( ptjet );
	  hAllFlav_Lam_PosTag_JetEta->Fill( etajet );
	}
	if ( GamP ) {
	  hAllFlav_Gam_PosTag_NTracks->Fill( ntagtracks );
	  hAllFlav_Gam_PosTag_JetPt->Fill( ptjet );
	  hAllFlav_Gam_PosTag_JetEta->Fill( etajet );
	}
	if ( FakP ) {
	  hAllFlav_Fak_PosTag_NTracks->Fill( ntagtracks );
	  hAllFlav_Fak_PosTag_JetPt->Fill( ptjet );
	  hAllFlav_Fak_PosTag_JetEta->Fill( etajet );
	}
      }
      if ( varneg > 0 ) hAllFlav_Tagger->Fill(-varneg );
      if ( varpos > 0 ) {
	hAllFlav_Tagger->Fill( varpos );
	hAllFlav_Veto_Tagger->Fill( varpos );
      }
      
      if ( varneg > 0 ) {
	if      ( catN == 1 ) hAllFlav_Tagger_Gam->Fill(-varneg );
	else if ( catN == 2 ) hAllFlav_Tagger_K0s->Fill(-varneg );
	else if ( catN == 3 ) hAllFlav_Tagger_Lam->Fill(-varneg );
	else if ( catN == 4 ) hAllFlav_Tagger_Bwd->Fill(-varneg );
	else if ( catN == 5 ) hAllFlav_Tagger_Cwd->Fill(-varneg );
	else if ( catN == 6 ) hAllFlav_Tagger_Tau->Fill(-varneg );
	else if ( catN == 7 ) hAllFlav_Tagger_Int->Fill(-varneg );
	else if ( catN == 8 ) hAllFlav_Tagger_Fak->Fill(-varneg );
	else if ( catN == 9 ) hAllFlav_Tagger_Bad->Fill(-varneg );
	else             hAllFlav_Tagger_Oth->Fill(-varneg );
      }
      if ( varpos > 0 ) {
	if      ( catP == 1 ) hAllFlav_Tagger_Gam->Fill( varpos );
	else if ( catP == 2 ) hAllFlav_Tagger_K0s->Fill( varpos );
	else if ( catP == 3 ) hAllFlav_Tagger_Lam->Fill( varpos );
	else if ( catP == 4 ) hAllFlav_Tagger_Bwd->Fill( varpos );
	else if ( catP == 5 ) hAllFlav_Tagger_Cwd->Fill( varpos );
	else if ( catP == 6 ) hAllFlav_Tagger_Tau->Fill( varpos );
	else if ( catP == 7 ) hAllFlav_Tagger_Int->Fill( varpos );
	else if ( catP == 8 ) hAllFlav_Tagger_Fak->Fill( varpos );
	else if ( catP == 9 ) hAllFlav_Tagger_Bad->Fill( varpos );
	else             hAllFlav_Tagger_Oth->Fill( varpos );
      }
      
      if ( flavour == 1 || flavour == 21 ) { // light q+gluon
	if ( TagNeg ) {
	  hLightFlav_NegTag_NTracks->Fill( ntagtracks );
	  hLightFlav_NegTag_JetPt->Fill( ptjet );
	  hLightFlav_NegTag_JetEta->Fill( etajet );
	}
	if ( varpos > 0 ) hLightFlav_PosTag->Fill( varpos );
	if ( TagPos ) {
	  hLightFlav_PosTag_NTracks->Fill( ntagtracks );
	  hLightFlav_PosTag_JetPt->Fill( ptjet );
	  hLightFlav_PosTag_JetEta->Fill( etajet );
	  if ( K0sP ) {
	    hLightFlav_K0s_PosTag_NTracks->Fill( ntagtracks );
	    hLightFlav_K0s_PosTag_JetPt->Fill( ptjet );
	    hLightFlav_K0s_PosTag_JetEta->Fill( etajet );
	  }
	  if ( LamP ) {
	    hLightFlav_Lam_PosTag_NTracks->Fill( ntagtracks );
	    hLightFlav_Lam_PosTag_JetPt->Fill( ptjet );
	    hLightFlav_Lam_PosTag_JetEta->Fill( etajet );
	  }
	  if ( GamP ) {
	    hLightFlav_Gam_PosTag_NTracks->Fill( ntagtracks );
	    hLightFlav_Gam_PosTag_JetPt->Fill( ptjet );
	    hLightFlav_Gam_PosTag_JetEta->Fill( etajet );
	  }
	  if ( FakP ) {
	    hLightFlav_Fak_PosTag_NTracks->Fill( ntagtracks );
	    hLightFlav_Fak_PosTag_JetPt->Fill( ptjet );
	    hLightFlav_Fak_PosTag_JetEta->Fill( etajet );
	  }
	}
	if ( varneg > 0 ) hLightFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) {
	  hLightFlav_Tagger->Fill( varpos );
	  hLightFlav_Veto_Tagger->Fill( varpos );
	  if ( K0sP ) hLightFlav_K0s_Veto_Tagger->Fill( varpos );
	  if ( LamP ) hLightFlav_Lam_Veto_Tagger->Fill( varpos );
	  if ( GamP ) hLightFlav_Gam_Veto_Tagger->Fill( varpos );
	  if ( FakP ) hLightFlav_Fak_Veto_Tagger->Fill( varpos );
	}
	if ( varneg > 0 ) {
	  if     ( catN == 1 ) hLightFlav_Tagger_Gam->Fill(-varneg );
	  else if ( catN == 2 ) hLightFlav_Tagger_K0s->Fill(-varneg );
	  else if ( catN == 3 ) hLightFlav_Tagger_Lam->Fill(-varneg );
	  else if ( catN == 4 ) hLightFlav_Tagger_Bwd->Fill(-varneg );
	  else if ( catN == 5 ) hLightFlav_Tagger_Cwd->Fill(-varneg );
	  else if ( catN == 6 ) hLightFlav_Tagger_Tau->Fill(-varneg );
	  else if ( catN == 7 ) hLightFlav_Tagger_Int->Fill(-varneg );
	  else if ( catN == 8 ) hLightFlav_Tagger_Fak->Fill(-varneg );
	  else if ( catN == 9 ) hLightFlav_Tagger_Bad->Fill(-varneg );
	  else               hLightFlav_Tagger_Oth->Fill(-varneg );
	}
	if ( varpos > 0 ) {
	  if     ( catP == 1 ) hLightFlav_Tagger_Gam->Fill( varpos );
	  else if ( catP == 2 ) hLightFlav_Tagger_K0s->Fill( varpos );
	  else if ( catP == 4 ) hLightFlav_Tagger_Bwd->Fill( varpos );
	  else if ( catP == 3 ) hLightFlav_Tagger_Lam->Fill( varpos );
	  else if ( catP == 4 ) hLightFlav_Tagger_Bwd->Fill( varpos );
	  else if ( catP == 5 ) hLightFlav_Tagger_Cwd->Fill( varpos );
	  else if ( catP == 6 ) hLightFlav_Tagger_Tau->Fill( varpos );
	  else if ( catP == 7 ) hLightFlav_Tagger_Int->Fill( varpos );
	  else if ( catP == 8 ) hLightFlav_Tagger_Fak->Fill( varpos );
	  else if ( catP == 9 ) hLightFlav_Tagger_Bad->Fill( varpos );
	  else               hLightFlav_Tagger_Oth->Fill( varpos );
	}
      } // light q+gluon
      
      if ( flavour == 21 ) { // gluon jets
	if ( varpos > 0 ) hGluonFlav_PosTag->Fill( varpos );
	if ( TagPos ) {
	  hGluonFlav_PosTag_NTracks->Fill( ntagtracks );
	  hGluonFlav_PosTag_JetPt->Fill( ptjet );
	  hGluonFlav_PosTag_JetEta->Fill( etajet );
	}
	if ( varneg > 0 ) hGluonFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hGluonFlav_Tagger->Fill( varpos );
	if ( varpos > 0 ) hGluonFlav_Veto_Tagger->Fill( varpos );
      } // gluon jets
      
      else if ( flavour == 1 ) { // uds jets
	if ( varpos > 0 ) hUDSFlav_PosTag->Fill( varpos );
	if ( TagPos ) {
	  hUDSFlav_PosTag_NTracks->Fill( ntagtracks );
	  hUDSFlav_PosTag_JetPt->Fill( ptjet );
	  hUDSFlav_PosTag_JetEta->Fill( etajet );
	}
	if ( varneg > 0 ) hUDSFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hUDSFlav_Tagger->Fill( varpos );
	if ( varpos > 0 ) hUDSFlav_Veto_Tagger->Fill( varpos );
      } // uds jets
      
      else if ( flavour == 4 ) { // c jets
	if ( varpos > 0 ) hCFlav_PosTag->Fill( varpos );
	if ( TagPos ) {
	  hCFlav_PosTag_NTracks->Fill( ntagtracks );
	  hCFlav_PosTag_JetPt->Fill( ptjet );
	  hCFlav_PosTag_JetEta->Fill( etajet );
	}
	if ( varneg > 0 ) hCFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hCFlav_Tagger->Fill( varpos );
	if ( varpos > 0 ) hCFlav_Veto_Tagger->Fill( varpos );
      } // c jets
      
      else if ( flavour == 5 ) { // b jets
	if ( varpos > 0 ) hBFlav_PosTag->Fill( varpos );
	if ( TagPos ) {
	  hBFlav_PosTag_NTracks->Fill( ntagtracks );
	  hBFlav_PosTag_JetPt->Fill( ptjet );
	  hBFlav_PosTag_JetEta->Fill( etajet );
	}
	if ( varneg > 0 ) hBFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hBFlav_Tagger->Fill( varpos );
	if ( varpos > 0 ) hBFlav_Veto_Tagger->Fill( varpos );
      } // b jets
    }
    //*********************************
    // Veto on Positive Tag
    //$$
    if ( !Veto ) continue;
    //$$
    
    hData_Veto_NTracks->Fill( ntagtracks ) ;
    hData_Veto_JetPt->Fill( ptjet );
    hData_Veto_JetEta->Fill( etajet );
    if ( varneg > 0 ) hData_Veto_Tagger->Fill(-varneg );
    
    if(isData_!=0){
      hAllFlav_Veto_NTracks->Fill( ntagtracks ) ;
      hAllFlav_Veto_JetPt->Fill( ptjet );
      hAllFlav_Veto_JetEta->Fill( etajet );
      hAllFlav_Veto_Flavour->Fill( flavour );
      if ( varneg > 0 ) hAllFlav_Veto_Tagger->Fill(-varneg );
      
      if (flavour == 1 || flavour == 21  ) {
	hLightFlav_Veto_NTracks->Fill( ntagtracks );
	hLightFlav_Veto_JetPt->Fill( ptjet );
	hLightFlav_Veto_JetEta->Fill( etajet );
	if ( varneg > 0 ) {
	  hLightFlav_Veto_Tagger->Fill(-varneg );
	  if ( K0sN ) hLightFlav_K0s_Veto_Tagger->Fill(-varneg );
	  if ( LamN ) hLightFlav_Lam_Veto_Tagger->Fill(-varneg );
	  if ( GamN ) hLightFlav_Gam_Veto_Tagger->Fill(-varneg );
	  if ( FakN ) hLightFlav_Fak_Veto_Tagger->Fill(-varneg );
	}
      }
      
      if (flavour == 21) {
	hGluonFlav_Veto_NTracks->Fill( ntagtracks );
	hGluonFlav_Veto_JetPt->Fill( ptjet );
	hGluonFlav_Veto_JetEta->Fill( etajet );
	if ( varneg > 0 ) hGluonFlav_Veto_Tagger->Fill(-varneg );
      }
      else if (flavour == 1) {
	hUDSFlav_Veto_NTracks->Fill( ntagtracks );
	hUDSFlav_Veto_JetPt->Fill( ptjet );
	hUDSFlav_Veto_JetEta->Fill( etajet );
	if ( varneg > 0 ) hUDSFlav_Veto_Tagger->Fill(-varneg );
      }
      else if (flavour == 4) {
	hCFlav_Veto_NTracks->Fill( ntagtracks );
	hCFlav_Veto_JetPt->Fill( ptjet );
	hCFlav_Veto_JetEta->Fill( etajet );
	if ( varneg > 0 ) hCFlav_Veto_Tagger->Fill(-varneg );
      }
      else if (flavour == 5) {
	hBFlav_Veto_NTracks->Fill( ntagtracks );
	hBFlav_Veto_JetPt->Fill( ptjet );
	hBFlav_Veto_JetEta->Fill( etajet );
	if ( varneg > 0 ) hBFlav_Veto_Tagger->Fill(-varneg );
      }
    }
    //*********************************
    // Veto and Negative Tag
    //$$
    if ( !TagNeg ) continue;
    //$$
    
    hData_Veto_NegTag_NTracks->Fill( ntagtracks ) ;
    hData_Veto_NegTag_JetPt->Fill( ptjet );
    hData_Veto_NegTag_JetEta->Fill( etajet );
    
    if(isData_!=0){
      hAllFlav_Veto_NegTag_NTracks->Fill( ntagtracks ) ;
      hAllFlav_Veto_NegTag_JetPt->Fill( ptjet );
      hAllFlav_Veto_NegTag_JetEta->Fill( etajet );
      if ( K0sN ) {
	hAllFlav_K0s_Veto_NegTag_NTracks->Fill( ntagtracks );
	hAllFlav_K0s_Veto_NegTag_JetPt->Fill( ptjet );
	hAllFlav_K0s_Veto_NegTag_JetEta->Fill( etajet );
      }
      if ( LamN ) {
	hAllFlav_Lam_Veto_NegTag_NTracks->Fill( ntagtracks );
	hAllFlav_Lam_Veto_NegTag_JetPt->Fill( ptjet );
	hAllFlav_Lam_Veto_NegTag_JetEta->Fill( etajet );
      }
      if ( GamN ) {
	hAllFlav_Gam_Veto_NegTag_NTracks->Fill( ntagtracks );
	hAllFlav_Gam_Veto_NegTag_JetPt->Fill( ptjet );
	hAllFlav_Gam_Veto_NegTag_JetEta->Fill( etajet );
      }
      if ( FakN ) {
	hAllFlav_Fak_Veto_NegTag_NTracks->Fill( ntagtracks );
	hAllFlav_Fak_Veto_NegTag_JetPt->Fill( ptjet );
	hAllFlav_Fak_Veto_NegTag_JetEta->Fill( etajet );
      }
      
      if (flavour == 1 || flavour == 21  ) {
	hLightFlav_Veto_NegTag_NTracks->Fill( ntagtracks );
	hLightFlav_Veto_NegTag_JetPt->Fill( ptjet );
	hLightFlav_Veto_NegTag_JetEta->Fill( etajet );
      }
      
      if (flavour == 21) {
	hGluonFlav_Veto_NegTag_NTracks->Fill( ntagtracks );
	hGluonFlav_Veto_NegTag_JetPt->Fill( ptjet );
	hGluonFlav_Veto_NegTag_JetEta->Fill( etajet );
      }
      else if (flavour == 1) {
	hUDSFlav_Veto_NegTag_NTracks->Fill( ntagtracks );
	hUDSFlav_Veto_NegTag_JetPt->Fill( ptjet );
	hUDSFlav_Veto_NegTag_JetEta->Fill( etajet );
      }
      else if (flavour == 4) {
	hCFlav_Veto_NegTag_NTracks->Fill( ntagtracks );
	hCFlav_Veto_NegTag_JetPt->Fill( ptjet );
	hCFlav_Veto_NegTag_JetEta->Fill( etajet );
      }
      else if (flavour == 5) {
	hBFlav_Veto_NegTag_NTracks->Fill( ntagtracks );
	hBFlav_Veto_NegTag_JetPt->Fill( ptjet );
	hBFlav_Veto_NegTag_JetEta->Fill( etajet );
      }
    }
    
  } // end loop on jet
  
  hData_All_NJets->Fill( njets );
  hData_NJets->Fill( numjet );
  
  if(isData_!=0){
    hAllFlav_All_NJets->Fill( njets );
    hAllFlav_NJets->Fill( numjet );
  }
}


float 
MistagAnalyzer::calculPtRel()
{
  float threturn = 0;
  
  
  return threturn;

}



// ------------ method called once each job just before starting event loop  ------------
void
MistagAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MistagAnalyzer::endJob() {
}
reco::JetFlavour MistagAnalyzer::getMatchedParton(const reco::CaloJet &jet)
{
  reco::JetFlavour jetFlavour;
  
  if (flavourMatchOptionf == "fastMC") {
    
    //    jetFlavour.underlyingParton4Vec(jet.p4());

    // DISABLED

    
  } else if (flavourMatchOptionf == "genParticle") {
    for ( JetFlavourMatchingCollection::const_iterator j  = theJetPartonMapf->begin();
	  j != theJetPartonMapf->end();
	  j ++ ) {
      RefToBase<Jet> aJet  = (*j).first;       //      const JetFlavour aFlav = (*j).second;
      if ( fabs(aJet->phi() - jet.phi()) < 1.e-5 && fabs(aJet->eta() - jet.eta())< 1.e-5 ){
	// matched
		jetFlavour = reco::JetFlavour (aJet->p4(), math::XYZPoint(0,0,0), (*j).second.getFlavour());
		//	jetFlavour.flavour((*j).second.getFlavour());
      }
    }
    
    return jetFlavour;
  }
  
  return jetFlavour;
}









int
MistagAnalyzer::TaggedJet(reco::CaloJet calojet, edm::Handle<reco::JetTagCollection > jetTags ) {
  double small = 1.e-5;
  int result = -1;
  
  
  //std::string moduleLabel = (jetTags).provenance()->moduleLabel();
  //std::string processname = (jetTags).provenance()->processName();
  for (size_t t = 0; t < jetTags->size(); t++) {
    edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].first;
    if (jet_p.isNull()) {
      continue;
    }
    if (DeltaR<reco::Candidate>()( calojet, *jet_p ) < small) {
      result = (int) t;
    }
  }
  return result;
}





//define this as a plug-in
DEFINE_FWK_MODULE(MistagAnalyzer);


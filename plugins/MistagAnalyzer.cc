
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
  
  //svComputer = new CombinedSVComputer(iConfig.getParameter<edm::ParameterSet>("csvComputerPSet"));
  //svComputer(iConfig.getParameter<edm::ParameterSet>("csvComputerPSet"));
 
  useTrackHistory_       = iConfig.getParameter<bool>("useTrackHistory");
  jetCorrector_          = iConfig.getParameter<std::string>("jetCorrector");
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
  
///////////////
// Some Histograms

  hData_All_NJets       = new TH1F("hData_All_NJets","nb. of jets",21,-0.5,20.5);
  hData_All_NTracks     = new TH1F("hData_All_NTracks","nb. of tracks",20,0.5,20.5);
  hData_All_JetPt       = new TH1F("hData_All_JetPt","pt(jet)",21,20.,230.);
  hData_All_JetEta      = new TH1F("hData_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_NJets           = new TH1F("hData_NJets","nb. of jets",21,-0.5,20.5);
  hData_NTracks         = new TH1F("hData_NTracks","nb. of tracks",20,0.5,20.5);
  hData_JetPt           = new TH1F("hData_JetPt","pt(jet)",21,20.,230.);
  hData_JetEta          = new TH1F("hData_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_Tagger          = new TH1F("hData_Tagger","Tagger",100,-50.,50.);
  hData_Tagger_TCHE     = new TH1F("hData_Tagger_TCHE","Tagger_TCHE",100,-50.,50.);
  hData_Tagger_TCHP     = new TH1F("hData_Tagger_TCHP","Tagger_TCHP",100,-50.,50.);
  hData_Tagger_JP       = new TH1F("hData_Tagger_JP","Tagger_JP",100,-50.,50.);
  hData_Tagger_SSV      = new TH1F("hData_Tagger_SSV","Tagger_SSV",100,-50.,50.);
  hData_Tagger_CSV      = new TH1F("hData_Tagger_CSV","Tagger_CSV",100,-50.,50.);
  hData_Tagger_MU       = new TH1F("hData_Tagger_MU","Tagger_MU",100,-50.,50.);
  
  hAllFlav_Flavour         = new TH1F("hAllFlav_Flavour","Flavour",22,-0.5,21.5);
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
  
  hLightFlav_Tagger        = new TH1F("hLightFlav_Tagger","Tagger",100,-50.,50.);
  hGluonFlav_Tagger        = new TH1F("hGluonFlav_Tagger","Tagger",100,-50.,50.);
  hUDSFlav_Tagger          = new TH1F("hUDSFlav_Tagger","Tagger",100,-50.,50.);
  hCFlav_Tagger            = new TH1F("hCFlav_Tagger","Tagger",100,-50.,50.);
  hBFlav_Tagger            = new TH1F("hBFlav_Tagger","Tagger",100,-50.,50.);
 
  nTuplesJets =  new TNtuple("Jets","All Jets",
  "Event_Njets:Event_Trig:Jet_Number:Jet_Flavour:Jet_Ntracks:Jet_Pt:Jet_Eta:Jet_Phi:TagN_TC1:TagP_TC1:TagN_TC2:TagP_TC2:TagN_TC3:TagP_TC3:TagN_JP:TagP_JP:Tag_JP:TagN_SSV:TagP_SSV:TagN_CSV:TagP_CSV:TagN_MU:TagP_MU:Cat_TC1:Cat_TC2:Cat_TC3:Cat_JP:Cat_SSV:Cat_MU:Muon_hits:Muon_chi2:Muon_d0:Muon_Pt:Muon_PtRel");
  }

 
MistagAnalyzer::~MistagAnalyzer()
 
{
  rootFile_->cd();

  hData_All_NJets       ->Write();
  hData_All_NTracks     ->Write();
  hData_All_JetPt       ->Write();
  hData_All_JetEta      ->Write();
  hData_NJets           ->Write();
  hData_NTracks         ->Write();
  hData_JetPt           ->Write();
  hData_JetEta          ->Write();
  hData_Tagger          ->Write();
  hData_Tagger_TCHE     ->Write();
  hData_Tagger_TCHP     ->Write();
  hData_Tagger_JP       ->Write();
  hData_Tagger_SSV      ->Write();
  hData_Tagger_CSV      ->Write();
  hData_Tagger_MU       ->Write();
  
  hAllFlav_Flavour         ->Write();
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
 
  hLightFlav_Tagger        ->Write();
  hGluonFlav_Tagger        ->Write();
  hUDSFlav_Tagger          ->Write();
  hCFlav_Tagger            ->Write();
  hBFlav_Tagger            ->Write();
 
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
  if(useTrackHistory_) classifier_.newEvent(iEvent, iSetup);
  
  const JetCorrector *acorrector = JetCorrector::getJetCorrector(jetCorrector_,iSetup);
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
  edm::Handle<reco::JetTagCollection> jetTags_CombinedSvtx;
  iEvent.getByLabel(combinedSvtxModuleName_, jetTags_CombinedSvtx);
  edm::Handle<reco::JetTagCollection> jetTags_negCombinedSvtx;
  iEvent.getByLabel(combinedSvtxNegModuleName_, jetTags_negCombinedSvtx);
  
  edm::Handle<reco::JetTagCollection> jetTags_Svtx;
  iEvent.getByLabel(svtxModuleName_, jetTags_Svtx);
  edm::Handle<reco::JetTagCollection> jetTags_negSvtx;
  iEvent.getByLabel(svtxNegModuleName_, jetTags_negSvtx);
  
  //------------------------------------------------------
  //Soft muon tagger
  //------------------------------------------------------
  
  edm::Handle<reco::JetTagCollection> jetTags_softMu;
  iEvent.getByLabel(softMuonModuleName_, jetTags_softMu);
  edm::Handle<reco::JetTagCollection> jetTags_softMuneg;
  iEvent.getByLabel(softMuonNegModuleName_, jetTags_softMuneg);
  
  edm::Handle<reco::SoftLeptonTagInfoCollection> tagInos_softmuon;
  iEvent.getByLabel(softMuonTagInfoName_, tagInos_softmuon);
  

  //------------------------------------------------------
  //Trigger info
  //------------------------------------------------------
  
  TriggerResults tr;
  Handle<TriggerResults> h_trigRes;
  iEvent.getByLabel(InputTag("TriggerResults::HLT8E29"), h_trigRes);
  tr = *h_trigRes;
  int BitTrigger = 0;

  if ( TriggerInfo_ ) { 
    vector<string> triggerList;
    Service<service::TriggerNamesService> tns;
    bool foundNames = tns->getTrigPaths(tr,triggerList);
    if (!foundNames) cout << "Could not get trigger names!\n";
    if (tr.size()!=triggerList.size()) cout << "ERROR: length of names and paths not the same: " << triggerList.size() << "," << tr.size() << endl;
// dump trigger list at first event
//   std::cout << "*********************" << std::endl;
    for (unsigned int i=0; i< tr.size(); i++) {
    if ( !tr[i].accept() == 1 ) continue;
      //std::cout << "trigger list " << triggerList[i] << std::endl;
// Trigger table 10**31
//       if( triggerList[i] == "HLT_Jet110"                     ) BitTrigger +=1 ;
//       if( triggerList[i] == "HLT_Jet180"                     ) BitTrigger +=2 ;   
//       if( triggerList[i] == "HLT_Jet250"                     ) BitTrigger +=4 ;   
//       if( triggerList[i] == "HLT_DiJetAve30"                 ) BitTrigger +=10 ;   
//       if( triggerList[i] == "HLT_DiJetAve50"                 ) BitTrigger +=20 ;   
//       if( triggerList[i] == "HLT_DiJetAve70"                 ) BitTrigger +=40 ;   
//       if( triggerList[i] == "HLT_TripleJet85"                ) BitTrigger +=100 ;   
//       if( triggerList[i] == "HLT_QuadJet30"                  ) BitTrigger +=200 ;   
//       if( triggerList[i] == "HLT_QuadJet60"                  ) BitTrigger +=400 ;   
//       if( triggerList[i] == "HLT_BTagMu_DoubleJet60_Relaxed" ) BitTrigger +=1000 ;	      
//       if( triggerList[i] == "HLT_BTagMu_TripleJet40_Relaxed" ) BitTrigger +=2000 ;         
//       if( triggerList[i] == "HLT_BTagMu_QuadJet30_Relaxed" )   BitTrigger +=4000 ;         
//       if( triggerList[i] == "HLT_Ele15_SW_L1R" ) BitTrigger +=10000 ;		
//       if( triggerList[i] == "HLT_Ele15_LW_L1R" ) BitTrigger +=20000 ;		
//       if( triggerList[i] == "HLT_Mu9" )          BitTrigger +=40000 ;		
// Trigger table 10**29
      if( triggerList[i] == "HLT_Jet15U"  		     ) BitTrigger +=1 ;
      if( triggerList[i] == "HLT_Jet30U"                     ) BitTrigger +=2 ;   
      if( triggerList[i] == "HLT_Jet50U"                     ) BitTrigger +=4 ;   
      if( triggerList[i] == "HLT_DiJetAve15U"                ) BitTrigger +=10 ;   
      if( triggerList[i] == "HLT_DiJetAve30U"                ) BitTrigger +=20 ;   
      if( triggerList[i] == "HLT_QuadJet15U"                 ) BitTrigger +=100 ;   
      if( triggerList[i] == "HLT_HT100U"                     ) BitTrigger +=200 ;   
      if( triggerList[i] == "HLT_BTagMu_Jet10U"              ) BitTrigger +=1000 ;	      
//      std::cout << "trigger list " << triggerList[i]  <<  endl;
    }
  }
  
  if (flavourMatchOptionf == "fastMC" && isData_) {
    iEvent.getByLabel(flavourSourcef, jetMC);
    for(JetFlavourMatchingCollection::const_iterator iter =
	  jetMC->begin(); iter != jetMC->end(); iter++)
      flavoursMapf.insert(std::pair <const edm::RefToBase<reco::Jet>, unsigned int>((*iter).first, (unsigned int)((*iter).second).getFlavour()));
  } else if (flavourMatchOptionf == "genParticle") {
    iEvent.getByLabel (flavourSourcef, theJetPartonMapf);
  }
  
  int cap0, cap1, cap2, cap3, cap4, cap5, cap6, cap7, cap8;
  int can0, can1, can2, can3, can4, can5, can6, can7, can8;

  //*********************************
  // Loop over the jets
  
  CaloJetCollection::const_iterator jet;
  
  int Njets = 0;
  for ( jet = recoJets.begin(); jet != recoJets.end(); ++jet ) {
    double JES     =  acorrector->correction(*jet, iEvent, iSetup);
    double jetpt =   (*jet).pt()  ;
    double jeteta=   (*jet).eta() ;
    if ( (jetpt * JES) <= minJetPt_ || std::fabs( jeteta ) >= maxJetEta_ ) continue;
    Njets++; 
//    cout << "Jet " <<  Njets << " " << jetpt * JES <<  endl;
  }
    
  int numjet = 0;
  for ( jet = recoJets.begin(); jet != recoJets.end(); ++jet ) {

    int ith_tagged = -1;      
    float Flavour  ;    
    int flavour    ; 
    if ( isData_ != 0 ) {
      //Flavour  = getMatchedParton(*jet).flavour();    
      //flavour    = int(TMath::Abs( Flavour )); 
      flavour = abs(getMatchedParton(*jet).getFlavour());   
      if ( flavour >= 1 && flavour <= 3 ) flavour = 1;  
      Flavour = flavour;
      //Jet_Flavour = flavour;
    }

    double JES     =  acorrector->correction(*jet, iEvent, iSetup);

    double jetpt =   (*jet).pt()  ;
    double jeteta=   (*jet).eta() ;
    double ptjet = jetpt* JES;
    
//     bool TagPos = false;
//     bool TagNeg = false;
//     bool Veto = false;

    float etajet = TMath::Abs( (*jet).eta());
    float phijet = (*jet).phi();
    if (phijet < 0.) phijet += 2*TMath::Pi();
    
    //*****************************************************************
    //Taggers 
    //*****************************************************************
   
    ith_tagged = this->TaggedJet(*jet,jetTags_JP);
    
    int ntagtracks = (*tagInfo)[ith_tagged].probabilities(0).size();
    float Proba  = (*jetTags_JP)[ith_tagged].second;
    
    ith_tagged = this->TaggedJet(*jet,jetTags_PosJP);
    float ProbaP = (*jetTags_PosJP)[ith_tagged].second;
    
    ith_tagged = this->TaggedJet(*jet,jetTags_NegJP);
    float ProbaN = (*jetTags_NegJP)[ith_tagged].second;
    
    ith_tagged              = this->TaggedJet(*jet,jetTags_CombinedSvtx);
    float   CombinedSvtx    = (*jetTags_CombinedSvtx)[ith_tagged].second;
    ith_tagged              = this->TaggedJet(*jet,jetTags_negCombinedSvtx);
    float   CombinedSvtxN   = (*jetTags_negCombinedSvtx)[ith_tagged].second;
    
    //reco::TaggingVariableList variables =
    //                   svComputer((*tagInfo)[ith_tagged], (*tagInfoSVx)[ith_tagged]);
    //int vertexCat = variables[reco::btau::vertexCategory];
    //int vertexCat = -1;

    ith_tagged            = this->TaggedJet(*jet,jetTags_Svtx);
    float   Svtx          = (*jetTags_Svtx)[ith_tagged].second;
    ith_tagged            = this->TaggedJet(*jet,jetTags_negSvtx);
    float   SvtxN         = (*jetTags_negSvtx)[ith_tagged].second;
    
    float   SoftM   = 0;
    float   SoftMN  = 0;
    
    ith_tagged            = this->TaggedJet(*jet,jetTags_softMu);
    if( (*jetTags_softMu)[ith_tagged].second     > -100000 )   SoftM     = (*jetTags_softMu)[ith_tagged].second;
    ith_tagged            = this->TaggedJet(*jet,jetTags_softMuneg);
    if( (*jetTags_softMuneg)[ith_tagged].second  > -100000 )  SoftMN     = ((*jetTags_softMuneg)[ith_tagged].second);
    if(SoftMN > 0) SoftMN = -1*SoftMN;
    
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
    
    //*****************************************************************
    //get track histories for 1st, 2nd and 3rd track (TC)
    //*****************************************************************
    int Category1 = 0,  Category2 = 0, Category3 = 0;
    
    // Track history
    if (useTrackHistory_ && indexes.size()!=0 && isData_!=0) {
//$$ ????
      flags1P[TrackCategories::Conversion] ;
//$$ ????
      if ( flags1P[TrackCategories::BWeakDecay] )  Category1 += int(pow(10., -1 + 1)); 
      if ( flags1P[TrackCategories::CWeakDecay] )  Category1 += int(pow(10., -1 + 2)); 
      if ( flags1P[TrackCategories::TauDecay] )    Category1 += int(pow(10., -1 + 3)); 
      if ( flags1P[TrackCategories::Conversion] )  Category1 += int(pow(10., -1 + 4)); 
      if ( flags1P[TrackCategories::KsDecay] )     Category1 += int(pow(10., -1 + 5)); 
      if ( flags1P[TrackCategories::LambdaDecay] ) Category1 += int(pow(10., -1 + 6)); 
      if ( flags1P[TrackCategories::Interaction] ) Category1 += int(pow(10., -1 + 7)); 
      if ( flags1P[TrackCategories::Fake] )	   Category1 += int(pow(10., -1 + 8)); 
      if ( flags1P[TrackCategories::Bad] )	   Category1 += int(pow(10., -1 + 9)); 
      if ( idSize > 1 ) {
        if ( flags2P[TrackCategories::BWeakDecay] )  Category2 += int(pow(10., -1 + 1)); 
        if ( flags2P[TrackCategories::CWeakDecay] )  Category2 += int(pow(10., -1 + 2)); 
        if ( flags2P[TrackCategories::TauDecay] )    Category2 += int(pow(10., -1 + 3)); 
        if ( flags2P[TrackCategories::Conversion] )  Category2 += int(pow(10., -1 + 4)); 
        if ( flags2P[TrackCategories::KsDecay] )     Category2 += int(pow(10., -1 + 5)); 
        if ( flags2P[TrackCategories::LambdaDecay] ) Category2 += int(pow(10., -1 + 6)); 
        if ( flags2P[TrackCategories::Interaction] ) Category2 += int(pow(10., -1 + 7)); 
        if ( flags2P[TrackCategories::Fake] )	     Category2 += int(pow(10., -1 + 8)); 
        if ( flags2P[TrackCategories::Bad] )	     Category2 += int(pow(10., -1 + 9)); 
      }
      if ( idSize > 2 ) {
        if ( flags3P[TrackCategories::BWeakDecay] )  Category3 += int(pow(10., -1 + 1)); 
        if ( flags3P[TrackCategories::CWeakDecay] )  Category3 += int(pow(10., -1 + 2)); 
        if ( flags3P[TrackCategories::TauDecay] )    Category3 += int(pow(10., -1 + 3)); 
        if ( flags3P[TrackCategories::Conversion] )  Category3 += int(pow(10., -1 + 4)); 
        if ( flags3P[TrackCategories::KsDecay] )     Category3 += int(pow(10., -1 + 5)); 
        if ( flags3P[TrackCategories::LambdaDecay] ) Category3 += int(pow(10., -1 + 6)); 
        if ( flags3P[TrackCategories::Interaction] ) Category3 += int(pow(10., -1 + 7)); 
        if ( flags3P[TrackCategories::Fake] )	     Category3 += int(pow(10., -1 + 8)); 
        if ( flags3P[TrackCategories::Bad] )	     Category3 += int(pow(10., -1 + 9)); 
      }
      if ( flags1N[TrackCategories::BWeakDecay] )  Category1 += 2*int(pow(10., -1 + 1)); 
      if ( flags1N[TrackCategories::CWeakDecay] )  Category1 += 2*int(pow(10., -1 + 2)); 
      if ( flags1N[TrackCategories::TauDecay] )    Category1 += 2*int(pow(10., -1 + 3)); 
      if ( flags1N[TrackCategories::Conversion] )  Category1 += 2*int(pow(10., -1 + 4)); 
      if ( flags1N[TrackCategories::KsDecay] )     Category1 += 2*int(pow(10., -1 + 5)); 
      if ( flags1N[TrackCategories::LambdaDecay] ) Category1 += 2*int(pow(10., -1 + 6)); 
      if ( flags1N[TrackCategories::Interaction] ) Category1 += 2*int(pow(10., -1 + 7)); 
      if ( flags1N[TrackCategories::Fake] )	   Category1 += 2*int(pow(10., -1 + 8)); 
      if ( flags1N[TrackCategories::Bad] )	   Category1 += 2*int(pow(10., -1 + 9)); 
      if ( idSize > 1 ) {
        if ( flags2N[TrackCategories::BWeakDecay] )  Category2 += 2*int(pow(10., -1 + 1)); 
        if ( flags2N[TrackCategories::CWeakDecay] )  Category2 += 2*int(pow(10., -1 + 2)); 
        if ( flags2N[TrackCategories::TauDecay] )    Category2 += 2*int(pow(10., -1 + 3)); 
        if ( flags2N[TrackCategories::Conversion] )  Category2 += 2*int(pow(10., -1 + 4)); 
        if ( flags2N[TrackCategories::KsDecay] )     Category2 += 2*int(pow(10., -1 + 5)); 
        if ( flags2N[TrackCategories::LambdaDecay] ) Category2 += 2*int(pow(10., -1 + 6)); 
        if ( flags2N[TrackCategories::Interaction] ) Category2 += 2*int(pow(10., -1 + 7)); 
        if ( flags2N[TrackCategories::Fake] )	     Category2 += 2*int(pow(10., -1 + 8)); 
        if ( flags2N[TrackCategories::Bad] )	     Category2 += 2*int(pow(10., -1 + 9)); 
      }
      if ( idSize > 2 ) {
        if ( flags3N[TrackCategories::BWeakDecay] )  Category3 += 2*int(pow(10., -1 + 1)); 
        if ( flags3N[TrackCategories::CWeakDecay] )  Category3 += 2*int(pow(10., -1 + 2)); 
        if ( flags3N[TrackCategories::TauDecay] )    Category3 += 2*int(pow(10., -1 + 3)); 
        if ( flags3N[TrackCategories::Conversion] )  Category3 += 2*int(pow(10., -1 + 4)); 
        if ( flags3N[TrackCategories::KsDecay] )     Category3 += 2*int(pow(10., -1 + 5)); 
        if ( flags3N[TrackCategories::LambdaDecay] ) Category3 += 2*int(pow(10., -1 + 6)); 
        if ( flags3N[TrackCategories::Interaction] ) Category3 += 2*int(pow(10., -1 + 7)); 
        if ( flags3N[TrackCategories::Fake] )	     Category3 += 2*int(pow(10., -1 + 8)); 
        if ( flags3N[TrackCategories::Bad] )	     Category3 += 2*int(pow(10., -1 + 9)); 
      }
    }
      
    //*****************************************************************
    //get track histories of tracks in jets (for Jet Proba) 
    //*****************************************************************
    int CategoryJet = 0;
    
    if(useTrackHistory_ && isData_!=0){
      ith_tagged = this->TaggedJet(*jet,jetTags_PosJP);
      TrackRefVector jetProbTracks( (*tagInfo)[ith_tagged].selectedTracks() );  
      
      cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0; 
      can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0; 

      for(unsigned int i=0; i<jetProbTracks.size(); i++){
	reco::TrackIPTagInfo::TrackIPData ip = ((*tagInfo)[ith_tagged].impactParameterData())[i];

	if(ip.ip3d.significance()> 0) {
	  TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
	  if ( theFlag[TrackCategories::BWeakDecay] )	cap0 = 1; 
	  if ( theFlag[TrackCategories::CWeakDecay] )	cap1 = 1; 
	  if ( theFlag[TrackCategories::TauDecay] )	cap2 = 1;  
	  if ( theFlag[TrackCategories::Conversion] )   cap3 = 1; 
	  if ( theFlag[TrackCategories::KsDecay] )	cap4 = 1; 
	  if ( theFlag[TrackCategories::LambdaDecay] )  cap5 = 1; 
	  if ( theFlag[TrackCategories::Interaction] )  cap6 = 1; 
	  if ( theFlag[TrackCategories::Fake] ) 	cap7 = 1; 
	  if ( theFlag[TrackCategories::Bad] )  	cap8 = 1; 
	}
	else {
	  TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
	  if ( theFlag[TrackCategories::BWeakDecay] )	can0 = 2; 
	  if ( theFlag[TrackCategories::CWeakDecay] )	can1 = 2; 
	  if ( theFlag[TrackCategories::TauDecay] )	can2 = 2;  
	  if ( theFlag[TrackCategories::Conversion] )   can3 = 2; 
	  if ( theFlag[TrackCategories::KsDecay] )	can4 = 2; 
	  if ( theFlag[TrackCategories::LambdaDecay] )  can5 = 2; 
	  if ( theFlag[TrackCategories::Interaction] )  can6 = 2; 
	  if ( theFlag[TrackCategories::Fake] ) 	can7 = 2; 
	  if ( theFlag[TrackCategories::Bad] )  	can8 = 2; 
	}
      }
      CategoryJet =  cap0+can0           + (cap1+can1)*10    + (cap2+can2)*100 
                  + (cap3+can3)*1000     + (cap4+can4)*10000 
                  + (cap5+can5)*100000   + (cap6+can6)*1000000 
		  + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }
    
    //*****************************************************************
    //get track histories  associated to sec. vertex (for simple SV)
    //*****************************************************************
    int CategorySVx = 0;
    
    if(useTrackHistory_ && isData_!=0) {
      ith_tagged =    this->TaggedJet(*jet,jetTags_Svtx);
      TrackRefVector  svxPostracks( (*tagInfoSVx)[ith_tagged].vertexTracks(0) );
      
      cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0; 
      can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0; 

      for(unsigned int i=0; i<svxPostracks.size(); i++){
	TrackCategories::Flags theFlag = classifier_.evaluate( svxPostracks[i] ).flags();
	if ( theFlag[TrackCategories::BWeakDecay] )   cap0 = 1; 
	if ( theFlag[TrackCategories::CWeakDecay] )   cap1 = 1; 
	if ( theFlag[TrackCategories::TauDecay] )     cap2 = 1;  
	if ( theFlag[TrackCategories::Conversion] )   cap3 = 1; 
	if ( theFlag[TrackCategories::KsDecay] )      cap4 = 1; 
	if ( theFlag[TrackCategories::LambdaDecay] )  cap5 = 1; 
	if ( theFlag[TrackCategories::Interaction] )  cap6 = 1; 
	if ( theFlag[TrackCategories::Fake] )	      cap7 = 1; 
	if ( theFlag[TrackCategories::Bad] )	      cap8 = 1; 
      }
      ith_tagged  =    this->TaggedJet(*jet,jetTags_negSvtx);
      TrackRefVector   svxNegtracks( (*tagInfoNegSVx)[ith_tagged].vertexTracks(0) );
      for(unsigned int i=0; i<svxNegtracks.size(); i++){
	TrackCategories::Flags theFlag = classifier_.evaluate( svxNegtracks[i] ).flags();
	if ( theFlag[TrackCategories::BWeakDecay] )   can0 = 2; 
	if ( theFlag[TrackCategories::CWeakDecay] )   can1 = 2; 
	if ( theFlag[TrackCategories::TauDecay] )     can2 = 2;  
	if ( theFlag[TrackCategories::Conversion] )   can3 = 2; 
	if ( theFlag[TrackCategories::KsDecay] )      can4 = 2; 
	if ( theFlag[TrackCategories::LambdaDecay] )  can5 = 2; 
	if ( theFlag[TrackCategories::Interaction] )  can6 = 2; 
	if ( theFlag[TrackCategories::Fake] )	      can7 = 2; 
	if ( theFlag[TrackCategories::Bad] )	      can8 = 2; 
      }
      CategorySVx =  cap0+can0           + (cap1+can1)*10    + (cap2+can2)*100 
                  + (cap3+can3)*1000     + (cap4+can4)*10000 
                  + (cap5+can5)*100000   + (cap6+can6)*1000000 
		  + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }
    
    float mu_NHits_tracker  = -10000;
    float mu_Chi2           = -10000;
    float mu_d0             = -10000;
    float mu_pT             = -10000;
    float ptRel             = -10000;
    //*****************************************************************
    //get track histories of the muon (SoftMuon tagger)
    //*****************************************************************
    int CategoryMuon = 0;

    ith_tagged = this->TaggedJet(*jet,jetTags_softMu);
    if( (*tagInos_softmuon)[ith_tagged].leptons()!=0 ) {
      ptRel =  calculPtRel( (*(*tagInos_softmuon)[ith_tagged].lepton(0)), *jet, JES);
      mu_NHits_tracker  = (*tagInos_softmuon)[ith_tagged].lepton(0)->hitPattern().numberOfValidHits();
      mu_Chi2           = (*tagInos_softmuon)[ith_tagged].lepton(0)->normalizedChi2()                ;
      mu_d0             = (*tagInos_softmuon)[ith_tagged].lepton(0)->d0()                            ;
      mu_pT             = (*tagInos_softmuon)[ith_tagged].lepton(0)->pt()                            ;
    }
    if(  SoftM > 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ){
      if ( useTrackHistory_ && isData_!=0 ) {     
        TrackCategories::Flags theFlagP = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
        if ( theFlagP[TrackCategories::BWeakDecay] )   CategoryMuon += int(pow(10., -1 + 1)); 
        if ( theFlagP[TrackCategories::CWeakDecay] )   CategoryMuon += int(pow(10., -1 + 2)); 
        if ( theFlagP[TrackCategories::TauDecay] )     CategoryMuon += int(pow(10., -1 + 3));  
        if ( theFlagP[TrackCategories::Conversion] )   CategoryMuon += int(pow(10., -1 + 4)); 
        if ( theFlagP[TrackCategories::KsDecay] )      CategoryMuon += int(pow(10., -1 + 5)); 
        if ( theFlagP[TrackCategories::LambdaDecay] )  CategoryMuon += int(pow(10., -1 + 6)); 
        if ( theFlagP[TrackCategories::Interaction] )  CategoryMuon += int(pow(10., -1 + 7)); 
        if ( theFlagP[TrackCategories::Fake] )	       CategoryMuon += int(pow(10., -1 + 8)); 
        if ( theFlagP[TrackCategories::Bad] )	       CategoryMuon += int(pow(10., -1 + 9)); 
      }
    }

    ith_tagged = this->TaggedJet(*jet,jetTags_softMuneg);
     //std::cout << "SoftMN " << SoftMN << std::endl;
    if ( SoftMN < 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ) {
      ptRel =  calculPtRel( (*(*tagInos_softmuon)[ith_tagged].lepton(0)), *jet, JES);
      mu_NHits_tracker  = (*tagInos_softmuon)[ith_tagged].lepton(0)->hitPattern().numberOfValidHits();
      mu_Chi2           = (*tagInos_softmuon)[ith_tagged].lepton(0)->normalizedChi2()                ;
      mu_d0             = (*tagInos_softmuon)[ith_tagged].lepton(0)->d0()                            ;
      mu_pT             = (*tagInos_softmuon)[ith_tagged].lepton(0)->pt()                            ;
      if ( useTrackHistory_ && isData_!=0 ) {     
        TrackCategories::Flags theFlagN = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
        if ( theFlagN[TrackCategories::BWeakDecay] )   CategoryMuon += 2*int(pow(10., -1 + 1)); 
        if ( theFlagN[TrackCategories::CWeakDecay] )   CategoryMuon += 2*int(pow(10., -1 + 2)); 
        if ( theFlagN[TrackCategories::TauDecay] )     CategoryMuon += 2*int(pow(10., -1 + 3));  
        if ( theFlagN[TrackCategories::Conversion] )   CategoryMuon += 2*int(pow(10., -1 + 4)); 
        if ( theFlagN[TrackCategories::KsDecay] )      CategoryMuon += 2*int(pow(10., -1 + 5)); 
        if ( theFlagN[TrackCategories::LambdaDecay] )  CategoryMuon += 2*int(pow(10., -1 + 6)); 
        if ( theFlagN[TrackCategories::Interaction] )  CategoryMuon += 2*int(pow(10., -1 + 7)); 
        if ( theFlagN[TrackCategories::Fake] )	       CategoryMuon += 2*int(pow(10., -1 + 8)); 
        if ( theFlagN[TrackCategories::Bad] )	       CategoryMuon += 2*int(pow(10., -1 + 9)); 
      }
    }
        
    //*********************************
    // Jet selection
    
    if ( (jetpt * JES) <= minJetPt_ || std::fabs( jeteta ) >= maxJetEta_ ) continue;
    numjet++;
    
    float Ijet =       numjet;
    float Ntagtracks = ntagtracks;
    
    
    //*****************************************************************
    // the ntuple
    //*****************************************************************
     //Event_Njets:Event_Trig:
     //Jet_Number:Jet_Flavour:Jet_Ntracks:Jet_Pt:Jet_Eta:Jet_Phi:
     //TagN_TC1:TagP_TC1:TagN_TC2:TagP_TC2:TagN_TC3:TagP_TC3:TagN_JP:TagP_JP :Tag_JP:
     //TagN_SSV:TagP_SSV:TagN_CSV:TagP_CSV:TagN_MU:TagP_MU:
     //Cat_TC1:Cat_TC2:Cat_TC3:Cat_JP:Cat_SSV:Cat_MU:
     //Muon_hits:Muon_chi2:Muon_d0:Muon_Pt:Muon_PtRel");
 
    float jet_input[34] =
      {Njets, BitTrigger, 
       Ijet, Flavour, Ntagtracks, ptjet, (*jet).eta(), phijet, 
       Ip1N, Ip1P, Ip2N, Ip2P, Ip3N, Ip3P, ProbaN, ProbaP, Proba,
       SvtxN, Svtx, CombinedSvtxN, CombinedSvtx, SoftMN, SoftM, 
       Category1, Category2, Category3, CategoryJet, CategorySVx, CategoryMuon,
       mu_NHits_tracker,mu_Chi2,mu_d0,mu_pT,ptRel
      };
    
    
    nTuplesJets->Fill(jet_input);
    
    //*****************************************************************
    //define positive and negative tags
    //*****************************************************************
    float varpos = -1000.;
    float varneg = -1000.;
    
    if ( selTagger_ == 0 ) {       // jet proba
      if ( ProbaP > 0 ) varpos = 20.*ProbaP;
      if ( ProbaN > 0 ) varneg = 20.*ProbaN;
//       if ( ProbaP > tagCut_ ) TagPos = true;
//       if ( ProbaN > tagCut_ ) TagNeg = true;
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
      if  (Svtx > 1)   varpos =  10.*Svtx - 10.;
      if  (SvtxN < -1) varneg = -10.*SvtxN - 10.;
//       if ( Svtx  > tagCut_ ) TagPos = true;
//       if (-SvtxN > tagCut_ ) TagNeg = true; 
    }
    else if ( selTagger_ == 5 ) {    // SV combined
      if  (CombinedSvtx > 0)  varpos = 50.*CombinedSvtx;
      if  (CombinedSvtxN > 0) varneg = 50.*CombinedSvtxN;
//       if ( CombinedSvtx  > tagCut_ ) TagPos = true;
//       if (-CombinedSvtxN > tagCut_ ) TagNeg = true;   
    }
    else if ( selTagger_ == 6 ) {    // soft muon ptrel
      if  ( SoftM  > 0) varpos =  5*SoftM;
      if  ( SoftMN < 0) varneg = -5*SoftMN;
//       if  (  SoftM > tagCut_ ) TagPos = true;
//       if  (  SoftMN> tagCut_ ) TagNeg = true;   
    }
    
//     if ( selTagger_ >= 1 && selTagger_ <= 3) {
//       if ( varpos > tagCut_ ) TagPos = true;
//       if ( varneg > tagCut_ ) TagNeg = true;
//     }
    
//     // veto on positive tag
//     if ( Ip1P < vetoPos_ ) Veto = true;
    
    //***************************
    
    int cat1P = 0, cat2P = 0, cat3P = 0, catP = 0;
    int cat1N = 0, cat2N = 0, cat3N = 0, catN = 0;

    // Track history
    if (useTrackHistory_ && indexes.size()!=0 && isData_!=0) {
      flags1P[TrackCategories::Conversion] ;
           if ( flags1P[TrackCategories::BWeakDecay] )  cat1P = 1; 
      else if ( flags1P[TrackCategories::CWeakDecay] )  cat1P = 2; 
      else if ( flags1P[TrackCategories::TauDecay] )    cat1P = 3; 
      else if ( flags1P[TrackCategories::Conversion] )  cat1P = 4; 
      else if ( flags1P[TrackCategories::KsDecay] )     cat1P = 5; 
      else if ( flags1P[TrackCategories::LambdaDecay] ) cat1P = 6; 
      else if ( flags1P[TrackCategories::Interaction] ) cat1P = 7; 
      else if ( flags1P[TrackCategories::Fake] )        cat1P = 8; 
      else if ( flags1P[TrackCategories::Bad] )         cat1P = 9; 
      if(idSize > 1){
	     if ( flags2P[TrackCategories::BWeakDecay] )  cat2P = 1;
	else if ( flags2P[TrackCategories::CWeakDecay] )  cat2P = 2;
	else if ( flags2P[TrackCategories::TauDecay] )    cat2P = 3;
	else if ( flags2P[TrackCategories::Conversion] )  cat2P = 4;
	else if ( flags2P[TrackCategories::KsDecay] )     cat2P = 5;
	else if ( flags2P[TrackCategories::LambdaDecay] ) cat2P = 6;
	else if ( flags2P[TrackCategories::Interaction] ) cat2P = 7;
	else if ( flags2P[TrackCategories::Fake] )        cat2P = 8;
	else if ( flags2P[TrackCategories::Bad] )         cat2P = 9;
      }
      if(idSize > 2){
	     if ( flags3P[TrackCategories::BWeakDecay] )  cat3P = 1;
	else if ( flags3P[TrackCategories::CWeakDecay] )  cat3P = 2;
	else if ( flags3P[TrackCategories::TauDecay] )    cat3P = 3;
	else if ( flags3P[TrackCategories::Conversion] )  cat3P = 4;
	else if ( flags3P[TrackCategories::KsDecay] )     cat3P = 5;
	else if ( flags3P[TrackCategories::LambdaDecay] ) cat3P = 6;
	else if ( flags3P[TrackCategories::Interaction] ) cat3P = 7;
	else if ( flags3P[TrackCategories::Fake] )        cat3P = 8;
	else if ( flags3P[TrackCategories::Bad] )         cat3P = 9;
      }
           if ( flags1N[TrackCategories::BWeakDecay] )  cat1N = 1;
      else if ( flags1N[TrackCategories::CWeakDecay] )  cat1N = 2;
      else if ( flags1N[TrackCategories::TauDecay] )    cat1N = 3;
      else if ( flags1N[TrackCategories::Conversion] )  cat1N = 4;
      else if ( flags1N[TrackCategories::KsDecay] )     cat1N = 5;
      else if ( flags1N[TrackCategories::LambdaDecay] ) cat1N = 6;
      else if ( flags1N[TrackCategories::Interaction] ) cat1N = 7;
      else if ( flags1N[TrackCategories::Fake] )        cat1N = 8;
      else if ( flags1N[TrackCategories::Bad] )         cat1N = 9;
      if(idSize > 1){
	     if ( flags2N[TrackCategories::BWeakDecay] )  cat2N = 1;
	else if ( flags2N[TrackCategories::CWeakDecay] )  cat2N = 2;
	else if ( flags2N[TrackCategories::TauDecay] )    cat2N = 3;
	else if ( flags2N[TrackCategories::Conversion] )  cat2N = 4;
	else if ( flags2N[TrackCategories::KsDecay] )     cat2N = 5;
	else if ( flags2N[TrackCategories::LambdaDecay] ) cat2N = 6;
	else if ( flags2N[TrackCategories::Interaction] ) cat2N = 7;
	else if ( flags2N[TrackCategories::Fake] )        cat2N = 8;
	else if ( flags2N[TrackCategories::Bad] )         cat2N = 9;
      }
      if(idSize > 2){
	     if ( flags3N[TrackCategories::BWeakDecay] )  cat3N = 1;
	else if ( flags3N[TrackCategories::CWeakDecay] )  cat3N = 2;
	else if ( flags3N[TrackCategories::TauDecay] )    cat3N = 3;
	else if ( flags3N[TrackCategories::Conversion] )  cat3N = 4;
	else if ( flags3N[TrackCategories::KsDecay] )     cat3N = 5;
	else if ( flags3N[TrackCategories::LambdaDecay] ) cat3N = 6;
	else if ( flags3N[TrackCategories::Interaction] ) cat3N = 7;
	else if ( flags3N[TrackCategories::Fake] )        cat3N = 8;
	else if ( flags3N[TrackCategories::Bad] )         cat3N = 9;
      }
      
      if ( selTagger_ == 0 || selTagger_ == 3 ) {
	catP = cat3P;
	catN = cat3N;
      }
      else if ( selTagger_ == 2 ) {
	catP = cat2P;
	catN = cat2N;
      }
      else if ( selTagger_ == 1 ) {
	catP = cat1P;
	catN = cat1N;
      }
    }
    
    //*****************************************************************
    //fill the histograms
    //*****************************************************************

    hData_All_NTracks->Fill( ntagtracks ) ;
    hData_All_JetPt->Fill( ptjet );
    hData_All_JetEta->Fill( etajet );
    
    //*********************************
    // Taggability
    //$$
    if ( ntagtracks < ntrackMin_ ) continue;
    //$$

    hData_NTracks->Fill( ntagtracks ) ;
    hData_JetPt->Fill( ptjet );
    hData_JetEta->Fill( etajet );
    
    if ( isData_ != 0 ) {
      hAllFlav_Flavour->Fill( flavour );
    }

    //*********************************
    // Tagging
    
    if ( varneg > 0 ) hData_Tagger->Fill(-varneg );
    if ( varpos > 0 ) hData_Tagger->Fill( varpos );

    if ( Ip2P > 0 )   hData_Tagger_TCHE->Fill( Ip2P );
    if ( Ip2N > 0 )   hData_Tagger_TCHE->Fill(-Ip2N );
    if ( Ip3P > 0 )   hData_Tagger_TCHP->Fill( Ip3P );
    if ( Ip3N > 0 )   hData_Tagger_TCHP->Fill(-Ip3N );
    if ( ProbaP > 0 ) hData_Tagger_JP->Fill( 20.*ProbaP );
    if ( ProbaN > 0 ) hData_Tagger_JP->Fill(-20.*ProbaN );
    if ( Svtx  >  1 ) hData_Tagger_SSV->Fill( 10.*Svtx - 10 );
    if ( SvtxN < -1 ) hData_Tagger_SSV->Fill( 10*SvtxN + 10 );
    if  (CombinedSvtx > 0 )  hData_Tagger_CSV->Fill( 50.*CombinedSvtx );
    if ( CombinedSvtxN > 0 ) hData_Tagger_CSV->Fill(-50.*CombinedSvtxN );
    if  ( SoftM  > 0)  hData_Tagger_MU->Fill( 5*SoftM );
    if  ( SoftMN < 0 ) hData_Tagger_MU->Fill( 5*SoftMN );

    if ( isData_ != 0 ) {
      if ( varneg > 0 ) hAllFlav_Tagger->Fill(-varneg );
      if ( varpos > 0 ) hAllFlav_Tagger->Fill( varpos );
      
      if ( varneg > 0 ) {
	if      ( catN == 1 ) hAllFlav_Tagger_Bwd->Fill(-varneg );
	else if ( catN == 2 ) hAllFlav_Tagger_Cwd->Fill(-varneg );
	else if ( catN == 3 ) hAllFlav_Tagger_Tau->Fill(-varneg );
	else if ( catN == 4 ) hAllFlav_Tagger_Gam->Fill(-varneg );
	else if ( catN == 5 ) hAllFlav_Tagger_K0s->Fill(-varneg );
	else if ( catN == 6 ) hAllFlav_Tagger_Lam->Fill(-varneg );
	else if ( catN == 7 ) hAllFlav_Tagger_Int->Fill(-varneg );
	else if ( catN == 8 ) hAllFlav_Tagger_Fak->Fill(-varneg );
	else if ( catN == 9 ) hAllFlav_Tagger_Bad->Fill(-varneg );
	else                  hAllFlav_Tagger_Oth->Fill(-varneg );
      }
      if ( varpos > 0 ) {
	if      ( catP == 1 ) hAllFlav_Tagger_Bwd->Fill( varpos );
	else if ( catP == 2 ) hAllFlav_Tagger_Cwd->Fill( varpos );
	else if ( catP == 3 ) hAllFlav_Tagger_Tau->Fill( varpos );
	else if ( catP == 4 ) hAllFlav_Tagger_Gam->Fill( varpos );
	else if ( catP == 5 ) hAllFlav_Tagger_K0s->Fill( varpos );
	else if ( catP == 6 ) hAllFlav_Tagger_Lam->Fill( varpos );
	else if ( catP == 7 ) hAllFlav_Tagger_Int->Fill( varpos );
	else if ( catP == 8 ) hAllFlav_Tagger_Fak->Fill( varpos );
	else if ( catP == 9 ) hAllFlav_Tagger_Bad->Fill( varpos );
	else                  hAllFlav_Tagger_Oth->Fill( varpos );
      }
      
      if ( flavour == 1 || flavour == 21 ) { // light q+gluon
	if ( varneg > 0 ) hLightFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hLightFlav_Tagger->Fill( varpos );
      } // light q+gluon
      
      if ( flavour == 21 ) { // gluon jets
	if ( varneg > 0 ) hGluonFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hGluonFlav_Tagger->Fill( varpos );
      } // gluon jets
      
      else if ( flavour == 1 ) { // uds jets
	if ( varneg > 0 ) hUDSFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hUDSFlav_Tagger->Fill( varpos );
      } // uds jets
      
      else if ( flavour == 4 ) { // c jets
	if ( varneg > 0 ) hCFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hCFlav_Tagger->Fill( varpos );
      } // c jets
      
      else if ( flavour == 5 ) { // b jets
	if ( varneg > 0 ) hBFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hBFlav_Tagger->Fill( varpos );
      } // b jets
    }
    
  } // end loop on jet
  
  hData_All_NJets->Fill( Njets );
  hData_NJets->Fill( numjet );
}


float 
MistagAnalyzer::calculPtRel(reco::Track theMuon, reco::Jet theJet, double JES )
{
  double pmu = TMath::Sqrt( theMuon.px()*theMuon.px() + theMuon.py()*theMuon.py()  + theMuon.pz()*theMuon.pz() );

  double jetpx = theJet.px() * JES +theMuon.px();
  double jetpy = theJet.py() * JES +theMuon.py();
  double jetpz = theJet.pz() * JES +theMuon.pz();
  double jetp = TMath::Sqrt(jetpx*jetpx + jetpy*jetpy + jetpz*jetpz);

  double ptrel  = ( jetpx * theMuon.px()  + jetpy * theMuon.py() + jetpz * theMuon.pz() ) / jetp;
  ptrel = TMath::Sqrt( pmu * pmu  - ptrel * ptrel );

  return ptrel;
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
    
    //jetFlavour.underlyingParton4Vec(jet.p4());
    
  } else if (flavourMatchOptionf == "hepMC") {
    
    //jetFlavour = jetFlavourIdentifier_.identifyBasedOnPartons(jet);
    
  } else if (flavourMatchOptionf == "genParticle") {
    for ( JetFlavourMatchingCollection::const_iterator j  = theJetPartonMapf->begin();
	  j != theJetPartonMapf->end();
	  j ++ ) {
      RefToBase<Jet> aJet  = (*j).first;       //      const JetFlavour aFlav = (*j).second;
      if ( fabs(aJet->phi() - jet.phi()) < 1.e-5 && fabs(aJet->eta() - jet.eta())< 1.e-5 ){
	// matched
	jetFlavour = reco::JetFlavour (aJet->p4(), math::XYZPoint(0,0,0), (*j).second.getFlavour());
		
	//jetFlavour.flavour((*j).second.getFlavour());
      }
    }
    
    return jetFlavour;
  }
  
  return jetFlavour;
}

int MistagAnalyzer::TaggedJet(reco::CaloJet calojet, edm::Handle<reco::JetTagCollection > jetTags ) {
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
    std::string moduleLabel = (jetTags).provenance()->moduleLabel();
    std::string processname = (jetTags).provenance()->processName();
  }
  return result;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MistagAnalyzer);


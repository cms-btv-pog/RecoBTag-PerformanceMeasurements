
#include "RecoBTag/PerformanceMeasurements/interface/MistagAnalyzer.h"

MistagAnalyzer::MistagAnalyzer(const edm::ParameterSet& iConfig): classifier_(iConfig)  
{
    
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
  
  produceJetProbaTree_ = iConfig.getParameter<bool>("produceJetProbaTree");
  
  
  
  
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

  //fs->make<TH1F>(

  hData_All_NJets       = fs->make<TH1F>("hData_All_NJets","nb. of jets",21,-0.5,20.5);
  hData_All_NTracks     = fs->make<TH1F>("hData_All_NTracks","nb. of tracks",20,0.5,20.5);
  hData_All_JetPt       = fs->make<TH1F>("hData_All_JetPt","pt(jet)",21,20.,230.);
  hData_All_JetEta      = fs->make<TH1F>("hData_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_NJets           = fs->make<TH1F>("hData_NJets","nb. of jets",21,-0.5,20.5);
  hData_NTracks         = fs->make<TH1F>("hData_NTracks","nb. of tracks",20,0.5,20.5);
  hData_JetPt           = fs->make<TH1F>("hData_JetPt","pt(jet)",21,20.,230.);
  hData_JetEta          = fs->make<TH1F>("hData_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_Tagger          = fs->make<TH1F>("hData_Tagger","Tagger",100,-50.,50.);
  hData_Tagger_TCHE     = fs->make<TH1F>("hData_Tagger_TCHE","Tagger_TCHE",100,-50.,50.);
  hData_Tagger_TCHP     = fs->make<TH1F>("hData_Tagger_TCHP","Tagger_TCHP",100,-50.,50.);
  hData_Tagger_JP       = fs->make<TH1F>("hData_Tagger_JP","Tagger_JP",100,-50.,50.);
  hData_Tagger_SSV      = fs->make<TH1F>("hData_Tagger_SSV","Tagger_SSV",100,-50.,50.);
  hData_Tagger_CSV      = fs->make<TH1F>("hData_Tagger_CSV","Tagger_CSV",100,-50.,50.);
  hData_Tagger_MU       = fs->make<TH1F>("hData_Tagger_MU","Tagger_MU",100,-50.,50.);
  
  hAllFlav_Flavour         = fs->make<TH1F>("hAllFlav_Flavour","Flavour",22,-0.5,21.5);
  hAllFlav_Tagger          = fs->make<TH1F>("hAllFlav_Tagger","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Gam      = fs->make<TH1F>("hAllFlav_Tagger_Gam","Tagger",100,-50.,50.);
  hAllFlav_Tagger_K0s      = fs->make<TH1F>("hAllFlav_Tagger_K0s","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Lam      = fs->make<TH1F>("hAllFlav_Tagger_Lam","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Bwd      = fs->make<TH1F>("hAllFlav_Tagger_Bwd","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Cwd      = fs->make<TH1F>("hAllFlav_Tagger_Cwd","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Tau      = fs->make<TH1F>("hAllFlav_Tagger_Tau","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Int      = fs->make<TH1F>("hAllFlav_Tagger_Int","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Fak      = fs->make<TH1F>("hAllFlav_Tagger_Fak","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Bad      = fs->make<TH1F>("hAllFlav_Tagger_Bad","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Oth      = fs->make<TH1F>("hAllFlav_Tagger_Oth","Tagger",100,-50.,50.);
  
  hLightFlav_Tagger        = fs->make<TH1F>("hLightFlav_Tagger","Tagger",100,-50.,50.);
  hGluonFlav_Tagger        = fs->make<TH1F>("hGluonFlav_Tagger","Tagger",100,-50.,50.);
  hUDSFlav_Tagger          = fs->make<TH1F>("hUDSFlav_Tagger","Tagger",100,-50.,50.);
  hCFlav_Tagger            = fs->make<TH1F>("hCFlav_Tagger","Tagger",100,-50.,50.);
  hBFlav_Tagger            = fs->make<TH1F>("hBFlav_Tagger","Tagger",100,-50.,50.);
  
  
  
  
  
  smalltree    = fs->make<TTree>("ttree", "ttree");
  
  
  smalltree->Branch("BitTrigger"          ,&BitTrigger          ,"BitTrigger/I");
  smalltree->Branch("nSelJets"            ,&nSelJets            ,"nSelJets/I");
  smalltree->Branch("nJetCand"            ,&nJetCand            ,"nJetCand/I");
  
  
  
  if(produceJetProbaTree_){
     
    //--------------------------------------
    //tracks information 
    //--------------------------------------
    
    
    
    smalltree->Branch("nTrackCand"                        ,&nTrackCand                       ,"nTrackCand/I");
    smalltree->Branch("TrackCand_p"                       ,TrackCand_p                       ,"TrackCand_p[nTrackCand]/F");
    smalltree->Branch("TrackCand_dz"                      ,TrackCand_dz                      ,"TrackCand_dz[nTrackCand]/F");
    smalltree->Branch("TrackCand_d0"                      ,TrackCand_d0                      ,"TrackCand_d0[nTrackCand]/F");
    smalltree->Branch("TrackCand_pt"                      ,TrackCand_pt                      ,"TrackCand_pt[nTrackCand]/F");
    smalltree->Branch("TrackCand_eta"                     ,TrackCand_eta                     ,"TrackCand_eta[nTrackCand]/F");
    smalltree->Branch("TrackCand_e"                       ,TrackCand_e                       ,"TrackCand_e[nTrackCand]/F");
    smalltree->Branch("TrackCand_phi"                     ,TrackCand_phi                     ,"TrackCand_phi[nTrackCand]/F");
    
    
    smalltree->Branch("TrackCand_charge"                  ,TrackCand_charge,                  "TrackCand_charge[nTrackCand]/F");
    smalltree->Branch("TrackCand_chi2"                    ,TrackCand_chi2,                    "TrackCand_chi2[nTrackCand]/F");
    smalltree->Branch("TrackCand_Normchi2"                ,TrackCand_Normchi2,                "TrackCand_Normchi2[nTrackCand]/F");
    smalltree->Branch("TrackCand_IPsignificance"          ,TrackCand_IPsignificance,          "TrackCand_IPsignificance[nTrackCand]/F");
    smalltree->Branch("TrackCand_TransverseIPsignificance",TrackCand_TransverseIPsignificance,"TrackCand_TransverseIPsignificance[nTrackCand]/F");
    smalltree->Branch("TrackCand_TransverseIP"            ,TrackCand_TransverseIP            ,"TrackCand_TransverseIP[nTrackCand]/F");
    smalltree->Branch("TrackCand_Proba"                   ,TrackCand_Proba                   ,"TrackCand_[nTrackCand]Proba/F");
    smalltree->Branch("TrackCand_nHitPixel"               ,TrackCand_nHitPixel,               "TrackCand_nHitPixel[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitTOB"                 ,TrackCand_nHitTOB,                 "TrackCand_nHitTOB[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitTIB"                 ,TrackCand_nHitTIB,                 "TrackCand_nHitTIB[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitTID"                 ,TrackCand_nHitTID,                 "TrackCand_nHitTID[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitTEC"                 ,TrackCand_nHitTEC,                 "TrackCand_nHitTEC[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitTracker"             ,TrackCand_nHitTracker,             "TrackCand_nHitTracker[nTrackCand]/I");
    //smalltree->Branch("TrackCand_isHitL1"                 ,TrackCand_isHitL1,                 "TrackCand_isHitL1[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitPixelEC"             ,TrackCand_nHitPixelEC,             "TrackCand_nHitPixelEC[nTrackCand]/I");
    smalltree->Branch("TrackCand_nHitPixelBL"             ,TrackCand_nHitPixelBL,             "TrackCand_nHitPixelBL[nTrackCand]/I");
    
    
    smalltree->Branch("TrackCand_nHitPixelBL"             ,TrackCand_nHitPixelBL,             "TrackCand_nHitPixelBL[nTrackCand]/I");
    smalltree->Branch("JetCand_nFirstTrack",               JetCand_nFirstTrack               ,"JetCand_nFirstTrack[nJetCand]/I");
    smalltree->Branch("JetCand_nLastTrack",                JetCand_nLastTrack                ,"JetCand_nLastTrack[nJetCand]/I"); 
    smalltree->Branch("Ntagtracks",                        Ntagtracks                        ,"Ntagtracks[nJetCand]/I"); 
    
    
    //--------------------------------------
    //secondary vertex information 
    //--------------------------------------
    smalltree->Branch("nSecondaryV"            ,&nSecondaryV           ,"nSecondaryV/I");
    smalltree->Branch("SecondaryV_x"           ,SecondaryV_x           ,"SecondaryV_x[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_y"           ,SecondaryV_y           ,"SecondaryV_y[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_z"           ,SecondaryV_z           ,"SecondaryV_z[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_ex"          ,SecondaryV_ex          ,"SecondaryV_ex[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_ey"          ,SecondaryV_ey          ,"SecondaryV_ey[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_ez"          ,SecondaryV_ez          ,"SecondaryV_ez[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_chi2"        ,SecondaryV_chi2        ,"SecondaryV_chi2[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_ndf"         ,SecondaryV_ndf         ,"SecondaryV_ndf[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_flightDistance"         ,SecondaryV_flightDistance         ,"SecondaryV_flightDistance[nSecondaryV]/F");
    smalltree->Branch("SecondaryV_flightDistanceError"    ,SecondaryV_flightDistanceError    ,"SecondaryV_flightDistanceError[nSecondaryV]/F");

  }
    

  
    
  
  smalltree->Branch("JetCand_pt",               JetCand_pt               ,"JetCand_pt[nJetCand]/F");
  smalltree->Branch("JetCand_jes",              JetCand_jes              ,"JetCand_jes[nJetCand]/F");
  smalltree->Branch("JetCand_eta",              JetCand_eta              ,"JetCand_eta[nJetCand]/F");
  smalltree->Branch("JetCand_phi",              JetCand_phi              ,"JetCand_phi[nJetCand]/F");
  smalltree->Branch("JetCand_multiplicity",     JetCand_multiplicity     ,"JetCand_multiplicity[nJetCand]/I");
  smalltree->Branch("JetCand_flavor",           JetCand_flavor           ,"JetCand_flavor[nJetCand]/I");
  smalltree->Branch("JetCand_nFirstTrack",      JetCand_nFirstTrack      ,"JetCand_nFirstTrack[nJetCand]/I");
  smalltree->Branch("JetCand_nLastTrack",       JetCand_nLastTrack       ,"JetCand_nLastTrack[nJetCand]/I"); 
  smalltree->Branch("JetCand_nFirstSV"         ,JetCand_nFirstSV         ,"JetCand_nFirstSV[nJetCand]/I");
  smalltree->Branch("JetCand_nLastSV"          ,JetCand_nLastSV          ,"JetCand_nLastSV[nJetCand]/I");
  smalltree->Branch("JetCand_Ip1N",             JetCand_Ip1N             ,"JetCand_Ip1N[nJetCand]/F");
  smalltree->Branch("JetCand_Ip1P",             JetCand_Ip1P             ,"JetCand_Ip1P[nJetCand]/F");
  smalltree->Branch("JetCand_Ip2N",             JetCand_Ip2N             ,"JetCand_Ip2N[nJetCand]/F");
  smalltree->Branch("JetCand_Ip2P",             JetCand_Ip2P             ,"JetCand_Ip2P[nJetCand]/F");
  smalltree->Branch("JetCand_Ip3N",             JetCand_Ip3N             ,"JetCand_Ip3N[nJetCand]/F");
  smalltree->Branch("JetCand_Ip3P",             JetCand_Ip3P             ,"JetCand_Ip3P[nJetCand]/F");
  smalltree->Branch("JetCand_ProbaN",           JetCand_ProbaN           ,"JetCand_ProbaN[nJetCand]/F");
  smalltree->Branch("JetCand_ProbaP",           JetCand_ProbaP           ,"JetCand_ProbaP[nJetCand]/F");
  smalltree->Branch("JetCand_Proba",            JetCand_Proba            ,"JetCand_Proba[nJetCand]/F");
  smalltree->Branch("JetCand_SvtxN",            JetCand_SvtxN            ,"JetCand_SvtxN[nJetCand]/F");
  smalltree->Branch("JetCand_Svtx",             JetCand_Svtx             ,"JetCand_Svtx[nJetCand]/F");
  smalltree->Branch("JetCand_CombinedSvtxN",    JetCand_CombinedSvtxN    ,"JetCand_CombinedSvtxN[nJetCand]/F");
  smalltree->Branch("JetCand_CombinedSvtx",     JetCand_CombinedSvtx     ,"JetCand_CombinedSvtx[nJetCand]/F");
  smalltree->Branch("JetCand_SoftMN",           JetCand_SoftMN           ,"JetCand_SoftMN[nJetCand]/F");
  smalltree->Branch("JetCand_SoftM",            JetCand_SoftM            ,"JetCand_SoftM[nJetCand]/F");
  smalltree->Branch("JetCand_Category1",        JetCand_Category1        ,"JetCand_Category1[nJetCand]/F");
  smalltree->Branch("JetCand_Category2",        JetCand_Category2        ,"JetCand_Category2[nJetCand]/F");
  smalltree->Branch("JetCand_Category3",        JetCand_Category3        ,"JetCand_Category3[nJetCand]/F");
  smalltree->Branch("JetCand_CategoryJet",      JetCand_CategoryJet      ,"JetCand_CategoryJet[nJetCand]/F");
  smalltree->Branch("JetCand_CategorySVx",      JetCand_CategorySVx      ,"JetCand_CategorySVx[nJetCand]/F");
  smalltree->Branch("JetCand_CategoryMuon",     JetCand_CategoryMuon     ,"JetCand_CategoryMuon[nJetCand]/F");
  smalltree->Branch("JetCand_mu_NHits_tracker", JetCand_mu_NHits_tracker ,"JetCand_mu_NHits_tracker[nJetCand]/F");
  smalltree->Branch("JetCand_mu_Chi2",          JetCand_mu_Chi2          ,"JetCand_mu_Chi2[nJetCand]/F");
  smalltree->Branch("JetCand_mu_pT",            JetCand_mu_pT            ,"JetCand_mu_pT[nJetCand]/F");
  smalltree->Branch("JetCand_mu_d0",            JetCand_mu_d0            ,"JetCand_mu_d0[nJetCand]/F");
  smalltree->Branch("JetCand_ptRel",            JetCand_ptRel            ,"JetCand_ptRel[nJetCand]/F");
  
  
  /*Njets, BitTrigger, 
       Ijet, Flavour, Ntagtracks, ptjet, (*jet).eta(), phijet, 
       Ip1N, Ip1P, Ip2N, Ip2P, Ip3N, Ip3P, ProbaN, ProbaP, Proba,
       SvtxN, Svtx, CombinedSvtxN, CombinedSvtx, SoftMN, SoftM, 
       Category1, Category2, Category3, CategoryJet, CategorySVx, CategoryMuon,
       mu_NHits_tracker,mu_Chi2,mu_d0,mu_pT,ptRel*/

}

 
MistagAnalyzer::~MistagAnalyzer()
 
{
 

  /*hData_All_NJets       ->Write();
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
  smalltree->Write();*/
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
  

  nTrackCand  = 0;
  nSecondaryV = 0; 
    
    
    
  //------------------------------------------------------
  //Trigger info
  //------------------------------------------------------
  
  TriggerResults tr;
  Handle<TriggerResults> h_trigRes;
  //iEvent.getByLabel(InputTag("TriggerResults::HLT8E29"), h_trigRes);
  iEvent.getByLabel(InputTag("TriggerResults::HLT"), h_trigRes);
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
  nJetCand = 0;
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
    
    
    JetCand_flavor[nJetCand] = Flavour;
    JetCand_eta[nJetCand]    = (*jet).eta();
    JetCand_phi[nJetCand]    = (*jet).phi();
    JetCand_pt[nJetCand]     = (*jet).pt();
    
    
    
    
    double JES     =  acorrector->correction(*jet, iEvent, iSetup);
    JetCand_jes[nJetCand]     = JES;

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
    JetCand_nFirstTrack[nJetCand]  = nTrackCand;
    const edm::RefVector<reco::TrackCollection> & assotracks((*tagInfo)[ith_tagged].selectedTracks());
    if(produceJetProbaTree_){
      
      for(unsigned int trackIdx=0; trackIdx < assotracks.size(); trackIdx++){
	//(*tagInfo)[ith_tagged].probability(trackIdx,0);
	
	
	TrackCand_p[nTrackCand]              = (assotracks[trackIdx])->p();
	TrackCand_d0[nTrackCand]             = (assotracks[trackIdx])->d0();
	TrackCand_dz[nTrackCand]             = (assotracks[trackIdx])->dz();
	TrackCand_pt[nTrackCand]             = (assotracks[trackIdx])->pt();
	TrackCand_phi[nTrackCand]            = (assotracks[trackIdx])->phi();
	TrackCand_eta[nTrackCand]            = (assotracks[trackIdx])->eta();
	TrackCand_charge[nTrackCand]         = (assotracks[trackIdx])->charge();
	TrackCand_chi2[nTrackCand]           = (assotracks[trackIdx])->chi2();    
	TrackCand_Normchi2[nTrackCand]       = (assotracks[trackIdx])->normalizedChi2();   
	TrackCand_zIP[nTrackCand]            = (assotracks[trackIdx])->dz();//-(*pv).z();        
	
	
	TrackCand_nHitTracker[nTrackCand]          = (assotracks[trackIdx])->hitPattern().numberOfValidStripHits();
	TrackCand_nHitPixel[nTrackCand]            = (assotracks[trackIdx])->hitPattern().numberOfValidPixelHits();
	
	
	
	TrackCand_nHitTOB[nTrackCand] = (assotracks[trackIdx])->hitPattern().numberOfValidStripTOBHits();
	TrackCand_nHitTIB[nTrackCand] = (assotracks[trackIdx])->hitPattern().numberOfValidStripTIBHits();
	TrackCand_nHitTEC[nTrackCand] = (assotracks[trackIdx])->hitPattern().numberOfValidStripTECHits();
	TrackCand_nHitTID[nTrackCand] = (assotracks[trackIdx])->hitPattern().numberOfValidStripTIDHits();
	
	
	TrackCand_nHitPixelEC[nTrackCand] = (assotracks[trackIdx])->hitPattern().numberOfValidPixelEndcapHits();
	TrackCand_nHitPixelBL[nTrackCand] = (assotracks[trackIdx])->hitPattern().numberOfValidPixelBarrelHits();
       
	
	nTrackCand++;
	
      }
    }
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    JetCand_nLastTrack[nJetCand]   = nTrackCand;
    
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
    if( (*jetTags_softMuneg)[ith_tagged].second  > -100000 )   SoftMN    = ((*jetTags_softMuneg)[ith_tagged].second);
    if(SoftMN > 0) SoftMN = -1*SoftMN;
    
    
              
    JetCand_ProbaN[nJetCand]        = ProbaN;         
    JetCand_ProbaP[nJetCand]        = ProbaP;          
    JetCand_Proba[nJetCand]         = Proba;         
    JetCand_SvtxN[nJetCand]         = SvtxN;        
    JetCand_Svtx[nJetCand]          = Svtx;         
    JetCand_CombinedSvtxN[nJetCand] = CombinedSvtxN;
    JetCand_CombinedSvtx[nJetCand]  = CombinedSvtx;  
    JetCand_SoftMN[nJetCand]        = SoftMN;         
    JetCand_SoftM[nJetCand]         = SoftM;          
 
    
    
    
    
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
    
    JetCand_Ip1P[nJetCand]   = -1000;
    JetCand_Ip1N[nJetCand]  =  1000;
    int index = 0;
    for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
	 itipdata != ipdata.end(); itipdata++){
      double ip3D = (*itipdata).ip3d.significance() ;
      //negatif tracks
      if( ip3D > JetCand_Ip1P[nJetCand]  ) JetCand_Ip1P[nJetCand] = ip3D;         
      if( ip3D < JetCand_Ip1N[nJetCand]  ) JetCand_Ip1N[nJetCand] = ip3D;
      index++;
    }
    
    JetCand_Ip1N[nJetCand] = -JetCand_Ip1N[nJetCand];
    ith_tagged = this->TaggedJet(*jet,jetTags_TCHighEff);
    JetCand_Ip2P[nJetCand]   = (*jetTags_TCHighEff)[ith_tagged].second;
    ith_tagged = this->TaggedJet(*jet,jetTags_TCHighPur);
    JetCand_Ip3P[nJetCand]   = (*jetTags_TCHighPur)[ith_tagged].second;
    ith_tagged = this->TaggedJet(*jet,jetTags_NegTCHighEff);
    JetCand_Ip2N[nJetCand]   = (*jetTags_NegTCHighEff)[ith_tagged].second;
    ith_tagged = this->TaggedJet(*jet,jetTags_NegTCHighPur);
    JetCand_Ip3N[nJetCand]   = (*jetTags_NegTCHighPur)[ith_tagged].second;
    
    //*****************************************************************
    //get track histories for 1st, 2nd and 3rd track (TC)
    //*****************************************************************
    JetCand_Category1[nJetCand] = 0;
    JetCand_Category2[nJetCand] = 0;
    JetCand_Category3[nJetCand] = 0;
    
    // Track history
    if (useTrackHistory_ && indexes.size()!=0 && isData_!=0) {
      //$$ ????
      flags1P[TrackCategories::Conversion] ;
      //$$ ????
      if ( flags1P[TrackCategories::BWeakDecay] )  JetCand_Category1[nJetCand] += int(pow(10., -1 + 1)); 
      if ( flags1P[TrackCategories::CWeakDecay] )  JetCand_Category1[nJetCand] += int(pow(10., -1 + 2)); 
      if ( flags1P[TrackCategories::TauDecay] )    JetCand_Category1[nJetCand] += int(pow(10., -1 + 3)); 
      if ( flags1P[TrackCategories::Conversion] )  JetCand_Category1[nJetCand] += int(pow(10., -1 + 4)); 
      if ( flags1P[TrackCategories::KsDecay] )     JetCand_Category1[nJetCand] += int(pow(10., -1 + 5)); 
      if ( flags1P[TrackCategories::LambdaDecay] ) JetCand_Category1[nJetCand] += int(pow(10., -1 + 6)); 
      if ( flags1P[TrackCategories::Interaction] ) JetCand_Category1[nJetCand] += int(pow(10., -1 + 7)); 
      if ( flags1P[TrackCategories::Fake] )	   JetCand_Category1[nJetCand] += int(pow(10., -1 + 8)); 
      if ( flags1P[TrackCategories::Bad] )	   JetCand_Category1[nJetCand] += int(pow(10., -1 + 9)); 
      if ( idSize > 1 ) {
        if ( flags2P[TrackCategories::BWeakDecay] )  JetCand_Category2[nJetCand] += int(pow(10., -1 + 1)); 
        if ( flags2P[TrackCategories::CWeakDecay] )  JetCand_Category2[nJetCand] += int(pow(10., -1 + 2)); 
        if ( flags2P[TrackCategories::TauDecay] )    JetCand_Category2[nJetCand] += int(pow(10., -1 + 3)); 
        if ( flags2P[TrackCategories::Conversion] )  JetCand_Category2[nJetCand] += int(pow(10., -1 + 4)); 
        if ( flags2P[TrackCategories::KsDecay] )     JetCand_Category2[nJetCand] += int(pow(10., -1 + 5)); 
        if ( flags2P[TrackCategories::LambdaDecay] ) JetCand_Category2[nJetCand] += int(pow(10., -1 + 6)); 
        if ( flags2P[TrackCategories::Interaction] ) JetCand_Category2[nJetCand] += int(pow(10., -1 + 7)); 
        if ( flags2P[TrackCategories::Fake] )	     JetCand_Category2[nJetCand] += int(pow(10., -1 + 8)); 
        if ( flags2P[TrackCategories::Bad] )	     JetCand_Category2[nJetCand] += int(pow(10., -1 + 9)); 
      }
      if ( idSize > 2 ) {
        if ( flags3P[TrackCategories::BWeakDecay] )  JetCand_Category3[nJetCand] += int(pow(10., -1 + 1)); 
        if ( flags3P[TrackCategories::CWeakDecay] )  JetCand_Category3[nJetCand] += int(pow(10., -1 + 2)); 
        if ( flags3P[TrackCategories::TauDecay] )    JetCand_Category3[nJetCand] += int(pow(10., -1 + 3)); 
        if ( flags3P[TrackCategories::Conversion] )  JetCand_Category3[nJetCand] += int(pow(10., -1 + 4)); 
        if ( flags3P[TrackCategories::KsDecay] )     JetCand_Category3[nJetCand] += int(pow(10., -1 + 5)); 
        if ( flags3P[TrackCategories::LambdaDecay] ) JetCand_Category3[nJetCand] += int(pow(10., -1 + 6)); 
        if ( flags3P[TrackCategories::Interaction] ) JetCand_Category3[nJetCand] += int(pow(10., -1 + 7)); 
        if ( flags3P[TrackCategories::Fake] )	     JetCand_Category3[nJetCand] += int(pow(10., -1 + 8)); 
        if ( flags3P[TrackCategories::Bad] )	     JetCand_Category3[nJetCand] += int(pow(10., -1 + 9)); 
      }
      if ( flags1N[TrackCategories::BWeakDecay] )  JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 1)); 
      if ( flags1N[TrackCategories::CWeakDecay] )  JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 2)); 
      if ( flags1N[TrackCategories::TauDecay] )    JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 3)); 
      if ( flags1N[TrackCategories::Conversion] )  JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 4)); 
      if ( flags1N[TrackCategories::KsDecay] )     JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 5)); 
      if ( flags1N[TrackCategories::LambdaDecay] ) JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 6)); 
      if ( flags1N[TrackCategories::Interaction] ) JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 7)); 
      if ( flags1N[TrackCategories::Fake] )	   JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 8)); 
      if ( flags1N[TrackCategories::Bad] )	   JetCand_Category1[nJetCand] += 2*int(pow(10., -1 + 9)); 
      if ( idSize > 1 ) {
        if ( flags2N[TrackCategories::BWeakDecay] )  JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 1)); 
        if ( flags2N[TrackCategories::CWeakDecay] )  JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 2)); 
        if ( flags2N[TrackCategories::TauDecay] )    JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 3)); 
        if ( flags2N[TrackCategories::Conversion] )  JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 4)); 
        if ( flags2N[TrackCategories::KsDecay] )     JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 5)); 
        if ( flags2N[TrackCategories::LambdaDecay] ) JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 6)); 
        if ( flags2N[TrackCategories::Interaction] ) JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 7)); 
        if ( flags2N[TrackCategories::Fake] )	     JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 8)); 
        if ( flags2N[TrackCategories::Bad] )	     JetCand_Category2[nJetCand] += 2*int(pow(10., -1 + 9)); 
      }
      if ( idSize > 2 ) {
        if ( flags3N[TrackCategories::BWeakDecay] )  JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 1)); 
        if ( flags3N[TrackCategories::CWeakDecay] )  JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 2)); 
        if ( flags3N[TrackCategories::TauDecay] )    JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 3)); 
        if ( flags3N[TrackCategories::Conversion] )  JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 4)); 
        if ( flags3N[TrackCategories::KsDecay] )     JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 5)); 
        if ( flags3N[TrackCategories::LambdaDecay] ) JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 6)); 
        if ( flags3N[TrackCategories::Interaction] ) JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 7)); 
        if ( flags3N[TrackCategories::Fake] )	     JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 8)); 
        if ( flags3N[TrackCategories::Bad] )	     JetCand_Category3[nJetCand] += 2*int(pow(10., -1 + 9)); 
      }
    }
    
    
    
    
    
    
    
    
    
    
    
    //*****************************************************************
    //get track histories of tracks in jets (for Jet Proba) 
    //*****************************************************************
    JetCand_CategoryJet[nJetCand] = 0;
    
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
      JetCand_CategoryJet[nJetCand] =  cap0+can0           + (cap1+can1)*10    + (cap2+can2)*100 
                  + (cap3+can3)*1000     + (cap4+can4)*10000 
                  + (cap5+can5)*100000   + (cap6+can6)*1000000 
		  + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }
    
    //*****************************************************************
    //get track histories  associated to sec. vertex (for simple SV)
    //*****************************************************************
    JetCand_CategorySVx[nJetCand] = 0;
    
    JetCand_nFirstSV[nJetCand]  = nSecondaryV;
    
    
    
    ith_tagged =    this->TaggedJet(*jet,jetTags_Svtx);
    TrackRefVector  svxPostracks( (*tagInfoSVx)[ith_tagged].vertexTracks(0) );
    
    if(produceJetProbaTree_){
      
      
      SecondaryV_x[nSecondaryV]    = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).x();
	SecondaryV_y[nSecondaryV]    = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).y();
	SecondaryV_z[nSecondaryV]    = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).z();
	SecondaryV_ex[nSecondaryV]   = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).xError();
	SecondaryV_ey[nSecondaryV]   = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).yError();
	SecondaryV_ez[nSecondaryV]   = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).zError();
	SecondaryV_chi2[nSecondaryV] = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).chi2();
	SecondaryV_ndf[nSecondaryV]  = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).ndof();
	
	SecondaryV_flightDistance[nSecondaryV]      = (*tagInfoSVx)[ith_tagged].flightDistance(0).value();
	SecondaryV_flightDistanceError[nSecondaryV] = (*tagInfoSVx)[ith_tagged].flightDistance(0).error();
	
	nSecondaryV++;
    }
    
    
    
    
    
    
    if(useTrackHistory_ && isData_!=0) {
      
      
      
      
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
      JetCand_CategorySVx[nJetCand] =  cap0+can0           + (cap1+can1)*10    + (cap2+can2)*100 
                  + (cap3+can3)*1000     + (cap4+can4)*10000 
                  + (cap5+can5)*100000   + (cap6+can6)*1000000 
		  + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }
    
    
    JetCand_nLastSV[nJetCand]   = nSecondaryV;
    JetCand_mu_NHits_tracker[nJetCand]  = -10000;
    JetCand_mu_Chi2[nJetCand]           = -10000;
    JetCand_mu_d0[nJetCand]             = -10000;
    JetCand_mu_pT[nJetCand]             = -10000;
    JetCand_ptRel[nJetCand]             = -10000;
    
    
    //*****************************************************************
    //get track histories of the muon (SoftMuon tagger)
    //*****************************************************************
    JetCand_CategoryMuon[nJetCand] = 0;

    ith_tagged = this->TaggedJet(*jet,jetTags_softMu);
    if( (*tagInos_softmuon)[ith_tagged].leptons()!=0 ) {
      JetCand_ptRel[nJetCand]             = calculPtRel( (*(*tagInos_softmuon)[ith_tagged].lepton(0)), *jet, JES);
      JetCand_mu_NHits_tracker[nJetCand]  = (*tagInos_softmuon)[ith_tagged].lepton(0)->hitPattern().numberOfValidHits();
      JetCand_mu_Chi2[nJetCand]           = (*tagInos_softmuon)[ith_tagged].lepton(0)->normalizedChi2()                ;
      JetCand_mu_d0[nJetCand]             = (*tagInos_softmuon)[ith_tagged].lepton(0)->d0()                            ;
      JetCand_mu_pT[nJetCand]             = (*tagInos_softmuon)[ith_tagged].lepton(0)->pt()                            ;
    }
    if(  SoftM > 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ){
      if ( useTrackHistory_ && isData_!=0 ) {     
        TrackCategories::Flags theFlagP = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
        if ( theFlagP[TrackCategories::BWeakDecay] )   JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 1)); 
        if ( theFlagP[TrackCategories::CWeakDecay] )   JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 2)); 
        if ( theFlagP[TrackCategories::TauDecay] )     JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 3));  
        if ( theFlagP[TrackCategories::Conversion] )   JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 4)); 
        if ( theFlagP[TrackCategories::KsDecay] )      JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 5)); 
        if ( theFlagP[TrackCategories::LambdaDecay] )  JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 6)); 
        if ( theFlagP[TrackCategories::Interaction] )  JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 7)); 
        if ( theFlagP[TrackCategories::Fake] )	       JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 8)); 
        if ( theFlagP[TrackCategories::Bad] )	       JetCand_CategoryMuon[nJetCand] += int(pow(10., -1 + 9)); 
      }
    }

    ith_tagged = this->TaggedJet(*jet,jetTags_softMuneg);
     //std::cout << "SoftMN " << SoftMN << std::endl;
    if ( SoftMN < 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ) {
      JetCand_ptRel[nJetCand] =  calculPtRel( (*(*tagInos_softmuon)[ith_tagged].lepton(0)), *jet, JES);
      JetCand_mu_NHits_tracker[nJetCand]  = (*tagInos_softmuon)[ith_tagged].lepton(0)->hitPattern().numberOfValidHits();
      JetCand_mu_Chi2[nJetCand]           = (*tagInos_softmuon)[ith_tagged].lepton(0)->normalizedChi2()                ;
      JetCand_mu_d0[nJetCand]             = (*tagInos_softmuon)[ith_tagged].lepton(0)->d0()                            ;
      JetCand_mu_pT[nJetCand]             = (*tagInos_softmuon)[ith_tagged].lepton(0)->pt()                            ;
      if ( useTrackHistory_ && isData_!=0 ) {     
        TrackCategories::Flags theFlagN = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
        if ( theFlagN[TrackCategories::BWeakDecay] )   JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 1)); 
        if ( theFlagN[TrackCategories::CWeakDecay] )   JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 2)); 
        if ( theFlagN[TrackCategories::TauDecay] )     JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 3));  
        if ( theFlagN[TrackCategories::Conversion] )   JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 4)); 
        if ( theFlagN[TrackCategories::KsDecay] )      JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 5)); 
        if ( theFlagN[TrackCategories::LambdaDecay] )  JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 6)); 
        if ( theFlagN[TrackCategories::Interaction] )  JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 7)); 
        if ( theFlagN[TrackCategories::Fake] )	       JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 8)); 
        if ( theFlagN[TrackCategories::Bad] )	       JetCand_CategoryMuon[nJetCand] += 2*int(pow(10., -1 + 9)); 
      }
    }
        
    //*********************************
    // Jet selection
    
    if ( (jetpt * JES) <= minJetPt_ || std::fabs( jeteta ) >= maxJetEta_ ) continue;
    numjet++;
    nSelJets    =       numjet;
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
 
    /*float jet_input[34] =
      {Njets, BitTrigger, 
       Ijet, Flavour, Ntagtracks, ptjet, (*jet).eta(), phijet, 
       Ip1N, Ip1P, Ip2N, Ip2P, Ip3N, Ip3P, ProbaN, ProbaP, Proba,
       SvtxN, Svtx, CombinedSvtxN, CombinedSvtx, SoftMN, SoftM, 
       Category1, Category2, Category3, CategoryJet, CategorySVx, CategoryMuon,
       mu_NHits_tracker,mu_Chi2,mu_d0,mu_pT,ptRel
       };*/
    
    
    // nTuplesJets->Fill(jet_input);
    
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
      if ( JetCand_Ip1P[nJetCand] > 0. ) varpos = JetCand_Ip1P[nJetCand];
      if ( JetCand_Ip1N[nJetCand] > 0. ) varneg = JetCand_Ip1N[nJetCand];
    }
    else if ( selTagger_ == 2 ) {     // TC High Eff.
      if ( JetCand_Ip2P[nJetCand] > 0. ) varpos = JetCand_Ip2P[nJetCand];
      if ( JetCand_Ip2N[nJetCand] > 0. ) varneg = JetCand_Ip2N[nJetCand];
    }
    else if ( selTagger_ == 3 ) {     // TC High Pure.
      if ( JetCand_Ip3P[nJetCand] > 0. ) varpos = JetCand_Ip3P[nJetCand];
      if ( JetCand_Ip3N[nJetCand] > 0. ) varneg = JetCand_Ip3N[nJetCand];
    }
    else if ( selTagger_ == 4 ) {    // SV simple
      if  (JetCand_Svtx[nJetCand] > 1)   varpos =  10.*JetCand_Svtx[nJetCand] - 10.;
      if  (JetCand_SvtxN[nJetCand] < -1) varneg = -10.*JetCand_SvtxN[nJetCand] - 10.;
      //       if ( Svtx  > tagCut_ ) TagPos = true;
      //       if (-SvtxN > tagCut_ ) TagNeg = true; 
    }
    else if ( selTagger_ == 5 ) {    // SV combined
      if  (JetCand_CombinedSvtx[nJetCand] > 0)  varpos = 50.*JetCand_CombinedSvtx[nJetCand];
      if  (JetCand_CombinedSvtxN[nJetCand] > 0) varneg = 50.*JetCand_CombinedSvtxN[nJetCand];
      //       if ( CombinedSvtx  > tagCut_ ) TagPos = true;
      //       if (-CombinedSvtxN > tagCut_ ) TagNeg = true;   
    }
    else if ( selTagger_ == 6 ) {    // soft muon ptrel
      if  ( JetCand_SoftM[nJetCand]  > 0) varpos =  5*JetCand_SoftM[nJetCand];
      if  ( JetCand_SoftMN[nJetCand] < 0) varneg = -5*JetCand_SoftMN[nJetCand];
      //       if  (  SoftM > tagCut_ ) TagPos = true;
      //       if  (  SoftMN> tagCut_ ) TagNeg = true;   
    }
    
    //     if ( selTagger_ >= 1 && selTagger_ <= 3) {
    //       if ( varpos > tagCut_ ) TagPos = true;
    //       if ( varneg > tagCut_ ) TagNeg = true;
    //     }
    
    //     // veto on positive tag
    //     if ( JetCand_Ip1P[nJetCand] < vetoPos_ ) Veto = true;
    
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

    if ( JetCand_Ip2P[nJetCand] > 0 )   hData_Tagger_TCHE->Fill( JetCand_Ip2P[nJetCand] );
    if ( JetCand_Ip2N[nJetCand] > 0 )   hData_Tagger_TCHE->Fill(-JetCand_Ip2N[nJetCand] );
    if ( JetCand_Ip3P[nJetCand] > 0 )   hData_Tagger_TCHP->Fill( JetCand_Ip3P[nJetCand] );
    if ( JetCand_Ip3N[nJetCand] > 0 )   hData_Tagger_TCHP->Fill(-JetCand_Ip3N[nJetCand] );
    if ( JetCand_ProbaP[nJetCand] > 0 ) hData_Tagger_JP->Fill( 20.*JetCand_ProbaP[nJetCand] );
    if ( JetCand_ProbaN[nJetCand] > 0 ) hData_Tagger_JP->Fill(-20.*JetCand_ProbaN[nJetCand] );
    if ( JetCand_Svtx[nJetCand]  >  1 ) hData_Tagger_SSV->Fill( 10.*JetCand_Svtx[nJetCand] - 10 );
    if ( JetCand_SvtxN[nJetCand] < -1 ) hData_Tagger_SSV->Fill( 10*JetCand_SvtxN[nJetCand] + 10 );
    if ( JetCand_CombinedSvtx[nJetCand] > 0  )  hData_Tagger_CSV->Fill( 50.*JetCand_CombinedSvtx[nJetCand] );
    if ( JetCand_CombinedSvtxN[nJetCand] > 0 )  hData_Tagger_CSV->Fill(-50.*JetCand_CombinedSvtxN[nJetCand] );
    if ( JetCand_SoftM[nJetCand]  > 0        )  hData_Tagger_MU->Fill( 5*JetCand_SoftM[nJetCand] );
    if ( JetCand_SoftMN[nJetCand] < 0        )  hData_Tagger_MU->Fill( 5*JetCand_SoftMN[nJetCand] );

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
    
    
    
    nJetCand++;
    
  }// end loop on jet
  smalltree->Fill();
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


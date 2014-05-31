#include "RecoBTag/PerformanceMeasurements/interface/BTagAnalyzer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

PFJetIDSelectionFunctor pfjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE );
PFJetIDSelectionFunctor pfjetIDTight( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT );

pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

BTagAnalyzer::BTagAnalyzer(const edm::ParameterSet& iConfig):
  classifier_(iConfig) ,
  pv(0),
  PFJet80(0),
  thelepton1(0.,0.,0.,0.),
  thelepton2(0.,0.,0.,0.),
  computer(0),
  nsubjettinessCalculator(Njettiness::onepass_kt_axes, NsubParameters(1.0, 0.8, 0.8)),
  cap0(0),
  cap1(0),
  cap2(0),
  cap3(0),
  cap4(0),
  cap5(0),
  cap6(0),
  cap7(0),
  cap8(0),
  can0(0),
  can1(0),
  can2(0),
  can3(0),
  can4(0),
  can5(0),
  can6(0),
  can7(0),
  can8(0),
  hadronizerType_(0)
{
  //now do what ever initialization you need
  std::string module_type  = iConfig.getParameter<std::string>("@module_type");
  std::string module_label = iConfig.getParameter<std::string>("@module_label");
  cout << "Constructing " << module_type << ":" << module_label << endl;

  // Parameters
  runSubJets_ = iConfig.getParameter<bool>("runSubJets");
  allowJetSkipping_ = iConfig.getParameter<bool>("allowJetSkipping");
  minJetPt_  = iConfig.getParameter<double>("MinPt");
  maxJetEta_ = iConfig.getParameter<double>("MaxEta");

  selTagger_ = iConfig.getParameter<int>("selTagger");

  useSelectedTracks_    = iConfig.getParameter<bool> ("useSelectedTracks");
  useTrackHistory_      = iConfig.getParameter<bool> ("useTrackHistory");
  produceJetTrackTree_  = iConfig.getParameter<bool> ("produceJetTrackTree");
  produceAllTrackTree_  = iConfig.getParameter<bool> ("produceAllTrackTree");
  producePtRelTemplate_ = iConfig.getParameter<bool> ("producePtRelTemplate");

  use_ttbar_filter_     = iConfig.getParameter<bool> ("use_ttbar_filter");
  channel_              = iConfig.getParameter<edm::InputTag> ("channel");

  // Modules
  primaryVertexColl_   = iConfig.getParameter<edm::InputTag>("primaryVertexColl");
  tracksColl_          = iConfig.getParameter<edm::InputTag>("tracksColl");

  src_                 = iConfig.getParameter<edm::InputTag>("src");

  JetCollectionTag_ = iConfig.getParameter<edm::InputTag>("Jets");
  FatJetCollectionTag_ = iConfig.getParameter<edm::InputTag>("FatJets");
  PrunedFatJetCollectionTag_ = iConfig.getParameter<edm::InputTag>("PrunedFatJets");

  trackCHEBJetTags_    = iConfig.getParameter<std::string>("trackCHEBJetTags");
  trackCNegHEBJetTags_ = iConfig.getParameter<std::string>("trackCNegHEBJetTags");

  trackCHPBJetTags_    = iConfig.getParameter<std::string>("trackCHPBJetTags");
  trackCNegHPBJetTags_ = iConfig.getParameter<std::string>("trackCNegHPBJetTags");

  jetPBJetTags_        = iConfig.getParameter<std::string>("jetPBJetTags");
  jetPNegBJetTags_     = iConfig.getParameter<std::string>("jetPNegBJetTags");
  jetPPosBJetTags_     = iConfig.getParameter<std::string>("jetPPosBJetTags");

  jetBPBJetTags_        = iConfig.getParameter<std::string>("jetBPBJetTags");
  jetBPNegBJetTags_     = iConfig.getParameter<std::string>("jetBPNegBJetTags");
  jetBPPosBJetTags_     = iConfig.getParameter<std::string>("jetBPPosBJetTags");

  simpleSVHighEffBJetTags_      = iConfig.getParameter<std::string>("simpleSVHighEffBJetTags");
  simpleSVNegHighEffBJetTags_   = iConfig.getParameter<std::string>("simpleSVNegHighEffBJetTags");
  simpleSVHighPurBJetTags_      = iConfig.getParameter<std::string>("simpleSVHighPurBJetTags");
  simpleSVNegHighPurBJetTags_   = iConfig.getParameter<std::string>("simpleSVNegHighPurBJetTags");

  combinedSVBJetTags_     = iConfig.getParameter<std::string>("combinedSVBJetTags");
  combinedSVNegBJetTags_  = iConfig.getParameter<std::string>("combinedSVNegBJetTags");
  combinedSVPosBJetTags_  = iConfig.getParameter<std::string>("combinedSVPosBJetTags");

  combinedSVRetrainedBJetTags_     = iConfig.getParameter<std::string>("combinedSVRetrainedBJetTags");
  combinedSVRetrainedNegBJetTags_  = iConfig.getParameter<std::string>("combinedSVRetrainedNegBJetTags");
  combinedSVRetrainedPosBJetTags_  = iConfig.getParameter<std::string>("combinedSVRetrainedPosBJetTags");

  combinedCSVJPBJetTags_     = iConfig.getParameter<std::string>("combinedCSVJPBJetTags");
  combinedCSVJPNegBJetTags_  = iConfig.getParameter<std::string>("combinedCSVJPNegBJetTags");
  combinedCSVJPPosBJetTags_  = iConfig.getParameter<std::string>("combinedCSVJPPosBJetTags");

  combinedCSVSLBJetTags_     = iConfig.getParameter<std::string>("combinedCSVSLBJetTags");
  combinedCSVSLNegBJetTags_  = iConfig.getParameter<std::string>("combinedCSVSLNegBJetTags");
  combinedCSVSLPosBJetTags_  = iConfig.getParameter<std::string>("combinedCSVSLPosBJetTags");

  combinedCSVJPSLBJetTags_     = iConfig.getParameter<std::string>("combinedCSVJPSLBJetTags");
  combinedCSVJPSLNegBJetTags_  = iConfig.getParameter<std::string>("combinedCSVJPSLNegBJetTags");
  combinedCSVJPSLPosBJetTags_  = iConfig.getParameter<std::string>("combinedCSVJPSLPosBJetTags");

  simpleIVFSVHighPurBJetTags_ = iConfig.getParameter<std::string>("simpleIVFSVHighPurBJetTags");
  simpleIVFSVHighEffBJetTags_ = iConfig.getParameter<std::string>("simpleIVFSVHighEffBJetTags");
  doubleIVFSVHighEffBJetTags_ = iConfig.getParameter<std::string>("doubleIVFSVHighEffBJetTags");
  combinedIVFSVBJetTags_      = iConfig.getParameter<std::string>("combinedIVFSVBJetTags");
  combinedIVFSVPosBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVPosBJetTags");

  //softMuonBJetTags_       = iConfig.getParameter<std::string>("softMuonBJetTags");
  //softMuonNegBJetTags_    = iConfig.getParameter<std::string>("softMuonNegBJetTags");
  //softMuonTagInfoName_    = iConfig.getParameter<std::string>("softMuonTagInfoName");

  softPFMuonBJetTags_       = iConfig.getParameter<std::string>("softPFMuonBJetTags");
  softPFMuonNegBJetTags_    = iConfig.getParameter<std::string>("softPFMuonNegBJetTags");
  softPFMuonPosBJetTags_    = iConfig.getParameter<std::string>("softPFMuonPosBJetTags");

  softPFElectronBJetTags_       = iConfig.getParameter<std::string>("softPFElectronBJetTags");
  softPFElectronNegBJetTags_    = iConfig.getParameter<std::string>("softPFElectronNegBJetTags");
  softPFElectronPosBJetTags_    = iConfig.getParameter<std::string>("softPFElectronPosBJetTags");

  softPFMuonTagInfos_      = iConfig.getParameter<std::string>("softPFMuonTagInfos");
  softPFElectronTagInfos_  = iConfig.getParameter<std::string>("softPFElectronTagInfos");

  muonCollectionName_       = iConfig.getParameter<edm::InputTag>("muonCollectionName");
  patMuonCollectionName_    = iConfig.getParameter<edm::InputTag>("patMuonCollectionName");

  triggerTable_             = iConfig.getParameter<edm::InputTag>("triggerTable");

  SVComputer_               = iConfig.getParameter<edm::InputTag>( "svComputer");

  triggerPathNames_        = iConfig.getParameter<std::vector<std::string> >("TriggerPathNames");
  ttbarTriggerPathNames_   = iConfig.getParameter<std::vector<std::string> >("TTbarTriggerPathNames");
  PFJet80TriggerPathNames_ = iConfig.getParameter<std::vector<std::string> >("PFJet80TriggerPathNames");

  ///////////////
  // TTree

  smalltree = fs->make<TTree>("ttree", "ttree");

  //--------------------------------------
  // event information
  //--------------------------------------
  EventInfo.RegisterTree(smalltree);
  if ( runSubJets_ )          EventInfo.RegisterPatMuonTree(smalltree);
  if ( use_ttbar_filter_ )    EventInfo.RegisterTTbarTree(smalltree);
  if ( produceJetTrackTree_ ) EventInfo.RegisterJetTrackTree(smalltree);
  if ( produceAllTrackTree_ ) EventInfo.RegisterAllTrackTree(smalltree);

  //--------------------------------------
  // jet information
  //--------------------------------------
  JetInfo[0].RegisterTree(smalltree,(runSubJets_ ? "JetInfo" : ""));
  if ( runSubJets_ )          JetInfo[0].RegisterSubJetSpecificTree(smalltree,(runSubJets_ ? "JetInfo" : ""));
  if ( produceJetTrackTree_ ) JetInfo[0].RegisterJetTrackTree(smalltree,(runSubJets_ ? "JetInfo" : ""));
  if ( runSubJets_ ) {
    JetInfo[1].RegisterTree(smalltree,"FatJetInfo");
    JetInfo[1].RegisterFatJetSpecificTree(smalltree,"FatJetInfo");
    if ( produceJetTrackTree_ ) JetInfo[1].RegisterJetTrackTree(smalltree,"FatJetInfo");
  }

  //// Book Histograms
  Histos[0] = new BookHistograms(fs->mkdir( "HistJets" )) ;
  if (runSubJets_) Histos[1] = new BookHistograms(fs->mkdir( "HistFatJets" )) ;

  cout << module_type << ":" << module_label << " constructed" << endl;
}


BTagAnalyzer::~BTagAnalyzer()
{
}


static std::vector<std::size_t> sortedIndexes(const std::vector<TrackIPTagInfo::TrackIPData >& values)
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
void BTagAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  //------------------------------------------------------
  //Event informations
  //------------------------------------------------------
  EventInfo.Run = iEvent.id().run();
  isData_ = iEvent.isRealData();

  if ( !isData_ && EventInfo.Run > 0 ) EventInfo.Run = -EventInfo.Run;

  EventInfo.Evt  = iEvent.id().event();
  EventInfo.LumiBlock  = iEvent.luminosityBlock();

  // Tag Jets
  if ( useTrackHistory_ ) classifier_.newEvent(iEvent, iSetup);

  edm::Handle <PatJetCollection> jetsColl;
  iEvent.getByLabel (JetCollectionTag_, jetsColl);

  edm::Handle <PatJetCollection> fatjetsColl;
  edm::Handle <PatJetCollection> prunedfatjetsColl;
  if (runSubJets_) {
    iEvent.getByLabel(FatJetCollectionTag_, fatjetsColl) ;
    iEvent.getByLabel(PrunedFatJetCollectionTag_, prunedfatjetsColl) ;
  }

  JetToJetMap fatJetToPrunedFatJetMap;
  if( runSubJets_ )
  {
    for(PatJetCollection::const_iterator it = fatjetsColl->begin(); it != fatjetsColl->end(); ++it)
    {
      PatJetCollection::const_iterator prunedJetMatch;
      bool prunedJetMatchFound = false;
      float dR = 0.8;
      for(PatJetCollection::const_iterator pjIt = prunedfatjetsColl->begin(); pjIt != prunedfatjetsColl->end(); ++pjIt)
      {
        float dR_temp = reco::deltaR( it->p4(), pjIt->p4() );
        if( dR_temp < dR )
        {
          prunedJetMatchFound = true;
          dR = dR_temp;
          prunedJetMatch = pjIt;
        }
      }
      if( !prunedJetMatchFound ) edm::LogError("NoMatchingGroomedJet") << "Matching pruned jet not found."; // This should never happen but just in case
      fatJetToPrunedFatJetMap[&(*it)] = &(*prunedJetMatch);
    }
  }

  //------------------------------------------------------
  // Determine hadronizer type (done only once per job)
  //------------------------------------------------------
  if( !isData_ && hadronizerType_ == 0 )
  {
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(src_, genEvtInfoProduct);

    std::string moduleName = "";
    const edm::Provenance& prov = iEvent.getProvenance(genEvtInfoProduct.id());
    if( genEvtInfoProduct.isValid() )
      moduleName = prov.moduleName();

    if( moduleName.find("Pythia8")!=std::string::npos )
      hadronizerType_ |= ( 1 << 1 ); // set the 2nd bit
    else // assuming Pythia6
      hadronizerType_ |= ( 1 << 0 ); // set the 1st bit
  }

  //------------------------------------------------------
  // MC informations
  //------------------------------------------------------
  EventInfo.pthat      = -1.;
  EventInfo.nPUtrue    = -1.;
  EventInfo.nPU        = 0;
  EventInfo.ncQuarks   = 0;
  EventInfo.nbQuarks   = 0;
  EventInfo.nBHadrons  = 0;
  EventInfo.nDHadrons  = 0;
  EventInfo.nDaughters = 0;
  EventInfo.nGenV0     = 0;
  EventInfo.nGenlep    = 0;
  EventInfo.nGenquark  = 0;
  EventInfo.nTrkAll    = 0;
  EventInfo.nPatMuon   = 0;
  EventInfo.mcweight   = 1.;

  bool AreBHadrons = false;


  //---------------------------- Start MC info ---------------------------------------//
  if ( !isData_ ) {
    // EventInfo.pthat
    edm::Handle<GenEventInfoProduct> geninfos;
    iEvent.getByLabel( "generator",geninfos );
    EventInfo.mcweight=geninfos->weight();
    if (geninfos->binningValues().size()>0) EventInfo.pthat = geninfos->binningValues()[0];
    // pileup

    edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator ipu;
    for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
      for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
        EventInfo.PU_bunch[EventInfo.nPU]      =  ipu->getBunchCrossing();
        EventInfo.PU_z[EventInfo.nPU]          = (ipu->getPU_zpositions())[i];
        EventInfo.PU_sumpT_low[EventInfo.nPU]  = (ipu->getPU_sumpT_lowpT())[i];
        EventInfo.PU_sumpT_high[EventInfo.nPU] = (ipu->getPU_sumpT_highpT())[i];
        EventInfo.PU_ntrks_low[EventInfo.nPU]  = (ipu->getPU_ntrks_lowpT())[i];
        EventInfo.PU_ntrks_high[EventInfo.nPU] = (ipu->getPU_ntrks_highpT())[i];
        ++EventInfo.nPU;
      }
      EventInfo.nPUtrue = ipu->getTrueNumInteractions();
      if(EventInfo.nPU==0) EventInfo.nPU = ipu->getPU_NumInteractions(); // needed in case getPU_zpositions() is empty
    }

  //------------------------------------------------------
  // generated particles
  //------------------------------------------------------
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel ("genParticles", genParticles);

    EventInfo.GenPVz = -1000.;

    for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      int ID = abs(genIt.pdgId());
      unsigned int nDaughters = genIt.numberOfDaughters();

      if ( isHardProcess(genIt.status()) ) EventInfo.GenPVz = genIt.vz();

      // prompt parton
      if ( isHardProcess(genIt.status()) && genIt.p4().pt() > 10. ) {
        EventInfo.Genquark_pT[EventInfo.nGenquark]     = genIt.p4().pt();
        EventInfo.Genquark_eta[EventInfo.nGenquark]    = genIt.p4().eta();
        EventInfo.Genquark_phi[EventInfo.nGenquark]    = genIt.p4().phi();
        EventInfo.Genquark_pdgID[EventInfo.nGenquark]  = genIt.pdgId();
        if (genIt.numberOfMothers()>0) EventInfo.Genquark_mother[EventInfo.nGenquark] = genIt.mother()->pdgId();   // protection
        else EventInfo.Genquark_mother[EventInfo.nGenquark] = -99999;
        //  cout << " status " << genIt.status() << " pdgId " << genIt.pdgId()  << " pT " << genIt.p4().pt() << " mother " << moth->pdgId() << endl;
        ++EventInfo.nGenquark;
      }

      // b and c quarks from the end of parton showering and before hadronization
      if ( ID == 4 || ID == 5 ) {
        if( nDaughters > 0 ) {
          int nparton_daughters = 0;
          for (unsigned int d=0; d<nDaughters; ++d) {
            int daughterID = abs(genIt.daughter(d)->pdgId());
            if( (daughterID == 1 || daughterID == 2 || daughterID == 3 ||
                 daughterID == 4 || daughterID == 5 || daughterID == 6 || daughterID == 21))
              nparton_daughters++;
          }
          if( nparton_daughters == 0 ) {
            if ( ID == 5 ) {
              EventInfo.bQuark_pT[EventInfo.nbQuarks]  = genIt.p4().pt();
              EventInfo.bQuark_eta[EventInfo.nbQuarks] = genIt.p4().eta();
              EventInfo.bQuark_phi[EventInfo.nbQuarks] = genIt.p4().phi();
              EventInfo.bQuark_pdgID[EventInfo.nbQuarks] = genIt.pdgId();
              EventInfo.bQuark_status[EventInfo.nbQuarks] = genIt.status();
              EventInfo.bQuark_fromGSP[EventInfo.nbQuarks] = isFromGSP(&genIt);
              ++EventInfo.nbQuarks;
            }
            if ( ID == 4 ) {
              EventInfo.cQuark_pT[EventInfo.ncQuarks]  = genIt.p4().pt();
              EventInfo.cQuark_eta[EventInfo.ncQuarks] = genIt.p4().eta();
              EventInfo.cQuark_phi[EventInfo.ncQuarks] = genIt.p4().phi();
              EventInfo.cQuark_pdgID[EventInfo.ncQuarks] = genIt.pdgId();
              EventInfo.cQuark_status[EventInfo.ncQuarks] = genIt.status();
              EventInfo.cQuark_fromGSP[EventInfo.ncQuarks] = isFromGSP(&genIt);
              ++EventInfo.ncQuarks;
            }
          }
        }
      }

      if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) AreBHadrons = true;
//       // Primary b Hadrons
//       if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
//         //  cout << " pdgId " << genIt.pdgId()  << " pT " << genIt.p4().pt() << " mother " << mother->pdgId() << endl;
//         EventInfo.BHadron_pT[EventInfo.nBHadrons]    = genIt.p4().pt();
//         EventInfo.BHadron_eta[EventInfo.nBHadrons]   = genIt.p4().eta();
//         EventInfo.BHadron_phi[EventInfo.nBHadrons]   = genIt.p4().phi();
//         EventInfo.BHadron_mass[EventInfo.nBHadrons]  = genIt.mass();
//         EventInfo.BHadron_pdgID[EventInfo.nBHadrons] = genIt.pdgId();
//         EventInfo.BHadron_status[EventInfo.nBHadrons] = genIt.status();
//         EventInfo.BHadron_mother[EventInfo.nBHadrons] = genIt.mother()->pdgId();
//         // check if any of the daughters is also B hadron
//         int hasBHadronDaughter = 0;
//         for (unsigned int d=0; d<nDaughters; ++d) {
//           int daughterID = abs(genIt.daughter(d)->pdgId());
//           if ( (daughterID/100)%10 == 5 || (daughterID/1000)%10 == 5 ) { hasBHadronDaughter = 1; break; }
//         }
//         EventInfo.BHadron_hasBdaughter[EventInfo.nBHadrons] = hasBHadronDaughter;
//         ++EventInfo.nBHadrons;
//       }

      // Final c Hadrons
      if ( (ID/100)%10 == 4 || (ID/1000)%10 == 4 ) {
        // check if daughter is not a DHadron
	bool hasNoDHadronAsDaughter = true;
	for (unsigned int d=0; d<nDaughters; d++) {
	  const Candidate* daughter = genIt.daughter(d);
	  int IDdaughter = abs(daughter->pdgId());
	  bool hasDHadron = (IDdaughter/100)%10 == 4 || (IDdaughter/1000)%10 == 4;
	  if ( hasDHadron ) hasNoDHadronAsDaughter = false;
	}

	// variables
	if ( hasNoDHadronAsDaughter ) {
	  EventInfo.DHadron_pT[EventInfo.nDHadrons]    = genIt.p4().pt();
	  EventInfo.DHadron_eta[EventInfo.nDHadrons]   = genIt.p4().eta();
	  EventInfo.DHadron_phi[EventInfo.nDHadrons]   = genIt.p4().phi();
	  EventInfo.DHadron_mass[EventInfo.nDHadrons]  = genIt.mass();
	  EventInfo.DHadron_pdgID[EventInfo.nDHadrons] = genIt.pdgId();

	  int nCharged = 0;
	  // Loop on daughters
	  int numberChargedDaughters = 0;
	  for (unsigned int d=0; d<nDaughters; d++) {
	    const Candidate* daughter = genIt.daughter(d);
	    int IDdaughter = abs(daughter->pdgId());
	    EventInfo.DHadron_DaughtersPdgID[EventInfo.nDaughters] = IDdaughter;
	    // take vertex of first daughter for SV
	    if ( d == 0 ) {
	      EventInfo.DHadron_SVx[EventInfo.nDHadrons] = daughter->vx();
   	      EventInfo.DHadron_SVy[EventInfo.nDHadrons] = daughter->vy();
    	      EventInfo.DHadron_SVz[EventInfo.nDHadrons] = daughter->vz();
	    }
	    EventInfo.nDaughters++;
	    if ( daughter->charge() != 0 ) numberChargedDaughters++;
	    // charged track multiplicity
	    unsigned int nDaughters2 = daughter->numberOfDaughters();
	    if (nDaughters2 == 0 && daughter->charge() != 0
				 && daughter->p4().pt()>1.) nCharged++;
	    for (unsigned int d2=0; d2<nDaughters2; d2++) {
	      const Candidate* daughter2 = daughter->daughter(d2);
	      unsigned int nDaughters3 = daughter2->numberOfDaughters();
	      if (nDaughters3 == 0 && daughter2->charge() != 0
	        		   && daughter2->p4().pt()>1.) nCharged++;
	      for (unsigned int d3=0; d3<nDaughters3; d3++) {
		const Candidate* daughter3 = daughter2->daughter(d3);
		unsigned int nDaughters4 = daughter3->numberOfDaughters();
		if (nDaughters4 == 0 && daughter3->charge() != 0
				     && daughter3->p4().pt()>1.) nCharged++;
		for (unsigned int d4=0; d4<nDaughters4; d4++) {
		  const Candidate* daughter4 = daughter3->daughter(d4);
		  unsigned int nDaughters5 = daughter4->numberOfDaughters();
		  if (nDaughters5 == 0 && daughter4->charge() != 0
				       && daughter4->p4().pt()>1.) nCharged++;
		} // D -> X -> X -> X -> X
	      } // D -> X -> X -> X
	    } // D -> X -> X
	  } // loop on daughters from final c hadron

	  EventInfo.DHadron_nDaughters[EventInfo.nDHadrons] = nDaughters;
	  EventInfo.DHadron_nChargedDaughters[EventInfo.nDHadrons] = numberChargedDaughters;
	  EventInfo.DHadron_nCharged[EventInfo.nDHadrons] = nCharged;
// cout << " D hadron " << EventInfo.nDHadrons << "  id " << EventInfo.DHadron_pdgID[EventInfo.nDHadrons] << "  pT " << EventInfo.DHadron_pT[EventInfo.nDHadrons] << "  eta " << EventInfo.DHadron_eta[EventInfo.nDHadrons] << "  phi " << EventInfo.DHadron_phi[EventInfo.nDHadrons] << endl;
	  EventInfo.nDHadrons++;
	} // final c hadron
      } // c hadron

      // Leptons
      if ( (ID == 11 || ID == 13) && genIt.p4().pt() > 3. ) {
       if (genIt.numberOfMothers()>0) {  // protection for herwig ttbar mc
        const Candidate * moth1 = genIt.mother();
        if ( moth1->pdgId() != genIt.pdgId() ) {
          EventInfo.Genlep_pT[EventInfo.nGenlep]    = genIt.p4().pt();
          EventInfo.Genlep_eta[EventInfo.nGenlep]   = genIt.p4().eta();
          EventInfo.Genlep_phi[EventInfo.nGenlep]   = genIt.p4().phi();
          EventInfo.Genlep_pdgID[EventInfo.nGenlep] = genIt.pdgId();
          EventInfo.Genlep_status[EventInfo.nGenlep] = genIt.status();
          int ID1 = abs( moth1->pdgId() );
          int ID2 = -9999;
          int ID3 = -9999;
          int ID4 = -9999;
          int isWZ = 0, istau = 0, isB = 0, isD = 0;
          if (moth1->numberOfMothers()>0) {   // protection for herwig ttbar mc
            const Candidate * moth2 = moth1->mother();
            ID2 = abs( moth2->pdgId() );
            if (moth2->numberOfMothers()>0) {
              const Candidate * moth3 = moth2->mother();
              ID3 = abs( moth3->pdgId() );
              if (moth3->numberOfMothers()>0) {
                const Candidate * moth4 = moth3->mother();
                ID4 = abs( moth4->pdgId() );
              }
            }
          }
          if ( ID1 == 15 ) istau = 1;
          if ( ID1 == 24 || ID2 == 24 || ID3 == 24 || ID4 == 24 ||
              ID1 == 25 || ID2 == 25 || ID3 == 25 || ID4 == 25 ) isWZ = 1;
          if ( ID1 == 411 || ID1 == 421 || ID1 == 431 ||
              ID1 ==4112 || ID1 ==4122 || ID1 ==4132 ||
              ID1 ==4212 || ID1 ==4222 || ID1 ==4232 || ID1 ==4332 ) isD = 1;
          if ( ID2 == 411 || ID2 == 421 || ID2 == 431 ||
              ID2 ==4112 || ID2 ==4122 || ID2 ==4132 ||
              ID2 ==4212 || ID2 ==4222 || ID2 ==4232 || ID2 ==4332 ) isD = 1;
          if ( ID1 == 511 || ID1 == 521 || ID1 == 531 || ID1 == 541 ||
              ID1 ==5112 || ID1 ==5122 || ID1 ==5132 ||
              ID1 ==5212 || ID1 ==5222 || ID1 ==5232 || ID1 ==5332 ) isB = 1;
          if ( ID2 == 511 || ID2 == 521 || ID2 == 531 || ID2 == 541 ||
              ID2 ==5112 || ID2 ==5122 || ID2 ==5132 ||
              ID2 ==5212 || ID2 ==5222 || ID2 ==5232 || ID2 ==5332 ) isB = 1;
          if ( ID3 == 511 || ID3 == 521 || ID3 == 531 || ID3 == 541 ||
              ID3 ==5112 || ID3 ==5122 || ID3 ==5132 ||
              ID3 ==5212 || ID3 ==5222 || ID3 ==5232 || ID3 ==5332 ) isB = 1;
          if ( ID4 == 511 || ID4 == 521 || ID4 == 531 || ID4 == 541 ||
              ID4 ==5112 || ID4 ==5122 || ID4 ==5132 ||
              ID4 ==5212 || ID4 ==5222 || ID4 ==5232 || ID4 ==5332 ) isB = 1;
          if ( isB+isD != 0 ) EventInfo.Genlep_mother[EventInfo.nGenlep] = 5*isB + 4*isD;
          else                EventInfo.Genlep_mother[EventInfo.nGenlep] = 10*istau + 100*isWZ;
          //  cout << " lepton " << EventInfo.nGenlep << " pdgID " << EventInfo.Genlep_pdgID[EventInfo.nGenlep]
          //       << " moth1 " << moth1->pdgId() << " moth2 " << moth2->pdgId()
          //       << " moth3 " << moth3->pdgId() << " moth4 " << moth4->pdgId()
          //       << " pT " << EventInfo.Genlep_pT[EventInfo.nGenlep]
          //       << " from " << EventInfo.Genlep_mother[EventInfo.nGenlep] << endl;
          ++EventInfo.nGenlep;
        }
       }
      }

      // V0: K0s and Lambda
      if ( ID == 310 || ID == 3122 ) {
	 int nCharged = 0;
	 float SVx = 0., SVy = 0., SVz = 0.;
	 for (unsigned int d=0; d<nDaughters; d++) {
	   const Candidate* daughter = genIt.daughter(d);
	   if ( daughter->charge() != 0 && daughter->pt() > 1. ) {
	     nCharged++;
   	     SVx = daughter->vx();
   	     SVy = daughter->vy();
   	     SVz = daughter->vz();
	   }
	 }
	 bool inAcceptance = false;
	 if ( nCharged > 0 ) {
	   float Radius = TMath::Sqrt( SVx*SVx + SVy*SVy );
	   float AbsEta = abs( genIt.p4().eta() );
	   if ( (AbsEta < 2.0 && Radius < 7.3) ||
	        (AbsEta < 2.5 && Radius < 4.4) ||
	        (AbsEta > 1.8 && AbsEta < 2.5 && abs(SVz) < 34.5 ) )
	     inAcceptance = true;
	 }
	 if ( inAcceptance ) {
	   EventInfo.GenV0_pT[EventInfo.nGenV0]    = genIt.p4().pt();
	   EventInfo.GenV0_eta[EventInfo.nGenV0]   = genIt.p4().eta();
	   EventInfo.GenV0_phi[EventInfo.nGenV0]   = genIt.p4().phi();
	   EventInfo.GenV0_pdgID[EventInfo.nGenV0] = genIt.pdgId();
	   EventInfo.GenV0_SVx[EventInfo.nGenV0] = SVx;
	   EventInfo.GenV0_SVy[EventInfo.nGenV0] = SVy;
   	   EventInfo.GenV0_SVz[EventInfo.nGenV0] = SVz;
	   EventInfo.GenV0_nCharged[EventInfo.nGenV0] = nCharged;
	   EventInfo.nGenV0++;
	 }
      }

    } // end loop on generated particles

  //------------------------------------------------------
  // generated particles: looking for b hadrons only
  //------------------------------------------------------
    if ( AreBHadrons ) {
    for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      int ID = abs(genIt.pdgId());
      unsigned int nDaughters = genIt.numberOfDaughters();

      // b Hadrons
      if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
        //  cout << " pdgId " << genIt.pdgId()  << " pT " << genIt.p4().pt() << " mother " << mother->pdgId() << endl;
        EventInfo.BHadron_pT[EventInfo.nBHadrons]    = genIt.p4().pt();
        EventInfo.BHadron_eta[EventInfo.nBHadrons]   = genIt.p4().eta();
        EventInfo.BHadron_phi[EventInfo.nBHadrons]   = genIt.p4().phi();
        EventInfo.BHadron_mass[EventInfo.nBHadrons]  = genIt.mass();
        EventInfo.BHadron_pdgID[EventInfo.nBHadrons] = genIt.pdgId();
        if (genIt.numberOfMothers()>0) EventInfo.BHadron_mother[EventInfo.nBHadrons] = genIt.mother()->pdgId();
        else EventInfo.BHadron_mother[EventInfo.nBHadrons] = -99999;
        // check if any of the daughters is also B hadron
        int hasBHadronDaughter = 0;
        for (unsigned int d=0; d<nDaughters; ++d) {
          int daughterID = abs(genIt.daughter(d)->pdgId());
          if ( (daughterID/100)%10 == 5 || (daughterID/1000)%10 == 5 ) { hasBHadronDaughter = 1; break; }
        }
        EventInfo.BHadron_hasBdaughter[EventInfo.nBHadrons] = hasBHadronDaughter;
	int idD1 = -1, idD2 = -1;
	int nCharged = 0;
	// take vertex of first daughter for SV
  	EventInfo.BHadron_SVx[EventInfo.nBHadrons] = genIt.daughter(0)->vx();
   	EventInfo.BHadron_SVy[EventInfo.nBHadrons] = genIt.daughter(0)->vy();
    	EventInfo.BHadron_SVz[EventInfo.nBHadrons] = genIt.daughter(0)->vz();

	if ( !hasBHadronDaughter ) {
	  // Loop on daughters: B -> X
	  for (unsigned int d=0; d<nDaughters; d++) {
	    const Candidate* daughter = genIt.daughter(d);
	    int IDdaughter = abs(daughter->pdgId());
	    bool isDHadron = (IDdaughter/100)%10 == 4 || (IDdaughter/1000)%10 == 4;
	    bool hasNoDHadronAsDaughter = true;
	    bool isDHadron2, hasNoDHadronAsDaughter2, isDHadron3;
	    unsigned int nDaughters2 = daughter->numberOfDaughters();

            // B -> D
	    if ( isDHadron ) {
	      for (unsigned int d2=0; d2<nDaughters2; d2++) {
	        const Candidate* daughter2 = daughter->daughter(d2);
	        int IDdaughter2 = abs(daughter2->pdgId());
	        isDHadron2 = (IDdaughter2/100)%10 == 4 || (IDdaughter2/1000)%10 == 4;
	        if ( isDHadron2 ) hasNoDHadronAsDaughter = false;
	      }
	      if ( hasNoDHadronAsDaughter ) {
	        for (int c=0; c<EventInfo.nDHadrons; c++) {
	          if ( abs( daughter->p4().pt()
		                   - EventInfo.DHadron_pT[c] ) < 0.0001 ) {
		    if ( idD1 < 0 ) idD1 = c;
		    else if ( idD2 < 0 ) idD2 = c;
// cout << " D hadron from B " << c << "  id " << daughter->pdgId() << "  pT " << daughter->p4().pt()  << "  eta " << daughter->p4().eta() << "  phi " << daughter->p4().phi() << endl;
		  }
	        }
	      }
	    }

	    if ( !isDHadron || !hasNoDHadronAsDaughter ) {
	      if (nDaughters2 == 0 && daughter->charge() != 0
	  	    		   && daughter->p4().pt() > 1.) nCharged++;
	      // B -> X -> X
	      for (unsigned int d2=0; d2<nDaughters2; d2++) {
	        const Candidate* daughter2 = daughter->daughter(d2);
	        int IDdaughter2 = abs(daughter2->pdgId());
	        isDHadron2 = (IDdaughter2/100)%10 == 4 || (IDdaughter2/1000)%10 == 4;
	        hasNoDHadronAsDaughter2 = true;
	        unsigned int nDaughters3 = daughter2->numberOfDaughters();

	        // B -> D*(*) -> D
	        if ( isDHadron2 ) {
	          for (unsigned int d3=0; d3<nDaughters3; d3++) {
	            const Candidate* daughter3 = daughter2->daughter(d3);
	            int IDdaughter3 = abs(daughter3->pdgId());
	            isDHadron3 = (IDdaughter3/100)%10 == 4 || (IDdaughter3/1000)%10 == 4;
	            if ( isDHadron3 ) hasNoDHadronAsDaughter2 = false;
	          }
	          if ( hasNoDHadronAsDaughter2 ) {
	            for (int c=0; c<EventInfo.nDHadrons; c++) {
	              if ( abs( daughter2->p4().pt()
		                       - EventInfo.DHadron_pT[c] ) < 0.0001 ) {
		        if ( idD1 < 0 ) idD1 = c;
		        else if ( idD2 < 0 ) idD2 = c;
// cout << " D hadron from B->D*(*) " << c << "  id " << daughter2->pdgId() << "  pT " << daughter2->p4().pt() << "  eta " << daughter2->p4().eta() << "  phi " << daughter2->p4().phi() << endl;
		      }
	            }
	          }
	        }

	        if ( !isDHadron2 || !hasNoDHadronAsDaughter2 ) {
	          if (nDaughters3 == 0 && daughter2->charge() != 0
	        		       && daughter2->p4().pt() > 1.) nCharged++;
	          // B -> X -> X -> X
	          for (unsigned int d3=0; d3<nDaughters3; d3++) {
	            const Candidate* daughter3 = daughter2->daughter(d3);
	            int IDdaughter3 = abs(daughter3->pdgId());
	            isDHadron3 = (IDdaughter3/100)%10 == 4 || (IDdaughter3/1000)%10 == 4;
	            unsigned int nDaughters4 = daughter3->numberOfDaughters();

	            // B -> D** -> D* -> D
	            if ( isDHadron3 ) {
	              for (int c=0; c<EventInfo.nDHadrons; c++) {
	                if ( abs( daughter3->p4().pt()
		                         - EventInfo.DHadron_pT[c] ) < 0.0001 ) {
		          if ( idD1 < 0 ) idD1 = c;
		          else if ( idD2 < 0 ) idD2 = c;
// cout << " D hadron from B->D**->D* " << c << "  id " << daughter3->pdgId() << "  pT " << daughter3->p4().pt() << "  eta " << daughter3->p4().eta() << "  phi " << daughter3->p4().phi() << endl;
		        }
	              }
	            }
                    else {
	              if (nDaughters4 == 0 && daughter3->charge() != 0
	        	  	           && daughter3->p4().pt() > 1.) nCharged++;
	              // B -> X -> X -> X -> X
	              for (unsigned int d4=0; d4<nDaughters4; d4++) {
	                const Candidate* daughter4 = daughter3->daughter(d4);
	                unsigned int nDaughters5 = daughter4->numberOfDaughters();
	                if (nDaughters5 == 0 && daughter4->charge() != 0
	        	   	             && daughter4->p4().pt()>1.) nCharged++;
	              }
	            }
	          } // B -> X -> X -> X
	        }
	      } // B -> X -> X
	    }
	  } // loop on daughters from final b hadron
// cout << " B hadron " << EventInfo.nBHadrons << "  id " << EventInfo.BHadron_pdgID[EventInfo.nBHadrons] << "  pT " << EventInfo.BHadron_pT[EventInfo.nBHadrons] << "  eta " << EventInfo.BHadron_eta[EventInfo.nBHadrons]<< "  phi " << EventInfo.BHadron_phi[EventInfo.nBHadrons] << "  ->D " << idD1 << " " << idD2 << "  nCh " << nCharged << endl;
	} // final b hadron

    	EventInfo.BHadron_nCharged[EventInfo.nBHadrons] = nCharged;
    	EventInfo.BHadron_DHadron1[EventInfo.nBHadrons] = idD1;
    	EventInfo.BHadron_DHadron2[EventInfo.nBHadrons] = idD2;
 	EventInfo.nBHadrons++;
      } // b hadron
    } // end loop on generated particles
    }

  //------------------------------------------------------
  // simulated PV
  //------------------------------------------------------
    if ( useTrackHistory_ ) {
      edm::Handle<edm::HepMCProduct> theHepMCProduct;
      iEvent.getByLabel("generator",theHepMCProduct);
      std::vector<simPrimaryVertex> simpv;
      simpv=getSimPVs(theHepMCProduct);
      //       cout << "simpv.size() " << simpv.size() << endl;
    }
  }
  //---------------------------- End MC info ---------------------------------------//
  //   cout << "EventInfo.Evt:" <<EventInfo.Evt << endl;
  //   cout << "EventInfo.pthat:" <<EventInfo.pthat << endl;


  //------------------------------------------------------
  // ttbar information
  //------------------------------------------------------
  if (use_ttbar_filter_) {
    edm::Handle<int> pIn;
    iEvent.getByLabel(channel_, pIn);
    EventInfo.ttbar_chan=*pIn;
    edm::Handle<vector< double  >> pIn2;
    iEvent.getByLabel(channel_, pIn2);
    if (pIn2->size()==8) {
      EventInfo.lepton1_pT=(*pIn2)[0];
      EventInfo.lepton1_eta=(*pIn2)[1];
      EventInfo.lepton1_phi=(*pIn2)[2];
      EventInfo.lepton2_pT=(*pIn2)[3];
      EventInfo.lepton2_eta=(*pIn2)[4];
      EventInfo.lepton2_phi=(*pIn2)[5];
      EventInfo.met=(*pIn2)[6];
      EventInfo.mll=(*pIn2)[7];
    }
    else {
      EventInfo.lepton1_pT=-1;
      EventInfo.lepton1_eta=-1;
      EventInfo.lepton1_phi=-1;
      EventInfo.lepton2_pT=-1;
      EventInfo.lepton2_eta=-1;
      EventInfo.lepton2_phi=-1;
      EventInfo.met=-1;
      EventInfo.mll=-1;
    }
  }

  //------------------------------------------------------
  // PAT Muons
  //------------------------------------------------------
  edm::Handle<std::vector<pat::Muon> >  patMuonsHandle;
  if( runSubJets_ )
  {
    iEvent.getByLabel(patMuonCollectionName_,patMuonsHandle);

    for( std::vector<pat::Muon>::const_iterator it = patMuonsHandle->begin(); it != patMuonsHandle->end(); ++it )
    {
      if( !it->isGlobalMuon() ) continue;

      EventInfo.PatMuon_isGlobal[EventInfo.nPatMuon] = 1;
      EventInfo.PatMuon_isPF[EventInfo.nPatMuon]     = it->isPFMuon();
      EventInfo.PatMuon_nTkHit[EventInfo.nPatMuon]   = it->innerTrack()->hitPattern().numberOfValidHits();
      EventInfo.PatMuon_nPixHit[EventInfo.nPatMuon]  = it->innerTrack()->hitPattern().numberOfValidPixelHits();
      EventInfo.PatMuon_nOutHit[EventInfo.nPatMuon]  = it->innerTrack()->trackerExpectedHitsOuter().numberOfHits();
      EventInfo.PatMuon_nMuHit[EventInfo.nPatMuon]   = it->outerTrack()->hitPattern().numberOfValidMuonHits();
      EventInfo.PatMuon_nMatched[EventInfo.nPatMuon] = it->numberOfMatches();
      EventInfo.PatMuon_chi2[EventInfo.nPatMuon]     = it->globalTrack()->normalizedChi2();
      EventInfo.PatMuon_chi2Tk[EventInfo.nPatMuon]   = it->innerTrack()->normalizedChi2();
      EventInfo.PatMuon_pt[EventInfo.nPatMuon]       = it->pt();
      EventInfo.PatMuon_eta[EventInfo.nPatMuon]      = it->eta();
      EventInfo.PatMuon_phi[EventInfo.nPatMuon]      = it->phi();
      EventInfo.PatMuon_vz[EventInfo.nPatMuon]       = it->vz();
      EventInfo.PatMuon_IP[EventInfo.nPatMuon]       = it->dB(pat::Muon::PV3D);
      EventInfo.PatMuon_IPsig[EventInfo.nPatMuon]    = (it->dB(pat::Muon::PV3D))/(it->edB(pat::Muon::PV3D));
      EventInfo.PatMuon_IP2D[EventInfo.nPatMuon]     = it->dB(pat::Muon::PV2D);
      EventInfo.PatMuon_IP2Dsig[EventInfo.nPatMuon]  = (it->dB(pat::Muon::PV2D))/(it->edB(pat::Muon::PV2D));

      ++EventInfo.nPatMuon;
    }
  }


  //------------------------------------------------------
  // Muons
  //------------------------------------------------------
  edm::Handle<edm::View<reco::Muon> >  muonsHandle;
  iEvent.getByLabel(muonCollectionName_,muonsHandle);
  muons = *muonsHandle;

  //----------------------------------------
  // Transient track for IP calculation
  //----------------------------------------
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  //------------------
  // Primary vertex
  //------------------
  iEvent.getByLabel(primaryVertexColl_,primaryVertex);

  //bool newvertex = false;

  bool pvFound = (primaryVertex->size() != 0);
  if ( pvFound ) {
    pv = &(*primaryVertex->begin());
  }
  else {
    reco::Vertex::Error e;
    e(0,0)=0.0015*0.0015;
    e(1,1)=0.0015*0.0015;
    e(2,2)=15.*15.;
    reco::Vertex::Point p(0,0,0);
    pv=  new reco::Vertex(p,e,1,1,1);
    //newvertex = true;
  }
  //   GlobalPoint Pv_point = GlobalPoint((*pv).x(), (*pv).y(), (*pv).z());
  EventInfo.PVz = (*primaryVertex)[0].z();
  EventInfo.PVez = (*primaryVertex)[0].zError();

  EventInfo.nPV=0;
  for (unsigned int i = 0; i< primaryVertex->size() ; ++i) {
    EventInfo.PV_x[EventInfo.nPV]      = (*primaryVertex)[i].x();
    EventInfo.PV_y[EventInfo.nPV]      = (*primaryVertex)[i].y();
    EventInfo.PV_z[EventInfo.nPV]      = (*primaryVertex)[i].z();
    EventInfo.PV_ex[EventInfo.nPV]     = (*primaryVertex)[i].xError();
    EventInfo.PV_ey[EventInfo.nPV]     = (*primaryVertex)[i].yError();
    EventInfo.PV_ez[EventInfo.nPV]     = (*primaryVertex)[i].zError();
    EventInfo.PV_chi2[EventInfo.nPV]   = (*primaryVertex)[i].normalizedChi2();
    EventInfo.PV_ndf[EventInfo.nPV]    = (*primaryVertex)[i].ndof();
    EventInfo.PV_isgood[EventInfo.nPV] = (*primaryVertex)[i].isValid();
    EventInfo.PV_isfake[EventInfo.nPV] = (*primaryVertex)[i].isFake();

    ++EventInfo.nPV;
  }

  //------------------------------------------------------
  // Trigger info
  //------------------------------------------------------

  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByLabel(triggerTable_, trigRes);

  EventInfo.nBitTrigger = int(triggerPathNames_.size()/32)+1;
  for(int i=0; i<EventInfo.nBitTrigger; ++i) EventInfo.BitTrigger[i] = 0;
  if (use_ttbar_filter_) EventInfo.trig_ttbar = 0;

  std::vector<std::string> triggerList;
  edm::Service<edm::service::TriggerNamesService> tns;
  bool foundNames = tns->getTrigPaths(*trigRes,triggerList);
  if ( !foundNames ) edm::LogError("TriggerNamesNotFound") << "Could not get trigger names!";
  if ( trigRes->size() != triggerList.size() ) edm::LogError("TriggerPathLengthMismatch") << "Length of names and paths not the same: "
    << triggerList.size() << "," << trigRes->size() ;

  processTrig(trigRes, triggerList);

  PFJet80 = false;
  for(std::vector<std::string>::const_iterator itTrigPathNames = PFJet80TriggerPathNames_.begin();
      itTrigPathNames != PFJet80TriggerPathNames_.end(); ++itTrigPathNames)
  {
    std::vector<std::string>::const_iterator it = std::find(triggerPathNames_.begin(), triggerPathNames_.end(), *itTrigPathNames);
    if( it != triggerPathNames_.end() )
    {
      int triggerIdx = ( it - triggerPathNames_.begin() );
      int bitIdx = int(triggerIdx/32);
      if ( EventInfo.BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) PFJet80 = true;
    }
  }

  //------------- added by Camille-----------------------------------------------------------//
  edm::ESHandle<JetTagComputer> computerHandle;
  iSetup.get<JetTagComputerRecord>().get( SVComputer_.label(), computerHandle );

  computer = dynamic_cast<const GenericMVAJetTagComputer*>( computerHandle.product() );

  computer->passEventSetup(iSetup);
  //------------- end added-----------------------------------------------------------//

  if (use_ttbar_filter_) {
    thelepton1.SetPtEtaPhiM(EventInfo.lepton1_pT, EventInfo.lepton1_eta, EventInfo.lepton1_phi, 0.);
    thelepton2.SetPtEtaPhiM(EventInfo.lepton2_pT, EventInfo.lepton2_eta, EventInfo.lepton2_phi, 0.);
  }

  //------------------------------------------------------
  // All tracks info
  //------------------------------------------------------
  edm::Handle<reco::TrackCollection> tracksHandle;
  if( produceAllTrackTree_ )
  {
    iEvent.getByLabel(tracksColl_,tracksHandle);

    for(reco::TrackCollection::const_iterator trk = tracksHandle->begin(); trk != tracksHandle->end(); ++trk)
    {
      EventInfo.TrkAll_d0[EventInfo.nTrkAll]        = trk->d0();
      EventInfo.TrkAll_dz[EventInfo.nTrkAll]        = trk->dz();
      EventInfo.TrkAll_p[EventInfo.nTrkAll]         = trk->p();
      EventInfo.TrkAll_pt[EventInfo.nTrkAll]        = trk->pt();
      EventInfo.TrkAll_eta[EventInfo.nTrkAll]       = trk->eta();
      EventInfo.TrkAll_phi[EventInfo.nTrkAll]       = trk->phi();
      EventInfo.TrkAll_chi2[EventInfo.nTrkAll]      = trk->normalizedChi2();
      EventInfo.TrkAll_charge[EventInfo.nTrkAll]    = trk->charge();
      EventInfo.TrkAll_nHitAll[EventInfo.nTrkAll]   = trk->numberOfValidHits();
      EventInfo.TrkAll_nHitPixel[EventInfo.nTrkAll] = trk->hitPattern().numberOfValidPixelHits();
      EventInfo.TrkAll_nHitStrip[EventInfo.nTrkAll] = trk->hitPattern().numberOfValidStripHits();
      EventInfo.TrkAll_nHitTIB[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTIBHits();
      EventInfo.TrkAll_nHitTID[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTIDHits();
      EventInfo.TrkAll_nHitTOB[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTOBHits();
      EventInfo.TrkAll_nHitTEC[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidStripTECHits();
      EventInfo.TrkAll_nHitPXB[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidPixelBarrelHits();
      EventInfo.TrkAll_nHitPXF[EventInfo.nTrkAll]   = trk->hitPattern().numberOfValidPixelEndcapHits();
      EventInfo.TrkAll_isHitL1[EventInfo.nTrkAll]   = trk->hitPattern().hasValidHitInFirstPixelBarrel();
      EventInfo.TrkAll_nSiLayers[EventInfo.nTrkAll] = trk->hitPattern().trackerLayersWithMeasurement();
      EventInfo.TrkAll_nPxLayers[EventInfo.nTrkAll] = trk->hitPattern().pixelLayersWithMeasurement();

      reco::TrackRef trkRef(tracksHandle, trk - tracksHandle->begin() );

      setTracksPV(trkRef, primaryVertex,
                  EventInfo.TrkAll_PV[EventInfo.nTrkAll],
                  EventInfo.TrkAll_PVweight[EventInfo.nTrkAll]);
      ++EventInfo.nTrkAll;
    }
  }
  //------------------------------------------------------

  //------------------------------------------------------
  // Jet info
  //------------------------------------------------------
  int iJetColl = 0 ;
  //// Do jets
  processJets(jetsColl, fatjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;
  if (runSubJets_) {
    iJetColl = 1 ;
    processJets(fatjetsColl, jetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;
  }
  //------------------------------------------------------

  //// Fill TTree
  if ( EventInfo.BitTrigger > 0 || EventInfo.Run < 0 ) {
    smalltree->Fill();
  }

  return;
}


void BTagAnalyzer::processTrig(const edm::Handle<edm::TriggerResults>& trigRes, const std::vector<std::string>& triggerList)
{
  for (unsigned int i = 0; i < trigRes->size(); ++i) {

    if ( !trigRes->at(i).accept() ) continue;


    for (std::vector<std::string>::const_iterator itTrigPathNames = triggerPathNames_.begin();
        itTrigPathNames != triggerPathNames_.end(); ++itTrigPathNames)
    {
      int triggerIdx = ( itTrigPathNames - triggerPathNames_.begin() );
      int bitIdx = int(triggerIdx/32);
      if ( NameCompatible(*itTrigPathNames,triggerList[i]) ) EventInfo.BitTrigger[bitIdx] |= ( 1 << (triggerIdx - bitIdx*32) );
    }

    if (use_ttbar_filter_)
    {
      for (std::vector<std::string>::const_iterator itTrigPathNames = ttbarTriggerPathNames_.begin();
          itTrigPathNames != ttbarTriggerPathNames_.end(); ++itTrigPathNames)
      {
        if ( NameCompatible(*itTrigPathNames,triggerList[i]) ) EventInfo.trig_ttbar |= ( 1 << ( itTrigPathNames - ttbarTriggerPathNames_.begin() ) );
      }
    } //// if use_ttbar_filter_
  } //// Loop over trigger names

  return;
}


void BTagAnalyzer::processJets(const edm::Handle<PatJetCollection>& jetsColl, const edm::Handle<PatJetCollection>& jetsColl2,
    const edm::Event& iEvent, const edm::EventSetup& iSetup, const JetToJetMap& fatJetToPrunedFatJetMap, const int iJetColl)
{

  int numjet = 0;
  JetInfo[iJetColl].nMuon = 0;
  JetInfo[iJetColl].nPFElectron = 0;
  JetInfo[iJetColl].nPFMuon = 0;

  JetInfo[iJetColl].nJet = 0;
  JetInfo[iJetColl].nTrack = 0;
  JetInfo[iJetColl].nTrkInc = 0;
  JetInfo[iJetColl].nSV = 0;

  //// Loop over the jets
  for ( PatJetCollection::const_iterator pjet = jetsColl->begin(); pjet != jetsColl->end(); ++pjet ) {

    double ptjet  = pjet->pt()  ;
    double jeteta = pjet->eta() ;
    double jetphi = pjet->phi() ;

    if ( allowJetSkipping_ && ( ptjet < minJetPt_ || std::fabs( jeteta ) > maxJetEta_ ) ) continue;

    //// overlap removal with lepton from ttbar selection
    if (use_ttbar_filter_) {
      TLorentzVector thejet;
      thejet.SetPtEtaPhiM(ptjet, jeteta, jetphi, 0.);
      double deltaR1 = thejet.DeltaR(thelepton1);
      double deltaR2 = thejet.DeltaR(thelepton2);
      if (EventInfo.ttbar_chan>=0 && (deltaR1 < 0.5 || deltaR2 < 0.5)) continue;
    }
    //// end of removal

    int flavour  =-1  ;
    if ( !isData_ ) {
      flavour = abs( pjet->partonFlavour() );
      if ( flavour >= 1 && flavour <= 3 ) flavour = 1;
    }

    JetInfo[iJetColl].Jet_flavour[JetInfo[iJetColl].nJet] = flavour;
    JetInfo[iJetColl].Jet_eta[JetInfo[iJetColl].nJet]     = pjet->eta();
    JetInfo[iJetColl].Jet_phi[JetInfo[iJetColl].nJet]     = pjet->phi();
    JetInfo[iJetColl].Jet_pt[JetInfo[iJetColl].nJet]      = pjet->pt();
    JetInfo[iJetColl].Jet_mass[JetInfo[iJetColl].nJet]    = pjet->mass();
    JetInfo[iJetColl].Jet_genpt[JetInfo[iJetColl].nJet]   = ( pjet->genJet()!=0 ? pjet->genJet()->pt() : -1. );

    retpf.set(false);
    JetInfo[iJetColl].Jet_looseID[JetInfo[iJetColl].nJet]  = ( pfjetIDLoose( *pjet, retpf ) ? 1 : 0 );
    retpf.set(false);
    JetInfo[iJetColl].Jet_tightID[JetInfo[iJetColl].nJet]  = ( pfjetIDTight( *pjet, retpf ) ? 1 : 0 );

    JetInfo[iJetColl].Jet_jes[JetInfo[iJetColl].nJet]      = pjet->pt()/pjet->correctedJet("Uncorrected").pt();
    JetInfo[iJetColl].Jet_residual[JetInfo[iJetColl].nJet] = pjet->pt()/pjet->correctedJet("L3Absolute").pt();

    if( runSubJets_ && iJetColl == 0 )
    {
      int fatjetIdx=-1;
      for( PatJetCollection::const_iterator jIt = jetsColl2->begin(); jIt != jetsColl2->end(); ++jIt )
      {
        if( &(*pjet) == fatJetToPrunedFatJetMap.find(&(*jIt))->second->daughter(0) ||
            &(*pjet) == fatJetToPrunedFatJetMap.find(&(*jIt))->second->daughter(1) )
        {
          fatjetIdx = int( jIt - jetsColl2->begin() );
          break;
        }
      }
      JetInfo[iJetColl].Jet_FatJetIdx[JetInfo[iJetColl].nJet] = fatjetIdx;

      if( ptjet==0. ) // special treatment for pT=0 subjets
      {
        ++numjet;
        ++JetInfo[iJetColl].nJet;
        continue;
      }
    }

    int subjet1Idx = -1, subjet2Idx = -1;

    if ( runSubJets_ && iJetColl == 1 )
    {
      // N-subjettiness
      std::vector<fastjet::PseudoJet> fjConstituents;
      std::vector<edm::Ptr<reco::PFCandidate> > constituents = pjet->getPFConstituents();
      std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator m;
      for ( m = constituents.begin(); m != constituents.end(); ++m )
      {
        reco::PFCandidatePtr constit = *m;
        if (constit->pt() == 0)
        {
          edm::LogWarning("NullTransverseMomentum") << "dropping input candidate with pt=0";
          continue;
        }
        fjConstituents.push_back(fastjet::PseudoJet(constit->px(),constit->py(),constit->pz(),constit->energy()));
        fjConstituents.back().set_user_index(m - constituents.begin());
      }

      JetInfo[iJetColl].Jet_tau1[JetInfo[iJetColl].nJet] = nsubjettinessCalculator.getTau(1,fjConstituents);
      JetInfo[iJetColl].Jet_tau2[JetInfo[iJetColl].nJet] = nsubjettinessCalculator.getTau(2,fjConstituents);

      JetInfo[iJetColl].Jet_ptPruned[JetInfo[iJetColl].nJet]   = fatJetToPrunedFatJetMap.find(&(*pjet))->second->pt();
      JetInfo[iJetColl].Jet_jesPruned[JetInfo[iJetColl].nJet]  = fatJetToPrunedFatJetMap.find(&(*pjet))->second->pt()/fatJetToPrunedFatJetMap.find(&(*pjet))->second->correctedJet("Uncorrected").pt();
      JetInfo[iJetColl].Jet_etaPruned[JetInfo[iJetColl].nJet]  = fatJetToPrunedFatJetMap.find(&(*pjet))->second->eta();
      JetInfo[iJetColl].Jet_phiPruned[JetInfo[iJetColl].nJet]  = fatJetToPrunedFatJetMap.find(&(*pjet))->second->phi();
      JetInfo[iJetColl].Jet_massPruned[JetInfo[iJetColl].nJet] = fatJetToPrunedFatJetMap.find(&(*pjet))->second->mass();

      for( PatJetCollection::const_iterator jIt = jetsColl2->begin(); jIt != jetsColl2->end(); ++jIt )
      {
        if( &(*jIt) == fatJetToPrunedFatJetMap.find(&(*pjet))->second->daughter(0) ) subjet1Idx = int( jIt - jetsColl2->begin() );
        if( &(*jIt) == fatJetToPrunedFatJetMap.find(&(*pjet))->second->daughter(1) ) subjet2Idx = int( jIt - jetsColl2->begin() );
        if( subjet1Idx>=0 && subjet2Idx>=0 ) break;
      }
      JetInfo[iJetColl].Jet_SubJet1Idx[JetInfo[iJetColl].nJet] = subjet1Idx;
      JetInfo[iJetColl].Jet_SubJet2Idx[JetInfo[iJetColl].nJet] = subjet2Idx;

      int nsubjettracks = 0, nsharedsubjettracks = 0;

      if( subjet1Idx>=0 && subjet2Idx>=0 ) // protection for pathological cases of pruned fat jets with only one constituent which results in an undefined subjet 2 index
      {
        for(int sj=0; sj<2; ++sj)
        {
          int subjetIdx = (sj==0 ? subjet1Idx : subjet2Idx); // subjet index
          int compSubjetIdx = (sj==0 ? subjet2Idx : subjet1Idx); // companion subjet index
          int nTracks = ( jetsColl2->at(subjetIdx).hasTagInfo("impactParameter") ? jetsColl2->at(subjetIdx).tagInfoTrackIP("impactParameter")->selectedTracks().size() : 0 );

          for(int t=0; t<nTracks; ++t)
          {
            if( reco::deltaR( jetsColl2->at(subjetIdx).tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->eta(), jetsColl2->at(subjetIdx).tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->phi(), jetsColl2->at(subjetIdx).eta(), jetsColl2->at(subjetIdx).phi() ) < 0.3 )
            {
              ++nsubjettracks;
              if( reco::deltaR( jetsColl2->at(subjetIdx).tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->eta(), jetsColl2->at(subjetIdx).tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->phi(), jetsColl2->at(compSubjetIdx).eta(), jetsColl2->at(compSubjetIdx).phi() ) < 0.3 )
              {
                if(sj==0) ++nsharedsubjettracks;
              }
            }
          }
        }
      }

      JetInfo[iJetColl].Jet_nsubjettracks[JetInfo[iJetColl].nJet] = nsubjettracks-nsharedsubjettracks;
      JetInfo[iJetColl].Jet_nsharedsubjettracks[JetInfo[iJetColl].nJet] = nsharedsubjettracks;
    }

    float etajet = TMath::Abs( pjet->eta() );
    float phijet = pjet->phi();
    if (phijet < 0.) phijet += 2*TMath::Pi();

    //*****************************************************************
    // Taggers
    //*****************************************************************

    // Loop on Selected Tracks

    const reco::TrackRefVector & selectedTracks( pjet->tagInfoTrackIP("impactParameter")->selectedTracks() );
    const reco::TrackRefVector & tracks( pjet->tagInfoTrackIP("impactParameter")->tracks() );

    int ntagtracks = 0;
    if (useSelectedTracks_) ntagtracks = selectedTracks.size();
    else ntagtracks = tracks.size();

    JetInfo[iJetColl].Jet_ntracks[JetInfo[iJetColl].nJet] = ntagtracks;

    JetInfo[iJetColl].Jet_nFirstTrack[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nTrack;
    JetInfo[iJetColl].Jet_nFirstTrkInc[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrkInc;
    int k=0;

    int nseltracks = 0;
    int nsharedtracks = 0;

    unsigned int trackSize = selectedTracks.size();
    if ( !useSelectedTracks_ ) trackSize = tracks.size();

    if ( allowJetSkipping_ && trackSize==0 ) continue;

    for (unsigned int itt=0; itt < trackSize; ++itt)
    {
      reco::Track  ptrack;
      reco::TrackRef  ptrackRef;
      if(useSelectedTracks_) 
      {
        ptrack = *(selectedTracks[itt]);
        ptrackRef = selectedTracks[itt];
      }
      else
      {
        ptrack = *(tracks[itt]);
        ptrackRef = tracks[itt];
      }

      TransientTrack transientTrack = trackBuilder->build(ptrack);
      GlobalVector direction(pjet->px(), pjet->py(), pjet->pz());

      //--------------------------------
      Double_t decayLength=-1;
      TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), *pv, direction, transientTrack.field());
      if (closest.isValid())
        decayLength =  (closest.globalPosition()-   RecoVertex::convertPos(pv->position())).mag();
      else
        decayLength = -1;

      Double_t distJetAxis =  IPTools::jetTrackDistance(transientTrack, direction, *pv).second.value();

      JetInfo[iJetColl].Track_dist[JetInfo[iJetColl].nTrack]     = distJetAxis;
      JetInfo[iJetColl].Track_length[JetInfo[iJetColl].nTrack]   = decayLength;

      JetInfo[iJetColl].Track_dxy[JetInfo[iJetColl].nTrack]      = ptrack.dxy(pv->position());
      JetInfo[iJetColl].Track_dz[JetInfo[iJetColl].nTrack]       = ptrack.dz(pv->position());
      JetInfo[iJetColl].Track_zIP[JetInfo[iJetColl].nTrack]      = ptrack.dz()-(*pv).z();

      float deta = ptrack.eta() - JetInfo[iJetColl].Jet_eta[JetInfo[iJetColl].nJet];
      float dphi = ptrack.phi() - JetInfo[iJetColl].Jet_phi[JetInfo[iJetColl].nJet];

      if ( dphi > TMath::Pi() ) dphi = 2.*TMath::Pi() - dphi;
      float deltaR = TMath::Sqrt(deta*deta + dphi*dphi);

      bool pass_cut_trk = false;
      if (std::fabs(distJetAxis) < 0.07 && decayLength < 5.0) pass_cut_trk = true;

      if (std::fabs(distJetAxis) < 0.07 && decayLength < 5.0
                                        && deltaR < 0.3) nseltracks++;

      if ( runSubJets_ && iJetColl == 1 && pass_cut_trk && subjet1Idx >= 0 && subjet2Idx >= 0 ) {

        deta = ptrack.eta() - JetInfo[0].Jet_eta[subjet1Idx];
        dphi = ptrack.phi() - JetInfo[0].Jet_phi[subjet1Idx];
        if ( dphi > TMath::Pi() ) dphi = 2.*TMath::Pi() - dphi;
        float dR1 = TMath::Sqrt(deta*deta + dphi*dphi);

        deta = ptrack.eta() - JetInfo[0].Jet_eta[subjet2Idx];
        dphi = ptrack.phi() - JetInfo[0].Jet_phi[subjet2Idx];
        if ( dphi > TMath::Pi() ) dphi = 2.*TMath::Pi() - dphi;
        float dR2 = TMath::Sqrt(deta*deta + dphi*dphi);

        if ( dR1 < 0.3 && dR2 < 0.3 ) nsharedtracks++;
      }

      // track selection
      if ( (useSelectedTracks_ && pass_cut_trk) ||  !useSelectedTracks_) {

        if ( useSelectedTracks_ ) {
          JetInfo[iJetColl].Track_IP2D[JetInfo[iJetColl].nTrack]     = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.value();
          JetInfo[iJetColl].Track_IP2Dsig[JetInfo[iJetColl].nTrack]  = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.significance();
          JetInfo[iJetColl].Track_IP[JetInfo[iJetColl].nTrack]       = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.value();
          JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack]    = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.significance();
          JetInfo[iJetColl].Track_IP2Derr[JetInfo[iJetColl].nTrack]  = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.error();
          JetInfo[iJetColl].Track_IPerr[JetInfo[iJetColl].nTrack]    = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.error();
          JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]    = pjet->tagInfoTrackIP("impactParameter")->probabilities(0)[k];
        }
        else {
          Measurement1D ip2d    = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pv).second;
          Measurement1D ip3d    = IPTools::signedImpactParameter3D(trackBuilder->build(ptrack), direction, *pv).second;
          //Measurement1D ip2dsig = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pv).first;
          //Measurement1D ip3dsig = IPTools::signedImpactParameter3D(trackBuilder->build(ptrack), direction, *pv).first;

          JetInfo[iJetColl].Track_IP2D[JetInfo[iJetColl].nTrack]     = (ip2d.value());
          JetInfo[iJetColl].Track_IP2Dsig[JetInfo[iJetColl].nTrack]  = (ip2d.value())/(ip2d.error());
          JetInfo[iJetColl].Track_IP[JetInfo[iJetColl].nTrack]       = (ip3d.value());
          JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack]    = (ip3d.value())/(ip3d.error());
          JetInfo[iJetColl].Track_IP2Derr[JetInfo[iJetColl].nTrack]  = (ip2d.error());
          JetInfo[iJetColl].Track_IPerr[JetInfo[iJetColl].nTrack]    = (ip3d.error());
          JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]    = -1000;
        }

        JetInfo[iJetColl].Track_p[JetInfo[iJetColl].nTrack]        = ptrack.p();
        JetInfo[iJetColl].Track_pt[JetInfo[iJetColl].nTrack]       = ptrack.pt();
        JetInfo[iJetColl].Track_eta[JetInfo[iJetColl].nTrack]      = ptrack.eta();
        JetInfo[iJetColl].Track_phi[JetInfo[iJetColl].nTrack]      = ptrack.phi();
        JetInfo[iJetColl].Track_chi2[JetInfo[iJetColl].nTrack]     = ptrack.normalizedChi2();
        JetInfo[iJetColl].Track_charge[JetInfo[iJetColl].nTrack]   = ptrack.charge();

        JetInfo[iJetColl].Track_nHitAll[JetInfo[iJetColl].nTrack]  = ptrack.numberOfValidHits();
        JetInfo[iJetColl].Track_nHitPixel[JetInfo[iJetColl].nTrack]= ptrack.hitPattern().numberOfValidPixelHits();
        JetInfo[iJetColl].Track_nHitStrip[JetInfo[iJetColl].nTrack]= ptrack.hitPattern().numberOfValidStripHits();
        JetInfo[iJetColl].Track_nHitTIB[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTIBHits();
        JetInfo[iJetColl].Track_nHitTID[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTIDHits();
        JetInfo[iJetColl].Track_nHitTOB[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTOBHits();
        JetInfo[iJetColl].Track_nHitTEC[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidStripTECHits();
        JetInfo[iJetColl].Track_nHitPXB[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidPixelBarrelHits();
        JetInfo[iJetColl].Track_nHitPXF[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().numberOfValidPixelEndcapHits();
        JetInfo[iJetColl].Track_isHitL1[JetInfo[iJetColl].nTrack]  = ptrack.hitPattern().hasValidHitInFirstPixelBarrel();

        setTracksPV(ptrackRef, primaryVertex,
                    JetInfo[iJetColl].Track_PV[JetInfo[iJetColl].nTrack],
                    JetInfo[iJetColl].Track_PVweight[JetInfo[iJetColl].nTrack]);

        if( pjet->hasTagInfo("secondaryVertex") )
        {
          setTracksSV(ptrackRef, pjet->tagInfoSecondaryVertex("secondaryVertex"),
                      JetInfo[iJetColl].Track_isfromSV[JetInfo[iJetColl].nTrack],
                      JetInfo[iJetColl].Track_SV[JetInfo[iJetColl].nTrack],
                      JetInfo[iJetColl].Track_SVweight[JetInfo[iJetColl].nTrack]);
        }
        else
        {
          JetInfo[iJetColl].Track_isfromSV[JetInfo[iJetColl].nTrack] = 0;
          JetInfo[iJetColl].Track_SV[JetInfo[iJetColl].nTrack] = -1;
          JetInfo[iJetColl].Track_SVweight[JetInfo[iJetColl].nTrack] = 0.;
        }

        JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] = 0;

        if ( useTrackHistory_ && !isData_ ) {
          TrackCategories::Flags theFlag ;
          if(useSelectedTracks_) theFlag  = classifier_.evaluate( pjet->tagInfoTrackIP("impactParameter")->selectedTracks()[k] ).flags();
          else                     theFlag  = classifier_.evaluate( pjet->tagInfoTrackIP("impactParameter")->tracks()[k] ).flags();

          if ( theFlag[TrackCategories::BWeakDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 1);
          if ( theFlag[TrackCategories::CWeakDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 2);
          if ( theFlag[TrackCategories::TauDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 3);
          if ( theFlag[TrackCategories::ConversionsProcess] )JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 4);
          if ( theFlag[TrackCategories::KsDecay] )	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 5);
          if ( theFlag[TrackCategories::LambdaDecay] )       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 6);
          if ( theFlag[TrackCategories::HadronicProcess] )   JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 7);
          if ( theFlag[TrackCategories::Fake] ) 	       JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 8);
          if ( theFlag[TrackCategories::SharedInnerHits] )   JetInfo[iJetColl].Track_history[JetInfo[iJetColl].nTrack] += pow(10, -1 + 9);
        }

        JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = -1;

        TLorentzVector track4P, jet4P;
        track4P.SetPtEtaPhiM(ptrack.pt(), ptrack.eta(), ptrack.phi(), 0. );
        jet4P.SetPtEtaPhiM( pjet->pt(), pjet->eta(), pjet->phi(), pjet->energy() );


        if ( jet4P.DeltaR(track4P) < 0.3 && std::fabs(distJetAxis) < 0.07 && decayLength < 5.0 ) {

          if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) Histos[iJetColl]->TrackProbaNeg->Fill(-JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]);
          if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 && PFJet80 ) Histos[iJetColl]->TrackProbJet80->Fill(-JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack]);

          if ( findCat( &ptrack, *&cat0 ) ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 0;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat0->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat0->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat0->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat0->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat1 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack]  = 1;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat1->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat1->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat1->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat1->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat2 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 2;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat2->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat2->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat2->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat2->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat3 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 3;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat3->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat3->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat3->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat3->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat4 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 3 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 4;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat4->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat4->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat4->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat4->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat5 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 3 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 4 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 5;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat5->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat5->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat5->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat5->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat6 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 3 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 4 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 5 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 6;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat6->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat6->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat6->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat6->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat7 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 3 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 4 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 5 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 6 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 7;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat7->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat7->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat7->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat7->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat8 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 3 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 4 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 5 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 6 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 7 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 8;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat8->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat8->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat8->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat8->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }

          if ( findCat( &ptrack, *&cat9 ) && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 0 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 1 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 2 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 3 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack]!= 4 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 5 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 6 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 7 && JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] != 8 ) {
            JetInfo[iJetColl].Track_category[JetInfo[iJetColl].nTrack] = 9;
            if ( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] < 0 ) {
              Histos[iJetColl]->IPSign_cat9->Fill( JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
              Histos[iJetColl]->TrackProbaNeg_Cat9->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              if ( PFJet80 ) {
                Histos[iJetColl]->IPSign_cat9->Fill( -JetInfo[iJetColl].Track_IPsig[JetInfo[iJetColl].nTrack] );
                Histos[iJetColl]->TrackProbJet80_Cat9->Fill( -JetInfo[iJetColl].Track_Proba[JetInfo[iJetColl].nTrack] );
              }
            }
          }
        }

        ++JetInfo[iJetColl].nTrack;
      }
      if ( useSelectedTracks_ ) JetInfo[iJetColl].Jet_ntracks[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nTrack-JetInfo[iJetColl].Jet_nFirstTrack[JetInfo[iJetColl].nJet];

      if ( producePtRelTemplate_ ) {
        if ( ptrack.quality(reco::TrackBase::highPurity)
            // Remove the tracks that are pixel-less
            && ptrack.algo()!=(reco::TrackBase::iter4)
            && ptrack.algo()!=(reco::TrackBase::iter5)
            && deltaR < 0.4
            && ptrack.pt() > 5.
            && ptrack.numberOfValidHits() >= 11
            && ptrack.hitPattern().numberOfValidPixelHits() >= 2
            && ptrack.normalizedChi2() < 10
            && ptrack.trackerExpectedHitsOuter().numberOfHits() <= 2
            && ptrack.dz()-(*pv).z() < 1. ) {

          if ( iJetColl != 1 ) {
            if ( useSelectedTracks_ ) {
              JetInfo[iJetColl].TrkInc_IP[JetInfo[iJetColl].nTrkInc]    = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.value();
              JetInfo[iJetColl].TrkInc_IPsig[JetInfo[iJetColl].nTrkInc] = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.significance();
            }
            JetInfo[iJetColl].TrkInc_pt[JetInfo[iJetColl].nTrkInc]    = ptrack.pt();
            JetInfo[iJetColl].TrkInc_eta[JetInfo[iJetColl].nTrkInc]   = ptrack.eta();
            JetInfo[iJetColl].TrkInc_phi[JetInfo[iJetColl].nTrkInc]   = ptrack.phi();
            JetInfo[iJetColl].TrkInc_ptrel[JetInfo[iJetColl].nTrkInc] = calculPtRel( ptrack , *pjet);

            ++JetInfo[iJetColl].nTrkInc;
          }
        }
      }
      ++k;

    } //// end loop on tracks

    JetInfo[iJetColl].Jet_nseltracks[JetInfo[iJetColl].nJet] = nseltracks;

    if ( runSubJets_ && iJetColl == 1 )
      JetInfo[iJetColl].Jet_nsharedtracks[JetInfo[iJetColl].nJet] = nsharedtracks;

    JetInfo[iJetColl].Jet_nLastTrack[JetInfo[iJetColl].nJet]   = JetInfo[iJetColl].nTrack;
    JetInfo[iJetColl].Jet_nLastTrkInc[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nTrkInc;

    // b-tagger discriminants
    float Proba  = pjet->bDiscriminator(jetPBJetTags_.c_str());
    float ProbaN = pjet->bDiscriminator(jetPNegBJetTags_.c_str());
    float ProbaP = pjet->bDiscriminator(jetPPosBJetTags_.c_str());

    float Bprob  = pjet->bDiscriminator(jetBPBJetTags_.c_str());
    float BprobN = pjet->bDiscriminator(jetBPNegBJetTags_.c_str());
    float BprobP = pjet->bDiscriminator(jetBPPosBJetTags_.c_str());

    float CombinedSvtx  = pjet->bDiscriminator(combinedSVBJetTags_.c_str());
    float CombinedSvtxN = pjet->bDiscriminator(combinedSVNegBJetTags_.c_str());
    float CombinedSvtxP = pjet->bDiscriminator(combinedSVPosBJetTags_.c_str());

    std::vector<const reco::BaseTagInfo*>  baseTagInfos_test;
    JetTagComputer::TagInfoHelper helper_test(baseTagInfos_test);
    baseTagInfos_test.push_back( pjet->tagInfoTrackIP("impactParameter") );
    baseTagInfos_test.push_back( pjet->tagInfoSecondaryVertex("secondaryVertex") );
    TaggingVariableList vars_test = computer->taggingVariables(helper_test);
    float JetVtxCategory = -9999.;
    if(vars_test.checkTag(reco::btau::vertexCategory)) JetVtxCategory = ( vars_test.get(reco::btau::vertexCategory) );

    float RetCombinedSvtx  = pjet->bDiscriminator(combinedSVRetrainedBJetTags_.c_str());
    float RetCombinedSvtxN = pjet->bDiscriminator(combinedSVRetrainedNegBJetTags_.c_str());
    float RetCombinedSvtxP = pjet->bDiscriminator(combinedSVRetrainedPosBJetTags_.c_str());

    float CombinedCSVJP  = pjet->bDiscriminator(combinedCSVJPBJetTags_.c_str());
    float CombinedCSVJPN = pjet->bDiscriminator(combinedCSVJPNegBJetTags_.c_str());
    float CombinedCSVJPP = pjet->bDiscriminator(combinedCSVJPPosBJetTags_.c_str());

    float CombinedCSVSL  = pjet->bDiscriminator(combinedCSVSLBJetTags_.c_str());
    float CombinedCSVSLN = pjet->bDiscriminator(combinedCSVSLNegBJetTags_.c_str());
    float CombinedCSVSLP = pjet->bDiscriminator(combinedCSVSLPosBJetTags_.c_str());

    float CombinedCSVJPSL  = pjet->bDiscriminator(combinedCSVJPSLBJetTags_.c_str());
    float CombinedCSVJPSLN = pjet->bDiscriminator(combinedCSVJPSLNegBJetTags_.c_str());
    float CombinedCSVJPSLP = pjet->bDiscriminator(combinedCSVJPSLPosBJetTags_.c_str());

    float SimpleIVF_HP    = pjet->bDiscriminator(simpleIVFSVHighPurBJetTags_.c_str());
    float SimpleIVF_HE    = pjet->bDiscriminator(simpleIVFSVHighEffBJetTags_.c_str());
    float DoubleIVF_HE    = pjet->bDiscriminator(doubleIVFSVHighEffBJetTags_.c_str());
    float CombinedIVF     = pjet->bDiscriminator(combinedIVFSVBJetTags_.c_str());
    float CombinedIVF_P   = pjet->bDiscriminator(combinedIVFSVPosBJetTags_.c_str());

    float Svtx    = pjet->bDiscriminator(simpleSVHighEffBJetTags_.c_str());
    float SvtxN   = pjet->bDiscriminator(simpleSVNegHighEffBJetTags_.c_str());
    float SvtxHP  = pjet->bDiscriminator(simpleSVHighPurBJetTags_.c_str());
    float SvtxNHP = pjet->bDiscriminator(simpleSVNegHighPurBJetTags_.c_str());

    float SoftM  = pjet->bDiscriminator(softPFMuonBJetTags_.c_str());
    float SoftMN = pjet->bDiscriminator(softPFMuonNegBJetTags_.c_str());
    float SoftMP = pjet->bDiscriminator(softPFMuonPosBJetTags_.c_str());

    // PFMuon information
    for (unsigned int leptIdx = 0; leptIdx < (pjet->hasTagInfo(softPFMuonTagInfos_.c_str()) ? pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->leptons() : 0); ++leptIdx) {

      JetInfo[iJetColl].PFMuon_IdxJet[JetInfo[iJetColl].nPFMuon]    = JetInfo[iJetColl].nJet;
      JetInfo[iJetColl].PFMuon_pt[JetInfo[iJetColl].nPFMuon]        = pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->lepton(leptIdx)->pt();
      JetInfo[iJetColl].PFMuon_eta[JetInfo[iJetColl].nPFMuon]       = pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->lepton(leptIdx)->eta();
      JetInfo[iJetColl].PFMuon_phi[JetInfo[iJetColl].nPFMuon]       = pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->lepton(leptIdx)->phi();
      JetInfo[iJetColl].PFMuon_ptrel[JetInfo[iJetColl].nPFMuon]     = calculPtRel( *(pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->lepton(leptIdx)), *pjet );
      JetInfo[iJetColl].PFMuon_ratio[JetInfo[iJetColl].nPFMuon]     = (pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->properties(leptIdx).ratio);
      JetInfo[iJetColl].PFMuon_ratioRel[JetInfo[iJetColl].nPFMuon]  = (pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->properties(leptIdx).ratioRel);
      JetInfo[iJetColl].PFMuon_deltaR[JetInfo[iJetColl].nPFMuon]    = (pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->properties(leptIdx).deltaR);
      JetInfo[iJetColl].PFMuon_IP[JetInfo[iJetColl].nPFMuon]        = (pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->properties(leptIdx).sip3d);
      JetInfo[iJetColl].PFMuon_IP2D[JetInfo[iJetColl].nPFMuon]      = (pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->properties(leptIdx).sip2d);

      JetInfo[iJetColl].PFMuon_GoodQuality[JetInfo[iJetColl].nPFMuon] = 0;
      int muIdx = matchMuon( pjet->tagInfoSoftLepton(softPFMuonTagInfos_.c_str())->lepton(leptIdx), muons );
      if ( muIdx != -1 && muons[muIdx].isGlobalMuon() == 1 ) {

        JetInfo[iJetColl].PFMuon_GoodQuality[JetInfo[iJetColl].nPFMuon] = 1;

        if (muons[muIdx].outerTrack()->hitPattern().numberOfValidMuonHits()>0 &&
            muons[muIdx].numberOfMatches()>1 && muons[muIdx].innerTrack()->hitPattern().numberOfValidHits()>10 &&
            muons[muIdx].innerTrack()->hitPattern().numberOfValidPixelHits()>1 &&
            muons[muIdx].innerTrack()->trackerExpectedHitsOuter().numberOfHits()<3 &&
            muons[muIdx].globalTrack()->normalizedChi2()<10. && muons[muIdx].innerTrack()->normalizedChi2()<10.)
          JetInfo[iJetColl].PFMuon_GoodQuality[JetInfo[iJetColl].nPFMuon] = 2;

      }
      ++JetInfo[iJetColl].nPFMuon;
    }

    // Old soft muon information for perfomance study safety
    for (unsigned int leptIdx = 0; leptIdx < pjet->tagInfoSoftLepton("softMuon")->leptons(); ++leptIdx){

      int muIdx = matchMuon( pjet->tagInfoSoftLepton("softMuon")->lepton(leptIdx), muons );
      if ( muIdx != -1 && muons[muIdx].isGlobalMuon() ) {
        JetInfo[iJetColl].Muon_IdxJet[JetInfo[iJetColl].nMuon]   = JetInfo[iJetColl].nJet;
        JetInfo[iJetColl].Muon_ptrel[JetInfo[iJetColl].nMuon]    = calculPtRel( *(pjet->tagInfoSoftLepton("softMuon")->lepton(leptIdx)), *pjet );
        JetInfo[iJetColl].Muon_nTkHit[JetInfo[iJetColl].nMuon]   = muons[muIdx].innerTrack()->hitPattern().numberOfValidHits();
        JetInfo[iJetColl].Muon_nPixHit[JetInfo[iJetColl].nMuon]  = muons[muIdx].innerTrack()->hitPattern().numberOfValidPixelHits();
        JetInfo[iJetColl].Muon_nOutHit[JetInfo[iJetColl].nMuon]  = muons[muIdx].innerTrack()->trackerExpectedHitsOuter().numberOfHits();
        JetInfo[iJetColl].Muon_nMuHit[JetInfo[iJetColl].nMuon]   = muons[muIdx].outerTrack()->hitPattern().numberOfValidMuonHits();
        JetInfo[iJetColl].Muon_chi2[JetInfo[iJetColl].nMuon]     = muons[muIdx].globalTrack()->normalizedChi2();
        JetInfo[iJetColl].Muon_chi2Tk[JetInfo[iJetColl].nMuon]   = muons[muIdx].innerTrack()->normalizedChi2();
        JetInfo[iJetColl].Muon_pt[JetInfo[iJetColl].nMuon]       = muons[muIdx].pt();
        JetInfo[iJetColl].Muon_eta[JetInfo[iJetColl].nMuon]      = muons[muIdx].eta();
        JetInfo[iJetColl].Muon_phi[JetInfo[iJetColl].nMuon]      = muons[muIdx].phi();
        JetInfo[iJetColl].Muon_ratio[JetInfo[iJetColl].nMuon]    = (pjet->tagInfoSoftLepton("softMuon")->properties(leptIdx).ratio);
        JetInfo[iJetColl].Muon_ratioRel[JetInfo[iJetColl].nMuon] = (pjet->tagInfoSoftLepton("softMuon")->properties(leptIdx).ratioRel);
        JetInfo[iJetColl].Muon_deltaR[JetInfo[iJetColl].nMuon]   = (pjet->tagInfoSoftLepton("softMuon")->properties(leptIdx).deltaR);

        JetInfo[iJetColl].Muon_isGlobal[JetInfo[iJetColl].nMuon] = 1;
        JetInfo[iJetColl].Muon_nMatched[JetInfo[iJetColl].nMuon] = muons[muIdx].numberOfMatches() ;
        JetInfo[iJetColl].Muon_vz[JetInfo[iJetColl].nMuon]       = muons[muIdx].vz();
        int mutkid = getMuonTk(muons[muIdx].innerTrack()->pt(), iJetColl);
        JetInfo[iJetColl].Muon_TrackIdx[JetInfo[iJetColl].nMuon] = ( mutkid>=0 ? mutkid: -1 );

        JetInfo[iJetColl].Muon_Proba[JetInfo[iJetColl].nMuon] = -100.;
        if ( mutkid >= 0 ) {
          JetInfo[iJetColl].Muon_Proba[JetInfo[iJetColl].nMuon] = JetInfo[iJetColl].Track_Proba[mutkid];
        }

        JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] = 0;
        if ( useTrackHistory_ && !isData_ ) {
          TrackCategories::Flags theFlagP = classifier_.evaluate( pjet->tagInfoSoftLepton("softMuon")->lepton(leptIdx) ).flags();
          if ( theFlagP[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 1));
          if ( theFlagP[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 2));
          if ( theFlagP[TrackCategories::TauDecay] )           JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 3));
          if ( theFlagP[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 4));
          if ( theFlagP[TrackCategories::KsDecay] )            JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 5));
          if ( theFlagP[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 6));
          if ( theFlagP[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 7));
          if ( theFlagP[TrackCategories::Fake] )               JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 8));
          if ( theFlagP[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Muon_hist[JetInfo[iJetColl].nMuon] += int(pow(10., -1 + 9));
        }

        //---------------------------------
        // calculate IP/s of muons' tracks
        //---------------------------------
        reco::TransientTrack tt = trackBuilder->build(muons[muIdx].innerTrack());
        GlobalVector direction(pjet->px(), pjet->py(), pjet->pz());
        Measurement1D ip   = IPTools::signedImpactParameter3D(tt, direction, (*primaryVertex)[0]).second;
        Measurement1D ip2d = IPTools::signedTransverseImpactParameter(tt, direction, (*primaryVertex)[0]).second;

        JetInfo[iJetColl].Muon_IP[JetInfo[iJetColl].nMuon]   = ip.value();
        JetInfo[iJetColl].Muon_IP2D[JetInfo[iJetColl].nMuon] = ip2d.value();
        JetInfo[iJetColl].Muon_IPsig[JetInfo[iJetColl].nMuon]   = (ip.value())/(ip.error());
        JetInfo[iJetColl].Muon_IP2Dsig[JetInfo[iJetColl].nMuon] = (ip2d.value())/(ip2d.error());

        ++JetInfo[iJetColl].nMuon;
      }
    }

    float SoftE  = pjet->bDiscriminator(softPFElectronBJetTags_.c_str());
    float SoftEN = pjet->bDiscriminator(softPFElectronNegBJetTags_.c_str());
    float SoftEP = pjet->bDiscriminator(softPFElectronPosBJetTags_.c_str());

    // PFElectron information
    for (unsigned int leptIdx = 0; leptIdx < (pjet->hasTagInfo(softPFElectronTagInfos_.c_str()) ? pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->leptons() : 0); ++leptIdx) {

      JetInfo[iJetColl].PFElectron_IdxJet[JetInfo[iJetColl].nPFElectron]    = JetInfo[iJetColl].nJet;
      JetInfo[iJetColl].PFElectron_pt[JetInfo[iJetColl].nPFElectron]        = pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->lepton(leptIdx)->pt();
      JetInfo[iJetColl].PFElectron_eta[JetInfo[iJetColl].nPFElectron]       = pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->lepton(leptIdx)->eta();
      JetInfo[iJetColl].PFElectron_phi[JetInfo[iJetColl].nPFElectron]       = pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->lepton(leptIdx)->phi();
      JetInfo[iJetColl].PFElectron_ptrel[JetInfo[iJetColl].nPFElectron]     = calculPtRel( *(pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->lepton(leptIdx)), *pjet );
      JetInfo[iJetColl].PFElectron_ratio[JetInfo[iJetColl].nPFElectron]     = (pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->properties(leptIdx).ratio);
      JetInfo[iJetColl].PFElectron_ratioRel[JetInfo[iJetColl].nPFElectron]  = (pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->properties(leptIdx).ratioRel);
      JetInfo[iJetColl].PFElectron_deltaR[JetInfo[iJetColl].nPFElectron]    = (pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->properties(leptIdx).deltaR);
      JetInfo[iJetColl].PFElectron_IP[JetInfo[iJetColl].nPFElectron]        = (pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->properties(leptIdx).sip3d);
      JetInfo[iJetColl].PFElectron_IP2D[JetInfo[iJetColl].nPFElectron]      = (pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->properties(leptIdx).sip2d);
      JetInfo[iJetColl].PFElectron_mva_e_pi[JetInfo[iJetColl].nPFElectron]  = -999.;//pjet->tagInfoSoftLepton(softPFElectronTagInfos_.c_str())->properties(leptIdx).mva_e_pi;

      ++JetInfo[iJetColl].nPFElectron;
    }

    // Jet information
    JetInfo[iJetColl].Jet_ProbaN[JetInfo[iJetColl].nJet]   = ProbaN;
    JetInfo[iJetColl].Jet_ProbaP[JetInfo[iJetColl].nJet]   = ProbaP;
    JetInfo[iJetColl].Jet_Proba[JetInfo[iJetColl].nJet]    = Proba;
    JetInfo[iJetColl].Jet_BprobN[JetInfo[iJetColl].nJet]   = BprobN;
    JetInfo[iJetColl].Jet_BprobP[JetInfo[iJetColl].nJet]   = BprobP;
    JetInfo[iJetColl].Jet_Bprob[JetInfo[iJetColl].nJet]    = Bprob;
    JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet]     = SvtxN;
    JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet]      = Svtx;
    JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet]   = SvtxNHP;
    JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet]    = SvtxHP;
    JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] = CombinedSvtxN;
    JetInfo[iJetColl].Jet_CombSvxP[JetInfo[iJetColl].nJet] = CombinedSvtxP;
    JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet]  = CombinedSvtx;
    JetInfo[iJetColl].Jet_RetCombSvxN[JetInfo[iJetColl].nJet] = RetCombinedSvtxN;
    JetInfo[iJetColl].Jet_RetCombSvxP[JetInfo[iJetColl].nJet] = RetCombinedSvtxP;
    JetInfo[iJetColl].Jet_RetCombSvx[JetInfo[iJetColl].nJet]  = RetCombinedSvtx;
    JetInfo[iJetColl].Jet_CombCSVJP_N[JetInfo[iJetColl].nJet] = CombinedCSVJPN;
    JetInfo[iJetColl].Jet_CombCSVJP_P[JetInfo[iJetColl].nJet] = CombinedCSVJPP;
    JetInfo[iJetColl].Jet_CombCSVJP[JetInfo[iJetColl].nJet]  = CombinedCSVJP;
    JetInfo[iJetColl].Jet_CombCSVSL_N[JetInfo[iJetColl].nJet] = CombinedCSVSLN;
    JetInfo[iJetColl].Jet_CombCSVSL_P[JetInfo[iJetColl].nJet] = CombinedCSVSLP;
    JetInfo[iJetColl].Jet_CombCSVSL[JetInfo[iJetColl].nJet]  = CombinedCSVSL;
    JetInfo[iJetColl].Jet_CombCSVJPSL_N[JetInfo[iJetColl].nJet] = CombinedCSVJPSLN;
    JetInfo[iJetColl].Jet_CombCSVJPSL_P[JetInfo[iJetColl].nJet] = CombinedCSVJPSLP;
    JetInfo[iJetColl].Jet_CombCSVJPSL[JetInfo[iJetColl].nJet]  = CombinedCSVJPSL;
    JetInfo[iJetColl].Jet_SimpIVF_HP[JetInfo[iJetColl].nJet]  = SimpleIVF_HP;
    JetInfo[iJetColl].Jet_SimpIVF_HE[JetInfo[iJetColl].nJet]  = SimpleIVF_HE;
    JetInfo[iJetColl].Jet_DoubIVF_HE[JetInfo[iJetColl].nJet]  = DoubleIVF_HE;
    JetInfo[iJetColl].Jet_CombIVF[JetInfo[iJetColl].nJet]     = CombinedIVF;
    JetInfo[iJetColl].Jet_CombIVF_P[JetInfo[iJetColl].nJet]   = CombinedIVF_P;
    JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet]  = SoftMN;
    JetInfo[iJetColl].Jet_SoftMuP[JetInfo[iJetColl].nJet]  = SoftMP;
    JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet]   = SoftM;
    JetInfo[iJetColl].Jet_SoftElN[JetInfo[iJetColl].nJet]  = SoftEN;
    JetInfo[iJetColl].Jet_SoftElP[JetInfo[iJetColl].nJet]  = SoftEP;
    JetInfo[iJetColl].Jet_SoftEl[JetInfo[iJetColl].nJet]   = SoftE;

    JetInfo[iJetColl].Jet_VtxCat[JetInfo[iJetColl].nJet] = JetVtxCategory;

    std::vector<TrackIPTagInfo::TrackIPData>  ipdata = pjet->tagInfoTrackIP("impactParameter")->impactParameterData();
    std::vector<std::size_t> indexes( sortedIndexes(ipdata) );

    TrackCategories::Flags flags1P;
    TrackCategories::Flags flags2P;
    TrackCategories::Flags flags3P;
    TrackCategories::Flags flags1N;
    TrackCategories::Flags flags2N;
    TrackCategories::Flags flags3N;
    int idSize = 0;

    if ( useTrackHistory_ && indexes.size() != 0 && !isData_ ) {

      idSize = indexes.size();
      flags1P = classifier_.evaluate( selectedTracks[indexes[0]] ).flags();
      if(idSize > 1) flags2P = classifier_.evaluate( selectedTracks[indexes[1]] ).flags();
      if(idSize > 2) flags3P = classifier_.evaluate( selectedTracks[indexes[2]] ).flags();
      flags1N = classifier_.evaluate( selectedTracks[indexes[idSize-1]] ).flags();
      if(idSize > 1) flags2N = classifier_.evaluate( selectedTracks[indexes[idSize-2]] ).flags();
      if(idSize > 2) flags3N = classifier_.evaluate( selectedTracks[indexes[idSize-3]] ).flags();
    }

    JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCHEBJetTags_.c_str());
    JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCNegHEBJetTags_.c_str());

    JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCHPBJetTags_.c_str());
    JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet]   = pjet->bDiscriminator(trackCNegHPBJetTags_.c_str());

    //*****************************************************************
    //get track histories for 1st, 2nd and 3rd track (TC)
    //*****************************************************************
    JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] = 0;

    // Track Histosry
    if ( useTrackHistory_ && indexes.size()!=0 && !isData_ ) {
      if ( flags1P[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 1));
      if ( flags1P[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 2));
      if ( flags1P[TrackCategories::TauDecay] )           JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 3));
      if ( flags1P[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 4));
      if ( flags1P[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 5));
      if ( flags1P[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 6));
      if ( flags1P[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 7));
      if ( flags1P[TrackCategories::Fake] )	          JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 8));
      if ( flags1P[TrackCategories::SharedInnerHits] )	  JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 9));
      if ( idSize > 1 ) {
        if ( flags2P[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 1));
        if ( flags2P[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 2));
        if ( flags2P[TrackCategories::TauDecay] )           JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 3));
        if ( flags2P[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 4));
        if ( flags2P[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 5));
        if ( flags2P[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 6));
        if ( flags2P[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 7));
        if ( flags2P[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 8));
        if ( flags2P[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 9));
      }
      if ( idSize > 2 ) {
        if ( flags3P[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 1));
        if ( flags3P[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 2));
        if ( flags3P[TrackCategories::TauDecay] )           JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 3));
        if ( flags3P[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 4));
        if ( flags3P[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 5));
        if ( flags3P[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 6));
        if ( flags3P[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 7));
        if ( flags3P[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 8));
        if ( flags3P[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += int(pow(10., -1 + 9));
      }
      if ( flags1N[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 1));
      if ( flags1N[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 2));
      if ( flags1N[TrackCategories::TauDecay] )           JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 3));
      if ( flags1N[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 4));
      if ( flags1N[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 5));
      if ( flags1N[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 6));
      if ( flags1N[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 7));
      if ( flags1N[TrackCategories::Fake] )	          JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 8));
      if ( flags1N[TrackCategories::SharedInnerHits] )	  JetInfo[iJetColl].Jet_hist1[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 9));
      if ( idSize > 1 ) {
        if ( flags2N[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 1));
        if ( flags2N[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 2));
        if ( flags2N[TrackCategories::TauDecay] )           JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 3));
        if ( flags2N[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 4));
        if ( flags2N[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 5));
        if ( flags2N[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 6));
        if ( flags2N[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 7));
        if ( flags2N[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 8));
        if ( flags2N[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist2[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 9));
      }
      if ( idSize > 2 ) {
        if ( flags3N[TrackCategories::BWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 1));
        if ( flags3N[TrackCategories::CWeakDecay] )         JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 2));
        if ( flags3N[TrackCategories::TauDecay] )           JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 3));
        if ( flags3N[TrackCategories::ConversionsProcess] ) JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 4));
        if ( flags3N[TrackCategories::KsDecay] )            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 5));
        if ( flags3N[TrackCategories::LambdaDecay] )        JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 6));
        if ( flags3N[TrackCategories::HadronicProcess] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 7));
        if ( flags3N[TrackCategories::Fake] )	            JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 8));
        if ( flags3N[TrackCategories::SharedInnerHits] )    JetInfo[iJetColl].Jet_hist3[JetInfo[iJetColl].nJet] += 2*int(pow(10., -1 + 9));
      }
    }

    //*****************************************************************
    //get track Histosries of tracks in jets (for Jet Proba)
    //*****************************************************************
    JetInfo[iJetColl].Jet_histJet[JetInfo[iJetColl].nJet] = 0;

    if(useTrackHistory_ && !isData_){
      TrackRefVector jetProbTracks( pjet->tagInfoTrackIP("impactParameter")->selectedTracks() );

      cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0;
      can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0;

      for (unsigned int i=0; i<jetProbTracks.size(); ++i) {
        reco::TrackIPTagInfo::TrackIPData ip = (pjet->tagInfoTrackIP("impactParameter")->impactParameterData())[i];

        if ( ip.ip3d.significance() > 0 ) {
          TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
          if ( theFlag[TrackCategories::BWeakDecay] )	      cap0 = 1;
          if ( theFlag[TrackCategories::CWeakDecay] )	      cap1 = 1;
          if ( theFlag[TrackCategories::TauDecay] )	      cap2 = 1;
          if ( theFlag[TrackCategories::ConversionsProcess] ) cap3 = 1;
          if ( theFlag[TrackCategories::KsDecay] )	      cap4 = 1;
          if ( theFlag[TrackCategories::LambdaDecay] )        cap5 = 1;
          if ( theFlag[TrackCategories::HadronicProcess] )    cap6 = 1;
          if ( theFlag[TrackCategories::Fake] ) 	      cap7 = 1;
          if ( theFlag[TrackCategories::SharedInnerHits] )    cap8 = 1;
        }
        else {
          TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
          if ( theFlag[TrackCategories::BWeakDecay] )	      can0 = 2;
          if ( theFlag[TrackCategories::CWeakDecay] )	      can1 = 2;
          if ( theFlag[TrackCategories::TauDecay] )	      can2 = 2;
          if ( theFlag[TrackCategories::ConversionsProcess] ) can3 = 2;
          if ( theFlag[TrackCategories::KsDecay] )	      can4 = 2;
          if ( theFlag[TrackCategories::LambdaDecay] )        can5 = 2;
          if ( theFlag[TrackCategories::HadronicProcess] )    can6 = 2;
          if ( theFlag[TrackCategories::Fake] ) 	      can7 = 2;
          if ( theFlag[TrackCategories::SharedInnerHits] )    can8 = 2;
        }
      }
      JetInfo[iJetColl].Jet_histJet[JetInfo[iJetColl].nJet] =  cap0+can0 + (cap1+can1)*10 + (cap2+can2)*100
        + (cap3+can3)*1000     + (cap4+can4)*10000
        + (cap5+can5)*100000   + (cap6+can6)*1000000
        + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }

    //*****************************************************************
    //get track Histories  associated to sec. vertex (for simple SV)
    //*****************************************************************
    JetInfo[iJetColl].Jet_histSvx[JetInfo[iJetColl].nJet] = 0;
    JetInfo[iJetColl].Jet_nFirstSV[JetInfo[iJetColl].nJet]  = JetInfo[iJetColl].nSV;
    JetInfo[iJetColl].Jet_SV_multi[JetInfo[iJetColl].nJet]  = pjet->tagInfoSecondaryVertex("secondaryVertex")->nVertices();

    std::vector<const reco::BaseTagInfo*>  baseTagInfos;
    JetTagComputer::TagInfoHelper helper(baseTagInfos);
    baseTagInfos.push_back( pjet->tagInfoTrackIP("impactParameter") );
    baseTagInfos.push_back( pjet->tagInfoSecondaryVertex("secondaryVertex") );

    TaggingVariableList vars = computer->taggingVariables(helper);

    float SVmass = 0.;
    if ( JetInfo[iJetColl].Jet_SV_multi[JetInfo[iJetColl].nJet] > 0) {
      SVmass = ( vars.checkTag(reco::btau::vertexMass) ? vars.get(reco::btau::vertexMass) : -9999 );
    }
    JetInfo[iJetColl].Jet_SvxMass[JetInfo[iJetColl].nJet] = SVmass;

    // if secondary vertices present
    if ( JetInfo[iJetColl].Jet_SV_multi[JetInfo[iJetColl].nJet] > 0 )
    {

      JetInfo[iJetColl].SV_x[JetInfo[iJetColl].nSV]    = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).x();
      JetInfo[iJetColl].SV_y[JetInfo[iJetColl].nSV]    = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).y();
      JetInfo[iJetColl].SV_z[JetInfo[iJetColl].nSV]    = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).z();
      JetInfo[iJetColl].SV_ex[JetInfo[iJetColl].nSV]   = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).xError();
      JetInfo[iJetColl].SV_ey[JetInfo[iJetColl].nSV]   = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).yError();
      JetInfo[iJetColl].SV_ez[JetInfo[iJetColl].nSV]   = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).zError();
      JetInfo[iJetColl].SV_chi2[JetInfo[iJetColl].nSV] = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).chi2();
      JetInfo[iJetColl].SV_ndf[JetInfo[iJetColl].nSV]  = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).ndof();

      //cout << "pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).ndof() " << pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).ndof() << endl;

      JetInfo[iJetColl].SV_flight[JetInfo[iJetColl].nSV]	  = pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDistance(0).value();
      JetInfo[iJetColl].SV_flightErr[JetInfo[iJetColl].nSV]   = pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDistance(0).error();
      JetInfo[iJetColl].SV_flight2D[JetInfo[iJetColl].nSV]	  = pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDistance(0, true).value();
      JetInfo[iJetColl].SV_flight2DErr[JetInfo[iJetColl].nSV] = pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDistance(0, true).error();
      JetInfo[iJetColl].SV_nTrk[JetInfo[iJetColl].nSV]        = pjet->tagInfoSecondaryVertex("secondaryVertex")->nVertexTracks();
      JetInfo[iJetColl].SV_nTrk_firstVxt[JetInfo[iJetColl].nSV] =  pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0).tracksSize();

      // ------------------------added by Camille ---------------------------------------------------------------------------------------//

      edm::RefToBase<Jet> jet = pjet->tagInfoTrackIP("impactParameter")->jet();

      if(vars.checkTag(reco::btau::vertexEnergyRatio)) JetInfo[iJetColl].SV_energy_ratio[JetInfo[iJetColl].nSV] = ( vars.get(reco::btau::vertexEnergyRatio) );
      else JetInfo[iJetColl].SV_energy_ratio[JetInfo[iJetColl].nSV] = ( -9999 );

      if(vars.checkTag(reco::btau::vertexJetDeltaR)) JetInfo[iJetColl].SV_deltaR_jet[JetInfo[iJetColl].nSV] = (  vars.get(reco::btau::vertexJetDeltaR) );
      else JetInfo[iJetColl].SV_deltaR_jet[JetInfo[iJetColl].nSV] = ( -9999 );

      if(vars.checkTag(reco::btau::trackSip3dSigAboveCharm) ) JetInfo[iJetColl].SV_aboveC[JetInfo[iJetColl].nSV] = (  vars.get(reco::btau::trackSip3dSigAboveCharm ));
      else JetInfo[iJetColl].SV_aboveC[JetInfo[iJetColl].nSV] = (  -9999 );

      if(vars.checkTag(reco::btau::vertexMass)) JetInfo[iJetColl].SV_mass[JetInfo[iJetColl].nSV] = ( vars.get(reco::btau::vertexMass));
      else  JetInfo[iJetColl].SV_mass[JetInfo[iJetColl].nSV] = ( -9999 );

      Int_t totcharge=0;
      TrackKinematics vertexKinematics;

      const Vertex &vertex = pjet->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(0);

      Bool_t hasRefittedTracks = vertex.hasRefittedTracks();

      TrackRefVector vertexTracks = pjet->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks(0);
      for(TrackRefVector::const_iterator track = vertexTracks.begin();
          track != vertexTracks.end(); ++track) {
        Double_t w = pjet->tagInfoSecondaryVertex("secondaryVertex")->trackWeight(0, *track);
        if (w < 0.5)
          continue;
        if (hasRefittedTracks) {
          Track actualTrack = vertex.refittedTrack(*track);
          vertexKinematics.add(actualTrack, w);
          totcharge+=actualTrack.charge();
        }
        else {
          vertexKinematics.add(**track, w);
          const reco::Track& mytrack = **track;
          totcharge+=mytrack.charge();
        }
      }

      math::XYZVector jetDir = jet->momentum().Unit();

      Bool_t useTrackWeights = true;
      math::XYZTLorentzVector vertexSum = useTrackWeights
        ? vertexKinematics.weightedVectorSum()
        : vertexKinematics.vectorSum();

      math::XYZTLorentzVector flightDir( pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDirection(0).x(), pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDirection(0).y(), pjet->tagInfoSecondaryVertex("secondaryVertex")->flightDirection(0).z(), 0  );
      JetInfo[iJetColl].SV_deltaR_sum_jet[JetInfo[iJetColl].nSV] = ( Geom::deltaR(vertexSum, jetDir) );
      JetInfo[iJetColl].SV_deltaR_sum_dir[JetInfo[iJetColl].nSV] = ( Geom::deltaR( flightDir, vertexSum ) );
      JetInfo[iJetColl].SV_vtx_pt[JetInfo[iJetColl].nSV] = vertex.p4().pt();
      JetInfo[iJetColl].SV_vtx_eta[JetInfo[iJetColl].nSV] = vertex.p4().eta();
      JetInfo[iJetColl].SV_vtx_phi[JetInfo[iJetColl].nSV] = vertex.p4().phi();

      math::XYZVector jetVertex = (math::XYZVector(jet->vx(),jet->vy(),jet->vz()));

      Line::PositionType pos(GlobalPoint(vertex.p4().x(),vertex.p4().y(),vertex.p4().z()));
      Line::DirectionType dir(GlobalVector(flightDir.px(),flightDir.py(),flightDir.pz()));
      Line trackline(pos,dir);
      // get the Jet  line
      Line::PositionType pos2(GlobalPoint(jetVertex.x(),jetVertex.y(),jetVertex.z()));
      Line::DirectionType dir2(GlobalVector(jetDir.x(),jetDir.y(),jetDir.z()));
      Line jetline(pos2,dir2);
      // now compute the distance between the two lines
      JetInfo[iJetColl].SV_vtxDistJetAxis[JetInfo[iJetColl].nSV] = (jetline.distance(trackline)).mag();

      // total charge at the secondary vertex
      JetInfo[iJetColl].SV_totCharge[JetInfo[iJetColl].nSV]=totcharge;
      // ------------------------end added ---------------------------------------------------------------------------------------//

      ++JetInfo[iJetColl].nSV;

    } //// if secondary vertices present
    JetInfo[iJetColl].Jet_nLastSV[JetInfo[iJetColl].nJet] = JetInfo[iJetColl].nSV;


    cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0;
    can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0;

    TrackRefVector svxPostracks( pjet->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks(0) );

    //*****************************************************************
    // for Mistag studies
    //*****************************************************************
    if ( useTrackHistory_ && !isData_ ) {
      for (unsigned int i=0; i<svxPostracks.size(); ++i) {
        TrackCategories::Flags theFlag = classifier_.evaluate( svxPostracks[i] ).flags();
        if ( theFlag[TrackCategories::BWeakDecay] )         cap0 = 1;
        if ( theFlag[TrackCategories::CWeakDecay] )         cap1 = 1;
        if ( theFlag[TrackCategories::TauDecay] )           cap2 = 1;
        if ( theFlag[TrackCategories::ConversionsProcess] ) cap3 = 1;
        if ( theFlag[TrackCategories::KsDecay] )            cap4 = 1;
        if ( theFlag[TrackCategories::LambdaDecay] )        cap5 = 1;
        if ( theFlag[TrackCategories::HadronicProcess] )    cap6 = 1;
        if ( theFlag[TrackCategories::Fake] )	            cap7 = 1;
        if ( theFlag[TrackCategories::SharedInnerHits] )    cap8 = 1;
      }
    }

    TrackRefVector svxNegtracks( pjet->tagInfoSecondaryVertex("secondaryVertexNegative")->vertexTracks(0) );

    if ( useTrackHistory_ && !isData_ ) {
      for (unsigned int i=0; i<svxNegtracks.size(); ++i) {
        TrackCategories::Flags theFlag = classifier_.evaluate( svxNegtracks[i] ).flags();
        if ( theFlag[TrackCategories::BWeakDecay] )         can0 = 2;
        if ( theFlag[TrackCategories::CWeakDecay] )         can1 = 2;
        if ( theFlag[TrackCategories::TauDecay] )           can2 = 2;
        if ( theFlag[TrackCategories::ConversionsProcess] ) can3 = 2;
        if ( theFlag[TrackCategories::KsDecay] )            can4 = 2;
        if ( theFlag[TrackCategories::LambdaDecay] )        can5 = 2;
        if ( theFlag[TrackCategories::HadronicProcess] )    can6 = 2;
        if ( theFlag[TrackCategories::Fake] )	            can7 = 2;
        if ( theFlag[TrackCategories::SharedInnerHits] )    can8 = 2;
      }
      JetInfo[iJetColl].Jet_histSvx[JetInfo[iJetColl].nJet] =  cap0+can0 + (cap1+can1)*10 + (cap2+can2)*100
        + (cap3+can3)*1000     + (cap4+can4)*10000
        + (cap5+can5)*100000   + (cap6+can6)*1000000
        + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }

    //*********************************
    // Jet analysis

    Histos[iJetColl]->hData_All_NTracks->Fill( ntagtracks ) ;
    Histos[iJetColl]->hData_All_JetPt->Fill( ptjet );
    Histos[iJetColl]->hData_All_JetEta->Fill( etajet );

    //*****************************************************************
    //define positive and negative tags
    //*****************************************************************
    float varpos = -1000.;
    float varneg = -1000.;

    if ( selTagger_ == 0 ) {       // jet proba
      if ( ProbaP > 0 ) varpos = 20.*ProbaP;
      if ( ProbaN > 0 ) varneg = 20.*ProbaN;
    }
    else if ( selTagger_ == 2 ) {     // TC High Eff.
      if ( JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet] > 0. ) varpos = JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet];
      if ( JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet] > 0. ) varneg = JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 3 ) {     // TC High Pure.
      if ( JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet] > 0. ) varpos = JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet];
      if ( JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet] > 0. ) varneg = JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 4 ) {    // SSVHE
      if  (JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet] > 1)   varpos =  10.*JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet] - 10.;
      if  (JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] < -1) varneg = -10.*JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] - 10.;
    }
    else if ( selTagger_ == 5 ) {    // SV combined
      if  (JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet] > 0)  varpos = 50.*JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet];
      if  (JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] > 0) varneg = 50.*JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 6 ) {    // soft muon ptrel
      if  ( JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet]  > 0) varpos =  5*JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet];
      if  ( JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet] < 0) varneg = -5*JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet];
    }
    else if ( selTagger_ == 7 ) {    // SSVHP
      if  (JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet] > 1)   varpos =  10.*JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet] - 10.;
      if  (JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] < -1) varneg = -10.*JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] - 10.;
    }

    Histos[iJetColl]->hData_NTracks->Fill( ntagtracks ) ;
    Histos[iJetColl]->hData_JetPt->Fill( ptjet );
    Histos[iJetColl]->hData_JetEta->Fill( etajet );

    //*********************************
    // Tagging

    if ( varneg > 0 ) Histos[iJetColl]->hData_Tagger->Fill(-varneg );
    if ( varpos > 0 ) Histos[iJetColl]->hData_Tagger->Fill( varpos );

    if ( JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHE->Fill( JetInfo[iJetColl].Jet_Ip2P[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHE->Fill(-JetInfo[iJetColl].Jet_Ip2N[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHP->Fill( JetInfo[iJetColl].Jet_Ip3P[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->    hData_Tagger_TCHP->Fill(-JetInfo[iJetColl].Jet_Ip3N[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_ProbaP[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->  hData_Tagger_JP->Fill( 10*JetInfo[iJetColl].Jet_ProbaP[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_ProbaN[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->  hData_Tagger_JP->Fill(-10*JetInfo[iJetColl].Jet_ProbaN[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet]  >  1 ) Histos[iJetColl]->   hData_Tagger_SSVHE->Fill( 5*JetInfo[iJetColl].Jet_Svx[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] < -1 ) Histos[iJetColl]->   hData_Tagger_SSVHE->Fill( 5*JetInfo[iJetColl].Jet_SvxN[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet]  >  1 ) Histos[iJetColl]-> hData_Tagger_SSVHP->Fill( 5*JetInfo[iJetColl].Jet_SvxHP[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] < -1 ) Histos[iJetColl]-> hData_Tagger_SSVHP->Fill( 5*JetInfo[iJetColl].Jet_SvxNHP[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet] > 0  ) Histos[iJetColl]->hData_Tagger_CSV->Fill( 25*JetInfo[iJetColl].Jet_CombSvx[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] > 0 ) Histos[iJetColl]->hData_Tagger_CSV->Fill(-25*JetInfo[iJetColl].Jet_CombSvxN[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet]  > 0 ) Histos[iJetColl]-> hData_Tagger_MU->Fill( 2.5*JetInfo[iJetColl].Jet_SoftMu[JetInfo[iJetColl].nJet] );
    if ( JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet] < 0 ) Histos[iJetColl]-> hData_Tagger_MU->Fill( 2.5*JetInfo[iJetColl].Jet_SoftMuN[JetInfo[iJetColl].nJet] );


    if ( !isData_ )
    {
      Histos[iJetColl]->hAllFlav_Flavour->Fill( flavour );
      if ( varneg > 0 ) Histos[iJetColl]->hAllFlav_Tagger->Fill(-varneg );
      if ( varpos > 0 ) Histos[iJetColl]->hAllFlav_Tagger->Fill( varpos );
      if ( flavour == 1 || flavour == 21 ) { // light q+gluon
        if ( varneg > 0 ) Histos[iJetColl]->hLightFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hLightFlav_Tagger->Fill( varpos );
      } // light q+gluon
      if ( flavour == 21 ) { // gluon jets
        if ( varneg > 0 ) Histos[iJetColl]->hGluonFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hGluonFlav_Tagger->Fill( varpos );
      } // gluon jets
      else if ( flavour == 1 ) { // uds jets
        if ( varneg > 0 ) Histos[iJetColl]->hUDSFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hUDSFlav_Tagger->Fill( varpos );
      } // uds jets
      else if ( flavour == 4 ) { // c jets
        if ( varneg > 0 ) Histos[iJetColl]->hCFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hCFlav_Tagger->Fill( varpos );
      } // c jets
      else if ( flavour == 5 ) { // b jets
        if ( varneg > 0 ) Histos[iJetColl]->hBFlav_Tagger->Fill(-varneg );
        if ( varpos > 0 ) Histos[iJetColl]->hBFlav_Tagger->Fill( varpos );
      } // b jets
    }

    // Track History
    if ( useTrackHistory_ && indexes.size()!=0 && !isData_ ) {

      int cat1P = 0, cat2P = 0, cat3P = 0, catP = 0;
      int cat1N = 0, cat2N = 0, cat3N = 0, catN = 0;
      flags1P[TrackCategories::ConversionsProcess] ;
      if ( flags1P[TrackCategories::BWeakDecay] )              cat1P = 1;
      else if ( flags1P[TrackCategories::CWeakDecay] )         cat1P = 2;
      else if ( flags1P[TrackCategories::TauDecay] )           cat1P = 3;
      else if ( flags1P[TrackCategories::ConversionsProcess] ) cat1P = 4;
      else if ( flags1P[TrackCategories::KsDecay] )            cat1P = 5;
      else if ( flags1P[TrackCategories::LambdaDecay] )        cat1P = 6;
      else if ( flags1P[TrackCategories::HadronicProcess] )    cat1P = 7;
      else if ( flags1P[TrackCategories::Fake] )               cat1P = 8;
      else if ( flags1P[TrackCategories::SharedInnerHits] )    cat1P = 9;
      if(idSize > 1){
        if ( flags2P[TrackCategories::BWeakDecay] )              cat2P = 1;
        else if ( flags2P[TrackCategories::CWeakDecay] )         cat2P = 2;
        else if ( flags2P[TrackCategories::TauDecay] )           cat2P = 3;
        else if ( flags2P[TrackCategories::ConversionsProcess] ) cat2P = 4;
        else if ( flags2P[TrackCategories::KsDecay] )            cat2P = 5;
        else if ( flags2P[TrackCategories::LambdaDecay] )        cat2P = 6;
        else if ( flags2P[TrackCategories::HadronicProcess] )    cat2P = 7;
        else if ( flags2P[TrackCategories::Fake] )               cat2P = 8;
        else if ( flags2P[TrackCategories::SharedInnerHits] )    cat2P = 9;
      }
      if(idSize > 2){
        if ( flags3P[TrackCategories::BWeakDecay] )              cat3P = 1;
        else if ( flags3P[TrackCategories::CWeakDecay] )         cat3P = 2;
        else if ( flags3P[TrackCategories::TauDecay] )           cat3P = 3;
        else if ( flags3P[TrackCategories::ConversionsProcess] ) cat3P = 4;
        else if ( flags3P[TrackCategories::KsDecay] )            cat3P = 5;
        else if ( flags3P[TrackCategories::LambdaDecay] )        cat3P = 6;
        else if ( flags3P[TrackCategories::HadronicProcess] )    cat3P = 7;
        else if ( flags3P[TrackCategories::Fake] )               cat3P = 8;
        else if ( flags3P[TrackCategories::SharedInnerHits] )    cat3P = 9;
      }
      if ( flags1N[TrackCategories::BWeakDecay] )              cat1N = 1;
      else if ( flags1N[TrackCategories::CWeakDecay] )         cat1N = 2;
      else if ( flags1N[TrackCategories::TauDecay] )           cat1N = 3;
      else if ( flags1N[TrackCategories::ConversionsProcess] ) cat1N = 4;
      else if ( flags1N[TrackCategories::KsDecay] )            cat1N = 5;
      else if ( flags1N[TrackCategories::LambdaDecay] )        cat1N = 6;
      else if ( flags1N[TrackCategories::HadronicProcess] )    cat1N = 7;
      else if ( flags1N[TrackCategories::Fake] )               cat1N = 8;
      else if ( flags1N[TrackCategories::SharedInnerHits] )    cat1N = 9;
      if(idSize > 1){
        if ( flags2N[TrackCategories::BWeakDecay] )              cat2N = 1;
        else if ( flags2N[TrackCategories::CWeakDecay] )         cat2N = 2;
        else if ( flags2N[TrackCategories::TauDecay] )           cat2N = 3;
        else if ( flags2N[TrackCategories::ConversionsProcess] ) cat2N = 4;
        else if ( flags2N[TrackCategories::KsDecay] )            cat2N = 5;
        else if ( flags2N[TrackCategories::LambdaDecay] )        cat2N = 6;
        else if ( flags2N[TrackCategories::HadronicProcess] )    cat2N = 7;
        else if ( flags2N[TrackCategories::Fake] )               cat2N = 8;
        else if ( flags2N[TrackCategories::SharedInnerHits] )    cat2N = 9;
      }
      if(idSize > 2){
        if ( flags3N[TrackCategories::BWeakDecay] )              cat3N = 1;
        else if ( flags3N[TrackCategories::CWeakDecay] )         cat3N = 2;
        else if ( flags3N[TrackCategories::TauDecay] )           cat3N = 3;
        else if ( flags3N[TrackCategories::ConversionsProcess] ) cat3N = 4;
        else if ( flags3N[TrackCategories::KsDecay] )            cat3N = 5;
        else if ( flags3N[TrackCategories::LambdaDecay] )        cat3N = 6;
        else if ( flags3N[TrackCategories::HadronicProcess] )    cat3N = 7;
        else if ( flags3N[TrackCategories::Fake] )               cat3N = 8;
        else if ( flags3N[TrackCategories::SharedInnerHits] )    cat3N = 9;
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

      if ( varneg > 0 ) {
        if      ( catN == 1 ) Histos[iJetColl]->hAllFlav_Tagger_Bwd->Fill(-varneg );
        else if ( catN == 2 ) Histos[iJetColl]->hAllFlav_Tagger_Cwd->Fill(-varneg );
        else if ( catN == 3 ) Histos[iJetColl]->hAllFlav_Tagger_Tau->Fill(-varneg );
        else if ( catN == 4 ) Histos[iJetColl]->hAllFlav_Tagger_Gam->Fill(-varneg );
        else if ( catN == 5 ) Histos[iJetColl]->hAllFlav_Tagger_K0s->Fill(-varneg );
        else if ( catN == 6 ) Histos[iJetColl]->hAllFlav_Tagger_Lam->Fill(-varneg );
        else if ( catN == 7 ) Histos[iJetColl]->hAllFlav_Tagger_Int->Fill(-varneg );
        else if ( catN == 8 ) Histos[iJetColl]->hAllFlav_Tagger_Fak->Fill(-varneg );
        else if ( catN == 9 ) Histos[iJetColl]->hAllFlav_Tagger_Bad->Fill(-varneg );
        else                  Histos[iJetColl]->hAllFlav_Tagger_Oth->Fill(-varneg );
      }
      if ( varpos > 0 ) {
        if      ( catP == 1 ) Histos[iJetColl]->hAllFlav_Tagger_Bwd->Fill( varpos );
        else if ( catP == 2 ) Histos[iJetColl]->hAllFlav_Tagger_Cwd->Fill( varpos );
        else if ( catP == 3 ) Histos[iJetColl]->hAllFlav_Tagger_Tau->Fill( varpos );
        else if ( catP == 4 ) Histos[iJetColl]->hAllFlav_Tagger_Gam->Fill( varpos );
        else if ( catP == 5 ) Histos[iJetColl]->hAllFlav_Tagger_K0s->Fill( varpos );
        else if ( catP == 6 ) Histos[iJetColl]->hAllFlav_Tagger_Lam->Fill( varpos );
        else if ( catP == 7 ) Histos[iJetColl]->hAllFlav_Tagger_Int->Fill( varpos );
        else if ( catP == 8 ) Histos[iJetColl]->hAllFlav_Tagger_Fak->Fill( varpos );
        else if ( catP == 9 ) Histos[iJetColl]->hAllFlav_Tagger_Bad->Fill( varpos );
        else                  Histos[iJetColl]->hAllFlav_Tagger_Oth->Fill( varpos );
      }
    }

    ++numjet;
    ++JetInfo[iJetColl].nJet;
  } // end loop on jet

  Histos[iJetColl]->hData_All_NJets->Fill( JetInfo[iJetColl].nJet );
  Histos[iJetColl]->hData_NJets->Fill( numjet );

  return;
} // BTagAnalyzer:: processJets


float BTagAnalyzer::calculPtRel(const reco::Track& theMuon, const pat::Jet& theJet )
{
  double pmu = TMath::Sqrt( theMuon.px()*theMuon.px() + theMuon.py()*theMuon.py()  + theMuon.pz()*theMuon.pz() );

  double jetpx = 0 ;
  double jetpy = 0 ;
  double jetpz = 0 ;

  if( theJet.isCaloJet() ){
    jetpx = theJet.px() + theMuon.px();
    jetpy = theJet.py() + theMuon.py();
    jetpz = theJet.pz() + theMuon.pz();
  }else{
    jetpx = theJet.px();
    jetpy = theJet.py();
    jetpz = theJet.pz();
  }

  double jetp = TMath::Sqrt(jetpx*jetpx + jetpy*jetpy + jetpz*jetpz);

  double ptrel  = ( jetpx * theMuon.px()  + jetpy * theMuon.py() + jetpz * theMuon.pz() ) / jetp;
  ptrel = TMath::Sqrt( pmu * pmu  - ptrel * ptrel );

  return ptrel;
}


void BTagAnalyzer::setTracksPV( const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight )
{
  iPV = -1;
  PVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::VertexCollection::const_iterator IV;
  typedef reco::Vertex::trackRef_iterator IT;

  for(IV iv=pvHandle->begin(); iv!=pvHandle->end(); ++iv)
  {
    const reco::Vertex & vtx = *iv;
    // loop over tracks in vertices
    for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
    {
      const reco::TrackBaseRef & baseRef = *it;
      // one of the tracks in the vertex is the same as the track considered in the function
      if( baseRef == trackBaseRef )
      {
        float w = vtx.trackWeight(baseRef);
        // select the vertex for which the track has the highest weight
        if( w > PVweight )
        {
          PVweight = w;
          iPV = ( iv - pvHandle->begin() );
          break;
        }
      }
    }
  }

}


void BTagAnalyzer::setTracksSV( const reco::TrackRef & trackRef, const reco::SecondaryVertexTagInfo * svTagInfo, int & isFromSV, int & iSV, float & SVweight )
{
  isFromSV = 0;
  iSV = -1;
  SVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::Vertex::trackRef_iterator IT;

  size_t nSV = svTagInfo->nVertices();
  for(size_t iv=0; iv<nSV; ++iv)
  {
    const reco::Vertex & vtx = svTagInfo->secondaryVertex(iv);
    // loop over tracks in vertices
    for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
    {
      const reco::TrackBaseRef & baseRef = *it;
      // one of the tracks in the vertex is the same as the track considered in the function
      if( baseRef == trackBaseRef )
      {
        float w = vtx.trackWeight(baseRef);
        // select the vertex for which the track has the highest weight
        if( w > SVweight )
        {
          SVweight = w;
          isFromSV = 1;
          iSV = iv;
          break;
        }
      }
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void BTagAnalyzer::beginJob() {
  //cat 0
  cat0.etaMax = 2.5;
  cat0.etaMin = 0;
  cat0.nHitsMax= 50;
  cat0.nHitsMin= 8;
  cat0.nPixelHitsMax = 1;
  cat0.nPixelHitsMin = 1;
  cat0.pMax= 5000;
  cat0.pMin= 0;
  cat0.chiMin= 0;
  cat0.chiMax= 5;
  cat0.withFirstPixel = 0;

  //cat 1
  cat1.etaMax        = 2.5;
  cat1.etaMin        = 0;
  cat1.nHitsMax      = 50;
  cat1.nHitsMin      = 8;
  cat1.nPixelHitsMax = 8;
  cat1.nPixelHitsMin = 2;
  cat1.pMax          = 5000;
  cat1.pMin          = 0;
  cat1.chiMin        = 2.5;
  cat1.chiMax        = 5;
  cat1.withFirstPixel = 0;

  //cat 2
  cat2.etaMax        = 0.8;
  cat2.etaMin        = 0;
  cat2.nHitsMax      = 50;
  cat2.nHitsMin      = 8;
  cat2.nPixelHitsMax = 8;
  cat2.nPixelHitsMin = 3;
  cat2.pMax          = 8;
  cat2.pMin          = 0;
  cat2.chiMin        = 0;
  cat2.chiMax        = 2.5;
  cat2.withFirstPixel = 0;

  //cat 3
  cat3.etaMax        = 1.6;
  cat3.etaMin        = 0.8;
  cat3.nHitsMax      = 50;
  cat3.nHitsMin      = 8;
  cat3.nPixelHitsMax = 8;
  cat3.nPixelHitsMin = 3;
  cat3.pMax          = 8;
  cat3.pMin          = 0;
  cat3.chiMin        = 0;
  cat3.chiMax        = 2.5;
  cat3.withFirstPixel = 0;

  //cat 4
  cat4.etaMax        = 2.5;
  cat4.etaMin        = 1.6;
  cat4.nHitsMax      = 50;
  cat4.nHitsMin      = 8;
  cat4.nPixelHitsMax = 8;
  cat4.nPixelHitsMin = 3;
  cat4.pMax          = 8;
  cat4.pMin          = 0;
  cat4.chiMin        = 0;
  cat4.chiMax        = 2.5;
  cat4.withFirstPixel = 0;

  //cat 5
  cat5.etaMax        = 2.5;
  cat5.etaMin        = 0;
  cat5.nHitsMax      = 50;
  cat5.nHitsMin      = 8;
  cat5.nPixelHitsMax = 8;
  cat5.nPixelHitsMin = 2;
  cat5.pMax          = 8;
  cat5.pMin          = 0;
  cat5.chiMin        = 0;
  cat5.chiMax        = 2.5;
  cat5.withFirstPixel = 0;

  //cat 6
  cat6.etaMax        = 0.8;
  cat6.etaMin        = 0;
  cat6.nHitsMax      = 50;
  cat6.nHitsMin      = 8;
  cat6.nPixelHitsMax = 8;
  cat6.nPixelHitsMin = 3;
  cat6.pMax          = 5000;
  cat6.pMin          = 8;
  cat6.chiMin        = 0;
  cat6.chiMax        = 2.5;
  cat6.withFirstPixel = 0;

  //cat 7
  cat7.etaMax        = 1.6;
  cat7.etaMin        = 0.8;
  cat7.nHitsMax      = 50;
  cat7.nHitsMin      = 8;
  cat7.nPixelHitsMax = 8;
  cat7.nPixelHitsMin = 3;
  cat7.pMax          = 5000;
  cat7.pMin          = 8;
  cat7.chiMin        = 0;
  cat7.chiMax        = 2.5;
  cat7.withFirstPixel = 0;

  //cat 8
  cat8.etaMax        = 2.5;
  cat8.etaMin        = 1.6;
  cat8.nHitsMax      = 50;
  cat8.nHitsMin      = 8;
  cat8.nPixelHitsMax = 8;
  cat8.nPixelHitsMin = 3;
  cat8.pMax          = 5000;
  cat8.pMin          = 8;
  cat8.chiMin        = 0;
  cat8.chiMax        = 2.5;
  cat8.withFirstPixel = 0;

  //cat 9
  cat9.etaMax        = 2.5;
  cat9.etaMin        = 0;
  cat9.nHitsMax      = 50;
  cat9.nHitsMin      = 8;
  cat9.nPixelHitsMax = 2;
  cat9.nPixelHitsMin = 2;
  cat9.pMax          = 5000;
  cat9.pMin          = 8;
  cat9.chiMin        = 0;
  cat9.chiMax        = 2.5;
  cat9.withFirstPixel = 0;
}


// ------------ method called once each job just after ending the event loop  ------------
void BTagAnalyzer::endJob() {
}


std::vector< float > BTagAnalyzer::getTrackProbabilies(std::vector< float > v, const int ipType){

  std::vector< float > vectTrackProba;

  for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); ++q) {
    //positives and negatives tracks
    double p3d = 0;
    if(ipType == 0){
      if ( *q>=0){p3d= (*q)/2.;}else{p3d=1.+(*q)/2.;}
      //if(-log(p3d)> 5) p3d=exp(-5.0);
      vectTrackProba.push_back(p3d);
    }

    //positives tracks only
    if(ipType == 1 && *q >=0 ){
      vectTrackProba.push_back(*q);
    }

    //negatives tracks only
    if(ipType == 2 && *q <0){
      vectTrackProba.push_back(-(*q));
    }
  }

  return vectTrackProba;
}


double BTagAnalyzer::calculProbability(std::vector< float > v) {

  int ngoodtracks=v.size();
  double SumJet=0.;
  double m_minTrackProb = 0.005;
  for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); ++q) {
    SumJet += (*q>m_minTrackProb)?log(*q):log(m_minTrackProb);
  }

  double ProbJet;
  double Loginvlog = 0;

  if ( SumJet < 0. ) {
    if ( ngoodtracks >= 2 ) {
      Loginvlog = log(-SumJet);
    }
    double Prob = 1.;
    double lfact = 1.;
    for (int l=1; l!=ngoodtracks; ++l) {
      lfact *= l;
      Prob += exp( l*Loginvlog-log(1.*lfact) );
    }
    double LogProb = log(Prob);
    ProbJet = std::min(exp( std::max(LogProb+SumJet,-30.) ), 1.);
  } else {
    ProbJet = 1.;
  }
  return ProbJet;
}


bool BTagAnalyzer::findCat(const reco::Track* track, CategoryFinder& d) {
  //int numcat=-1;

  double   p = track->p();
  double eta = track->eta();
  double chi = track->normalizedChi2();
  int   nhit = track->numberOfValidHits();
  int   npix = track->hitPattern().numberOfValidPixelHits();

  bool result = ( p > d.pMin  &&  p  < d.pMax		  &&
      fabs(eta) > d.etaMin    &&  fabs(eta) < d.etaMax    &&
      nhit >= d.nHitsMin      &&  nhit <= d.nHitsMax	  &&
      npix >= d.nPixelHitsMin &&  npix <= d.nPixelHitsMax &&
      chi >= d.chiMin         &&  chi <= d.chiMax );

  return result;
}


int BTagAnalyzer::matchMuon(const edm::RefToBase<reco::Track>& theMuon, const edm::View<reco::Muon>& muons ){
  double small = 1.e-3;
  int matchedMu = -1;
  for(unsigned int i=0; i<muons.size(); ++i){
    double muonpt = -10000;
    if( muons[i].isGlobalMuon() )                               muonpt = muons[i].globalTrack()->pt() ;
    if(!muons[i].isGlobalMuon() &&  muons[i].isTrackerMuon())   muonpt = muons[i].innerTrack()->pt() ;

    if ( fabs(theMuon->pt() - muonpt )  < small  ){ matchedMu = i; }

  }
  return matchedMu;
}


std::vector<BTagAnalyzer::simPrimaryVertex> BTagAnalyzer::getSimPVs(const edm::Handle<edm::HepMCProduct>& evtMC){

  std::vector<BTagAnalyzer::simPrimaryVertex> simpv;
  const HepMC::GenEvent* evt=evtMC->GetEvent();
  if (evt) {

    for ( HepMC::GenEvent::vertex_const_iterator vitr= evt->vertices_begin();
        vitr != evt->vertices_end(); ++vitr ) { // loop for vertex ...
      HepMC::FourVector pos = (*vitr)->position();
      //HepLorentzVector pos = (*vitr)->position();

      bool hasMotherVertex=false;

      for ( HepMC::GenVertex::particle_iterator
          mother  = (*vitr)->particles_begin(HepMC::parents);
          mother != (*vitr)->particles_end(HepMC::parents);
          ++mother ) {
        HepMC::GenVertex * mv=(*mother)->production_vertex();
        if (mv) {
          hasMotherVertex=true;
          break; //if verbose_, print all particles of gen vertices
        }
      }

      if (hasMotherVertex) {continue;}

      // could be a new vertex, check  all primaries found so far to avoid multiple entries
      const double mm=0.1;
      simPrimaryVertex sv(pos.x()*mm,pos.y()*mm,pos.z()*mm);
      simPrimaryVertex *vp=NULL;  // will become non-NULL if a vertex is found and then point to it
      for (std::vector<simPrimaryVertex>::iterator v0=simpv.begin();
          v0!=simpv.end(); ++v0) {
        if ( (fabs(sv.x-v0->x)<1e-5) && (fabs(sv.y-v0->y)<1e-5)
            && (fabs(sv.z-v0->z)<1e-5) ) {
          vp=&(*v0);
          break;
        }
      }

      if (!vp) {	  // this is a new vertex
        simpv.push_back(sv);
        vp=&simpv.back();
      }
      vp->genVertex.push_back((*vitr)->barcode());
    }
  }
  return simpv;
}


int BTagAnalyzer::getMuonTk(double pt, const int iJetColl) {
  int idxTk = -1;
  for (int itk = 0; itk < JetInfo[iJetColl].nTrack ; ++itk) {
    if ( fabs(pt-JetInfo[iJetColl].Track_pt[itk]) < 1e-5 ) idxTk = itk;
  }
  return idxTk;
}


int BTagAnalyzer::isFromGSP(const reco::Candidate* c)
{
  int isFromGSP = 0;

  if( c->numberOfMothers() == 1 ) {
    const reco::Candidate* dau = c;
    const reco::Candidate* mom = c->mother();
    while( dau->numberOfMothers() == 1 && !( isHardProcess(mom->status()) && (abs(mom->pdgId())==4 || abs(mom->pdgId())==5) ) ) {
      if( abs(mom->pdgId())==21 )
      {
        isFromGSP = 1;
        break;
      }
      dau = mom;
      mom = dau->mother();
    }
  }

  return isFromGSP;
}


bool BTagAnalyzer::isHardProcess(const int status)
{
  // if Pythia8
  if( hadronizerType_ & (1 << 1) )
  {
    if( status>=21 && status<=29 )
      return true;
  }
  else // assuming Pythia6
  {
    if( status==3 )
      return true;
  }

  return false;
}

// -------------------------------------------------------------------------
// NameCompatible
// -------------------------------------------------------------------------
bool BTagAnalyzer::NameCompatible(const std::string& pattern, const std::string& name)
{
  const boost::regex regexp(edm::glob2reg(pattern));

  return boost::regex_match(name,regexp);
}


//define this as a plug-in
DEFINE_FWK_MODULE(BTagAnalyzer);



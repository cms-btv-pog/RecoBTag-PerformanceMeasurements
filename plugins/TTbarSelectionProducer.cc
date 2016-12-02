#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" 
#include "RecoBTag/PerformanceMeasurements/interface/TTbarSelectionProducer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace edm;



TTbarSelectionProducer::TTbarSelectionProducer(const edm::ParameterSet& iConfig) :
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerColl"))), 
  prunedGenParticleCollectionName_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  generatorevt_(consumes<GenEventInfoProduct>(edm::InputTag("generator",""))),
  generatorlhe_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  RecoHBHENoiseFilter_(consumes<bool>(iConfig.getParameter<edm::InputTag>("RecoHBHENoiseFilter"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxColl"))),
  bsToken_(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot",""))),
  electronToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronColl"))),
  conversionsToken_(mayConsume< reco::ConversionCollection >(iConfig.getParameter<edm::InputTag>("conversions"))),
  electronIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronIdMap"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonColl"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetColl"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metColl")))
{
  verbose_           = iConfig.getParameter<int > ("verbose");
  triggerBitsProc_    = iConfig.getParameter<edm::InputTag>("triggerColl").process();
  
  trigNamesToSel_     = iConfig.getParameter<std::vector<std::string> >("trigNamesToSel");
  trigChannels_       = iConfig.getParameter<std::vector<int> >("trigChannels");
  doTrigSel_          = iConfig.getParameter<bool>("doTrigSel");

  //MET filters
  std::vector<edm::InputTag> metFiltersInputs=iConfig.getParameter<std::vector<edm::InputTag> >("metFilters");
  for(size_t i=0; i<metFiltersInputs.size(); i++)
    metFilters_.push_back(consumes<edm::TriggerResults>(metFiltersInputs[i]));
  metFiltersToApply_ = iConfig.getParameter<std::vector<std::string> >("metFiltersToApply");

  //Configuration for electrons
  electron_cut_pt_   = iConfig.getParameter<double>        ("electron_cut_pt");
  electron_cut_eta_  = iConfig.getParameter<double>        ("electron_cut_eta");
  electron_cut_iso_  = iConfig.getParameter<double>        ("electron_cut_iso");
  
  //Configuration for muons
  muon_cut_pt_   = iConfig.getParameter<double>        ("muon_cut_pt");
  muon_cut_eta_  = iConfig.getParameter<double>        ("muon_cut_eta");
  muon_cut_iso_  = iConfig.getParameter<double>        ("muon_cut_iso");
  
  //Configuration for jets
  jet_cut_pt_   = iConfig.getParameter<double>        ("jet_cut_pt");
  jet_cut_eta_  = iConfig.getParameter<double>        ("jet_cut_eta");
  
  //Configuration for met
  met_cut_   = iConfig.getParameter<double>        ("met_cut");
  
  //produce
  produces<int>("topChannel");
  produces<int>("topTrigger");
  produces<int>("topMETFilter");
  produces<std::vector<reco::GenParticle> >();
  produces<std::vector<pat::Electron> >();
  produces<std::vector<pat::Muon> >();
  produces<std::vector<pat::Jet> >();
  produces<std::vector<pat::MET> >();

  //save some control histograms
  edm::Service<TFileService> fs;
  histos_["wgtcounter"] = fs->make<TH1F>("wgtcounter",";Weight;Weight sum",500,0,500);
  histos_["triggerpaths"] = fs->make<TH1F>("triggerpaths",";Trigger paths;Channel",trigChannels_.size(),0,trigChannels_.size());
  for(size_t i=0; i<trigChannels_.size(); i++)
    {
      histos_["triggerpaths"]->GetXaxis()->SetBinLabel(i+1,trigNamesToSel_[i].c_str());
      histos_["triggerpaths"]->SetBinContent(i+1,trigChannels_[i]);
    }

  histos_["metfilters"] = fs->make<TH1F>("metfilters",";MET filter;",metFiltersToApply_.size()+1,0,metFiltersToApply_.size()+1);
  histos_["metfilters"]->GetXaxis()->SetBinLabel(1,"HBHENoiseRECO");
  for(size_t i=0; i<metFiltersToApply_.size(); i++) histos_["metfilters"]->GetXaxis()->SetBinLabel(i+2,metFiltersToApply_[i].c_str());

  std::string ch[]={"inc","ee","mumu","emu","e","mu"};
  for(size_t i=0; i<sizeof(ch)/sizeof(std::string); i++)
    {
      histos_["cutflow_"+ch[i]] = fs->make<TH1F>(("cutflow_"+ch[i]).c_str(),"Selection level", 4, 0, 4);
      histos_["m_"+ch[i]]       = fs->make<TH1F>(("m_"+ch[i]).c_str(),      "Invariant mass [GeV]",200,0.,1000);
      histos_["met_"+ch[i]]     = fs->make<TH1F>(("met_"+ch[i]).c_str(),    "Missing transverse energy [GeV]", 100,0., 500);
      histos_["njets_"+ch[i]]   = fs->make<TH1F>(("njets_"+ch[i]).c_str(),  "Jet multiplicity",10,0.,10.);
    }
  for(std::map<std::string, TH1F *>::iterator it =histos_.begin();
      it!=histos_.end();
      it++)
    it->second->Sumw2();
}


TTbarSelectionProducer::~TTbarSelectionProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TTbarSelectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //keep track of generator level weights for normalization
   if(!iEvent.isRealData())
     {
       edm::Handle<GenEventInfoProduct> evt;
       iEvent.getByToken(generatorevt_,evt);
       float w0(1.0);
       if(evt.isValid()) w0=evt->weight();
       histos_["wgtcounter"]->Fill(0.,w0);
       
       edm::Handle<LHEEventProduct> evet;
       iEvent.getByToken(generatorlhe_,evet);
       if(evet.isValid())
	 {
	   float asdd=evet->originalXWGTUP();
	   for(unsigned int i=0  ; i<evet->weights().size();i++){
	     float asdde=evet->weights()[i].wgt;
	     float wi=w0*asdde/asdd;
	     histos_["wgtcounter"]->Fill(i+1,wi);
	   }
	 }
     }

   if(verbose_>5) std::cout << "in TTbarSelectionProducer::produce " << std::endl;

   //check trigger
   int trigWord(0);
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_, triggerBits);
   bool changedConfig = false;
   if (!hltConfig.init(iEvent.getRun(), iSetup, triggerBitsProc_, changedConfig))
     {
       std::cout << "Initialization of HLTConfigProvider failed!!" << std::endl;
       return;
     }   
   std::vector<bool> passTriggers(trigNamesToSel_.size(),false);
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) 
     {
       if(!triggerBits->accept(i)) continue;
       std::string trigName=hltConfig.triggerNames()[i];
       int tctr(0);
       for(std::vector<std::string>::iterator trigIt=trigNamesToSel_.begin();
	   trigIt!=trigNamesToSel_.end();
	   trigIt++, tctr++)
	 {
	   if(trigName.find(*trigIt)==std::string::npos) continue;
	   trigWord |= (1<<tctr);
	   passTriggers[tctr];
	 }
     }
   if(verbose_>5) std::cout << "Trigger word is: " << trigWord << std::endl;
   
   //disable trigger selection in MC, only trigger word will be stored
   bool isData=iEvent.isRealData();
   if(!isData) doTrigSel_=false;
   
   //MET filter decisions
   int metfilterWord(0);
   for(size_t im=0; im<metFilters_.size(); im++)
     {
       edm::Handle< edm::TriggerResults> metFilterBits;
       iEvent.getByToken(metFilters_[im],metFilterBits);
       if(!metFilterBits.isValid()) continue;

       //check accepted filters
       const edm::TriggerNames &metFilterNames = iEvent.triggerNames(*metFilterBits);
       for (unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
	 std::string name=metFilterNames.triggerName(i);
	 bool accept(metFilterBits->accept(i));
	 for(size_t j=0; j<metFiltersToApply_.size(); j++)
	   {
	     if(name.find(metFiltersToApply_[j])==std::string::npos) continue;
	     //std::cout << name << std::endl;
	     metfilterWord |= ((!accept) << (j+1));
	   }
       }
     }

   //for early runs
   if(isData)
     {
       edm::Handle<bool> HBHENoiseFilterResultHandle;
       iEvent.getByToken(RecoHBHENoiseFilter_, HBHENoiseFilterResultHandle);
       bool result(false);
       if(HBHENoiseFilterResultHandle.isValid()) result=*HBHENoiseFilterResultHandle;
       metfilterWord |= result;
     }
   
   //std::cout << metfilterWord << std::endl;

   // ----------------
   // Primary vertex
   //------------------ 
   edm::Handle<reco::VertexCollection> primaryVertices;
   iEvent.getByToken(vtxToken_, primaryVertices);
   const reco::Vertex &pVtx = *(primaryVertices->begin());

   //------------------------------------------
   //get beam spot
   //------------------------------------------
   Handle<reco::BeamSpot> bsHandle;
   iEvent.getByToken(bsToken_, bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   //------------------------------------------
   //Selection of muons
   //------------------------------------------
   edm::Handle<pat::MuonCollection> muHa;
   iEvent.getByToken(muonToken_, muHa);
   std::vector<pat::Muon> selMuons, selNonIsoMuons;
   for (const pat::Muon &mu : *muHa) 
     { 
       bool passKin( mu.pt() > muon_cut_pt_  && fabs(mu.eta()) < muon_cut_eta_ );
       if(!passKin) continue;

       //cf. https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
       //bool isMedium(muon::isMediumMuon(mu));
       bool isTight(muon::isTightMuon(mu,pVtx));
       bool passID(isTight);
       if(!passID) continue;

       double nhIso   = mu.neutralHadronIso();
       double puchIso = mu.puChargedHadronIso();
       double chIso   = mu.chargedHadronIso() ;
       double gIso    = mu.photonIso() ;
       double relIso  = (TMath::Max(Float_t(nhIso+gIso-0.5*puchIso),Float_t(0.))+chIso)/mu.pt();
       bool passIso( relIso < muon_cut_iso_ );
       if(!passIso)  continue;
       selMuons.push_back(mu);
     }
   if(verbose_>5) std::cout << "\t Selected  " << selMuons.size() << " muons" << std::endl;


   //------------------------------------------
   //Selection of electrons
   //------------------------------------------
   edm::Handle<edm::View<pat::Electron> > elHa;
   iEvent.getByToken(electronToken_, elHa);
   std::vector<pat::Electron> selElectrons;
   edm::Handle<reco::ConversionCollection> convHa;
   iEvent.getByToken(conversionsToken_, convHa);
   edm::Handle<edm::ValueMap<bool> > eIDHa;
   iEvent.getByToken(electronIdMapToken_ ,eIDHa);
   for (size_t i = 0; i < elHa->size(); ++i)
     {
       const auto el = elHa->ptrAt(i);

       bool passKin(el->pt() > electron_cut_pt_ && 
		    fabs(el->superCluster()->eta()) < electron_cut_eta_ && 
		    (el->isEB() || el->isEE()));
       if(!passKin) continue;

       // Conversion rejection
       bool passConvVeto = !ConversionTools::hasMatchedConversion(*el,convHa,beamspot.position());

       //cut-based electron id+iso
       //cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
       bool passElectronID = (*eIDHa)[el];
       bool passID( passConvVeto && passElectronID);
       if(!passID) continue;

       selElectrons.push_back(*el);
     }
   if(verbose_>5) std::cout << "\t Selected " << selElectrons.size() << " electrons" << std::endl;
   
   
   //------------------------------------------
   //Selection of jet
   //------------------------------------------
   edm::Handle<pat::JetCollection> jetHa;
   iEvent.getByToken(jetToken_, jetHa);
   std::vector<pat::Jet> selJets;
   for (const pat::Jet &j : *jetHa) 
     {
       bool passKin( j.pt() > jet_cut_pt_  && fabs(j.eta()) < jet_cut_eta_ );
       if(!passKin) continue;
	   
       // check overlap with electron and muon
       float minDR(99999.);
       for(size_t ilep=0; ilep<selMuons.size(); ilep++)
	 {
	   double dR = deltaR(j,selMuons[ilep]);
	   if(dR < minDR) minDR=dR;
	 }
       for(size_t ilep=0; ilep<selElectrons.size(); ilep++)
	 {
	   double dR = deltaR(j,selElectrons[ilep]);
	   if(dR < minDR) minDR=dR;
	 }
       bool hasOverlap(minDR<0.4);
       if(hasOverlap) continue;
       
       selJets.push_back(j);
     }
   if(verbose_>5) std::cout << "\t Selected "<< selJets.size() << " jets" << std::endl;
   

   edm::Handle<pat::METCollection> metHa;
   iEvent.getByToken(metToken_, metHa);
   const pat::MET &met = metHa->front();
   std::vector<pat::MET> selMETs(1,met);

   //assign channel
   float mll(0);
   int chsel=AssignChannel(selElectrons, selMuons, trigWord);
   if(verbose_>5) cout << "\t Channel assigned after lepton selection " << chsel << endl;

   bool passLepSel(chsel!=0);
   std::vector<std::string> tags(1,"inc");
   if(abs(chsel)==11*11) 
     {
       tags.push_back("ee");
       mll=(selElectrons[0].p4()+selElectrons[1].p4()).mass();
     }
   if(abs(chsel)==11*13) 
     {
       tags.push_back("emu");
       mll=(selElectrons[0].p4()+selMuons[0].p4()).mass();
     }
   if(abs(chsel)==13*13)
     {
       tags.push_back("mumu");
       mll=(selMuons[0].p4()+selMuons[1].p4()).mass();
     }
   if(abs(chsel)==13)    tags.push_back("mu");
   if(abs(chsel)==11)    tags.push_back("e");

   //jet selection
   bool passJetSel(false);
   if((abs(chsel)==13 || abs(chsel)==11) && selJets.size()>=4) passJetSel=true;
   if(abs(chsel)>13 && selJets.size()>=1)                      passJetSel=true;
   if(verbose_>5) std::cout << "\t Pass jet selection-" << passJetSel << std::endl;

   //MET selection
   bool passMetSel(true);
   if(abs(chsel)==11*11 || abs(chsel)==13*13) passMetSel=(met.pt()>met_cut_);
   if(verbose_>5) std::cout << "\t Pass met selection-" << passMetSel << std::endl;
   
   //fill control histos
   for(size_t i=0; i<tags.size(); i++)
     {
       histos_["cutflow_"+tags[i]]->Fill(0);
       if(passLepSel                            ) histos_["cutflow_"+tags[i]]->Fill(1);
       if(passLepSel && passJetSel              ) histos_["cutflow_"+tags[i]]->Fill(2);
       if(passLepSel && passJetSel && passMetSel) histos_["cutflow_"+tags[i]]->Fill(3);
       if(              passJetSel && passMetSel) histos_["m_"+tags[i]]->Fill(mll);
       if(passLepSel && passJetSel              ) histos_["met_"+tags[i]]->Fill(mll);
       if(passLepSel &&               passMetSel) histos_["njets_"+tags[i]]->Fill(selJets.size());
     }

   //save summary of selected objects into event
   if(!passLepSel || !passJetSel || !passMetSel) 
     {
       chsel=0;
       selElectrons.clear();
       selMuons.clear();
       selJets.clear();
       selMETs.clear();
       if(verbose_>5) std::cout << "\t Event is *not* selected " << std::endl;
     }
   else if(verbose_>5) std::cout << "\t Event is selected " << std::endl;


   //select particles from the hard process for MC
   std::vector<reco::GenParticle> selGen;
   if(!iEvent.isRealData())
     {
       edm::Handle<reco::GenParticleCollection> gpHa;
       iEvent.getByToken(prunedGenParticleCollectionName_,gpHa);
       int genChannel(1);
       for (const reco::GenParticle &g : *gpHa)
	 {
	   if(!g.isHardProcess()) continue;
	   if(abs(g.pdgId())==11 || abs(g.pdgId())==13) genChannel*=g.pdgId();
	   selGen.push_back(g);
	 }
       if(verbose_>5) std::cout << "\t gen level channel is " << genChannel << std::endl;
     }
   
   std::auto_ptr<int> trigWordOut( new int(trigWord) );
   iEvent.put(trigWordOut,"topTrigger");
   std::auto_ptr<int> metfilterWordOut( new int(metfilterWord) );
   iEvent.put(metfilterWordOut,"topMETFilter");
   std::auto_ptr<int > chOut( new int (chsel) );
   iEvent.put(chOut,"topChannel");
   std::auto_ptr< vector<reco::GenParticle> > genColl( new vector<reco::GenParticle>(selGen) );
   iEvent.put(genColl);
   auto_ptr<vector<pat::Electron> > eleColl( new vector<pat::Electron>(selElectrons) );
   iEvent.put( eleColl );
   auto_ptr<vector<pat::Muon> > muColl( new vector<pat::Muon>(selMuons) );
   iEvent.put( muColl );
   auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet>(selJets) );
   iEvent.put( jetColl );
   auto_ptr<vector<pat::MET> > metColl( new vector<pat::MET>(selMETs) );
   iEvent.put( metColl );
}

// ------------ method called once each job just before starting event loop  ------------
void
TTbarSelectionProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TTbarSelectionProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
TTbarSelectionProducer::beginRun(const edm::Run & iRun, edm::EventSetup const & iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TTbarSelectionProducer::endRun(const edm::Run & iRun, edm::EventSetup const & iSetup)
{
  try{
    edm::Service<TFileService> fs;
    
    edm::Handle<LHERunInfoProduct> lheruninfo;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByLabel( "externalLHEProducer", lheruninfo );
    
    LHERunInfoProduct myLHERunInfoProduct = *(lheruninfo.product());
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); 
	 iter!=myLHERunInfoProduct.headers_end(); 
	 iter++)
      {
	std::string tag("generator");
	if(iter->tag()!="") tag+="_"+iter->tag();
		
	std::vector<std::string> lines = iter->lines();
	std::vector<std::string> prunedLines;
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) 
	  {
	    if(lines.at(iLine)=="") continue;
	    if(lines.at(iLine).find("weightgroup")!=std::string::npos) continue;
	    prunedLines.push_back( lines.at(iLine) );
	  }
	
	if(histos_.find(tag)==histos_.end()) 
	  histos_[tag]=fs->make<TH1F>(tag.c_str(),tag.c_str(),prunedLines.size(),0,prunedLines.size());
	for (unsigned int iLine = 0; iLine<prunedLines.size(); iLine++) 
	  histos_[tag]->GetXaxis()->SetBinLabel(iLine+1,prunedLines.at(iLine).c_str());  
      }
  }
  catch(...){
    std::cout << "Failed to retrieve LHERunInfoProduct" << std::endl;
  }
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TTbarSelectionProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TTbarSelectionProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTbarSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



int
TTbarSelectionProducer::AssignChannel(std::vector<pat::Electron> &selElectrons,
				      std::vector<pat::Muon> &selMuons,
				      int trigWord)
{
  int chsel(0);

  //check which triggers fired
  bool triggerSingleMu(false), triggerSingleEle(false), triggerDoubleMu(false), triggerMuEG(false), triggerDoubleEle(false);
  if(doTrigSel_)
    {
      for(size_t i=0; i<trigChannels_.size(); i++)
	{
	  bool hasTrigger((trigWord>>i) & 0x1);
	  triggerSingleMu  |= ( abs(trigChannels_[i])==13 && hasTrigger);
	  triggerSingleEle |= ( abs(trigChannels_[i])==11 && hasTrigger);
	  triggerDoubleMu  |= ( abs(trigChannels_[i])==13*13 && hasTrigger);
	  triggerMuEG      |= ( abs(trigChannels_[i])==11*13 && hasTrigger);
	  triggerDoubleEle |= ( abs(trigChannels_[i])==11*11 && hasTrigger);
	}
    }
  else
    {
      triggerSingleMu=true;
      triggerSingleEle=true;
      triggerDoubleMu=true; 
      triggerMuEG=true; 
      triggerDoubleEle=true;
    }

  //assign dilepton channel
  if(selMuons.size()==1 && selElectrons.size()==0 && triggerSingleMu) chsel=selMuons[0].charge()*(-13);
  if(selMuons.size()==0 && selElectrons.size()==1 && triggerSingleEle) chsel=selElectrons[0].charge()*(-11);
  if(selMuons.size()+selElectrons.size()>=2)
    {
      float sumPt_ee(0.);
      int iee[2]={-1,-1};
      for(size_t i=0; i<selElectrons.size(); i++)
	for(size_t j=i+1; j<selElectrons.size(); j++)
	  {
	    float sumPt=selElectrons[i].pt()+selElectrons[j].pt();
	    if(sumPt<sumPt_ee) continue;
	    sumPt_ee=sumPt;
	    iee[0]=i; iee[1]=j;
	  }

      float sumPt_mm(0.);
      int imm[2]={-1,-1};
      for(size_t i=0; i<selMuons.size(); i++)
	for(size_t j=i+1; j<selMuons.size(); j++)
	  {
	    float sumPt=selMuons[i].pt()+selMuons[j].pt();
	    if(sumPt<sumPt_mm) continue;
	    sumPt_mm=sumPt;
	    imm[0]=i; imm[1]=j;
	  }

      float sumPt_em(0.);
      int iem[2]={-1,-1};
      for(size_t i=0; i<selMuons.size(); i++)
	for(size_t j=0; j<selElectrons.size(); j++)
	  {
	    float sumPt=selMuons[i].pt()+selElectrons[j].pt();
	    if(sumPt<sumPt_em) continue;
	    sumPt_em=sumPt;
	    iem[0]=i; iem[1]=j;
	  }

      if(sumPt_em>=sumPt_ee && sumPt_em>=sumPt_mm && triggerMuEG)
	{
	  std::vector<pat::Muon> newSelMuons(1,selMuons[iem[0]]);             selMuons=newSelMuons;
	  std::vector<pat::Electron> newSelElectrons(1,selElectrons[iem[1]]); selElectrons=newSelElectrons;
	  chsel=selElectrons[0].charge()*(-11)*selMuons[0].charge()*(-13);
	}
      else if(sumPt_mm>=sumPt_ee && sumPt_mm>=sumPt_em && triggerDoubleMu)
	{
	  selElectrons.clear();
	  std::vector<pat::Muon> newSelMuons(2); newSelMuons[0]=selMuons[imm[0]]; newSelMuons[1]=selMuons[imm[1]]; selMuons=newSelMuons;
	  chsel=selMuons[0].charge()*(-13)*selMuons[1].charge()*(-13);
	}
      else if(sumPt_ee>=sumPt_mm && sumPt_ee>=sumPt_em && triggerDoubleEle)
	{
	  selMuons.clear();
	  std::vector<pat::Electron> newSelElectrons(2); newSelElectrons[0]=selElectrons[iee[0]]; newSelElectrons[1]=selElectrons[iee[1]]; selElectrons=newSelElectrons;
	  chsel=selElectrons[0].charge()*(-11)*selElectrons[1].charge()*(-11);
	}
    }
  
  return chsel;
}




//define this as a plug-in
DEFINE_FWK_MODULE(TTbarSelectionProducer);

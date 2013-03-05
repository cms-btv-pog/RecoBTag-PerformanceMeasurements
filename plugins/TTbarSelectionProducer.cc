

#include "RecoBTag/PerformanceMeasurements/interface/TTbarSelectionProducer.h"

using namespace std;
using namespace edm;
 


TTbarSelectionProducer::TTbarSelectionProducer(const edm::ParameterSet& iConfig)
{
   //register your products
   
   isData_            = iConfig.getParameter<bool > ("isData");
   
   //Configuration for electrons
   
   electronColl_      = iConfig.getParameter<edm::InputTag> ("electronColl");
   electron_cut_pt_   = iConfig.getParameter<double>        ("electron_cut_pt");
   electron_cut_eta_  = iConfig.getParameter<double>        ("electron_cut_eta");
   electron_cut_iso_  = iConfig.getParameter<double>        ("electron_cut_iso");
   
   //Configuration for muons
   
   
   muonColl_      = iConfig.getParameter<edm::InputTag> ("muonColl");
   muon_cut_pt_   = iConfig.getParameter<double>        ("muon_cut_pt");
   muon_cut_eta_  = iConfig.getParameter<double>        ("muon_cut_eta");
   muon_cut_iso_  = iConfig.getParameter<double>        ("muon_cut_iso");
   
   //Configuration for jets 
   
   jetColl_      = iConfig.getParameter<edm::InputTag> ("jetColl");
   jet_cut_pt_   = iConfig.getParameter<double>        ("jet_cut_pt");
   jet_cut_eta_  = iConfig.getParameter<double>        ("jet_cut_eta");
   
   //Configuration for met 
   
   
   metColl_   = iConfig.getParameter<edm::InputTag> ("metColl");
   met_cut_   = iConfig.getParameter<double>        ("met_cut");
   
   
   trackColl_ = iConfig.getParameter<edm::InputTag> ("trackColl");
   
   
   //doBeamSpot_      = iConfig.getParameter<bool>	    ("doBeamSpot");  
   //beamSpotProducer_ = iConfig.getParameter<edm::InputTag>  ("beamSpotProducer");

   
   
   produces<int>();
//   produces<vector<TLorentzVector>>();
   produces<vector<double>>();
//   produces<double>();
   
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
   
   
   // some histograms  
   hcheck_cutflow        = fs->make<TH1F>("hcheck_cutflow","Selection level", 8, -0.5, 7.5);
   hcheck_m_ee           = fs->make<TH1F>("hcheck_m_ee","M_{e e}",200,0.,1000);
   hcheck_m_emu          = fs->make<TH1F>("hcheck_m_emu","M_{e #mu}",200,0.,1000);
   hcheck_m_mumu         = fs->make<TH1F>("hcheck_m_mumu","M_{#mu #mu}",200,0.,1000);
   hcheck_met_ee         = fs->make<TH1F>("hcheck_met_ee","MET (ee channel)", 100,0., 500);
   hcheck_met_emu        = fs->make<TH1F>("hcheck_met_emu","MET (e #mu channel)", 100,0., 500);
   hcheck_met_mumu       = fs->make<TH1F>("hcheck_met_mumu","MET (#mu #mu channel)", 100,0., 500);

   
   
   
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
  
   //std::cout << "in TTbarSelectionProducer::produce " << std::endl;

   
   //------------------------------------------
   //get bField
   //------------------------------------------
   
   float bField;
   
   if(isData_){
   
     edm::Handle<DcsStatusCollection> dcsHandle;
     iEvent.getByLabel("scalersRawToDigi", dcsHandle);
     
     float currentToBFieldScaleFactor = 2.09237036221512717e-04;
     float current = (*dcsHandle)[0].magnetCurrent();
     bField = current*currentToBFieldScaleFactor;
     
   }else{
     
      edm::ESHandle<MagneticField> magneticField;
      iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
      bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
   }

   //std::cout << "get bFiled : done " << std::endl;
   
    //------------------------------------------
   //get tracks, used for conversion
   //------------------------------------------
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(trackColl_, tracks);
   
   //std::cout << "get tracks : done " << std::endl;
   
   //------------------------------------------
   //get beam spot
   //------------------------------------------
   const reco::BeamSpot* bs = 0;
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByType(recoBeamSpotHandle);
   bs = recoBeamSpotHandle.product();
  
   //std::cout << "get beam spot : done " << std::endl;
   
   int channel = -1;
   int ind_cutflow=0;
   std::vector< TLorentzVector  > TheLeptons;
   double Mll_2 = -1;
   double themet = -1.;

   
   //------------------------------------------
   //Selection of muons
   //------------------------------------------
   
   std::vector< TLorentzVector  > p4Muon;
   std::vector< int  > chargeMuon;
   
   
   edm::Handle< std::vector<pat::Muon> >  muHa;
   iEvent.getByLabel(muonColl_, muHa);

   for (vector < pat::Muon >::const_iterator it = muHa->begin (); it != muHa->end (); it++){
     
     const pat::Muon * patmuon = &*it;
     if ( patmuon->pt() > muon_cut_pt_  && fabs(patmuon->eta()) < muon_cut_eta_ ){
       
       bool passmuonID = false;
       
       if(
         patmuon->isGlobalMuon() &&
         patmuon->isTrackerMuon()&&
	 patmuon->globalTrack()->normalizedChi2() < 10 &&
         patmuon->innerTrack()->numberOfValidHits() > 10 &&
         patmuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0  &&
         patmuon->innerTrack()->dxy(*bs)  < 0.02
	 ) passmuonID = true;
	 
	 
         double neutralHadronIso= patmuon->neutralHadronIso();
         double chargedHadronIso= patmuon->chargedHadronIso() ;
         double photonIso= patmuon->photonIso() ;
         
         double relIso =  ( neutralHadronIso+ chargedHadronIso + photonIso)/patmuon->pt() ;
    
         bool passIso = false;
	 if( relIso < muon_cut_iso_ ) passIso = true;
	 
         //std::cout << " muon " << patmuon->pt() << " " << fabs(patmuon->eta()) << " " << passmuonID << " " << passIso << std::endl ; 

	 if(passIso && passmuonID){
	   
           TLorentzVector themuon;
	   themuon.SetPtEtaPhiM(patmuon->pt(), patmuon->eta(),  patmuon->phi(), 0.);
	 
	   p4Muon.push_back(themuon);
	   chargeMuon.push_back(patmuon->charge());
	 }
     }
   } 
   
   //std::cout << "get muons : done " << chargeMuon.size() << std::endl;
   
   
   //------------------------------------------
   //Selection of electrons
   //------------------------------------------
   
   
   std::vector< TLorentzVector  > p4Elec;
   std::vector< int  > chargeElec;
   
   edm::Handle< std::vector<pat::Electron> >  elHa;
   iEvent.getByLabel(electronColl_, elHa);

   for (vector < pat::Electron >::const_iterator it = elHa->begin (); it != elHa->end (); it++){
     
     const pat::Electron * patelec = &*it;
     
     if ( patelec->pt() > electron_cut_pt_  && fabs(patelec->eta()) < electron_cut_eta_ ){
       
       
       double theta = 2*atan(exp(-1*patelec->superCluster()->eta()));
       double ET_SC = patelec->superCluster()->energy()*sin(theta);

       bool passTrackCut = false;
       if(   patelec->ecalDrivenSeed() 
          && patelec->gsfTrack().isNonnull()  
          && fabs(patelec->gsfTrack()->dxy(*bs)) < 0.04 
	  && ET_SC > 15 ) passTrackCut = true;
       
       
       // --------------------- Conversion ---------------------
       
       
       bool passConvReject = false;

/*
       ConversionFinder convFinder;
       ConversionInfo convInfo = convFinder.getConversionInfo(*patelec, tracks, bField);

       //std::cout << " conversion " << convInfo.dist() << " " << convInfo.dcot() << " " <<  patelec->gsfTrack()-> trackerExpectedHitsInner().numberOfLostHits() << std::endl;

       // fill informations 
       if( (convInfo.dist() >= 0.02 ||  convInfo.dcot() >= 0.02 ) && 
        patelec->gsfTrack()-> trackerExpectedHitsInner().numberOfLostHits() <2) passConvReject = true;;
*/
        if (patelec->passConversionVeto() && patelec->gsfTrack()-> trackerExpectedHitsInner().numberOfLostHits() <1) passConvReject = true;
       
       
       // --------------------- eID ---------------------
       bool passeID = false;
       const std::vector< std::pair<std::string,float> > patids  = patelec->electronIDs();
       for (unsigned int i=0;i<patids.size();i++){
	 //if(patids[i].first == "mvaTrigV0" && patids[i].second > 0 && patids[i].second < 1) passeID = true;
	 if(patids[i].first == "mvaTrigV0" && patids[i].second > 0.5 ) passeID = true;
         //std::cout << " id " << patids[i].first << std::endl;
       }
       
       double neutralHadronIso = patelec->neutralHadronIso();
       double chargedHadronIso = patelec->chargedHadronIso() ;
       double photonIso = patelec->photonIso() ;
       double relIso =  ( neutralHadronIso+ chargedHadronIso + photonIso)/patelec->pt() ;
       
       bool passIso = false;
       if(relIso < electron_cut_iso_) passIso = true;

       //std::cout << " electron " << patelec->pt() << " " << fabs(patelec->eta()) << " " << passTrackCut << " " << passConvReject << " " << passIso << " " << passeID << std::endl;
       
       if(passIso && passeID && passConvReject &&passTrackCut ){
       
          TLorentzVector theelectron;
	  theelectron.SetPtEtaPhiM(patelec->pt(), patelec->eta(),  patelec->phi(), 0.);
	 
	  p4Elec.push_back(theelectron);
	  chargeElec.push_back(patelec->charge());
       }
     }
   } 
	
   
   //std::cout << "get electrons : done " << chargeElec.size() << std::endl;
   
   
   //------------------------------------------
   //Selection of jet
   //------------------------------------------
   
   std::vector< TLorentzVector  > p4Jet;
   
   
   edm::Handle< std::vector<pat::Jet> >  jetHa;
   iEvent.getByLabel(jetColl_, jetHa);
   
   
   int nSelJets=0;
   
   for (vector < pat::Jet >::const_iterator it = jetHa->begin (); it != jetHa->end (); it++){
     
     const pat::Jet * patjet = &*it;
     if ( patjet->pt() > jet_cut_pt_  && fabs(patjet->eta()) < jet_cut_eta_ ){
       
       // check overlap with electron and muon 
      TLorentzVector thejet;
      thejet.SetPtEtaPhiM(patjet->pt(), patjet->eta(),  patjet->phi(), 0.);
      double deltaRmu = 10000;
      double deltaRel = 10000;

      for(unsigned int imu=0; imu< p4Muon.size(); imu++)
      {
        double deltaR = thejet.DeltaR(p4Muon[imu]);
        if(deltaR < deltaRmu) deltaRmu = deltaR;
      }

      for(unsigned int iel=0; iel< p4Elec.size(); iel++)
      {
        double deltaR = thejet.DeltaR(p4Elec[iel]);
        if(deltaR < deltaRel) deltaRel = deltaR;
      }

      if( deltaRmu > 0.5  && deltaRel > 0.5) {
        nSelJets++;
	p4Jet.push_back(thejet);
      }
       
       
     }
     
     
   } 
   
   //std::cout << "get jets : done " << nSelJets << std::endl;
   
   edm::Handle< std::vector<pat::MET> >  metHa;
   iEvent.getByLabel(metColl_, metHa);

   /*for (vector < pat::MET >::const_iterator it = metHa->begin (); it != metHa->end (); it++){
     
     const pat::MET * patmet = &*it;
     if ( patmet->pt() > met_cut_  ){
       std::cout << "in loop on met" << std::endl;
     }
   } */
  
   //std::cout << "get met : done " << std::endl;
  
  bool passSel = false;
  if ((p4Elec.size()+p4Muon.size()) >=1) ind_cutflow++;
  if( (p4Elec.size()+p4Muon.size()) >=2){
    
    ind_cutflow++;  
    int idxLept1 = -1;
    int idxLept2 = -1;
    //std::cout << " Lepton : " << p4Elec.size()+p4Muon.size() << std::endl;
    GetLeptonPair(p4Elec, p4Muon, chargeElec, chargeMuon, idxLept1, idxLept2, channel);
    
    if(channel >=0){
    
      ind_cutflow++;
      double Mll = -1;
      if(channel == 0) Mll = (p4Elec[idxLept1]+p4Elec[idxLept2]).M(); 
      if(channel == 1) Mll = (p4Muon[idxLept1]+p4Muon[idxLept2]).M();
      if(channel == 2) Mll = (p4Elec[idxLept1]+p4Muon[idxLept2]).M();

      if(channel == 0)  {
         TheLeptons.push_back(p4Elec[idxLept1]);
         TheLeptons.push_back(p4Elec[idxLept2]);
      }
      else if (channel == 1) {
         TheLeptons.push_back(p4Muon[idxLept1]);
         TheLeptons.push_back(p4Muon[idxLept2]);
      }
      else if (channel == 2) {
         TheLeptons.push_back(p4Elec[idxLept1]);
         TheLeptons.push_back(p4Muon[idxLept2]);
      }


      Mll_2 = (TheLeptons[0]+TheLeptons[1]).M();
      if (fabs(Mll_2-Mll)>0.01) std::cout << " Mll_2 " << Mll_2 << " Mll " << Mll << std::endl;
      //std::cout << " Lepton pair : " << Mll << std::endl;
      
      const pat::MET    *met = 0;
      met = &(metHa->front());
      themet = sqrt(pow(met->px(), 2) + pow(met->px(), 2) );
      
      
      if( Mll > 20 && ( channel ==2 || ( channel <=1 && (Mll < 76 || Mll > 116 )  ) ) ) {
        ind_cutflow++;
	if (nSelJets >=2) {
         ind_cutflow++;
         if (themet >  met_cut_  ||  channel ==2) {
             passSel = true;
             ind_cutflow++;
             if (channel==0) {
                 hcheck_m_ee->Fill(Mll)   ;
                 hcheck_met_ee->Fill(themet) ;
             }
             else if (channel==1) {
                 hcheck_m_mumu->Fill(Mll)   ;
                 hcheck_met_mumu->Fill(themet) ;
             }
             else if (channel==2) {
                 hcheck_m_emu->Fill(Mll)   ;
                 hcheck_met_emu->Fill(themet) ;
/*
                 std::cout << " TTbarProducer " << std::endl; 
                 cout << " lepton 1 " << TheLeptons[0].Pt() << " " << TheLeptons[0].Eta() << endl;
                 cout << " lepton 2 " << TheLeptons[1].Pt() << " " << TheLeptons[1].Eta() << endl;
                 for(unsigned int ij=0; ij< p4Jet.size(); ij++) {
                   cout << " jet " << ij << "   " << p4Jet[ij].Pt() << " " << p4Jet[ij].Eta() << endl;
                 }
                 cout << " " << endl;
*/
             }
         }
       }
     }

	 
	 
    }
    
    
  }
   //std::cout << "get selection : done" << passSel << std::endl;
  
   if(!passSel) channel = -1;
   hcheck_cutflow->Fill(0)        ;
   if (ind_cutflow>0) {
    for (int ii=1; ii<=ind_cutflow; ii++) {
     hcheck_cutflow->Fill(ii)        ;
    }
   }
   //std::cout << " selection " << passSel << " channel " << channel << std::endl;
   
  
  std::auto_ptr<int > pOut( new int (channel) );
  iEvent.put(pOut);

  vector<double> thelep_and_met;
  if (channel>=0) {
   thelep_and_met.push_back(TheLeptons[0].Pt());
   thelep_and_met.push_back(TheLeptons[0].Eta());
   thelep_and_met.push_back(TheLeptons[0].Phi());
   thelep_and_met.push_back(TheLeptons[1].Pt());
   thelep_and_met.push_back(TheLeptons[1].Eta());
   thelep_and_met.push_back(TheLeptons[1].Phi());
   thelep_and_met.push_back(themet);
   thelep_and_met.push_back(Mll_2);
  }
  std::auto_ptr<std::vector<double>> pOut2 (new std::vector<double> (thelep_and_met) );
  iEvent.put(pOut2);
/*
  std::auto_ptr<std::vector<TLorentzVector>> pOut2 (new std::vector<TLorentzVector> (TheLeptons) );
  iEvent.put(pOut2);

  std::auto_ptr<double > pOut3( new double (themet) );
  iEvent.put(pOut3);
*/
  
  //std::cout << "print output : done " << std::endl;


 
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
TTbarSelectionProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TTbarSelectionProducer::endRun(edm::Run&, edm::EventSetup const&)
{
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



void
TTbarSelectionProducer::GetLeptonPair(
                          std::vector<TLorentzVector> elec_in, std::vector<TLorentzVector> muon_in, 
                          std::vector<int> elec_charge, std::vector<int> muon_charge, 
			  int &idxLept1, int &idxLept2, int &thechannel){
			  
			  
  float sum_pT_ee = 0.;
  bool pass_elec = false;
  int ie1 = -1;
  int ie2 = -1;
  if (elec_in.size () >= 2) {
    for (unsigned int i = 0; i < elec_in.size (); i++) {
      for (unsigned int j = i + 1; j < elec_in.size (); j++) {
	if (pass_elec)
	  continue;
	if ( elec_charge[i] != elec_charge[j] ){
	  pass_elec = true;
	  sum_pT_ee = elec_in[i].Pt () + elec_in[j].Pt ();
	  ie1 = i;
	  ie2 = j;
	}
      }
    }
  }

  float sum_pT_mumu = 0.;
  bool pass_muon = false;
  int imu1 = -1;
  int imu2 = -1;
  if (muon_in.size () >= 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = i + 1; j < muon_in.size (); j++) {
	if (pass_muon)
	  continue;
	if ( muon_charge[i] != muon_charge[j] ){
	  pass_muon = true;
	  sum_pT_mumu = muon_in[i].Pt () + muon_in[j].Pt ();
	  imu1 = i;
	  imu2 = j;
	}
      }
    }
  }


  float sum_pT_emu_start = 0.;
  float sum_pT_emu = 0.;
  int je1 = -1;
  int jmu2 = -1;
  if (muon_in.size () >= 1 && elec_in.size () >= 1) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = 0; j < elec_in.size (); j++) {
	if ( (muon_charge[i] != elec_charge[j]) ){
	  sum_pT_emu = muon_in[i].Pt () + elec_in[j].Pt ();
	  if (sum_pT_emu > sum_pT_emu_start) {
	    sum_pT_emu_start = sum_pT_emu;
	    je1 = j;
	    jmu2 = i;
	  }
	}
      }
    }
  }


  float sum[3] = { sum_pT_ee, sum_pT_mumu, sum_pT_emu };
  int sortedIndices[3];
  TMath::Sort (3, sum, sortedIndices);
  if (sortedIndices[0] == 0 && sum_pT_ee != 0.) {
    idxLept1 = ie1;
    idxLept2 = ie2;
    thechannel = 0;
  }
  else if (sortedIndices[0] == 1 && sum_pT_mumu != 0.) {
    idxLept1 = imu1;
    idxLept2 = imu2;
    thechannel = 1;
  }
  else if (sortedIndices[0] == 2 && sum_pT_emu != 0.) {
    idxLept1 = je1;
    idxLept2 = jmu2;
    thechannel = 2;
  }
  
  
  


		  
}




//define this as a plug-in
DEFINE_FWK_MODULE(TTbarSelectionProducer);

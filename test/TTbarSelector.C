#include "TTbarSelector.h"

using namespace std;

bool TTbarSelector::passTTbarSelection(bool isData, vector<TLorentzVector> theLeptColl, vector<Int_t> theLeptIds, vector< pair< TLorentzVector, Float_t> > theJetColl, Int_t ttbar_trigWord, Float_t ttbar_w[250], int ttbar_nw, TH1F* wgtcounter, TString syst, bool computeEvtWgtOnly)
{

        bool applyTriggerEff_ = false;
        bool applyLepSelEff_ = false;

        if(!isData) applyTriggerEff_ = true;
        if(!isData) applyLepSelEff_ = true;

        // Check of µµ/µe/ee channel
        Int_t theChannel = 1;
        for (unsigned short int ih = 0; ih < theLeptColl.size(); ih++)
        {
                theChannel *= theLeptIds[ih];
        }

        // Apply trigger efficiency
        Float_t trigWgtLo(1.0), trigWgtNom(1.0), trigWgtHi(1.0);
        if(applyTriggerEff_)
	{
	       pair<float,float> eff = getTriggerEfficiency(theChannel);
	       trigWgtLo       *= eff.first-eff.second;
	       trigWgtNom      *= eff.first;
	       trigWgtHi       *= eff.first+eff.second;
	}

        // Apply lepton selection efficiency
        Float_t lepSelEffLo(1.0), lepSelEffNom(1.0), lepSelEffHi(1.0);
        if(applyLepSelEff_)
	{
	        for(size_t il=0; il< theLeptColl.size(); il++)
	        {
	                pair<float,float> lepSF = getLeptonSelectionEfficiencyScaleFactor(theLeptIds[il], theLeptColl[il].Pt(), theLeptColl[il].Eta() );
	                lepSelEffLo  *= (lepSF.first-lepSF.second);
	                lepSelEffNom *= lepSF.first;
	                lepSelEffHi  *= (lepSF.first+lepSF.second);
	        }
	}

        Float_t qcdScaleLo(1.0),qcdScaleHi(1.0);
        if(!isData && ttbar_nw>17)
        {
                qcdScaleLo = ttbar_w[9]*(wgtcounter->GetBinContent(10)/wgtcounter->GetBinContent(1));
                qcdScaleHi = ttbar_w[5]*(wgtcounter->GetBinContent(6)/wgtcounter->GetBinContent(1));
        }

        // nominal event weight
        evWgt = 1.0;
        if(!isData)
        {
                if(ttbar_nw != 0) evWgt *= trigWgtNom*lepSelEffNom*ttbar_w[0];
                else              evWgt *= trigWgtNom*lepSelEffNom;

                //weights for systematic uncertainties
                if     (syst == "trig__minus")                    evWgt *= trigWgtLo/trigWgtNom;
                else if(syst == "trig__plus")                     evWgt *= trigWgtHi/trigWgtNom;
                else if(syst == "lept__minus")                    evWgt *= lepSelEffLo/lepSelEffNom;
                else if(syst == "lept__plus")                     evWgt *= lepSelEffHi/lepSelEffNom;
                else if(syst == "scale1__plus"  && ttbar_nw > 0)  evWgt *= qcdScaleHi/ttbar_w[0];
                else if(syst == "scale1__minus" && ttbar_nw > 0)  evWgt *= qcdScaleLo/ttbar_w[0];
        }

        // if we do not want to compute the selection but only the weight calculation (for PU reweighting)
        if(computeEvtWgtOnly) return true;

        // Dilepton cut
        if ( theLeptColl.size()  != 2) return false;

        if (theChannel != -13*11) return false;

        int nPassPt=0; // Additional offline pt cut (tight)
        for( unsigned short int ilep = 0; ilep < theLeptColl.size(); ilep++ ){ if( theLeptColl[ilep].Pt() > 25 ) nPassPt++; }
        if( nPassPt != 2 ) return false;

        // Pass trigger
        bool passTrigger_ = passTrigger(isData, theChannel, ttbar_trigWord);
        if(!passTrigger_) return false;


        // Jet cleaning
        theSelJetColl.clear();

        for (unsigned short int ijet = 0; ijet < theJetColl.size(); ijet++)
        {
                // Jet Cleaning
                if(theLeptColl[0].DeltaR(theJetColl[ijet].first) < 0.4 || theLeptColl[1].DeltaR(theJetColl[ijet].first) < 0.4 ) continue;

        	// JES 
	        vector<float> jesSF(3,1.0);
	        theJECuncertainty_->setJetEta(fabs(theJetColl[ijet].first.Eta()));
	        theJECuncertainty_->setJetPt(theJetColl[ijet].first.Pt());
	        float unc = theJECuncertainty_->getUncertainty(true);
	        jesSF[1]=(1.+fabs(unc));
	        jesSF[2]=(1.-fabs(unc));
	 
                // JER 
	        vector<float> jerSF= getJetResolutionScales(jesSF[0]*theJetColl[ijet].first.Pt(), theJetColl[ijet].first.Eta(), theJetColl[ijet].second);

                // Save a copy of the jet
	        TLorentzVector initialJet(theJetColl[ijet].first);
	        TLorentzVector nominalJet(initialJet*jesSF[0]*jerSF[0]);

	        // apply energy shifts according to systematic variation
	        TLorentzVector variedJet = nominalJet;
	                
	        if(syst == "jes__plus")  variedJet *= jesSF[1]/jesSF[0];
	        if(syst == "jes__minus") variedJet *= jesSF[2]/jesSF[0];
	        if(syst == "jer__minus") variedJet *= jerSF[1]/jerSF[0];
	        if(syst == "jer__plus")  variedJet *= jerSF[2]/jerSF[0];

	      
	        //check if can be selected for this variation
	        if(variedJet.Pt() < 20 || TMath::Abs(variedJet.Eta()) > 2.4) continue;

                theSelJetColl.push_back(make_pair(variedJet, ijet));

        }

        // Dijet cut
        if ( theSelJetColl.size() < 2) return false;

        // dilepton invariant mass
        TLorentzVector dilepton = theLeptColl[0] + theLeptColl[1];
        mll_ = dilepton.M();
        if(mll_ < 12) return false;


        return true;                              
}


bool TTbarSelector::passTrigger(bool isData, Int_t ttbar_chan, Int_t ttbar_trigWord)
{
    vector<std::pair<Int_t,Int_t> > trigBits;
    if(isData)
    {
        if(ttbar_chan == -13*11) trigBits.push_back(make_pair( 0,ttbar_chan)); //Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL
        if(ttbar_chan == -13*11) trigBits.push_back(make_pair( 1,ttbar_chan)); //Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
        if(ttbar_chan == -13*11) trigBits.push_back(make_pair( 2,ttbar_chan)); //Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
    }            
    else
    {
        trigBits.push_back(make_pair( 0,-13*11)); //Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL
        trigBits.push_back(make_pair( 1,-13*11)); //Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
        trigBits.push_back(make_pair( 2,-13*11)); //Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
    }


    bool hasTrigger;
    if( trigBits.size()==0) hasTrigger = true;
    else hasTrigger = false;

    if( isData )
    {
        for(unsigned int ibit=0; ibit<trigBits.size(); ibit++)
        {
            if(trigBits[ibit].second != ttbar_chan) continue;
            hasTrigger |= ((ttbar_trigWord>>trigBits[ibit].first) & 1);
        }
    }else{
        hasTrigger=true;
    }

    return hasTrigger;

}

//////////////////////////
//  CTag Selection    ////
/////////////////////////
bool TTbarSelector::passSemiLepTTbarSelection(bool isData, vector<TLorentzVector> theLeptColl, vector<Int_t> theLeptIds, vector< pair< TLorentzVector, Float_t> > theJetColl, Int_t ttbar_trigWord, Float_t ttbar_w[250], int ttbar_nw, TH1F* wgtcounter, TString syst, bool computeEvtWgtOnly)
{

        bool applyTriggerEff_ = false;
        bool applyLepSelEff_ = false;

        if(!isData) applyTriggerEff_ = true;
        if(!isData) applyLepSelEff_ = true;

        // Check of mu/e channel
        Int_t theChannel = 1;
        for (unsigned short int ih = 0; ih < theLeptColl.size(); ih++)
        {
                theChannel = theLeptIds[ih];
        }

        // Apply trigger efficiency
        Float_t trigWgtLo(1.0), trigWgtNom(1.0), trigWgtHi(1.0);
        if(applyTriggerEff_)
	{
	       pair<float,float> eff = getTriggerEfficiency(theChannel);
	       trigWgtLo       *= eff.first-eff.second;
	       trigWgtNom      *= eff.first;
	       trigWgtHi       *= eff.first+eff.second;
	}

        // Apply lepton selection efficiency
        Float_t lepSelEffLo(1.0), lepSelEffNom(1.0), lepSelEffHi(1.0);
        if(applyLepSelEff_)
	{
	        for(size_t il=0; il< theLeptColl.size(); il++)
	        {
	                pair<float,float> lepSF = getLeptonSelectionEfficiencyScaleFactor(theLeptIds[il], theLeptColl[il].Pt(), theLeptColl[il].Eta() );
	                lepSelEffLo  *= (lepSF.first-lepSF.second);
	                lepSelEffNom *= lepSF.first;
	                lepSelEffHi  *= (lepSF.first+lepSF.second);
	        }
	}

        Float_t qcdScaleLo(1.0),qcdScaleHi(1.0);
        if(!isData && ttbar_nw>17)
        {
                qcdScaleLo = ttbar_w[9]*(wgtcounter->GetBinContent(10)/wgtcounter->GetBinContent(1));
                qcdScaleHi = ttbar_w[5]*(wgtcounter->GetBinContent(6)/wgtcounter->GetBinContent(1));
        }

        // nominal event weight
        evWgt = 1.0;
        if(!isData)
        {
                if(ttbar_nw != 0) evWgt *= trigWgtNom*lepSelEffNom*ttbar_w[0];
                else              evWgt *= trigWgtNom*lepSelEffNom;

                //weights for systematic uncertainties
                if     (syst == "trig__minus")                    evWgt *= trigWgtLo/trigWgtNom;
                else if(syst == "trig__plus")                     evWgt *= trigWgtHi/trigWgtNom;
                else if(syst == "lept__minus")                    evWgt *= lepSelEffLo/lepSelEffNom;
                else if(syst == "lept__plus")                     evWgt *= lepSelEffHi/lepSelEffNom;
                else if(syst == "scale1__plus"  && ttbar_nw > 0)  evWgt *= qcdScaleHi/ttbar_w[0];
                else if(syst == "scale1__minus" && ttbar_nw > 0)  evWgt *= qcdScaleLo/ttbar_w[0];
        }

        // if we do not want to compute the selection but only the weight calculation (for PU reweighting)
        if(computeEvtWgtOnly) return true;

        // SingleLepton cut
        if ( theLeptColl.size()  != 1) return false;

        int nPassPt=0; // Additional offline pt cut (tight)

        if (abs(theChannel) != 13 ) return false; //SingleMu
        for( unsigned short int ilep = 0; ilep < theLeptColl.size(); ilep++ ){ if( theLeptColl[ilep].Pt() > 25 ) nPassPt++; }
        if( nPassPt != 1 ) return false;

        //if (abs(theChannel) != 11 ) return false; //SingleElectron
        //for( unsigned short int ilep = 0; ilep < theLeptColl.size(); ilep++ ){ if( theLeptColl[ilep].Pt() > 40 ) nPassPt++; }
        //if( nPassPt != 1 ) return false;       

        // Pass trigger
        bool passSingleTrigger_ = passSingleTrigger(isData, theChannel, ttbar_trigWord);
        if(!passSingleTrigger_) return false;


        // Jet cleaning
        theSelJetColl.clear();

        for (unsigned short int ijet = 0; ijet < theJetColl.size(); ijet++)
        {
                // Jet Cleaning
                if(theLeptColl[0].DeltaR(theJetColl[ijet].first) < 0.4 ) continue;

        	// JES 
	        vector<float> jesSF(3,1.0);
	        theJECuncertainty_->setJetEta(fabs(theJetColl[ijet].first.Eta()));
	        theJECuncertainty_->setJetPt(theJetColl[ijet].first.Pt());
	        float unc = theJECuncertainty_->getUncertainty(true);
	        jesSF[1]=(1.+fabs(unc));
	        jesSF[2]=(1.-fabs(unc));
	 
                // JER 
	        vector<float> jerSF= getJetResolutionScales(jesSF[0]*theJetColl[ijet].first.Pt(), theJetColl[ijet].first.Eta(), theJetColl[ijet].second);

                // Save a copy of the jet
	        TLorentzVector initialJet(theJetColl[ijet].first);
	        TLorentzVector nominalJet(initialJet*jesSF[0]*jerSF[0]);

	        // apply energy shifts according to systematic variation
	        TLorentzVector variedJet = nominalJet;
	                
	        if(syst == "jes__plus")  variedJet *= jesSF[1]/jesSF[0];
	        if(syst == "jes__minus") variedJet *= jesSF[2]/jesSF[0];
	        if(syst == "jer__minus") variedJet *= jerSF[1]/jerSF[0];
	        if(syst == "jer__plus")  variedJet *= jerSF[2]/jerSF[0];

	      
	        //check if can be selected for this variation
	        if(variedJet.Pt() < 25 || TMath::Abs(variedJet.Eta()) > 2.4) continue;

                theSelJetColl.push_back(make_pair(variedJet, ijet));

        }

        // Four jets cut
        if ( theSelJetColl.size() < 4) return false;

        return true;                              
}

bool TTbarSelector::passSingleTrigger(bool isData, Int_t ttbar_chan, Int_t ttbar_trigWord)
{
    vector<std::pair<Int_t,Int_t> > trigBits;
    if(isData)
    {
        if(ttbar_chan == -11) trigBits.push_back(make_pair( 5,ttbar_chan)); //Ele35_WPLoose_Gsf_v
	if(ttbar_chan ==  11) trigBits.push_back(make_pair( 5,ttbar_chan)); //Ele35_WPLoose_Gsf_v
        if(ttbar_chan == -13) trigBits.push_back(make_pair( 3,ttbar_chan)); //IsoMu22_v
        if(ttbar_chan ==  13) trigBits.push_back(make_pair( 3,ttbar_chan)); //IsoMu22_v
//        if(ttbar_chan == -11) trigBits.push_back(make_pair( 6,ttbar_chan)); //Ele45_WPLoose_Gsf_v
//        if(ttbar_chan ==  11) trigBits.push_back(make_pair( 6,ttbar_chan)); //Ele45_WPLoose_Gsf_v
//        if(ttbar_chan == -13) trigBits.push_back(make_pair( 4,ttbar_chan)); //IsoMu24_v
//        if(ttbar_chan ==  13) trigBits.push_back(make_pair( 4,ttbar_chan)); //IsoMu24_v
    }
    else
    {
        trigBits.push_back(make_pair( 5,-11));
        trigBits.push_back(make_pair( 5, 11));
        trigBits.push_back(make_pair( 3,-13));
        trigBits.push_back(make_pair( 3, 13));
//	  trigBits.push_back(make_pair( 6,-11));
//        trigBits.push_back(make_pair( 6, 11));
//        trigBits.push_back(make_pair( 3,-13));
//        trigBits.push_back(make_pair( 3, 13));
    }


    bool hasTrigger;
    if( trigBits.size()==0) hasTrigger = true;
    else hasTrigger = false;

   if( isData )
    {
        for(unsigned int ibit=0; ibit<trigBits.size(); ibit++)
        {   
            if(trigBits[ibit].second != ttbar_chan) continue;
            hasTrigger |= ((ttbar_trigWord>>trigBits[ibit].first) & 1);
        }
    }else{
        hasTrigger=true;
    }

    return hasTrigger;

}

//Sources                                                                                                                                                          
// CMS AN 022/2015 v15                                                                                                                                     
// https://indico.cern.ch/event/434078/#preview:1614815  
pair<float,float> TTbarSelector::getTriggerEfficiency(int channel)
{
  pair<float,float>res(1.0,0.0);
  //if(channel == -11*13) { res.first=0.91; res.second=0.05; }
  //if(channel == -11*13) { res.first=1.0; res.second=0.05; }
  if(channel == -11*13) { res.first=0.901; res.second=0.015; } //  ICHEP2016 dataset
  //if(channel == -11*11) { res.first=0.95; res.second=0.05; }
  if(channel == -11*11) { res.first=1.0; res.second=0.05; }
  //if(channel == -13*13) { res.first=0.99; res.second=0.001; } // this is a single muon trigger
  if(channel == -13*13) { res.first=1.0; res.second=0.05; } // this is a single muon trigger

 return res;
}

//Sources
// CMS AN 022/2015 v15
pair<float,float> TTbarSelector::getLeptonSelectionEfficiencyScaleFactor(int id,float pt,float eta)
{
  pair<float,float>res(1.0,0.0);
 
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
vector<float> TTbarSelector::getJetResolutionScales(float pt, float eta, float genjpt)
{
  vector<float> res(3,1.0);

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


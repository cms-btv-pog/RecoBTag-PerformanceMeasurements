/**_________________________________________________________________
   class:   BTagLeptonEvent.cc
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagLeptonEvent.cc,v 1.2 2007/09/27 03:23:01 yumiceva Exp $

________________________________________________________________**/


#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"

ClassImp(BTagLeptonEvent)

// ROOT

//_______________________________________________________________
BTagLeptonEvent::BTagLeptonEvent() {

	this->Reset();
	
}

//_______________________________________________________________
BTagLeptonEvent::~BTagLeptonEvent() {
}

//_______________________________________________________________
void BTagLeptonEvent::Reset() {

	pdgid.clear();
	pt.clear();
    eta.clear();
    phi.clear();
	e.clear();
    charge.clear();
    //p.clear();
    trkchi2.clear();
    trkndof.clear();
	chi2.clear();
	ndof.clear();
	SArechits.clear();
    trkrechits.clear();
	d0.clear();
	d0sigma.clear();
	//vx.clear();
	//vy.clear();
	//vz.clear();
	jet_deltaR.clear();
	jet_ptrel.clear();
	
	mc_pt.clear();
	mc_eta.clear();
	mc_phi.clear();
	mc_charge.clear();
	//mc_p.clear();
	mc_pdgid.clear();
	//mc_vx.clear();
	//mc_vy.clear();
	//mc_vz.clear();
    mc_mother_pdgid.clear();

}

/**_________________________________________________________________
   class:   BTagEvent.cc
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagEvent.cc,v 1.1 2007/09/24 18:26:47 yumiceva Exp $

________________________________________________________________**/


#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"

ClassImp(BTagEvent)

// ROOT

//_______________________________________________________________
BTagEvent::BTagEvent() {

	this->Reset();
	
}

//_______________________________________________________________
BTagEvent::~BTagEvent() {
}

//_______________________________________________________________
void BTagEvent::Reset() {

	event   = -1;
	run     = -1;

		
	// mc

	// reco
	njets       = -1;
	nmuons      = -1;
	
	nvertices     = -1;
	ngenjets    = -1;

	jet_p.clear();
	jet_pt.clear();
	jet_eta.clear();
	jet_phi.clear();
	jet_e.clear();
	jet_et.clear();
	jet_ntrks.clear();
	jet_vx.clear();
	jet_vy.clear();
	jet_vz.clear();
	jet_flavour_alg.clear();
	jet_flavour_phy.clear();
	jet_isbtagged.clear();
	jet_hasLepton.clear();
	
	jetcorrection.clear();

	genjet_p.clear();
	genjet_pt.clear();
	genjet_eta.clear();
	genjet_phi.clear();
	genjet_et.clear();
	genjet_vx.clear();
	genjet_vy.clear();
	genjet_vz.clear();

	btag_TrkCounting_disc2D_1trk.clear();
	btag_TrkCounting_disc2D_2trk.clear();
	btag_TrkCounting_disc2D_3trk.clear();
	btag_TrkCounting_disc3D_1trk.clear();
	btag_TrkCounting_disc3D_2trk.clear();
	btag_TrkCounting_disc3D_3trk.clear();

	btag_JetProb_disc2D.clear();
	btag_JetProb_disc3D.clear();

	btag_NegTag_disc2D_1trk.clear();
	btag_NegTag_disc2D_2trk.clear();
	btag_NegTag_disc2D_3trk.clear();
	btag_NegTag_disc3D_1trk.clear();
	btag_NegTag_disc3D_2trk.clear();
	btag_NegTag_disc3D_3trk.clear();

	lepton.clear();
	
}



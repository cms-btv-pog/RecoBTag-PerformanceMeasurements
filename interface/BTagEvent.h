#ifndef BTagEvent_h
#define BTagEvent_h

/**_________________________________________________________________
   class:   BTagEvent.h
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagEvent.h,v 1.2 2007/09/27 03:23:00 yumiceva Exp $

________________________________________________________________**/


// C++ standard

// ROOT
#include "TObject.h"
#include "TMatrixDSym.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"


class BTagEvent : public TObject {

  public:

	BTagEvent();
	~BTagEvent();

	void Reset();
	
	Int_t event; // event number
	Int_t run;   // run number

		
	// mc

	// reco
	Int_t njets;    // number of jets
	Int_t nmuons;   // number of muons
	Int_t nvertices;  // number of vertices
	Int_t ngenjets; // number of generated jets

	//std::vector< float > jet_p;
	std::vector< float > jet_pt;
	std::vector< float > jet_eta;
	std::vector< float > jet_phi;
	std::vector< float > jet_e;
	std::vector< float > jet_et;
	//std::vector< float > jet_vx;
	//std::vector< float > jet_vy;
	//std::vector< float > jet_vz;
	std::vector< int > jet_ntrks;

	std::vector< int > jet_flavour_phy;
	std::vector< int > jet_flavour_alg;
	//std::vector< int > jet_isbtagged;
	std::vector< int > jet_hasLepton;
	std::vector< float > jetcorrection;
	
	//std::vector< float > genjet_p;
	std::vector< float > genjet_pt;
	std::vector< float > genjet_eta;
	std::vector< float > genjet_phi;
	std::vector< float > genjet_et;
	//std::vector< float > genjet_vx;
	//std::vector< float > genjet_vy;
	//std::vector< float > genjet_vz;

	std::vector< float > btag_TrkCounting_disc2D_1trk;
	std::vector< float > btag_TrkCounting_disc2D_2trk;
	std::vector< float > btag_TrkCounting_disc2D_3trk;
	std::vector< float > btag_TrkCounting_disc3D_1trk;
	std::vector< float > btag_TrkCounting_disc3D_2trk;
	std::vector< float > btag_TrkCounting_disc3D_3trk;

	//std::vector< float > btag_JetProb_disc2D;
	std::vector< float > btag_JetProb_disc3D;

	//std::vector< float > btag_NegTag_disc2D_1trk;
	//std::vector< float > btag_NegTag_disc2D_2trk;
	//std::vector< float > btag_NegTag_disc2D_3trk;
	std::vector< float > btag_NegTag_disc3D_1trk;
	std::vector< float > btag_NegTag_disc3D_2trk;
	std::vector< float > btag_NegTag_disc3D_3trk;
		
	std::vector< BTagLeptonEvent > lepton;
	
	ClassDef(BTagEvent,1);

};

#endif

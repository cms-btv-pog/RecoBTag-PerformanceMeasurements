#ifndef BTagEvent_h
#define BTagEvent_h

/**_________________________________________________________________
   class:   BTagEvent.h
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagEvent.h,v 1.1 2007/09/24 18:26:47 yumiceva Exp $

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

	std::vector< double > jet_p;
	std::vector< double > jet_pt;
	std::vector< double > jet_eta;
	std::vector< double > jet_phi;
	std::vector< double > jet_e;
	std::vector< double > jet_et;
	std::vector< double > jet_vx;
	std::vector< double > jet_vy;
	std::vector< double > jet_vz;
	std::vector< int > jet_ntrks;

	std::vector< int > jet_flavour_phy;
	std::vector< int > jet_flavour_alg;
	std::vector< int > jet_isbtagged;
	std::vector< int > jet_hasLepton;
	std::vector< double > jetcorrection;
	
	std::vector< double > genjet_p;
	std::vector< double > genjet_pt;
	std::vector< double > genjet_eta;
	std::vector< double > genjet_phi;
	std::vector< double > genjet_et;
	std::vector< double > genjet_vx;
	std::vector< double > genjet_vy;
	std::vector< double > genjet_vz;

	std::vector< double > btag_TrkCounting_disc2D_1trk;
	std::vector< double > btag_TrkCounting_disc2D_2trk;
	std::vector< double > btag_TrkCounting_disc2D_3trk;
	std::vector< double > btag_TrkCounting_disc3D_1trk;
	std::vector< double > btag_TrkCounting_disc3D_2trk;
	std::vector< double > btag_TrkCounting_disc3D_3trk;

	std::vector< double > btag_JetProb_disc2D;
	std::vector< double > btag_JetProb_disc3D;

	std::vector< double > btag_NegTag_disc2D_1trk;
	std::vector< double > btag_NegTag_disc2D_2trk;
	std::vector< double > btag_NegTag_disc2D_3trk;
	std::vector< double > btag_NegTag_disc3D_1trk;
	std::vector< double > btag_NegTag_disc3D_2trk;
	std::vector< double > btag_NegTag_disc3D_3trk;
		
	std::vector< BTagLeptonEvent > lepton;
	
	ClassDef(BTagEvent,1);

};

#endif

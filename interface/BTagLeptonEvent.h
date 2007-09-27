#ifndef BTagLeptonEvent_h
#define BTagLeptonEvent_h

/**_________________________________________________________________
   class:   BTagLeptonEvent.h
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagLeptonEvent.h,v 1.1 2007/09/24 18:26:47 yumiceva Exp $

________________________________________________________________**/

#include "TObject.h"

class BTagLeptonEvent : public TObject {

  public:

	BTagLeptonEvent();
	~BTagLeptonEvent();

	void Reset();

	std::vector< int > pdgid;
	std::vector< double > pt;
	std::vector< double > eta;
	std::vector< double > phi;
	std::vector< double > e;
	std::vector< double > charge;
	std::vector< double > p;
	std::vector< double > trkchi2;
	std::vector< double > trkndof;
	std::vector< double > chi2;
	std::vector< double > ndof;
	std::vector< int > SArechits;
	std::vector< int > trkrechits;
	std::vector< double > d0;
	std::vector< double > d0sigma;
	std::vector< double > vx;
	std::vector< double > vy;
	std::vector< double > vz;
	std::vector< double > jet_deltaR;
	std::vector< double > jet_ptrel;

	std::vector< double > mc_pt;
	std::vector< double > mc_eta;
	std::vector< double > mc_phi;
	std::vector< double > mc_charge;
	std::vector< double > mc_p;
	std::vector< int > mc_pdgid;
	std::vector< double > mc_vx;
	std::vector< double > mc_vy;
	std::vector< double > mc_vz;
	std::vector< int > mc_mother_pdgid;

	ClassDef(BTagLeptonEvent,1);
};

#endif

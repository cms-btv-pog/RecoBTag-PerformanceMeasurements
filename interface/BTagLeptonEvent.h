#ifndef BTagLeptonEvent_h
#define BTagLeptonEvent_h

/**_________________________________________________________________
   class:   BTagLeptonEvent.h
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagLeptonEvent.h,v 1.2 2007/09/27 03:23:00 yumiceva Exp $

________________________________________________________________**/

#include "TObject.h"

class BTagLeptonEvent : public TObject {

  public:

	BTagLeptonEvent();
	~BTagLeptonEvent();

	void Reset();

	std::vector< int > pdgid;
	std::vector< float > pt;
	std::vector< float > eta;
	std::vector< float > phi;
	std::vector< float > e;
	std::vector< int > charge;
	//std::vector< float > p;
	std::vector< float > trkchi2;
	std::vector< float > trkndof;
	std::vector< float > chi2;
	std::vector< float > ndof;
	std::vector< int > SArechits;
	std::vector< int > trkrechits;
	std::vector< float > d0;
	std::vector< float > d0sigma;
	//std::vector< float > vx;
	//std::vector< float > vy;
	//std::vector< float > vz;
	std::vector< float > jet_deltaR;
	std::vector< float > jet_ptrel;

	std::vector< float > mc_pt;
	std::vector< float > mc_eta;
	std::vector< float > mc_phi;
	std::vector< float > mc_charge;
	//std::vector< float > mc_p;
	std::vector< int > mc_pdgid;
	//std::vector< float > mc_vx;
	//std::vector< float > mc_vy;
	//std::vector< float > mc_vz;
	std::vector< int > mc_mother_pdgid;

	ClassDef(BTagLeptonEvent,1);
};

#endif

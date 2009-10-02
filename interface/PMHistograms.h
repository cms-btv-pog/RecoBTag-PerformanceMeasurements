//
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      PM_Histograms
//
/**\class PerformanceMeasurements/PM_Histograms

 Description:

         Author: Francisco Yumiceva, Fermilab
*/
//
// $Id: PM_Histogrmas,v 1.1 2009/09/22 02:45:37 yumiceva Exp $
//
//


#ifndef PMHistograms_H
#define PMHistograms_H

#include "PhysicsTools/FWLite/interface/TH1Store.h"
#include "TH2F.h"
#include "TString.h"
#include "TLorentzVector.h"

class PMHistograms
{

  public:

		//PMHistograms() {}
	PMHistograms( TH1Store *hstore ) {
		fstore = hstore;
	}
	~PMHistograms() { }

	void Add();
	void FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
					int JetFlavor, bool tagged=false);
	
  private:

	TH1Store *fstore;

};

#endif


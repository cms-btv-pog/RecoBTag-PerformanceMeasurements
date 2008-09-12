/**_________________________________________________________________
   class:   BTagLeptonEvent.cc
   package: RecoBTag/PerformanceMeasurements


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagLeptonEvent.cc,v 1.6 2008/08/11 06:06:22 bazterra Exp $

________________________________________________________________**/

#ifdef NOSCRAMV
#include "BTagLeptonEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"
#endif

ClassImp(BTagLeptonEvent)

// ROOT

//_______________________________________________________________
void BTagLeptonEvent::Reset()
{
    e.clear();
    pdgid.clear();
    chi2.clear();
    ndof.clear();
    SArechits.clear();
    //vx.clear();
    //vy.clear();
    //vz.clear();

    mc_pt.clear();
    mc_eta.clear();
    mc_phi.clear();
    mc_charge.clear();
    mc_e.clear();
    mc_pdgid.clear();
    //mc_vx.clear();
    //mc_vy.clear();
    //mc_vz.clear();
    mc_mother_pdgid.clear();
}

#ifndef BTagLeptonEvent_h
#define BTagLeptonEvent_h

/**_________________________________________________________________
   class:   BTagLeptonEvent.h
   package: RecoBTag/PerformanceMeasurements


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagLeptonEvent.h,v 1.6 2008/08/11 06:06:21 bazterra Exp $

________________________________________________________________**/

#ifdef NOSCRAMV
#include "BTagBaseTrackEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagBaseTrackEvent.h"
#endif

class BTagLeptonEvent : public BTagBaseTrackEvent
{

public:

    BTagLeptonEvent() : BTagBaseTrackEvent()
    {
        Reset();
    }
    ~BTagLeptonEvent() {}

    virtual void Reset();

    std::vector< float > e;
    std::vector< float > chi2;
    std::vector< float > ndof;
    std::vector< int > SArechits;
    //std::vector< float > vx;
    //std::vector< float > vy;
    //std::vector< float > vz;

    std::vector< int > pdgid;
    std::vector< float > mc_pt;
    std::vector< float > mc_eta;
    std::vector< float > mc_phi;
    std::vector< float > mc_charge;
    std::vector< float > mc_e;
    std::vector< int > mc_pdgid;
    //std::vector< float > mc_vx;
    //std::vector< float > mc_vy;
    //std::vector< float > mc_vz;
    std::vector< int > mc_mother_pdgid;

    ClassDef(BTagLeptonEvent,1);
};

#endif

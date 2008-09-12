#ifndef BTagTrackEvent_h
#define BTagTrackEvent_h

/**_________________________________________________________________
   class:   BTagTrackEvent.h
   package: RecoBTag/PerformanceMeasurements


 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.h,v 1.7 2008/09/09 17:04:50 bazterra Exp $
________________________________________________________________**/

#ifdef NOSCRAMV
#include "BTagBaseTrackEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagBaseTrackEvent.h"
#endif

class BTagTrackEvent : public BTagBaseTrackEvent
{

public:

    BTagTrackEvent() : BTagBaseTrackEvent()
    {
        Reset();
    }
    ~BTagTrackEvent() {}

    virtual void Reset();

    std::vector<float> ip2d, ip3d, dta;
    std::vector<float> ip2dSigma, ip3dSigma;

    ClassDef(BTagTrackEvent,1);
};

#endif

/**_________________________________________________________________
   class:   BTagTrackEvent.cc
   package: RecoBTag/PerformanceMeasurements


 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.cc,v 1.5 2008/09/08 19:06:24 bazterra Exp $

________________________________________________________________**/

#ifdef NOSCRAMV
#include "BTagTrackEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagTrackEvent.h"
#endif

ClassImp(BTagTrackEvent)

// ROOT

//_______________________________________________________________
void BTagTrackEvent::Reset()
{
    ip2d.clear();
    ip3d.clear();
    dta.clear();
    ip2dSigma.clear();
    ip3dSigma.clear();
}

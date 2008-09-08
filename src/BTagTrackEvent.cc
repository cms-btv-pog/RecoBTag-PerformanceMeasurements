/**_________________________________________________________________
   class:   BTagTrackEvent.cc
   package: RecoBTag/PerformanceMeasurements


 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.cc,v 1.4 2008/08/27 20:13:14 bazterra Exp $

________________________________________________________________**/


#include "RecoBTag/PerformanceMeasurements/interface/BTagTrackEvent.h"

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

/**_________________________________________________________________
   class:   BTagTrackEvent.cc
   package: RecoBTag/PerformanceMeasurements


 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.cc,v 1.3 2008/08/11 06:06:22 bazterra Exp $

________________________________________________________________**/


#include "RecoBTag/PerformanceMeasurements/interface/BTagTrackEvent.h"

ClassImp(BTagTrackEvent)

// ROOT

//_______________________________________________________________
void BTagTrackEvent::Reset()
{
    //	is.clear();
    ip2d.clear();
    ip3d.clear();
    dta.clear();
    ip2dSigma.clear();
    ip3dSigma.clear();
    is.clear();
}

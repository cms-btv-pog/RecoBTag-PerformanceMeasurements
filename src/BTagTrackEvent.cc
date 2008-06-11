/**_________________________________________________________________
   class:   BTagTrackEvent.cc
   package: RecoBTag/PerformanceMeasurements
   

 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.cc,v 1.1 2008/03/11 21:09:39 bazterra Exp $

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
	sdl.clear();
	dta.clear();
	ip2dSigma.clear();
	ip3dSigma.clear();
	sdlSigma.clear();
	dtaSigma.clear();
}

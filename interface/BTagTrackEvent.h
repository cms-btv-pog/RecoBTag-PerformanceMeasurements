#ifndef BTagTrackEvent_h
#define BTagTrackEvent_h

/**_________________________________________________________________
   class:   BTagTrackEvent.h
   package: RecoBTag/PerformanceMeasurements
   

 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.h,v 1.2 2008/03/17 02:23:03 bazterra Exp $
________________________________________________________________**/

#include "RecoBTag/PerformanceMeasurements/interface/BTagBaseTrackEvent.h"

class BTagTrackEvent : public BTagBaseTrackEvent {

 public:

	BTagTrackEvent() : BTagBaseTrackEvent() { Reset(); }
	~BTagTrackEvent() {} 

	virtual void Reset();

	std::vector<float> ip2d, ip3d, sdl, dta;
	std::vector<float> ip2dSigma, ip3dSigma, sdlSigma, dtaSigma;

	ClassDef(BTagTrackEvent,1);
};

#endif

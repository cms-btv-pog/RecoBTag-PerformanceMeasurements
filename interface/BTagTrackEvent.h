#ifndef BTagTrackEvent_h
#define BTagTrackEvent_h

/**_________________________________________________________________
   class:   BTagTrackEvent.h
   package: RecoBTag/PerformanceMeasurements
   

 author: Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagTrackEvent.h,v 1.4 2007/09/28 23:13:14 bazterra Exp $
________________________________________________________________**/

#include "RecoBTag/PerformanceMeasurements/interface/BTagBaseTrackEvent.h"

class BTagTrackEvent : public BTagBaseTrackEvent {

 public:

	BTagTrackEvent() : BTagBaseTrackEvent() { Reset(); }
	~BTagTrackEvent() {} 

	virtual void Reset();

	std::vector<float> ip2d, ip3d, sdl, dta;
	std::vector<float> ip2dSigma, ip3dSigma, sdlSigma, dtaSigma;

	enum Category {
		Fake = 0,
		Bad,
		SignalEvent,
 		PV,
 		SV,
 		TV,
 		Displaced,
 		Ks,
 		Lambda,
 		PhotonConversion,
 		Up,
 		Down,
 		Strange,
 		Charm,
 		Bottom,
 		Light,
 		Unknown
  	};

    typedef std::vector<bool> Flags;

	std::vector<Flags> is;

	ClassDef(BTagTrackEvent,1);
};

#endif

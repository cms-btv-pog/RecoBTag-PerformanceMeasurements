
#ifdef NOSCRAMV
#include "BTagEvent.h"
#include "BTagTrackEvent.h"
#include "BTagLeptonEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagTrackEvent.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"
#endif
 
#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
 
#pragma link C++ class BTagEvent+;
#pragma link C++ class BTagBaseTrackEvent+;
#pragma link C++ class BTagTrackEvent+;
#pragma link C++ class BTagLeptonEvent+;

#endif

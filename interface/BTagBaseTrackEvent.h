#ifndef BTagBaseTrackEvent_h
#define BTagBaseTrackEvent_h

/**_______________________________________________________________________________
   class:   BTagBaseTrackEvent.h
   package: RecoBTag/PerformanceMeasurements


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
         Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagBaseTrackEvent.h,v 1.3 2008/06/11 08:46:58 tboccali Exp $
________________________________________________________________________________**/

#include <vector>

#include "TObject.h"

// This class
class BTagBaseTrackEvent : public TObject
{

public:

    BTagBaseTrackEvent()
    {
        Reset();
    }
    ~BTagBaseTrackEvent() {}

    virtual void Reset();

    std::vector< float > pt;
    std::vector< float > eta;
    std::vector< float > phi;
    std::vector< int > charge;
    std::vector< float > trkchi2;
    std::vector< float > trkndof;
    std::vector< int > trkrechits;
    std::vector< float > d0;
    std::vector< float > d0sigma;

    std::vector< float > jet_deltaR;
    std::vector< float > jet_ptrel;

    ClassDef(BTagBaseTrackEvent,1);
};

#endif

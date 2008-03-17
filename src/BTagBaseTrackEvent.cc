/**_________________________________________________________________
   class:   BTagBaseTrackEvent.cc
   package: RecoBTag/PerformanceMeasurements
   

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
         Victor E. Bazterra, UIC (baites@fnal.gov)

 version $Id: BTagBaseTrackEvent.cc,v 1.1 2008/03/11 21:09:39 bazterra Exp $

________________________________________________________________**/

#include "RecoBTag/PerformanceMeasurements/interface/BTagBaseTrackEvent.h"

ClassImp(BTagBaseTrackEvent)

//_______________________________________________________________
void BTagBaseTrackEvent::Reset() 
{
    pt.clear();
    eta.clear();
    phi.clear();
    charge.clear();
    trkchi2.clear();
    trkndof.clear();
    trkrechits.clear();
    d0.clear();
    d0sigma.clear();
    jet_deltaR.clear();
    jet_ptrel.clear();	
}

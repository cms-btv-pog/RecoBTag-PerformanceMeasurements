#include "RecoBTag/PerformanceMeasurements/interface/TopJetSelection.h"

#include<iostream>

using namespace std;

TopJetSelection::TopJetSelection()
{
  // cut parameters
  etaMin     = 0.0;
  etaMax     = 2.4;

  etMin = 0.0;
  etMax = 9999.9;
}
// 'global' event selection based on basic variables


bool TopJetSelection::operator() (const TopJet & jet) const
{

  bool accept = true;

  // temporary fudge to correct for double loop error
  //  jetPartonMomentum /= 2.0;

  if ( fabs(jet.eta()) < etaMin  ||
       fabs(jet.eta()) > etaMax  ) return false;

  if ( jet.et() < etMin ||
       jet.et() > etMax ) return false;

  return accept;
}

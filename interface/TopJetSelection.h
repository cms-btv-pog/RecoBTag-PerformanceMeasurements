#ifndef TopJetSelection_H
#define TopJetSelection_H

#include "AnalysisDataFormats/TopObjects/interface/TopJet.h" 

/** \class TopJetSelection
 *
 *  Decide if jet and associated parton satisfy desired kinematic cuts.
 *
 */

class TopJetSelection {

 public:
  TopJetSelection();
  /// Returns true if jet and associated parton satisfy kinematic cuts.
  bool operator() (const TopJet & jet) const;

  /// Set cut parameters
  void setEtaMin ( double d ) { etaMin  = d ; } 
  void setEtaMax ( double d ) { etaMax  = d ; } 
  void setEtMin  ( double d ) { etMin = d ; } 
  void setEtMax  ( double d ) { etMax = d ; } 

 protected:

  // eta range 
  double etaMin ;   // these are meant as |eta| !!
  double etaMax ;

  // et range 
  double etMin ;   // these are meant as |eta| !!
  double etMax ;

  
  
} ;

#endif

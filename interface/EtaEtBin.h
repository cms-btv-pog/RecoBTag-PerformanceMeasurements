#ifndef EtaEtBin_H
#define EtaEtBin_H

#include "TString.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

/** \class EtaEtBin
 *
 *  Decide if jet/parton lie within desired rapidity/et range.
 *
 */

class EtaEtBin {

 public:

  EtaEtBin(bool etaActive_ , double etaMin_ , double etaMax_ ,
	    bool etActive_  , double etMin_  , double etMax_ ) ;
  EtaEtBin ( const EtaEtBin & other);

  ~EtaEtBin () {} ;

  /// String describes rapidity/et range.
  TString getDescriptionString () const;

  /// Get rapidity/et ranges and check whether rapidity/et cuts are active.
  bool   getEtaActive () const { return etaActive ; }
  double getEtaMin    () const { return etaMin    ; }
  double getEtaMax    () const { return etaMax    ; }

  bool   getEtActive () const { return etActive ; }
  double getEtMin    () const { return etMin    ; }
  double getEtMax    () const { return etMax    ; }


  /// Check if jet/parton are within rapidity/et cuts.
  bool inBin(const double & eta , const double & et) const;
  bool inBin(const pat::Jet & jet) const;

  bool operator==(const EtaEtBin & other) const;
  EtaEtBin* clone() const;


 private:

  // definition of the bin

  bool   etaActive ; // should cuts be applied?
  double etaMin ;
  double etaMax ;
  
  bool   etActive ; // should cuts be applied?
  double etMin ;
  double etMax ;
  

  // description string as built from bin definition
  
  
} ;


#endif

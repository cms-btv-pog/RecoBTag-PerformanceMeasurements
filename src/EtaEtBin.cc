#include "RecoBTag/PerformanceMeasurements/interface/EtaEtBin.h"

EtaEtBin::EtaEtBin ( bool etaActive_ , double etaMin_ , double etaMax_ ,
		     bool etActive_  , double etMin_  , double etMax_ )
  : etaActive ( etaActive_ ) , etaMin ( etaMin_ ) , etaMax ( etaMax_ ) ,
    etActive  (  etActive_ ) , etMin  (  etMin_ ) , etMax  (  etMax_ )   {
}

bool EtaEtBin::operator==(const EtaEtBin & other) const
{
  return (other.getDescriptionString()==this->getDescriptionString());
}

EtaEtBin::EtaEtBin ( const EtaEtBin & other)
{
  etaActive = other.getEtaActive ();
  etaMin = other.getEtaMin ();
  etaMax = other.getEtaMax ();
 
  etActive  = other.getEtActive ();
  etMin  = other.getEtMin ();
  etMax  = other.getEtMax ();
}

EtaEtBin* EtaEtBin::clone() const
{
  return new EtaEtBin(*this);
}



TString EtaEtBin::getDescriptionString() const
{
  // create string only from the active parts
  TString descr ( "" );

  if ( etaActive ) {
    descr += "_ETA_";
    descr += etaMin;
    descr += "-";
    descr += etaMax;
  }

  if ( etActive ) {
    descr += "_ET_";
    descr += etMin;
    descr += "-";
    descr += etMax;
  }
  if (!(etaActive||etActive)) descr="_GLOBAL";
  // remove blanks which are introduced when adding doubles
  descr.ReplaceAll ( " " , "" );

  return descr;
}

bool EtaEtBin::inBin(const pat::Jet & jet) const
{
  return inBin(jet.eta(), jet.et());
}

bool EtaEtBin::inBin (const double & eta , const double & et ) const {
  bool inEta = true;
  //
  if ( etaActive ) {
    if ( fabs(eta) < etaMin ) inEta = false;
    if ( fabs(eta) > etaMax ) inEta = false;
  }

  bool inEt = true;
  //
  if ( etActive ) {
    if ( et < etMin ) inEt = false;
    if ( et > etMax ) inEta = false;
  }

  return ( inEta && inEt );
}

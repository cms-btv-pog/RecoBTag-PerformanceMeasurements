#ifndef IPSPlotter_h
#define IPSPlotter_h

#include "TH1F.h"

#include "RecoBTag/PerformanceMeasurements/test/Plotters/BasePlotter.h"

class IPSPlotter : public BasePlotter {

 public:

  IPSPlotter(TString filename="") : BasePlotter(filename) {}

  void Write();

  ClassDef(IPSPlotter,1);

 protected:

  void Book(); 
  void Fill(BTagEvent*);

 private:

  TH1F * ips_; 
  TH1F * ipsB_;
  TH1F * ipsC_;
  TH1F * ipsKs_;
  TH1F * ipsElse_;
  TH1F * ipsLight_;
  TH1F * ipsLambda_;
  TH1F * ipsDisplaced_;
  TH1F * ipsConversion_;
  
};

#endif

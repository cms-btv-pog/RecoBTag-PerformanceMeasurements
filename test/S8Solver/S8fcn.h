#ifndef S8fcn_h
#define S8fcn_h

#include "Minuit2/FCNBase.h" 

using namespace ROOT::Minuit2;

class S8fcn : public FCNBase {

 public:
  // cache the current data

  void SetData(std::vector < double> data){    
    fdata = data;
  };
  void SetCorr(double kb, double beta, double kcl, double alpha, double delta, double gamma) {
    fkb = kb;
    fbeta = beta;
    fkcl = kcl;
    falpha = alpha;
	fdelta = delta;
	fgamma = gamma;
  }
  
  virtual double operator() (const std::vector<double>&) const;
  virtual double Up() const {return 1.;}
  
 private:

  std::vector < double> fdata;
  double fkb;
  double fbeta;
  double fkcl;
  double falpha;
  double fdelta;
  double fgamma;
};

#endif



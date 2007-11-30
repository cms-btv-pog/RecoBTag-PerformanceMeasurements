#ifndef S8FitSolver_h
#define S8FitSolver_h

#include "S8fcn.h"

#include "Minuit2/VariableMetricMinimizer.h"
#include "TString.h"
#include<map>

using namespace ROOT::Minuit2;

class S8FitSolver {

  public:
	S8FitSolver();
	void Solve(std::map<TString,double> input);
	std::map<TString,double> GetSolution(){ return fsolution;};
	std::map<TString,double> GetSolutionErr(){ return fsolutionErr;};
	void Verbose(bool option) { fVerbose = option; };
	void Init(std::map< TString, double > input);
	
  private:
	bool fVerbose;
	ModularFunctionMinimizer* theFitter;
	S8fcn *thePDF;
	double ff_minimum;
	std::map<TString,double> fsolution;
	std::map<TString,double> fsolutionErr;
	std::map< TString, double > finitmap;
	double falpha;
	double fbeta;
	double fkb;
	double fkcl;
	//double data[8];
        //ClassDef(S8FitSolver,1);
};

#endif

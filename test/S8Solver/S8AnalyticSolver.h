#ifndef S8AnalyticSolver_h
#define S8AnalyticSolver_h

#include "TString.h"
#include<map>

class S8AnalyticSolver {

  public:
	S8AnalyticSolver();
	//~S8AnalyticSolver();
	void Solve(std::map<TString,double> input);
	std::map<TString,double> GetSolution(){ return fsolution;};
	std::map<TString,double> GetSolutionErr(){ return fsolutionErr;};
	void Verbose(bool option) { fVerbose = option; };
	double sqr(double x) {
		return x*x;
	};

  private:
	bool fVerbose;
	std::map<TString,double> fsolution;
	std::map<TString,double> fsolutionErr;

	//ClassDef(S8AnalyticSolver,1);
};

#endif

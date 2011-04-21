////////////////////////////////////////////////
//  System 8 solver                           //
//  Author : Benoit Clement                   //
//           benoit.clement@ires.in2p3.fr     //
//  Date : 01 Dec 2003                        //
//  modified P.Van Hove	2011                  //
////////////////////////////////////////////////

#ifndef System8Solver_h
#define System8Solver_h

#include "TNamed.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include <string>
#include <vector>
#include <map>
#include <Math/RootFinder.h>
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

class TH1;
class NumericInputGroup;

class System8
{
public:
  System8() ;
  ~System8() {;}
  int Solve();

  void   SetInput    (double, double, double, double, double, double, double);
  void   SetCorrelation(double, double, double, double, double, double, double, double);  
  int    GetNbSolutions()    {return fSolutions.size();}
  double GetSolution(int idx, int i, int j)  {fIndex=idx; return UKN(i,j); }  
  void   SetNpt(int n)      {fNpt=n;}

  TGraph* MakeGraph(int idx, int npt=1000);
  
private:  
  void SetResult(double v)     {fFa=v;}

  int fIndex;
  int fNpt;  
  double fq[3], fQ[3], fQQQ ;
  double kc[2][4] ;
  
  double UKN(int i,int j); // system unknowns   

  // Internal system computing. See corresponding ps file for details
  std::vector<std::pair<int, double> > fSolutions;
  
  double fFa  ; // "a" sample fraction 
  double fSign; // "a" sample fraction 
  long double fScale[2]; // Scaling to avoid some numerical problems
           
  long double Q(int i, int j) {if(j==0) return fq[i]; if(j==1) return fQ[i]; return fQQQ;}

  long double U()  {return fFa     ;}
  long double V()  {return U()-1.  ;}

  long double A(int i) {return  U()*Q(i,1)-kc[0][i]*Q(i,0)*Q((i+1)%3,0);}
  long double B(int i) {return -V()*kc[0][i]*Q(i,0)                    ;}
  long double C(int i) {return  V()*kc[0][i]*Q((i+1)%3,0)              ;}
  long double D(int i) {return  V()*(V()*kc[0][i]-U()*kc[1][i])       ;}

  long double M() {return C(0)*C(1)*D(2)+D(0)*A(1)*D(2)+C(0)*D(1)*B(2)+D(0)*B(1)*B(2) ;}
  long double N() {return C(0)*C(1)*C(2)+D(0)*A(1)*C(2)+C(0)*D(1)*A(2)+D(0)*B(1)*A(2)-A(0)*C(1)*D(2)-B(0)*A(1)*D(2)-A(0)*D(1)*B(2)-B(0)*B(1)*B(2) ;}
  long double P() {return A(0)*C(1)*C(2)+B(0)*A(1)*C(2)+A(0)*D(1)*A(2)+B(0)*B(1)*A(2) ;}
  long double DELTA() {return N()*N()+4*M()*P();}
  long double R()      {return kc[0][3]*Q(0,0)*Q(1,0)*Q(2,0) - Q(0,2)*U()*U();}
  long double S(int i) {return kc[0][3]*V()*Q((i+1)%3,0)*Q((i+2)%3,0)        ;}
  long double T(int i) {return kc[0][3]*V()*V()*Q((i+2)%3,0)                 ;}
  long double W()      {return V()*(kc[0][3]*V()*V() - kc[1][3]*U()*U())     ;}   
   
  long double Gamma(int i);
  void EvalScale();
  long double ZeroFunction(double x) ;// for wrapper
};

//////////////////////////////////////////////////////////////////////////////////////////////
class SolutionSelector
{
public: 
  SolutionSelector() ;
  ~SolutionSelector() ;
  void SetInitialOrder(Int_t n, Int_t a = 1);
  void SetInitialOrder(std::string s1, std::string s2, std::string s3);
  void ResetOrder();
  double GetSolution(int i, int j)      {return fSolution[i*4+j];} 
  double GetSolution(int i)             {return fSolution[i];}   
  double* GetSolution()  {return fSolution;}
  bool FindSolution(System8*);
 
private:
  int    fOrder[4]    ;
  double fSolution[8];
};

//////////////////////////////////////////////////////////////////////////////////////////////
class System8Solver
{
public :

  enum ErrorType  {NONE, SYST, STAT, ALL};
  System8Solver() ;
  System8Solver(double n  , double n1 , double n2 , double n3,
                double n12, double n23, double n31, double n123) ;
  System8Solver(double* n) ;
  ~System8Solver();

  void Reset();
    void setInput(const NumericInputGroup &);
  void SetInput(double n  , double n1 , double n2 , double n3,
                double n12, double n23, double n31, double n123) ;
  void SetInput(double* n) { SetInput(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]);}

  void SetCorr(double c12a, double c23a, double c31a,double c123a,
                double c12b, double c23b, double c31b,double c123b) ;
  void SetCorr(double* c) {SetCorr(c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);}

  void SetCorrError(double c12a, double c23a, double c31a, double c123a,
                    double c12b, double c23b, double c31b, double c123b) ;
  void SetCorrError(double* c) {SetCorrError(c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);}
  
  void SetCovMatrixSign(TMatrixDSym CovMatrixSign);
  void SetCovMatrixBack(TMatrixDSym CovMatrixBack);
  void SetCovMatrix(TMatrixDSym CovMatrixSign, TMatrixDSym CovMatrixBack) {SetCovMatrixSign(CovMatrixSign); SetCovMatrixBack(CovMatrixBack);} 
  void SetCovMatrix(TMatrixDSym CovMatrix);

  void SetBinNumber(Int_t i){fBin = i;}
  void SetError( ErrorType err, int niter ) {kError = err ; fIter = niter;} 
  void SetPrecision( Int_t pr ) {fNpt = pr ;}
  void SetNbErrorIteration( Int_t pr ) {fIter = pr ;}
  void SetInitialOrder(Int_t n, Int_t a = 1)    {fSelector.SetInitialOrder(n,a);}
  void SetInitialOrder(std::string s1, std::string s2, std::string s3)  {fSelector.SetInitialOrder(s1,s2,s3);}
  void ResetOrder(){fSelector.ResetOrder();}
  Int_t Solve() ;
  double GetResult  (Int_t n) ;
  double GetErrorInf(Int_t n, std::string opt="") ;
  double GetErrorSup(Int_t n, std::string opt="") ;
  double GetResult  (std::string label) ;
  double GetErrorInf(std::string label, std::string opt="") ;
  double GetErrorSup(std::string label, std::string opt="") ;
  Int_t  GetNmiss(){return (fNmiss_Syst+fNmiss_Stat);}
  Int_t  GetNmissStat(){if(kError==ALL || kError==STAT) return (fNmiss_Stat); return 0;}
  Int_t  GetNmissSyst(){if(kError==ALL || kError==SYST) return (fNmiss_Syst); return 0;}
  TMatrixDSym GetCovMatrix(){return fCovmatrix;}
  void DrawKappa();

    TH1 *result(const int &);
    double getError(const int &);
    double getCentralValue(const int &);

private :

  System8 fSystem;
  SolutionSelector fSelector;

  TRandom3* fRndm;
  double fCorr[2][4]     ;
  TMatrixDSym fCovmatrix; 
  TH1F* Kappa;

  double fInput[8]    ;
  double fIndep[8]    ;
  
  ErrorType   kError  ;

  double* fResult;
  double fErrorinf_Stat[8] ;
  double fErrorsup_Stat[8] ;
  double fErrorinf_Syst[8] ;
  double fErrorsup_Syst[8] ;

  Int_t fBin;
  Int_t fNpt;
  Int_t fIter;
  Int_t fNmiss_Syst;
  Int_t fNmiss_Stat;


  void ComputeErrors() ;
    void fitErrors();
  void MakeSystem(double *shift);
  int SolveSystem(double *res);

  double GetValue(double* tab, std::string label) ;

    private:
        TH1    *_result[8];
        double  _result_errors[8];
        double _central_values[8];
};

void demo(Double_t* inputs, Double_t* corr, Int_t k, Double_t *ans);

#endif

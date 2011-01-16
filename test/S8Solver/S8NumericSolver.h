////////////////////////////////////////////////
//  System 8 solver                           //
//  Author : Benoit Clement                   //
//           benoit.clement@ires.in2p3.fr     //
//  Date : 01 Dec 2003                        //
// modified D.Gele 2007                       //
////////////////////////////////////////////////

#ifndef S8NumericSolver_S8NumericSolver_h
#define S8NumericSolver_S8NumericSolver_h


#include "TNamed.h"
#include "TRandom3.h"

#include <string>

// some root includes
#include <utility>

// user include files
#include <iostream>
#include <math.h>
#include <vector>
#include <map>

class TH1;
class NumericInput;

class S8NumericSolver : public TNamed
{
    private :
        TRandom3* fRndm;
        Double_t fCorr[2][3]     ;
        Double_t fCorrerr[2][3]  ;

        Double_t kc[2][3] ;

        Double_t fInput[8]    ;
        Double_t fIndep[8]    ;
        Double_t q[3], Q[3], QQQ ;

        Int_t    fAsym[2]    ;
        Int_t    kError      ;

        Double_t fResult[8];
        std::map< int, double > fMapResult;
        std::map< int, double > fMapErrorInf_Stat;
        std::map< int, double > fMapErrorSup_Stat;

        bool fForceSol;
        Int_t fpickSol;

        Double_t fNb;
        Double_t fErrorinf_Stat[8] ;
        Double_t fErrorsup_Stat[8] ;
        Double_t fErrorinf_Syst[8] ;
        Double_t fErrorsup_Syst[8] ;

        Int_t fNpt;
        Int_t fIter;

        double fAveRes;
        bool fverbose;
        bool fAveResSetup;

        void ComputeErrors() ;
        void fitErrors();

        void MakeSystem(Double_t *shift);
        int SolveSystem(Double_t *res);
        Int_t FindSolution(Double_t* res, int n);
        Double_t GetValue(Double_t* tab, std::string label) ;

        Double_t E(int i,int j);

        Double_t U(int i = 0) {return q[i];}
        Double_t V(int i = 0) {i=0; return W()-1.;}
        Double_t W(int i = 0) {i=0; return fNb;}

        Double_t A(int i = 2) {return W()*Q[i]-kc[1][i]*U(i)*U((i+1)%3) 	;}
        Double_t B(int i = 2) {return -kc[1][i]*U(i)*V() 			;}
        Double_t C(int i = 2) {return V()*U((i+1)%3)*kc[1][i]       	;}
        Double_t D(int i = 2) {return V()*(V()*kc[1][i]-W()*kc[0][i])	;}

        Double_t G() {return A(1)*C(2)+B(1)*A(2) ;}
        Double_t H() {return A(1)*D(2)+B(1)*B(2) ;}
        Double_t I() {return C(1)*C(2)+D(1)*A(2) ;}
        Double_t J() {return C(1)*D(2)+D(1)*B(2) ;}

        Double_t K() {return (A(0)*I()+B(0)*G()) ;}
        Double_t L() {return (A(0)*J()+B(0)*H()) ;}
        Double_t M() {return (C(0)*I()+D(0)*G()) ;}
        Double_t N() {return (C(0)*J()+D(0)*H()) ;}

        Double_t P() {return M()-L() ;}

        Double_t X() ;
        Double_t Y() ;
        Double_t Z() ;
        Double_t T() ;

    public :
        S8NumericSolver() ;
        S8NumericSolver(std::string name) ;
        S8NumericSolver(std::string name, Double_t n, Double_t n1, Double_t n2, Double_t n3,
        Double_t n12, Double_t n23, Double_t n31, Double_t n123) ;
        S8NumericSolver(std::string name, Double_t* n) ;
        ~S8NumericSolver();

        void Reset();
        void SetInput(Double_t n, Double_t n1, Double_t n2, Double_t n3,
        Double_t n12, Double_t n23, Double_t n31, Double_t n123) ;

        void setInput(const NumericInputGroup &);
        void SetCorr(Double_t c12a, Double_t c23a, Double_t c31a,
        Double_t c12b, Double_t c23b, Double_t c31b) ;

        void SetCorrError(Double_t c12a, Double_t c23a, Double_t c31a,
        Double_t c12b, Double_t c23b, Double_t c31b) ;

        void SetError( Int_t b );
        void SetPrecision( Int_t pr ) {fNpt = pr ;}
        void SetNbErrorIteration( Int_t pr ) {fIter = pr ;}

        void SetSolution( Int_t asol)
        {
            fForceSol = true;
            fpickSol = asol;
        }

        void SetInitialOrder(Int_t n, Int_t a = 1);
        void SetInitialOrder(std::string s1, std::string s2, std::string s3);

        Int_t Solve() ;

        void SetAverageRes( double n ) { fAveRes = n; fAveResSetup = true; }

        Double_t GetResultVec(Int_t n); 

        Double_t GetResult  (Int_t n) ;
        Double_t GetErrorInf(Int_t n, std::string opt="") ;
        Double_t GetErrorSup(Int_t n, std::string opt="") ;
        Double_t GetErrorInfVec(Int_t n);
        Double_t GetErrorSupVec(Int_t n);
        Double_t GetResult  (std::string label) ;
        Double_t GetErrorInf(std::string label, std::string opt="") ;
        Double_t GetErrorSup(std::string label, std::string opt="") ;

        Double_t ZeroFunctionFb(Double_t x) ;

        TH1 *result(const int &);
        double getError(const int &);

    private:
        TH1    *_result[8];
        double  _result_errors[8];
};

#endif

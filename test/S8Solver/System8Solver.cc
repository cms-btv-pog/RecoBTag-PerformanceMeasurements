/////////////////////////////////////////////////////////////////////////
//
//  System is a method to determine selection effiencies and
//  signal/background fractions on data with minimum monte-carlo input
//
//  It requires :
//     - 1 data sample containing 2 signal components
//     - 3 selection criteria (cuts on sensitive variable) with
//       various efficiencies on the 2 signals.
//
//  The solved system is :
//   / 1    = na + nb
//  |  q1   = ea1*na + eb1*nb
//  |  q2   = ea2*na + eb2*nb
//  /  q3   = ea3*na + eb3*nb
//  \\  q12  = ca12*ea1*ea2*na + cb12*eb1*eb2*nb
//  |  q23  = ca23*ea2*ea3*na + cb23*eb2*eb3*nb
//  |  q31  = ca31*ea3*ea1*na + cb31*eb3*eb1*nb
//   \\ q123 = ca123*ea1*ea2*ea3*na + cb123*eb2*eb3*nb
//
//   Input :
//     qxxx        : fractions of data surviving selection(s) xxx
//     caxy (cbxy) : correlation factor between cuts x and y on signal a (b)
//
//   Output : 
//     na (nb)     : fraction signal a (b)
//     eax (ebx)   : efficiency of cut x on signal a (b)
//
/////////////////////////////////////////////////////////////////////////

#include <boost/lexical_cast.hpp>

#include <TF1.h>
#include <TH1F.h>

#include "System8Solver.h"
#include "TNtupleD.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "S8NumericInput.h"

#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

using boost::lexical_cast;

/////////////////////////////////////////////////////////////////////////
System8Solver::System8Solver()
/////////////////////////////////////////////////////////////////////////
{
// Generic constructor
this->Reset();
}

/////////////////////////////////////////////////////////////////////////
System8Solver::System8Solver(double n, double n1, double n2 , double n3,
                             double n12, double n23, double n31, double n123)
/////////////////////////////////////////////////////////////////////////
{
// Constructor with input initialization :
// n    : total numer of events in data sample
// nx   : number of events passing cut x
// nxy  : number of events passing cuts x and y
// n123 : number of events passing all 3 cuts.
this->Reset();
this->SetInput(n,n1,n2,n3,n12,n23,n31,n123);
}

/////////////////////////////////////////////////////////////////////////
System8Solver::System8Solver(double* n)
/////////////////////////////////////////////////////////////////////////
{
this->Reset();
this->SetInput(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]);
}

System8Solver::~System8Solver()
{
    for(int i = 0; 8 > i; ++i)
        delete *(_result + i);
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::Reset()
/////////////////////////////////////////////////////////////////////////
{
// Reseting inputs, correlations, results and errors
kError = NONE;
fNmiss_Syst = 0;
fNmiss_Stat = 0;
fRndm = new TRandom3(12345);
Kappa = new TH1F("name","name",500,0.8,1.2);
fCovmatrix.ResizeTo(8,8);

    for(int i = 0; 8 > i; ++i)
    {
        *(_result + i) = new TH1F((string("res") + lexical_cast<string>(i)).c_str(),
                                  (string("Result") + lexical_cast<string>(i)).c_str(),
                                  200, 0, 2.);
        *(_result_errors + i) = 0;
        *(_central_values + i) = 0;
    }

for(int i=0 ; i<8 ; ++i)
  {
  fInput[i]  = 0. ;
  fIndep[i]  = 0. ; 
  fErrorinf_Stat[i] = 0. ;
  fErrorsup_Stat[i] = 0.;
  fErrorinf_Syst[i] = 0. ;
  fErrorsup_Syst[i] = 0. ;

  fCorr[i/4][i%4]    = 1. ;
  for(int j=0;j<8;j++) fCovmatrix(i,j) = 0.;
  }
fBin = 0;
fSystem.SetNpt(5000);
fResult = fSelector.GetSolution();
}

void System8Solver::setInput(const NumericInputGroup &inputGroup)
{
    const NumericInput &input = inputGroup.input;

    SetInput(input.n.all.first,
             input.n.mu.first,
             input.n.tag.first,
             input.p.all.first,
             input.n.muTag.first,
             input.p.tag.first,
             input.p.mu.first,
             input.p.muTag.first);

    const Coefficients &coef = inputGroup.coefficients;
        
    SetCorr(coef.kappaB.first,
            coef.beta.first,
            coef.delta.first,
            coef.kappaB123.first,
            coef.kappaCL.first,
            coef.alpha.first,
            coef.gamma.first,
            coef.kappaCL123.first);

    SetCorrError(sqrt(coef.kappaB.second),
                 sqrt(coef.beta.second),
                 sqrt(coef.delta.second),
                 sqrt(coef.kappaB123.second),
                 sqrt(coef.kappaCL.second),
                 sqrt(coef.alpha.second),
                 sqrt(coef.gamma.second),
                 sqrt(coef.kappaCL123.second));
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::SetInput(double n, double n1, double n2 , double n3,
                             double n12, double n23, double n31, double n123)
/////////////////////////////////////////////////////////////////////////
{
// Set the 8 input values - cf constructor
fInput[0] = n    ;
fInput[1] = n1   ;
fInput[2] = n2   ;
fInput[3] = n3   ;
fInput[4] = n12  ;
fInput[5] = n23  ;
fInput[6] = n31  ;
fInput[7] = n123 ;

fIndep[0] = n123                                         ;
fIndep[1] = n23 - fIndep[0]                              ;
fIndep[2] = n31 - fIndep[0]                              ; 
fIndep[3] = n3  - fIndep[0] - fIndep[1] - fIndep[2]      ;
fIndep[4] = n12 - fIndep[0]                              ;
fIndep[5] = n2  - n12       - fIndep[1]                  ;
fIndep[6] = n1  - n12       - fIndep[2]                  ;
fIndep[7] = n   - fIndep[6] - fIndep[5] - fIndep[4] - n3 ;
for(int i=0; i<8;++i) fIndep[i] = sqrt(fIndep[i]);
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::SetCorr(double c12a, double c23a, double c31a, double c123a,
                            double c12b, double c23b, double c31b, double c123b)
/////////////////////////////////////////////////////////////////////////
{
// Set the correlation factors
// 3 first parameters : correlation for signal a.
// 3 last parameters  : correlation for signal b.
fCorr[0][0] = c12a  ;
fCorr[1][0] = c12b  ;
fCorr[0][1] = c23a  ;
fCorr[1][1] = c23b  ;
fCorr[0][2] = c31a  ;
fCorr[1][2] = c31b  ;
fCorr[0][3] = c123a ;
fCorr[1][3] = c123b ;
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::SetCorrError(double c12a, double c23a, double c31a,double c123a,
                            double c12b, double c23b, double c31b, double c123b)
/////////////////////////////////////////////////////////////////////////
{
// Set the errors on correlation factors
// Same order than SetCorr
fCovmatrix(0,0) =  pow(c12a,2)  ;
fCovmatrix(1,1) =  pow(c23a,2)  ;
fCovmatrix(2,2) =  pow(c31a,2)  ;
fCovmatrix(3,3) =  pow(c123a,2) ;
fCovmatrix(4,4) =  pow(c12b,2)  ;
fCovmatrix(5,5) =  pow(c23b,2)  ;
fCovmatrix(6,6) =  pow(c31b,2)  ;
fCovmatrix(7,7) =  pow(c123b,2) ;

}

/////////////////////////////////////////////////////////////////////////
void System8Solver::SetCovMatrix(TMatrixDSym CovMatrix){
/////////////////////////////////////////////////////////////////////////
  if(CovMatrix.IsSymmetric())   fCovmatrix = CovMatrix;
  else { std::cout<<"The setted covariance matrix is not symmetric"<<std::endl; exit(1);}
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::SetCovMatrixSign(TMatrixDSym CovMatrixSign)
/////////////////////////////////////////////////////////////////////////
{
  // Set cov Matrix for signal
  // Diagnonal == error² 
  // Of course covariance must be sup than variance² mean ii > ij ; all j
 if(CovMatrixSign.IsSymmetric())  for(int i = 0; i<4 ; i++) for(int j = 0; j<4;j++) fCovmatrix(i,j) = CovMatrixSign(i,j);
 else {std::cout<<"Signal covariance matrix is not symmetric"<<std::endl; exit(1);}
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::SetCovMatrixBack(TMatrixDSym CovMatrixBack)
/////////////////////////////////////////////////////////////////////////
{
// Set cov Matrix for background
  if(CovMatrixBack.IsSymmetric())  for(int i = 0; i<4 ; i++) for(int j = 0; j<4;j++) fCovmatrix(i+4,j+4) = CovMatrixBack(i,j) ;
 else {std::cout<<"Background covariance matrix is not symmetric"<<std::endl; exit(1);}
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::MakeSystem(double *shift)
/////////////////////////////////////////////////////////////////////////
{
double N = fInput[0]+shift[0];
double q[7], k[8];
for(int i=0;i<7;++i) q[i] = (fInput[i+1]+shift[i+1])/N ;
fSystem.SetInput(q[0],q[1],q[2],q[3],q[4],q[5],q[6]);
 for(int i=0;i<8;++i) {k[i] = fCorr[i/4][i%4] + shift[i+8]; if(i==0) Kappa->Fill(k[i]);}
fSystem.SetCorrelation(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]);
}

/////////////////////////////////////////////////////////////////////////
int System8Solver::Solve()
/////////////////////////////////////////////////////////////////////////
{
double shift[16]={0};
MakeSystem(shift);
fSystem.Solve();
if(!fSelector.FindSolution(&fSystem)) { fNmiss_Stat = 666666; return 0; }
// bool truc =  fSelector.FindSolution(&fSystem); //                             ICCIICICICICIIC   ICI
//
    ComputeErrors();
    fitErrors();

return 1;
}

/////////////////////////////////////////////////////////////////////////
void System8Solver::ComputeErrors()
/////////////////////////////////////////////////////////////////////////
{
printf("Start computing errors -> mode : %d\n",kError);
if(kError==NONE) return;
TFile toto(Form("./result/System8/System8Output_Bin%i.root",fBin),"recreate");
TNtupleD stat("stat","stat","e1:e2:e3:e4:e5:e6:e7:e8");
TNtupleD syst("syst","syst","e1:e2:e3:e4:e5:e6:e7:e8");
double central[8]={0};
for(int i=0;i<8;++i)
  {
  central[i] = fResult[i];
  fErrorsup_Stat[i] = 0.;
  fErrorinf_Stat[i] = 0.;
  fErrorsup_Syst[i] = 0.;
  fErrorinf_Syst[i] = 0.;
  }
 double Ninf[8] = {0};
 double Nsup[8] = {0};
double shift[16] = {0.};
if(kError==STAT || kError==ALL) // Stat error
  {
  printf("   Computing stat errors\n");
  for(int i=0;i<8;++i) Ninf[i] = Nsup[i] = 0;
  for(int i=8;i<16;++i) shift[i] = 0.;
  for(int i=0;i<fIter;++i)
    {
   if(i%100==0) std::cout<<"Stat iteration number "<<i<<" / "<<fIter<<std::endl;
    double w[8] = {0};
    for(int j=0;j<8;++j) w[j] = fRndm->Gaus(0.,1.)*fIndep[j];
    shift[0] = w[0]+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7];
    shift[1] = w[0] + w[2] + w[4] + w[6];
    shift[2] = w[0] + w[1] + w[4] + w[5];
    shift[3] = w[0] + w[1] + w[2] + w[3];
    shift[4] = w[0] + w[4];
    shift[5] = w[0] + w[1];
    shift[6] = w[0] + w[2];
    shift[7] = w[0] ;
    MakeSystem(shift);
    fSystem.Solve();
    
    if(!fSelector.FindSolution(&fSystem)){fNmiss_Stat++; continue;}
    stat.Fill(fResult);
    for(int j=0;j<8;++j)
      {
      if(fResult[j]<central[j]) {fErrorinf_Stat[j]+= pow(fResult[j]-central[j],2); Ninf[j]++;}
      if(fResult[j]>central[j]) {fErrorsup_Stat[j]+= pow(fResult[j]-central[j],2); Nsup[j]++;}


        // Test if value is within X% from the central value
        //
        // if (fabs(fResult[j] - central[j]) >= .2 * central[j])
        //    continue;

        // Idea: only values from the Central solution should be
        //       kept. Therefore add something like test if new
        //       value is within, say 10% from the central value.
        //
        _result[j]->Fill(fResult[j]);

        //if(fResult[j]<central[j] && fResult[j]>0) {fErrorinf_Stat[j]+= pow(fResult[j]-central[j],2); Ninf++;}
        //if(fResult[j]>central[j] && fResult[j]<1) {fErrorsup_Stat[j]+= pow(fResult[j]-central[j],2); Nsup++;}
      }
    }
  for(int i=0;i<8;++i)
    {
    fErrorinf_Stat[i]/=Ninf[i];
    fErrorsup_Stat[i]/=Nsup[i];
    }
  }
if(kError==SYST || kError==ALL) // Syst error
  {
  printf("   Computing syst errors\n");
  for(int i=0;i<8;++i) Ninf[i] = Nsup[i] = 0;
  double Rot[8][8];
  TVectorD EigenValues; 
  TMatrixD Cht = fCovmatrix.EigenVectors(EigenValues);   
  for(int i =0; i<8;i++) for(int j =0; j<8;j++) Rot[i][j] = Cht(i,j); 
  for(int i=0;i<8;++i) shift[i] = 0.;
  for(int i=0;i<fIter;++i)
    {
    if(i%100==0) std::cout<<"Syst iteration number "<<i<<" / "<<fIter<<std::endl;
    double RotKErr[8] = {0};
    for(Int_t j=0;j<8;j++) RotKErr[j] = fRndm->Gaus(0.,1.)*TMath::Sqrt(EigenValues(j)); 
    for(Int_t j=0;j<8;++j){
      shift[j+8] = 0; 
      for(Int_t k =0;k<8;k++) shift[j+8] += Rot[j][k]*RotKErr[k]; 
    }
    MakeSystem(shift);
    fSystem.Solve();
    if(!fSelector.FindSolution(&fSystem)){ fNmiss_Syst++; continue;}
    syst.Fill(fResult);
    for(int j=0;j<8;++j)
      {
      if(fResult[j]<central[j]) {fErrorinf_Syst[j]+= pow(fResult[j]-central[j],2); Ninf[j]++;}
      if(fResult[j]>central[j]) {fErrorsup_Syst[j]+= pow(fResult[j]-central[j],2); Nsup[j]++;}
      }
    }
  for(int i=0;i<8;++i)
    {
    fErrorinf_Syst[i]/=Ninf[i];
    fErrorsup_Syst[i]/=Nsup[i];
    }
   }
for(int i=0;i<8;++i) fResult[i] = central[i];
toto.Write();
toto.Close();
}

void System8Solver::fitErrors()
{
    for(int i = 0; 8 > i; ++i)
    {
        TH1 *result = _result[i];
        if (!result->GetEntries())
        {
            cerr << "Error plot is empty for result: " << i << endl;

            continue;
        }

        cout << " Fit result: " << i << endl;

        result->Fit("gaus", "0+");

        double error = _result[i]->GetFunction("gaus")->GetParameter(2);
        _result_errors[i] = error;

        double central_value = _result[i]->GetFunction("gaus")->GetParameter(1);
        _central_values[i] = central_value;

        cout << " central[" << i << "]: " << central_value
            << "   sigma[" << i << "]: " << error << endl;
        cout << endl;
    }
}

double System8Solver::getError(const int &id)
{
    if (0 > id ||
        7 < id)

        return -1;

    return _result_errors[id];
}

double System8Solver::getCentralValue(const int &id)
{
    if (0 > id ||
        7 < id)

        return -1;

    return _central_values[id];
}

TH1 *System8Solver::result(const int &id)
{
    if (0 > id ||
        7 < id)

        return 0;

    return _result[id];
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetValue(double* tab, std::string label)
/////////////////////////////////////////////////////////////////////////
{
if(label == "na")  return tab[0];
if(label == "nb")  return tab[1];
if(label == "ea1") return tab[2];
if(label == "eb1") return tab[3];
if(label == "ea2") return tab[4];
if(label == "eb2") return tab[5];
if(label == "ea3") return tab[6];
if(label == "eb3") return tab[7];
return 0;
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetResult(Int_t n)
/////////////////////////////////////////////////////////////////////////
{
// return the system 8 results :
// n = 0 -> na
// n = 1 -> nb
// n = 2 -> ea1
// n = 3 -> eb1
// n = 4 -> ea2
// n = 5 -> eb2
// n = 6 -> ea3
// n = 7 -> eb3

if(n>=0 && n<8) return fResult[n];
return 0;
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetErrorSup(Int_t n ,std::string opt)
/////////////////////////////////////////////////////////////////////////
{
double err = 0.;
if(n>=0 && n<8 && kError!=NONE)
  {
  if(opt=="All" || opt=="" || opt == "Stat") err += fErrorsup_Stat[n];
  if(opt=="All" || opt=="" || opt == "Syst") err += fErrorsup_Syst[n];
  }
return sqrt(err);
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetErrorInf(Int_t n ,std::string opt)
/////////////////////////////////////////////////////////////////////////
{
double err = 0.;
if(n>=0 && n<8 && kError!=NONE)
  {
  if(opt=="All" || opt=="" || opt == "Stat") err += fErrorinf_Stat[n];
  if(opt=="All" || opt=="" || opt == "Syst") err += fErrorinf_Syst[n];
  }
return sqrt(err);
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetResult(std::string label)
/////////////////////////////////////////////////////////////////////////
{
// return the system 8 results :
// label must be one of those
// "na", "nb", "ea1", "eb1", "ea2", "eb2", "ea3", "eb3"
return this->GetValue(fResult, label);
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetErrorSup(std::string label,std::string opt)
/////////////////////////////////////////////////////////////////////////
{
double err = 0.;
if(opt=="All" || opt=="" || opt == "Stat") err+= this->GetValue(fErrorsup_Stat, label);
if(opt=="All" || opt=="" || opt == "Syst") err+= this->GetValue(fErrorsup_Syst, label);
return sqrt(err);
}

/////////////////////////////////////////////////////////////////////////
double System8Solver::GetErrorInf(std::string label,std::string opt)
/////////////////////////////////////////////////////////////////////////
{
double err = 0.;
if(opt=="All" || opt=="" || opt == "Stat") err+= this->GetValue(fErrorinf_Stat, label);
if(opt=="All" || opt=="" || opt == "Syst") err+= this->GetValue(fErrorinf_Syst, label);
return sqrt(err);
}


/////////////////////////////////////////////////////////////////////////
void System8Solver::DrawKappa()
/////////////////////////////////////////////////////////////////////////
{
  new TCanvas;
  Kappa->Draw();

}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
SolutionSelector::SolutionSelector()
/////////////////////////////////////////////////////////////////////////
{
for(int i=0; i<4;++i) fOrder[i]=0;
for(int j=0; j<8;++j) fSolution[j]=0;
}

/////////////////////////////////////////////////////////////////////////
SolutionSelector::~SolutionSelector()
/////////////////////////////////////////////////////////////////////////
{;}

/////////////////////////////////////////////////////////////////////////
void SolutionSelector::ResetOrder()
/////////////////////////////////////////////////////////////////////////
{
for(int i=0; i<4;++i) SetInitialOrder(i,0);
}


/////////////////////////////////////////////////////////////////////////
void SolutionSelector::SetInitialOrder(int n, int a)
/////////////////////////////////////////////////////////////////////////
{
if(a) fOrder[n]=a/abs(a);
else  fOrder[n]=0;
}

/////////////////////////////////////////////////////////////////////////
void SolutionSelector::SetInitialOrder(std::string s1, std::string s2, std::string s3)
/////////////////////////////////////////////////////////////////////////
{
// Define criteria used to resolve ambiguity in system solution.
// Force a parameter to be larger to its
// symmetric.
// n=0 -> na, nb
// n=1 -> ea1, eb1
// n=2 -> ea2, eb2
// n=3 -> ea3, eb3
// If not called, default is na > nb.
char s = s1[1];
char t = s3[1];
char u = s1[0];
char v = s3[0];
int n=-1 , a=0;
if(u=='n' && v=='n') n = 0;
if(u=='e' && v=='e')
  {
  if(s1[2]==s3[2]) n = ((int) s1[2]) - 48;
  else
    {
    printf("\nIgnore SetInitialOrder(\"%s\", \"%s\", \"%s\") : Wrong parameters.\n",s1.c_str(),s2.c_str(),s3.c_str());
    return;
    }
  }
if     (s=='a' && t=='b' && s2=="<") a = -1;
else if(s=='a' && t=='b' && s2==">") a = 1;
else if(s=='b' && t=='a' && s2=="<") a = 1;
else if(s=='b' && t=='a' && s2==">") a = -1;
if(a==0 || n==-1)
  {
  printf("\nIgnore SetInitialOrder(\"%s\", \"%s\", \"%s\") : Wrong parameters.\n",s1.c_str(),s2.c_str(),s3.c_str());
  return;
  }
SetInitialOrder(n,a);
}

/////////////////////////////////////////////////////////////////////////
bool SolutionSelector::FindSolution(System8* sys)
/////////////////////////////////////////////////////////////////////////
{
// number of solution
int nsol = sys->GetNbSolutions();
double tmp[8];
int n=0;
for(int i=0; i<nsol;++i)
  {
  bool good = true;
  for(int j=0; j<8 && good;++j) 
    {
    tmp[j]=sys->GetSolution(i,j%2,j/2);// sorted na,nb,ea1,eb1...
    if(tmp[j]>1 || tmp[j]<0) good = false;
    if(j%2==1 && fOrder[j/2]!=0) good = good && ( tmp[ fOrder[j/2] < 0 ? j-1 : j ] < tmp[ fOrder[j/2] < 0 ? j : j-1] );
    }
  if(good) // potential good solution
    {
      for(int j=0; j<8; ++j) {fSolution[j] = tmp[j]; /* std::cout<<Form(" ---------->tmp   %i   ",n)<<tmp[j]<<std::endl;*/}
      n++;
    }
  }
 if(n>1) std::cout << "Warning : several solutions compatible with the selection criteria  " << n << std::endl;
if(n==0){std::cout << "Warning : no solution compatible with the selection criteria" << std::endl; }
return n==1;
}


/////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////
System8::System8()
{
// number of points depends on efficiencies : if fa or fb close to 1 or 0 --> Increase the number of point 
  fNpt   = 2000;     
}

/////////////////////////////////////////////////////////////////////////
void   System8::SetInput(double v1, double v2, double v3, double v4, double v5, double v6, double v7)
/////////////////////////////////////////////////////////////////////////
{
fq[0] = v1;
fq[1] = v2;
fq[2] = v3;
fQ[0] = v4;
fQ[1] = v5;
fQ[2] = v6;
fQQQ  = v7;
}
/////////////////////////////////////////////////////////////////////////
void   System8::SetCorrelation(double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8)
/////////////////////////////////////////////////////////////////////////
{
kc[0][0]=v1;
kc[0][1]=v2;
kc[0][2]=v3;
kc[0][3]=v4;
kc[1][0]=v5;
kc[1][1]=v6;
kc[1][2]=v7;
kc[1][3]=v8;
}

/////////////////////////////////////////////////////////////////////////
long double System8::Gamma(int i)
/////////////////////////////////////////////////////////////////////////
{
long double res = 0;
if(i==0)
  {
  res-=R() *B(0)*C(2);
  res+=S(1)*A(0)*C(2);
  res-=S(2)*B(0)*A(2);
  res+=T(1)*A(0)*A(2);
  }
else if(i==1)
  {
  res-=R() *B(0)*D(2);
  res+=R() *D(0)*C(2);
  res-=S(0)*B(0)*C(2);
  res+=S(1)*A(0)*D(2);
  res-=S(1)*C(0)*C(2);
  res-=S(2)*B(0)*B(2);
  res+=S(2)*D(0)*A(2);
  res+=T(0)*A(0)*C(2);
  res+=T(1)*A(0)*B(2);
  res-=T(1)*C(0)*A(2);
  res-=T(2)*B(0)*A(2);
  res+=W() *A(0)*A(2);  
  }
else if(i==2)
  {
  res+=R() *D(0)*D(2);
  res-=S(0)*B(0)*D(2);
  res+=S(0)*D(0)*C(2);
  res-=S(1)*C(0)*D(2);
  res+=S(2)*D(0)*B(2);
  res+=T(0)*A(0)*D(2);   
  res-=T(0)*C(0)*C(2);
  res-=T(1)*C(0)*B(2);  
  res-=T(2)*B(0)*B(2);
  res+=T(2)*D(0)*A(2);
  res+=W() *A(0)*B(2);
  res-=W() *C(0)*A(2);  
  }
else if(i==3)
  {
  res+=S(0)*D(0)*D(2);
  res-=T(0)*C(0)*D(2);  
  res+=T(2)*D(0)*B(2);
  res-=W() *C(0)*B(2);
  } 
 if(fSign == -1)     return fScale[0]*res;
 else if(fSign == 1) return fScale[1]*res;
 return 1e12*res;   
}

//////////////////////////////////////////////////////////////////////////////////////
double System8::UKN(int i,int j) // i=0..1 : sample , j=0..3 : tag
//////////////////////////////////////////////////////////////////////////////////////
{
fSign = fSolutions[fIndex].first;
fFa   = fSolutions[fIndex].second;

double e = -1;
if(j==0)
  {
  if(i==0) e =  U(); // f_a
  else     e = -V(); // f_b
  }
else
  {
  if(i==0) e = (Q(j-1,0)+V()*UKN(1,j))/U();// eps_a^(X/Y/Z)
  else
    {
    if     (j==1) //eps_b^(X) 
      {
      e = 0.5*(-N()+fSign*TMath::Sqrt(DELTA()))/M();
      }
    else if(j==2) //eps_b^(Y) 
      {
      double ex = UKN(1,1);
      e = (A(0)-C(0)*ex)/(D(0)*ex-B(0)) ;
      }
    else if(j==3) //eps_b^(Z) 
      {
      double ex = UKN(1,1);
      e = (A(2)+B(2)*ex)/(C(2)+D(2)*ex) ;
      }
    }
  }
return e;
}

//////////////////////////////////////////////////////////////////////////////////////
long double System8::ZeroFunction(double x)
//////////////////////////////////////////////////////////////////////////////////////
{
fFa = x;
long double d = DELTA();
if(d<0) return -100;
long double Eps = 0.5*(-N() +fSign*sqrt(d))/M();
long double A1 = Gamma(0)+Eps*(Gamma(1)+Eps*(Gamma(2)+Eps*Gamma(3)));
return A1;
}

//////////////////////////////////////////////////////////////////////////////////////
void System8::EvalScale()
//////////////////////////////////////////////////////////////////////////////////////
{
  fScale[0] = 1;
  fScale[1] = 1;
  long double Scale_temp[2]  = {0};
  double n    = 7;
  double step = 1/n; 
  for(int sign = -1; sign <=1 ; sign +=2 ){
    fSign = sign;
    for(int i = 1; i<n-1;i++){
      long double A1 =  ZeroFunction(i*step);
      if(Scale_temp[(sign+1)/2] < fabs(1/A1) ) Scale_temp[(sign+1)/2] = fabs(1/A1);
    }
  }
  fScale[0]= Scale_temp[0];
  fScale[1]= Scale_temp[1];
}

//////////////////////////////////////////////////////////////////////////////////////
int   System8::Solve()
//////////////////////////////////////////////////////////////////////////////////////
{
std::pair<int, double> tmp;
fSolutions.clear();
EvalScale();
for(int i=-1; i<=1; i+=2)
  {
  fSign = i;  
  tmp.first = i;
  long double f = 0.;
  long double fold = 0.;
  double x = 0.;
  double xold = 0.;
  double u = 1./fNpt;
  for(int i=0;i<=fNpt;++i)
    {
    x = i*u ;
    f = ZeroFunction(x);
    if(fabs(f)<1e-10) {continue;}
    if(fold*f < 0)
      {
      if((x-xold) == 0) continue;
      double aa = (f-fold)/(x-xold);
      double bb = f-aa*x;
      tmp.second = -bb/aa;
      fSolutions.push_back(tmp);
      }
    fold = f;
    xold = x;
    }
  }
return fSolutions.size();
}
//////////////////////////////////////////////////////////////////////////////////////
TGraph* System8::MakeGraph(int idx, int npt)
//////////////////////////////////////////////////////////////////////////////////////
{
double u=1./npt;
double x=0, f=0;
fSign=idx;
TGraph* gr = new TGraph;
for(int i=0;i<=npt;++i)
  {
  x = i*u ;
  f = ZeroFunction(x);
  gr->SetPoint(i,x,f);         
  }
return gr;  
}
  
void demo(Double_t* inputs, Double_t* corr, Int_t k, Double_t *ans)
{
System8Solver sol(inputs);
//sol.SetCorr(kappa_b,beta,delta,C123B,kappa_cl,alpha,chi,C123CL);
sol.SetCorr(corr[0], corr[1], corr[2], corr[3], corr[4], corr[5], corr[6], corr[7]);

sol.SetCorrError(0.,0.,0.,0.,0.,0.,0.,0.);
System8Solver::ErrorType myErrType = System8Solver::STAT;
sol.SetError(myErrType, 10000);
//$$
sol.SetNbErrorIteration(1000);
//$$
sol.SetInitialOrder(1,1);
 if(!sol.Solve())
 {
 sol.SetInitialOrder(0,-1);   // OK en principe...
 if(!sol.Solve())
   {
   sol.SetInitialOrder(3,1);   // OK en principe...
   if(!sol.Solve())
     {
     sol.SetInitialOrder(2,1);   // OK en principe...
     if(!sol.Solve())
       cout << "Arrggg : Impossible de déterminer une solution physique" << endl;
     }
   }
 }

// Final result :

 if( k == 0 )
 {
   cout << endl;
   cout << "Output values : " <<endl;
   cout << "nb/n  : " << sol.GetResult("na")  << "  + " << sol.GetErrorSup("na")  << "  - " << sol.GetErrorInf("na") << endl;
   cout << "ncl/n : " << sol.GetResult("nb")  << "  + " << sol.GetErrorSup("nb")  << "  - " << sol.GetErrorInf("nb") << endl;
   cout << "eopp  : " << sol.GetResult("ea3") << "  + " << sol.GetErrorSup("ea3") << "  - " << sol.GetErrorInf("ea3") << endl;
   cout << "ropp  : " << sol.GetResult("eb3") << "  + " << sol.GetErrorSup("eb3") << "  - " << sol.GetErrorInf("eb3") << endl;
   cout << "emu   : " << sol.GetResult("ea1") << "  + " << sol.GetErrorSup("ea1") << "  - " << sol.GetErrorInf("ea1") << endl;
   cout << "rmu   : " << sol.GetResult("eb1") << "  + " << sol.GetErrorSup("eb1") << "  - " << sol.GetErrorInf("eb1") << endl;
   cout << "etag  : " << sol.GetResult("ea2") << "  + " << sol.GetErrorSup("ea2") << "  - " << sol.GetErrorInf("ea2") << endl;
   cout << "rtag  : " << sol.GetResult("eb2") << "  + " << sol.GetErrorSup("eb2") << "  - " << sol.GetErrorInf("eb2")  << endl;
   cout << endl;
  
   cout << "-------------------------------------" << endl;
   cout << "nb/n  : " << sol.GetResult("na")  << "  + " << sol.GetErrorSup("na")  << "  - " << sol.GetErrorInf("na") << endl;
   cout << "emu   : " << sol.GetResult("ea1") << "  + " << sol.GetErrorSup("ea1") << "  - " << sol.GetErrorInf("ea1") << endl;
   cout << "etag  : " << sol.GetResult("ea2") << "  + " << sol.GetErrorSup("ea2") << "  - " << sol.GetErrorInf("ea2") << endl;
   cout << "eopp  : " << sol.GetResult("ea3") << "  + " << sol.GetErrorSup("ea3") << "  - " << sol.GetErrorInf("ea3") << endl;
   cout << "-------------------------------------" << endl;

   cout << "  yfrac[" << k << "] = " << sol.GetResult("na") 
       << "; eyfrac[" << k << "] = " << (sol.GetErrorSup("na")+sol.GetErrorInf("na"))/2. << ";"  << endl;
   cout << "  zfrac[" << k << "] = " << sol.GetResult("nb") 
       << "; ezfrac[" << k << "] = " << (sol.GetErrorSup("nb")+sol.GetErrorInf("nb"))/2. << ";"  << endl;
   cout << "   yopp[" << k << "] = " << sol.GetResult("ea3") 
       << ";  eyopp[" << k << "] = " << (sol.GetErrorSup("ea3")+sol.GetErrorInf("ea3"))/2. << ";"  << endl;
   cout << "   zopp[" << k << "] = " << sol.GetResult("eb3") 
       << ";  ezopp[" << k << "] = " << (sol.GetErrorSup("eb3")+sol.GetErrorInf("eb3"))/2. << ";"  << endl;
   cout << "    ymu[" << k << "] = " << sol.GetResult("ea1") 
       << ";   eymu[" << k << "] = " << (sol.GetErrorSup("ea1")+sol.GetErrorInf("ea1"))/2. << ";"  << endl;
   cout << "    zmu[" << k << "] = " << sol.GetResult("eb1") 
       << ";   ezmu[" << k << "] = " << (sol.GetErrorSup("eb1")+sol.GetErrorInf("eb1"))/2. << ";"  << endl;
   cout << "   ytag[" << k << "] = " << sol.GetResult("ea2") 
       << ";  eytag[" << k << "] = " << (sol.GetErrorSup("ea2")+sol.GetErrorInf("ea2"))/2. << ";"  << endl;
   cout << "   ztag[" << k << "] = " << sol.GetResult("eb2") 
       << ";  eztag[" << k << "] = " << (sol.GetErrorSup("eb2")+sol.GetErrorInf("eb2"))/2. << ";"  << endl;
   cout << endl;
 }
 ans[0] = sol.GetResult("ea2");
 ans[1] = sol.GetErrorSup("ea2");
 ans[2] = sol.GetErrorInf("ea2");
 return ;
}


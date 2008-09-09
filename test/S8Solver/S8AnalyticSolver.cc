
#include "S8AnalyticSolver.h"

#include <math.h>

#include<iostream>

//______________________________________________________________________________
// Function to solve system8 analytically
// Lars Sonnenschein 17/August/2003
//
// Francisco Yumiceva modified Oct. 2007
//______________________________________________________________________________

//ClassImp(S8AnalyticSolver)

S8AnalyticSolver::S8AnalyticSolver() {
	fVerbose = false;
}


void S8AnalyticSolver::Solve(std::map< TString, double > input) {

	
	double data[8];
	
	data[0] = input["n"];
	data[1] = input["p"];
	data[2] = input["nMu"];
	data[3] = input["pMu"];
	data[4] = input["nTag"];
	data[5] = input["pTag"];
	data[6] = input["nMuTag"];
	data[7] = input["pMuTag"];
	
	//const int dim=8;
  double n = data[0];
  double p = data[1];
  double nmu = data[2];
  double pmu = data[3];
  double nsvx = data[4];
  double psvx = data[5];
  double nall = data[6];
  double pall = data[7];

  if (fVerbose) {
	  std::cout << "S8AnalyticSolver input parameters:" << std::endl;
	  std::cout << "n=" << n << std::endl;
	  std::cout << "p=" << p << std::endl;
	  std::cout << "nmu=" << nmu << std::endl;
	  std::cout << "pmu=" << pmu << std::endl;
	  std::cout << "nsvx=" << nsvx << std::endl;
	  std::cout << "psvx=" << psvx << std::endl;
	  std::cout << "nall=" << nall << std::endl;
	  std::cout << "pall=" << pall << std::endl;
  }
   
  double denom = sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx);
  double denomX =   pmu*nsvx+pall*n-nall*p-psvx*nmu;
  if (fVerbose) std::cout << denom << " " << denomX << std::endl;

  double Esvx1 = 1./2./(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu+sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx));

  double Esvx2 = 1./2./(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu-sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx));

  if (fVerbose) std::cout << " ESVX " << Esvx1 << " " << Esvx2<< std::endl;
  //error of first Esvx:
			 
  double a = sqr(-1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu-sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*pmu+1./2./(-p*nmu+n*pmu)*(pall-1./2.*1./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pall*pmu*nsvx+2.*sqr(pall)*n-2.*pall*nall*p-2.*psvx*pall*nmu+4.*pmu*nall*psvx)))*n; 

    a += sqr(1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu+sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*nmu+1./2./(-p*nmu+n*pmu)*(-nall+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*nall-2.*pall*n*nall+2.*sqr(nall)*p-2.*nall*psvx*nmu+4.*nmu*nsvx*pall)))*p;

    a += sqr(1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu+sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*p+1./2./(-p*nmu+n*pmu)*(-psvx+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*nsvx*psvx*pmu-2.*psvx*n*pall-2.*nall*p*psvx+2.*sqr(psvx)*nmu+4.*nsvx*p*pall)))*nmu*(1-nmu/n);

    a += sqr(-1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu+sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*n+1./2./(-p*nmu+n*pmu)*(nsvx+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(2.*pmu*sqr(nsvx)-2.*nsvx*pall*n-2.*nsvx*nall*p-2.*nsvx*psvx*nmu+4.*n*nall*psvx)))*pmu*(1-pmu/p);

    a += sqr(1./2./(-p*nmu+n*pmu)*(pmu+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(2.*sqr(pmu)*nsvx-2.*pmu*pall*n-2.*pmu*nall*p-2.*pmu*psvx*nmu+4.*p*nmu*pall)))*nsvx*(1-nsvx/n);
    
    a += sqr(1./2./(-p*nmu+n*pmu)*(-nmu+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*nmu-2.*pall*n*nmu-2.*nall*p*nmu+2.*psvx*sqr(nmu)+4.*n*pmu*nall)))*psvx*(1-psvx/p);
    
    a += sqr(1./2./(-p*nmu+n*pmu)*(-p+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*p-2.*pall*n*p+2.*nall*sqr(p)-2.*p*psvx*nmu+4.*n*pmu*psvx)))*nall*(1-nall/n);
    
    a += sqr(1./2./(-p*nmu+n*pmu)*(n+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*n+2.*pall*sqr(n)-2.*n*nall*p-2.*n*psvx*nmu+4.*p*nmu*nsvx)))*pall*(1-pall/p);


  double Esvx1err = sqrt(a);

  

  double b = sqr(-1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu+sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*pmu+1./2./(-p*nmu+n*pmu)*(pall-1./2.*1./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pall*pmu*nsvx+2.*sqr(pall)*n-2.*pall*nall*p-2.*psvx*pall*nmu+4.*pmu*nall*psvx)))*n;


  b += sqr(1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu-sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*nmu+1./2./(-p*nmu+n*pmu)*(-nall+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*nall-2.*pall*n*nall+2.*sqr(nall)*p-2.*nall*psvx*nmu+4.*nmu*nsvx*pall)))*p;
  
  b += sqr(1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu-sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*p+1./2./(-p*nmu+n*pmu)*(-psvx+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*nsvx*psvx*pmu-2.*psvx*n*pall-2.*nall*p*psvx+2.*sqr(psvx)*nmu+4.*nsvx*p*pall)))*nmu*(1-nmu/n);

  b += sqr(-1./2./sqr(-p*nmu+n*pmu)*(pmu*nsvx+pall*n-nall*p-psvx*nmu-sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx))*n+1./2./(-p*nmu+n*pmu)*(nsvx+1./2./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(2.*pmu*sqr(nsvx)-2.*nsvx*pall*n-2.*nsvx*nall*p-2.*nsvx*psvx*nmu+4.*n*nall*psvx)))*pmu*(1-pmu/p);

  b += sqr(1./2./(-p*nmu+n*pmu)*(pmu-1./2.*1./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(2.*sqr(pmu)*nsvx-2.*pmu*pall*n-2.*pmu*nall*p-2.*pmu*psvx*nmu+4.*p*nmu*pall)))*nsvx*(1-nsvx/n);

  b += sqr(1./2./(-p*nmu+n*pmu)*(-nmu-1./2.*1./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*nmu-2.*pall*n*nmu-2.*nall*p*nmu+2.*psvx*sqr(nmu)+4.*n*pmu*nall)))*psvx*(1-psvx/p);
    
  b += sqr(1./2./(-p*nmu+n*pmu)*(-p-1./2.*1./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*p-2.*pall*n*p+2.*nall*sqr(p)-2.*p*psvx*nmu+4.*n*pmu*psvx)))*nall*(1-nall/n);

  b += sqr(1./2./(-p*nmu+n*pmu)*(n-1./2.*1./sqrt(sqr(pmu)*sqr(nsvx)-2.*pmu*nsvx*pall*n-2.*pmu*nsvx*nall*p-2.*pmu*nsvx*psvx*nmu+sqr(pall)*sqr(n)-2.*pall*n*nall*p-2.*pall*n*psvx*nmu+sqr(nall)*sqr(p)-2.*nall*p*psvx*nmu+sqr(psvx)*sqr(nmu)+4.*p*nmu*nsvx*pall+4.*n*pmu*nall*psvx)*(-2.*pmu*nsvx*n+2.*pall*sqr(n)-2.*n*nall*p-2.*n*psvx*nmu+4.*p*nmu*nsvx)))*pall*(1-pall/p);


  //error of second Esvx:
  double Esvx2err = sqrt(b);

  if (fVerbose) std::cout << " ESVX " << Esvx1 << " " << Esvx2<< std::endl;
  double Emu1 = -Esvx1*(-nall*pmu+nmu*pall)/(nall*psvx+nsvx*Esvx1*pmu-psvx*Esvx1*nmu-nsvx*pall);
  if (fVerbose) std::cout << " Emu1 " << Emu1 << std::endl;
  double Emu2 = -Esvx2*(-nall*pmu+nmu*pall)/(nall*psvx+nsvx*Esvx2*pmu-psvx*Esvx2*nmu-nsvx*pall);
  if (fVerbose) std::cout << " Emu2 " << Emu2 << std::endl;

  double Emu1err = sqr((-nall*pmu+pall*nmu)/(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)-Esvx1*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*(-nsvx*pmu+psvx*nmu))*sqr(Esvx1);
  Emu1err += sqr(-Esvx1*pmu/(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)+Esvx1*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*psvx)*nall*(1-nall/n);
  Emu1err += sqr(-Esvx1*nall/(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)+sqr(Esvx1)*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*nsvx)*pmu*(1-pmu/p);
  Emu1err += sqr(Esvx1*nmu/(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)-Esvx1*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*nsvx)*pall*(1-pall/p);
  Emu1err += sqr(Esvx1*pall/(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)-sqr(Esvx1)*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*psvx)*nmu*(1-nmu/n);
  Emu1err += sqr(-Esvx1*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*(-nall+Esvx1*nmu))*psvx*(1-psvx/p);
  Emu1err += sqr(-Esvx1*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx1*nsvx*pmu+Esvx1*psvx*nmu+nsvx*pall)*(-Esvx1*pmu+pall))*nsvx*(1-nsvx/n);
  Emu1err = sqrt(Emu1err);

  double Emu2err = sqr((-nall*pmu+pall*nmu)/(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)-Esvx2*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*(-nsvx*pmu+psvx*nmu))*sqr(Esvx2);
  Emu2err += sqr(-Esvx2*pmu/(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)+Esvx2*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*psvx)*nall*(1-nall/n);
  Emu2err += sqr(-Esvx2*nall/(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)+sqr(Esvx2)*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*nsvx)*pmu*(1-pmu/p);
  Emu2err += sqr(Esvx2*nmu/(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)-Esvx2*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*nsvx)*pall*(1-pall/p);
  Emu2err += sqr(Esvx2*pall/(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)-sqr(Esvx2)*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*psvx)*nmu*(1-nmu/n);
  Emu2err += sqr(-Esvx2*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*(-nall+Esvx2*nmu))*psvx*(1-psvx/p);
  Emu2err += sqr(-Esvx2*(-nall*pmu+pall*nmu)/sqr(-nall*psvx-Esvx2*nsvx*pmu+Esvx2*psvx*nmu+nsvx*pall)*(-Esvx2*pmu+pall))*nsvx*(1-nsvx/n);
  Emu2err = sqrt(Emu1err);




  double Rsvx1 = (-nall*Esvx1*pmu+nall*pall+Esvx1*Emu1*pmu*nsvx-Emu1*nsvx*pall)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu);
  double Rsvx2 = (-nall*Esvx2*pmu+nall*pall+Esvx2*Emu2*pmu*nsvx-Emu2*nsvx*pall)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu);


  double Rsvx1err = sqr(-(Esvx1*pmu-pall)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)+(nall*Esvx1*pmu-nall*pall-Esvx1*Emu1*pmu*nsvx+Emu1*pall*nsvx)/sqr(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)*(-Emu1*p+pmu))*nall*(1-nall/n);
  Rsvx1err += sqr(-(nall*pmu-Emu1*pmu*nsvx)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)+(nall*Esvx1*pmu-nall*pall-Esvx1*Emu1*pmu*nsvx+Emu1*pall*nsvx)/sqr(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)*(nmu*Emu1*p-nmu*pmu))*sqr(Esvx1err);
  Rsvx1err += sqr(-(nall*Esvx1-Esvx1*Emu1*nsvx)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)+(nall*Esvx1*pmu-nall*pall-Esvx1*Emu1*pmu*nsvx+Emu1*pall*nsvx)/sqr(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)*(nall-Esvx1*nmu))*pmu*(1-pmu/p);
  Rsvx1err += sqr(-(-nall+Emu1*nsvx)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu))*pall;
  Rsvx1err += sqr(-(-Esvx1*nsvx*pmu+nsvx*pall)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)+(nall*Esvx1*pmu-nall*pall-Esvx1*Emu1*pmu*nsvx+Emu1*pall*nsvx)/sqr(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)*(-nall*p+Esvx1*p*nmu))*sqr(Emu1err);
  Rsvx1err += sqr(-(-pmu*Emu1*Esvx1+Emu1*pall)/(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu))*nsvx*(1-nsvx/n);
  Rsvx1err += sqr((nall*Esvx1*pmu-nall*pall-Esvx1*Emu1*pmu*nsvx+Emu1*pall*nsvx)/sqr(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)*(-nall*Emu1+Esvx1*nmu*Emu1))*p;
  Rsvx1err += sqr((nall*Esvx1*pmu-nall*pall-Esvx1*Emu1*pmu*nsvx+Emu1*pall*nsvx)/sqr(-nall*Emu1*p+nall*pmu+Esvx1*nmu*Emu1*p-Esvx1*nmu*pmu)*(Emu1*Esvx1*p-Esvx1*pmu))*nmu*(1-nmu/n);
  Rsvx1err = sqrt(Rsvx1err);


  double Rsvx2err = sqr(-(Esvx2*pmu-pall)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)+(nall*Esvx2*pmu-nall*pall-Esvx2*Emu2*pmu*nsvx+Emu2*pall*nsvx)/sqr(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)*(-Emu2*p+pmu))*nall*(1-nall/n);
  Rsvx2err += sqr(-(nall*pmu-Emu2*pmu*nsvx)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)+(nall*Esvx2*pmu-nall*pall-Esvx2*Emu2*pmu*nsvx+Emu2*pall*nsvx)/sqr(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)*(nmu*Emu2*p-nmu*pmu))*sqr(Esvx2err);
  Rsvx2err += sqr(-(nall*Esvx2-Esvx2*Emu2*nsvx)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)+(nall*Esvx2*pmu-nall*pall-Esvx2*Emu2*pmu*nsvx+Emu2*pall*nsvx)/sqr(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)*(nall-Esvx2*nmu))*pmu*(1-pmu/p);
  Rsvx2err += sqr(-(-nall+Emu2*nsvx)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu))*pall;
  Rsvx2err += sqr(-(-Esvx2*nsvx*pmu+nsvx*pall)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)+(nall*Esvx2*pmu-nall*pall-Esvx2*Emu2*pmu*nsvx+Emu2*pall*nsvx)/sqr(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)*(-nall*p+Esvx2*p*nmu))*sqr(Emu2err);
  Rsvx2err += sqr(-(-pmu*Emu2*Esvx2+Emu2*pall)/(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu))*nsvx*(1-nsvx/n);
  Rsvx2err += sqr((nall*Esvx2*pmu-nall*pall-Esvx2*Emu2*pmu*nsvx+Emu2*pall*nsvx)/sqr(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)*(-nall*Emu2+Esvx2*nmu*Emu2))*p;
  Rsvx2err += sqr((nall*Esvx2*pmu-nall*pall-Esvx2*Emu2*pmu*nsvx+Emu2*pall*nsvx)/sqr(-nall*Emu2*p+nall*pmu+Esvx2*nmu*Emu2*p-Esvx2*nmu*pmu)*(Emu2*Esvx2*p-Esvx2*pmu))*nmu*(1-nmu/n);
  Rsvx2err = sqrt(Rsvx2err);



  double Rmu1 = Rsvx1*Emu1*(-nall+Esvx1*nmu)/(-nall*Esvx1+Emu1*Esvx1*nsvx-Rsvx1*Emu1*nsvx+Esvx1*Rsvx1*nmu);
  // double Rmu2 = Rsvx2*Emu2*(-nall+Esvx2*nmu)/(-nall*Esvx2+Emu2*Esvx2*nsvx-Rsvx2*Emu2*nsvx+Esvx2*Rsvx2*nmu);


  double Rmu1err=sqr(Emu1*(-nall+Esvx1*nmu)/(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)-Rsvx1*Emu1*(-nall+Esvx1*nmu)/sqr(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)*(Esvx1*nmu-Emu1*nsvx))*sqr(Rsvx1err);
  Rmu1err += sqr(Rsvx1*(-nall+Esvx1*nmu)/(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)-Rsvx1*Emu1*(-nall+Esvx1*nmu)/sqr(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)*(Esvx1*nsvx-Rsvx1*nsvx))*sqr(Emu1err);
  Rmu1err += sqr(-Rsvx1*Emu1/(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)+Rsvx1*Emu1*(-nall+Esvx1*nmu)/sqr(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)*Esvx1)*nall*(1-nall/n);
  Rmu1err += sqr(Rsvx1*Emu1*nmu/(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)-Rsvx1*Emu1*(-nall+Esvx1*nmu)/sqr(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)*(-nall+Emu1*nsvx+nmu*Rsvx1))*sqr(Esvx1err);
  Rmu1err += sqr(Rsvx1*Emu1*Esvx1/(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)-sqr(Rsvx1)*Emu1*(-nall+Esvx1*nmu)/sqr(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)*Esvx1)*nmu*(1-nmu/n);
  Rmu1err += sqr(-Rsvx1*Emu1*(-nall+Esvx1*nmu)/sqr(-nall*Esvx1+Emu1*Esvx1*nsvx+Esvx1*nmu*Rsvx1-Rsvx1*Emu1*nsvx)*(Emu1*Esvx1-Rsvx1*Emu1))*nsvx*(1-nsvx/n);
  Rmu1err = sqrt(Rmu1err);
  
  double Rmu2err=sqr(Emu2*(-nall+Esvx2*nmu)/(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)-Rsvx2*Emu2*(-nall+Esvx2*nmu)/sqr(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)*(Esvx2*nmu-Emu2*nsvx))*sqr(Rsvx2err);
  Rmu2err += sqr(Rsvx2*(-nall+Esvx2*nmu)/(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)-Rsvx2*Emu2*(-nall+Esvx2*nmu)/sqr(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)*(Esvx2*nsvx-Rsvx2*nsvx))*sqr(Emu2err);
  Rmu2err += sqr(-Rsvx2*Emu2/(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)+Rsvx2*Emu2*(-nall+Esvx2*nmu)/sqr(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)*Esvx2)*nall*(1-nall/n);
  Rmu2err += sqr(Rsvx2*Emu2*nmu/(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)-Rsvx2*Emu2*(-nall+Esvx2*nmu)/sqr(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)*(-nall+Emu2*nsvx+nmu*Rsvx2))*sqr(Esvx2err);
  Rmu2err += sqr(Rsvx2*Emu2*Esvx2/(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)-sqr(Rsvx2)*Emu2*(-nall+Esvx2*nmu)/sqr(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)*Esvx2)*nmu*(1-nmu/n);
  Rmu2err += sqr(-Rsvx2*Emu2*(-nall+Esvx2*nmu)/sqr(-nall*Esvx2+Emu2*Esvx2*nsvx+Esvx2*nmu*Rsvx2-Rsvx2*Emu2*nsvx)*(Emu2*Esvx2-Rsvx2*Emu2))*nsvx*(1-nsvx/n);
  Rmu2err=sqrt(Rmu2err);


  double nb1 = -(-nall+Rsvx1*nmu)/Emu1/(Esvx1-Rsvx1);
  // double nb2 = -(-nall+Rsvx2*nmu)/Emu2/(Esvx2-Rsvx2);

  double nb1err=sqrt(1./sqr(Emu1*(Esvx1-Rsvx1))*nall*(1-nall/n) + sqr((nmu*Emu1*(Esvx1-Rsvx1)-(-nall+Rsvx1*nmu)*Emu1)/sqr(Emu1*(Esvx1-Rsvx1)))*sqr(Rsvx1err) + 1./sqr(Emu1*(Esvx1-Rsvx1))*nmu*(1-nmu/n) + sqr((-nall+Rsvx1*nmu)/(sqr(Emu1)*(Esvx1-Rsvx1)))*sqr(Emu1err) + sqr((-nall+Rsvx1*nmu)*Emu1/sqr(Emu1*Esvx1+Emu1*Rsvx1))*sqr(Esvx1err));
  // double nb2err=sqrt(1./sqr(Emu2*(Esvx2-Rsvx2))*nall*(1-nall/n) + sqr((nmu*Emu2*(Esvx2-Rsvx2)-(-nall+Rsvx2*nmu)*Emu2)/sqr(Emu2*(Esvx2-Rsvx2)))*sqr(Rsvx2err) + 1./sqr(Emu2*(Esvx2-Rsvx2))*nmu*(1-nmu/n) + sqr((-nall+Rsvx2*nmu)/(sqr(Emu2)*(Esvx2-Rsvx2)))*sqr(Emu2err) + sqr((-nall+Rsvx2*nmu)*Emu2/sqr(Emu2*Esvx2+Emu2*Rsvx2))*sqr(Esvx2err));


  double pb1 = -(-pall+Rmu1*Rsvx1*p)/(Emu1*Esvx1-Rmu1*Rsvx1);
  // double pb2 = -(-pall+Rmu2*Rsvx2*p)/(Emu2*Esvx2-Rmu2*Rsvx2);

  double pb1err=sqrt(sqr(1./(Emu1*Esvx1-Rmu1*Rsvx1))*pall*(1-pall/p) + sqr((-Rsvx1*p*(Emu1*Esvx1-Rmu1*Rsvx1)+(pall-Rmu1*Rsvx1*p)*Rsvx1)/sqr(Emu1*Esvx1-Rmu1*Rsvx1))*sqr(Rmu1) + sqr((-Rmu1*p*(Emu1*Esvx1-Rmu1*Rsvx1)+(pall-Rmu1*Rsvx1*p)*Rmu1)/sqr(Emu1*Esvx1-Rmu1*Rsvx1))*sqr(Rsvx1) + sqr((Rmu1*Rsvx1)/(Emu1*Esvx1-Rmu1*Rsvx1))*p + sqr((-pall+Rmu1*Rsvx1*p)*Esvx1/sqr(Emu1*Esvx1-Rmu1*Rsvx1))*sqr(Emu1err) + sqr((-pall+Rmu1*Rsvx1*p)*Emu1/sqr(Emu1*Esvx1-Rmu1*Rsvx1))*sqr(Esvx1err));

  // double pb2err=sqrt(sqr(1./(Emu2*Esvx2-Rmu2*Rsvx2))*pall*(1-pall/p) + sqr((-Rsvx2*p*(Emu2*Esvx2-Rmu2*Rsvx2)+(pall-Rmu2*Rsvx2*p)*Rsvx2)/sqr(Emu2*Esvx2-Rmu2*Rsvx2))*sqr(Rmu2) + sqr((-Rmu2*p*(Emu2*Esvx2-Rmu2*Rsvx2)+(pall-Rmu2*Rsvx2*p)*Rmu2)/sqr(Emu2*Esvx2-Rmu2*Rsvx2))*sqr(Rsvx2) + sqr((Rmu2*Rsvx2)/(Emu2*Esvx2-Rmu2*Rsvx2))*p + sqr((-pall+Rmu2*Rsvx2*p)*Esvx2/sqr(Emu2*Esvx2-Rmu2*Rsvx2))*sqr(Emu2err) + sqr((-pall+Rmu2*Rsvx2*p)*Emu2/sqr(Emu2*Esvx2-Rmu2*Rsvx2))*sqr(Esvx2err));


  double nq1 = n-nb1;
  // double nq2 = n-nb2;

  double nq1err=sqrt(n+sqr(nb1err));
  // double nq2err=sqrt(n+sqr(nb2err));



  double pq1 = p-pb1;
  // double pq2 = p-pb2;

  double pq1err=sqrt(p+sqr(pb1err));
  // double pq2err=sqrt(p+sqr(pb2err));

  //std::cout << "First solution:" << std::endl;
  //std::cout << "Esvx=" << Esvx1 << " Emu=" << Emu1 << " Rsvx=" << Rsvx1 << " Rmu=" << Rmu1 
  //     << " nb=" << nb1 << " pb=" << pb1 << " nq=" << nq1 << " pq=" << pq1 << std::endl;
  //std::cout << "Second solution:" << std::endl;
  //std::cout << "Esvx=" << Esvx2 << " Emu=" << Emu2 << " Rsvx=" << Rsvx2 << " Rmu=" << Rmu2 
  //     << " nb=" << nb2 << " pb=" << pb2 << " nq=" << nq2 << " pq=" << pq2 << std::endl;

  if (fVerbose) {
	  std::cout << "S8AnalyticSolver output parameters:" << std::endl;
	  std::cout << "Esvx=" << Esvx1 << "+/-" << Esvx1err << std::endl;
	  std::cout << "Emu=" << Emu1 << "+/-" << Emu1err << std::endl; 
	  std::cout << "Rsvx=" << Rsvx1 << "+/-" << Rsvx1err << std::endl;
	  std::cout << "Rmu=" << Rmu1 << "+/-" << Rmu1err << std::endl;
	  std::cout << "nb=" << nb1 << "+/-" << nb1err << std::endl;
	  std::cout << "pb=" << pb1 << "+/-" << pb1err << std::endl;
	  std::cout << "nq=" << nq1 << "+/-" << nq1err << std::endl;
	  std::cout << "pq=" << pq1 << "+/-" << pq1err << std::endl;
  }

  fsolution["effTag_b"] = Esvx1;
  fsolution["effTag_cl"] = Rsvx1;
  fsolution["effMu_b"] = Emu1;
  fsolution["effMu_cl"] = Rmu1;
  fsolution["n_b"] = nb1;
  fsolution["n_cl"] = nq1;
  fsolution["p_b"] = pb1;
  fsolution["p_cl"] = pq1;

  fsolutionErr["effTag_b"] = Esvx1err;
  fsolutionErr["effTag_cl"] = Rsvx1err;
  fsolutionErr["effMu_b"] = Emu1err;
  fsolutionErr["effMu_cl"] = Rmu1err;
  fsolutionErr["n_b"] = nb1err;
  fsolutionErr["n_cl"] = nq1err;
  fsolutionErr["p_b"] = pb1err;
  fsolutionErr["p_cl"] = pq1err;


	
}









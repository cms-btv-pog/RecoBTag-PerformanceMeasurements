
#include "S8fcn.h"

#include "TMath.h"

#include<iostream>

//______________________________________________________________________
double S8fcn::operator() (const std::vector<double>& xval) const {
  
  double f = 0.0;

  Double_t c[8] = {0} , bmu[8] = {0}, cmu[8] = {0};
  bmu[0] = xval[0] + xval[1]; // 0
  bmu[1] = xval[2] + xval[3]; // 1
  bmu[2] = xval[4]*xval[0] + xval[5]*xval[1]; //nmu 2
  bmu[3] = fdelta * xval[4]*xval[2] + fgamma * xval[5]*xval[3]; //pmu 3
  bmu[4] = xval[6]*xval[0] + xval[7]*xval[1]; //nsvx 4
  bmu[5] = fbeta * xval[6]*xval[2] + falpha * xval[7]*xval[3]; //psvx 5
  bmu[6] = fkb * xval[4]*xval[6]*xval[0] + fkcl * xval[5]*xval[7]*xval[1]; //nall 6
  bmu[7] = fdelta* fbeta * fkb * xval[4]*xval[6]*xval[2] + fgamma * falpha * fkcl * xval[5]*xval[7]*xval[3]; //nall 7 

  
  c[0] = fdata[0]-fdata[2]-fdata[4]+fdata[6];
  c[1] = fdata[1]-fdata[3]-fdata[5]+fdata[7];
  c[2] = fdata[2]-fdata[6];
  c[3] = fdata[3]-fdata[7];  
  c[4] = fdata[4]-fdata[6];  
  c[5] = fdata[5]-fdata[7];  
  c[6] = fdata[6];  
  c[7] = fdata[7];
  
  cmu[0] = bmu[0]-bmu[2]-bmu[4]+bmu[6];
  cmu[1] = bmu[1]-bmu[3]-bmu[5]+bmu[7]; 
  cmu[2] = bmu[2]-bmu[6];
  cmu[3] = bmu[3]-bmu[7]; 
  cmu[4] = bmu[4]-bmu[6];
  cmu[5] = bmu[5]-bmu[7]; 
  cmu[6] = bmu[6];
  cmu[7] = bmu[7];

  for(int i=0; i<8; i++) f += c[i]*TMath::Log(cmu[i])-cmu[i];
  f = -2.0*f;
  

  // new approach:
  //for (int i=0; i<8; i++) {
//	  f += (fdata[i] - bmu[i]);
  //}
  
  //std::cout << "data0= " << fdata[0] << " fbeta= " << fbeta << " f= " << f << std::endl;
  return f;
}


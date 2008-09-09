#include "S8FitSolver.h"

#include "TMath.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameterState.h"

#include<iostream>
#include<cmath>

//______________________________________________________________________________
// Function to solve system8 analytically
// Lars Sonnenschein 17/August/2003
//
// Francisco Yumiceva modified Oct. 2007
//______________________________________________________________________________

//ClassImp(S8FitSolver)

using namespace ROOT::Minuit2;


S8FitSolver::S8FitSolver() {
        fVerbose = false;
	thePDF = new S8fcn();
	theFitter    = new VariableMetricMinimizer();
}

void S8FitSolver::Init(std::map< TString, double > input) {
	finitmap = input;
}
void S8FitSolver::Solve(std::map< TString, double > input) {

  std::vector<double> data;
  //std::cout << " begins Solve" << std::endl;
  
  data.push_back(input["n"]);
  data.push_back(input["p"]);
  data.push_back(input["nMu"]);
  data.push_back(input["pMu"]);
  data.push_back(input["nTag"]);
  data.push_back(input["pTag"]);
  data.push_back(input["nMuTag"]);
  data.push_back(input["pMuTag"]);
				 
  //std::cout << " got inputs" << std::endl;
  thePDF->SetData(data);
  //std::cout << "pdf fed with data" << std::endl;
  thePDF->SetCorr(input["kappa_b"],input["beta"],input["kappa_cl"],input["alpha"],input["delta"],input["gamma"]);
  //std::cout << "pdf configured" << std::endl;

  // check if init values are posite otherwise zero
  for( std::map<TString,double>::const_iterator im = finitmap.begin(); im != finitmap.end(); ++ im) {
	  if (im->second < 0 ) finitmap[im->first] = 0.;
	  if (isnan(im->second) || std::isinf(im->second) ) finitmap[im->first] = 0.;
	  if ( ( (im->first == "effMu_b") || (im->first == "effTag_b") ) &&
		   ( im->second > 1 || im->second <0) ) finitmap[im->first] = 0.6;
	  if ( ( (im->first == "effMu_cl") || (im->first == "effTag_cl")) &&
		   ( im->second > 1 || im->second <0) ) finitmap[im->first] = 0.1;
  }
  
  MnUserParameters upar;
  upar.Add("nb", finitmap["n_b"],  sqrt(finitmap["n_b"])/10., finitmap["n_b"] - sqrt(finitmap["n_b"]),finitmap["n_b"] + sqrt(finitmap["n_b"]));
  upar.Add("nq", finitmap["n_cl"], sqrt(finitmap["n_cl"])/10.,finitmap["n_cl"] - sqrt(finitmap["n_cl"]),finitmap["n_cl"] + sqrt(finitmap["n_cl"]));
  upar.Add("pb", finitmap["p_b"],   sqrt(finitmap["p_b"])/10.,finitmap["p_b"] - sqrt(finitmap["p_b"]),finitmap["p_b"] + sqrt(finitmap["p_b"]));
  upar.Add("pq", finitmap["p_cl"], sqrt(finitmap["p_cl"])/10.,finitmap["p_cl"] - sqrt(finitmap["p_cl"]),finitmap["p_cl"] + sqrt(finitmap["p_cl"]));
  upar.Add("Emu", finitmap["effMu_b"], 0.001, finitmap["effMu_b"] , finitmap["effMu_b"]+0.1); 
  upar.Add("Rmu", finitmap["effMu_cl"], 0.001, finitmap["effMu_cl"] , finitmap["effMu_cl"]+0.2); 
  upar.Add("Esvx", finitmap["effTag_b"], 0.001, finitmap["effTag_b"]+0.05 , finitmap["effTag_b"]+0.2);
  upar.Add("Rsvx", finitmap["effTag_cl"], 0.001,finitmap["effTag_cl"]+0.05 , finitmap["effTag_cl"]+0.2);
  
  MnMigrad migrad(*thePDF, upar);
  
  FunctionMinimum fmin = migrad(20000,0.001);//10000,0.01);//migrad(maxFCN,tolerance);
  
  if(!fmin.IsValid()) {
    //try with higher strategy
    //std::cout<<"FM is invalid, try with strategy = 2."<<std::endl;
    MnMigrad migrad2(*thePDF, fmin.UserState(), MnStrategy(2));
    fmin = migrad2();
  } 
  if(!fmin.IsValid()) { 
    //try with higher strategy 
    //std::cout<<"FM is invalid, try with strategy = 2."<<std::endl; 
    MnMigrad migrad3(*thePDF, fmin.UserState(), MnStrategy(2)); 
    fmin = migrad3(); 
  }  

  ff_minimum = fmin.Fval();

  std::cout << fmin << std::endl;
  //std::cout << ff_minimum << std::endl;
  
  //if ( !isnan(ff_minimum) && !isinf(ff_minimum) ) {
  //if ( fmin.IsValid() ) {
  if ( true ) {
    fsolution["effTag_b"] = fmin.UserState().Value("Esvx");//fmin.Parameters().Vec()(6);
	fsolution["effTag_cl"] = fmin.UserState().Value("Rsvx");//fmin.Parameters().Vec()(7);  
	fsolution["effMu_b"] = fmin.UserState().Value("Emu");//fmin.Parameters().Vec()(4);  
	fsolution["effMu_cl"] =fmin.UserState().Value("Rmu");// fmin.Parameters().Vec()(5);  
	fsolution["n_b"] = fmin.UserState().Value("nb");//fmin.Parameters().Vec()(0); 
	fsolution["n_cl"] = fmin.UserState().Value("nq");//fmin.Parameters().Vec()(1);  
	fsolution["p_b"] = fmin.UserState().Value("pb");//fmin.Parameters().Vec()(2);
	fsolution["p_cl"] = fmin.UserState().Value("pq");//fmin.Parameters().Vec()(3); 
  
	  fsolutionErr["effTag_b"] = TMath::Sqrt(fmin.UserState().Error("Esvx"));//fmin.Error().Matrix()(6,6);  
	  fsolutionErr["effTag_cl"] = TMath::Sqrt(fmin.UserState().Error("Rsvx"));//fmin.Error().Matrix()(7,7);
	  fsolutionErr["effMu_b"] = TMath::Sqrt(fmin.UserState().Error("Emu"));//fmin.Error().Matrix()(4,4);
	  fsolutionErr["effMu_cl"] = TMath::Sqrt(fmin.UserState().Error("Rmu"));//fmin.Error().Matrix()(5,5);
	  fsolutionErr["n_b"] = TMath::Sqrt(fmin.UserState().Error("nb"));//fmin.Error().Matrix()(0,0);
	  fsolutionErr["n_cl"] = TMath::Sqrt(fmin.UserState().Error("nq"));//fmin.Error().Matrix()(1,1);
	  fsolutionErr["p_b"] = TMath::Sqrt(fmin.UserState().Error("pb"));//fmin.Error().Matrix()(2,2);
	  fsolutionErr["p_cl"] = TMath::Sqrt(fmin.UserState().Error("pq"));//fmin.Error().Matrix()(3,3);
  } else {
	  fsolution["effTag_b"] = 0.;
	  fsolution["effTag_cl"] = 0.;
	  fsolution["effMu_b"] = 0.;
	  fsolution["effMu_cl"] = 0.;
	  fsolution["n_b"] = 0.;
	  fsolution["n_cl"] = 0.;
	  fsolution["p_b"] = 0.;
	  fsolution["p_cl"] = 0.;
  
	  fsolutionErr["effTag_b"] = 0.;
	  fsolutionErr["effTag_cl"] =0.;
	  fsolutionErr["effMu_b"] = 0.;
	  fsolutionErr["effMu_cl"] =0.;
	  fsolutionErr["n_b"] = 0.;
	  fsolutionErr["n_cl"] =0.;
	  fsolutionErr["p_b"] = 0.;
	  fsolutionErr["p_cl"] =0.;
	  
  }
  
}


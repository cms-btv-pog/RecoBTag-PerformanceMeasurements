#ifndef S8bPerformance_h
#define S8bPerformance_h
/** \class S8bPerformance
 *
 * Analyze ROOT files produced by analyzer and create plots
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: S8bPerformance.h,v 1.1 2007/10/17 03:24:36 yumiceva Exp $
 *
 */

#include "TGraph.h"

class S8bPerformance {

  public:
	S8bPerformance() {};
	void Set(std::string name) {
		fname = name;
		fNcuts = 40;
		b_all = c_all = udsg_all = 0;
		for( int i=0; i<fNcuts; ++i) {
		  b_tagged[i] = 0;
		  c_tagged[i] = 0;
		  udsg_tagged[i] = 0;
		}
		//h1["b"+fname] = new TH1F("b"+TString(fname),"b efficiency",40,0.,0.8);
		//h1["c"+fname] = new TH1F("c"+TString(fname),"c efficiency",40,0.,0.8);
		//h1["udsg"+fname] = new TH1F("udsg"+TString(fname),"udsg efficiency",40,0.,0.8);
		//g1["b"+fname].Set(fNcuts);//		b_vsbgraph.Set(fNcuts);
		//g1["c"+fname].Set(fNcuts);//		b_vscgraph.Set(fNcuts);
		//g1["udsg"+fname].Set(fNcuts);//		b_vsudsggraph.Set(fNcuts);
	};
	void SetMinDiscriminator( double value) { fMinDisc = value; };
	void SetMaxDiscriminator( double value) { fMaxDisc = value; };
	void SetNcuts( int value ) { 
	  fNcuts = value;
	  for( int i=0; i<fNcuts; ++i) {
	    b_tagged[i] = 0;
	    c_tagged[i] = 0;
	    udsg_tagged[i] = 0;
	  }
	};
	void Add( double discriminator, int flavor ) {

		fdisc = discriminator;
		fflavor = flavor;
		double binwidth = (fMaxDisc-fMinDisc)/((double) fNcuts);
		if (fflavor==5) b_all++;
		if (fflavor==4) c_all++;
		if ((fflavor>0 && fflavor<4)||(fflavor==21)) udsg_all++;
		
		for( int i=0; i<fNcuts; ++i ) {
			double cut = i*binwidth + fMinDisc;
			if (cut > fMaxDisc) cut = fMaxDisc;

			if ( fdisc > cut ) {
				if ( fflavor==5 ) b_tagged[i]++;
				if ( fflavor==4 ) c_tagged[i]++;
				if ((fflavor>0 && fflavor<4)||(fflavor==21)) udsg_tagged[i]++;
			}
		}
	};

	void Eval() {
		int ip = 0;
		for(std::map<int,double>::const_iterator im=b_tagged.begin(); im!=b_tagged.end(); ++im) {

			double eff = (im->second)/((double)b_all);
			double err = sqrt(eff*(1-eff)/b_all);
			b_eff[ip] = eff;
			b_effErr[ip] = err;
			//g1["b"+fname].SetPoint(ip,eff,eff);
			//h1["b"+fname]->Fill(eff);
			ip++;
		}
		ip=0;
		for(std::map<int,double>::const_iterator im=c_tagged.begin(); im!=c_tagged.end(); ++im) {
		  
		  double eff = (im->second)/((double)c_all);
		  double err = sqrt(eff*(1-eff)/c_all);
		  c_eff[ip] = eff;
		  c_effErr[ip] = err;
		  //g1["c"+fname].SetPoint(ip,b_tagged[ip],eff);
			//h1["c"+fname]->Fill(eff);
			ip++;
		}
		ip=0;
		for(std::map<int,double>::const_iterator im=udsg_tagged.begin(); im!=udsg_tagged.end(); ++im) {

			double eff = (im->second)/((double)udsg_all);
			double err = sqrt(eff*(1-eff)/udsg_all);
			udsg_eff[ip] = eff;
			udsg_effErr[ip] = err;
			//g1["udsg"+fname].SetPoint(ip,b_tagged[ip],eff);
			//h1["udsg"+fname]->Fill(eff);
			ip++;
		}
				
	};
	std::map< int,double> GetMap(TString option="b") {
	  if (option=="b") return b_eff;
	  if (option=="c") return c_eff;
	  if (option=="udsg") return udsg_eff;
	  if (option=="bErr") return b_effErr;
	  if (option=="cErr") return c_effErr;
	  if (option=="udsgErr") return udsg_effErr;
	  //return NULL;
	};

	Double_t *GetArray(TString option="b") {
	  std::map<int,double> amap = GetMap(option);
	  Double_t *array;//[(int)amap.size()];
	  int ip=0;
	  for(std::map<int,double>::const_iterator im=amap.begin(); im!=amap.end(); ++im) {
	    array[ip] = im->second;
	  }
	  return array;
	  
	};

  private:
	double fdisc;
	int    fflavor;
	double fMinDisc;
	double fMaxDisc;
	int    fNcuts;
	std::string fname;
	//std::map< std::string, TH1* > h1;
	//std::map< std::string, TGraph > g1;
	std::map< int, double > b_eff;
	std::map< int, double > c_eff;
	std::map< int, double > udsg_eff;
	std::map< int, double > b_effErr;
	std::map< int, double > c_effErr;
	std::map< int, double > udsg_effErr;

	std::map< int, double > b_tagged;
	std::map< int, double > c_tagged;
	std::map< int, double > udsg_tagged;
	double b_all;
	double c_all;
	double udsg_all;
	
	
	
};

#endif


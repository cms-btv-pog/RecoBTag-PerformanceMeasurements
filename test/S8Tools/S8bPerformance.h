#ifndef S8bPerformance_h
#define S8bPerformance_h
/** \class S8bPerformance
 *
 * Analyze ROOT files produced by analyzer and create plots
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: S8bPerformance.h,v 1.4 2007/10/16 18:27:09 yumiceva Exp $
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
		//h1["b"+fname] = new TH1F("b"+TString(fname),"b efficiency",40,0.,0.8);
		//h1["c"+fname] = new TH1F("c"+TString(fname),"c efficiency",40,0.,0.8);
		//h1["udsg"+fname] = new TH1F("udsg"+TString(fname),"udsg efficiency",40,0.,0.8);
		g1["b"+fname].Set(fNcuts);//		b_vsbgraph.Set(fNcuts);
		g1["c"+fname].Set(fNcuts);//		b_vscgraph.Set(fNcuts);
		g1["udsg"+fname].Set(fNcuts);//		b_vsudsggraph.Set(fNcuts);
	};
	void SetMinDiscriminator( double value) { fMinDisc = value; };
	void SetMaxDiscriminator( double value) { fMaxDisc = value; };
	void SetNcuts( int value ) { fNcuts = value; };
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

	std::map<std::string, TGraph> GetMap() {

		int ip = 0;
		for(std::map<int,double>::const_iterator i=b_tagged.begin(); i!=b_tagged.end(); ++i) {

			double eff = (i->second)/((double)b_all);
			g1["b"+fname].SetPoint(ip,eff,eff);
			//h1["b"+fname]->Fill(eff);
			ip++;
		}
		ip=0;
		for(std::map<int,double>::const_iterator i=c_tagged.begin(); i!=c_tagged.end(); ++i) {

			double eff = (i->second)/((double)c_all);
			g1["c"+fname].SetPoint(ip,b_tagged[ip],eff);
			//h1["c"+fname]->Fill(eff);
			ip++;
		}
		ip=0;
		for(std::map<int,double>::const_iterator i=udsg_tagged.begin(); i!=udsg_tagged.end(); ++i) {

			double eff = (i->second)/((double)udsg_all);
			g1["udsg"+fname].SetPoint(ip,b_tagged[ip],eff);
			//h1["udsg"+fname]->Fill(eff);
			ip++;
		}
		//b_vsbgraph.Print();
		//b_vscgraph.Print();
		//b_vsudsggraph.Print();
		
		//h1["b"+fname] = new TH1F("b"+TString(fname),"b efficiency",fNcuts,0.,0.8);;//b_vsbgraph.GetHistogram();
		//h1["c"+fname] = b_vscgraph.GetHistogram();
		//h1["udsg"+fname] = b_vsudsggraph.GetHistogram();
		
		return g1;
		
	};
	
  private:
	double fdisc;
	int    fflavor;
	double fMinDisc;
	double fMaxDisc;
	int    fNcuts;
	std::string fname;
	//std::map< std::string, TH1* > h1;
	std::map< std::string, TGraph > g1;
	std::map< int, double > b_tagged;
	std::map< int, double > c_tagged;
	std::map< int, double > udsg_tagged;
	double b_all;
	double c_all;
	double udsg_all;
	//TGraph b_vsbgraph;
	//TGraph b_vscgraph;
	//TGraph b_vsudsggraph;
	
	
};

#endif


#ifndef S8Plotter_h
#define S8Plotter_h
/** \class S8Plotter
 *  
 * Analyze ROOT files produced by analyzer and create plots
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: S8Plotter.h,v 1.15 2008/09/12 16:41:32 bazterra Exp $
 *
 */

#include "TROOT.h"
#include "TAxis.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"

#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>

#ifdef NOSCRAMV
#include "BTagEvent.h"
#include "BTagLeptonEvent.h"
#include "S8bPerformance.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"
#include "S8bPerformance.h"
#endif


class S8Plotter {

  public:

	S8Plotter();
	S8Plotter(TString filename);
	virtual ~S8Plotter();

	virtual Int_t  GetEntry(Long64_t entry);
	virtual void     Init();
	virtual void     Loop();
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Show(Long64_t entry = -1);
	virtual Int_t    Cut(Long64_t entry);
	virtual Bool_t   Notify();
	//virtual std::string itoa(int i);
	//virtual std::string ftoa(float i);
	void Add(TString filename);
	void Verbose(bool option=true) { fVerbose = option; };
	void Book();
	void SetTagger(TString name) { ftagger = name; };
	void SetTaggerLevel(TString name) { flevel = name; };
	void SetAwayTagger(TString name) { fAwaytagger = name; };
	void SetAwayTaggerLevel(TString name) { fAwaylevel = name; };
	void SampleName(TString sample) {
		fsamplename = sample;
	};
	//void SetPtMin(double value) { fMinPt = value; };
	//void SetPtMax(double value) { fMaxPt = value; };
	void FillHistos(std::string type, TLorentzVector p4Jet, double ptrel, int flavor, bool tagged);
	void Print(std::string extension="png",std::string tag="") {
		if ( tag != "" ) tag = "_"+tag;
		
		for(std::map<std::string,TCanvas*>::const_iterator icv=cv_map.begin(); icv!=cv_map.end(); ++icv){

			std::string tmpname = icv->first;
			TCanvas *acv = icv->second;
			acv->Print(TString(tmpname+tag+"."+extension));
		}		
	};
	
	void SaveToFile(TString filename="S8_plots.root") {

		//if (!foutfile) {
			
		foutfile = new TFile(filename,"RECREATE");
		for(std::map<std::string,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih){
			TH1 *htemp = ih->second;
			htemp->Write();
		}
		for(std::map<std::string,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih){
			TH2 *htemp = ih->second;
			htemp->Write();
		}
		//for(std::map<std::string,TGraph >::const_iterator ih=g1.begin(); ih!=g1.end(); ++ih){
		//TGraph gtemp = ih->second;
		//gtemp.Write();
		//}
		
		foutfile->Write();
		foutfile->Close();
		std::cout << " File " << filename << " saved." << std::endl;
	};
	//bool IsGoodMuon() {
	//};
	void Clean() {
		for(std::map<std::string,TH1* >::const_iterator ih=h1.begin(); ih!=h1.end(); ++ih){
			delete ih->second;
		}
		for(std::map<std::string,TH2* >::const_iterator ih=h2.begin(); ih!=h2.end(); ++ih){
			delete ih->second;
		}
		delete gTC2_b;
		delete gTC2_c;
		delete gTC2_udsg;
		delete gTC3_b;
		delete gTC3_c;
		delete gTC3_udsg;
		delete gTP_b;
		delete gTP_c;
		delete gTP_udsg;
		//delete multiTC2;
		delete multiTC3;
		delete multiTP;
	}

	double EffErr(double n, double p) {
		double w= n/p;
		return sqrt(w*(1-w)/p);
	}
	void SetAxis(TString name, TAxis axis) {
	  if (name == "Pt") {
		  fJetPtAxis.Set(axis.GetNbins(),(axis.GetXbins())->GetArray());
		  fMinPt = fJetPtAxis.GetXmin();
		  fMaxPt = fJetPtAxis.GetXmax();
	  }
	  if (name == "Eta") {
		  fJetEtaAxis.Set(axis.GetNbins(),(axis.GetXbins())->GetArray());
		  fMinEta = fJetEtaAxis.GetXmin();
		  fMaxEta = fJetEtaAxis.GetXmax();
	  }
	  if (name == "CorrPt") {
		  fCorrPtAxis.Set(axis.GetNbins(),(axis.GetXbins())->GetArray());
		  fMinCorrPt = fCorrPtAxis.GetXmin();
		  fMaxCorrPt = fCorrPtAxis.GetXmax();
	  }
	  if (name == "CorrEta") {
		  fCorrEtaAxis.Set(axis.GetNbins(),(axis.GetXbins())->GetArray());
		  fMinCorrEta = fCorrEtaAxis.GetXmin();
		  fMaxCorrEta = fCorrEtaAxis.GetXmax();
	  }
	}
	void PrintInfo();
	void PrintBins(TString option="Pt") {
	  if (option =="Pt" || option =="Eta") {
	    const Double_t *ax = 0;
	    int nbins = 0;
	    if (option=="Pt") { 
	      ax = (fJetPtAxis.GetXbins())->GetArray();
	      nbins = fJetPtAxis.GetNbins();
	    }
	    if (option=="Eta") {
	      ax = (fJetEtaAxis.GetXbins())->GetArray();
	      nbins = fJetEtaAxis.GetNbins();
	    }
		if (option=="CorrPt") { 
	      ax = (fCorrPtAxis.GetXbins())->GetArray();
	      nbins = fCorrPtAxis.GetNbins();
	    }
		if (option=="CorrEta") {
	      ax = (fCorrEtaAxis.GetXbins())->GetArray();
	      nbins = fCorrEtaAxis.GetNbins();
	    }
		
	    std::cout << option << " binning =" << std::endl;
	    std::cout << "{";
	    for( int i=0; i<nbins; ++i) {
	      std::cout << ax[i];
	      if (i!=nbins-1) std::cout << ",";
	    }
	    std::cout << "}" << std::endl;
	  } else { std::cout << " don't know about " << option << std::endl;}
	}

  private:
	
	BTagEvent *fS8evt;
	
	TChain            *fChain;
	Int_t             fCurrent;
	TFile            *ffile;
	TFile            *foutfile;
	TString           fsamplename;
	bool              fVerbose;
	TString           ftagger;
	TString           flevel;
	TString           fAwaytagger;
	TString           fAwaylevel;
	double            fMinPt;
	double            fMaxPt;
	double            fMinEta;
	double            fMaxEta;
	double            fMinCorrPt;
	double            fMaxCorrPt;
	double            fMinCorrEta;
	double            fMaxCorrEta;
	
	TAxis             fJetPtAxis;
	TAxis             fJetEtaAxis;
	TAxis             fCorrPtAxis;
	TAxis             fCorrEtaAxis;
	
	std::map<std::string, TCanvas*> cv_map;
	std::map<std::string, TH1*> h1;
	std::map<std::string, TH2*> h2;
	std::map<std::string, TGraph> gg1;
	
	std::vector < std::string > quark_label;
	std::map<std::string, int>  quark_color;
	std::map<std::string, std::string > cut_label;
	std::map< TString, float > fTrackCountingMap; // discriminator cut and level
	std::map< TString, float > fTrackProbabilityMap; // discriminator cut and level
	std::map< TString, float > fbTaggerMap; // point to the selected tagger
	std::map< TString, float > fbAwayTaggerMap; // point to the selected tagger
	
	S8bPerformance fperformanceTC2;
	S8bPerformance fperformanceTC3;
	S8bPerformance fperformanceTP;
	
	TGraphErrors *gTC2_b;	
	TGraphErrors *gTC2_c;	
	TGraphErrors *gTC2_udsg;	
	TGraphErrors *gTC3_b;	
	TGraphErrors *gTC3_c;	
	TGraphErrors *gTC3_udsg;	
	TGraphErrors *gTP_b;	
	TGraphErrors *gTP_c;	
	TGraphErrors *gTP_udsg;	
	TMultiGraph *multiTC2;
	TMultiGraph *multiTC3;
	TMultiGraph *multiTP;

};

#endif

#ifdef S8Plotter_cxx
S8Plotter::S8Plotter() {

	S8Plotter("");
	
}

S8Plotter::S8Plotter(TString filename)
{

	// defaults
	fVerbose = false;
	ftagger = "TrackCounting";
	flevel  = "Loose";
	fAwaytagger = "TrackCounting";
	fAwaylevel = "Loose";

	fTrackCountingMap["Loose"]  = 2.0; // use TC2:high eff.
	fTrackCountingMap["Medium"] = 4.2; // use TC2:high eff.
	fTrackCountingMap["Tight"]  = 4.1;// use TC3:high purity
 
	fTrackProbabilityMap["Loose"] = 0.24;
	fTrackProbabilityMap["Medium"] = 0.49;
	fTrackProbabilityMap["Tight"] = 0.74;

	// default Binning
//	const int nptarray = 11;
//	const int nptarray = 6;
//	const int nptarray = 2;
	const int nptarray = 4;
	const int netaarray = 10;
	const int ncorrptarray = 5;
	const int ncorretaarray = 5;
//	Double_t jetptbins[nptarray] = {20.,30., 40., 50., 60., 70., 80, 90., 100., 120., 140., 230.};
//{20.,30., 40., 50., 60., 70., 80, 90., 100., 120., 140., 160., 180., 230.};
//	Double_t jetptbins[nptarray] = {30., 40., 50., 60., 70., 80, 90., 100., 120., 140., 230.};
//	Double_t jetptbins[nptarray] = {30., 50., 70., 90., 120., 230.};
	Double_t jetptbins[nptarray] = {30., 70., 120.,230.};
//	Double_t jetptbins[nptarray] = {30., 230.};

	Double_t jetetabins[netaarray] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5};
	Double_t corrptbins[ncorrptarray] = {20.,40.,60.,80.,230.};
	Double_t corretabins[ncorrptarray] = {0.,0.5,1.,1.5,2.5};
	fJetPtAxis.Set(nptarray-1, jetptbins );
	fJetEtaAxis.Set(netaarray-1, jetetabins );
	fCorrPtAxis.Set(ncorrptarray-1, corrptbins );
	fCorrEtaAxis.Set(ncorretaarray-1, corretabins );
	
	fMinPt = fJetPtAxis.GetXmin();
	fMaxPt = fJetPtAxis.GetXmax();
	fMinEta = fJetEtaAxis.GetXmin();
	fMaxEta = fJetEtaAxis.GetXmax();
	fMinCorrPt = fCorrPtAxis.GetXmin();
	fMaxCorrPt = fCorrPtAxis.GetXmax();
	fMinCorrEta = fCorrEtaAxis.GetXmin();
	fMaxCorrEta = fCorrEtaAxis.GetXmax();
	
	// eff plots for all discriminator values
	fperformanceTC2.Set("TC2");
	fperformanceTC3.Set("TC3");
	fperformanceTP.Set("TP");
	fperformanceTC2.SetMinDiscriminator(-1);
	fperformanceTC2.SetMaxDiscriminator(15);
	fperformanceTC3.SetMinDiscriminator(-1);
	fperformanceTC3.SetMaxDiscriminator(15);
	fperformanceTP.SetMinDiscriminator(0);
	fperformanceTP.SetMaxDiscriminator(1);
	
	
	// chain
	fChain = new TChain("summary");
	if ( filename != "" ) fChain->Add(filename);

	// event container
	fS8evt = new BTagEvent();
			
}
void S8Plotter::Add(TString filename)
{
	fChain->Add(filename);
}

S8Plotter::~S8Plotter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

std::string itoa(int i) {
  char temp[20];
  sprintf(temp,"%d",i);
  return((std::string)temp);
}

std::string ftoa(float i, const char *format="0") {
  char temp[20];
  if ( format == "0" ) sprintf(temp,"%.2f",i);
  else sprintf(temp,format,i);
  
  return((std::string)temp);
}

std::string GetStringFlavour(int i) {
	if ( i == 5 ) { return "_b"; }
	else if ( i == 4 ) { return "_c"; }
	else if ( i>0 && i<=3 ) { return "_udsg"; }
	else if ( i == 21 ) { return "_udsg"; }
	else if ( i == 0 ) { return "ambiguous";}
	else { 
		//cout << " not defined flavour " << i << endl;
		return "";
	}
}

void NormalizeHistograms(std::vector< TH1* > hist) {

	for ( std::vector< TH1* >::size_type ihist = 0; ihist != hist.size() ; ++ihist ) {
		hist[ihist]->Scale(1./hist[ihist]->Integral());
	}
	
}

Int_t S8Plotter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t S8Plotter::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void S8Plotter::Init()
{
	// check b-tagging
	if ( ftagger == "TrackCounting" ) 		  fbTaggerMap = fTrackCountingMap;
	else if ( ftagger == "TrackProbability" ) fbTaggerMap = fTrackProbabilityMap;
	else { 
	  std::cout << " No b tagger " << ftagger << " available, options are: TrackCounting, TrackProbability" << std::endl;
	  exit(1);//gApplication->Terminate(); 
	}
	if ( fAwaytagger == "TrackCounting" ) 		  fbAwayTaggerMap = fTrackCountingMap;
	else if ( fAwaytagger == "TrackProbability" ) fbAwayTaggerMap = fTrackProbabilityMap;
	else { 
	  std::cout << " No b tagger, for away jet, " << ftagger << " available, options are: TrackCounting, TrackProbability" << std::endl;
	  exit(1);//gApplication->Terminate(); 
	}

	
	//if (!tree) return;

	//fChain = tree;
	fCurrent = -1;
	
	fChain->SetBranchAddress("s8.",&fS8evt);
	//fChain->SetBranchAddress("summaryKVF.",&fBTagSummary[1]);
	//fChain->SetBranchAddress("summaryTKF.",&fBTagSummary[2]);

	//cout << "Init done" << endl;
}

void S8Plotter::PrintInfo() {

	std::cout << " Tagger: " << ftagger << std::endl;
	std::cout << " Level:  " << flevel << " (discriminator>" << fbTaggerMap[flevel] << ")" << std::endl;
	std::cout << " Tagger for away-jet: " << fAwaytagger << std::endl;
	std::cout << " Level of Tagger for away-jet:  " << fAwaylevel << " (discriminator>" << fbTaggerMap[fAwaylevel] << ")" << std::endl;
	
}


Bool_t S8Plotter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void S8Plotter::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t S8Plotter::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef S8Plotter_cxx

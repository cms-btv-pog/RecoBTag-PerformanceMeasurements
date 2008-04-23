#ifndef S8Solver_h
#define S8Solver_h

#include "TROOT.h"
#include "TNamed.h"
#include "TString.h"
#include "TArrayD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

#include <string>
#include <map>


class S8Solver {

  public:
        S8Solver();
		//S8Solver(std::string name);
		virtual ~S8Solver(){};
		void LoadHistos();
		void SetName(TString value) {fthename= value;}
		void SetData( TString filename ) { finputFile = new TFile(filename); }
		void SetCorrData( TString filename ) { finputCorrFile = new TFile(filename); }
		void SetPtrelCut ( double value ) { fminPtrel = value; }
		void SetPtrelMaxCut ( double value ) { fMaxPtrel = value; }
		void SetPtFits() { fcategory = "pT"; }
		void SetEtaFits() { fcategory = "Eta"; }
		void SetRebin( bool option) { frebin = option; }
		void SetMinimum( double value ) { fMin = value; }
		void SetMaximum( double value ) { fMax = value; }
		void RecalculateFactors( bool option ) { frecalculateFactors = option;}
		void SetAlphaConstant(bool option) { fAlphaConst = option; }
		void SetBetaConstant(bool option) { fBetaConst = option; }
		void SetKappabConstant(bool option) { fKappabConst = option; }
		void SetKappaclConstant(bool option) { fKappaclConst = option; }
		void SetAlphaFactor(double value) { fAlphaf = value; }
		void SetBetaFactor(double value) { fBetaf = value; }
		void SetKappabFactor(double value) { fKappabf = value; }
		void SetKappaclFactor(double value) { fKappaclf = value; }
		void SetDeltaConstant(bool option) { fDeltaConst = option; }
		void SetGammaConstant(bool option) { fGammaConst = option; }
		void SetMethod( TString option ) { fmethod = option; };
		void Solve();
        void Verbose(bool option) { fVerbose = option; };
		void PrintData(TString option="");
		void Draw(int maxNbins=0);
		void Print(TString extension="eps");
		void Save(TString filename="solver.root") {

			TFile *ofile = new TFile(filename,"RECREATE");

			// true efficiency
			feffTag_b->Write();
			feffTag_cl->Write();
			feffmu_b->Write();
			feffmu_cl->Write();
			// input
			fnHisto->Write();
			fpHisto->Write();
			fnHistoMu->Write();
			fpHistoMu->Write();
			fnHistoSvx->Write();
			fpHistoSvx->Write();
			fnHistoAll->Write();
			fpHistoAll->Write();
			fh_kb->Write();
			fh_kcl->Write();
			fh_alpha->Write();
			fh_beta->Write();
			fh_delta->Write();
			fh_gamma->Write();
			// output
			geffTag_b->Write();
			gS8effTag_b->Write();
			geffmu_b->Write();
			gS8effmu_b->Write();
				
			//ofile->Write();
			ofile->Close();
			delete ofile;
		};
		
		
  protected:
		void GetInput();
		
  private:
        bool fVerbose;
		TFile *finputFile;
		TFile *finputCorrFile;
		TString fcategory;
		TString fmethod;
		TString fthename;
		double fMin;
		double fMax;
		bool fAlphaConst;
		bool fBetaConst;
		bool fKappabConst;
		bool fKappaclConst;
		bool fDeltaConst;
		bool fGammaConst;
		bool fisCorrFile;
		double fAlphaf;
		double fBetaf;
		double fKappabf;
		double fKappaclf;
		bool frebin;
		bool frecalculateFactors;
		double fminPtrel;
		double fMaxPtrel;
		std::map< TString, double > TotalInput;
		std::map< int, std::map< TString, double> > BinnedInput;
		std::map< TString, double > fTotalSolution;
		std::map< TString, double > fTotalSolutionErr;
		std::map< int, std::map< TString, double> > fBinnedSolution;
		std::map< int, std::map< TString, double> > fBinnedSolutionErr;

		std::map<TString, TCanvas*> cv_map;
		
		TH2F* fnHistoBase;
		TH2F* fpHistoBase;
		TH2F* fnSvxHistoBase;
		TH2F* fpSvxHistoBase;
		TH1D* fnHisto;
		TH1D* fpHisto;
		TH1D* fnHistoMu;
		TH1D* fpHistoMu;
		TH1D* fnHistoSvx;
		TH1D* fpHistoSvx;
		TH1D* fnHistoAll;
		TH1D* fpHistoAll;
		TH1D* fh_kb;
		TH1D* fh_kcl;
		TH1D* fh_alpha;
		TH1D* fh_beta;
		TH1D* fh_delta;
		TH1D* fh_gamma;
		TH1D* feffTag_b;
		TH1D* feffTag_cl;
		TH1D* feffmu_b;
		TH1D* feffmu_cl;
		TH1D* feffTagmu_b;
		TH1D* feffTagmu_cl;

		TH1D *halljets_b;  
		TH1D *halljets_cl;  
		TH1D *htagjets_b ;  
		TH1D *htagjets_cl;  
		TH1D *halljets_b_ptrel;  
		TH1D *halljets_cl_ptrel;  
		TH1D *htagjets_b_ptrel;  
		TH1D *htagjets_cl_ptrel;  
		TH1D *halloppjets_b;  
		TH1D *halloppjets_cl;  
		TH1D *halloppjets_b_ptrel;  
		TH1D *halloppjets_cl_ptrel; 
		TH1D *htagoppjets_b_ptrel;  
		TH1D *htagoppjets_cl_ptrel; 
		TH1D *htagoppjets_b;  
		TH1D *htagoppjets_cl; 
		
		TGraphErrors *geffTag_b;
		TGraphErrors *gS8effTag_b;
		TGraphErrors *geffmu_b;
		TGraphErrors *gS8effmu_b;
		TGraphErrors *ginput_n;
		TGraphErrors *ginput_nmu;
		TGraphErrors *ginput_ntag;
		TGraphErrors *ginput_ntagmu;
		TGraphErrors *ginput_p;
		TGraphErrors *ginput_pmu;
		TGraphErrors *ginput_ptag;
		TGraphErrors *ginput_ptagmu; 
		
		ClassDef(S8Solver,1);


};

#endif

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
#include <utility>

#include "S8NumericInput.h"

class S8Solver
{
  public:
        S8Solver();
		//S8Solver(std::string name);
		virtual ~S8Solver(){};
		void Clear();
		void LoadHistos();
		void SetName(TString value) {fthename= value;}
		void SetData( TString filename ) { finputFile = new TFile(filename); }
		void SetCorrData( TString filename ) { finputCorrFile = new TFile(filename); }
		void SetPtrelCut ( double value ) { fminPtrel = value; }
		void SetPtrelMaxCut ( double value ) { fMaxPtrel = value; }
		void SetPtFits() { fcategory = "pT"; }
		void SetEtaFits() { fcategory = "eta"; }
		void SetRebin( int nbins ) { frebin = true; fnbins = nbins; }
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
		void SetDeltaFactor(double value) { fDeltaf = value; }
		void SetGammaFactor(double value) { fGammaf = value; }
		void SetDeltaConstant(bool option) { fDeltaConst = option; }
		void SetGammaConstant(bool option) { fGammaConst = option; }
		void SetMethod( TString option ) { fmethod = option; }
		void UseMCTrue(bool option) { fusemctrue = option; }
		void Solve();
		void SetSolution( int bin, int solution) {
			// bin 0 corresponds to average solution
				fPickSolutionMap[bin] = solution;
		}

        void setDoBinnedSolution(bool flag) { _doBinnedSolution = flag; }


        void Verbose(bool option) { fVerbose = option; };
		void PrintData(TString option="");
		void DumpTable(std::string filename="table.txt");
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
			geffTag_cl->Write();
			gS8effTag_cl->Write();
			geffmu_cl->Write();
			gS8effmu_cl->Write();
				
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
        bool fusemctrue;

		double fAlphaf;
		double fBetaf;
		double fKappabf;
		double fKappaclf;
		double fDeltaf;
		double fGammaf;
		bool frebin; int fnbins;
		bool frecalculateFactors;
        bool _doBinnedSolution;
		double fminPtrel;
		double fMaxPtrel;
		std::map< TString, double > TotalInput;
        std::map< TString, double > TotalInputErr;
		std::map< int, std::map< TString, double> > BinnedInput;
		std::map< TString, double > fTotalSolution;
		std::map< TString, double > fTotalSolutionErr;
		std::map< int, int > fPickSolutionMap;
		std::map< int, std::map< TString, double> > fBinnedSolution;
		std::map< int, std::map< TString, double> > fBinnedSolutionErr;

		std::map<TString, TCanvas*> cv_map;

		//int fPickBin; int fPicknSol;

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
		
		TH1D *b_halljets_ptrel ;       
		TH1D *cl_halljets_ptrel ;      
		TH1D *b_halljets_tagged;       
		TH1D *cl_halljets_tagged ;     
		TH1D *b_halljets_ptreltagged ; 
		TH1D *cl_halljets_ptreltagged ;
		TH1D *b_halloppjets_tagged ;   
		TH1D *cl_halloppjets_tagged;   
		TH1D *b_halloppjets_ptrel ;    
		TH1D *cl_halloppjets_ptrel;    
		TH1D *b_halloppjets_ptreltagged ; 
		TH1D *cl_halloppjets_ptreltagged ;


		TGraphErrors *geffTag_b;
		TGraphErrors *gS8effTag_b;
		TGraphErrors *geffmu_b;
		TGraphErrors *gS8effmu_b;
		TGraphErrors *geffTag_cl;
		TGraphErrors *gS8effTag_cl;
		TGraphErrors *geffmu_cl;
		TGraphErrors *gS8effmu_cl;
		TGraphErrors *ginput_n;
		TGraphErrors *ginput_nmu;
		TGraphErrors *ginput_ntag;
		TGraphErrors *ginput_ntagmu;
		TGraphErrors *ginput_p;
		TGraphErrors *ginput_pmu;
		TGraphErrors *ginput_ptag;
		TGraphErrors *ginput_ptagmu; 

        TH1 *_averageResults[8];

        void calculateCorrCoefficients(NumericInput &numericInput);

        NumericInput              _totalInput;
        std::vector<NumericInput> _binnedInput;
		
		ClassDef(S8Solver,1);
};

Measurement calculateEfficiency(const TH1 *, const TH1 *);
Measurement calculateEfficiency(const Measurement &, const Measurement &);

Measurement integrate(const TH1 *);
Measurement measurementFromBin(const TH1 *, const int &);

#endif

#ifndef PtrelSolver_h
#define PtrelSolver_h

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"

#include "math.h"
#include "fstream.h"
#include "iomanip.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "vector.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"

#define PT_BASE  100
#define ETA_BASE 1

Double_t pdf(Double_t *xx, Double_t *par);
Double_t pdf1(Double_t *xx, Double_t *par);
Double_t combined_pdf(Double_t *xx, Double_t *par);

static double pthatbining[100];
static double etabining[100];

class PtrelSolver : public TObject {

 protected:

  Int_t pthat_bins;
  Int_t eta_bins;
  

  Double_t x_min, x_max;
  Double_t fit_min, fit_max;
  Int_t hist_bins;

  Int_t pt_threshold;
  Int_t pt_sumbins;
  Int_t eta_threshold;
  Int_t eta_sumbins;


  // define cuts
  TCut muon, jet;

  TCanvas *c1;
  std::fstream file;
  TF1     *func;
  TH1F    *hist;
  TFile   *inputfile;
  TFile   *mistagfile;
  TFile   *outfile;

  TLatex  *label;


  std::vector<std::vector<double> > *eff_table;
  std::vector<std::vector<double> > *pdfs_b;
  std::vector<std::vector<double> > *pdfs_c;
  TObjArray                          combined_pdfs;


  std::vector<std::vector<double> > *pdfs_b_tag;
  std::vector<std::vector<double> > *pdfs_c_tag;


  TObjArray                          combined_pdfs_tag;
  TObjArray                          combined_pdfs_sys;
  TObjArray                          combined_pdfs_sys_tag;

  TObjArray data, mc;  // store histograms.

  void init();


  int    getBin(double xx, double *binning);
  int    index( double pt, double eta, double others = 0);
  int    index( int ptnum, int etanum, int others=0);
  double effErr(double N1, double N1_err, double N2, double N2_err);


  TF1   *getAPdf(int ii);
  TF1   *getAPdf(int pt_bin, int eta_bin);
  TF1   *getPdfByIndex(TObjArray *list, int jj, const char *tag=0);
  TF1   *getPdfByIndex(int ii); // return the pdf before tagging
  TF1   *getPdfByIndex(int ii, const char *tag); // other methods
  TF1   *getTaggedPdfByIndex(int ii); // return the pdf after tagging


  // calculate the efficiency with a fit to the ptrel distributions.
  // there are several cases:
  // a) used the same templates before/after tagging.
  // b) used different templates before/after tagging.
  void effCal(TH1F *hist, TH1F *hist_tag, TF1  *pdf, std::vector<double> *eff);
  void effCal(TH1F *hist, TH1F *hist_tag, TF1  *pdf, TF1  *pdf_tag,  std::vector<double> *eff);
  void effCal(TH1F *hist, TH1F *hist_tag, TF1  *pdf, std::vector<double> *eff, const char *rootfilename);


 public:

  // Default constructor.
  PtrelSolver();
  PtrelSolver(double fitmin, double fitmax);
  PtrelSolver(double fitmin, double fitmax, int ptbins, int etabins);

  // Destructor.
  virtual ~PtrelSolver();

  
  void   setPtAverage( int threshold, int sum);
  void   setEtaAverage(int threshold, int sum);







  void readTemplates(const char *filename, std::vector<std::vector<double> > *parameters);
  void buildPdfs(TObjArray *combined, std::vector<std::vector<double> > *b, std::vector<std::vector<double> > *c, const char *tag = 0);
  TF1 *buildAPdf(int ii, std::vector<std::vector<double> > *b, std::vector<std::vector<double> > *c);


  void initPdfs(const char *b_pdf = "b_flavor.data", const char *c_pdf = "c_flavor.data");
  void initPdfs(const char *b_pdf,  const char *c_pdf, TObjArray *combined, const char *tag);


  void makeTemplates(int flavor, const char *inputfilename, const char *pdffilename, const char *rootfilename, const char *tag  ="", const char *pthist = "npT", const char *etahist="nEta", bool latex = false);

  void makeTemplates(int flavor = 5, const char *filename = "b_flavor.dat", const char *rootfilename = "b_pdfs.root" );

  

  TH1F *getData(double pt_min, double pt_max, double eta_min, double eta_max, TCut cut);
  void getMCEff(double pt_min, double pt_max, double eta_min, double eta_max, TCut tagger, double &eff, double &err);


  // build efficiecny basing give data & templates. 
  void estimate(TCut tagger, 
		std::vector<std::vector<double> > *eff, 
		std::vector<std::vector<double> > *b, 
		std::vector<std::vector<double> > *c,
		const char *rootfilename = "fit_result.root",
		const char *filename = "fit_result.dat");



  void makePlot(const char *epsname, TH1 *hist, TF1 *pdf);

  void makeEffTable(std::vector<std::vector<double> > *eff, const char *filename="eff.table");
  void makeEffHists(std::vector<std::vector<double> > *eff, const char *tag ="run1");
  void makeHistEPS(TObjArray &data, TObjArray &mc);
  void test(bool sys= false);
  void make(bool sys = false);
  void makettbar();


  TH1F *getMCeff(TFile *file, const char *hist); 
  void measure(const char *inputfilename, const char *outputfilename, const char *tag, const char *pthist, const char *etahist, bool sys);

  void measure(const char *inputfilename, const char *tag, const char *pthist="npT", const char *etahist= "nEta", bool sys= false);

  void counting(const char *tag, const char *outfilename, const char *inputfilename, const char *mistagfilename = 0);

  void Fit(TH1F *data, TF1 *pdf, std::vector<double> *num, std::vector<double> *num_err);
  


  TGraphErrors *checkLinearity2(TF1 *pdf, int steps, int total_events= 1000, bool verbose = true);
  TGraphErrors *checkPurity2(TF1 *pdf, int steps, int num_bs, bool verbose);


  ClassDef(PtrelSolver, 1)
};

#endif


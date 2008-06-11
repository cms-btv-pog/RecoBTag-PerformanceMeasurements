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
#include <fstream>
#include <iomanip>
#include "stdio.h"
#include "string"
#include "stdlib.h"
#include "vector"
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
  void effCal(TH1F *hist, TH1F *hist_tag, TF1  *pdf, TF1 *pdf_tag, std::vector<double> *eff, const char *rootfilename);

  void Fit(TH1F *data, TF1 *pdf, std::vector<double> *num, std::vector<double> *num_err);

  bool locateFile(const char *file); 

 public:

  // Default constructor.
  PtrelSolver();
  PtrelSolver(double fitmin, double fitmax);
  PtrelSolver(double fitmin, double fitmax, int ptbins, int etabins);

  // Destructor.
  virtual ~PtrelSolver();

  
  void   setPtAverage( int threshold, int sum);
  void   setEtaAverage(int threshold, int sum);
  TH1F  *getMCeff(TFile *file, const char *hist,     const char *tagger); 
  TH1F  *getMCeff(TFile *file, const char *dir, const char *hist, const char *tagger); 
  TH1F  *getMCeff(const char *filename, const char *dir, const char *hist, const char *tagger); 

  void counting(
    const char *,
    const char *,
    const char *
  );


  void counting(
    const char *, 
    const char *,
    const char *,
    const char *,
    const char *    
  );

  void makeAllTemplatesPerTag(const char *inputfilename, const char *dir, const char *sampletag, const char *tag, const char *outputdir, const char *versiontag, bool sys = false);
  void makeTemplates(const char *flavor, const char *sampletag, bool sum, const char *inputfilename, const char *dir, const char *tag, const char *thehist, int pdfbase, const char *outputdir, const char *versiontag="", bool sys=false, bool latex=false);


  //  void makeTemplates(int flavor, const char *inputfilename, const char *pdffilename, const char *rootfilename, const char *tag  ="", const char *pthist = "npT", const char *etahist="nEta", bool latex = false);

  void readTemplates(const char *filename, std::vector<std::vector<double> > *parameters);
  void buildPdfs(TObjArray *combined, std::vector<std::vector<double> > *b, std::vector<std::vector<double> > *c, const char *tag = 0);
  TF1 *buildAPdf(int ii, std::vector<std::vector<double> > *b, std::vector<std::vector<double> > *c);
  void makeAllPerFile(const char *datafile, const char *outputdir, const char *versiontag);


  void initPdfs(const char *b_pdf,  const char *c_pdf, TObjArray  *combined, const char *tag);
  void initPdfs(const char *b_pdf,  const char *c_pdf, const char *pdftag);


  void makePlot(const char *epsname, TH1 *hist, TF1 *pdf);

  void makeEffTable(std::vector<std::vector<double> > *eff, const char *filename="eff.table");
  void makeEffHists(std::vector<std::vector<double> > *eff, const char *tag ="run1");
  void makeHistEPS(TObjArray &data, TObjArray &mc);



  // for fitting
  void measure(const char *sampletag, const char *inputfilename, const char *dir, const char *outfilename, const char *tag, const char *thehistname, bool sys=false, const char *mcfilename=0, const char *mcdir=0);


  
  void measureByCounting(
    const char *,
    const char *,
    const char *,
    const char *,
    const char *,
    const char *
  );


  //  pdf data file: sampletag + favor+"_templates_" + tag + versiontag.
  bool initPdfsByTag(const char * sampletag, const char * tag, const char * pdfdir, const char *versiontag, bool sys);
  bool initPdfsByTag(const char * directory, const char * tag, const char * versiontag);

  void allEff(const char *inputfilename, const char *dir, const char *outfilename, const char *sampletag, const char *pdfdir, const char *tagger, const char *versiontag="", bool sys=false);


  // calibration
  TGraphErrors *checkLinearity2(TF1 *pdf, int steps, int total_events= 1000, bool verbose = true);
  TGraphErrors *checkPurity2(TF1 *pdf, int steps, int num_bs, bool verbose);

  ClassDef(PtrelSolver, 1)
};

#endif


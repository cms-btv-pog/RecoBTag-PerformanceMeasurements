
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2F.h"

#include "math.h"
#include "iostream.h"
#include "fstream.h"
#include "iomanip.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "vector.h"

#include "PtrelSolver.h"
#include "analysis.h"

/**************************************************************************
 *
 * templates: d * [x^a * e( bx^2) +c]
 * d is the normalization factor. 
 *
 **************************************************************************/
Double_t pdf(Double_t *xx, Double_t *par) {


  Double_t x = xx[0];

  return par[3] * (pow(x, par[0]) * exp(par[1] * x * x) + par[2]);
}


// this is for c flavor
Double_t pdf1(Double_t *xx, Double_t *par) {


  Double_t x = xx[0];

  return par[4] * (pow(x, par[0]) * exp(par[1] * pow(x, par[3])) + par[2]);
}



// make combined functions, assuming that each has parameter size 4
Double_t combined_pdf(Double_t *xx, Double_t *par) {

  Double_t f= 0;

  f += pdf1(xx, &par[0]); 
  f += pdf1(xx, &par[5]); 

  return f;
}
/**************************************************************************/

ClassImp(PtrelSolver)


PtrelSolver::PtrelSolver() { init();}
PtrelSolver::PtrelSolver(double fitmin, double fitmax) {

  init();

  fit_min = fitmin;
  fit_max = fitmax;
}
PtrelSolver::PtrelSolver(double fitmin, double fitmax, int ptbins, int etabins) {

  init();


  fit_min = fitmin;
  fit_max = fitmax;

  pthat_bins = ptbins;
  eta_bins   = etabins;
}
PtrelSolver::~PtrelSolver() {}



void PtrelSolver::init() {

  x_min = 0; 
  x_max = 5;

  fit_min = 0;
  fit_max = 5;

  hist_bins  = 50;
  pthat_bins = 0;
  eta_bins   = 0;

  pt_threshold = 1;
  pt_sumbins = 1;

  eta_threshold = 1;
  eta_sumbins = 1;


  eff_table= new std::vector<std::vector<double> >;
  pdfs_b   = new std::vector<std::vector<double> >;
  pdfs_c   = new std::vector<std::vector<double> >;


  pdfs_b_tag   = new std::vector<std::vector<double> >;
  pdfs_c_tag   = new std::vector<std::vector<double> >;

  c1       = new TCanvas("c1", "", 600, 600);


  label = new TLatex;
  label->SetNDC();
}


void PtrelSolver::setPtAverage(int threshold, int sum) {

  pt_threshold = threshold;
  pt_sumbins = sum;
}

void PtrelSolver::setEtaAverage(int threshold, int sum) {

  eta_threshold = threshold;
  eta_sumbins = sum;
}


/**************************************************************************
 *
 * access pdf index, etc.
 * binning start with 1, 2, ..., n
 *
 **************************************************************************/
int  PtrelSolver::getBin(double xx, double *binning) {

  int nbin = 0;
  while (xx >= binning[nbin]  ) {

    nbin ++;
  }
  return nbin;
}

int PtrelSolver::index(double pt, double eta, double others) {

  int ptnum = getBin(pt, pthatbining);
  int etanum = getBin(eta, etabining);

  return ptnum*PT_BASE + etanum*ETA_BASE;
}
int PtrelSolver::index(int ptnum, int etanum, int others) { return ptnum*PT_BASE + etanum*ETA_BASE; }

// calculate error of N2/N1
double PtrelSolver::effErr(double N1, double N1_err, double N2, double N2_err) {

  return (1.0*N2/N1) * sqrt(pow(1.0*N1_err/N1, 2) + pow(1.0*N2_err/N2, 2));
}



/**************************************************************************
 * 
 *  access pdfs
 *  pdf is expected to be arranged as "pt num * etabins + eta num"
 *         index starts from 0.
 **************************************************************************/
TF1 *PtrelSolver::getAPdf(int ii) {

  if (combined_pdfs.GetLast() < ii) return 0;
  return (TF1 *) combined_pdfs.At(ii);
}
TF1 *PtrelSolver::getAPdf(int pt_bin, int eta_bin) {

  int index = (pt_bin-1) *eta_bins + (eta_bin -1);
  return this->getAPdf(index);
}

TF1 *PtrelSolver::getPdfByIndex(TObjArray *list, int jj, const char *tag) {
  if (!list) return 0;

  // build the pdf index name
  char index[100];
  if (tag && strlen(tag) >0 ) sprintf(index, "%s_%d", tag, jj);
  else sprintf(index, "%d", jj);


  for (int ii = 0; ii <= list->GetLast(); ii ++) {

    TString pdfname(  ((TF1 *)list->At(ii))->GetName()  );
    if (!pdfname.CompareTo(index)) return (TF1 *) list->At(ii);
  }

  return 0;
}


TF1 *PtrelSolver::getPdfByIndex(int jj) {this->getPdfByIndex( &combined_pdfs, jj);}

TF1 *PtrelSolver::getTaggedPdfByIndex(int jj) {

  return  this->getPdfByIndex( &combined_pdfs_tag, jj, "tag");
}

TF1 *PtrelSolver::getPdfByIndex(int jj, const char *tag) {


  TObjArray *list = 0;

  if ( !strcmp(tag, "") )        list = &combined_pdfs;
  if ( !strcmp(tag, "tag") )     list = &combined_pdfs_tag;
  if ( !strcmp(tag, "sys") )     list = &combined_pdfs_sys;
  if ( !strcmp(tag, "sys_tag") ) list = &combined_pdfs_sys_tag;

  return this->getPdfByIndex(list, jj, tag);
}



/**************************************************************************
 *
 * measure the efficciencies with given input data file
 *
 **************************************************************************/
void PtrelSolver::effCal(TH1F *hist,
			 TH1F *hist_tag, 
			 TF1  *pdf,
			 std::vector<double> *eff) {

  if (!eff || !pdf || !hist || !hist_tag) return;
  eff->clear();
  x_min = hist->GetXaxis()->GetXmin();
  x_max = hist->GetXaxis()->GetXmax();


  std::vector<double>  Nb, Nb_err;
  std::vector<double>  Nb_tag, Nb_tag_err;


  Fit(hist, pdf, &Nb, &Nb_err);
  Fit(hist_tag, pdf, &Nb_tag, &Nb_tag_err);


  (*eff).push_back( 1.0*Nb_tag[0]/Nb[0] );
  (*eff).push_back( effErr(Nb[0], Nb_err[0], Nb_tag[0], Nb_tag_err[0]));

  return;
}


void PtrelSolver::effCal(TH1F *hist,
			 TH1F *hist_tag, 
			 TF1  *pdf,
			 TF1 *pdf_tag,
			 std::vector<double> *eff) {
  
  if (!eff || !pdf || !hist || !hist_tag || !pdf_tag) return;
  eff->clear();


  std::vector<double>  Nb, Nb_err;
  std::vector<double>  Nb_tag, Nb_tag_err;


  Fit(hist, pdf, &Nb, &Nb_err);
  Fit(hist_tag, pdf_tag, &Nb_tag, &Nb_tag_err);


  (*eff).push_back( 1.0*Nb_tag[0]/Nb[0] );
  (*eff).push_back( effErr(Nb[0], Nb_err[0], Nb_tag[0], Nb_tag_err[0]));

  return;
}


//
//  make results from root histograms
void PtrelSolver::effCal(TH1F *hist,
			 TH1F *hist_tag, 
			 TF1  *pdf,
			 std::vector<double> *eff,
			 const char *rootfilename){

  if (!eff || !pdf || !hist || !hist_tag || !rootfilename) return;
  TFile *rootfile = new TFile(rootfilename, "UPDATE");


  std::vector<double>  Nb, Nb_err;
  std::vector<double>  Nb_tag, Nb_tag_err;
  char hist_name[100];


  //  a) before applying tagging
  Fit(hist, pdf, &Nb, &Nb_err);
  sprintf(hist_name, "%s_before", pdf->GetName());
  hist->SetName(hist_name);
  rootfile->cd();
  hist->Write();
  c1->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->Draw("PE");
  pdf->Draw("SAME"); pdf->SetLineWidth(2);
  pdf->SetLineColor(kBlack);

  TF1 sig(*pdf);
  sig.SetParameter("Nc", 0);
  sig.SetLineColor(kRed);sig.SetLineWidth(2);
  sig.Draw("SAME");


  TF1 bkg(*pdf);
  bkg.SetParameter("Nb", 0);
  bkg.SetLineColor(kBlue);sig.SetLineWidth(2);
  bkg.Draw("SAME");
  TString epsname(hist_name);
  epsname.Append(".eps");
  c1->SaveAs(epsname.Data());
  /**********************************************************************/


  //b) after applying tagging
  Fit(hist_tag, pdf, &Nb_tag, &Nb_tag_err);
  sprintf(hist_name, "%s_tag", pdf->GetName());
  hist->SetName(hist_name);
  rootfile->cd();
  hist->Write();
  c1->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->Draw("PE");
  pdf->Draw("SAME"); pdf->SetLineWidth(2);
  pdf->SetLineColor(kBlack);

  TF1 sig_tag(*pdf);
  sig_tag.SetParameter("Nc", 0);
  sig_tag.SetLineColor(kRed);sig_tag.SetLineWidth(2);
  sig_tag.Draw("SAME");


  TF1 bkg_tag(*pdf);
  bkg_tag.SetParameter("Nb", 0);
  bkg_tag.SetLineColor(kBlue);sig_tag.SetLineWidth(2);
  bkg_tag.Draw("SAME");
  TString epsname_tag(hist_name);
  epsname_tag.Append(".eps");
  c1->SaveAs(epsname_tag.Data());
  /**********************************************************************/



  (*eff).push_back( 1.0*Nb_tag[0]/Nb[0] );
  (*eff).push_back( effErr(Nb[0], Nb_err[0], Nb_tag[0], Nb_tag_err[0]));


  rootfile->cd();
  rootfile->Write();
  rootfile->Close();
}





TH1F *PtrelSolver::getMCeff(TFile *file, const char *hist) {
  if (!file || !hist) return 0;

  TString raw("b_"); raw.Append(hist);
  TString tag("b_"); tag.Append(hist);
  tag.Insert(3, "cmb");

  TString eff(hist); eff.Append("_eff_mc");

  TH2F *raw_hist = (TH2F *)file->Get(raw.Data());
  TH2F *tag_hist = (TH2F *)file->Get(tag.Data());


  TH1F *h1 = (TH1F *)tag_hist->ProjectionX("total");
  TH1F *h2 = (TH1F *)raw_hist->ProjectionX("tagged");

  h1->Sumw2(); h2->Sumw2();

  h1->Divide(h2);
  h1->SetName(eff.Data());
  return h1;
}

void PtrelSolver::measure(const char *inputfilename, const char *tag, const char *pthist, const char *etahist, bool sys) {

  TString outfile(tag);
  outfile.Append(".root");

  measure(inputfilename, outfile.Data(), tag, pthist, etahist, sys);
}


void PtrelSolver::measure(const char *inputfilename, const char *outfilename, const char *tag, const char *pthist, const char *etahist, bool sys) {
  



  TH1F *hist=0,  *hist_tag=0, *eff_sys=0;
  TH2F *hist2=0, *hist2_tag = 0;
  TH1F *eff =0;
  TH1F *mc_eff;
  TH1F *ratio;
  TF1  *thePdf= 0;
  TF1  *thePdf_tag= 0;
  TF1  *thePdf_sys= 0;
  TF1  *thePdf_sys_tag= 0;
  int   pdf_index;
  int   nbins;
  std::vector<double> result, result_sys;

  TString eff_name_sys;
  TString eff_name; 
  TString tag_name; // to access the tagged histogram

  //  input and output files.
  inputfile = new TFile(inputfilename, "READ");
  outfile = new TFile(outfilename, "RECREATE");
  char data_name[100];


  c1->cd();

  // 1. pt dependence
  // duplicate the binning for efficiency table
  hist2 = (TH2F *)inputfile->Get(pthist);
  tag_name.Append(pthist); tag_name.Insert(1, "cmb");
  hist2_tag = (TH2F *)inputfile->Get(tag_name.Data());
  std::cout << "information: " << "accessing histogram ... " << pthist << std::endl;
  std::cout << "information: " << "accessing histogram ... " << tag_name.Data() << std::endl;
  if (!hist2 || !hist2_tag) {

    std::cout << "information: error in access histograms " << std::endl;
    return;
  }

  eff_name.Append("eff_pt_");
  eff_name.Append(tag);
  eff = (TH1F *)hist2->ProjectionX(eff_name.Data());
  eff_name_sys.Append(eff_name.Data());
  eff_name_sys.Append("_sys");
  eff_sys = (TH1F *)hist2->ProjectionX(eff_name_sys.Data());


  nbins = eff->GetNbinsX();
  int total_num_plots = nbins * 2;
  total_num_plots = (int) (sqrt(total_num_plots) + 1);
  TCanvas *c2 = new TCanvas("c2", "", 900, 900);
  c2->Divide(total_num_plots, total_num_plots);

  for (int ii = 0; ii <= nbins; ii++) {


    c1->cd();
    //    std::cout << "processing bin " << ii << std::endl; 
    thePdf = getPdfByIndex( index(ii, 0) );
    thePdf_tag = getTaggedPdfByIndex( index(ii, 0) );
    if (!thePdf || !thePdf_tag) {
      std::cout << "can't locate the correct PDF for this data set ..." 
		<< std::endl;
      break;
    } else {


      //      std::cout << "information: proper pdfs are accessed " << std::endl;
    }


    // before tag
    sprintf(data_name, "ptbin_%d_%s", ii, tag);
    if (ii !=0) {

      hist = (TH1F *)hist2->ProjectionY(data_name, ii, ii);
    } else 
      hist = (TH1F *)hist2->ProjectionY(data_name, 1, nbins);


    // after tag
    sprintf(data_name, "ptbin_%d_%s_tag", ii, tag);
    if (ii !=0) {
     
      hist_tag = (TH1F *)hist2_tag->ProjectionY(data_name, ii, ii);
    }    else 
      hist_tag = (TH1F *)hist2_tag->ProjectionY(data_name, 1, nbins);



    //    std::cout << " reach here" << std::endl;
    effCal(hist, hist_tag, thePdf, thePdf_tag,  &result);

    if (ii == 3) {
      formatHist1(hist, "p_{T}^{rel} (GeV)", "Events");
      formatHist1(hist_tag, "p_{T}^{rel} (GeV)", "Events");
      makePlot("pt40_50.eps", hist, thePdf);
      makePlot("pt40_50_tag.eps", hist_tag, thePdf_tag);

    }

    c2->cd(2*ii+1); hist->Draw("HIST"); thePdf->Draw("SAME");
    c2->cd(2*ii+2); hist_tag->Draw("HIST"); thePdf_tag->Draw("SAME");
    c1->cd();


    // estimate systematic from some fitting
    if (sys) {

      thePdf_sys     = getPdfByIndex( index(ii, 0), "sys" );
      thePdf_sys_tag = getPdfByIndex( index(ii, 0), "sys_tag" );
      if (thePdf_sys && thePdf_sys_tag) {
	effCal(hist, hist_tag, thePdf_sys, thePdf_sys_tag,  &result_sys);
	result_sys[1] = sqrt(pow(result[1], 2) + pow(result_sys[0] - result[0] , 2));


	c2->cd(2*ii+1); hist->Draw("HIST"); thePdf_sys->Draw("SAME");
	c2->cd(2*ii+2); hist_tag->Draw("HIST"); thePdf_sys_tag->Draw("SAME");
	c1->cd();


	eff_sys->SetBinContent(ii, result[0]);
	eff_sys->SetBinError(ii, result_sys[1]); // only the error is different.
      } else {
	
	std::cout << "information: no systematics applied !" << std::endl;
      }
    }


    if (ii == 0) {
      std::cout <<  "tagger - " << tag << " average efficiency " 
		<< result[0] <<"+/-" 
		<< result[1] <<"(stat.)+/-";

      if (sys)  std::cout <<result_sys[1] << "(total)" << std::endl;
      continue;
    }



    eff->SetBinContent(ii, result[0]);
    eff->SetBinError(ii, result[1]);


    outfile->cd();
    hist->Write();
    hist_tag->Write();
  }



  TString mc_eff_name(eff_name); 
  inputfile->cd();

  mc_eff = getMCeff(inputfile, pthist);

  hist = Divide(eff, mc_eff);

  outfile->cd();
  eff->Write();
  eff_sys->Write();
  mc_eff->SetLineColor(kRed);
  mc_eff->Write();
  hist->SetLineColor(kBlue);
  hist->Write();



  // make histograms.
  //
  TString title;
  TLegend *leg = 0;
  const char *epsname= 0;

  eff->SetMarkerStyle(20);
  eff->SetMinimum(0);
  eff->SetMaximum(1.35);
  formatHist1(eff, "p_{T} (GeV)", "Eff.");
  eff_sys->SetLineColor(kBlack);
  eff_sys->SetMarkerColor(kBlack);

  epsname = saveAsEPS(c1, eff, "PE1");

  mc_eff->SetMarkerSize(0);
  leg = new TLegend(0.25, 0.77, 0.6, 0.92);
  title.Append(tag);
  leg->SetHeader(title.Data());
  leg->AddEntry(eff, "P_{T}^{rel} Fit", "pl");
  leg->AddEntry(mc_eff, "MC Truth", "l");
  leg->SetFillColor(0);
  leg->SetFillColor(0);
  leg->Draw();
  c1->cd();
  mc_eff->Draw("SAMEHISTE1");
  if (sys) eff_sys->Draw("SAMEPE1");
  c1->SaveAs(epsname);

  TString big_plot(tag);
  big_plot.Append("_pt.eps");
  c2->SaveAs(big_plot.Data());

  // 2. eta dependence
  // duplicate the binning for efficiency table
  tag_name.Resize(0);
  tag_name.Append(etahist); tag_name.Insert(1, "cmb");
  std::cout << "information: " << "accessing histogram ... " << etahist << std::endl;
  std::cout << "information: " << "accessing histogram ... " << tag_name.Data() << std::endl;

  hist2 = (TH2F *)inputfile->Get(etahist);
  hist2_tag = (TH2F *)inputfile->Get(tag_name.Data());
  eff_name.Resize(0);
  eff_name.Append("eff_eta_");
  eff_name.Append(tag);
  eff = (TH1F *)hist2->ProjectionX(eff_name.Data());

  eff_name_sys.Resize(0);
  eff_name_sys.Append(eff_name);
  eff_name_sys.Append("_sys");
  eff_sys = (TH1F *)hist2->ProjectionX(eff_name_sys.Data());


  nbins = eff->GetNbinsX();
  big_plot.Resize(0);
  big_plot.Append(tag);
  big_plot.Append("_eta.eps");
  total_num_plots = nbins * 2;
  total_num_plots = (int) (sqrt(total_num_plots) + 1);
  c2 = new TCanvas("c2", "", 900, 900);
  c2->Divide(total_num_plots, total_num_plots);

  //  std::cout << "total operating points " << nbins << std::endl;
  for (int ii = 1; ii <= nbins; ii++) {

    c1->cd();
    thePdf = getPdfByIndex( index(0, ii) );
    thePdf_tag = getTaggedPdfByIndex( index(0, ii) );
    if (!thePdf || !thePdf_tag) {
      std::cout << "can't locate the correct PDF for this data set ..." 
		<< std::endl;
      break;
    }
    //    std::cout << "fitting data with PDF " << index(0, ii) << std::endl;

    // before tag
    sprintf(data_name, "etabin_%d_%s", ii, tag);
    hist = (TH1F *)hist2->ProjectionY(data_name, ii, ii);


    // after tag
    sprintf(data_name, "etabin_%d_%s_tag", ii, tag);
    hist_tag =(TH1F *) hist2_tag->ProjectionY(data_name, ii, ii);
    

    effCal(hist, hist_tag, thePdf, thePdf_tag,  &result);
    c2->cd(2*(ii-1)+1); hist->Draw("HIST"); thePdf->Draw("SAME");
    c2->cd(2*(ii-1)+2); hist_tag->Draw("HIST"); thePdf_tag->Draw("SAME");
    c1->cd();


    if (sys) {
      thePdf_sys     = getPdfByIndex( index(0, ii), "sys" );
      thePdf_sys_tag = getPdfByIndex( index(0, ii), "sys_tag" );
      if (thePdf_sys && thePdf_sys_tag) {
	effCal(hist, hist_tag, thePdf_sys, thePdf_sys_tag,  &result_sys);
	result_sys[1] = sqrt(pow(result[1], 2) + pow(result_sys[0] - result[0] , 2));


	eff_sys->SetBinContent(ii, result[0]);
	eff_sys->SetBinError(ii, result_sys[1]);


	c2->cd(2*ii+1); hist->Draw("HIST"); thePdf_sys->Draw("SAME");
	c2->cd(2*ii+2); hist_tag->Draw("HIST"); thePdf_sys_tag->Draw("SAME");
	c1->cd();


      } else {

	std::cout << "information: no systematics applied !" << std::endl;
      }
    }

    eff->SetBinContent(ii, result[0]);
    eff->SetBinError(ii, result[1]);



    outfile->cd();
    hist->Write();
    hist_tag->Write();
  }

  c2->SaveAs(big_plot.Data());
  delete c2;
  c1->cd();



  mc_eff = getMCeff(inputfile, etahist);
  hist = Divide(eff, mc_eff);

  mc_eff->SetLineColor(kRed);
  hist->SetLineColor(kBlue);

  outfile->cd();
  eff->Write();
  eff_sys->Write();
  mc_eff->Write();
  hist->Write();


  // make plots
  eff->SetMarkerStyle(20);
  eff->SetMinimum(0);
  eff->SetMaximum(1.35);
  eff_sys->SetLineColor(kBlack);
  eff_sys->SetMarkerColor(kBlack);

  mc_eff->SetMarkerSize(0);

  formatHist1(eff, "|#eta|", "Eff.");
  epsname = saveAsEPS(c1, eff, "PE1");
  leg = new TLegend(0.55, 0.75, 0.9, 0.90);
  title.Resize(0);
  title.Append(tag);
  leg->SetHeader(title.Data());
  leg->AddEntry(eff, "P_{T}^{rel} Fit ", "pl");
  leg->AddEntry(mc_eff, "MC Truth", "l");
  leg->SetFillColor(0);
  leg->Draw();
  c1->cd();
  mc_eff->Draw("SAMEHISTE1");
  if (sys) eff_sys->Draw("SAMEPE1");
  // hist->Draw("SAMEHISTE");
  c1->SaveAs(epsname);



  // save output
  outfile->Write();
  outfile->Close();
}


// Uses the counting method to measure efficiencies
void PtrelSolver::counting(
  const char *tag,
  const char *outfilename,
  const char *inputfilename,
  const char *mistagfilename
)
{
  TString baseName, name; 

  // output file with results
  outfile = new TFile(outfilename, "RECREATE");

  // input file with data.
  inputfile = new TFile(inputfilename, "READ");

  // input file with mistag data.
  if (mistagfilename)
    mistagfile = new TFile(mistagfilename, "READ");
 
  char data_name[100];

  // pt dependence

  // pT dependency for muon-in-jet sample
  TH2F * hist2d = (TH2F*) inputfile->Get("npT");

  // pT dependency for muon-in-jet-away-jet-tagged sample
  TH2F * hist2d_tag = (TH2F*) inputfile->Get("ppT");

  // pT dependency for muon-in-jet-tagged-away-jet-tagged sample
  TH2F * hist2d_two_tag = (TH2F*) inputfile->Get("pcmbpT");
    
  // caculation of mistag (extracted from Daniel's histograms)
  TH1F * mistag = (TH1F*) mistagfile->Get("Hj132");
  mistag->Add((TH1*) mistagfile->Get("Hj432"));
  TH1F * denominator = (TH1F*) mistagfile->Get("Hj102");
  denominator->Add((TH1*) mistagfile->Get("Hj402"));  
  mistag->Divide(mistag, denominator, 1., 1., "B");
   
  // histogram names
  baseName = "eff_counting_pt_";
  baseName.Append(tag);
  
  TH1F * eff = (TH1F*) hist2d_tag->ProjectionX(baseName.Data());
  Int_t nbins = eff->GetNbinsX();

  name = baseName;
  name.Append("_plus");
  TH1F * effPlus = (TH1F*) hist2d_tag->ProjectionX(name.Data());

  name = baseName;
  name.Append("_minus");
  TH1F * effMinus = (TH1F*) hist2d_tag->ProjectionX(name.Data());

  name = baseName;
  name.Append("_ncl_mc");
  TH2F * hist2d_cl = (TH2F*) inputfile->Get("cl_npT");
  TH1F * ncl_mc = (TH1F*) hist2d_cl->ProjectionX(name.Data(), -1, -1, "e");

  name = baseName;
  name.Append("_nb_tag_mc");
  TH2F * hist2d_tag_b = (TH2F*) inputfile->Get("b_ppT");
  TH1F * nb_tag_mc = (TH1F*) hist2d_tag_b->ProjectionX(name.Data(), -1, -1, "e");

  name = baseName;
  name.Append("_mtg_mc");
  TH2F * hist2d_tag_cl = (TH2F*) inputfile->Get("cl_ppT");
  TH1F * mtg_mc = (TH1F*) hist2d_tag_cl->ProjectionX(name.Data(), -1, -1, "e");
  mtg_mc->Divide(mtg_mc, ncl_mc, 1., 1., "e");

  name = baseName;
  name.Append("_nb_two_tag_mc");
  TH2F * hist2d_two_tag_b = (TH2F*) inputfile->Get("b_pcmbpT");
  TH1F * nb_two_tag_mc = (TH1F*) hist2d_two_tag_b->ProjectionX(name.Data(), -1, -1, "e");

  TH1F * hist;
  TH1F * hist_tag;
  TH1F * hist_two_tag;

  TF1 * thePdf1;
  TF1 * thePdf2;
   
  std::vector<double>  n, n_err;
  std::vector<double>  nt, nt_err;
  
  Double_t Nst, Mtg, MtgErr;
  Double_t Eff, EffMtgPlus, EffMtgMinus; 
  Double_t EffStatE2, EffSystE2Plus, EffSystE2Minus;
        
  for (int ii = 1; ii <= nbins; ii++)
  {
    std::cout << "processing bin " << ii << std::endl; 

    // get the right template
    thePdf1 = getPdfByIndex( index(ii, 0) );
    thePdf2 = getTaggedPdfByIndex( index(ii, 0) );

    // if (!thePdf1 || !thePdf2 || !thePdfSys1 || !thePdfSys2)
    if (!thePdf1 || !thePdf2)
    {
      std::cout << "can't locate the correct PDF for this data set ..." << std::endl;
      break;
    }

    // before tag
    sprintf(data_name, "ptbin_%d_%s", ii, tag);
    hist = (TH1F*) hist2d->ProjectionY(data_name, ii, ii);

    // after tagging away jet
    sprintf(data_name, "ptbin_%d_%s_tag", ii, tag);
    hist_tag = (TH1F*) hist2d_tag->ProjectionY(data_name, ii, ii);
    
    // after tagging mu-jet and away jet
    sprintf(data_name, "ptbin_%d_%s_two_tag", ii, tag);
    hist_two_tag = (TH1F*) hist2d_two_tag->ProjectionY(data_name, ii, ii);

    // get the b fraction in the muon-in-jet-away-jet-tagged
    Fit(hist, thePdf1, &n, &n_err);
          
    // get the b fraction in the muon-in-jet-tagged-away-jet-tagged
    Fit(hist_two_tag, thePdf2, &nt, &nt_err);
        
    Nst = (Double_t) hist_tag->GetEntries();

    // red the measure lights mistag rate or get it from MC
    if (mistagfilename)
    {
      Mtg = mistag->GetBinContent(ii); 
      MtgErr = mistag->GetBinError(ii); 
    }
    else
    {
      Mtg = mtg_mc->GetBinContent(ii);
      MtgErr = mtg_mc->GetBinError(ii);
    }
    
    Eff         = nt[0] / ( Nst - n[1] * Mtg );
    EffMtgPlus  = nt[0] / ( Nst - n[1] * Mtg * 1.05 );
    EffMtgMinus = nt[0] / ( Nst - n[1] * Mtg * 0.95 );

    EffStatE2 = pow(nt_err[0] / ( Nst - n[1] *  Mtg ), 2)
              + pow(nt[0] * n[1] * MtgErr / pow(Nst - Mtg * n[1], 2), 2)
              + pow(nt[0] * Mtg  * n_err[1] / pow(Nst - Mtg * n[1], 2), 2);
              
    EffSystE2Plus  = pow( Eff - EffMtgPlus, 2);
    EffSystE2Minus = pow( Eff - EffMtgMinus, 2);

    std::cout << "muon-in-jet (# measure c, # truth c)                        : (" << n[1] << ", "  << ncl_mc->GetBinContent(ii) << ")" << std::endl;
    std::cout << "muon-in-jet-awat-jet-tagged (# measure b, # truth b)        : (" << (Nst - n[1] * Mtg) << ", " << nb_tag_mc->GetBinContent(ii) << ")" << std::endl;
    std::cout << "muon-in-jet-tagged-awat-jet-tagged (# measure b, # truth b) : (" << nt[0] << ", " << nb_two_tag_mc->GetBinContent(ii) << ")" << std::endl;
    
    std::cout << "Eff             : " << Eff << std::endl;
    std::cout << "EffStatErr      : " << sqrt(EffStatE2) << std::endl;
    std::cout << "EffSystErrPlus  : " << sqrt(EffSystE2Plus) << std::endl;
    std::cout << "EffSystErrMinus : " << sqrt(EffSystE2Minus) << std::endl;    
    std::cout << "Mistag          : " << Mtg << " pm " << MtgErr << std::endl;
    
    eff->SetBinContent(ii, Eff);
    eff->SetBinError(ii, sqrt(EffStatE2));

    effPlus->SetBinContent(ii, Eff+sqrt(EffStatE2+EffSystE2Plus));
    effPlus->SetBinError(ii, 0.0);
    effMinus->SetBinContent(ii, Eff-sqrt(EffStatE2+EffSystE2Minus));
    effMinus->SetBinError(ii, 0.0);

    outfile->cd();
    hist->Write();
    hist_tag->Write();
    hist_two_tag->Write();
  }

  inputfile->cd();
  
  hist = (TH1F *)inputfile->Get("eff_TaggedBothJets_b");  
  TH1F * mc_eff = new TH1F(*hist);

  name = baseName;
  name.Append("_mc");
  mc_eff->SetName(name.Data());
  mc_eff->SetLineColor(kRed);
 
  outfile->cd();
  eff->Write();
  mc_eff->Write();
  effPlus->Write();
  effMinus->Write();

  // make histograms.
  //
  TString title;
  TLegend *leg = 0;
  const char *epsname= 0;

  formatHist1(effPlus, "p_{T} (GeV)", "Eff.");
  effPlus->GetYaxis()->SetRangeUser(0.0,1.4);
  effPlus->SetFillColor(kBlack);
  effPlus->SetFillStyle(3005);
  effPlus->SetLineColor(kWhite);
  
  epsname = saveAsEPS(c1, effPlus, "BAR");
  effMinus->SetLineColor(kWhite);
  effMinus->Draw("BAR SAME");
  eff->Draw("SAME");
  mc_eff->Draw("SAME E");
    
  mc_eff->SetMarkerSize(0);  
  leg = new TLegend(0.25, 0.8, 0.6, 0.9);
  leg->AddEntry(eff, "Counting", "pl");
  leg->AddEntry(mc_eff, "MC truth", "l");
  leg->SetFillColor(0);
  leg->Draw();

  gPad->RedrawAxis();

  c1->cd();
  c1->SaveAs(epsname);

  // Calculation of the over all efficiency ---------------------------
 
  // MC Mistag
  name = baseName;
  name.Append("_total_mtg_mc");
  mtg_mc = (TH1F*) hist2d_tag_cl->ProjectionX(name.Data(),-1,-1,"e");
 
  // get ptrel of all the bins
  thePdf1 = getPdfByIndex( index(0, 0) );  
  thePdf2 = getTaggedPdfByIndex( index(0, 0) );
  
  // before tag
  sprintf(data_name, "ptbin_%d_%s", 0,tag);
  hist = (TH1F*) hist2d->ProjectionY(data_name);

  // after tagging away jet
  sprintf(data_name, "ptbin_%d_%s_tag", 0, tag);
  hist_tag = (TH1F*) hist2d_tag->ProjectionY(data_name);
    
  // after tagging mu-jet and away jet
  sprintf(data_name, "ptbin_%d_%s_two_tag", 0, tag);
  hist_two_tag = (TH1F*) hist2d_two_tag->ProjectionY(data_name);

  // get the b fraction in the muon-in-jet-away-jet-tagged
  Fit(hist, thePdf1, &n, &n_err);
    
  // get the b fraction in the muon-in-jet-tagged-away-jet-tagged
  Fit(hist_two_tag, thePdf2, &nt, &nt_err);

  // red the measure lights mistag rate or get it from MC
  if (mistagfilename)
  { 
  	// FIXME
    Mtg = 0.0; 
    MtgErr = 0.0; 
  }
  else
  {
    Mtg = mtg_mc->GetEntries()/ncl_mc->GetEntries();  	
    MtgErr = 0;
  }
   
  Nst = (Double_t) hist_tag->GetEntries();

  Eff = nt[0] / ( Nst - n[1] * Mtg );
  EffMtgPlus  = nt[0] / ( Nst - n[1] * Mtg * 1.05 );
  EffMtgMinus = nt[0] / ( Nst - n[1] * Mtg * 0.95 );

  EffStatE2 = pow(nt_err[0] / ( Nst - n[1] *  Mtg ), 2)
            + pow(nt[0] * n[1] * MtgErr   / pow(Nst - Mtg * n[1], 2), 2)
            + pow(nt[0] * Mtg  * n_err[1] / pow(Nst - Mtg * n[1], 2), 2);
            
  EffSystE2Plus  = pow( Eff - EffMtgPlus, 2);
  EffSystE2Minus = pow( Eff - EffMtgMinus, 2);

  std::cout << std::endl << "===============================================================================" << std::endl;

  std::cout << "muon-in-jet (# measure c, # truth c)                        : (" << n[1] << ", " << ncl_mc->GetEntries() << ")" << std::endl;
  std::cout << "muon-in-jet-awat-jet-tagged (# measure b, # truth b)        : (" << (Nst - n[1] * Mtg) << ", " << nb_tag_mc->GetEntries() << ")" << std::endl;
  std::cout << "muon-in-jet-tagged-awat-jet-tagged (# measure b, # truth b) : (" << nt[0] << ", " << nb_two_tag_mc->GetEntries() << ")" << std::endl;
    
  std::cout << "Total Eff             : " << Eff << std::endl;
  std::cout << "Total EffStatErr      : " << sqrt(EffStatE2) << std::endl;
  std::cout << "Total EffSystErrPlus  : " << sqrt(EffSystE2Plus) << std::endl;
  std::cout << "Total EffSystErrMinus : " << sqrt(EffSystE2Minus) << std::endl;   
  std::cout << "Total EffErrPlus      : " << sqrt(EffStatE2+EffSystE2Plus) << std::endl;   
  std::cout << "Total EffErrMinus     : " << sqrt(EffStatE2+EffSystE2Minus) << std::endl;      
  std::cout << "Total Mistag          : " << Mtg << " mc: " << MtgErr << std::endl;

  std::cout << "===============================================================================" << std::endl << std::endl;
    
  // --------------------------------------------

  // Eta dependence

  // eta dependency for muon-in-jet sample
  hist2d = (TH2F*) inputfile->Get("nEta");

  // eta dependency for muon-in-jet-away-jet-tagged sample
  hist2d_tag = (TH2F*) inputfile->Get("pEta");

  // eta dependency for muon-in-jet-tagged-away-jet-tagged sample
  hist2d_two_tag = (TH2F*) inputfile->Get("pcmbEta");
    
  // caculation of mistag (extracted from Daniel's histograms)
  mistag = (TH1F*) mistagfile->Get("Hj133");
  mistag->Add((TH1*) mistagfile->Get("Hj433"));
  denominator = (TH1F*) mistagfile->Get("Hj103");
  denominator->Add((TH1*) mistagfile->Get("Hj403"));  
  mistag->Divide(mistag, denominator, 1., 1., "B");
   
  // histogram names
  baseName = "eff_counting_eta_";
  baseName.Append(tag);
  
  eff = (TH1F*) hist2d_tag->ProjectionX(baseName.Data());
  nbins = eff->GetNbinsX();
  
  name = baseName;
  name.Append("_plus");
  effPlus = (TH1F*) hist2d_tag->ProjectionX(name.Data());

  name = baseName;
  name.Append("_minus");
  effMinus = (TH1F*) hist2d_tag->ProjectionX(name.Data());

  name = baseName;
  name.Append("_ncl_mc");
  hist2d_cl = (TH2F*) inputfile->Get("cl_nEta");
  ncl_mc = (TH1F*) hist2d_cl->ProjectionX(name.Data(), -1, -1, "e");

  name = baseName;
  name.Append("_nb_tag_mc");
  hist2d_tag_b = (TH2F*) inputfile->Get("b_pEta");
  nb_tag_mc = (TH1F*) hist2d_tag_b->ProjectionX(name.Data(), -1, -1, "e");

  name = baseName;
  name.Append("_mtg_mc");
  hist2d_tag_cl = (TH2F*) inputfile->Get("cl_pEta");
  mtg_mc = (TH1F*) hist2d_tag_cl->ProjectionX(name.Data(), -1, -1, "e");
  mtg_mc->Divide(mtg_mc, ncl_mc, 1., 1., "e");

  name = baseName;
  name.Append("_nb_two_tag_mc");
  hist2d_two_tag_b = (TH2F*) inputfile->Get("b_pcmbEta");
  nb_two_tag_mc = (TH1F*) hist2d_two_tag_b->ProjectionX(name.Data(), -1, -1,"e");
  
  for (int ii = 1; ii <= nbins; ii++)
  {
    std::cout << "processing bin " << ii << std::endl; 

    // get the right template
    thePdf1 = getPdfByIndex( index(0, ii) );
    thePdf2 = getTaggedPdfByIndex( index(0, ii) );

    if (!thePdf1 || !thePdf2)
    {
      std::cout << "can't locate the correct PDF for this data set ..." << std::endl;
      break;
    }

    // before tag
    sprintf(data_name, "etabin_%d_%s", ii, tag);
    hist = (TH1F*) hist2d->ProjectionY(data_name, ii, ii);

    // after tagging away jet
    sprintf(data_name, "etabin_%d_%s_tag", ii, tag);
    hist_tag = (TH1F*) hist2d_tag->ProjectionY(data_name, ii, ii);
    
    // after tagging mu-jet and away jet
    sprintf(data_name, "etabin_%d_%s_two_tag", ii, tag);
    hist_two_tag = (TH1F*) hist2d_two_tag->ProjectionY(data_name, ii, ii);

    // get the b fraction in the muon-in-jet-away-jet-tagged
    Fit(hist, thePdf1, &n, &n_err);
    
    // get the b fraction in the muon-in-jet-tagged-away-jet-tagged
    Fit(hist_two_tag, thePdf2, &nt, &nt_err);
   
    Nst = (Double_t) hist_tag->GetEntries();

    // red the measure lights mistag rate or get it from MC
    if (mistagfilename)
    {
      Mtg = mistag->GetBinContent(ii); 
      MtgErr = mistag->GetBinError(ii); 
    }
    else
    {
      Mtg = mtg_mc->GetBinContent(ii);
      MtgErr = mtg_mc->GetBinError(ii);
    }

    Eff         = nt[0] / ( Nst - n[1] * Mtg );
    EffMtgPlus  = nt[0] / ( Nst - n[1] * Mtg * 1.05 );
    EffMtgMinus = nt[0] / ( Nst - n[1] * Mtg * 0.95 );

    EffStatE2 = pow(nt_err[0] / ( Nst - n[1] *  Mtg ), 2)
              + pow(nt[0] * n[1] * MtgErr   / pow(Nst - Mtg * n[1], 2), 2)
              + pow(nt[0] * Mtg  * n_err[1] / pow(Nst - Mtg * n[1], 2), 2);
              
    EffSystE2Plus  = pow( Eff - EffMtgPlus, 2);
    EffSystE2Minus = pow( Eff - EffMtgMinus, 2);

    std::cout << "muon-in-jet (# measure c, # truth c)                        : (" << n[1] << ", " << ncl_mc->GetBinContent(ii) << ")" << std::endl;
    std::cout << "muon-in-jet-awat-jet-tagged (# measure b, # truth b)        : (" << (Nst - n[1] * Mtg) << ", " << nb_tag_mc->GetBinContent(ii) << ")" << std::endl;
    std::cout << "muon-in-jet-tagged-awat-jet-tagged (# measure b, # truth b) : (" << nt[0] << "," << nb_two_tag_mc->GetBinContent(ii) << ")" << std::endl;
        
    std::cout << "Eff             : " << Eff << std::endl;
    std::cout << "EffStatErr      : " << sqrt(EffStatE2) << std::endl;
    std::cout << "EffSystErrPlus  : " << sqrt(EffSystE2Plus) << std::endl;
    std::cout << "EffSystErrMinus : " << sqrt(EffSystE2Minus) << std::endl;    
    std::cout << "Mistag : " << Mtg << " mc: " << MtgErr << std::endl;
    
    eff->SetBinContent(ii, Eff);
    eff->SetBinError(ii, sqrt(EffStatE2));

    effPlus->SetBinContent(ii, Eff+sqrt(EffStatE2+EffSystE2Plus));
    effPlus->SetBinError(ii, 0.0);
    effMinus->SetBinContent(ii, Eff-sqrt(EffStatE2+EffSystE2Minus));
    effMinus->SetBinError(ii, 0.0);

    outfile->cd();
    hist->Write();
    hist_tag->Write();
    hist_two_tag->Write();    
  }

  inputfile->cd();
  hist = (TH1F *)inputfile->Get("eff_TaggedBothJets_eta_b");
  mc_eff = new TH1F(*hist);

  name = baseName;
  name.Append("_mc");
  mc_eff->SetName(name.Data());
  mc_eff->SetLineColor(kRed);

  outfile->cd();
  eff->Write();
  mc_eff->Write();
  effPlus->Write();
  effMinus->Write();
  
  formatHist1(effPlus, "|#eta|", "Eff.");
  effPlus->GetYaxis()->SetRangeUser(0.0,1.4);
  effPlus->SetFillColor(kBlack);
  effPlus->SetFillStyle(3005);
  effPlus->SetLineColor(kWhite);
  
  epsname = saveAsEPS(c1, effPlus, "BAR");
  effMinus->SetLineColor(kWhite);
  effMinus->Draw("BAR SAME");
  eff->Draw("SAME");
  mc_eff->Draw("SAME E");
  
  mc_eff->SetMarkerSize(0);  
  leg = new TLegend(0.25, 0.8, 0.6, 0.9);
  leg->AddEntry(eff, "Counting", "pl");
  leg->AddEntry(mc_eff, "MC truth", "l");
  leg->SetFillColor(0);
  leg->Draw();

  gPad->RedrawAxis();

  c1->cd();
  c1->SaveAs(epsname);
    
  // save output
  outfile->Write();
  outfile->Close();
}


void PtrelSolver::makePlot(const char *epsname, TH1 *hist, TF1 *pdf) {

  if (!epsname || !hist || !pdf) return;

  c1->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->Draw("PE");
  pdf->Draw("SAME"); pdf->SetLineWidth(2);
  pdf->SetLineColor(kBlack);

  TF1 sig_tag(*pdf);
  sig_tag.SetParameter("Nc", 0);
  sig_tag.SetLineColor(kRed);sig_tag.SetLineWidth(2);
  sig_tag.Draw("SAME");


  TF1 bkg_tag(*pdf);
  bkg_tag.SetParameter("Nb", 0);
  bkg_tag.SetLineColor(kBlue);sig_tag.SetLineWidth(2);
  bkg_tag.Draw("SAME");

  c1->SaveAs(epsname);
}



void PtrelSolver::estimate(TCut tagger, 
		      std::vector<std::vector<double> > *eff, 
		      std::vector<std::vector<double> > *b, 
		      std::vector<std::vector<double> > *c,
		      const char *rootfilename,
		      const char *filename) {

  if (!eff || !b || !c || !rootfilename) return;
  TFile *rootfile = new TFile(rootfilename, "RECREATE");

  file.open(filename, ios::trunc | ios::out); // dat file to save fit resutls.
  eff->clear();


  std::vector<double>  Nb, Nb_err;
  std::vector<double>  Nb_tag, Nb_tag_err;


  char hist_name[100];
  for (unsigned int ii = 0; ii < (*b).size(); ii ++) {
    std::vector<double> xx;
    (*eff).push_back(xx); // add one entry

    (*eff)[ii].push_back( (*b)[ii][0] ); file<< setw(5) << (*b)[ii][0] << " ";


    //    (*eff)[ii].push_back( (*b)[ii][1] ); file<< setw(5) << (*b)[ii][1] << " ";
    // (*eff)[ii].push_back( (*b)[ii][2] ); file<< setw(5) << (*b)[ii][2] << " ";
    // (*eff)[ii].push_back( (*b)[ii][3] ); file<< setw(5) << (*b)[ii][3] << " ";
    // build the binning


    func = buildAPdf(ii, b, c);


    // data efficiency
    hist = getData( (*b)[ii][0], (*b)[ii][1], (*b)[ii][2], (*b)[ii][3], muon );
    Fit(hist, func, &Nb, &Nb_err);
    file<< setw(5) << hist->GetEntries() << " ";
    file<< setw(5) << Nb[0] << " ";
    file<< setw(5) << Nb_err[0] << " ";
    file<< setw(5) << Nb[1] << " ";
    file<< setw(5) << Nb_err[1] << " ";

    /**********************************************************************/
    //  save the results, make some plots
    sprintf(hist_name, "data_pt%0.3f_eta%0.3f", (*b)[ii][0],  (*b)[ii][2]);
    hist->SetName(hist_name);
    rootfile->cd();
    hist->Write();
    c1->cd();
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
    hist->Draw("PE");
    func->Draw("SAME"); func->SetLineWidth(2);
    func->SetLineColor(kBlack);

    TF1 sig(*func);
    sig.SetParameter("Nc", 0);
    sig.SetLineColor(kRed);sig.SetLineWidth(2);
    sig.Draw("SAME");


    TF1 bkg(*func);
    bkg.SetParameter("Nb", 0);
    bkg.SetLineColor(kBlue);sig.SetLineWidth(2);
    bkg.Draw("SAME");
    TString epsname(hist_name);
    epsname.Append(".eps");
    c1->SaveAs(epsname.Data());
    /**********************************************************************/



    hist = getData( (*b)[ii][0], (*b)[ii][1], (*b)[ii][2], (*b)[ii][3], muon && tagger );
    Fit(hist, func, &Nb_tag, &Nb_tag_err);
    file<< setw(5) << hist->GetEntries() << " ";
    file<< setw(5) << Nb_tag[0] << " ";
    file<< setw(5) << Nb_tag_err[0] << " ";
    file<< setw(5) << Nb_tag[1] << " ";
    file<< setw(5) << Nb_tag_err[1] << std::endl;

    /**********************************************************************/
    //  save the results, make some plots
    sprintf(hist_name, "data_pt%0.3f_eta%0.3f_tag", (*b)[ii][0],  (*b)[ii][2]);
    hist->SetName(hist_name);
    rootfile->cd();
    hist->Write();
    c1->cd();
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
    hist->Draw("PE");
    func->Draw("SAME"); func->SetLineWidth(2);
    func->SetLineColor(kBlack);

    TF1 sig_tag(*func);
    sig_tag.SetParameter("Nc", 0);
    sig_tag.SetLineColor(kRed);sig_tag.SetLineWidth(2);
    sig_tag.Draw("SAME");


    TF1 bkg_tag(*func);
    bkg_tag.SetParameter("Nb", 0);
    bkg_tag.SetLineColor(kBlue);sig_tag.SetLineWidth(2);
    bkg_tag.Draw("SAME");
    TString epsname_tag(hist_name);
    epsname_tag.Append(".eps");
    c1->SaveAs(epsname_tag.Data());
    /**********************************************************************/



    (*eff)[ii].push_back( 1.0*Nb_tag[0]/Nb[0] );
    (*eff)[ii].push_back( effErr(Nb[0], Nb_err[0], Nb_tag[0], Nb_tag_err[0]));



    //  MC efficiency
    double mc_eff, mc_eff_err;

    this->getMCEff( (*b)[ii][0], (*b)[ii][1], (*b)[ii][2], (*b)[ii][3], tagger, mc_eff, mc_eff_err);


    (*eff)[ii].push_back( mc_eff);
    (*eff)[ii].push_back( mc_eff_err);
  }

  rootfile->Write();
  rootfile->Close();

  file.close();
}




void PtrelSolver::getMCEff(double pt_min, double pt_max, double eta_min, double eta_max, TCut tagger, double &eff, double &err) {



  TH1F *data= new TH1F("data", "", hist_bins, x_min, x_max);

  char eta_cut[100];
  sprintf(eta_cut, "abs(jetEta) >= %f && abs(jetEta) < %f", eta_min, eta_max);
  TChain ntuple("ntuple");

  if (fabs(pt_min -30) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb30_50_b.root");
    // ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc30_50_b.root");
  }

  if (fabs(pt_min -50) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb50_80_b.root");
    //ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc50_80_b.root");
  }


  if (fabs(pt_min -80) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb80_120_b.root");
    //ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc80_120_b1.root");
  }

  if (fabs(pt_min -120) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb120_170_b.root");
    //ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc120_170_b.root");
  }

  ntuple.Project(data->GetName(), "pTrel", eta_cut && muon);


  double events = data->GetEntries();

  ntuple.Project(data->GetName(), "pTrel", eta_cut && muon&& tagger);

  eff = data->GetEntries()*1.0/events;
  err = sqrt(data->GetEntries())/events;
}


TH1F *PtrelSolver::getData(double pt_min, double pt_max, double eta_min, double eta_max, TCut cut) {


  TH1F *data= new TH1F("data", "", hist_bins, x_min, x_max);

  char eta_cut[100];
  sprintf(eta_cut, "abs(jetEta) >= %f && abs(jetEta) < %f", eta_min, eta_max);
  TChain ntuple("ntuple");

  if (fabs(pt_min -30) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb30_50_b.root");
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc30_50_b.root");
  }

  if (fabs(pt_min -50) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb50_80_b.root");
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc50_80_b.root");
  }


  if (fabs(pt_min -80) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb80_120_b.root");
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc80_120_b1.root");
  }

  if (fabs(pt_min -120) < 0.001) {
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/bb120_170_b.root");
    ntuple.Add("/uscms_data/d1/ptan/data/Rel131/cc120_170_b.root");
  }

  ntuple.Project(data->GetName(), "pTrel", cut && eta_cut);

  data->GetXaxis()->SetTitle("P_{Trel} (GeV)");
  data->GetYaxis()->SetTitle("Events");

  return data;
}


void PtrelSolver::makeEffTable(std::vector<std::vector<double> > *eff, const char *filename) {
  if (!eff) return;

  file.open(filename, ios::trunc | ios::out);
  file << setw(10) << "#bin_index" << " "
       << setw(10) << "eff_data" << " "
       << setw(10) << "error" << " "
       << setw(10) << "eff_MC" << " "
       << setw(10) << "error" << " "
       << std::endl;

  for (unsigned int ii = 0; ii < (*eff).size(); ii ++) {

    for (unsigned int jj = 0; jj < (*eff)[ii].size(); jj ++) {

      file <<  setw(10) << (*eff)[ii][jj] << " ";
    }
    file << std::endl;
  }

  file.close();
}


//
// ********************************************
// this needs further work on it
//
void PtrelSolver::makeEffHists(std::vector<std::vector<double> > *eff, const char *tag) {
  if (!eff || !tag) return;

  // assuming pt, eta binned in increasing oder.   
  // first) build the binning of pt, eta distributions. 
  int num_pt_bins, num_eta_bins;
  double pt_bin[100], eta_bin[100];

  num_pt_bins = 0;
  num_eta_bins = 0;

  // make the binning for pt.
  for (unsigned int ii = 0; ii < (*eff).size(); ii ++) {

    if (num_pt_bins == 0) {
      pt_bin[num_pt_bins] = (*eff)[ii][0]; 
    }

    if ((*eff)[ii][0] > pt_bin[num_pt_bins])  {
      num_pt_bins ++;
      pt_bin[num_pt_bins] = (*eff)[ii][0];
    }
    if ((*eff)[ii][1] > pt_bin[num_pt_bins])  {
      num_pt_bins ++;
      pt_bin[num_pt_bins] = (*eff)[ii][1];
    }
  }


  // make the binning for eta. 
  for (unsigned int ii = 0; ii < (*eff).size(); ii ++) {

    if (num_eta_bins == 0) {
      eta_bin[num_eta_bins] = (*eff)[ii][2]; 
    }

    if ((*eff)[ii][2] > eta_bin[num_eta_bins])  { 
      num_eta_bins ++;
      eta_bin[num_eta_bins] = (*eff)[ii][2]; 
    }
    if ((*eff)[ii][3] > eta_bin[num_eta_bins])  {
      num_eta_bins ++;
      eta_bin[num_eta_bins] = (*eff)[ii][3];
    }
  }

  for (int ii = 0; ii <= num_pt_bins; ii ++) {

    std::cout << setw(10) << pt_bin[ii] <<"\t";
  }
  std::cout << std::endl;


  for (int ii = 0; ii <= num_eta_bins; ii ++) {

    std::cout << setw(10) << eta_bin[ii] <<"\t";
  }
  std::cout << std::endl;



  //  2.1) make pt dependence for every eta bin.
  TH1F *hist_data, *hist_mc;
  int eff_data_col = 4;
  int eff_mc_col = 6;

  for (int ii = 0; ii < num_eta_bins; ii ++) {
    char hist_name[30];

    // data histogram
    sprintf(hist_name, "%s_pt_eta%.2f_data", tag, eta_bin[ii]);
    std::cout << "histogram name " << hist_name << std::endl;
    hist_data = new TH1F(hist_name, hist_name, num_pt_bins, pt_bin);

    // mc histogram
    sprintf(hist_name, "%s_pt_eta%.2f_mc", tag, eta_bin[ii]);
    std::cout << "histogram name " << hist_name << std::endl;
    hist_mc = new TH1F(hist_name, hist_name, num_pt_bins, pt_bin);

    for (unsigned int jj = 0; jj < (*eff).size(); jj ++) {

      // loop over all entry in efficiency table and select
      // entry in the same eta bin. 
      if (fabs( (*eff)[jj][2] - eta_bin[ii])> 0.01) continue;


      int theBin = this->getBin( (*eff)[jj][0], pt_bin);
      hist_data->SetBinContent( theBin, (*eff)[jj][eff_data_col]); 
      hist_data->SetBinError(   theBin, (*eff)[jj][eff_data_col+1]); 

      hist_mc->SetBinContent( theBin, (*eff)[jj][eff_mc_col]); 
      hist_mc->SetBinError(   theBin, (*eff)[jj][eff_mc_col+1]); 
    }


    hist_data->GetXaxis()->SetTitle("P_{t} (GeV)");
    hist_data->GetYaxis()->SetTitle("Events");
    hist_mc->GetXaxis()->SetTitle("P_{t} (GeV)");
    hist_mc->GetYaxis()->SetTitle("Events");    
    data.AddLast(hist_data);
    mc.AddLast(hist_mc);
  }



  //  2.2) make eta dependence for every pt bin.
  for (int ii = 0; ii < num_pt_bins; ii ++) {
    char hist_name[30];

    // data histogram
    sprintf(hist_name, "%s_eta_pt%4f_data", tag, pt_bin[ii]);
    std::cout << "histogram name " << hist_name << std::endl;
    hist_data = new TH1F(hist_name, hist_name, num_eta_bins, eta_bin);

    // mc histogram
    sprintf(hist_name, "%s_eta_pt%4f_mc", tag, pt_bin[ii]);
    std::cout << "histogram name " << hist_name << std::endl;
    hist_mc = new TH1F(hist_name, hist_name, num_eta_bins, eta_bin);

    for (unsigned int jj = 0; jj < (*eff).size(); jj ++) {

      // loop over all entry in efficiency table and select
      // entry in the same pt. 
      if (fabs( (*eff)[jj][0] - pt_bin[ii])> 0.01) continue;


      int theBin = this->getBin( (*eff)[jj][2], eta_bin);
      hist_data->SetBinContent( theBin, (*eff)[jj][eff_data_col]); 
      hist_data->SetBinError(   theBin, (*eff)[jj][eff_data_col+1]); 

      hist_mc->SetBinContent( theBin, (*eff)[jj][eff_mc_col]); 
      hist_mc->SetBinError(   theBin, (*eff)[jj][eff_mc_col+1]); 

    }
    //    hist_data->SetMA(1.1);
    //hist_data->SetMinimum(0);
    //hist_mc->SetMinimum(0);
    hist_data->GetXaxis()->SetTitle("#eta");
    hist_data->GetYaxis()->SetTitle("Events");
    hist_mc->GetXaxis()->SetTitle("#eta");
    hist_mc->GetYaxis()->SetTitle("Events");
    data.AddLast(hist_data);
    mc.AddLast(hist_mc);
  }


  makeHistEPS(data, mc);
  return;
}


void PtrelSolver::makeHistEPS(TObjArray &data, TObjArray &mc) {
 
  if (data.GetLast() != mc.GetLast()) return;

  c1 = new TCanvas("c1", "", 600, 600);
  TH1F *hist_data, *hist_mc;
  for (int ii = 0; ii < data.GetLast()+1; ii ++) {

    c1->cd();
    hist_data = (TH1F *) data.At(ii);
    hist_mc = (TH1F *) mc.At(ii);
    if (!hist_data || !hist_mc) continue;


    hist_data->SetMinimum(0);
    hist_data->SetMaximum(hist_data->GetMaximum() *1.4);

    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerSize(1.2);
    hist_data->SetMarkerColor(kBlack);
    hist_data->SetLineColor(kBlack);

    hist_mc->SetMarkerStyle(21);
    hist_mc->SetMarkerSize(1.2);
    hist_mc->SetMarkerColor(kRed);
    hist_mc->SetLineColor(kRed);


    hist_data->Draw("PE");
    hist_mc->Draw("SAMEPE");

    TLegend *leg = new TLegend( 0.23, 0.2, 0.4, 0.4);

    leg->AddEntry(hist_data, "Data ", "pl");
    leg->AddEntry(hist_mc, "MC ", "pl");
    

    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->Draw();

    TString eps_name(hist_data->GetName());
    eps_name.Append(".eps");
    c1->SaveAs(eps_name.Data());
  }
}






//
// make the templates depending on the pt/eta bins
void PtrelSolver::makeTemplates(int flavor, const char *filename, const char *rootfilename) {

  if (!filename || !rootfilename) return;

  TFile *rootfile = new TFile(rootfilename, "RECREATE");
  if (!rootfile) return;


  // create output text file
  // and the header comments
  file.open(filename, ios::trunc | ios::out);
  file << setw(8) << "#pt_min" << " "
       << setw(8) << "pt_max" << " "
       << setw(8) << "eta_min" << " "
       << setw(8)<< "eta_max" << " "
       << setw(10) << "const_a" << " "
       << setw(10) << "error" << " "
       << setw(10) << "const_b" << " "
       << setw(10) << "error" << " "
       << setw(10)  << "const_c" << " "
       << setw(10)  << "error" << " "
       << setw(10)  << "const_d" << " "
       << setw(10)  << "error" << " "
       << endl;


  // template function
  switch (flavor) {

  case 5:
    func = new TF1("func", pdf, x_min, x_max, 4); break;
  case 4:  
    func = new TF1("func", pdf1, x_min, x_max, 5); break;
  default: break;
  } 
 if (!func) return;



  for (int ii = 0; ii < pthat_bins; ii ++) {

    // pt cut
    char pt_cut[100];
    sprintf(pt_cut, "(jetPt >= %f && jetPt < %f)", pthatbining[ii], pthatbining[ii+1]);


    //
    // load different data files
    // ?
    // 
    TChain *ntuple = new TChain("ntuple");
    switch (flavor) {

    case 5: // b flavor

      if (ii ==0) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/bb30_50_a.root");
      if (ii ==1) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/bb50_80_a.root");
      if (ii ==2) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/bb80_120_a.root");
      if (ii ==3) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/bb120_170_a.root");
      break;

    case 4: // c flavor
      if (ii ==0) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/cc30_50_a.root");
      if (ii ==1) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/cc50_80_a.root");
      if (ii ==2) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/cc80_120_a1.root");
      if (ii ==3) ntuple->Add("/uscms_data/d1/ptan/data/Rel131/cc120_170_a.root");

      break;

    default: break;
    }



    //
    // right now, no pt cut since it loads from different files. 
    // in same pt_bin, all events in all eta range used to make
    // templates;
    for (int jj = 0; jj < eta_bins; jj ++) {

      // eta cut
      char eta_cut[100];
      sprintf(eta_cut, "(abs(jetEta) >= %f && abs(jetEta) < %f)", etabining[jj], etabining[jj+1]);


      char hist_name[20];
      sprintf(hist_name, "hist_%d", ii*10 + jj);
      hist= new TH1F(hist_name, "", hist_bins, x_min, x_max);
      hist->GetXaxis()->SetTitle("p_{Trel} (GeV)");
      hist->GetYaxis()->SetTitle("Events/GeV");



      //  making template
      //  ntuple->Project(hist_name, "pTrel", muon && eta_cut);
      ntuple->Project(hist_name, "pTrel", muon);
      c1->cd();
      hist->Draw("PE");
      switch (flavor) {

      case 5:
	func->SetParameters(1.30774,-0.71646, 0.0375143, 100); break;
      case 4:
	func->SetParameters(1.30774,-0.71646, 0.0375143, 100, 2); break;
      default:
	break;
      }


      rootfile->cd();
      hist->Write(); 

      hist->Fit(func);


      TString eps_name(hist_name);
      switch (flavor) {

      case 5: eps_name.Append("_b"); break;
      case 4: eps_name.Append("_c"); break;
      default: 
	break;
      }
      eps_name.Append(".eps");
      c1->SaveAs(eps_name);

     

      file << setw(8)<< pthatbining[ii] << " "
	   << setw(8)  << pthatbining[ii+1] << " "
	   << setw(8)  << etabining[jj] << " "
	   << setw(8)	   << etabining[jj+1] << " ";
	
      for (int kk = 0; kk < func->GetNumberFreeParameters(); kk ++) {
	
	file << setw(10)   << func->GetParameter(kk) << " "
	     << setw(10)   << func->GetParError(kk) << " ";

      }
      file << endl;
    }
  }

  rootfile->Write();
  rootfile->Close();
  file.close();
}


/**************************************************************************
 *
 * make templates based on the 2-dimension histograms.
 *
 **************************************************************************/
void PtrelSolver::makeTemplates(int flavor, const char *inputfilename, const char *pdffilename, const char *rootfilename, const char *tag, const char *pthist, const char *etahist, bool latex) {


  if (!pdffilename || !rootfilename ) return;
  std::cout <<"information: " << "make templates for flavor " <<flavor 
	    << " with " << inputfilename
	    << " [" << tag << ", "
	    << pthist << ", "
	    << etahist << "]"
	    << std::endl;

  TFile *inputfile = new TFile(inputfilename, "READ");
  TFile *rootfile = new TFile(rootfilename, "RECREATE");
  if (!inputfile || !rootfile) return;

  file.open(pdffilename, ios::trunc | ios::out);
  

  file << "#";
  file << setw(8) << "pdf index" << " "
       << setw(10) << "const_a" << " "
       << setw(10) << "error" << " "
       << setw(10) << "const_b" << " "
       << setw(10) << "error" << " "
       << setw(10)  << "const_c" << " "
       << setw(10)  << "error" << " "
       << setw(10)  << "const_d" << " "
       << setw(10)  << "error" << " "
       << std::endl;

  if (!rootfile || !inputfile) return;


  char hist_name[100];
  TH1F *hist=0;
  TH2F *hist2=0;
  TF1  *thePdf= 0;
  int   pdf_index;
  int   nbins;
  TString histname;
  histname.Append(pthist);



  // 1) for pt dependence.
  if (flavor == 5) {
    histname.Insert(0, "b_");
    hist2 = (TH2F *)inputfile->Get(histname.Data());

  }  else {
    histname.Insert(0, "cl_");
    hist2 = (TH2F *)inputfile->Get(histname.Data());
  }
  if (!hist2) {

    std::cout << "can't access the data set  ... " << std::endl;
    rootfile->Close();
    inputfile->Close();
    file.close();
    return;
  }

  nbins = hist2->GetNbinsX();
  x_min = hist2->GetYaxis()->GetXmin();
  x_max = hist2->GetYaxis()->GetXmax();


  //  set the fit range for the histograms. 
  fit_min = x_min;
  if (flavor ==5) {
    fit_max = 4.0;
  } else fit_max = 3.;
  //    fit_max = 4;



  //  hist2->Draw();
  //  std::cout << "pTrel range ... [" << x_min << ", " 
  //	    << x_max << "] " << std::endl;
  int total_num =(int)( sqrt(nbins *1.));
  std::cout << "information: " << total_num << " ... " << std::endl;
  TCanvas *c2 = new TCanvas("c2", "", 900, 900);
  c2->Divide( total_num + 1, total_num+1); 

  for (int ii = 0; ii <= nbins; ii++) {

    c1->cd();
    //std::cout << "processing bin " << ii << " ... " << std::endl;
    pdf_index = index(ii, 0);

    sprintf(hist_name, "%s_flavor_%d_%d", tag, flavor, pdf_index);


    if (ii ==0) 
      hist = (TH1F *)hist2->ProjectionY(hist_name, 1, nbins);
    else if (ii >= pt_threshold && pt_sumbins) {

      int lowbin = 1;
      if (ii+1-pt_sumbins>lowbin) lowbin = ii+1-pt_sumbins;
      int highbin = lowbin+pt_sumbins-1;
      if (nbins <highbin) highbin = nbins;

      hist = (TH1F *)hist2->ProjectionY(hist_name, lowbin, highbin);
      std::cout << "information: projecting between [" 
		<< lowbin << ", " 
		<< highbin << "]" << std::endl;
    } else 
      hist = (TH1F *)hist2->ProjectionY(hist_name, ii, ii);


    // template thePdftion
    switch (flavor) {

    case 5:

      thePdf = new TF1("thePdf", pdf1, x_min, x_max, 5); 
      thePdf->SetParameters(1.30774,-0.51646, 0.00475143, 2.1, 800); break;
 
      
    case 4:  

      thePdf = new TF1("thePdf", pdf1, x_min, x_max, 5);
      thePdf->SetParameters(1.30774,-2.51646, 0.00475143, 1.1, 800); break;
      //      thePdf->SetParameters(2.6774,-7.1646, 0.00375143, 1.0, 100); break;
 
    default: break;
    } 
    if (!thePdf) {

      std::cout << "can't initialize the pdf for fitting ... " << std::endl;
      break;
    }


    formatHist1(hist, "p_{T}^{rel} GeV", "Events");
    hist->Fit(thePdf, "Q", "", fit_min, fit_max);
    rootfile->cd();
    hist->Write(); 
    

    TString eps_name(hist_name);
    eps_name.Append(".eps");
    c1->SaveAs(eps_name);


    char text_label[100];
    if (ii !=0) 
      sprintf(text_label, "[%3.f, %3.f] GeV", hist2->GetXaxis()->GetBinLowEdge(ii), hist2->GetXaxis()->GetBinLowEdge(ii) + hist2->GetXaxis()->GetBinWidth(ii) );
    else 
      sprintf(text_label, "[%3.f, %3.f] GeV", hist2->GetXaxis()->GetBinLowEdge(1), hist2->GetXaxis()->GetBinLowEdge(nbins) + hist2->GetXaxis()->GetBinWidth(nbins) );
      
    c2->cd(ii+1); hist->Draw(); label->DrawLatex(0.5, 0.75, text_label);
    c1->cd();


  

    // save fitted parameters into .dat file
    if (latex) {

      if (ii == 0)
	file <<"all pt bins&"  << "\t";
      else {

	file<< fixed<<setprecision(3) << "[" << hist2->GetXaxis()->GetBinLowEdge(ii)
	    << ", " << hist2->GetXaxis()->GetBinLowEdge(ii) + hist2->GetXaxis()->GetBinWidth(ii) 
	    <<"] &\t";
      }

      for (int kk = 0; kk < thePdf->GetNumberFreeParameters()-1; kk ++) {
	  
	file << scientific << setprecision(2) << thePdf->GetParameter(kk) << "$\\pm$"
	     << thePdf->GetParError(kk) << "& "
	     <<"\t";

      }
	  
      file << scientific << setprecision(2) << thePdf->GetParameter(thePdf->GetNumberFreeParameters() -1) 
	   << "$\\pm$"
	   << thePdf->GetParError(thePdf->GetNumberFreeParameters()-1 ) 
	   << "\\\\\\hline";

    } else {
      file << setw(8)<< pdf_index << " ";
      for (int kk = 0; kk < thePdf->GetNumberFreeParameters(); kk ++) {
	
	file << setw(10)   << thePdf->GetParameter(kk) << " "
	     << setw(10)   << thePdf->GetParError(kk) << " ";

      }
    }
    file << std::endl;
  }

  TString total_plot;
  total_plot.Append( (char) (flavor+'0'), 1);
  total_plot.Append("_");
  total_plot.Append(tag);
  total_plot.Append("_");
  total_plot.Append(pthist);
  total_plot.Append(".eps");
  c2->SaveAs(total_plot.Data());


  fit_max = 3.8;


  // 2) for eta dependence.
  histname.Resize(0);
  histname.Append(etahist);
  if (flavor == 5) {
    histname.Insert(0, "b_");
    hist2 = (TH2F *)inputfile->Get(histname.Data());

  }  else {
    histname.Insert(0, "cl_");
    hist2 = (TH2F *)inputfile->Get(histname.Data());
  }
  if (!hist2) {

    std::cout << "can't access the data set  ... " << std::endl;
    rootfile->Close();
    inputfile->Close();
    file.close();
    return;
  }

  nbins = hist2->GetNbinsX();
  x_min = hist2->GetYaxis()->GetXmin();
  x_max = hist2->GetYaxis()->GetXmax();



  total_num =(int)( sqrt(nbins *1.));
  c2 = new TCanvas("c2", "", 900, 900);
  c2->Divide( total_num + 1, total_num+1); 
 
  //  hist2->Draw();
  // std::cout << "pTrel range ... [" << x_min << ", " 
  //	    << x_max << "] " << std::endl;

  for (int ii = 1; ii <= nbins; ii++) {

    c1->cd();

    //    std::cout << "processing bin " << ii << " ... " << std::endl;
    pdf_index = index(0, ii);

    sprintf(hist_name, "%s_flavor_%d_%d", tag, flavor, pdf_index);

    if (ii >=eta_threshold && eta_sumbins) {
      
      int lowbin = 1;
      if (ii+1-eta_sumbins>lowbin) lowbin = ii+1-eta_sumbins;
      int highbin = lowbin + eta_sumbins-1;
      if (nbins <highbin) highbin = nbins;

      hist = (TH1F *)hist2->ProjectionY(hist_name, lowbin, highbin);
      std::cout << "information: projecting between [" 
		<< lowbin << ", " 
		<< highbin << "]" << std::endl;

    }  else 
      hist = (TH1F *)hist2->ProjectionY(hist_name, ii, ii);

    // template thePdftion
    switch (flavor) {

    case 5:

      thePdf = new TF1("thePdf", pdf1, x_min, x_max, 5); 
      thePdf->SetParameters(1.30774,-0.91646, 0.0175143, 2, 100); break;

    case 4:  

      thePdf = new TF1("thePdf", pdf1, x_min, x_max, 5);
      thePdf->SetParameters(1.30774,-0.71646, 0.0375143, 1.5, 100); break;

    default: break;
    } 
    if (!thePdf) {

      std::cout << "can't initialize the pdf for fitting ... " << std::endl;
      break;
    }

    formatHist1(hist, "p_{T}^{rel} GeV", "Events");
    hist->Fit(thePdf, "Q", "", fit_min, fit_max);
    rootfile->cd();
    hist->Write(); 
    
    hist->Draw();
    TString eps_name(hist_name);
    eps_name.Append(".eps");
    c1->SaveAs(eps_name);
  

    c2->cd(ii); hist->Draw();
    char text_label[100];
    sprintf(text_label, "[%1.2f, %1.2f]", hist2->GetXaxis()->GetBinLowEdge(ii), hist2->GetXaxis()->GetBinLowEdge(ii) + hist2->GetXaxis()->GetBinWidth(ii) );

    c2->cd(ii); hist->Draw(); label->DrawLatex(0.5, 0.75, text_label);

    c1->cd();
  
    // save fitted parameters into .dat file
    if (latex) {



      file<< fixed << setprecision(3) << "[" << hist2->GetXaxis()->GetBinLowEdge(ii)
	  << ", " << hist2->GetXaxis()->GetBinLowEdge(ii) + hist2->GetXaxis()->GetBinWidth(ii) 
	  <<"] &\t";
      

      for (int kk = 0; kk < thePdf->GetNumberFreeParameters()-1; kk ++) {
	
	file << scientific << setprecision(2) << thePdf->GetParameter(kk) << "$\\pm$"
	     << thePdf->GetParError(kk) << "& "
	     <<"\t";

      }
	  
      file << scientific << setprecision(2) << thePdf->GetParameter(thePdf->GetNumberFreeParameters() -1) 
	   << "$\\pm$"
	   << thePdf->GetParError(thePdf->GetNumberFreeParameters()-1 ) 
	   << "\\\\\\hline";
    } else {

      file << setw(8)<< pdf_index << " ";
      for (int kk = 0; kk < thePdf->GetNumberFreeParameters(); kk ++) {
	
	file << setw(10)   << thePdf->GetParameter(kk) << " "
	     << setw(10)   << thePdf->GetParError(kk) << " ";

      }
    }
    file << std::endl;
  }
  total_plot.Resize(0);
  total_plot.Append( (char) (flavor+'0'), 1);
  total_plot.Append("_");
  total_plot.Append(tag);
  total_plot.Append("_");
  total_plot.Append(etahist);
  total_plot.Append(".eps");
  c2->SaveAs(total_plot.Data());


  delete c2;
  rootfile->Write();
  rootfile->Close();
  file.close();
}


/**************************************************************************
 *
 *  read templates/build pdfs. 
 *
 *
 *  build the 2-components pdf
 *  the constant of b or c component is "Nb" and "Nc"
 *
 **************************************************************************/
void PtrelSolver::readTemplates(const char *filename, std::vector<std::vector<double> > *parameters) {

  if (!filename) return;
  parameters->clear();

  std::fstream file;
  char buff[256];
  file.open(filename, ios::in);
  int rows = 0, cols= 0;

  while (!file.eof()) {


    cols = 0;

    file.getline(buff, 256, '\n');
    //    std::cout << buff << std::endl;
    TString container(buff);
    if (container.Contains('#')) continue; //skip comments

    TString tmp("");
    for (int ii = 0 ; ii < container.Length(); ii ++) {
      

      if (container[ii] != ' ') tmp.Append(container[ii], 1);
      else {

	if (!tmp.Length()) continue; // skip space
	std::vector<double> newone;
	if (cols == 0) parameters->push_back(newone); // a new line

	(*parameters)[rows].push_back( atof(tmp.Data()) );
	tmp.Resize(0);
	cols ++;
	//std::cout << rows << " " << cols << std::endl;
      }
    }
    if (tmp.Length()) 	{
      (*parameters)[rows].push_back( atof(tmp.Data()));

      //      (*parameters)[rows][cols] = atof(tmp.Data());
    }   

    rows ++;
  }


  // show the decoded results. 
  // for (unsigned int ii = 0; ii < (*parameters).size() ; ii ++) {

  //    for (unsigned int jj = 0; jj < (*parameters)[ii].size(); jj ++) 
  //  std::cout << setw(10) <<(*parameters)[ii][jj] <<" "; 

  //    std::cout << std::endl;
  // }

  file.close();
}


TF1 *PtrelSolver::buildAPdf(int ii, 
		       std::vector<std::vector<double> > *b, 
		       std::vector<std::vector<double> > *c) {

  if ((unsigned int)ii > (*c).size()) {
      
    std::cout << "binning of parameter files is not matched !" << std::endl;
    return 0;
  }
  //  if ((*b)[ii].size() != 9 || (*c)[ii].size() != 11) {

  //    std::cout << "binning of parameter files is corrupted !" << std::endl;
  //  return 0;
  // }

  if ( (int)((*b)[ii][0]) != (int)((*c)[ii][0]) ) {
    
    std::cout << "pdf is being built with different pt eta bins!" << std::endl;
  }


  char pdf_name[50];
  sprintf(pdf_name, "pdf_index_%d", (int)((*b)[ii][0]));
  //  std::cout << "build pdf " << pdf_name << std::endl;


  // two components pdf with 8 parameters. 
  TF1 *ff=0;


  unsigned int shift = 1; // the cols for pt, eta, etc. 
  unsigned int total_params = ((*b)[ii].size() -shift)/2;
  unsigned int total_params_c = ((*c)[ii].size() -shift)/2;
  ff = new TF1(pdf_name, combined_pdf, x_min, x_max, total_params + total_params_c);


  for (unsigned int jj = 0; jj < total_params; jj ++) {

    if ((unsigned int) jj == (total_params -1) ) 
      func->SetParName(jj , "Nb");
    else       
      func->FixParameter(jj,(*b)[ii][jj*2 + shift]); 
  }
  

  for (unsigned int jj = 0; jj < total_params_c; jj ++) {
    
    if ((unsigned int) jj == (total_params_c -1) ) 
      func->SetParName(jj+total_params , "Nc"); 
    else 
      func->FixParameter(jj+total_params,(*c)[ii][jj*2 + shift]); 
  }


  return ff;
}


//  build the 2-components pdf
//  the constant of b or c component is "Nb" and "Nc"
void PtrelSolver::buildPdfs(TObjArray *combined, 
		       std::vector<std::vector<double> > *b, 
		       std::vector<std::vector<double> > *c,
		       const char *tag) {
  


  for (unsigned int ii = 0; ii < (*b).size(); ii ++) {
    if ((unsigned int)ii > (*c).size()) {
      
      std::cout << "binning of parameter files is not matched !" << std::endl;
      return;
    }

    if ( (int)((*b)[ii][0]) != (int)((*c)[ii][0]) ) {

      std::cout << "pdf is being built with different pt eta bins!" << std::endl;
    }


    char pdf_name[50];
    if (tag)    sprintf(pdf_name, "%s_%d", tag, (int)((*b)[ii][0]));
    else 
      sprintf(pdf_name, "%d", (int)((*b)[ii][0]));

    //    std::cout << "build pdf " << pdf_name << std::endl;


    unsigned int shift = 1; // the cols for pt, eta, etc. 
    unsigned int total_params = ((*b)[ii].size() -shift) /2;
    unsigned int total_params_c = ((*c)[ii].size() -shift) /2;


    // two components pdf
    // b flavor component first,
    // followed by c flavor
    func = new TF1(pdf_name, combined_pdf, x_min, x_max, total_params + total_params_c);


    for (unsigned int jj = 0; jj < total_params; jj ++) {

      if ((unsigned int) jj == (total_params -1) ) 
	func->SetParName(jj , "Nb");
      else       
	func->FixParameter(jj,(*b)[ii][jj*2 + shift]); 
    }


    for (unsigned int jj = 0; jj < total_params_c; jj ++) {

      if ((unsigned int) jj == (total_params_c -1) ) 
	func->SetParName(jj+total_params , "Nc"); 
      else 
	func->FixParameter(jj+total_params,(*c)[ii][jj*2 + shift]); 
    }
    
    combined->AddLast(func);
  }
}


void PtrelSolver::initPdfs(const char *b_pdf,  const char *c_pdf, TObjArray *combined, const char *tag) {

  std::vector<std::vector<double> > *b   = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *c   = new std::vector<std::vector<double> >;


  combined->Clear();

  this->readTemplates(b_pdf, b );
  this->readTemplates(c_pdf, c );

  this->buildPdfs(combined, b, c, tag);

  delete b;
  delete c;
}


void PtrelSolver::initPdfs(const char *b_pdf, const char *c_pdf) {


  pdfs_b->clear();
  pdfs_c->clear();
  combined_pdfs.Clear();

  this->readTemplates(b_pdf, pdfs_b );
  this->readTemplates(c_pdf, pdfs_c );

  this->buildPdfs(&combined_pdfs, pdfs_b, pdfs_c);
}


/**************************************************************************
 *
 *  main function to fit data given a pdf.
 *  fit a given data distribution in "data" at the range of "xlow, xhigh", 
 *  calculate fitted results
 *
 **************************************************************************/
void PtrelSolver::Fit(TH1F *data, 
		      TF1 *pdf, 
		      std::vector<double> *num, 
		      std::vector<double> *num_err) {
  
  if (!data || !pdf || !num || !num_err) return;


  num->clear();
  num_err->clear();


  x_min = data->GetXaxis()->GetXmin();
  x_max = data->GetXaxis()->GetXmax();

  //  hard coded range !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  fit_min = x_min;
  fit_max = 4.1;

  double tmp_b = pdf->GetParameter("Nb");
  double tmp_c = pdf->GetParameter("Nc");
  pdf->SetParameter("Nb", 1);
  pdf->SetParameter("Nc", 0);
  Double_t area_b = pdf->Integral(x_min, x_max);


  pdf->SetParameter("Nb", 0);
  pdf->SetParameter("Nc", 1);
  Double_t area_c = pdf->Integral(x_min, x_max);

  pdf->SetParameter("Nb", tmp_b);
  pdf->SetParameter("Nc", tmp_c);


  data->Fit(pdf, "Q", "", fit_min, fit_max);
  double total_events = data->GetEntries();


  (*num).push_back(pdf->GetParameter("Nb") * area_b);
  (*num).push_back(pdf->GetParameter("Nc") * area_c);
  (*num_err).push_back( pdf->GetParError(pdf->GetParNumber("Nb")) * area_b );
  (*num_err).push_back( pdf->GetParError(pdf->GetParNumber("Nc")) * area_c );

  double sum = 0;
  for (int unsigned ii = 0; ii <  (*num).size(); ii ++) sum += (*num)[ii];
  double scale = total_events / sum;

  for (unsigned int ii = 0; ii <  (*num).size(); ii ++) {

    (*num)[ii] = (*num)[ii]*scale;
    (*num_err)[ii] = (*num_err)[ii]*scale;
  }

  return;
}



void PtrelSolver::make(bool sys) {
  

  setPtAverage(1, 1);
  setEtaAverage(1, 1);


  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "b_flavor_notag.data", "b_pdf_notag.root", "no_tag", "ppT", "pEta");



  //makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "b_flavor_TCL.data", "b_pdf_TCL.root", "TCL", "pcmbpT", "pcmbEta");
  //makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPL_AwayTCL.root", "b_flavor_TPL.data", "b_pdf_TPL.root", "TPL", "pcmbpT", "pcmbEta");


  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCM_AwayTCL.root", "b_flavor_TCM.data", "b_pdf_TCM.root", "TCM", "pcmbpT", "pcmbEta");
  //makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPM_AwayTCL.root", "b_flavor_TPM.data", "b_pdf_TPM.root", "TPM", "pcmbpT", "pcmbEta");


  //  setPtAverage(1, 1);
  // setEtaAverage(1, 1);
  //makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCT_AwayTCL.root", "b_flavor_TCT.data", "b_pdf_TCT.root", "TCT", "pcmbpT", "pcmbEta");
  //makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPT_AwayTCL.root", "b_flavor_TPT.data", "b_pdf_TPT.root", "TPT", "pcmbpT", "pcmbEta");


  //  fit_max = 3;

  // c flavor
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "c_flavor_notag.data", "c_pdf_notag.root", "no_tag", "ppT", "pEta");


  //  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "c_flavor_TCL.data", "c_pdf_TCL.root", "TCL", "pcmbpT", "pcmbEta");
  //makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPL_AwayTCL.root", "c_flavor_TPL.data", "c_pdf_TPL.root", "TPL", "pcmbpT", "pcmbEta");


  setPtAverage(1, 1);
  setEtaAverage(1, 1);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCM_AwayTCL.root", "c_flavor_TCM.data", "c_pdf_TCM.root", "TCM", "pcmbpT", "pcmbEta");
  // makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPM_AwayTCL.root", "c_flavor_TPM.data", "c_pdf_TPM.root", "TPM", "pcmbpT", "pcmbEta");


  return;

  setPtAverage(1, 1);
  setEtaAverage(1, 1);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCT_AwayTCL.root", "c_flavor_TCT.data", "c_pdf_TCT.root", "TCT", "pcmbpT", "pcmbEta");
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPT_AwayTCL.root", "c_flavor_TPT.data", "c_pdf_TPT.root", "TPT", "pcmbpT", "pcmbEta");



  return;

  if (!sys) {
  
  // make all the templates for the fitting 
  setPtAverage(1, 1);
  setEtaAverage(1, 1);

  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "b_flavor_notag.data", "b_pdf_notag.root");



  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "b_flavor_TCL.data", "b_pdf_TCL.root", "TCL", "ncmbpT", "ncmbEta");
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPL_AwayTCL.root", "b_flavor_TPL.data", "b_pdf_TPL.root", "TPL", "ncmbpT", "ncmbEta");


  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCM_AwayTCL.root", "b_flavor_TCM.data", "b_pdf_TCM.root", "TCM", "ncmbpT", "ncmbEta");
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPM_AwayTCL.root", "b_flavor_TPM.data", "b_pdf_TPM.root", "TPM", "ncmbpT", "ncmbEta");


  setPtAverage(1, 1);
  setEtaAverage(1, 1);
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCT_AwayTCL.root", "b_flavor_TCT.data", "b_pdf_TCT.root", "TCT", "ncmbpT", "ncmbEta");
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPT_AwayTCL.root", "b_flavor_TPT.data", "b_pdf_TPT.root", "TPT", "ncmbpT", "ncmbEta");


  //  fit_max = 3;

  // c flavor
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "c_flavor_notag.data", "c_pdf_notag.root");


  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "c_flavor_TCL.data", "c_pdf_TCL.root", "TCL", "ncmbpT", "ncmbEta");
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPL_AwayTCL.root", "c_flavor_TPL.data", "c_pdf_TPL.root", "TPL", "ncmbpT", "ncmbEta");


  setPtAverage(1, 1);
  setEtaAverage(1, 1);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCM_AwayTCL.root", "c_flavor_TCM.data", "c_pdf_TCM.root", "TCM", "ncmbpT", "ncmbEta");
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPM_AwayTCL.root", "c_flavor_TPM.data", "c_pdf_TPM.root", "TPM", "ncmbpT", "ncmbEta");


  setPtAverage(1, 1);
  setEtaAverage(1, 1);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCT_AwayTCL.root", "c_flavor_TCT.data", "c_pdf_TCT.root", "TCT", "ncmbpT", "ncmbEta");
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPT_AwayTCL.root", "c_flavor_TPT.data", "c_pdf_TPT.root", "TPT", "ncmbpT", "ncmbEta");


  } else {


  //
  // make all the templates for systematics study
    setPtAverage(8, 2);
    setEtaAverage(1, 1);

    makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCL_AwayTCL_Mu10.root", "b_flavor_sys_notag.data", "b_pdf_sys_notag.root", "sys");

    makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCL_AwayTCL_Mu10.root", "b_flavor_sys_TCL.data", "b_pdf_sys_TCL.root", "sys_TCL", "ncmbpT", "ncmbEta");
    makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TPL_AwayTCL_Mu10.root", "b_flavor_sys_TPL.data", "b_pdf_sys_TPL.root", "sys_TPL", "ncmbpT", "ncmbEta");


    setPtAverage(1, 2);
    setEtaAverage(1, 2);
  makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCM_AwayTCL_Mu10.root", "b_flavor_sys_TCM.data", "b_pdf_sys_TCM.root", "sys_TCM", "ncmbpT", "ncmbEta");
  makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TPM_AwayTCL_Mu10.root", "b_flavor_sys_TPM.data", "b_pdf_sys_TPM.root", "sys_TPM", "ncmbpT", "ncmbEta");

  setPtAverage(1, 3);
  setEtaAverage(1, 3);
  makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCT_AwayTCL_Mu10.root", "b_flavor_sys_TCT.data", "b_pdf_sys_TCT.root", "sys_TCT", "ncmbpT", "ncmbEta");
  makeTemplates(5, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TPT_AwayTCL_Mu10.root", "b_flavor_sys_TPT.data", "b_pdf_sys_TPT.root", "sys_TPT", "ncmbpT", "ncmbEta");



  // c flavor
    setPtAverage(7, 2);
    setEtaAverage(1, 1);

  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCL_AwayTCL_Mu10.root", "c_flavor_sys_notag.data", "c_pdf_sys_notag.root", "sys");

  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCL_AwayTCL_Mu10.root", "c_flavor_sys_TCL.data", "c_pdf_sys_TCL.root", "sys_TCL", "ncmbpT", "ncmbEta");
  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TPL_AwayTCL_Mu10.root", "c_flavor_sys_TPL.data", "c_pdf_sys_TPL.root", "sys_TPL", "ncmbpT", "ncmbEta");


  setPtAverage(1, 3);
  setEtaAverage(1, 3);
  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCM_AwayTCL_Mu10.root", "c_flavor_sys_TCM.data", "c_pdf_sys_TCM.root", "sys_TCM", "ncmbpT", "ncmbEta");
  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TPM_AwayTCL_Mu10.root", "c_flavor_sys_TPM.data", "c_pdf_sys_TPM.root", "sys_TPM", "ncmbpT", "ncmbEta");


  setPtAverage(1, 7);
  setEtaAverage(1, 5);
  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TCT_AwayTCL_Mu10.root", "c_flavor_sys_TCT.data", "c_pdf_sys_TCT.root", "sys_TCT", "ncmbpT", "ncmbEta");
 
  makeTemplates(4, "/uscms/home/baites/work/CMSSW/1_3_3/src/RecoBTag/PerformanceMeasurements/test/MuBBbarCCbar_TPT_AwayTCL_Mu10.root", "c_flavor_sys_TPT.data", "c_pdf_sys_TPT.root", "sys_TPT", "ncmbpT", "ncmbEta");
  }
}


void PtrelSolver::makettbar() {

  //
  // make all the templates for systematics study
  setPtAverage(1, 2);
  setEtaAverage(1, 3);

  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCL_AwayTCL.root", "b_flavor_sys_notag.data", "b_pdf_sys_notag.root", "sys", "pcmbpT", "pcmbEta");

  return;

  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCL_AwayTCL.root", "b_flavor_sys_TCL.data", "b_pdf_sys_TCL.root", "sys_TCL", "pcmbpT", "pcmbEta");
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TPL_AwayTCL.root", "b_flavor_sys_TPL.data", "b_pdf_sys_TPL.root", "sys_TPL", "pcmbpT", "pcmbEta");

  setPtAverage(1, 4);
  setEtaAverage(1, 4);
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCM_AwayTCL.root", "b_flavor_sys_TCM.data", "b_pdf_sys_TCM.root", "sys_TCM", "pcmbpT", "pcmbEta");
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TPM_AwayTCL.root", "b_flavor_sys_TPM.data", "b_pdf_sys_TPM.root", "sys_TPM", "pcmbpT", "pcmbEta");

  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCT_AwayTCL.root", "b_flavor_sys_TCT.data", "b_pdf_sys_TCT.root", "sys_TCT", "pcmbpT", "pcmbEta");
  makeTemplates(5, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TPT_AwayTCL.root", "b_flavor_sys_TPT.data", "b_pdf_sys_TPT.root", "sys_TPT", "pcmbpT", "pcmbEta");



  // c flavor
  setPtAverage(1, 3);
  setEtaAverage(1, 3);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCL_AwayTCL.root", "c_flavor_sys_notag.data", "c_pdf_sys_notag.root", "sys", "pcmbpT", "pcmbEta");



  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCL_AwayTCL.root", "c_flavor_sys_TCL.data", "c_pdf_sys_TCL.root", "sys_TCL", "pcmbpT", "pcmbEta");
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TPL_AwayTCL.root", "c_flavor_sys_TPL.data", "c_pdf_sys_TPL.root", "sys_TPL", "pcmbpT", "pcmbEta");


  setPtAverage(1, 5);
  setEtaAverage(1, 5);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCM_AwayTCL.root", "c_flavor_sys_TCM.data", "c_pdf_sys_TCM.root", "sys_TCM", "pcmbpT", "pcmbEta");
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TPM_AwayTCL.root", "c_flavor_sys_TPM.data", "c_pdf_sys_TPM.root", "sys_TPM", "pcmbpT", "pcmbEta");


  setPtAverage(1, 14);
  setEtaAverage(1, 10);
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TCT_AwayTCL.root", "c_flavor_sys_TCT.data", "c_pdf_sys_TCT.root", "sys_TCT", "pcmbpT", "pcmbEta");
 
  makeTemplates(4, "/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_TTbar_TPT_AwayTCL.root", "c_flavor_sys_TPT.data", "c_pdf_sys_TPT.root", "sys_TPT", "ncmbpT", "ncmbEta");

}



void PtrelSolver::test(bool sys) {
  
  
  initPdfs("./tmp/b_flavor_notag.data", "./tmp/c_flavor_notag.data");
  initPdfs("./tmp/b_flavor_TCM.data", "./tmp/c_flavor_TCM.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_notag.data", "c_flavor_sys_notag.data", &combined_pdfs_sys, "sys");
  initPdfs("b_flavor_sys_TCL.data", "c_flavor_sys_TCL.data", &combined_pdfs_sys_tag, "sys_tag");
 

  //  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCL_AwayTCL.root", "TCLoose", "ppT", "pEta", sys);



  //    initPdfs("b_flavor_DTCM.data", "c_flavor_DTCM.data", &combined_pdfs_tag, "tag");
  //  initPdfs("b_flavor_sys_TCM.data", "c_flavor_sys_TCM.data", &combined_pdfs_sys_tag, "sys_tag");
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCM_AwayTCL.root", "TCMedium",  "ppT", "pEta", sys);
  //  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCT_AwayTCL.root", "TCTight", "ppT", "pEta", sys);



  return;
  initPdfs("b_flavor_DTPL.data", "c_flavor_DTPL.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TPL.data", "c_flavor_sys_TPL.data", &combined_pdfs_sys_tag, "sys_tag");
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TPL_AwayTCL.root", "TPLoose", "ppT", "pEta",sys);



  // making the measurement
  // 2) medium TCL tagger
  initPdfs("b_flavor_DTPM.data", "c_flavor_DTPM.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TPM.data", "c_flavor_sys_TPM.data", &combined_pdfs_sys_tag, "sys_tag");
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TPM_AwayTCL.root", "TPMedium","ppT", "pEta", sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TPT_AwayTCL.root", "TPTight","ppT", "pEta", sys);


  return;
  





  // making the measurement
  // 1) loose TCL tagger
  setPtAverage(1, 1);
  setEtaAverage(1, 1);

  initPdfs("b_flavor_notag.data", "c_flavor_notag.data");
  initPdfs("b_flavor_TCL.data", "c_flavor_TCL.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_notag.data", "c_flavor_sys_notag.data", &combined_pdfs_sys, "sys");
  initPdfs("b_flavor_sys_TCL.data", "c_flavor_sys_TCL.data", &combined_pdfs_sys_tag, "sys_tag");
  
  // measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCL_AwayTCL.root", "MuBBbarCCbar_TCL", "npT", "nEta", sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCL_AwayTCL.root", "TCLoose", "npT", "nEta", sys);


  // making the measurement
  // 2) medium TCL tagger

  initPdfs("b_flavor_TCM.data", "c_flavor_TCM.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TCM.data", "c_flavor_sys_TCM.data", &combined_pdfs_sys_tag, "sys_tag");
  //measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCM_AwayTCL.root", "MuBBbarCCbar_TCM",  "npT", "nEta", sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCM_AwayTCL.root", "TCMedium",  "npT", "nEta", sys);


  // making the measurement
  // 3) tight TC tagger
  initPdfs("b_flavor_TCT.data", "c_flavor_TCT.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TCT.data", "c_flavor_sys_TCT.data", &combined_pdfs_sys_tag, "sys_tag");

  // measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TCT_AwayTCL.root", "MuBBbarCCbar_TCT",  "npT", "nEta",sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCT_AwayTCL.root", "TCTight", "npT", "nEta", sys);


  initPdfs("b_flavor_TPL.data", "c_flavor_TPL.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TPL.data", "c_flavor_sys_TPL.data", &combined_pdfs_sys_tag, "sys_tag");

  //  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPL_AwayTCL.root", "MuBBbarCCbar_TPL", "npT", "nEta", sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TPL_AwayTCL.root", "TPLoose", "npT", "nEta",sys);



  // making the measurement
  // 2) medium TCL tagger
  initPdfs("b_flavor_TPM.data", "c_flavor_TPM.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TPM.data", "c_flavor_sys_TPM.data", &combined_pdfs_sys_tag, "sys_tag");
  // measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPM_AwayTCL.root", "MuBBbarCCbar_TPM","npT", "nEta", sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TPM_AwayTCL.root", "TPMedium","npT", "nEta", sys);



  // making the measurement
  // 3) tight TC tagger
  initPdfs("b_flavor_TPT.data", "c_flavor_TPT.data", &combined_pdfs_tag, "tag");
  initPdfs("b_flavor_sys_TPT.data", "c_flavor_sys_TPT.data", &combined_pdfs_sys_tag, "sys_tag");
  // measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_MuBBbarCCbar_TPT_AwayTCL.root", "MuBBbarCCbar_TPT","npT", "nEta", sys);
  measure("/uscms_data/d1/lpcbtag/yumiceva/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TPT_AwayTCL.root", "TPTight","npT", "nEta", sys);



}


/**************************************************************************
 *
 *  calibrate the fitting template, check for statistics bias, etc. 
 *  mainly try to see if we can extract the fraction of b-jets in a 
 *  steps:        points between 0, 1
 *  total_events: MC toy sample size
 *
 **************************************************************************/
TGraphErrors *PtrelSolver::checkPurity2(TF1 *pdf, int steps, int num_bs, bool verbose) {

  if (!pdf) return 0;


  double tmp_b = pdf->GetParameter("Nb");
  double tmp_c = pdf->GetParameter("Nc");
  pdf->SetParameter("Nb", 1);
  pdf->SetParameter("Nc", 0);
  Double_t area_b = pdf->Integral(x_min, x_max);


  pdf->SetParameter("Nb", 0);
  pdf->SetParameter("Nc", 1);
  Double_t area_c = pdf->Integral(x_min, x_max);

  pdf->SetParameter("Nb", tmp_b);
  pdf->SetParameter("Nc", tmp_c);

  if (verbose) {

      std::cout << " Sb = " << area_b << "\t";
      std::cout << " Sc = " << area_c << std::endl;
  }


  // generate MC toys
  TF1 signal, bkg;

  Double_t *purity     = new Double_t[steps+1];
  Double_t *purity_err = new Double_t[steps+1];


  Double_t *Nb     = new Double_t[steps+1];
  Double_t *Nc     = new Double_t[steps+1];
  Double_t *Nb_err = new Double_t[steps+1];
  Double_t *Nc_err = new Double_t[steps+1];


  Double_t *Nb_gen     = new Double_t[steps+1];
  Double_t *Nb_gen_err = new Double_t[steps+1];
  Double_t *Nc_gen     = new Double_t[steps+1];
  Double_t *Nc_gen_err = new Double_t[steps+1];



  TCanvas *c1, *c2, *c3;
  TLegend *leg;
  if (verbose) {
    
    c1 = new TCanvas("c1", "", 600, 600);
    c2 = new TCanvas("c2", "", 600, 600);
    c3 = new TCanvas("c3", "", 600, 600);

    leg = new TLegend(0.5, 0.5, 0.9, 0.8);
  }



  // generated data.
  // fit with template
  // check the linearity
  for (int ii = 0; ii < steps; ii ++) {


    //  purity is (ii+1)/steps
    purity[ii] = (ii+1)*1.0/steps;
    purity_err[ii] = 0;

    Nb[ii] = num_bs;
    Nc[ii] = (Int_t)(num_bs * steps*1.0 / (ii+1) - num_bs);

    Double_t total_events = Nb[ii] + Nc[ii];


    Nb_gen[ii] = Nb[ii];    
    Nb_gen_err[ii] = TMath::Sqrt(Nb[ii]);
    Nc_gen[ii] = Nc[ii];    
    Nc_gen_err[ii] = TMath::Sqrt(Nc[ii]);
    pdf->SetParameter("Nb", (Nb[ii]/total_events)  / area_b);
    pdf->SetParameter("Nc", (Nc[ii]/total_events)  / area_c);



    if (verbose) {

      std::cout << "Nb * Sb = " << Nb[ii] << "\t";
      std::cout << "Nc * Sc = " << Nc[ii] << std::endl;
      c1->cd();
      pdf->Draw();
      c1->SaveAs("pdf.eps");
      //c1->cd();
    }


    // generate data
    hist = new TH1F("data", "", hist_bins, x_min, x_max);
    hist->FillRandom(pdf->GetName(), (int)total_events);
    hist->Fit(pdf);

    Nb[ii] = pdf->GetParameter("Nb") * area_b;
    Nb_err[ii] = pdf->GetParError(pdf->GetParNumber("Nb")) * area_b;

    Nc[ii] = pdf->GetParameter("Nc") * area_c;
    Nc_err[ii] = pdf->GetParError(pdf->GetParNumber("Nc")) * area_c;

    double scale = total_events / ( Nb[ii] + Nc[ii]);
    // normalize
    Nb[ii]  = scale * Nb[ii];
    Nb_err[ii]  = scale * Nb_err[ii];
    Nc[ii]  = scale * Nc[ii];
    Nc_err[ii]  = scale * Nc_err[ii];


    if (verbose) {
    

      std::cout << "purity " << purity[ii] << "; Nb" << Nb[ii] << std::endl; 

      c2->cd();
      hist->GetYaxis()->SetLimits(-10, hist->GetMaximum() * 1.1);
      hist->GetXaxis()->SetTitle("pTrel (GeV)");
      hist->GetYaxis()->SetTitle("Events");
      hist->GetYaxis()->SetTitleOffset(1.4);
      hist->Draw("PE");
      
      


      //  make a plot
      pdf->SetLineColor(kBlack);
      pdf->Draw("SAME");


      bkg = TF1(*pdf);
      bkg.SetParameter("Nb", 0);
      bkg.SetLineColor(kBlue);
      bkg.Draw("SAME");


      signal = TF1(*pdf);
      signal.SetParameter("Nc", 0);
      signal.SetLineColor(kRed);
      signal.Draw("SAME");
    

      if (ii == steps/2) {

	leg->AddEntry(pdf, "Total Fit", "l");
	leg->AddEntry(&bkg, "c-jet ", "l");
	leg->AddEntry(&signal, "b-jet ", "l");
	
	c2->cd();
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->Draw();
	c2->SaveAs("example.eps");
      };
    }// end of plot in case verbose mode.


    Nb_err[ii] = Nb_err[ii]/Nb[ii];
    Nb[ii] = 1.0;
  }


  // make histograms
  TGraphErrors *result = new TGraphErrors(steps, purity, Nb, purity_err, Nb_err);

  if (verbose) {

    result->SetLineColor(kRed);
    result->SetMarkerColor(kRed);
    result->GetXaxis()->SetTitle("Purity of Sample");
    result->GetYaxis()->SetTitle("Fitted Num b jets");
    result->GetYaxis()->SetTitleOffset(1.4);
    c3->cd();
   
    result->Draw("AP");

    c3->SaveAs("purity.eps");
  }
  
  return result;
}


//  calibrate the fitting template, check for statistics bias, etc. 
//  mainly try to see if we can extract the fraction of b-jets in a 
//  steps:        points between 0, 1
//  total_events: MC toy sample size
TGraphErrors *PtrelSolver::checkLinearity2(TF1 *pdf, int steps, int total_events, bool verbose) {

  if (!pdf) return 0;


  double tmp_b = pdf->GetParameter("Nb");
  double tmp_c = pdf->GetParameter("Nc");
  pdf->SetParameter("Nb", 1);
  pdf->SetParameter("Nc", 0);
  Double_t area_b = pdf->Integral(x_min, x_max);


  pdf->SetParameter("Nb", 0);
  pdf->SetParameter("Nc", 1);
  Double_t area_c = pdf->Integral(x_min, x_max);

  pdf->SetParameter("Nb", tmp_b);
  pdf->SetParameter("Nc", tmp_c);

  if (verbose) {

      std::cout << " Sb = " << area_b << "\t";
      std::cout << " Sc = " << area_c << std::endl;
  }


  // generate MC toys
  TF1 signal, bkg;

  Double_t *Nb     = new Double_t[steps+1];
  Double_t *Nc     = new Double_t[steps+1];
  Double_t *Nb_err = new Double_t[steps+1];
  Double_t *Nc_err = new Double_t[steps+1];


  Double_t *Nb_gen     = new Double_t[steps+1];
  Double_t *Nb_gen_err = new Double_t[steps+1];
  Double_t *Nc_gen     = new Double_t[steps+1];
  Double_t *Nc_gen_err = new Double_t[steps+1];


  // make some plots for linearity check. 
  TCanvas *c1, *c2, *c3;
  TLegend *leg;
  if (verbose) {
    
    c1 = new TCanvas("c1", "", 600, 600);
    c2 = new TCanvas("c2", "", 600, 600);
    c3 = new TCanvas("c3", "", 600, 600);

    leg = new TLegend(0.5, 0.5, 0.9, 0.8);
  }



  // generated data.
  // fit with template
  // check the linearity
  for (int ii = 0; ii <= steps; ii ++) {


    
    Nb[ii] = 1.0 * ii/steps * total_events;
    Nc[ii] = total_events - Nb[ii];

    Nb_gen[ii] = Nb[ii];    
    Nb_gen_err[ii] = TMath::Sqrt(Nb[ii]);
    Nc_gen[ii] = Nc[ii];    
    Nc_gen_err[ii] = TMath::Sqrt(Nc[ii]);
    pdf->SetParameter("Nb", (Nb[ii]/total_events)  / area_b);
    pdf->SetParameter("Nc", (Nc[ii]/total_events)  / area_c);



    if (verbose) {

      std::cout << "Nb * Sb = " << Nb[ii] << "\t";
      std::cout << "Nc * Sc = " << Nc[ii] << std::endl;
      c1->cd();
      pdf->Draw();
      c1->SaveAs("pdf.eps");
      //c1->cd();
    }


    // generate data
    hist = new TH1F("data", "", hist_bins, x_min, x_max);
    hist->FillRandom(pdf->GetName(), total_events);
    //    hist->Sumw2();
    // hist->Scale(1.0 * total_events/hist->Integral(x_min, x_max));
    hist->Fit(pdf);

    Nb[ii] = pdf->GetParameter("Nb") * area_b;
    Nb_err[ii] = pdf->GetParError(pdf->GetParNumber("Nb")) * area_b;

    Nc[ii] = pdf->GetParameter("Nc") * area_c;
    Nc_err[ii] = pdf->GetParError(pdf->GetParNumber("Nc")) * area_c;

    double scale = total_events / ( Nb[ii] + Nc[ii]);
    // normalize
    Nb[ii]  = scale * Nb[ii];
    Nb_err[ii]  = scale * Nb_err[ii];
    Nc[ii]  = scale * Nc[ii];
    Nc_err[ii]  = scale * Nc_err[ii];


    if (verbose) {
    


      c2->cd();
      hist->GetYaxis()->SetLimits(-10, hist->GetMaximum() * 1.1);
      hist->GetXaxis()->SetTitle("pTrel (GeV)");
      hist->GetYaxis()->SetTitle("Events");
      hist->GetYaxis()->SetTitleOffset(1.4);
      hist->Draw("PE");
      
      


      //  make a plot
      pdf->SetLineColor(kBlack);
      pdf->Draw("SAME");


      bkg = TF1(*pdf);
      bkg.SetParameter("Nb", 0);
      bkg.SetLineColor(kBlue);
      bkg.Draw("SAME");


      signal = TF1(*pdf);
      signal.SetParameter("Nc", 0);
      signal.SetLineColor(kRed);
      signal.Draw("SAME");
    

      if (ii == steps/2) {

	leg->AddEntry(pdf, "Total Fit", "l");
	leg->AddEntry(&bkg, "c-jet ", "l");
	leg->AddEntry(&signal, "b-jet ", "l");
	
	c2->cd();
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->Draw();
	c2->SaveAs("example.eps");
      };
    }// end of plot in case verbose mode.
  }


  // make histograms
  TGraphErrors *result = new TGraphErrors(steps+1, Nb_gen, Nb, Nb_gen_err, Nb_err);

  if (verbose) {

    result->SetLineColor(kRed);
    result->SetMarkerColor(kRed);
    result->GetXaxis()->SetTitle("Expected Num b jets");
    result->GetYaxis()->SetTitle("Fitted Num b jets");
    c3->cd();
   
    result->Draw("AP");


    //    TGraphErrors *result1 = new TGraphErrors(steps+1, Nc_gen, Nc, Nc_gen_err, Nc_err);
    //result1->SetMarkerColor(kRed);
    //result1->SetLineColor(kBlack);
    //result1->Draw("SAMEAP");
    //c1->Update();

    c3->SaveAs("linear.eps");
  }
  
  return result;
}



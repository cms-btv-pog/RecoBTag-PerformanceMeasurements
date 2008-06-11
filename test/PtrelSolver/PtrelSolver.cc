#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2F.h"
#include "TROOT.h"


#include "math.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include <vector>

#include "PtrelSolver.h"
#include "analysis.h"

ClassImp(PtrelSolver)



using namespace std;


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

  gROOT->SetStyle("Plain");

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
  else sprintf(index, "_%d", jj);
  std::cout << "information: accessing pdf with name " << index << std::endl;

  for (int ii = 0; ii <= list->GetLast(); ii ++) {

    TString pdfname(  ((TF1 *)list->At(ii))->GetName()  );
    if (!pdfname.CompareTo(index)) return (TF1 *) list->At(ii);
  }

  return 0;
}


TF1 *PtrelSolver::getPdfByIndex(int jj) { return this->getPdfByIndex( &combined_pdfs, jj);}

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

  // normalized chi2
  (*eff).push_back( Nb[2] );
  (*eff).push_back( Nb_tag[2] );

  return;
}


//
//  make results from root histograms
void PtrelSolver::effCal(TH1F *hist, TH1F *hist_tag, TF1  *pdf, std::vector<double> *eff, const char *rootfilename){

  this->effCal(hist, hist_tag, pdf, pdf, eff, rootfilename);
}
void PtrelSolver::effCal(TH1F *hist, TH1F *hist_tag,  TF1  *pdf, TF1  *pdf_tag, std::vector<double> *eff, const char *rootfilename){

  if (!eff || !pdf || !pdf_tag || !hist || !hist_tag || !rootfilename) return;
  TFile *rootfile = new TFile(rootfilename, "UPDATE");


  std::vector<double>  Nb, Nb_err;
  std::vector<double>  Nb_tag, Nb_tag_err;
  char hist_name[100];


  //  a) before applying tagging
  Fit(hist, pdf, &Nb, &Nb_err);
  sprintf(hist_name, "%s", pdf->GetName());
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
  Fit(hist_tag, pdf_tag, &Nb_tag, &Nb_tag_err);
  sprintf(hist_name, "%s_tag", pdf_tag->GetName());
  hist->SetName(hist_name);
  rootfile->cd();
  hist->Write();
  c1->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->Draw("PE");
  pdf_tag->Draw("SAME"); pdf_tag->SetLineWidth(2);
  pdf_tag->SetLineColor(kBlack);

  TF1 sig_tag(*pdf_tag);
  sig_tag.SetParameter("Nc", 0);
  sig_tag.SetLineColor(kRed);sig_tag.SetLineWidth(2);
  sig_tag.Draw("SAME");


  TF1 bkg_tag(*pdf_tag);
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


TH1F *PtrelSolver::getMCeff(TFile *file, const char *dir, const char *hist, const char *tagger) {
  if (!file || !hist) return 0;

  TString raw(dir); raw.Append(hist); raw.Append("_b");

  TString tag(raw); tag += "_"; tag += tagger;
  tag.Insert(strlen(dir) + 1, "tag");



  TString eff(tagger); 
  eff += "_";
  eff += hist; // base name for the tagger & depending variables
  eff += "_mc";

  std::cout << "information: access MC truth " << raw.Data() << std::endl;
  std::cout << "information: access MC truth " << tag.Data() << std::endl;


  TH2F *raw_hist = (TH2F *)file->Get(raw.Data());
  TH2F *tag_hist = (TH2F *)file->Get(tag.Data());
  if (!raw_hist || !tag_hist) return 0;


  TString tmp(eff); tmp += "_after";
  TH1F *h1 = (TH1F *)tag_hist->ProjectionX(eff.Data());
  TH1F *h2 = (TH1F *)raw_hist->ProjectionX(tmp.Data());
  formatHist1(h1, h1->GetXaxis()->GetTitle(), "Efficiency");

  h1->Sumw2(); h2->Sumw2();
  h1->Divide(h2);
  h1->SetName(eff.Data());

  return h1;
}
TH1F *PtrelSolver::getMCeff(const char *filename, const char *dir, const char *hist, const char *tagger) {

  if (!filename) return 0;
  TFile *mceff = new TFile(filename, "READ");
  return this->getMCeff(mceff, dir, hist, tagger);
}



TH1F  *PtrelSolver::getMCeff(TFile *file, const char *hist, const char *tagger) {

  return this->getMCeff(file, "/Histograms/muon_in_jet/", hist, tagger);
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



/**************************************************************************
 *
 * make templates based on the 2-dimension histograms.
 *
 **************************************************************************/
/*
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




  //  1) for pt dependence
  histname.Append(pthist);
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


  int total_num =(int)( sqrt(nbins *1.));
  std::cout << "information: " << total_num << " ... " << std::endl;
  TCanvas *c2 = new TCanvas("c2", "", 900, 900);
  c2->Divide( total_num + 1, total_num+1); 

  for (int ii = 0; ii <= nbins; ii++) {

    c1->cd();
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
    std::cout << "information: " << hist->GetName() << " chi2/Ndof =" 
	      << thePdf->GetChisquare() << "//"
	      << thePdf->GetNDF()
	      << std::endl;
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
    std::cout << "information: " << hist->GetName() << " chi2/Ndof =" 
	      << thePdf->GetChisquare() << "//"
	      << thePdf->GetNDF()
	      << std::endl;
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
*/


// "b"
// "cl"
// thehistname+flavor
void PtrelSolver::makeTemplates(const char *flavor, const char *sampletag, bool sum, const char *inputfilename, const char *dir, const char *tag, const char *thehist, int pdfbase, const char *outputdir, const char *versiontag, bool sys, bool latex) {

  TString myhist(sampletag); myhist += "_"; myhist += thehist;
  const char *thehistname = myhist.Data();
 
  TString basename;
  basename += tag; basename += "_"; basename += thehistname;
  basename.Insert(0, flavor);

  if (!versiontag) return;
  std::cout <<"information: " << "make templates for flavor " << flavor 
	    << " with " << inputfilename
	    << " [" << tag << "] with indexing from "
	    << pdfbase 
	    << std::endl;

  // sampletag + flavor + "_templates_" +tag+ versiontag+"-sys"
  TString tmp(flavor); tmp.Insert(0, thehistname[0]); tmp += "_templates_"; tmp += tag;
  tmp.Insert(0, outputdir);
  tmp += versiontag;
  if (sys)  tmp += "-sys";
  file.open(tmp.Data(), ios::app | ios::out);

  tmp += ".root";
  TFile *inputfile = new TFile(inputfilename, "READ");
  TFile *rootfile  = new TFile(tmp.Data(),  "UPDATE");
  if (!inputfile || !rootfile) return;

  

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

  char  hist_name[100];
  TH1D *hist=0, *chi2=0;
  TH2D *hist2=0;
  TF1  *thePdf= 0;
  int   pdf_index;
  int   nbins;
  TString histname, chi2name;



  // based on the flavor and tag, determine the histgram used to build pdfs.
  histname += dir;
  histname += "/";
  histname += thehistname;
  histname += "_";
  histname += flavor;
  if (tag && strlen(tag) > 0)  {
    histname.Insert( strlen(dir) +2, "tag" );
    histname += "_";
    histname += tag;
  }

  std::cout << "information: building pdf from " 
	    << histname.Data() << std::endl;


  hist2 = (TH2D *) inputfile->Get(histname.Data());
  if (!hist2) {

    std::cout << "information: can't access the data set  ... " << std::endl;
    rootfile->Close();
    inputfile->Close();
    file.close();
    return;
  }

  nbins = hist2->GetNbinsX();
  x_min = hist2->GetYaxis()->GetXmin();
  x_max = hist2->GetYaxis()->GetXmax();
  chi2name += basename; chi2name += "_chi2";
  chi2 = new TH1D(chi2name.Data(), "", nbins+1, 0, nbins+1);
  formatHist1(chi2, "#chi^{2}", "Bin Index");

  int total_num =(int)( sqrt(nbins *1.));
  //  std::cout << "information: " << total_num << " ... " << std::endl;
  TCanvas *c2 = new TCanvas("c2", "", 900, 900);
  c2->Divide( total_num + 1, total_num+1); 

  // ii = 0 build the total pdf
  for (int ii = (!sum); ii <= nbins; ii++) {

    c1->cd();
    pdf_index = ii*pdfbase;

    sprintf(hist_name, "template_%s_%d", basename.Data(), pdf_index);
    
    if (ii ==0) 
      hist = hist2->ProjectionY(hist_name, 1, nbins, "e");
    else
      hist = hist2->ProjectionY(hist_name, ii, ii, "e");

    std::cout << hist2->GetSumw2N() << std::endl;
    std::cout << hist->GetSumw2N() << std::endl;

    // template thePdftion
    if (!strcmp(flavor, "b")) {

      thePdf = new TF1("thePdf", pdf1, x_min, x_max, 5); 
      thePdf->SetParameters(1.30774,-0.51646, 0.00475143, 2.1, 800);
    } else {

      thePdf = new TF1("thePdf", pdf1, x_min, x_max, 5);
      thePdf->SetParameters(1.30774,-2.51646, 0.00475143, 1.1, 800);
    } 


    if (!thePdf)  std::cout << "information: can't initialize the pdf for fitting ... " << std::endl;
    
    if (pdfbase == PT_BASE)  formatHist1(hist, "p_{T}^{rel} GeV", "Events");
    if (pdfbase == ETA_BASE) formatHist1(hist, "Pseudorapidity",  "Events");

    hist->Fit(thePdf, "Q", "", fit_min, fit_max);
    std::cout << "information: " << hist->GetName() << " chi2/Ndof =" 
	      << thePdf->GetChisquare() << "/"
	      << thePdf->GetNDF()
	      << std::endl;

     chi2->SetBinContent(ii+1,  thePdf->GetChisquare()*1.0/thePdf->GetNDF());
    rootfile->cd();
    hist->Write(); 
    

    char text_label[100];
    if (ii !=0) 
      sprintf(text_label, "[%3.f, %3.f]", hist2->GetXaxis()->GetBinLowEdge(ii), hist2->GetXaxis()->GetBinLowEdge(ii) + hist2->GetXaxis()->GetBinWidth(ii) );
    else 
      sprintf(text_label, "[%3.f, %3.f]", hist2->GetXaxis()->GetBinLowEdge(1), hist2->GetXaxis()->GetBinLowEdge(nbins) + hist2->GetXaxis()->GetBinWidth(nbins) );
      
    c2->cd(ii+1); hist->Draw("E"); label->DrawLatex(0.5, 0.75, text_label);
    // c1->cd();

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


  TString total_plot(basename);
  total_plot += ".eps";
  c2->SaveAs(total_plot.Data());
  delete c2;

  rootfile->cd();
  chi2->Write();
  rootfile->Write();
  rootfile->Close();
  file.close();
}

void PtrelSolver::makeAllTemplatesPerTag(const char *inputfilename, const char *dir, const char *sampletag, const char *tag, const char *outputdir, const char *versiontag, bool sys) {

  // clean the disk area
  TString command("rm -rf "); command += outputdir;
  command += "/";
  command += sampletag;
  command += "_templates_"; command += tag; 
  command += versiontag; 
  if (sys) command += "-sys"; 
  gSystem->Exec(command.Data());
  std::cout << "information:   " << command.Data() << std::endl;

  command += ".root";
  gSystem->Exec(command.Data());
  std::cout << "information:   " << command.Data() << std::endl;


  this->makeTemplates("b", sampletag, true, inputfilename, dir, tag, "pT", PT_BASE, outputdir, versiontag, sys);

  this->makeTemplates("cl", sampletag, true, inputfilename, dir, tag, "pT", PT_BASE, outputdir, versiontag, sys);


  this->makeTemplates("b", sampletag, false, inputfilename, dir, tag, "eta", ETA_BASE, outputdir, versiontag, sys);
  this->makeTemplates("cl", sampletag, false, inputfilename, dir, tag, "eta", ETA_BASE, outputdir, versiontag, sys);

}


//  pdf data file: sampletag + favor+"_templates_" + tag + versiontag.
bool   PtrelSolver::initPdfsByTag(const char *sampletag, const char *tag, const char *pdfdir, const char *versiontag, bool sys) {


  TString count(pdfdir); count += sampletag; 

  TString basename(pdfdir); basename += "/";
  basename += sampletag;
  basename += "_templates_";

  Int_t shift = count.Length() +1;

  TString b_pdf(basename), cl_pdf(basename);
  b_pdf.Insert(  shift, "b");  b_pdf  += versiontag; 
  if (sys) b_pdf  += "-sys";
  cl_pdf.Insert( shift, "cl"); cl_pdf += versiontag; 
  if (sys) cl_pdf += "-sys";

  if ( !locateFile(b_pdf.Data()) || !locateFile(cl_pdf.Data()) ) {
    std::cout << "information: pdf data file " << b_pdf.Data() << " & "
	      << cl_pdf.Data() << " dont exist ... "
	      << std::endl;
    return 0;
  }


  // initialize
  std::cout << "information: " << b_pdf.Data() << " & " << cl_pdf.Data() << std::endl;
  if (sys)  this->initPdfs(b_pdf.Data(), cl_pdf.Data(), "sys");
  else  this->initPdfs(b_pdf.Data(), cl_pdf.Data(), "");


  TString b_pdftag(basename), cl_pdftag(basename);
  b_pdftag.Insert(  shift, "b");  b_pdftag  += tag; b_pdftag  += versiontag;
  cl_pdftag.Insert( shift, "cl"); cl_pdftag += tag; cl_pdftag += versiontag; 
  if (sys) {

    b_pdftag  += "-sys";
    cl_pdftag += "-sys";
  }
  std::cout << "information: " << b_pdftag.Data() << " & " << cl_pdftag.Data() << std::endl;

  if ( !locateFile(b_pdftag.Data()) || !locateFile(cl_pdftag.Data()) ) {

    std::cout << "information: pdf data file " << b_pdftag.Data() << " & "
	      << cl_pdftag.Data() << " dont exist ... "
	      << std::endl;
    if (sys)  this->initPdfs(b_pdf.Data(), cl_pdf.Data(), "sys_tag");
    else  this->initPdfs(b_pdf.Data(), cl_pdf.Data(), "tag");
    std::cout << "information: initialize tag pdfs with untagged pdfs " << std::endl;

  } else {

    if (sys)  this->initPdfs(b_pdftag.Data(), cl_pdftag.Data(), "sys_tag");
    else  this->initPdfs(b_pdftag.Data(), cl_pdftag.Data(), "tag");
  }

  return true;
}

bool PtrelSolver::locateFile(const char *file) {

  TString cmd("ls "); cmd += file;
  Int_t status = gSystem->Exec(cmd.Data());

  if (!status) return true;
  return 0;
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

    file.getline(buff, 512, '\n');
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

  //   for (unsigned int jj = 0; jj < (*parameters)[ii].size(); jj ++) 
  //	std::cout << setw(10) <<(*parameters)[ii][jj] <<" "; 

  //     std::cout << std::endl;
  //   }

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

    std::cout << "information: build pdf --- " << pdf_name << std::endl;


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
void PtrelSolver::initPdfs(const char *b_pdf,  const char *c_pdf, const char *pdftag) {


  TObjArray *pdf = 0;

  if ( !strcmp(pdftag, "") )        pdf = &combined_pdfs;
  if ( !strcmp(pdftag, "tag") )     pdf = &combined_pdfs_tag;
  if ( !strcmp(pdftag, "sys") )     pdf = &combined_pdfs_sys;
  if ( !strcmp(pdftag, "sys_tag") ) pdf = &combined_pdfs_sys_tag;

  if (!pdf) {

    std::cout << "information: can't build the pdfs ... " << std::endl;
    return;
  }

  this->initPdfs(b_pdf, c_pdf, pdf, pdftag);
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


  //  double total_events = data->GetEntries();
  double total_events = data->Integral(x_min, x_max);
  std::cout << "total events" << data->GetEntries() 
	    << ": " << total_events << std::endl;

  double chi2 = 10000, b_frac = 0.1;
  while (chi2 >5 && b_frac < 1.0) {

    std::cout << "information: trying with b fraction of " << b_frac
	      << std::endl;
    pdf->SetParameter("Nb", total_events * b_frac /area_b);
    pdf->SetParameter("Nc", total_events * (1.0 -b_frac) /area_b );


    data->Fit(pdf, "Q", "", fit_min, fit_max);
    chi2 = pdf->GetChisquare()/pdf->GetNDF();
    b_frac += 0.1;
  }

  std::cout << "information: " << data->GetName() << " fitted with chi2/Ndof =" 
	    << pdf->GetChisquare() << "//"
	    << pdf->GetNDF()
	    << std::endl;

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
     
  // REMOVE -----   
  /*
  tmp_b = pdf->GetParameter("Nb");
  tmp_c = pdf->GetParameter("Nc");
  
  data->Draw("e1same");  

  TF1 * pdf2 = new TF1(*pdf);
  
  pdf->SetLineColor(kBlue);
  pdf->SetParameter("Nb", tmp_b);
  pdf->SetParameter("Nc", 0);  	
  pdf->Draw("same");

  pdf2->SetLineColor(kRed);
  pdf2->SetParameter("Nb", 0);
  pdf2->SetParameter("Nc", tmp_c);  	
  pdf2->Draw("same");
  */

  (*num).push_back( pdf->GetChisquare()*1.0/pdf->GetNDF() );

  return;
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






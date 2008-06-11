
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

// Uses the counting method to measure efficiencies
void PtrelSolver::counting(
  const char *inputfilename,
  const char *directory,
  const char *dependency,
  const char *tag,
  const char *outfilename
)
{
  // index base for pdf	
  int pdfbase = 0;
  
  if (!strcmp(dependency,"pT"))
    pdfbase = PT_BASE;
  else if (!strcmp(dependency,"eta"))
    pdfbase = ETA_BASE;
  else
  {
    std::cout << "Unknown dependency" << std::endl;
    return;
  }	
  
  // output file with results
  outfile = new TFile(outfilename, "UPDATE");

  // input file with data.
  inputfile = new TFile(inputfilename, "READ");
 
  char name[255], basename[255];
    
  // histogram with x dependency for muon-in-jet sample
  sprintf(name, "%s/n_%s", directory, dependency);
  TH2F * histogram2D = (TH2F*) inputfile->Get(name);

  // histogram with x dependency for muon-in-jet-away-jet-tagged sample
  sprintf(name, "%s/p_%s", directory, dependency);
  TH2F * histogram2D_0t = (TH2F*) inputfile->Get(name);

  // histogram with x dependency for muon-in-jet-tagged-away-jet-tagged sample
  sprintf(name, "%s/ptag_%s_%s", directory, dependency, tag);
  TH2F * histogram2D_tt = (TH2F*) inputfile->Get(name);
       
  // histogram base name
  sprintf(basename, "counting_%s_%s", dependency, tag);

  // histogram that contains the measure efficiency  
  TH1F * efficiency = (TH1F*) histogram2D_0t->ProjectionX(basename);
  Int_t nbins = efficiency->GetNbinsX();

  // histogram that contains the measure efficiency + systematic
  sprintf(name, "%s_plus", basename);
  TH1F * efficiencyPlus = (TH1F*) histogram2D_0t->ProjectionX(name);

  // histogram that contains the measure efficiency - systematic
  sprintf(name, "%s_minus", basename);
  TH1F * efficiencyMinus = (TH1F*) histogram2D_0t->ProjectionX(name);

  // histogram2D with x dependency for muon-in-jet sample from cl flavor
  sprintf(name, "%s/n_%s_cl", directory, dependency);
  TH2F * histogram2D_cl = (TH2F*) inputfile->Get(name);

  // 1D histogram with x dependency for muon-in-jet sample from cl flavor
  sprintf(name, "%s_n_cl", basename);
  TH1F * n_cl = (TH1F*) histogram2D_cl->ProjectionX(name, -1, -1, "e");

  // histogram with x dependency for muon-in-jet-away-jet-tagged sample
  sprintf(name, "%s/p_%s_b", directory, dependency);
  TH2F * histogram2D_0t_b = (TH2F*) inputfile->Get(name);

  // 1D histogram with x dependency for muon-in-jet-away-jet-tagged sample from b flavor  
  sprintf(name, "%s_p_b", basename);
  TH1F * p_b = (TH1F*) histogram2D_0t_b->ProjectionX(name, -1, -1, "e");

  // efficiency calculation from mc
  sprintf(name, "%s/ptag_%s_b_%s", directory, dependency, tag);  
  TH2F * histogram2D_tt_b = (TH2F*) inputfile->Get(name);  
  sprintf(name, "%s_ptag_b", basename);
  TH1F * ptag_b = (TH1F*) histogram2D_tt_b->ProjectionX(name, -1, -1, "e");
  sprintf(name, "%s_mc", basename);
  TH1F * mc = (TH1F*) histogram2D_tt_b->ProjectionX(name, -1, -1, "e");
  mc->Divide(mc, p_b, 1., 1., "e");

  // mistag rate calculation from mc
  sprintf(name, "%s/p_%s_cl", directory, dependency);
  TH2F * histogram2D_0t_cl = (TH2F*) inputfile->Get(name);
  sprintf(name, "%s_mistag", basename);
  TH1F * mistag = (TH1F*) histogram2D_0t_cl->ProjectionX(name, -1, -1, "e");
  mistag->Divide(mistag, n_cl, 1., 1., "e");

  // Auxiliaries histograms
  TH1F * histogram1D;
  TH1F * histogram1D_0t;
  TH1F * histogram1D_tt;

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
    thePdf1 = getPdfByIndex(pdfbase * ii, "");
    thePdf2 = getPdfByIndex(pdfbase * ii, "tag");

    if (!thePdf1 || !thePdf2)
    {
      std::cout << "can't locate the correct PDF for this data set ..." << std::endl;
      break;
    }

    // 1D histogram with x dependency for muon-in-jet sample
    sprintf(name, "%s_%s_%d", basename, tag, ii);
    histogram1D = (TH1F*) histogram2D->ProjectionY(name, ii, ii);

    // 1D histogram with x dependency for muon-in-jet-away-jet-tagged sample
    sprintf(name, "%s_%s_%d_0t", basename, tag, ii);
    histogram1D_0t = (TH1F*) histogram2D_0t->ProjectionY(name, ii, ii);
    
    // 1D histogram with x dependency for muon-in-jet-tagged-away-jet-tagged sample
    sprintf(name, "%s_%s_%d_tt", basename, tag, ii);
    histogram1D_tt = (TH1F*) histogram2D_tt->ProjectionY(name, ii, ii);

    // get the b fraction in the muon-in-jet
    Fit(histogram1D, thePdf1, &n, &n_err);
         
    // get the b fraction in the muon-in-jet-tagged-away-jet-tagged
    Fit(histogram1D_tt, thePdf2, &nt, &nt_err);
        
    Nst = (Double_t) histogram1D_0t->GetEntries();
    Mtg = mistag->GetBinContent(ii); 
    MtgErr = mistag->GetBinError(ii); 
    
    // efficiencies
    Eff         = nt[0] / ( Nst - n[1] * Mtg );
    EffMtgPlus  = nt[0] / ( Nst - n[1] * Mtg * 1.05 );
    EffMtgMinus = nt[0] / ( Nst - n[1] * Mtg * 0.95 );

    // efficiencies systematic and statistical errors
    EffStatE2 = pow(nt_err[0] / ( Nst - n[1] *  Mtg ), 2)
              + pow(nt[0] * n[1] * MtgErr / pow(Nst - Mtg * n[1], 2), 2)
              + pow(nt[0] * Mtg  * n_err[1] / pow(Nst - Mtg * n[1], 2), 2);
              
    EffSystE2Plus  = pow(Eff - EffMtgPlus, 2);
    EffSystE2Minus = pow(Eff - EffMtgMinus, 2);

    std::cout << "muon-in-jet (# measure c, # truth c)                        : (" << n[1] << ", " << n_cl->GetBinContent(ii) << ")" << std::endl;
    std::cout << "muon-in-jet-awat-jet-tagged (# measure b, # truth b)        : (" << (Nst - n[1] * Mtg) << ", " << p_b->GetBinContent(ii) << ")" << std::endl;
    std::cout << "muon-in-jet-tagged-awat-jet-tagged (# measure b, # truth b) : (" << nt[0] << ", " << ptag_b->GetBinContent(ii) << ")" << std::endl;
    
    std::cout << "Eff             : " << Eff << std::endl;
    std::cout << "EffStatErr      : " << sqrt(EffStatE2) << std::endl;
    std::cout << "EffSystErrPlus  : " << sqrt(EffSystE2Plus) << std::endl;
    std::cout << "EffSystErrMinus : " << sqrt(EffSystE2Minus) << std::endl;    
    std::cout << "Mistag          : " << Mtg << " pm " << MtgErr << std::endl;
    std::cout << std::endl;
    
    efficiency->SetBinContent(ii, Eff);
    efficiency->SetBinError(ii, sqrt(EffStatE2));

    efficiencyPlus->SetBinContent(ii, Eff+sqrt(EffStatE2+EffSystE2Plus));
    efficiencyPlus->SetBinError(ii, 0.0);
    efficiencyMinus->SetBinContent(ii, Eff-sqrt(EffStatE2+EffSystE2Minus));
    efficiencyMinus->SetBinError(ii, 0.0);

    // outfile->cd();
    // histogram1D->Write();
    // histogram1D_0t->Write();
    // histogram1D_tt->Write();
  }
    
  sprintf(name, "%s_mc", basename);  
  mc->SetName(name);
  mc->SetLineColor(kRed);
 
  outfile->cd();
  mc->Write();
  efficiency->Write();
  efficiencyPlus->Write();
  efficiencyMinus->Write();

  // make histograms.

  TString title;
  TLegend *leg = 0;
  const char *epsname= 0;

  if (!strcmp(dependency,"pT"))
    formatHist1(efficiencyPlus, "p_{T} (GeV)", "Efficiency");
  else if (!strcmp(dependency,"eta"))
    formatHist1(efficiencyPlus, "eta", "Efficiency");

  efficiencyPlus->GetYaxis()->SetRangeUser(0.0,1.4);
  efficiencyPlus->SetFillColor(kBlack);
  efficiencyPlus->SetFillStyle(3005);
  efficiencyPlus->SetLineColor(kWhite);
  
  epsname = saveAsEPS(c1, efficiencyPlus, "BAR");
  efficiencyMinus->SetLineColor(kWhite);
  efficiencyMinus->Draw("BAR SAME");
  efficiency->Draw("SAME");
  mc->Draw("SAME E");
    
  mc->SetMarkerSize(0);  
  leg = new TLegend(0.25, 0.8, 0.6, 0.9);
  leg->AddEntry(efficiency, "Counting", "pl");
  leg->AddEntry(mc, "MC truth", "l");
  leg->SetFillColor(0);
  leg->Draw();
  gPad->RedrawAxis();

  c1->cd();
  c1->SaveAs(epsname);

  outfile->Write();
  outfile->Close();
  
}


void PtrelSolver::counting(
  const char *inputfilename,
  const char *directory,
  const char *tag
)
{  
  // input file with data.
  inputfile = new TFile(inputfilename, "READ");
 
  char name[255], basename[255];
    
  // histogram with x dependency for muon-in-jet sample
  sprintf(name, "%s/n_pT", directory);
  TH2F * histogram2D = (TH2F*) inputfile->Get(name);

  // histogram with x dependency for muon-in-jet-away-jet-tagged sample
  sprintf(name, "%s/p_pT", directory);
  TH2F * histogram2D_0t = (TH2F*) inputfile->Get(name);

  // histogram with x dependency for muon-in-jet-tagged-away-jet-tagged sample
  sprintf(name, "%s/ptag_pT_%s", directory, tag);
  TH2F * histogram2D_tt = (TH2F*) inputfile->Get(name);

  // 1D histogram with x dependency for muon-in-jet sample
  sprintf(name, "%s_%s_%d", basename, tag, 0);  
  TH1F * histogram1D = (TH1F*) histogram2D->ProjectionY(name);

  // 1D histogram with x dependency for muon-in-jet-away-jet-tagged sample
  sprintf(name, "%s_%s_%d_0t", basename, tag, 0);
  TH1F * histogram1D_0t = (TH1F*) histogram2D_0t->ProjectionY(name);
    
  // 1D histogram with x dependency for muon-in-jet-tagged-away-jet-tagged sample
  sprintf(name, "%s_%s_%d_tt", basename, tag, 0);
  TH1F * histogram1D_tt = (TH1F*) histogram2D_tt->ProjectionY(name);

  // histogram with x dependency for muon-in-jet-away-jet-tagged sample
  sprintf(name, "%s/p_pT_b", directory);
  TH2F * histogram2D_0t_b = (TH2F*) inputfile->Get(name);

  // 1D histogram with x dependency for muon-in-jet-away-jet-tagged sample from b flavor  
  sprintf(name, "%s_p_b", basename);
  TH1F * p_b = (TH1F*) histogram2D_0t_b->ProjectionX(name, -1, -1, "e");

  // histogram2D with x dependency for muon-in-jet sample from cl flavor
  sprintf(name, "%s/n_pT_cl", directory);
  TH2F * histogram2D_cl = (TH2F*) inputfile->Get(name);

  // 1D histogram with x dependency for muon-in-jet sample from cl flavor
  sprintf(name, "%s_n_cl", basename);
  TH1F * n_cl = (TH1F*) histogram2D_cl->ProjectionX(name, -1, -1, "e");

  // mistag rate calculation from mc
  sprintf(name, "%s/p_pT_cl", directory);
  TH2F * histogram2D_0t_cl = (TH2F*) inputfile->Get(name);
  sprintf(name, "%s_p_cl", basename);
  TH1F * p_cl = (TH1F*) histogram2D_0t_cl->ProjectionX(name, -1, -1, "e");
  
  // b constent in ptag samples 
  sprintf(name, "%s/ptag_pT_b_%s", directory, tag);  
  TH2F * histogram2D_tt_b = (TH2F*) inputfile->Get(name);  
  sprintf(name, "%s_ptag_b", basename);
  TH1F * ptag_b = (TH1F*) histogram2D_tt_b->ProjectionX(name, -1, -1, "e");
   
  std::vector<double>  n, n_err;
  std::vector<double>  nt, nt_err;
  
  Double_t Nst, Mtg, MtgErr;
  Double_t Eff, EffMtgPlus, EffMtgMinus; 
  Double_t EffStatE2, EffSystE2Plus, EffSystE2Minus;

  // get the right template
  TF1 * thePdf1 = getPdfByIndex(0, "");
  TF1 * thePdf2 = getPdfByIndex(0, "tag");

  // get the b fraction in the muon-in-jet
  Fit(histogram1D, thePdf1, &n, &n_err);
    
  std::cout << n[0] << " " << n[1] << std::endl;  
    
  return;  
         
  // get the b fraction in the muon-in-jet-tagged-away-jet-tagged
  Fit(histogram1D_tt, thePdf2, &nt, &nt_err);

  Nst = (Double_t) histogram1D_0t->GetEntries();
  Mtg = p_cl->GetEntries()/n_cl->GetEntries();  	
  MtgErr = 0;
    
  // efficiencies
  Eff         = nt[0] / ( Nst - n[1] * Mtg );
  EffMtgPlus  = nt[0] / ( Nst - n[1] * Mtg * 1.05 );
  EffMtgMinus = nt[0] / ( Nst - n[1] * Mtg * 0.95 );

  // efficiencies systematic and statistical errors
  EffStatE2 = pow(nt_err[0] / ( Nst - n[1] *  Mtg ), 2) 
            + pow(nt[0] * n[1] * MtgErr / pow(Nst - Mtg * n[1], 2), 2)
            + pow(nt[0] * Mtg  * n_err[1] / pow(Nst - Mtg * n[1], 2), 2);
              
  EffSystE2Plus  = pow(Eff - EffMtgPlus, 2);
  EffSystE2Minus = pow(Eff - EffMtgMinus, 2);

  std::cout << std::endl << "===============================================================================" << std::endl;

  std::cout << "muon-in-jet (# measure c, # truth c)                        : (" << n[1] << ", " << n_cl->GetEntries() << ")" << std::endl;
  std::cout << "muon-in-jet-awat-jet-tagged (# measure b, # truth b)        : (" << (Nst - n[1] * Mtg) << ", " << p_b->GetEntries() << ")" << std::endl;
  std::cout << "muon-in-jet-tagged-awat-jet-tagged (# measure b, # truth b) : (" << nt[0] << ", " << ptag_b->GetEntries() << ")" << std::endl;
    
  std::cout << "Total Eff             : " << Eff << std::endl;
  std::cout << "Total EffStatErr      : " << sqrt(EffStatE2) << std::endl;
  std::cout << "Total EffSystErrPlus  : " << sqrt(EffSystE2Plus) << std::endl;
  std::cout << "Total EffSystErrMinus : " << sqrt(EffSystE2Minus) << std::endl;   
  std::cout << "Total EffErrPlus      : " << sqrt(EffStatE2+EffSystE2Plus) << std::endl;   
  std::cout << "Total EffErrMinus     : " << sqrt(EffStatE2+EffSystE2Minus) << std::endl;      
  std::cout << "Total Mistag          : " << Mtg << " mc: " << MtgErr << std::endl;

  std::cout << "===============================================================================" << std::endl << std::endl; 
}  


// Read the right pdf data file for the counting method
bool PtrelSolver::initPdfsByTag(const char * directory, const char * tag, const char * versiontag)
{
  char b_filename[256], cl_filename[256];

  sprintf(b_filename, "%s/nb_templates_%s", directory, versiontag);
  sprintf(cl_filename, "%s/ncl_templates_%s", directory, versiontag);

  if ( !locateFile(b_filename) || !locateFile(cl_filename) )
  {
    std::cout << "information: pdf data file " << b_filename << " & "
	          << cl_filename << " dont exist ... " << std::endl;
    return 0;
  }

  // initialize the first set of templates
  std::cout << "Information: " << b_filename << " & " << cl_filename << std::endl;
  initPdfs(b_filename, cl_filename, "");

  char btag_filename[256], cltag_filename[256];

//  sprintf(btag_filename, "%s/pb_templates_%s%s", directory, tag, versiontag);
//  sprintf(cltag_filename, "%s/pcl_templates_%s%s", directory, tag, versiontag);

  sprintf(btag_filename, "%s/pb_templates_%s%s", directory, tag, versiontag);
  sprintf(cltag_filename, "%s/pcl_templates_%s%s", directory, tag, versiontag);

  std::cout << btag_filename << std::endl;

  if ( !locateFile(btag_filename) || !locateFile(cltag_filename) )
  {
    std::cout << "information: pdf data file " << btag_filename << " & "
	          << cltag_filename << " dont exist ... " << std::endl;
    initPdfs(b_filename, cl_filename, "tag");
    std::cout << "information: initialize tag pdfs with untagged pdfs " << std::endl;
  }
  else
    initPdfs(btag_filename, cltag_filename, "tag");

  return true;
}


void PtrelSolver::measureByCounting(
  const char * inputfilename,
  const char * directory,
  const char * tag,
  const char * outfilename,
  const char * pdfdir,
  const char * versiontag
)
{
  if ( !initPdfsByTag(pdfdir, tag, versiontag) ) return;
  counting(inputfilename, directory, "pT", tag, outfilename); 
  counting(inputfilename, directory, "eta", tag, outfilename);
  counting(inputfilename, directory, tag);  
}




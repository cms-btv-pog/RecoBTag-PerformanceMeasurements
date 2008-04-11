
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

ClassImp(PtrelSolver)





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



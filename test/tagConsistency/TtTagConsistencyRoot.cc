/* -*- mode: c++ -*- */

#include "RecoBTag/PerformanceMeasurements/test/tagConsistency/TtTagConsistencyRoot.h"
#include "RecoBTag/PerformanceMeasurements/interface/RooGKCounter.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDCacheFile.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "AnalysisDataFormats/TopObjects/interface/TopLepton.h"
#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMET.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"




int TtTagConsistencyRoot::scan_plot( string data_list_file, string prefix,
				     double discr_min, double discr_max, double discr_step,
				     double mistag_min, double mistag_max )
{

  //gROOT->Reset();

  setTDRStyle();

  double xsec_tt0j  =   190.0;
  double xsec_tt1j  =   170.0;
  double xsec_tt2j  =   100.0;
  double xsec_tt3j  =    40.0;
  double xsec_tt4j  =    61.0;
  double xsec_ttbar = xsec_tt0j + xsec_tt1j + xsec_tt2j + xsec_tt3j + xsec_tt4j; 
  double xsec_w0j   = 30000.0;
  double xsec_w1j   =  8000.0;
  double xsec_w2j   =  2500.0;
  double xsec_w3j   =   722.0;
  double xsec_w4j   =   174.0;
  double xsec_w5j   =    45.0;

  double nEvents_tt0j = 418311.0;
  double nEvents_tt1j = 377989.0;
  double nEvents_tt2j = 211787.0;
  double nEvents_tt3j =  45569.0;
  double nEvents_tt4j =  20861.0;
  double nEvents_tt   =  nEvents_tt0j + nEvents_tt1j + nEvents_tt2j + nEvents_tt3j + nEvents_tt4j;
  double nEvents_ttbar = 56000.0;
  double nEvents_w0j = 174933.0;
  double nEvents_w1j = 169000.0;
  double nEvents_w2j = 180937.0;
  double nEvents_w3j =  68199.0;
  double nEvents_w4j =  72585.0;
  double nEvents_w5j =  55209.0;

  double nSelected_tt0j = 70779.0;
  double nSelected_tt1j = 95908.0;
  double nSelected_tt2j = 68102.0;
  double nSelected_tt3j = 16738.0;
  double nSelected_tt4j =  8621.0;
  double nSelected_ttbar = 1687;
  double nSelected_w0j =    44.0;
  double nSelected_w1j =   259.0;
  double nSelected_w2j =  1471.0;
  double nSelected_w3j =  3342.0;
  double nSelected_w4j = 11898.0;
  double nSelected_w5j = 12721.0;

  double discr[50];
  double discr_err[50];

  //double discr_min  =  2.0;
  //double discr_max  =  10.0;
  //double discr_step =  0.5;
  //
  //double discr_min  =  0.0;
  //double discr_max  =  1.0;
  //double discr_step =  0.05;

  int n = (int)( (discr_max - discr_min) / discr_step + 1.01 );

  TtTagConsistencyEfficiency b;

  double tt0j_weight = nEvents_tt/nEvents_tt0j * xsec_tt0j/xsec_ttbar ;
  double tt1j_weight = nEvents_tt/nEvents_tt1j * xsec_tt1j/xsec_ttbar ;
  double tt2j_weight = nEvents_tt/nEvents_tt2j * xsec_tt2j/xsec_ttbar ;
  double tt3j_weight = nEvents_tt/nEvents_tt3j * xsec_tt3j/xsec_ttbar ;
  double tt4j_weight = nEvents_tt/nEvents_tt4j * xsec_tt4j/xsec_ttbar ;

  int nFijk_sig = b . updateFijk2( "tagConsistency/tt0j.legacy.mc2.list", tt0j_weight, "signal" );
  nFijk_sig = b . updateFijk2( "tagConsistency/tt1j.legacy.mc2.list", tt1j_weight, "signal" );
  nFijk_sig = b . updateFijk2( "tagConsistency/tt2j.legacy.mc2.list", tt2j_weight, "signal" );
  nFijk_sig = b . updateFijk2( "tagConsistency/tt3j.legacy.mc2.list", tt3j_weight, "signal" );
  nFijk_sig = b . updateFijk2( "tagConsistency/tt4j.legacy.mc2.list", tt4j_weight, "signal" );

  /*
  double w0j_weight = nEvents_ttbar/nEvents_w0j * xsec_w0j/xsec_ttbar;
  double w1j_weight = nEvents_ttbar/nEvents_w1j * xsec_w1j/xsec_ttbar;
  double w2j_weight = nEvents_ttbar/nEvents_w2j * xsec_w2j/xsec_ttbar;
  double w3j_weight = nEvents_ttbar/nEvents_w3j * xsec_w3j/xsec_ttbar;
  double w4j_weight = nEvents_ttbar/nEvents_w4j * xsec_w4j/xsec_ttbar;
  double w5j_weight = nEvents_ttbar/nEvents_w5j * xsec_w5j/xsec_ttbar;
  */

  double w0j_weight = nEvents_w0j/nEvents_w0j * xsec_w0j/xsec_ttbar;
  double w1j_weight = nEvents_w0j/nEvents_w1j * xsec_w1j/xsec_ttbar;
  double w2j_weight = nEvents_w0j/nEvents_w2j * xsec_w2j/xsec_ttbar;
  double w3j_weight = nEvents_w0j/nEvents_w3j * xsec_w3j/xsec_ttbar;
  double w4j_weight = nEvents_w0j/nEvents_w4j * xsec_w4j/xsec_ttbar;
  double w5j_weight = nEvents_w0j/nEvents_w5j * xsec_w5j/xsec_ttbar;

  int nFijk_bg = b . updateFijk2( "tagConsistency/w0j.legacy.mc2.list", w0j_weight, "background" );
  nFijk_bg = b . updateFijk2( "tagConsistency/w1j.legacy.mc2.list", w1j_weight, "background" );
  nFijk_bg = b . updateFijk2( "tagConsistency/w2j.legacy.mc2.list", w2j_weight, "background" );
  nFijk_bg = b . updateFijk2( "tagConsistency/w3j.legacy.mc2.list", w3j_weight, "background" );
  nFijk_bg = b . updateFijk2( "tagConsistency/w4j.legacy.mc2.list", w4j_weight, "background" );
  nFijk_bg = b . updateFijk2( "tagConsistency/w5j.legacy.mc2.list", w5j_weight, "background" );

  //TtTagConsistencyEfficiency::tableContent tab_cont;
  TtTagConsistencyEfficiency::tableLine tab;

  b . cleanTable( tab );

  //b . updateTable( "./tagConsistency/legacy/fakeDataSet.tab", 4.0, 1.0, tab );
  //b . updateTableFromList( "./tagConsistency/trackCountingHighEffJetTags_csa07_chowder_topSemiLep_100pb.tab.list", 4.0, 1.0, tab );
  b . updateTableFromList( data_list_file, 0.5, 1.0, tab );

  //---- MC truth
  double goodB = tab . Nb;
  double goodC = tab . Nc;
  double goodL = tab . Nl;

  double taggedL[50];
  double taggedB[50];
  double taggedC[50];

  // ---> composing pure signal
  int n1[50];
  int n2[50];
  int n3[50];
  for ( int j = 0; j < n; j++ )
    {
      discr[j] = discr_min + j*discr_step;

      b . cleanTable( tab );

      //b . updateTable( "./tagConsistency/legacy/fakeDataSet.tab", discr[j], 1.0, tab );
      //b . updateTableFromList( "./tagConsistency/trackCountingHighEffJetTags_csa07_chowder_topSemiLep_100pb.tab.list", discr[j], 1.0, tab );
      b . updateTableFromList( data_list_file, discr[j], 1.0, tab );

      n1[j] = tab . N[1];
      n2[j] = tab . N[2];
      n3[j] = tab . N[3];

      taggedB[j] = tab . Nbtag;
      taggedC[j] = tab . Nctag;
      taggedL[j] = tab . Nltag;

      cout << n1[j] << "   " << n2[j] << "   " << n3[j] << endl;
      
    }

  double bgf = 0.1303; // fraction of BG in the fit
  //double bgf = ( w0j_weight * nSelected_w0j + w1j_weight * nSelected_w1j + w2j_weight * nSelected_w2j + w3j_weight * nSelected_w3j + w4j_weight * nSelected_w4j + w5j_weight * nSelected_w5j ) / nSelected_ttbar * nEvents_ttbar / nEvents_w0j;

  double epsb[50];
  double epsb_err[50];
  double epsc[50];
  double epsc_err[50];
  double epsl[50];
  double epsl_err[50];
  double mistag[50];
  double mistag_err[50];

  double epsb_true[50];
  double epsb_true_err[50];
  double epsc_true[50];
  double epsc_true_err[50];
  double epsl_true[50];
  double epsl_true_err[50];
  double mistag_true[50];
  double mistag_true_err[50];

  for (int i = 0; i<n; i++)
    {
      discr_err[i] = 0.0;

      epsl_true[i] = taggedL[i]/goodL;
      epsl_true_err[i] = epsl_true[i] / sqrt( taggedL[i] );
      epsl[i] = epsl_true[i];
      epsl_err[i] = epsl_true_err[i];

      epsb_true[i] = taggedB[i]/goodB;
      epsb_true_err[i] = epsb_true[i] / sqrt( taggedB[i] );

      epsc_true[i] = taggedC[i]/goodC;
      epsc_true_err[i] = epsc_true[i] / sqrt( taggedC[i] );


      mistag_true[i] = (taggedL[i]+taggedC[i])/(goodL+goodC+goodB);
      mistag_true_err[i] = mistag_true[i] / sqrt(taggedL[i]+taggedC[i]);

      //===> nominal fit, signal + bg
      //RooFitResult r = b . minTest2( 0, n1[i], n2[i], n3[i], 0, "msh", epsb_true[i], epsc_true[i], epsl[i], 560.0, bgf );

      TtTagConsistencyFitConfig _conf;
      _conf . n1 = n1[i];
      _conf . n2 = n2[i];
      _conf . n3 = n3[i];
      //_conf . bgf1 = 0.1303;
      //_conf . bgf2 = 0.1303;
      //_conf . bgf3 = 0.1303;
      //_conf . eb = epsb_true[i];
      //_conf . ec = epsc_true[i];
      //_conf . el = epsl[i];

      _conf . setConfig( n1[i], n2[i], n3[i], epsb_true[i], epsc_true[i], epsl[i] );

      //cout << "DEBUG: n1,2,3 = " << n1[i] << "   " << n2[i] << "   " << n3[i] << endl; 
      cout << "DEBUG: n1,2,3 = " << _conf.n1 << "   " << _conf.n2 << "   " << _conf.n3 << endl; 
      cout << "DEBUG: eb,ec,el = " << _conf.eb << "   " << _conf.ec << "   " << _conf.el << endl; 
      cout << "DEBUG: eb_max,ec_max,el_max = " << _conf.eb_max << "   " << _conf.ec_max << "   " << _conf.el_max << endl; 
      cout << "DEBUG: eb_min,ec_min,el_min = " << _conf.eb_min << "   " << _conf.ec_min << "   " << _conf.el_min << endl; 
      RooFitResult r = b . fit( _conf );

      RooRealVar * eb = dynamic_cast<RooRealVar*>(r.floatParsFinal().find("epsb"));
      //RooRealVar * ec = dynamic_cast<RooRealVar*>r.floatParsFinal().find("epsc");

      epsb[i] = eb -> getVal();
      epsb_err[i] = eb -> getError();
      //epsc[i] = ec -> getVal();
      //epsc_err[i] = ec -> getError();
    }

  //gStyle -> SetOptStat(0);

  TCanvas c1("c1","Tagging efficiencies",200,10,700,595);

  TH2F h( "h", "", 100, discr_min, discr_max, 100, 0.0, 1.0 );
  TH2F h2( "h2", "", 100, 0.0, 1.0, 100, mistag_min, mistag_max ); // eff vs mistag
  
  h . Draw();
  
  //c1->SetFillColor(42);
  c1 . SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  //c1->GetFrame()->SetBorderSize(12);

  TGraphErrors gr_epsb( n, discr, epsb, discr_err, epsb_err );
  //TGraphErrors gr_epsc( n, discr, epsc, discr_err, epsc_err );
  TGraphErrors gr_epsb_mistag_true( n, epsb, mistag_true, epsb_err, mistag_true_err );

  TGraphErrors gr_epsl( n, discr, epsl, discr_err, epsl_err );

  TGraphErrors gr_epsb_true( n, discr, epsb_true, discr_err, epsb_true_err );
  TGraphErrors gr_epsb_true_mistag_true( n, epsb_true, mistag_true, epsb_true_err, mistag_true_err );
  //TGraphErrors gr_epsc_true( n, discr, epsc_true, discr_err, epsc_true_err );

  gr_epsb . Draw("LP");
  //gr_epsc . Draw("LP");
  //gr_epsl . Draw("LP");

  gr_epsb_true . Draw("LP");
  //gr_epsc_true . Draw("LP");


  gr_epsb . SetLineWidth(2);
  gr_epsb . SetLineColor(4);
  gr_epsb . SetMarkerColor(4);
  gr_epsb . SetMarkerStyle(kOpenStar);
  gr_epsb . SetMarkerSize(1.5);

  //gr_epsc . SetLineWidth(2);
  //gr_epsc . SetLineColor(2);
  //gr_epsc . SetMarkerColor(2);

  gr_epsl . SetLineWidth(2);
  gr_epsl . SetLineColor(8);
  gr_epsl . SetMarkerColor(8);

  gr_epsb_true . SetLineWidth(2);
  gr_epsb_true . SetLineColor(1);
  gr_epsb_true . SetLineStyle(7);
  gr_epsb_true . SetMarkerColor(1);
  //gr_epsb_true . SetMarkerStyle(1);
  gr_epsb_true . SetMarkerStyle(kFullTriangleUp);
  gr_epsb_true . SetMarkerSize(1.5);

  /*
  gr_epsc_true . SetLineWidth(2);
  gr_epsc_true . SetLineColor(1);
  gr_epsc_true . SetLineStyle(7);
  gr_epsc_true . SetMarkerColor(1);
  gr_epsc_true . SetMarkerStyle(1);
  */

  c1 . Modified();

  Double_t xMin = discr_min;
  Double_t xMax = discr_max;
  //Double_t yMin = h -> GetYaxis() -> GetXmin();
  //Double_t yMax = h -> GetYaxis() -> GetXmax();
  Double_t yMin = 0;
  Double_t yMax = 1;

  TPaveText legend( xMin+0.4*(xMax-xMin),0.8*(yMax-yMin),xMin+0.95*(xMax-xMin),0.95*(yMax-yMin));
  legend . AddText("b-tagging efficiency");
  legend . AddText("preliminary");
  legend . SetBorderSize(0);
  legend . SetFillColor(0);
  legend . Draw();

  string buf = prefix + "_discr.eps";
  c1 . SaveAs( buf.c_str() );
  //c1 . SaveAs( "bTagEff.eps" );
  //c1 . SaveAs( "cTagEff.eps" );
  
  //c1 . SaveAs( "bcContour.eps" );

  //===> Draw eff vs. mistag
  h2.GetXaxis()->SetNdivisions(5);
  //h2.GetYaxis()->SetMoreLogLabels();
  h2.Draw();
  gr_epsb_mistag_true.Draw("P");
  gr_epsb_true_mistag_true.Draw("P");
  gr_epsb_mistag_true . SetLineWidth(2);
  gr_epsb_mistag_true . SetLineColor(4);
  gr_epsb_mistag_true . SetMarkerColor(4);
  gr_epsb_mistag_true . SetMarkerStyle(kOpenStar);
  gr_epsb_mistag_true . SetMarkerSize(1.5);
  gr_epsb_true_mistag_true . SetLineWidth(2);
  gr_epsb_true_mistag_true . SetLineColor(1);
  gr_epsb_true_mistag_true . SetLineStyle(7);
  gr_epsb_true_mistag_true . SetMarkerColor(1);
  gr_epsb_true_mistag_true . SetMarkerStyle(kFullTriangleUp);
  gr_epsb_true_mistag_true . SetMarkerSize(1.5);
  c1.SetLogy();
  c1 . Modified();
  buf = prefix + "_mistag.eps";
  c1 . SaveAs( buf.c_str() );

  return 0;
}






int TtTagConsistencyRoot::scan_plot_d( string data_list_file, string prefix,
				       double discr_min, double discr_max, double discr_step,
				       double mistag_min, double mistag_max,
				       double frac1, double frac2, double frac3,
				       bool use_qcd )
{

  //gROOT->Reset();

  setTDRStyle();

  double xsec_tt0j  =   190.0;
  double xsec_tt1j  =   170.0;
  double xsec_tt2j  =   100.0;
  double xsec_tt3j  =    40.0;
  double xsec_tt4j  =    61.0;
  double xsec_ttbar = xsec_tt0j + xsec_tt1j + xsec_tt2j + xsec_tt3j + xsec_tt4j; 
  double xsec_w0j   = 30000.0;
  double xsec_w1j   =  8000.0;
  double xsec_w2j   =  2500.0;
  double xsec_w3j   =   722.0;
  double xsec_w4j   =   174.0;
  double xsec_w5j   =    45.0;

  double nEvents_tt0j = 418311.0;
  double nEvents_tt1j = 377989.0;
  double nEvents_tt2j = 211787.0;
  double nEvents_tt3j =  45569.0;
  double nEvents_tt4j =  20861.0;
  double nEvents_tt   =  nEvents_tt0j + nEvents_tt1j + nEvents_tt2j + nEvents_tt3j + nEvents_tt4j;
  double nEvents_ttbar = 56000.0;
  double nEvents_w0j = 174933.0;
  double nEvents_w1j = 169000.0;
  double nEvents_w2j = 180937.0;
  double nEvents_w3j =  68199.0;
  double nEvents_w4j =  72585.0;
  double nEvents_w5j =  55209.0;

  double nSelected_tt0j = 70779.0;
  double nSelected_tt1j = 95908.0;
  double nSelected_tt2j = 68102.0;
  double nSelected_tt3j = 16738.0;
  double nSelected_tt4j =  8621.0;
  double nSelected_ttbar = 1687;
  double nSelected_w0j =    44.0;
  double nSelected_w1j =   259.0;
  double nSelected_w2j =  1471.0;
  double nSelected_w3j =  3342.0;
  double nSelected_w4j = 11898.0;
  double nSelected_w5j = 12721.0;

  double discr[50];
  double discr_err[50];

  //double discr_min  =  2.0;
  //double discr_max  =  10.0;
  //double discr_step =  0.5;
  //
  //double discr_min  =  0.0;
  //double discr_max  =  1.0;
  //double discr_step =  0.05;

  int n = (int)( (discr_max - discr_min) / discr_step + 1.01 );

  TtTagConsistencyEfficiency b;

  double tt0j_weight = nEvents_tt/nEvents_tt0j * xsec_tt0j/xsec_ttbar ;
  double tt1j_weight = nEvents_tt/nEvents_tt1j * xsec_tt1j/xsec_ttbar ;
  double tt2j_weight = nEvents_tt/nEvents_tt2j * xsec_tt2j/xsec_ttbar ;
  double tt3j_weight = nEvents_tt/nEvents_tt3j * xsec_tt3j/xsec_ttbar ;
  double tt4j_weight = nEvents_tt/nEvents_tt4j * xsec_tt4j/xsec_ttbar ;

  double nFijk_sig = b . updateFijk_d( "tagConsistency/tt0j.legacy.mc2.list", tt0j_weight, "signal" );
  nFijk_sig = b . updateFijk_d( "tagConsistency/tt1j.legacy.mc2.list", tt1j_weight, "signal" );
  nFijk_sig = b . updateFijk_d( "tagConsistency/tt2j.legacy.mc2.list", tt2j_weight, "signal" );
  nFijk_sig = b . updateFijk_d( "tagConsistency/tt3j.legacy.mc2.list", tt3j_weight, "signal" );
  nFijk_sig = b . updateFijk_d( "tagConsistency/tt4j.legacy.mc2.list", tt4j_weight, "signal" );

  /*
  double w0j_weight = nEvents_ttbar/nEvents_w0j * xsec_w0j/xsec_ttbar;
  double w1j_weight = nEvents_ttbar/nEvents_w1j * xsec_w1j/xsec_ttbar;
  double w2j_weight = nEvents_ttbar/nEvents_w2j * xsec_w2j/xsec_ttbar;
  double w3j_weight = nEvents_ttbar/nEvents_w3j * xsec_w3j/xsec_ttbar;
  double w4j_weight = nEvents_ttbar/nEvents_w4j * xsec_w4j/xsec_ttbar;
  double w5j_weight = nEvents_ttbar/nEvents_w5j * xsec_w5j/xsec_ttbar;
  */

  double w0j_weight = nEvents_w0j/nEvents_w0j * xsec_w0j/xsec_ttbar;
  double w1j_weight = nEvents_w0j/nEvents_w1j * xsec_w1j/xsec_ttbar;
  double w2j_weight = nEvents_w0j/nEvents_w2j * xsec_w2j/xsec_ttbar;
  double w3j_weight = nEvents_w0j/nEvents_w3j * xsec_w3j/xsec_ttbar;
  double w4j_weight = nEvents_w0j/nEvents_w4j * xsec_w4j/xsec_ttbar;
  double w5j_weight = nEvents_w0j/nEvents_w5j * xsec_w5j/xsec_ttbar;

  double nFijk_bg = b . updateFijk_d( "tagConsistency/w0j.legacy.mc2.list", w0j_weight, "background" );
  nFijk_bg = b . updateFijk_d( "tagConsistency/w1j.legacy.mc2.list", w1j_weight, "background" );
  nFijk_bg = b . updateFijk_d( "tagConsistency/w2j.legacy.mc2.list", w2j_weight, "background" );
  nFijk_bg = b . updateFijk_d( "tagConsistency/w3j.legacy.mc2.list", w3j_weight, "background" );
  nFijk_bg = b . updateFijk_d( "tagConsistency/w4j.legacy.mc2.list", w4j_weight, "background" );
  nFijk_bg = b . updateFijk_d( "tagConsistency/w5j.legacy.mc2.list", w5j_weight, "background" );

  if (use_qcd){
    cout << "training for QCD..." << endl;
    double qcd_weight = 1.0;
    nFijk_bg = b . updateFijk_d( "list_gumbo_mc2.list", qcd_weight, "background" );
  }

  TtTagConsistencyEfficiency::tableLine tab;

  tab.N[1]+=0.5;
  cout << "DEBUG1: tab.N[1],N[2],N[3] = " << tab.N[1] << "   " << tab.N[2] << "   " << tab.N[3] << endl; 

  b . cleanTable( tab );

  tab.N[1]+=0.5;
  cout << "DEBUG1: tab.N[1],N[2],N[3] = " << tab.N[1] << "   " << tab.N[2] << "   " << tab.N[3] << endl; 

  //b . updateTable( "./tagConsistency/legacy/fakeDataSet.tab", 4.0, 1.0, tab );
  //b . updateTableFromList( "./tagConsistency/trackCountingHighEffJetTags_csa07_chowder_topSemiLep_100pb.tab.list", 4.0, 1.0, tab );
  b . updateTableFromList_d( data_list_file, 0.5, 1.0, tab );

  //---- MC truth
  double goodB = tab . Nb;
  double goodC = tab . Nc;
  double goodL = tab . Nl;

  double taggedL[50];
  double taggedB[50];
  double taggedC[50];

  // ---> composing pure signal
  double n1[50];
  double n2[50];
  double n3[50];
  cout << "DEBUG1: n1,2,3 = " << n1[0] << "   " << n2[0] << "   " << n3[0] << endl; 
  for ( int j = 0; j < n; j++ )
    {
      discr[j] = discr_min + j*discr_step;

      b . cleanTable( tab );

      //b . updateTable( "./tagConsistency/legacy/fakeDataSet.tab", discr[j], 1.0, tab );
      //b . updateTableFromList( "./tagConsistency/trackCountingHighEffJetTags_csa07_chowder_topSemiLep_100pb.tab.list", discr[j], 1.0, tab );
      b . updateTableFromList_d( data_list_file, discr[j], 1.0, tab );

  cout << "DEBUG2: n1,2,3 = " << n1[0] << "   " << n2[0] << "   " << n3[0] << endl; 
      n1[j] = tab . N[1];
      n2[j] = tab . N[2];
      n3[j] = tab . N[3];
  cout << "DEBUG3: n1,2,3 = " << n1[0] << "   " << n2[0] << "   " << n3[0] << endl; 

      taggedB[j] = tab . Nbtag;
      taggedC[j] = tab . Nctag;
      taggedL[j] = tab . Nltag;

      cout << n1[j] << "   " << n2[j] << "   " << n3[j] << endl;
      
    }

  double bgf = 0.1303; // fraction of BG in the fit
  //double bgf = ( w0j_weight * nSelected_w0j + w1j_weight * nSelected_w1j + w2j_weight * nSelected_w2j + w3j_weight * nSelected_w3j + w4j_weight * nSelected_w4j + w5j_weight * nSelected_w5j ) / nSelected_ttbar * nEvents_ttbar / nEvents_w0j;

  double epsb[50];
  double epsb_err[50];
  double epsc[50];
  double epsc_err[50];
  double epsl[50];
  double epsl_err[50];
  double mistag[50];
  double mistag_err[50];

  double epsb_true[50];
  double epsb_true_err[50];
  double epsc_true[50];
  double epsc_true_err[50];
  double epsl_true[50];
  double epsl_true_err[50];
  double mistag_true[50];
  double mistag_true_err[50];

  bool first_fit = true;
  for (int i = 0; i<n; i++)
    {
      discr_err[i] = 0.0;

      epsl_true[i] = taggedL[i]/goodL;
      epsl_true_err[i] = epsl_true[i] / sqrt( taggedL[i] );
      epsl[i] = epsl_true[i];
      epsl_err[i] = epsl_true_err[i];

      epsb_true[i] = taggedB[i]/goodB;
      epsb_true_err[i] = epsb_true[i] / sqrt( taggedB[i] );

      epsc_true[i] = taggedC[i]/goodC;
      epsc_true_err[i] = epsc_true[i] / sqrt( taggedC[i] );


      mistag_true[i] = (taggedL[i]+taggedC[i])/(goodL+goodC+goodB);
      mistag_true_err[i] = mistag_true[i] / sqrt(taggedL[i]+taggedC[i]);

      TtTagConsistencyFitConfig _conf;
      _conf . n1_d = n1[i];
      _conf . n2_d = n2[i];
      _conf . n3_d = n3[i];
      //_conf . bgf1 = 0.217;
      //_conf . bgf1 = 0.3;
      //_conf . bgf2 = 0.024;
      //_conf . bgf3 = 0.006;
      _conf . bgf1 = frac1;
      _conf . bgf2 = frac2;
      _conf . bgf3 = frac3;
      _conf . eb = epsb_true[i];
      _conf . ec = epsc_true[i];
      _conf . el = epsl[i];

      //_conf . eb_min = epsb_true[i]-0.2;

      //cout << "DEBUG: n1,2,3 = " << n1[0] << "   " << n2[0] << "   " << n3[0] << endl; 
      //_conf . setConfig( n1[i], n2[i], n3[i], epsb_true[i], epsc_true[i], epsl[i] );

      //cout << "DEBUG: n1,2,3 = " << n1[0] << "   " << n2[0] << "   " << n3[0] << endl; 
      //cout << "DEBUG: n1,2,3 = " << n1[i] << "   " << n2[i] << "   " << n3[i] << endl; 
      //cout << "DEBUG: n1,2,3 = " << _conf.n1_d << "   " << _conf.n2_d << "   " << _conf.n3_d << endl; 
      //cout << "DEBUG: eb,ec,el = " << _conf.eb << "   " << _conf.ec << "   " << _conf.el << endl; 
      //cout << "DEBUG: eb_max,ec_max,el_max = " << _conf.eb_max << "   " << _conf.ec_max << "   " << _conf.el_max << endl; 
      //cout << "DEBUG: eb_min,ec_min,el_min = " << _conf.eb_min << "   " << _conf.ec_min << "   " << _conf.el_min << endl; 
      RooFitResult r = b . fit_d( _conf );
      //RooFitResult r = b . buildVar( _conf );

      // FIXME: development
      //if (first_fit){
      //	b . buildVar( _conf );
      //	first_fit = false;
      //}
      //RooFitResult r = b . fit_test( _conf );
      

      RooRealVar * eb = dynamic_cast<RooRealVar*>(r.floatParsFinal().find("epsb"));
      //RooRealVar * ec = dynamic_cast<RooRealVar*>r.floatParsFinal().find("epsc");

      epsb[i] = eb -> getVal();
      epsb_err[i] = eb -> getError();
      //epsc[i] = ec -> getVal();
      //epsc_err[i] = ec -> getError();
    }

  //gStyle -> SetOptStat(0);

  TCanvas c1("c1","Tagging efficiencies",700,595);

  TH2F h( "h", "", 100, discr_min, discr_max, 100, 0.0, 1.0 );
  TH2F h2( "h2", "", 100, 0.0, 1.0, 100, mistag_min, mistag_max ); // eff vs mistag
  
  h . Draw();
  
  //c1->SetFillColor(42);
  c1 . SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  //c1->GetFrame()->SetBorderSize(12);

  TGraphErrors gr_epsb( n, discr, epsb, discr_err, epsb_err );
  TGraphErrors gr_epsc( n, discr, epsc, discr_err, epsc_err );
  TGraphErrors gr_epsb_mistag_true( n, epsb, mistag_true, epsb_err, mistag_true_err );

  TGraphErrors gr_epsl( n, discr, epsl, discr_err, epsl_err );

  TGraphErrors gr_epsb_true( n, discr, epsb_true, discr_err, epsb_true_err );
  TGraphErrors gr_epsb_true_mistag_true( n, epsb_true, mistag_true, epsb_true_err, mistag_true_err );
  TGraphErrors gr_epsc_true( n, discr, epsc_true, discr_err, epsc_true_err );

  gr_epsb . Draw("LP");
  gr_epsc . Draw("LP");
  //gr_epsl . Draw("LP");

  gr_epsb_true . Draw("LP");
  gr_epsc_true . Draw("LP");


  gr_epsb . SetLineWidth(2);
  gr_epsb . SetLineColor(4);
  gr_epsb . SetMarkerColor(4);
  gr_epsb . SetMarkerStyle(kOpenStar);
  gr_epsb . SetMarkerSize(1.5);

  gr_epsc . SetLineWidth(2);
  gr_epsc . SetLineColor(2);
  gr_epsc . SetMarkerColor(2);

  gr_epsl . SetLineWidth(2);
  gr_epsl . SetLineColor(8);
  gr_epsl . SetMarkerColor(8);

  gr_epsb_true . SetLineWidth(2);
  gr_epsb_true . SetLineColor(1);
  gr_epsb_true . SetLineStyle(7);
  gr_epsb_true . SetMarkerColor(1);
  //gr_epsb_true . SetMarkerStyle(1);
  gr_epsb_true . SetMarkerStyle(kFullTriangleUp);
  gr_epsb_true . SetMarkerSize(1.5);

  gr_epsc_true . SetLineWidth(2);
  gr_epsc_true . SetLineColor(1);
  gr_epsc_true . SetLineStyle(7);
  gr_epsc_true . SetMarkerColor(1);
  gr_epsc_true . SetMarkerStyle(1);

  c1 . Modified();

  Double_t xMin = discr_min;
  Double_t xMax = discr_max;
  //Double_t yMin = h -> GetYaxis() -> GetXmin();
  //Double_t yMax = h -> GetYaxis() -> GetXmax();
  Double_t yMin = 0;
  Double_t yMax = 1;

  TPaveText legend( xMin+0.4*(xMax-xMin),0.8*(yMax-yMin),xMin+0.95*(xMax-xMin),0.95*(yMax-yMin));
  legend . AddText("b-tagging efficiency");
  legend . AddText("preliminary");
  legend . SetBorderSize(0);
  legend . SetFillColor(0);
  legend . Draw();

  string buf = prefix + "_discr.eps";
  c1 . SaveAs( buf.c_str() );
  //c1 . SaveAs( "bTagEff.eps" );
  //c1 . SaveAs( "cTagEff.eps" );
  
  //c1 . SaveAs( "bcContour.eps" );

  //===> Draw eff vs. mistag
  h2.GetXaxis()->SetNdivisions(5);
  //h2.GetYaxis()->SetMoreLogLabels();
  h2.Draw();
  gr_epsb_mistag_true.Draw("P");
  gr_epsb_true_mistag_true.Draw("P");
  gr_epsb_mistag_true . SetLineWidth(2);
  gr_epsb_mistag_true . SetLineColor(4);
  gr_epsb_mistag_true . SetMarkerColor(4);
  gr_epsb_mistag_true . SetMarkerStyle(kOpenStar);
  gr_epsb_mistag_true . SetMarkerSize(1.5);
  gr_epsb_true_mistag_true . SetLineWidth(2);
  gr_epsb_true_mistag_true . SetLineColor(1);
  gr_epsb_true_mistag_true . SetLineStyle(7);
  gr_epsb_true_mistag_true . SetMarkerColor(1);
  gr_epsb_true_mistag_true . SetMarkerStyle(kFullTriangleUp);
  gr_epsb_true_mistag_true . SetMarkerSize(1.5);
  c1.SetLogy();
  c1 . Modified();
  buf = prefix + "_mistag.eps";
  c1 . SaveAs( buf.c_str() );

  return 0;
}



void TtTagConsistencyRoot::setTDRStyle()
{
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

// For the Global title:

//  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}





void TtTagConsistencyRoot::fwlite_test( void )
{
  //cout << endl << "Beginning FWLite test..." << endl << endl;

  //AutoLibraryLoader::enable();
  
  //TFile file("./tagConsistency/test.root");
  //TFile file("/uscmst1b_scratch/lpc1/3DayLifetime/kukarzev/TQAFLayer1Output_50.root");
  //TDCacheFile * file = new TDCacheFile("/pnfs/cms/WAX/resilient/kukarzev/btag/csa07/chowder/topSemiLepMuon-100pb_2/TQAFLayer1Output_50.root");
  //TFile * file = TFile::Open("dcache:/pnfs/cms/WAX/resilient/kukarzev/btag/csa07/chowder/topSemiLepMuon-100pb_2/TQAFLayer1Output_50.root");
  TFile * file = TFile::Open("/uscmst1b_scratch/lpc1/3DayLifetime/kukarzev/TQAFLayer1Output_50.root");
  
  //fwlite::Event ev(&file);
  //fwlite::Event ev(file);

  /*
  cout << "starting event cycle..." << endl;
  RooGKCounter c(1,100);
  c.setPrintCount(true);
  for( ev.toBegin(); ! ev.atEnd(); ++ev) {
    //    fwlite::Handle<std::vector<TopJet> > caloJets;
    //    caloJets.getByLabel(ev,"selectedLayer1TopJets");

    c.count();

  fwlite::Handle< vector< TopJet > > caloJets;
  fwlite::Handle< vector< TopElectron > > electrons;
  fwlite::Handle< vector< TopMuon > > muons;
  //fwlite::Handle< vector< TopMET > > METs;
  fwlite::Handle< double> weightHandle;

  caloJets  . getByLabel( ev, "selectedLayer1TopJets" );
  electrons . getByLabel( ev, "selectedLayer1TopElectrons" );
  muons     . getByLabel( ev, "selectedLayer1TopMuons" );
  //METs      . getByLabel( ev, "selectedLayer1TopMETs" );


    //now can access data
  //std::cout <<" size "<<caloJets.ptr()->size()<<std::endl;
  int a = caloJets.ptr()->size();
  }
  cout << "done" <<std::endl;  
  */

  //delete file;
}


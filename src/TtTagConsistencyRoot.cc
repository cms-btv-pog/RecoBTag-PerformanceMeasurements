/* -*- mode: c++ -*- */

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyRoot.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"



int TtTagConsistencyRoot::scan_plot( void )
{

  //gROOT->Reset();

  //setTDRStyle();

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

  double discr_min  =  2.0;
  double discr_max  =  10.0;
  double discr_step =  0.5;

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
  //b . cleanTable( tab_cont );

  //b . updateTable( "./tagConsistency/csa07/Chowder-100pb-topSemiLepElectron/tag_consistency_trackCountingHighEffJetTags_83.tab", 4.0, 1.0, tab );
  //b . updateTable( "./tagConsistency/legacy/fakeDataSet.tab", 4.0, 1.0, tab );
  b . updateTableFromList( "./tagConsistency/trackCountingHighEffJetTags_csa07_chowder_topSemiLepElectron_100pb.tab.list", 4.0, 1.0, tab );

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
      //b . cleanTable( tab_cont );

      //b . updateTable( "./tagConsistency/legacy/fakeDataSet.tab", discr[j], 1.0, tab );
      b . updateTableFromList( "./tagConsistency/trackCountingHighEffJetTags_csa07_chowder_topSemiLepElectron_100pb.tab.list", discr[j], 1.0, tab );

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

  double epsb_true[50];
  double epsb_true_err[50];
  double epsc_true[50];
  double epsc_true_err[50];

  for (int i = 0; i<n; i++)
    {
      discr_err[i] = 0.0;

      epsl[i] = taggedL[i]/goodL;
      epsl_err[i] = epsl[i] / sqrt( taggedL[i] );

      epsb_true[i] = taggedB[i]/goodB;
      epsb_true_err[i] = epsb_true[i] / sqrt( taggedB[i] );

      epsc_true[i] = taggedC[i]/goodC;
      epsc_true_err[i] = epsc_true[i] / sqrt( taggedC[i] );

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



  h . Draw();
  
  //c1->SetFillColor(42);
  c1 . SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  //c1->GetFrame()->SetBorderSize(12);

  TGraphErrors gr_epsb( n, discr, epsb, discr_err, epsb_err );
  //TGraphErrors gr_epsc( n, discr, epsc, discr_err, epsc_err );

  TGraphErrors gr_epsl( n, discr, epsl, discr_err, epsl_err );

  TGraphErrors gr_epsb_true( n, discr, epsb_true, discr_err, epsb_true_err );
  //TGraphErrors gr_epsc_true( n, discr, epsc_true, discr_err, epsc_true_err );

  gr_epsb . Draw("LP");
  //gr_epsc . Draw("LP");
  //gr_epsl . Draw("LP");

  gr_epsb_true . Draw("LP");
  //gr_epsc_true . Draw("LP");


  gr_epsb . SetLineWidth(2);
  gr_epsb . SetLineColor(4);
  gr_epsb . SetMarkerColor(4);

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
  gr_epsb_true . SetMarkerStyle(1);

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

  c1 . SaveAs( "bTagEff_cFixed.eps" );
  //c1 . SaveAs( "bTagEff.eps" );
  //c1 . SaveAs( "cTagEff.eps" );
  
  //c1 . SaveAs( "bcContour.eps" );




}

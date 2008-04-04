// -*- C++ -*-
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      TtTagConsistencyEfficiency
// 
/**\class TtTagConsistencyEfficiency TtTagConsistencyEfficiency.cc RecoBTag/PerformanceMeasurements/src/TtTagConsistencyEfficiency.cc

 Description: Tag counting method for b-,c- and light tagging efficiency measurement with ttbar semileptonic events

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Fri Jun 29 14:53:10 CDT 2007
// $Id: TtTagConsistencyEfficiency.h,v 1.1.2.2 2008/03/21 21:09:17 kukartse Exp $
//
//

#ifndef PerformanceMeasurementsTtTagConsistencyEfficiency
#define PerformanceMeasurementsTtTagConsistencyEfficiency

#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <cstdio>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TH2F.h"

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyFitConfig.h"



using namespace std;



class TtTagConsistencyEfficiency
{

public:

  typedef struct _tableLine
  {
    _tableLine();
    double discr;
    double N[5];
    double Nb, Nc, Nl, Nx, Nbtag, Nctag, Nltag, Nxtag;    
  } tableLine;

  typedef struct _tableContent
  {
    int N[5];
    int Nb, Nc, Nl, Nbtag, Nctag, Nltag;    
  } tableContent;

  typedef struct var2
  {
    string name;
    string title;
    string formula[5];
  } BFormula;

  TtTagConsistencyEfficiency();
  ~TtTagConsistencyEfficiency(){};
  void cleanTable( tableContent & tab );
  void cleanTable( tableLine & tab );
  int readTable( const char * tabFileName, double discr );
  int readTable( vector<string> tabFileName, double discr );
  int readTable( vector<string> tabFileName, double discr, vector<double> weight, tableContent & _tab );
  int updateTable( vector<string> tabFileName, double discr, double weight, tableLine & _tab );
  int updateTable( string tabListFileName, double discr, double weight, tableContent & _tab );
  int updateTable( string tabFileName, double discr, double weight, tableLine & _tab );
  int updateTableFromList( string tabListFileName, double discr, double weight, tableLine & _tab );

  int getN( int index );
  int getNb( void );
  int getNc( void );
  int getNl( void );
  int getNbtag( void );
  int getNctag( void );
  int getNltag( void );

  RooFitResult minTest2( int n0, int n1, int n2, int n3, int n4, const char * option = "msh",  double eb = 0.65, double ec = 0.2, double el = 0.02, double xsec_ttbar = 830.0, double bgf = 0.2, double eb_min = 0.0, double eb_max = 1.0, double ec_min = 0.0, double ec_max = 1.0, double el_min = 0.0, double el_max = 1.0, double xsec_ttbar_min = 50.0, double xsec_ttbar_max = 1100.0, double bgf_min = 0.0, double bgf_max = 100.0 );

  RooFitResult fit( TtTagConsistencyFitConfig & config );
  RooFitResult fit_signal( TtTagConsistencyFitConfig & config, TH2F * _contour = NULL );
  RooFitResult fit_optimized( TtTagConsistencyFitConfig & config );

  TtTagConsistencyEfficiency::BFormula getFormula( int i, int j, int k, const char * suffix = "" );

  int sumFijk( const char * mcFileName );
  int readFijk( const char * mcFileName, string dataType, int minValue = 0 );
  int readFijk( vector<string> mcFileName, string dataType, int minValue = 0 );
  int readFijk( vector<string> mcFileName, vector<double> weight, string dataType, int minValue = 0 );
  int updateFijk( string mcListFileName, double weight, string dataType, int minValue = 0 );
  int updateFijk2( string mcListFileName, double weight, string dataType, int minValue = 0 );
  int test( void );

  int N[5];
  int Nb, Nc, Nl, Nbtag, Nctag, Nltag;
  int N_bg[5];
  int Nb_bg, Nc_bg, Nl_bg, Nbtag_bg, Nctag_bg, Nltag_bg;

  TH2F * contourBC;
  TH2F * contourBX;

  struct Fijk
  {
    Fijk();
    vector<int> i;
    vector<int> j;
    vector<int> k;
    //vector<int> value;
    std::map<int, std::map<int, std::map<int,int> > > value;
    int nEntries;
    int sum;
  };

  struct Fijk Fijk_sig;
  struct Fijk Fijk_bg;

protected:

  double _discr;
  bool inputFileRead;
  bool redo_fit_setup;

};

#endif

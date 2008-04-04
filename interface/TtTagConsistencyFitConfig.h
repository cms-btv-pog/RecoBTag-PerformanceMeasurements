// -*- C++ -*-
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      TtTagConsistencyFitConfig
// 
/**\class TtTagConsistencyFitConfig TtTagConsistencyFitConfig.cc RecoBTag/PerformanceMeasurements/src/TtTagConsistencyFitConfig.cc

 Description: Tag counting method for b-,c- and light tagging efficiency measurement with ttbar semileptonic events

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Fri Jun 29 14:53:10 CDT 2007
// $Id: TtTagConsistencyFitConfig.h,v 1.1.2.2 2008/03/21 21:09:17 kukartse Exp $
//
//

#ifndef PerformanceMeasurementsTtTagConsistencyFitConfig
#define PerformanceMeasurementsTtTagConsistencyFitConfig



class TtTagConsistencyFitConfig
{
public:
  
  TtTagConsistencyFitConfig();
  ~TtTagConsistencyFitConfig(){};

  void setConfig( int _n1, int _n2, int _n3, double _eb, double _ec, double _el );

  int n0, n1, n2, n3, n4;
  const char * option;
  double eb, ec, el;
  double eb_min, eb_max, ec_min, ec_max, el_min, el_max;
  double xsec_ttbar;
  double xsec_ttbar_min, xsec_ttbar_max;
  double bgf0, bgf1, bgf2, bgf3, bgf4;
  double bgf0_min, bgf0_max, bgf1_min, bgf1_max, bgf2_min, bgf2_max, bgf3_min, bgf3_max, bgf4_min, bgf4_max;
  double lumi;
  double efficiency;
};

#endif

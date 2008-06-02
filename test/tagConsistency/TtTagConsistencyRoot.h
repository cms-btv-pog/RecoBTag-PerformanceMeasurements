// -*- C++ -*-
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      TtTagConsistencyRoot
// 
/**\class TtTagConsistencyRoot TtTagConsistencyRoot.cc RecoBTag/PerformanceMeasurements/src/TtTagConsistencyRoot.cc

 Description: Tag counting method for b-,c- and light tagging efficiency measurement with ttbar semileptonic events

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Fri Jun 29 14:53:10 CDT 2007
// $Id: TtTagConsistencyRoot.h,v 1.1.2.1 2008/04/04 23:54:48 kukartse Exp $
//
//

#ifndef PerformanceMeasurementsTtTagConsistencyRoot
#define PerformanceMeasurementsTtTagConsistencyRoot


#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyEfficiency.h"



class TtTagConsistencyRoot
{

public:

  int scan_plot( string data_list_file, string prefix,
		 double discr_min=0.0, double discr_max=1.0, double discr_step=0.1,
		 double mistag_min=0.0, double mistag_max=1.0 );
  int scan_plot_d( string data_list_file, string prefix,
		   double discr_min=0.0, double discr_max=1.0, double discr_step=0.1,
		   double mistag_min=0.0, double mistag_max=1.0,
		   double frac1 = 0.3, double frac2 = 0.024, double frac3 = 0.006,
		   bool use_qcd = false );
  void setTDRStyle( void );
  void fwlite_test( void );

private:


};

#endif

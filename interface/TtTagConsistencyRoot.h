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
// $Id: TtTagConsistencyRoot.h,v 1.1.2.2 2008/03/21 21:09:17 kukartse Exp $
//
//

#ifndef PerformanceMeasurementsTtTagConsistencyRoot
#define PerformanceMeasurementsTtTagConsistencyRoot


#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyEfficiency.h"



class TtTagConsistencyRoot
{

public:

  int scan_plot( void );

private:


};

#endif

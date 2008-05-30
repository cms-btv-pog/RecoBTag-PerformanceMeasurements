/* -*- mode: c++ -*- */

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyFitConfig.h"
#include <iostream>

using namespace std;

TtTagConsistencyFitConfig::TtTagConsistencyFitConfig()
{
  n0 = 0;
  n1 = 0;
  n2 = 0;
  n3 = 0;
  n4 = 0;
  n0_d = 0.0;
  n1_d = 0.0;
  n2_d = 0.0;
  n3_d = 0.0;
  n4_d = 0.0;
  option = "msh";
  eb = 0.65;
  ec = 0.2;
  el = 0.01;
  xsec_ttbar = 561.0;  // pb
  xsec_ttbar_min = 0.0;
  xsec_ttbar_max = 100000.0;
  bgf0 = 3.64;
  bgf1 = 0.217;
  bgf2 = 0.024;
  bgf3 = 0.006;
  bgf4 = 0.002;
  bgf0_min = 0.0;
  bgf0_max = 5.0;
  bgf1_min = 0.0;
  bgf1_max = 5.0;
  bgf2_min = 0.0;
  bgf2_max = 5.0;
  bgf3_min = 0.0;
  bgf3_max = 5.0;
  bgf4_min = 0.0;
  bgf4_max = 5.0;
  eb_min = 0.0;
  eb_max = 1.0;
  ec_min = 0.0;
  ec_max = 1.0;
  el_min = 0.0;
  el_max = 1.0;
  
  lumi = 100.0; // 1/pb
  efficiency = 0.169;
}



void TtTagConsistencyFitConfig::dump()
{
  cout << endl << "=============> fit config dump " << endl << endl;
  cout << "n0 = " << n0 << endl;
  cout << "n1 = " << n1 << endl;
  cout << "n2 = " << n2 << endl;
  cout << "n3 = " << n3 << endl;
  cout << "n4 = " << n4 << endl;
  cout << "n0_d = " << n0_d << endl;
  cout << "n1_d = " << n1_d << endl;
  cout << "n2_d = " << n2_d << endl;
  cout << "n3_d = " << n3_d << endl;
  cout << "n4_d = " << n4_d << endl;
  cout << "eb_min, eb, eb_max = " << eb_min << "	" << eb << "	" << eb_max << endl;
  cout << "ec_min, ec, ec_max = " << ec_min << "	" << ec << "	" << ec_max << endl;
  cout << "el_min, el, el_max = " << el_min << "	" << el << "	" << el_max << endl;
  cout << "xsec_ttbar_min, xsec_ttbar, xsec_ttbar_max = " << xsec_ttbar_min << "	" << xsec_ttbar << "	" << xsec_ttbar_max << endl;
  cout << "bgf0_min, bgf0, bgf0_max = " << bgf0_min << "	" << bgf0 << "	" << bgf0_max << endl;
  cout << "bgf1_min, bgf1, bgf1_max = " << bgf1_min << "	" << bgf1 << "	" << bgf1_max << endl;
  cout << "bgf2_min, bgf2, bgf2_max = " << bgf2_min << "	" << bgf2 << "	" << bgf2_max << endl;
  cout << "bgf3_min, bgf3, bgf3_max = " << bgf3_min << "	" << bgf3 << "	" << bgf3_max << endl;
  cout << "bgf4_min, bgf4, bgf4_max = " << bgf4_min << "	" << bgf4 << "	" << bgf4_max << endl;
  cout << endl << "====================================================" << endl;
}



void TtTagConsistencyFitConfig::setConfig( int _n1, int _n2, int _n3, double _eb, double _ec, double _el )
  {
    n1=_n1;
    n2=_n2;
    n3=_n3;
    eb=_eb;
    ec=_ec;
    el=_el;

    eb_max=_eb+0.2;
    ec_max=_ec+0.2;
    el_max=_el+0.2;
    eb_min=_eb-0.2;
    ec_min=_ec-0.2;
    el_min=_el-0.2;

    bgf1 = 0.0;
    bgf2 = 0.0;
    bgf3 = 0.0;

    lumi = 100.0;

  }



void TtTagConsistencyFitConfig::setConfig( double _n1, double _n2, double _n3, double _eb, double _ec, double _el )
  {
    n1_d=_n1;
    n2_d=_n2;
    n3_d=_n3;
    eb=_eb;
    ec=_ec;
    el=_el;

    eb_max=_eb+0.2;
    ec_max=_ec+0.2;
    el_max=_el+0.2;
    eb_min=_eb-0.2;
    ec_min=_ec-0.2;
    el_min=_el-0.2;

    bgf1 = 0.0;
    bgf2 = 0.0;
    bgf3 = 0.0;

    lumi = 100.0;

  }

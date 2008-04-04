/* -*- mode: c++ -*- */

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyFitConfig.h"


TtTagConsistencyFitConfig::TtTagConsistencyFitConfig()
  {
    n0 = 0;
    n1 = 0;
    n2 = 0;
    n3 = 0;
    n4 = 0;
    option = "msh";
    eb = 0.65;
    ec = 0.2;
    el = 0.01;
    xsec_ttbar = 561.0;  // pb
    xsec_ttbar_min = 0.0;
    xsec_ttbar_max = 10000.0;
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

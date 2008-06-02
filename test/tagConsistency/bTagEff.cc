#include <iostream>
#include <boost/program_options.hpp>

#include "RecoBTag/PerformanceMeasurements/test/tagConsistency/TtTagConsistencyRoot.h"

namespace po = boost::program_options;
using namespace std;

int main( int argc, char ** argv )
{
  string prefix;
  double discr_min, discr_max, discr_step, mistag_min, mistag_max;
  double frac1, frac2, frac3;
  bool use_qcd;
  po::options_description general("General options");
  general.add_options()
    ("help", "produce help message")
    ("test", po::value<string>(), "print test string")
    ("scan-plot", "tag eff scan plot")
    ("scan-double", "tag eff scan plot (new version with double type entries)")
    ("data-list-file", po::value<string>(), "file containing the list of input .tab files")
    ("prefix", po::value<string>(&prefix)->default_value("btag"), "prefix for output files")
    ("min-discr", po::value<double>(&discr_min)->default_value(0.0), "low limit for the b-tagging discriminant")
    ("max-discr", po::value<double>(&discr_max)->default_value(1.0), "high limit for the b-tagging discriminant")
    ("discr-step", po::value<double>(&discr_step)->default_value(0.1), "b-tagging discriminant scan step")
    ("min-mistag", po::value<double>(&mistag_min)->default_value(0.0), "low limit for the mistagging rate")
    ("max-mistag", po::value<double>(&mistag_max)->default_value(1.0), "high limit for the mistagging rate")
    ("frac1", po::value<double>(&frac1)->default_value(0.3), "background fraction in bin 1")
    ("frac2", po::value<double>(&frac2)->default_value(0.024), "background fraction in bin 2")
    ("frac3", po::value<double>(&frac3)->default_value(0.006), "background fraction in bin 3")
    ("use-qcd", po::value<bool>(&use_qcd)->default_value(false), "consider QCD background (true/false)")
    ("fwlite", "FWLite test")
    ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
    cout << general << "\n";
    return 1;
  }
  
  if (vm.count("test")) {
    cout << "Test: " 
	 << vm["test"].as<string>() << ".\n";
  }


  if (vm.count("scan-plot") && vm.count("data-list-file")) {
    TtTagConsistencyRoot _root;
    _root . scan_plot( vm["data-list-file"].as<string>(), vm["prefix"].as<string>(),
		       vm["min-discr"].as<double>(), vm["max-discr"].as<double>(), vm["discr-step"].as<double>(),
		       vm["min-mistag"].as<double>(), vm["max-mistag"].as<double>());
    return 0;
  }
  
  if (vm.count("scan-double") && vm.count("data-list-file")) {
    TtTagConsistencyRoot _root;
    _root . scan_plot_d( vm["data-list-file"].as<string>(), vm["prefix"].as<string>(),
			 vm["min-discr"].as<double>(), vm["max-discr"].as<double>(), vm["discr-step"].as<double>(),
			 vm["min-mistag"].as<double>(), vm["max-mistag"].as<double>(),
			 vm["frac1"].as<double>(), vm["frac2"].as<double>(), vm["frac3"].as<double>(),
			 vm["use-qcd"].as<bool>() );
    return 0;
  }
  
  if (vm.count("fwlite")) {
    TtTagConsistencyRoot _root;
    _root . fwlite_test();
    return 0;
  }
  
  else {
    cout << "Test: Welcome to the desert of the real!\n";
  }

  return 0;
}

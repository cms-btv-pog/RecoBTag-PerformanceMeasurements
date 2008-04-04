#include <iostream>
#include <boost/program_options.hpp>

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyRoot.h"

namespace po = boost::program_options;
using namespace std;

int main( int argc, char ** argv )
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("test", po::value<string>(), "print test string")
    ("scan-plot", "tag eff scan plot")
    ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }
  
  if (vm.count("test")) {
    cout << "Test: " 
	 << vm["test"].as<string>() << ".\n";
  }

  if (vm.count("scan-plot")) {
    cout << desc << "\n";
    TtTagConsistencyRoot _root;
    _root . scan_plot();
    return 0;
  }
  
  else {
    cout << "Test: Welcome to the desert of the real!\n";
  }

  return 0;
}

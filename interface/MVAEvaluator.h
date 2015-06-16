#ifndef RecoBTag_PerformanceMeasurements_MVAEvaluator_h
#define RecoBTag_PerformanceMeasurements_MVAEvaluator_h

#include <memory>
#include <cassert>
#include <vector>
#include <string>
#include <map>

#include <TROOT.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


class MVAEvaluator {

  public:
    MVAEvaluator(const std::string & method, const std::string & weightFile) :
      mva_method_(weightFile), mva_weightFile_(weightFile) {};
    ~MVAEvaluator() {};

   void bookReader(const std::vector<std::string> & variables, const std::vector<std::string> & spectators);
   float evaluate(const std::map<std::string,float> & variables);

  private:
    std::unique_ptr<TMVA::Reader> mva_reader_;

    std::map<std::string,float> mva_variables_;
    std::map<std::string,float> mva_spectators_;

    std::string mva_method_;
    std::string mva_weightFile_;
};

#endif // RecoBTag_PerformanceMeasurements_MVAEvaluator_h


#include "RecoBTag/PerformanceMeasurements/interface/MVAEvaluator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Utils/interface/TMVAZipReader.h"


void MVAEvaluator::bookReader(const std::vector<std::string> & variables, const std::vector<std::string> & spectators)
{
  // initialize the TMVA reader
  mva_reader_.reset(new TMVA::Reader("Color:Silent:Error"));
  mva_reader_->SetVerbose(false);

  // add input variables
  for(std::vector<std::string>::const_iterator it = variables.begin(); it!=variables.end(); ++it)
  {
    mva_variables_.insert( std::pair<std::string,size_t>(*it,0.) );
    mva_reader_->AddVariable((*it).c_str(), &mva_variables_.at(*it));
  }

  // add spectator variables
  for(std::vector<std::string>::const_iterator it = spectators.begin(); it!=spectators.end(); ++it)
  {
    mva_spectators_.insert( std::pair<std::string,size_t>(*it,0.) );
    mva_reader_->AddSpectator(it->c_str(), &mva_spectators_.at(*it));
  }

  // load the TMVA weights
  reco::details::loadTMVAWeights(mva_reader_.get(), mva_method_.c_str(), mva_weightFile_.c_str());
}


float MVAEvaluator::evaluate(const std::map<std::string,float> & variables)
{
  assert( variables.size() >= mva_variables_.size() );

  // set the input variable values
  for(std::map<std::string,float>::iterator it = mva_variables_.begin(); it!=mva_variables_.end(); ++it)
  {
    if (variables.count(it->first)>0)
      it->second = variables.at(it->first);
    else
      edm::LogError("MissingInputVariable") << "Variable " << it->first << " is missing from the list of input variables.";
  }

  // evaluate the MVA
  float value = mva_reader_->EvaluateMVA(mva_method_.c_str());

  return value;
}

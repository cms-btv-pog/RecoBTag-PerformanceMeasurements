#ifndef VARIABLEPARSER_H
#define VARIABLEPARSER_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class VariableParser {

  public:

    VariableParser(bool isMC = false);

    std::vector<std::string> parseGroupsAndVariables(std::vector<edm::ParameterSet> groupSet, std::vector<edm::ParameterSet> variableSet);
    std::vector<std::string> getStoredVariables();
    void printStoredVariables();
    void printGroups(std::vector<edm::ParameterSet> groupSet);
    void printVariables(std::vector<edm::ParameterSet> variableSet);
    void saveStoredVariablesToFile(std::string filename="storedVariables.log");
    bool isToBeStored(std::string variableName);
    bool isMC(){return isMC_;}
    void resetStoredVariables(){storedVariables_.clear();}

  private:

    std::vector<std::string> storedVariables_;
    bool isMC_;

    void resolveVariableName(std::string fullName, std::string* variableName, std::string* prefix);
    void parseGroups(std::vector<edm::ParameterSet> groupSet, std::vector<edm::ParameterSet> variableSet);
    void parseVariables(std::vector<edm::ParameterSet> variableSet);
    void addVariable(edm::ParameterSet variable, bool force = false, std::string prefix = "");

};


#endif

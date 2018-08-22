#ifndef VARIABLEPARSER_H
#define VARIABLEPARSER_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <unordered_set>

class VariableParser {

  public:

    VariableParser(bool isMC = false);

    std::unordered_set<std::string> parseGroupsAndVariables(std::vector<edm::ParameterSet> groupSet, std::vector<edm::ParameterSet> variableSet);
    std::unordered_set<std::string> getStoredVariables();
    std::unordered_set<std::string> getRunOptions();
    void printStoredVariables();
    void printGroups(std::vector<edm::ParameterSet> groupSet);
    void printVariables(std::vector<edm::ParameterSet> variableSet);
    void printRunOptions();
    void saveStoredVariablesToFile(std::string filename="storedVariables.log");
    bool isToBeStored(std::string variableName);
    bool runOption(std::string runOption);
    bool isMC(){return isMC_;}
    void resetStoredVariables(){storedVariables_.clear();}
    void resetRunOptions(){runOptions_.clear();}

  private:

    std::unordered_set<std::string> storedVariables_;
    std::unordered_set<std::string> runOptions_;
    bool isMC_;

    void resolveVariableName(std::string fullName, std::string* variableName, std::string* prefix);
    void parseGroups(std::vector<edm::ParameterSet> groupSet, std::vector<edm::ParameterSet> variableSet);
    void parseVariables(std::vector<edm::ParameterSet> variableSet);
    void addVariable(edm::ParameterSet variable, bool force = false, std::string prefix = "");
    void addRunOption(std::string runOption);

};


#endif

#include "RecoBTag/PerformanceMeasurements/interface/VariableParser.h"


#include<iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>


using namespace std;

VariableParser::VariableParser(bool isMC){
  storedVariables_.clear();
  isMC_ = isMC;
}

unordered_set<string> VariableParser::parseGroupsAndVariables(vector<edm::ParameterSet> groupSet, vector<edm::ParameterSet> variableSet)
{
  parseGroups(groupSet, variableSet);
  parseVariables(variableSet);
  return storedVariables_;
}

unordered_set<string> VariableParser::getStoredVariables()
{
  return storedVariables_;
}

void VariableParser::printStoredVariables()
{
  cout << endl << "The following list of variables will be stored: " << endl;
  for(unordered_set<string>::iterator variable = storedVariables_.begin(); variable != storedVariables_.end(); ++variable)
  {
    cout << *variable << endl;
  }
  cout << endl;
}

void VariableParser::printGroups(vector<edm::ParameterSet> groupSet)
{
  for( vector<edm::ParameterSet>::iterator group = groupSet.begin(); group != groupSet.end(); ++group )
  {
    cout << left << setw(20) << "Group: " << group->getParameter<string>("group") << endl;
    cout << left << setw(20) << "Store: " << group->getParameter<bool>("store") << endl;
    cout << left << setw(20) << "Description: " << group->getParameter<string>("description") << endl;
    vector<string> variables = group->getParameter<vector<string>>("variables");
    cout << "Variables: " << endl;
    for( vector<string>::iterator variable = variables.begin(); variable != variables.end(); ++variable)
    {
      cout << setw(20) << " " << *variable << endl;
    }
    cout << endl;
  }
}

void VariableParser::printVariables(vector<edm::ParameterSet> variableSet)
{
  for( vector<edm::ParameterSet>::iterator variable = variableSet.begin(); variable != variableSet.end(); ++variable )
  {
    cout << left << setw(20) << "Variable: " << variable->getParameter<string>("variable") << endl;
    cout << left << setw(20) << "Store: " << variable->getParameter<bool>("store") << endl;
    cout << left << setw(20) << "MC only: " << variable->getParameter<bool>("mconly") << endl;
    cout << left << setw(20) << "Description: " << variable->getParameter<string>("description") << endl;
    vector<string> requires = variable->getParameter<vector<string>>("requires");
    cout << "Requires: " << endl;
    for( vector<string>::iterator require = requires.begin(); require != requires.end(); ++require)
    {
      cout << setw(20) << " " << *require << endl;
    }
    cout << endl;
  }

}

void VariableParser::saveStoredVariablesToFile(string filename)
{
  ofstream storedVariablesFile;
  storedVariablesFile.open(filename);
  for(unordered_set<string>::iterator variable = storedVariables_.begin(); variable != storedVariables_.end(); ++variable)
  {
    storedVariablesFile << *variable << endl;
  }
  storedVariablesFile.close();
  cout << "List of stored variables has been saved as " << filename << endl;
}

bool VariableParser::isToBeStored(string variableName)
{
  return (find(storedVariables_.begin(), storedVariables_.end(), variableName) != storedVariables_.end());
}

void VariableParser::resolveVariableName(string fullName, string* variableName, string* prefix)
{
  size_t pos = fullName.find(".");
  if(pos != string::npos)
  {
    *prefix = fullName.substr(0, pos+1);
    *variableName = fullName.substr(pos+1);
  }
  else
  {
    *variableName = fullName;
    *prefix = "";
  }
//   cout << "Prefix: " << *prefix << endl;
//   cout << "VariableName: " << *variableName << endl;
  return;
}

void VariableParser::parseGroups(vector<edm::ParameterSet> groupSet, vector<edm::ParameterSet> variableSet)
{
  for(vector<edm::ParameterSet>::iterator group = groupSet.begin(); group != groupSet.end(); ++group)
  {
    if(group->getParameter<bool>("store")) // group is to be stored
    {
      vector<string> variables = group->getParameter<vector<string>>("variables");
      for(vector<string>::iterator variable = variables.begin(); variable != variables.end(); ++variable)
      {
        bool found = false;
        if(!isToBeStored(*variable)) // variable is not contained yet
        {
          string variableName, prefix;
          resolveVariableName(*variable, &variableName, &prefix); // resolve variable name to find the corresponding ParameterSet
          for(vector<edm::ParameterSet>::iterator variablePSet = variableSet.begin(); variablePSet != variableSet.end(); ++variablePSet) //Get the corresponding set for dependencies
          {
            if(variableName == variablePSet->getParameter<string>("variable"))
            {
              addVariable(*variablePSet, true, prefix); // Add (force) the variable to the storedVariables_ unordered_set
              found = true;
              break;
            }
          }
          if(!found)
            cout << "WARNING: No PSet could be found for the variable " << *variable << " of group " << group->getParameter<string>("group") << ". Variable will not be stored. Please check for the PSet in variables_cfi.py" << endl;
        }
      }
    }
  }
}

void VariableParser::parseVariables(vector<edm::ParameterSet> variableSet)
{
  for(vector<edm::ParameterSet>::iterator variablePSet = variableSet.begin(); variablePSet != variableSet.end(); ++variablePSet)
  {
    addVariable(*variablePSet);
  }
}

void VariableParser::addVariable(edm::ParameterSet variablePSet, bool force, string prefix)
{
  string variableName = prefix + variablePSet.getParameter<string>("variable");
  bool store = variablePSet.getParameter<bool>("store") || force;
  bool mconly = variablePSet.getParameter<bool>("mconly");
  if(store && !( !isMC() && mconly )) // variable is to be stored (checking also MC/data setting)
  {
    if(!isToBeStored(variableName)) // variable is not contained yet
    {
      storedVariables_.insert(variableName);
      vector<string> requires = variablePSet.getParameter<vector<string>>("requires"); // Add required variables (currently does not support nested requirements... souldn't be too much of a problem)
      for(vector<string>::iterator require = requires.begin(); require != requires.end(); ++require)
      {
        if(!isToBeStored(prefix + *require))
          storedVariables_.insert(prefix + *require);
      }
    }
  }
}

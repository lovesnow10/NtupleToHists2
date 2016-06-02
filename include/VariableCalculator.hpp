#ifndef VARIABLECALCULATOR_HPP_
#define VARIABLECALCULATOR_HPP_
#include "TTree.h"
#include <string>
#include <iostream>
#include <map>

#include "tools.hpp"

class VariableCalculator {
private:
  TTree *mEvent;
  std::map<string, float> mVariablesMap;
public:
  VariableCalculator () {};
  virtual ~VariableCalculator () {};

  void CalculateVariables(TTree* event);
  float GetVarValue(string VarName);

};

#endif

#ifndef MAKEHISTS_HPP_
#define MAKEHISTS_HPP_

#include "HistStore.hpp"
#include "ConfigParser.hpp"
#include "ttbbNLO_syst.hpp"
#include "TTreeFormulaContainer.hpp"
#include "DSHandler.hpp"

#include "TFile.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "tools.hpp"

using namespace std;

class MakeHists {
private:
  TTree *mEvent;
  ConfigParser *mConfig;
  HistStore *hs;

  //[region][sample]
  std::map<string, std::map<string, float>> mRawYields;
  std::map<string, std::map<string, float>> mWeightedYields;

  void InitYields(DSHandler *ds);
  void FillYields();

public:
  MakeHists(){};
  virtual ~MakeHists(){};
  bool initialize(ConfigParser *config, DSHandler *ds);
  bool run(TTree *event, float norm, float weight_ttbb_Nominal,
           TTreeFormulaContainer *formulas);
  bool finalize(TFile *fFile);
};

#endif

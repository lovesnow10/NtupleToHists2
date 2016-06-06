/*
  A class inherit from AnalyseTemplate, loop over entries and generate hists
*/
#ifndef HPLUSRUN_HPP_
#define HPLUSRUN_HPP_

// My Headers
#include "HistStore.hpp"
#include "ttbbNLO_syst.hpp"
#include "ConfigParser.hpp"
#include "XSHelper.hpp"
#include "DSHandler.hpp"
#include "TTreeFormulaContainer.hpp"
#include "MakeHists.hpp"
#include "PrintInfo.hpp"
#include "tools.hpp"

#include "TChain.h"

// Common Headers
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>

using namespace std;

class HplusRun {
private:
  TChain *mWorker;
  TChain *mHelpWorker;
  ConfigParser *mConfig;
  XSHelper *mXSHelper;
  DSHandler *mDS;
  TTreeFormulaContainer *mFormulas;

  MakeHists *mMH;
  PrintInfo *mPI;

  TFile *fFile;

  long mMaxProcess;
  long mProcessed;

  std::map<string, bool> bControl;

  string mMCTree;
  string mDTTree;

  int mArgc;
  char **mArgv;

  // Cannot be copied
  HplusRun(HplusRun &hg);

public:
  HplusRun(int argv, char **argc);
  virtual ~HplusRun(){};

  bool initialize();
  bool run();
  bool finalize();
};

#endif

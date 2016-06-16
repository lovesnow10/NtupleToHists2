#ifndef TOOLS_HPP_
#define TOOLS_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"

using namespace std;

class Tools {
private:
  Tools(){};
  Tools(Tools &rhs){};
  ~Tools(){};
  Tools &operator=(Tools &rhs);

public:
  static Tools &Instance() {
    static Tools instance;
    return instance;
  }
  string GenName(string VarName, string Region, string Sample);
  float GetWeight(int DSID);
  string GetSampleType(int DSID);
  std::vector<string> GrabRootFiles(string path);
  string GetDataYear(int runNumber);

  template <typename T> T *GetTreeValue(TTree *fTree, string fVar) {
    T *fValue = 0x0;
    TLeaf *fLeaf = fTree->GetLeaf(fVar.c_str());
    if (fLeaf != 0x0) {
      fValue = (T *)fLeaf->GetValuePointer();
    }
    return fValue;
  }

  string GetValueType(TTree *fTree, string fVar) {
    string fType = "";
    TLeaf *fLeaf = fTree->GetLeaf(fVar.c_str());
    if (fLeaf != 0x0) {
      fType = fLeaf->GetTypeName();
    }
    return fType;
  }

  template <typename T> T GetSumUp(TTree *fTree, string fVar) {
    T sum = 0;
    long nentries = fTree->GetEntries();
    for (long ientry = 0; ientry < nentries; ++ientry) {
      fTree->GetEntry(ientry);
      sum = sum + *GetTreeValue<T>(fTree, fVar);
    }
    return sum;
  }
};
#endif

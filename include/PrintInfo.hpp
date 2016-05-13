#ifndef PRINTINFO_HPP_
#define PRINTINFO_HPP_

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <ctime>

#include "TTree.h"

#include "tools.hpp"

using namespace std;

class PrintInfo {
private:
  TTree *mEvent;
  ofstream mOut;
public:
  PrintInfo () {};
  virtual ~PrintInfo () {};
  bool initialize();
  bool run(TTree *event, std::map<string, float> weights, long ientry);
  bool finalize();
};


#endif

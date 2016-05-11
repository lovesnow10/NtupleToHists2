/*
  A simple parser convert TXT file to a Config map.
  Example like this:
  REGION 2jex2bex
  CUT nJets==2&&nBTags==2
  VAR pT_jet1
*/
#ifndef CONFIGPARSER_HPP_
#define CONFIGPARSER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

using namespace std;

class ConfigParser {
private:
  std::vector<string> mRegions;
  std::map<string, string> mRegionCutMap;
  std::map<string, std::map<string, std::vector<float>>> mRegionVarsMap;
  std::map<string, std::vector<string>> mCommonSettingMap;
  std::map<string, vector<string>> mRegionVarMap;

  string mConfigFileName;

  bool initialize();
  string FixString(string s);
  int GetNSpace(string s);

public:
  ConfigParser(string configfile);
  virtual ~ConfigParser(){};

  std::vector<string> GetRegions();
  string GetRegionCut(string RegionName);
  std::vector<float> GetVarBinning(std::string RegionName, std::string VarName);
  std::vector<string> GetCommonSetting(string settingname);
  vector<string> GetRegionVars(string RegionName);
};

#endif

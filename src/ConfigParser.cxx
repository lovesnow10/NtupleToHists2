#include "ConfigParser.hpp"
#include <sstream>
#include <algorithm>

ConfigParser::ConfigParser(string configfile) : mConfigFileName(configfile) {
  if (!initialize()) {
    std::cout << "Cannot Parse the File!" << std::endl;
    exit(-1);
  }
  std::cout << "ConfigParser:: "
            << "Initialization Successed!" << std::endl;
}

string ConfigParser::FixString(string s) {
  string ss = s.substr(0, s.find_first_of('#'));
  replace(ss.begin(), ss.end(), '"', ' ');
  if (ss == "")
    return ss;
  if (ss.find_first_not_of(' ')) {
    ss = ss.substr(ss.find_first_not_of(' '), ss.find_last_not_of(' '));
  } else {
    ss = ss.substr(ss.find_first_not_of(' '), ss.find_last_not_of(' ') + 1);
  }
  return ss;
}

int ConfigParser::GetNSpace(string s) { return count(s.begin(), s.end(), ' '); }

bool ConfigParser::initialize() {
  ifstream in(mConfigFileName.c_str());
  if (!in) {
    std::cout << "ConfigParser:: initialize:: "
              << "Cannot Open the File!" << std::endl;
    return false;
  }
  string line;
  string CurrentRegion;
  while (getline(in, line)) {
    line = FixString(line);
    if (line == "")
      continue;
    stringstream ss(line);
    int nSpace = this->GetNSpace(line);
    std::vector<string> tmpS;
    if (nSpace < 5) {
      tmpS.resize(6);
      for (int i = 0; i < nSpace + 1; i++) {
        ss >> tmpS[i];
      }
      for (int i = nSpace + 1; i < 6; i++) {
        tmpS[i] = "";
      }
    } else {
      tmpS.resize(nSpace + 1);
      for (int i = 0; i < nSpace + 1; i++) {
        ss >> tmpS[i];
      }
    }
    // Deal with COMMON
    if (tmpS[0] == "COMMON" && tmpS[1] != "") {
      for (int i = 2; i < nSpace + 1; i++) {
        mCommonSettingMap[tmpS[1]].push_back(tmpS[i]);
      }
      continue;
    } else if (tmpS[0] == "COMMON" && tmpS[1] == "") {
      std::cout << "ConfigParser:: initialize:: "
                << "Find COMMON But without Any Setting" << std::endl;
      return false;
    }
    // Deal with REGION
    if (tmpS[0] == "REGION" && tmpS[1] != "") {
      mRegions.push_back(tmpS[1]);
      CurrentRegion = tmpS[1];
      continue;
    } else if (tmpS[0] == "REGION" && tmpS[1] != "") {
      std::cout << "ConfigParser:: initialize:: "
                << "Find Region without Name!" << std::endl;
      return false;
    }
    // Deal with Vars
    if (tmpS[0] == "VAR" && tmpS[4] == "") {
      std::cout << "ConfigParser:: initialize:: "
                << "Find VAR with Lack of Args" << std::endl;
      return false;
    } else if (tmpS[0] == "VAR") {
      mRegionVarMap[CurrentRegion].push_back(tmpS[1]);
      std::vector<float> bins;
      bins.push_back(atof(tmpS[2].c_str()));
      bins.push_back(atof(tmpS[3].c_str()));
      bins.push_back(atof(tmpS[4].c_str()));
      mRegionVarsMap[CurrentRegion][tmpS[1]] = bins;
      continue;
    }
    // Deal with CUT
    if (tmpS[0] == "CUT" && tmpS[1] != "") {
      mRegionCutMap[CurrentRegion] = tmpS[1];
      continue;
    } else {
      std::cout << "ConfigParser:: initialize:: "
                << "Find CUT without Content!" << std::endl;
      return false;
    }
  }
  return true;
}

std::vector<string> ConfigParser::GetRegions() { return mRegions; }

string ConfigParser::GetRegionCut(string RegionName) {
  if (mRegionCutMap.find(RegionName) != mRegionCutMap.end()) {
    return mRegionCutMap.at(RegionName);
  } else {
    printf("ConfigParser:: GetRegions:: Cannot find REGION %s!!\n",
           RegionName.c_str());
    return "";
  }
}

std::vector<float> ConfigParser::GetVarBinning(std::string RegionName,
                                               std::string VarName) {
  std::vector<float> bins;
  bins.clear();
  if (mRegionVarsMap.find(RegionName) != mRegionVarsMap.end()) {
    std::map<string, std::vector<float>> var = mRegionVarsMap.at(RegionName);
    if (var.find(VarName) != var.end()) {
      bins = var.at(VarName);
    } else {
      printf(
          "ConfigParser:: GetVarBinning:: Cannot Find VAR %s in REGION %s!!\n",
          VarName.c_str(), RegionName.c_str());
    }
  } else {
    printf("ConfigParser:: GetVarBinning:: Cannot find REGION %s!!\n",
           RegionName.c_str());
  }
  return bins;
}

std::vector<string> ConfigParser::GetCommonSetting(string settingname) {
  std::vector<string> v;
  if (mCommonSettingMap.find(settingname) != mCommonSettingMap.end()) {
    v = mCommonSettingMap.at(settingname);
  } else {
    printf("ConfigParser:: GetCommonSetting:: Cannot find COMMON %s\n",
           settingname.c_str());
  }
  return v;
}

std::vector<string> ConfigParser::GetRegionVars(string RegionName) {
  std::vector<string> var;
  if (mRegionVarMap.find(RegionName) != mRegionVarMap.end()) {
    var = mRegionVarMap.at(RegionName);
  } else {
    printf("ConfigParser:: GetRegionVars:: Connot find REGION %s\n",
           RegionName.c_str());
  }
  return var;
}

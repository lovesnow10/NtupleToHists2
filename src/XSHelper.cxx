#include "XSHelper.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

XSHelper::XSHelper(string xsfile) {

  ifstream in(xsfile.c_str());
  string line;
  while (getline(in, line)) {
    if (line[0] == '#' || line.empty())
      continue;
    stringstream ss(line);
    string dsid, xsection, kfactor;
    ss >> dsid >> xsection >> kfactor;
    float xs = atof(xsection.c_str());
    float kf = atof(kfactor.c_str());
    m_XSMap[dsid] = xs * kf;
  }
}

float XSHelper::GetXS(const string key) {
  if (m_XSMap.find(key) == m_XSMap.end()) {
    cout << "Cannot find Cross Section for dsid " << key << endl;
    return 1.0;
  }
  return m_XSMap.at(key);
}

const float XSHelper::GetXS(const string key) const {
  if (m_XSMap.find(key) == m_XSMap.end()) {
    cout << "Cannot find Cross Section for dsid " << key << endl;
    return 0;
  }
  return m_XSMap.at(key);
}

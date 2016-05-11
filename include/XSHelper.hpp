#ifndef XSHELPER_HPP_
#define XSHELPER_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

class XSHelper {
public:
  XSHelper(string xsfile);
  ~XSHelper(){};

  float GetXS(const string key);
  const float GetXS(const string key) const;

private:
  XSHelper(XSHelper &xs){};
  map<string, float> m_XSMap;
};

#endif //#XSHELPER

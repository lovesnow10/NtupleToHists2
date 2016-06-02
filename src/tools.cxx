#include "tools.hpp"
#include "TList.h"
#include "TString.h"
#include "TSystemDirectory.h"

string Tools::GenName(string VarName, string Region, string Sample) {
  string hname = "hist";
  hname = hname + "_" + VarName + "_" + Region + "_" + Sample;
  return hname;
}

string Tools::GetSampleType(int DSID) {
  if (DSID == 0) {
    return "DATA";
  }
  if (DSID >= 361300 && DSID <= 361371) { // Wjets
    return "Fakes";
  }
  if (DSID >= 361372 && DSID <= 361491) {
    return "ZJets";
  }
  if (DSID >= 361063 && DSID <= 361086 &&
      !(DSID == 361070 || DSID == 361081 || DSID == 361082 || DSID == 361083)) {
    return "Dibosons";
  }
  if (DSID == 410015 || DSID == 410016) {
    return "SingleTop";
  }
  if (DSID == 410000 || DSID == 410009 || DSID == 410120 || DSID == 410121) {
    return "ttbar";
  }
  if (DSID >= 343365 && DSID <= 343367) {
    return "ttH";
  }
  if (DSID == 341541) {
    return "Hplus200";
  }
  if (DSID == 341542) {
    return "Hplus225";
  }
  if (DSID == 341543) {
    return "Hplus250";
  }
  if (DSID == 341544) {
    return "Hplus275";
  }
  if (DSID == 341545) {
    return "Hplus300";
  }
  if (DSID == 341546) {
    return "Hplus350";
  }
  if (DSID == 341547) {
    return "Hplus400";
  }
  if (DSID == 341548) {
    return "Hplus500";
  }
  if (DSID == 341549) {
    return "Hplus600";
  }
  if (DSID == 341550) {
    return "Hplus700";
  }
  if (DSID == 341551) {
    return "Hplus800";
  }
  if (DSID == 341552) {
    return "Hplus900";
  }
  if (DSID == 341553) {
    return "Hplus1000";
  }
  if (DSID == 341554) {
    return "Hplus1200";
  }
  if (DSID == 341555) {
    return "Hplus1400";
  }
  if (DSID == 341556) {
    return "Hplus1600";
  }
  if (DSID == 341557) {
    return "Hplus1800";
  }
  if (DSID == 341558) {
    return "Hplus2000";
  }
  if ((DSID >= 410066 && DSID <= 410068) || (DSID >= 410073 && DSID <= 410075) || (DSID >= 410111 && DSID <= 410116)) {
    return "ttV";
  }
  if (DSID == 410050)
  {
    return "tZ";
  }
  if (DSID == 410011 || DSID == 410012 || DSID == 410025 || DSID == 410026) {
    return "Fakes"; // SingleTop
  }
  if (DSID == 361070 || DSID == 361081 || DSID == 361082 || DSID == 361083) {
    return "Fakes"; // Dibosons
  }
  return "Unknown";
}

float Tools::GetWeight(int DSID) {
  float weight = 1.0;
  string type = GetSampleType(DSID);
  if (type == "ZJets") {
    int tmpID = DSID - 361371;
    if (DSID <= 361443) {
      int tmpType = tmpID % 3;
      if (tmpType == 1) {
        weight = 0.878;
      } else if (tmpType == 0) {
        weight = 1.135;
      }
    } else {
      int tmpType = tmpID % 2;
      if (tmpType == 1) {
        weight = 0.878;
      } else if (tmpType == 0) {
        weight = 1.135;
      }
    }
  }
  return weight;
}

std::vector<string> Tools::GrabRootFiles(string path) {
  TSystemDirectory dir(path.c_str(), path.c_str());
  TList *files = dir.GetListOfFiles();
  std::vector<string> rootfiles;
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file = (TSystemFile *)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.Contains(".root")) {
        rootfiles.push_back(path + "/" + fname.Data());
      }
    }
  }
  return rootfiles;
}

#include "DSHandler.hpp"
#include "TString.h"
#include "TCollection.h"
#include "TSystem.h"
#include <fstream>
#include <regex>

using namespace std;

DSHandler::DSHandler(string InputTxt) {
  ifstream in(InputTxt.c_str());
  if (!in) {
    printf("DSHandler:: Cannot Open File %s\n", InputTxt.c_str());
    exit(-1);
  }
  mPathsVec.clear();
  mSampleTypesVec.clear();
  mSampleTypesVec.push_back("Fakes");
  string line;
  int tmpID = 0;
  while (getline(in, line)) {
    mPathsVec.push_back(line);
    string ID = GetSampleID(line);
    if (ID != "") {
      string tmpType;
      if (ID.length() == 8) {
        tmpType = "DATA";
      } else {
        tmpType = Tools::Instance().GetSampleType(atoi(ID.c_str()));
      }
      if (find(mSampleTypesVec.begin(), mSampleTypesVec.end(), tmpType) ==
          mSampleTypesVec.end()) {
        mSampleTypesVec.push_back(tmpType);
      }
      mSampleTypeMap[tmpID++] = tmpType;
    }
  }
  this->Initialize();
}

void DSHandler::AddPath(string path) {
  mPathsVec.push_back(path);
  //  std::vector<string> tempVec = this->GrabRootFiles(path);
  //  mFilesMap.insert(make_pair(mTotalDS++, tempVec));
}

const std::vector<string> &DSHandler::GetFiles(int nDS) const {
  if (nDS < mTotalDS) {
    return mFilesMap.at(nDS);
  } else {
    printf("DSHandler:: GetFiles:: %i is Out-Of-Range!!\n", nDS);
    return tmpNone;
  }
}

bool DSHandler::Initialize() {
  mCurrentDS = 0;
  mFilesMap.clear();
  mTotalDS = mPathsVec.size();
  tmpNone.clear();
  int tmpID = 0;
  for (auto path : mPathsVec) {
    std::vector<string> tmpVec = GrabRootFiles(path);
    mFilesMap.insert(make_pair(tmpID++, tmpVec));
  }
  printf("DSHandler:: Initialize:: Initialization Succeed!\n");
  return true;
}

std::vector<string> DSHandler::GrabRootFiles(string path) {
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

void DSHandler::Clear() {
  mPathsVec.clear();
  mFilesMap.clear();
  mTotalDS = 0;
  mCurrentDS = 0;
  printf("DSHandler:: Clear:: Successfully Clear All!\n");
}

void DSHandler::DeletePath(string path) {
  auto _result = find(mPathsVec.begin(), mPathsVec.end(), path);
  if (_result == mPathsVec.end()) {
    printf("DSHandler:: DeletePath:: Cannot Find %s!!\n", path.c_str());
    return;
  }
  mPathsVec.erase(_result);
  this->Initialize();
  printf("DSHandler:: DeletePath:: Successfully Delete %s!\n", path.c_str());
}

const std::vector<string> &DSHandler::Next() {
  if (mCurrentDS < mTotalDS)
    return this->GetFiles(mCurrentDS++);
  else
    return tmpNone;
}

string DSHandler::GetSampleID(string path) {
  const regex pattern(".*\\.(\\d+)\\..*");
  std::match_results<std::string::const_iterator> result;
  bool valid = std::regex_match(path, result, pattern);
  if (!valid) {
    return "";
  } else {
    return result[1];
  }
}

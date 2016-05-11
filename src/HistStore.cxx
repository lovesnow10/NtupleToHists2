#include "HistStore.hpp"

using namespace std;

HistStore::HistStore() {}
void HistStore::AddHist(string hname, float nBins, float xlow, float xup) {
  TH1F *hist = new TH1F(hname.c_str(), hname.c_str(), nBins, xlow, xup);
  hist->Sumw2();
  hist->SetDirectory(0);
  mHistMap.insert(make_pair(hname, hist));
}

void HistStore::AddHist2D(string hname, float nBinsx, float xlow, float xup,
                          float nBinsy, float ylow, float yup) {
  TH2F *hist = new TH2F(hname.c_str(), hname.c_str(), nBinsx, xlow, xup, nBinsy,
                        ylow, yup);
  hist->SetDirectory(0);
  mHist2DMap.insert(make_pair(hname, hist));
}

TH1F *HistStore::GetHist(string hname) {
  if (mHistMap.find(hname) == mHistMap.end()) {
    printf("HistStore:: GetHist:: Cannot Find Histogram %s\n", hname.c_str());
    return nullptr;
  }
  return mHistMap.at(hname);
}

TH2F *HistStore::GetHist2D(string hname) {
  if (mHist2DMap.find(hname) == mHist2DMap.end()) {
    printf("HistStore:: GetHist2D:: Cannot Find Historam %s\n", hname.c_str());
    return nullptr;
  }
  return mHist2DMap.at(hname);
}

void HistStore::SaveAllHists(TFile *fFile) {
  mDir = (TDirectory *)fFile;
  for (auto hist : mHistMap) {
    hist.second->SetDirectory(mDir);
    hist.second->Write();
    printf("HistStore:: SaveAll:: Histogram %s Successfully Saved\n",
           hist.first.c_str());
  }
  for (auto hist : mHist2DMap) {
    hist.second->SetDirectory(mDir);
    hist.second->Write();
    printf("HistStore:: SaveAll:: Histogram %s Successfully Saved\n",
           hist.first.c_str());
  }
}

bool HistStore::HasHist(string hname) {
  if (mHistMap.find(hname) == mHistMap.end())
    return false;
  else
    return true;
}

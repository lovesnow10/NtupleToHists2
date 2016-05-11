#ifndef NNLOREWEIGHTER_HPP_
#define NNLOREWEIGHTER_HPP_

#include <string>
#include "TH1.h"
#include "TFile.h"

class NNLOReweighter {

  public:
    NNLOReweighter(int sampleID=0,std::string dirName="NNLOReweighter/data/");
    ~NNLOReweighter();
    void SetInputDirectory(std::string dirName="NNLOReweighter/data/");
    void SetInputFile(std::string fileName); // to specify a non-default file name
    void SetSampleID(int sampleID);
    void Init();
    float GetTtbarPtWeight(float ttbar_pt); // in MeV!!
    float GetTopPtWeight(float top_pt); // in MeV!!
    float GetTopPtAfterTtbarPtWeight(float top_pt); // in MeV!!
    float GetTtbarAndTopPtWeight(float ttbar_pt, float top_pt); // in MeV!!

  private:
    int m_sampleID;
    std::string m_weightsDirName;
    std::string m_weightsFileName;
    TFile *m_Weights_file;
    TH1* m_Hist_topPt;
    TH1* m_Hist_ttbarPt;
    TH1* m_Hist_topPtSeq;

};

#endif

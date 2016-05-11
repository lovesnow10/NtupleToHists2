#include "NNLOReweighter.hpp"
#include <iostream>

NNLOReweighter::NNLOReweighter(int sampleID, std::string m_weightsDirName) :
  m_sampleID(sampleID),
  m_weightsDirName(m_weightsDirName),
  m_weightsFileName(""),
  m_Weights_file(0),
  m_Hist_topPt(0),
  m_Hist_ttbarPt(0),
  m_Hist_topPtSeq(0)
{
}

NNLOReweighter::~NNLOReweighter(){
  m_Weights_file->Close();
  delete m_Weights_file;
}

void NNLOReweighter::SetInputDirectory(std::string dirName){
  m_weightsDirName = dirName;
}

void NNLOReweighter::SetInputFile(std::string fileName){
  m_weightsFileName = fileName;
}

void NNLOReweighter::SetSampleID(int sampleID){
  m_sampleID = sampleID;
}

void NNLOReweighter::Init(){
  // set the fileName according to sample ID
  if(m_sampleID!=0 && m_weightsFileName==""){
    if(m_sampleID==410000 || m_sampleID==410009 || m_sampleID==410120 || m_sampleID==410121) m_weightsFileName = m_weightsDirName+"/PowPyt6";
    if(m_sampleID==410004) m_weightsFileName = m_weightsDirName+"/PowHer";
    if(m_sampleID==410003) m_weightsFileName = m_weightsDirName+"/aMCHer";
    if(m_sampleID==410002) m_weightsFileName = m_weightsDirName+"/PowPytRadLo";
    if(m_sampleID==410001) m_weightsFileName = m_weightsDirName+"/PowPytRadHi";
    m_weightsFileName += "_TopPt_TTbarPt_TopPtSeq_rew.root";
  }
  if(m_weightsFileName==""){
    std::cout << "NNLOReweighter::WARNING: No input file nor valid sampleID specified. Not able to get reweighting." << std::endl;
    return;
  }
  //
  m_Weights_file = TFile::Open(m_weightsFileName.c_str());
  m_Hist_topPt    = (TH1*)m_Weights_file->Get("TopPt_rew"   );
  m_Hist_ttbarPt  = (TH1*)m_Weights_file->Get("TTbarPt_rew" );
  m_Hist_topPtSeq = (TH1*)m_Weights_file->Get("TopPtSeq_rew");
  m_Hist_topPt   ->SetDirectory(0);
  m_Hist_ttbarPt ->SetDirectory(0);
  m_Hist_topPtSeq->SetDirectory(0);
  m_Weights_file->Close();
}

float NNLOReweighter::GetTtbarPtWeight(float ttbar_pt){
  if(m_Hist_ttbarPt==0x0) return 1.;
  return m_Hist_ttbarPt ->GetBinContent( m_Hist_ttbarPt ->FindBin( ttbar_pt/1000. ) );
}

float NNLOReweighter::GetTopPtWeight(float top_pt){
  if(m_Hist_topPt==0x0) return 1.;
  return m_Hist_topPt->GetBinContent( m_Hist_topPt->FindBin( top_pt  /1000. ) );
}

float NNLOReweighter::GetTopPtAfterTtbarPtWeight(float top_pt){
  if(m_Hist_topPtSeq==0x0) return 1.;
  return m_Hist_topPtSeq->GetBinContent( m_Hist_topPtSeq->FindBin( top_pt  /1000. ) );
}

float NNLOReweighter::GetTtbarAndTopPtWeight(float ttbar_pt, float top_pt){
  if(m_Hist_ttbarPt==0x0 || m_Hist_topPtSeq==0x0) return 1.;
  return GetTtbarPtWeight(ttbar_pt)*GetTopPtAfterTtbarPtWeight(top_pt);
}

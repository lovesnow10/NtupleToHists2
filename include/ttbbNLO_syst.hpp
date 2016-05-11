#ifndef TTBBNLO_SYST_HPP_
#define TTBBNLO_SYST_HPP_

#include <map>
#include <string>
#include "TH2.h"
#include "TFile.h"

class HFSystDataMembers {

public:
  HFSystDataMembers();
  HFSystDataMembers(const HFSystDataMembers &dm);
  virtual ~HFSystDataMembers();

  // Declaration of needed variables as public data members
  int HF_Classification;

  float q1_pt;
  float q1_eta;
  float qq_pt;
  float qq_dr;
  float top_pt;
  float ttbar_pt;
};

class ttbbNLO_syst {

public:
  ttbbNLO_syst(std::string sampleID, std::string fileName,
               std::string fileName2);
  ~ttbbNLO_syst();
  void Init();

  std::map<std::string, float> GetttHFWeights(HFSystDataMembers *event);

  float ttbb_Nominal_rw();

  float ttbb_CSS_KIN_rw();
  float ttbb_MSTW_rw();
  float ttbb_NNPDF_rw();
  float ttbb_Q_CMMPS_rw();
  float ttbb_glosoft_rw();
  float ttbb_defaultX05_rw();
  float ttbb_defaultX2_rw();
  float ttbb_MPIup_rw();
  float ttbb_MPIdown_rw();

private:
  HFSystDataMembers *m_event;
  std::string m_sampleID;

  std::string weightsNormFileName;
  std::string weightsShapeFileName;
  TFile *Normsys_file;
  TFile *Shapesys_file;

  int m_HFClassification;
  float m_ttbar_pt;
  float m_top_pt;
  float m_var1;
  float m_var2;

  std::map<std::string, TH2F *> m_HistoMapForShapeRw;
  std::map<std::string, TH1D *> m_HistoMapForNormRw;

  std::map<std::string, float> m_MapForWeights;
};

#endif

//-------------------------------------//

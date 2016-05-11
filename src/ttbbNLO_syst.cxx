#include "ttbbNLO_syst.hpp"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

//-----------------------------------------------//

ttbbNLO_syst::ttbbNLO_syst(std::string sampleID, std::string fileName,
                           std::string fileName2)
    : m_event(0), m_sampleID(sampleID), weightsNormFileName(fileName),
      weightsShapeFileName(fileName2), Normsys_file(0), Shapesys_file(0),
      m_HFClassification(-999), m_ttbar_pt(-1.), m_top_pt(-1.), m_var1(-1.),
      m_var2(-1.) {

  m_HistoMapForShapeRw.clear();
  m_HistoMapForNormRw.clear();
}

ttbbNLO_syst::~ttbbNLO_syst() {

  Normsys_file->Close();
  Shapesys_file->Close();
  delete Normsys_file;
  delete Shapesys_file;
}

//--------//

std::map<std::string, float>
ttbbNLO_syst::GetttHFWeights(HFSystDataMembers *event) {

  m_event = event;

  // Retrieve HF info
  if (m_event->HF_Classification >= 100) { // tt + ≥1b
    m_HFClassification =
        (m_event->HF_Classification - (m_event->HF_Classification % 100)) / 100;
  } else if (m_event->HF_Classification > 0 &&
             m_event->HF_Classification < 100) { // tt + ≥1c
    m_HFClassification = -m_event->HF_Classification;
  } else if (m_event->HF_Classification == 0) { // tt light
    m_HFClassification = 0;
  } else if (m_event->HF_Classification == -1000) { // tt+b MPI
    m_HFClassification = -1000;
  } else { // other categories
    m_HFClassification = 0;
  }

  if (m_event->qq_pt != -99) {

    m_var1 = m_event->qq_pt / 1000.;
    m_var2 = m_event->qq_dr;

  } else if (m_event->qq_pt == -99) {
    m_var1 = m_event->q1_pt / 1000.;
    m_var2 = fabs(m_event->q1_eta);
  }

  m_top_pt = m_event->top_pt;
  m_ttbar_pt = m_event->ttbar_pt;

  m_MapForWeights.clear();

  if (m_HFClassification > 0) {

    m_MapForWeights["ttbb_Nominal_weight"] = ttbb_Nominal_rw();

    if (m_sampleID.find("410000") != std::string::npos ||
        m_sampleID.find("407009") != std::string::npos ||
        m_sampleID.find("407010") != std::string::npos ||
        m_sampleID.find("407011") != std::string::npos) {
      m_MapForWeights["ttbb_CSS_KIN_weight"] = ttbb_CSS_KIN_rw();
      m_MapForWeights["ttbb_MSTW_weight"] = ttbb_MSTW_rw();
      m_MapForWeights["ttbb_NNPDF_weight"] = ttbb_NNPDF_rw();
      m_MapForWeights["ttbb_Q_CMMPS_weight"] = ttbb_Q_CMMPS_rw();
      m_MapForWeights["ttbb_glosoft_weight"] = ttbb_glosoft_rw();
      m_MapForWeights["ttbb_defaultX05_weight"] = ttbb_defaultX05_rw();
      m_MapForWeights["ttbb_defaultX2_weight"] = ttbb_defaultX2_rw();
      m_MapForWeights["ttbb_MPIup_weight"] = ttbb_MPIup_rw();
      m_MapForWeights["ttbb_MPIdown_weight"] = ttbb_MPIdown_rw();
      m_MapForWeights["ttbb_MPIfactor_weight"] =
          m_MapForWeights["ttbb_Nominal_weight"];

    } else {
      m_MapForWeights["ttbb_CSS_KIN_weight"] = 1;
      m_MapForWeights["ttbb_MSTW_weight"] = 1;
      m_MapForWeights["ttbb_NNPDF_weight"] = 1;
      m_MapForWeights["ttbb_Q_CMMPS_weight"] = 1;
      m_MapForWeights["ttbb_glosoft_weight"] = 1;
      m_MapForWeights["ttbb_defaultX05_weight"] = 1;
      m_MapForWeights["ttbb_defaultX2_weight"] = 1;
      m_MapForWeights["ttbb_MPIup_weight"] = 1;
      m_MapForWeights["ttbb_MPIdown_weight"] = 1;
      m_MapForWeights["ttbb_MPIfactor_weight"] = 1;
    }
  } else {
    m_MapForWeights["ttbb_Nominal_weight"] = 1;
    m_MapForWeights["ttbb_CSS_KIN_weight"] = 1;
    m_MapForWeights["ttbb_MSTW_weight"] = 1;
    m_MapForWeights["ttbb_NNPDF_weight"] = 1;
    m_MapForWeights["ttbb_Q_CMMPS_weight"] = 1;
    m_MapForWeights["ttbb_glosoft_weight"] = 1;
    m_MapForWeights["ttbb_defaultX05_weight"] = 1;
    m_MapForWeights["ttbb_defaultX2_weight"] = 1;
    m_MapForWeights["ttbb_MPIup_weight"] = 1;
    m_MapForWeights["ttbb_MPIdown_weight"] = 1;
    if (m_HFClassification == -1000) {
      m_MapForWeights["ttbb_MPIfactor_weight"] = 1.5;
    } else {
      m_MapForWeights["ttbb_MPIfactor_weight"] = 1;
    }
  }
  return m_MapForWeights;
}

void ttbbNLO_syst::Init() {

  Shapesys_file = TFile::Open(weightsShapeFileName.c_str(), "READONLY");
  Normsys_file = TFile::Open(weightsNormFileName.c_str(), "READONLY");

  // Nominal
  m_HistoMapForShapeRw["Nominal_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Nominal_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Nominal_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Nominal_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Nominal_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["Nominal_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["Nominal_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["Nominal_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_Nominal_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["Nominal_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["Nominal_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["Nominal_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000");
  m_HistoMapForNormRw["Nominal_HFcateg"]->SetDirectory(0);

  // CSS_KIN
  m_HistoMapForShapeRw["CSS_KIN_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["CSS_KIN_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["CSS_KIN_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["CSS_KIN_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["CSS_KIN_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["CSS_KIN_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["CSS_KIN_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["CSS_KIN_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_CSSKIN_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["CSS_KIN_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["CSS_KIN_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["CSS_KIN_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_CSSKIN");
  m_HistoMapForNormRw["CSS_KIN_HFcateg"]->SetDirectory(0);

  // MSTW
  m_HistoMapForShapeRw["MSTW_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MSTW_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MSTW_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MSTW_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MSTW_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["MSTW_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["MSTW_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["MSTW_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MSTW_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["MSTW_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MSTW_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["MSTW_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_MSTW");
  m_HistoMapForNormRw["MSTW_HFcateg"]->SetDirectory(0);

  // NNPDF
  m_HistoMapForShapeRw["NNPDF_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["NNPDF_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["NNPDF_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["NNPDF_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["NNPDF_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["NNPDF_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["NNPDF_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["NNPDF_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_NNPDF_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["NNPDF_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["NNPDF_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["NNPDF_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_NNPDF");
  m_HistoMapForNormRw["NNPDF_HFcateg"]->SetDirectory(0);

  // Q_CMMPS
  m_HistoMapForShapeRw["Q_CMMPS_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Q_CMMPS_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Q_CMMPS_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Q_CMMPS_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["Q_CMMPS_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["Q_CMMPS_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["Q_CMMPS_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["Q_CMMPS_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_QCMMPS_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["Q_CMMPS_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["Q_CMMPS_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["Q_CMMPS_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_QCMMPS");
  m_HistoMapForNormRw["Q_CMMPS_HFcateg"]->SetDirectory(0);

  // glosoft
  m_HistoMapForShapeRw["glosoft_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["glosoft_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["glosoft_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["glosoft_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["glosoft_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["glosoft_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["glosoft_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["glosoft_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_glosoft_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["glosoft_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["glosoft_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["glosoft_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_glosoft");
  m_HistoMapForNormRw["glosoft_HFcateg"]->SetDirectory(0);

  // defaultX05
  m_HistoMapForShapeRw["defaultX05_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX05_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX05_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX05_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX05_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["defaultX05_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["defaultX05_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["defaultX05_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX05_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["defaultX05_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX05_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["defaultX05_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_X05");
  m_HistoMapForNormRw["defaultX05_HFcateg"]->SetDirectory(0);

  // defaultX2
  m_HistoMapForShapeRw["defaultX2_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX2_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX2_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX2_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["defaultX2_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["defaultX2_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["defaultX2_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["defaultX2_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_defaultX2_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["defaultX2_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["defaultX2_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["defaultX2_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_X2");
  m_HistoMapForNormRw["defaultX2_HFcateg"]->SetDirectory(0);

  // MPIup
  m_HistoMapForShapeRw["MPIup_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIup_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIup_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIup_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIup_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["MPIup_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["MPIup_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["MPIup_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIup_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["MPIup_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIup_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["MPIup_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_MPIup");
  m_HistoMapForNormRw["MPIup_HFcateg"]->SetDirectory(0);

  // MPIdown
  m_HistoMapForShapeRw["MPIdown_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIdown_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIdown_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIdown_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["MPIdown_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["MPIdown_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["MPIdown_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["MPIdown_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_MPIdown_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["MPIdown_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["MPIdown_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["MPIdown_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410000_MPIdown");
  m_HistoMapForNormRw["MPIdown_HFcateg"]->SetDirectory(0);

  // PowhegPythia_radHi
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radHi_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["PowhegPythia_radHi_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radHi_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["PowhegPythia_radHi_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radHi_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["PowhegPythia_radHi_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410001");
  m_HistoMapForNormRw["PowhegPythia_radHi_HFcateg"]->SetDirectory(0);

  // PowhegPythia_radLow
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radLow_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["PowhegPythia_radLow_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_radLow_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["PowhegPythia_radLow_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegPythia_radLow_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["PowhegPythia_radLow_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410002");
  m_HistoMapForNormRw["PowhegPythia_radLow_HFcateg"]->SetDirectory(0);

  // PowhegHerwigPP
  m_HistoMapForShapeRw["PowhegHerwigPP_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegHerwigPP_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegHerwigPP_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegHerwigPP_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["PowhegHerwigPP_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["PowhegHerwigPP_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["PowhegHerwigPP_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["PowhegHerwigPP_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_PowhegHerwigPP_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["PowhegHerwigPP_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["PowhegHerwigPP_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["PowhegHerwigPP_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410004");
  m_HistoMapForNormRw["PowhegHerwigPP_HFcateg"]->SetDirectory(0);

  // aMcAtNloHerwigPP
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_ttb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttB_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_ttB_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttbb_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_ttbb_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_tt3b_tt"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_tt3b_top_pt_ttbar_pt");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_ttb_q1_pt_q1_eta");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttB_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_ttB_q1_pt_q1_eta");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttbb_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_ttbb_qq_pt_qq_dr");
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_tt3b_bb"] =
      (TH2F *)Shapesys_file->Get("rwmap_aMcAtNloHerwigPP_tt3b_qq_pt_qq_dr");

  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttB_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttbb_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_tt3b_tt"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttB_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_ttbb_bb"]->SetDirectory(0);
  m_HistoMapForShapeRw["aMcAtNloHerwigPP_tt3b_bb"]->SetDirectory(0);

  m_HistoMapForNormRw["aMcAtNloHerwigPP_HFcateg"] =
      (TH1D *)Normsys_file->Get("h_410003");
  m_HistoMapForNormRw["aMcAtNloHerwigPP_HFcateg"]->SetDirectory(0);

  Shapesys_file->Close();
  Normsys_file->Close();
}

/// NOMINAL

float ttbbNLO_syst::ttbb_Nominal_rw() {

  std::string samplemap = "Nominal";

  if (m_sampleID.find("410001") != std::string::npos) {
    samplemap = "PowhegPythia_radHi";
  } else if (m_sampleID.find("410002") != std::string::npos) {
    samplemap = "PowhegPythia_radLow";
  } else if (m_sampleID.find("410003") != std::string::npos) {
    samplemap = "aMcAtNloHerwigPP";
  } else if (m_sampleID.find("410004") != std::string::npos) {
    samplemap = "PowhegHerwigPP";
  }

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw[samplemap + "_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw[samplemap + "_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_ttb_tt"]->FindBin(m_top_pt,
                                                             m_ttbar_pt));
    shape *= m_HistoMapForShapeRw[samplemap + "_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw[samplemap + "_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw[samplemap + "_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_ttB_tt"]->FindBin(m_top_pt,
                                                             m_ttbar_pt));
    shape *= m_HistoMapForShapeRw[samplemap + "_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw[samplemap + "_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw[samplemap + "_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_ttbb_tt"]->FindBin(m_top_pt,
                                                              m_ttbar_pt));
    shape *= m_HistoMapForShapeRw[samplemap + "_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw[samplemap + "_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw[samplemap + "_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_tt3b_tt"]->FindBin(m_top_pt,
                                                              m_ttbar_pt));
    shape *= m_HistoMapForShapeRw[samplemap + "_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw[samplemap + "_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------CSS_KIN

float ttbbNLO_syst::ttbb_CSS_KIN_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["CSS_KIN_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["CSS_KIN_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["CSS_KIN_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["CSS_KIN_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["CSS_KIN_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["CSS_KIN_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["CSS_KIN_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["CSS_KIN_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["CSS_KIN_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["CSS_KIN_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["CSS_KIN_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["CSS_KIN_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["CSS_KIN_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------MSTW

float ttbbNLO_syst::ttbb_MSTW_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["MSTW_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["MSTW_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MSTW_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["MSTW_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["MSTW_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MSTW_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["MSTW_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["MSTW_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MSTW_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["MSTW_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["MSTW_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MSTW_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MSTW_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------NNPDF

float ttbbNLO_syst::ttbb_NNPDF_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["NNPDF_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["NNPDF_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["NNPDF_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["NNPDF_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["NNPDF_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["NNPDF_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["NNPDF_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["NNPDF_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["NNPDF_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["NNPDF_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["NNPDF_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["NNPDF_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["NNPDF_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------Q_CMMPS

float ttbbNLO_syst::ttbb_Q_CMMPS_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["Q_CMMPS_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["Q_CMMPS_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["Q_CMMPS_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["Q_CMMPS_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["Q_CMMPS_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["Q_CMMPS_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["Q_CMMPS_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["Q_CMMPS_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["Q_CMMPS_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["Q_CMMPS_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["Q_CMMPS_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["Q_CMMPS_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["Q_CMMPS_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------glosoft

float ttbbNLO_syst::ttbb_glosoft_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["glosoft_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["glosoft_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["glosoft_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["glosoft_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["glosoft_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["glosoft_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["glosoft_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["glosoft_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["glosoft_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["glosoft_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["glosoft_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["glosoft_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["glosoft_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------defaultX05

float ttbbNLO_syst::ttbb_defaultX05_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["defaultX05_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["defaultX05_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_ttb_tt"]->FindBin(m_top_pt,
                                                           m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX05_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["defaultX05_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["defaultX05_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_ttB_tt"]->FindBin(m_top_pt,
                                                           m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX05_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["defaultX05_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["defaultX05_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_ttbb_tt"]->FindBin(m_top_pt,
                                                            m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX05_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["defaultX05_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["defaultX05_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_tt3b_tt"]->FindBin(m_top_pt,
                                                            m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX05_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX05_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------defaultX2

float ttbbNLO_syst::ttbb_defaultX2_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["defaultX2_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["defaultX2_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_ttb_tt"]->FindBin(m_top_pt,
                                                          m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX2_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["defaultX2_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["defaultX2_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_ttB_tt"]->FindBin(m_top_pt,
                                                          m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX2_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["defaultX2_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["defaultX2_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_ttbb_tt"]->FindBin(m_top_pt,
                                                           m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX2_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["defaultX2_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["defaultX2_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_tt3b_tt"]->FindBin(m_top_pt,
                                                           m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["defaultX2_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["defaultX2_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------MPIup

float ttbbNLO_syst::ttbb_MPIup_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["MPIup_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["MPIup_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIup_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["MPIup_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["MPIup_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIup_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["MPIup_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["MPIup_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIup_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["MPIup_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["MPIup_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIup_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIup_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//--------------------------------------------------MPIdown

float ttbbNLO_syst::ttbb_MPIdown_rw() {

  float norm = 0, shape = 0;

  switch (m_HFClassification) {
  case 10:
    norm = m_HistoMapForNormRw["MPIdown_HFcateg"]->GetBinContent(1);
    shape = m_HistoMapForShapeRw["MPIdown_ttb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_ttb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIdown_ttb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_ttb_bb"]->FindBin(m_var1, m_var2));
    break;
  case 1:
    norm = m_HistoMapForNormRw["MPIdown_HFcateg"]->GetBinContent(3);
    shape = m_HistoMapForShapeRw["MPIdown_ttB_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_ttB_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIdown_ttB_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_ttB_bb"]->FindBin(m_var1, m_var2));
    break;
  case 20:
    norm = m_HistoMapForNormRw["MPIdown_HFcateg"]->GetBinContent(2);
    shape = m_HistoMapForShapeRw["MPIdown_ttbb_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_ttbb_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIdown_ttbb_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_ttbb_bb"]->FindBin(m_var1, m_var2));
    break;
  default:
    norm = m_HistoMapForNormRw["MPIdown_HFcateg"]->GetBinContent(4);
    shape = m_HistoMapForShapeRw["MPIdown_tt3b_tt"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_tt3b_tt"]->FindBin(m_top_pt, m_ttbar_pt));
    shape *= m_HistoMapForShapeRw["MPIdown_tt3b_bb"]->GetBinContent(
        m_HistoMapForShapeRw["MPIdown_tt3b_bb"]->FindBin(m_var1, m_var2));
  }

  if (shape <= 0.)
    shape = 1.;
  return norm * shape;
}

//_________________________________________________________________________
//
HFSystDataMembers::HFSystDataMembers()
    : HF_Classification(-9), q1_pt(-99), q1_eta(-99), qq_pt(-99), qq_dr(-99),
      top_pt(-99), ttbar_pt(-99) {}

//_________________________________________________________________________
//
HFSystDataMembers::HFSystDataMembers(const HFSystDataMembers &dm) {

  HF_Classification = dm.HF_Classification;
  q1_pt = dm.q1_pt;
  q1_eta = dm.q1_eta;
  qq_pt = dm.qq_pt;
  qq_dr = dm.qq_dr;
  top_pt = dm.top_pt;
  ttbar_pt = dm.ttbar_pt;
}

//_________________________________________________________________________
//
HFSystDataMembers::~HFSystDataMembers() {}

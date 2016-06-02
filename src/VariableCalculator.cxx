#include "VariableCalculator.hpp"

float VariableCalculator::GetVarValue(string VarName) {
  if (mVariablesMap.find(VarName) == mVariablesMap.end()) {
    return -99999.0;
  } else {
    return mVariablesMap.at(VarName);
  }
}

void VariableCalculator::CalculateVariables(TTree *event) {
  mEvent = event;
  std::vector<float> jet_pt =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "jet_pt"));
  std::vector<float> jet_eta =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "jet_eta"));
  std::vector<float> jet_phi =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "jet_phi"));
  std::vector<float> jet_mv2c20 = *(
      Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "jet_mv2c20"));
  std::vector<float> el_pt =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "el_pt"));
  std::vector<float> el_eta =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "el_eta"));
  std::vector<float> el_phi =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "el_phi"));
  std::vector<float> mu_pt =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "mu_pt"));
  std::vector<float> mu_eta =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "mu_eta"));
  std::vector<float> mu_phi =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "mu_phi"));

  mVariablesMap["pT_jet1"] = jet_pt.at(0);
  mVariablesMap["pT_jet2"] = jet_pt.at(1);
  mVariablesMap["eta_jet1"] = jet_eta.at(0);
  mVariablesMap["eta_jet2"] = jet_eta.at(1);
  mVariablesMap["phi_jet1"] = jet_phi.at(0);
  mVariablesMap["phi_jet2"] = jet_phi.at(1);

  int nJets = *(Tools::Instance().GetTreeValue<int>(mEvent, "nJets"));
  for (int i = 0; i < nJets; i++) {
    float mv2c20 = jet_mv2c20.at(i);
    if (mv2c20 > -0.4434) {
      mVariablesMap["pT_bJet1"] = jet_pt.at(i);
      mVariablesMap["eta_bJet1"] = jet_eta.at(i);
      mVariablesMap["phi_bJet1"] = jet_phi.at(i);
      break;
    }
  }

  int nEl = *(Tools::Instance().GetTreeValue<int>(mEvent, "nElectrons"));
  int nMu = *(Tools::Instance().GetTreeValue<int>(mEvent, "nMuons"));
  if (nEl == 2 && nMu == 0) {
    mVariablesMap["pT_lep1"] = el_pt.at(0);
    mVariablesMap["pT_lep2"] = el_pt.at(1);
    mVariablesMap["eta_lep1"] = el_eta.at(0);
    mVariablesMap["eta_lep2"] = el_eta.at(1);
    mVariablesMap["phi_lep1"] = el_phi.at(0);
    mVariablesMap["phi_lep2"] = el_phi.at(1);
  } else if (nEl == 0 && nMu == 2) {
    mVariablesMap["pT_lep1"] = mu_pt.at(0);
    mVariablesMap["pT_lep2"] = mu_pt.at(1);
    mVariablesMap["eta_lep1"] = mu_eta.at(0);
    mVariablesMap["eta_lep2"] = mu_eta.at(1);
    mVariablesMap["phi_lep1"] = mu_phi.at(0);
    mVariablesMap["phi_lep2"] = mu_phi.at(1);
  } else if (nEl == 1 && nMu == 1) {
    float ept = el_pt.at(0);
    float mpt = mu_pt.at(0);
    if (ept > mpt) {
      mVariablesMap["pT_lep1"] = el_pt.at(0);
      mVariablesMap["pT_lep2"] = mu_pt.at(0);
      mVariablesMap["eta_lep1"] = el_eta.at(0);
      mVariablesMap["eta_lep2"] = mu_eta.at(0);
      mVariablesMap["phi_lep1"] = el_phi.at(0);
      mVariablesMap["phi_lep2"] = mu_eta.at(0);
    } else {
      mVariablesMap["pT_lep1"] = mu_pt.at(0);
      mVariablesMap["pT_lep2"] = el_pt.at(0);
      mVariablesMap["eta_lep1"] = mu_eta.at(0);
      mVariablesMap["eta_lep2"] = el_eta.at(0);
      mVariablesMap["phi_lep1"] = mu_phi.at(0);
      mVariablesMap["phi_lep2"] = el_phi.at(0);
    }
  }
}

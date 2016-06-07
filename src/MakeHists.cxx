#include "MakeHists.hpp"

using namespace std;

bool MakeHists::initialize(ConfigParser *config, DSHandler *ds) {
  mConfig = config;
  hs = new HistStore();
  calculator = new VariableCalculator();
  InitYields(ds);
  mTRFvariables.push_back("Mbb_MindR");
  mTRFvariables.push_back("dRbb_MaxPt");
  mTRFvariables.push_back("dRbb_MaxM");

  mVarToCalc.push_back("pT_jet1");
  mVarToCalc.push_back("pT_jet2");
  mVarToCalc.push_back("eta_jet1");
  mVarToCalc.push_back("eta_jet2");
  mVarToCalc.push_back("phi_jet1");
  mVarToCalc.push_back("phi_jet2");
  mVarToCalc.push_back("pT_bJet1");
  mVarToCalc.push_back("eta_bJet1");
  mVarToCalc.push_back("phi_bJet1");
  mVarToCalc.push_back("pT_lep1");
  mVarToCalc.push_back("pT_lep2");
  mVarToCalc.push_back("eta_lep1");
  mVarToCalc.push_back("eta_lep2");
  mVarToCalc.push_back("phi_lep1");
  mVarToCalc.push_back("phi_lep2");
  return true;
}

bool MakeHists::run(TTree *event, map<string, float> weights,
                    TTreeFormulaContainer *formulas,
                    std::map<string, bool> bControl) {
  mEvent = event;
  std::vector<string> mRegions;
  mRegions.clear();
  string mSample;
  bool isTRF = bControl.at("isTRF");
  bool doTTbarMerge = bControl.at("doTTbarMerge");
  //  bool doHFSplit = bControl.at("doHFSplit");
  int mcChannel =
      *(Tools::Instance().GetTreeValue<int>(mEvent, "mcChannelNumber"));
  if (mcChannel == 0)
    isTRF = false;

  // Fakes!
  bool isFake = false;
  mSample = Tools::Instance().GetSampleType(mcChannel);
  if (mcChannel == 0)
    isFake = false;
  if (mSample == "Fakes")
    isFake = true;

  // Get Weights!
  float mWeights = Tools::Instance().GetWeight(mcChannel);

  // No weight for DATA!
  if (mcChannel != 0) {
    if (mSample == "ttbar") {
      int HF_Classification =
          *(Tools::Instance().GetTreeValue<int>(mEvent, "HF_Classification"));
      if (HF_Classification == 0)
        mSample = "ttlight";
      else if (abs(HF_Classification) > 0 && abs(HF_Classification) < 100)
        mSample = "ttcc";
      else
        mSample = "ttbb";
    }

    float weight_mc =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_mc"));
    float weight_pileup =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_pileup"));
    float weight_jvt =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_jvt"));
    float weight_leptonSF =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_leptonSF"));
    float weight_bTagSF_77 =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_bTagSF_77"));
    if (mSample != "ttbb") {
      if (!isTRF)
        mWeights = mWeights * weight_mc * weight_pileup * weight_jvt *
                   weight_leptonSF * weight_bTagSF_77 *
                   weights["weight_ttbb_Nominal"] * weights["weight_NNLO"];
      else
        mWeights = mWeights * weight_mc * weight_pileup * weight_jvt *
                   weight_leptonSF * weights["weight_ttbb_Nominal"] *
                   weights["weight_NNLO"];
    } else {
      if (!isTRF)
        mWeights = mWeights * weight_mc * weight_pileup * weight_jvt *
                   weight_leptonSF * weight_bTagSF_77 *
                   weights["weight_ttbb_Nominal"];
      else
        mWeights = mWeights * weight_mc * weight_pileup * weight_jvt *
                   weight_leptonSF * weights["weight_ttbb_Nominal"];
    }
  }

  // Selectioins!
  for (int ientry = 0; ientry < formulas->GetEntries(); ++ientry) {
    TTreeFormula *formula = formulas->GetFormula(ientry);
    if (formula->EvalInstance()) {
      mRegions.push_back(formula->GetName());
    }
  }
  if (mRegions.empty())
    return false;

  std::map<string, float> mTRFweights;
  if (isTRF) {
    mTRFweights["2bex"] =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "trf_weight_77_2ex"));
    mTRFweights["3bex"] =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "trf_weight_77_3ex"));
    mTRFweights["4bin"] =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "trf_weight_77_4in"));
    std::vector<float> tmpTRFweights =
        *(Tools::Instance().GetTreeValue<vector<float>>(mEvent,
                                                        "trf_weight_77_in"));
    mTRFweights["2bin"] = tmpTRFweights[2];
    mTRFweights["3bin"] = tmpTRFweights[3];
  }
  // Fake leptons removal
  int nEl = *(Tools::Instance().GetTreeValue<int>(mEvent, "nElectrons"));
  int nMu = *(Tools::Instance().GetTreeValue<int>(mEvent, "nMuons"));
  if (mcChannel != 0 && mSample != "Fakes") {
    std::vector<int> el_true_type =
        *(Tools::Instance().GetTreeValue<std::vector<int>>(mEvent,
                                                           "el_true_type"));
    std::vector<int> el_true_origin =
        *(Tools::Instance().GetTreeValue<std::vector<int>>(mEvent,
                                                           "el_true_origin"));
    std::vector<int> mu_true_type =
        *(Tools::Instance().GetTreeValue<std::vector<int>>(mEvent,
                                                           "mu_true_type"));
    std::vector<int> mu_true_origin =
        *(Tools::Instance().GetTreeValue<std::vector<int>>(mEvent,
                                                           "mu_true_origin"));
    std::vector<int> el_true_originbkg =
        *(Tools::Instance().GetTreeValue<std::vector<int>>(
            mEvent, "el_true_originbkg"));
    int et1 = -1;
    int et2 = -1;
    int eo1 = -1;
    int eo2 = -1;
    int eb1 = -1;
    int eb2 = -1;

    if (nEl == 0 && nMu == 2) {
      et1 = mu_true_type.at(0);
      et2 = mu_true_type.at(1);
      eo1 = mu_true_origin.at(0);
      eo2 = mu_true_origin.at(1);
      if (!(et1 == 6 && et2 == 6 &&
            (eo1 == 10 || eo1 == 12 || eo1 == 13 || eo1 == 14) &&
            (eo2 == 10 || eo2 == 12 || eo2 == 13 || eo2 == 14)))
        mSample = "Fakes";
    } else if (nEl == 2 && nMu == 0) {
      et1 = el_true_type.at(0);
      et2 = el_true_type.at(1);
      eo1 = el_true_origin.at(0);
      eo2 = el_true_origin.at(1);
      eb1 = el_true_originbkg.at(0);
      eb2 = el_true_originbkg.at(1);

      bool pass =
          ((et1 == 2 && (eo1 == 10 || eo1 == 12 || eo1 == 13 || eo1 == 14)) ||
           (et1 == 4 && eo1 == 5 &&
            (eb1 == 10 || eb1 == 12 || eb1 == 13 || eb1 == 14)));
      pass =
          pass &&
          ((et2 == 2 && (eo2 == 10 || eo2 == 12 || eo2 == 13 || eo2 == 14)) ||
           (et2 == 4 && eo2 == 5 &&
            (eb2 == 10 || eb2 == 12 || eb2 == 13 || eb2 == 14)));
      if (!pass)
        mSample = "Fakes";
    } else if (nEl == 1 && nMu == 1) {
      et1 = el_true_type.at(0);
      et2 = mu_true_type.at(0);
      eo1 = el_true_origin.at(0);
      eo2 = mu_true_origin.at(0);
      eb1 = el_true_originbkg.at(0);

      bool pass =
          et2 == 6 && (eo2 == 10 || eo2 == 12 || eo2 == 13 || eo2 == 14);
      pass =
          pass &&
          ((et1 == 2 && (eo1 == 10 || eo1 == 12 || eo1 == 13 || eo1 == 14)) ||
           (et1 == 4 && eo1 == 5 &&
            (eb1 == 10 || eb1 == 12 || eb1 == 13 || eb1 == 14)));
      if (!pass)
        mSample = "Fakes";
    }
    // ttbar merge
    if (doTTbarMerge) {
      int TopHeavyFlavorFilterFlag = *(Tools::Instance().GetTreeValue<int>(
          mEvent, "TopHeavyFlavorFilterFlag"));
      bool truth_top_dilep_filter = *(Tools::Instance().GetTreeValue<bool>(
          mEvent, "truth_top_dilep_filter"));

      if (mcChannel == 410000 &&
          (TopHeavyFlavorFilterFlag == 5 || truth_top_dilep_filter == true))
        return false;
      if (mcChannel == 410120 && truth_top_dilep_filter == true)
        return false;
      if (mcChannel == 410009 && TopHeavyFlavorFilterFlag == 5)
        return false;
    }
  }
  // Calculate Variables
  /*  int nJets = *(Tools::Instance().GetTreeValue<int>(mEvent, "nJets"));
    //  int nBTags = mEvent->GetValue<int>("nBTags");
    if (mSample == "ttbar") {
      int HF_Classification =
          *(Tools::Instance().GetTreeValue<int>(mEvent, "HF_Classification"));
      if (HF_Classification == 0)
        mSample = "ttlight";
      else if (abs(HF_Classification) > 0 && abs(HF_Classification) < 100)
        mSample = "ttcc";
      else
        mSample = "ttbb";
    }

    float pT_jet1, pT_jet2, pT_bJet1, pT_lep1, pT_lep2;
    float eta_jet1, eta_jet2, eta_bJet1, eta_lep1, eta_lep2;
    float phi_jet1, phi_jet2, phi_bJet1, phi_lep1, phi_lep2;

    std::vector<float> jet_pt =
        *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent, "jet_pt"));
    std::vector<float> jet_eta =
        *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent,
    "jet_eta"));
    std::vector<float> jet_phi =
        *(Tools::Instance().GetTreeValue<std::vector<float>>(mEvent,
    "jet_phi"));
    std::vector<float> jet_mv2c20 = *(
        Tools::Instance().GetTreeValue<std::vector<float>>(mEvent,
    "jet_mv2c20"));
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

    pT_jet1 = jet_pt.at(0);
    pT_jet2 = jet_pt.at(1);
    eta_jet1 = jet_eta.at(0);
    eta_jet2 = jet_eta.at(1);
    phi_jet1 = jet_phi.at(0);
    phi_jet2 = jet_phi.at(1);
    for (int i = 0; i < nJets; i++) {
      float mv2c20 = jet_mv2c20.at(i);
      if (mv2c20 > -0.4434) {
        pT_bJet1 = jet_pt.at(i);
        eta_bJet1 = jet_eta.at(i);
        phi_bJet1 = jet_phi.at(i);
        break;
      }
    }
    // Setup leptons
    if (nEl == 2 && nMu == 0) {
      pT_lep1 = el_pt.at(0);
      pT_lep2 = el_pt.at(1);
      eta_lep1 = el_eta.at(0);
      eta_lep2 = el_eta.at(1);
      phi_lep1 = el_phi.at(0);
      phi_lep2 = el_phi.at(1);
    } else if (nEl == 0 && nMu == 2) {
      pT_lep1 = mu_pt.at(0);
      pT_lep2 = mu_pt.at(1);
      eta_lep1 = mu_eta.at(0);
      eta_lep2 = mu_eta.at(1);
      phi_lep1 = mu_phi.at(0);
      phi_lep2 = mu_phi.at(1);
    } else if (nEl == 1 && nMu == 1) {
      float ept = el_pt.at(0);
      float mpt = mu_pt.at(0);
      if (ept > mpt) {
        pT_lep1 = el_pt.at(0);
        pT_lep2 = mu_pt.at(0);
        eta_lep1 = el_eta.at(0);
        eta_lep2 = mu_eta.at(0);
        phi_lep1 = el_phi.at(0);
        phi_lep2 = mu_eta.at(0);
      } else {
        pT_lep1 = mu_pt.at(0);
        pT_lep2 = el_pt.at(0);
        eta_lep1 = mu_eta.at(0);
        eta_lep2 = el_eta.at(0);
        phi_lep1 = mu_phi.at(0);
        phi_lep2 = el_phi.at(0);
      }
    }*/
  calculator->CalculateVariables(mEvent);
  for (auto region : mRegions) {
    if (isTRF) {
      if (region.find("2bex") != string::npos) {
        mWeights = mWeights * mTRFweights["2bex"];
      }
      if (region.find("3bex") != string::npos) {
        mWeights = mWeights * mTRFweights["3bex"];
      }
      if (region.find("4bin") != string::npos) {
        mWeights = mWeights * mTRFweights["4bin"];
      }
    }
    mRawYields.at(region).at(mSample) += 1;
    mWeightedYields.at(region).at(mSample) += (weights["norm"] * mWeights);
    std::vector<string> vars = mConfig->GetRegionVars(region);
    for (auto var : vars) {
      float value;
      /*      if (var == "pT_jet1")
              value = pT_jet1;
            else if (var == "pT_jet2")
              value = pT_jet2;
            else if (var == "eta_jet1")
              value = eta_jet1;
            else if (var == "eta_jet2")
              value = eta_jet2;
            else if (var == "phi_jet1")
              value = phi_jet1;
            else if (var == "phi_jet2")
              value = phi_jet2;
            else if (var == "pT_bJet1")
              value = pT_bJet1;
            else if (var == "eta_bJet1")
              value = eta_bJet1;
            else if (var == "phi_bJet1")
              value = phi_bJet1;
            else if (var == "pT_lep1")
              value = pT_lep1;
            else if (var == "pT_lep2")
              value = pT_lep2;
            else if (var == "eta_lep1")
              value = eta_lep1;
            else if (var == "eta_lep2")
              value = eta_lep2;
            else if (var == "phi_lep1")
              value = phi_lep1;
            else if (var == "phi_lep2")
              value = phi_lep2;*/
      if (find(mVarToCalc.begin(), mVarToCalc.end(), var) != mVarToCalc.end()) {
        value = calculator->GetVarValue(var);
      } else {
        string valueType = Tools::Instance().GetValueType(mEvent, var);
        if (!isTRF) {
          if (valueType == "int" || valueType == "Int_t") {
            value = *(Tools::Instance().GetTreeValue<int>(mEvent, var));
          } else {
            value = *(Tools::Instance().GetTreeValue<float>(mEvent, var));
          }
        } else {
          if (valueType == "int" || valueType == "Int_t") {
            value = *(Tools::Instance().GetTreeValue<int>(mEvent, var));
          } else {
            if (find(mTRFvariables.begin(), mTRFvariables.end(), var) !=
                mTRFvariables.end()) {
              if (region.find("2bex") != string::npos) {
                string tmpVar = var + "_77_2ex";
                value =
                    *(Tools::Instance().GetTreeValue<float>(mEvent, tmpVar));
              } else if (region.find("3bex") != string::npos) {
                string tmpVar = var + "_77_3ex";
                value =
                    *(Tools::Instance().GetTreeValue<float>(mEvent, tmpVar));
              } else if (region.find("4bin") != string::npos) {
                string tmpVar = var + "_77_4in";
                value =
                    *(Tools::Instance().GetTreeValue<float>(mEvent, tmpVar));
              } else {
                value = *(Tools::Instance().GetTreeValue<float>(mEvent, var));
              }
            } else {
              value = *(Tools::Instance().GetTreeValue<float>(mEvent, var));
            }
          }
        }
      }
      string hname = Tools::Instance().GenName(var, region, mSample);
      if (hs->HasHist(hname)) {
        hs->GetHist(hname)->Fill(value, mWeights * weights["norm"]);
      } else {
        std::vector<float> bins = mConfig->GetVarBinning(region, var);
        hs->AddHist(hname, bins[0], bins[1], bins[2]);
        hs->GetHist(hname)->Fill(value, mWeights * weights["norm"]);
      }
    }
  }
  return true;
}

bool MakeHists::finalize(TFile *fFile) {
  FillYields();
  hs->SaveAllHists(fFile);
  printf("MakeHists:: finalize:: MakeHists has finished running!\n");
  return true;
}

void MakeHists::InitYields(DSHandler *ds) {
  std::vector<string> regions = mConfig->GetRegions();
  std::vector<string> samples = ds->GetAllTypes();
  int nSamples = samples.size();
  int nRegions = regions.size();
  if (find(samples.begin(), samples.end(), "ttbar") == samples.end()) {
    hs->AddHist2D("hist_raw_yields", nSamples, 0, nSamples, nRegions, 0,
                  nRegions);
    hs->AddHist2D("hist_weighted_yields", nSamples, 0, nSamples, nRegions, 0,
                  nRegions);
  } else {
    hs->AddHist2D("hist_raw_yields", nSamples + 2, 0, nSamples + 2, nRegions, 0,
                  nRegions);
    hs->AddHist2D("hist_weighted_yields", nSamples + 2, 0, nSamples + 2,
                  nRegions, 0, nRegions);
  }
  for (auto region : regions) {
    for (auto sample : samples) {
      if (sample != "ttbar") {
        mRawYields[region][sample] = 0;
        mWeightedYields[region][sample] = 0;
      } else {
        mRawYields[region]["ttlight"] = 0;
        mWeightedYields[region]["ttlight"] = 0;
        mRawYields[region]["ttcc"] = 0;
        mWeightedYields[region]["ttcc"] = 0;
        mRawYields[region]["ttbb"] = 0;
        mWeightedYields[region]["ttbb"] = 0;
      }
    }
  }
  int nx = 1, ny = 1;
  for (auto region : regions) {
    hs->GetHist2D("hist_raw_yields")
        ->GetYaxis()
        ->SetBinLabel(ny, region.c_str());
    hs->GetHist2D("hist_weighted_yields")
        ->GetYaxis()
        ->SetBinLabel(ny++, region.c_str());
  }
  for (auto sample : samples) {
    if (sample != "ttbar") {
      hs->GetHist2D("hist_raw_yields")
          ->GetXaxis()
          ->SetBinLabel(nx, sample.c_str());
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, sample.c_str());
    } else {
      hs->GetHist2D("hist_raw_yields")->GetXaxis()->SetBinLabel(nx, "ttlight");
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, "ttlight");
      hs->GetHist2D("hist_raw_yields")->GetXaxis()->SetBinLabel(nx, "ttcc");
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, "ttcc");
      hs->GetHist2D("hist_raw_yields")->GetXaxis()->SetBinLabel(nx, "ttbb");
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, "ttbb");
    }
  }
}

void MakeHists::FillYields() {
  int nx = hs->GetHist2D("hist_raw_yields")->GetNbinsX();
  int ny = hs->GetHist2D("hist_raw_yields")->GetNbinsY();
  for (int ix = 1; ix < nx + 1; ix++) {
    string xname =
        hs->GetHist2D("hist_raw_yields")->GetXaxis()->GetBinLabel(ix);
    for (int iy = 1; iy < ny + 1; iy++) {
      string yname =
          hs->GetHist2D("hist_raw_yields")->GetYaxis()->GetBinLabel(iy);
      string tmpName = yname + "_" + xname;
      printf("HistsGen:: FillYields:: Filling Yields %s\n", tmpName.c_str());
      float raw, weighted;
      raw = mRawYields.at(yname).at(xname);
      weighted = mWeightedYields.at(yname).at(xname);
      hs->GetHist2D("hist_raw_yields")->SetBinContent(ix, iy, raw);
      hs->GetHist2D("hist_weighted_yields")->SetBinContent(ix, iy, weighted);
    }
  }
}

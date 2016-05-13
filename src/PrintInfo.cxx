#include "PrintInfo.hpp"

bool PrintInfo::initialize() {
  time_t nowtime;
  nowtime = time(NULL);
  struct tm *local = localtime(&nowtime);
  char buf[80];
  strftime(buf, 80, "%Y%m%d_%H%M", local);
  string tmpName("Info");
  tmpName.append(buf);
  tmpName.append(".log");
  mOut = ofstream(tmpName.c_str());
  return true;
}

bool PrintInfo::run(TTree *event, std::map<string, float> weights,
                    long ientry) {
  mEvent = event;
  float norm = weights["norm"];
  float weight_ttbb_Nominal = weights["weight_ttbb_Nominal"];
  float weight_NNLO = weights["weight_NNLO"];

  int mcChannel =
      *(Tools::Instance().GetTreeValue<int>(mEvent, "mcChannelNumber"));
  int runNumber = *(Tools::Instance().GetTreeValue<int>(mEvent, "runNumber"));

  float weight_mc = 1.0;
  float weight_pileup = 1.0;
  float weight_leptonSF = 1.0;
  float weight_bTagSF_77 = 1.0;
  float weight_jvt = 1.0;
  //float truth_top_pt = 1.0;
  //float truth_ttbar_pt = 1.0;

  if (mcChannel != 0) {
    weight_mc = *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_mc"));
    weight_pileup =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_pileup"));
    weight_leptonSF =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_leptonSF"));
    weight_bTagSF_77 =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_bTagSF_77"));
    weight_jvt = *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_jvt"));

    //truth_top_pt =
        //*(Tools::Instance().GetTreeValue<float>(mEvent, "truth_top_pt"));
    //truth_ttbar_pt =
        //*(Tools::Instance().GetTreeValue<float>(mEvent, "truth_ttbar_pt"));
  }

  float mWeights = weight_mc * weight_pileup * weight_leptonSF * weight_bTagSF_77 *weight_jvt * weight_ttbb_Nominal * weight_NNLO * norm;

  mOut << setfill('=') << setw(80) << "" << endl;
  mOut << endl;
  mOut << setfill(' ') << setw(15)
       << "DSID: " << (mcChannel == 0 ? runNumber : mcChannel);
  mOut << setfill(' ') << setw(10) << "Entry:" << ientry << endl;
  mOut << endl;
  mOut << setfill('*') << setw(80) << "" << endl;
  mOut << setfill(' ') << setw(10) << "Normalization: " << norm << endl;
  mOut << setfill(' ') << setw(10) << "weight_mc: " << weight_mc << endl;
  mOut << setfill(' ') << setw(10) << "weight_pileup: " << weight_pileup
       << endl;
  mOut << setfill(' ') << setw(10) << "weight_jvt: " << weight_jvt << endl;
  mOut << setfill(' ') << setw(10) << "weight_leptonSF: " << weight_leptonSF
       << endl;
  mOut << setfill(' ') << setw(10) << "weight_bTagSF_77: " << weight_bTagSF_77
       << endl;
  mOut << setfill(' ') << setw(10)
       << "weight_ttbb_Nominal: " << weight_ttbb_Nominal << endl;
  mOut << setfill(' ') << setw(10) << "weight_NNLO: " << weight_NNLO << endl;
  mOut << setfill(' ') << setw(10) << "Weights: " << mWeights << endl;

  mOut << setfill('=') << setw(80) << "" << endl;
  mOut << endl;

  return true;
}

bool PrintInfo::finalize() {
  mOut.close();
  return true;
}

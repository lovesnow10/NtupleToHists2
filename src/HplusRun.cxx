#include "HplusRun.hpp"
#include "ttbbNLO_syst.hpp"
#include "NNLOReweighter.hpp"
#include "TROOT.h"

using namespace std;

HplusRun::HplusRun(int argc, char **argv) : mArgc(argc), mArgv(argv) {}

bool HplusRun::initialize() {
  if (mArgc < 3) {
    printf("HistGen:: initialize:: Given args %i Are Not Enough!", mArgc);
    return false;
  }
  mWorker = new TChain("nominal");
  mHelpWorker = new TChain("sumWeights");
  mConfig = new ConfigParser(mArgv[1]);
  mDS = new DSHandler(mArgv[2]);
  string XSFile = mConfig->GetCommonSetting("XSFile")[0];
  mXSHelper = new XSHelper(XSFile);

  if (!(mConfig->GetCommonSetting("MAX").empty())) mMaxProcess = atol(mConfig->GetCommonSetting("MAX")[0].c_str());
  else mMaxProcess = 0;
  mProcessed = 0;

  // initialize working classes
  mMH = new MakeHists();
  mMH->initialize(mConfig, mDS);

  mPI = new PrintInfo();
  mPI->initialize();

  printf("HplusRun:: initialize:: HistGen Initialization Succeeded!\n");
  return true;
}

bool HplusRun::run() {
  std::vector<string> files = mDS->Next();
  while (!files.empty()) {
    // Init general objs
    mWorker->Reset();
    mHelpWorker->Reset();
    for (auto file : files) {
      mWorker->Add(file.c_str());
      mHelpWorker->Add(file.c_str());
    }
    // Set up formulas
    std::vector<string> mRegions = mConfig->GetRegions();
    mFormulas = new TTreeFormulaContainer();
    for (auto region : mRegions) {
      string cut = mConfig->GetRegionCut(region);
      TTreeFormula *formula =
          new TTreeFormula(region.c_str(), cut.c_str(), mWorker);
      mFormulas->AddFormula(formula);
    }
    mWorker->SetNotify(mFormulas);

    long nentries = mWorker->GetEntries();
    if (nentries == 0)
      continue;
    long messageSlice = long(nentries / 6);
    if (messageSlice == 0)
      messageSlice++;
    mWorker->GetEntry(0);
    int runNumber =
        *(Tools::Instance().GetTreeValue<int>(mWorker, "runNumber"));
    int mcChannel =
        *(Tools::Instance().GetTreeValue<int>(mWorker, "mcChannelNumber"));
    // Calculate norm and init ttbbRW
    ttbbNLO_syst *ttbbRW;
    NNLOReweighter *nnloRW;
    float norm = 1.0;
    if (mcChannel != 0) {
      float totalEventsWeighted =
          Tools::Instance().GetSumUp<float>(mHelpWorker, "totalEventsWeighted");
      float xs = mXSHelper->GetXS(to_string(mcChannel));
      float lumi = atof((mConfig->GetCommonSetting("Luminosity")[0]).c_str());
      norm = xs * lumi / totalEventsWeighted;

      string NormFile = mConfig->GetCommonSetting("NormFile")[0];
      string ShapeFile = mConfig->GetCommonSetting("ShapeFile")[0];
      string NNLODir = mConfig->GetCommonSetting("NNLODir")[0];
      if (mcChannel == 410000 || mcChannel == 410009 || mcChannel == 410120 ||
          mcChannel == 410121) {
        ttbbRW = new ttbbNLO_syst("410000", NormFile, ShapeFile);
      } else {
        ttbbRW = new ttbbNLO_syst(to_string(mcChannel), NormFile, ShapeFile);
      }
      ttbbRW->Init();
      if (mcChannel == 410000 || mcChannel == 410009 || mcChannel == 410120 || mcChannel == 410121){
      nnloRW = new NNLOReweighter(mcChannel, NNLODir);
      nnloRW->Init();}
    }

    // Main loop!
    printf("Start Loop %i\n", mcChannel != 0 ? mcChannel : runNumber);
    for (long ientry = 0; ientry < nentries; ++ientry) {
      if ((mProcessed >= mMaxProcess) && (mMaxProcess > 0))
      {
        printf("HplusRun:: Run:: Reach Max Process Number %ld\n", mMaxProcess);
        break;
      }
      mWorker->GetEntry(ientry);
      if (ientry % messageSlice == 0) {
        printf("HplusRun:: Run:: ----------Now Processing %ld of %ld in "
               "%i----------\n",
               ientry + 1, nentries, mcChannel != 0 ? mcChannel : runNumber);
      }
      // Get ttbb weight
      float weight_ttbb_Nominal = 1.0;
      float weight_NNLO = 1.0;
      if (mcChannel != 0) {
        HFSystDataMembers *ttbb = new HFSystDataMembers();
        ttbb->HF_Classification = *(
            Tools::Instance().GetTreeValue<int>(mWorker, "HF_Classification"));
        ttbb->q1_eta =
            *(Tools::Instance().GetTreeValue<float>(mWorker, "q1_eta"));
        ttbb->q1_pt =
            *(Tools::Instance().GetTreeValue<float>(mWorker, "q1_pt"));
        ttbb->qq_dr =
            *(Tools::Instance().GetTreeValue<float>(mWorker, "qq_dr"));
        ttbb->qq_pt =
            *(Tools::Instance().GetTreeValue<float>(mWorker, "qq_pt"));
        ttbb->top_pt =
            *(Tools::Instance().GetTreeValue<float>(mWorker, "truth_top_pt")) *
            1e-3;
        ttbb->ttbar_pt = *(Tools::Instance().GetTreeValue<float>(
                             mWorker, "truth_ttbar_pt")) *
                         1e-3;

        weight_ttbb_Nominal =
            ttbbRW->GetttHFWeights(ttbb).at("ttbb_Nominal_weight");
        if (mcChannel == 410000 || mcChannel == 410009 || mcChannel == 410120 || mcChannel == 410121){
        float truth_top_pt = *(Tools::Instance().GetTreeValue<float>(mWorker, "truth_top_pt"));
        float truth_ttbar_pt = *(Tools::Instance().GetTreeValue<float>(mWorker, "truth_ttbar_pt"));
        weight_NNLO = nnloRW->GetTtbarAndTopPtWeight(truth_ttbar_pt, truth_top_pt);}

      }
      std::map<string, float> weights;
      weights["norm"] = norm;
      weights["weight_NNLO"] = weight_NNLO;
      weights["weight_ttbb_Nominal"] = weight_ttbb_Nominal;
      mMH->run(mWorker, weights, mFormulas);
      mPI->run(mWorker, weights, ientry);
      mProcessed++;
    }
    // Point to next DataSet!
    files = mDS->Next();
  }
  return true;
}

bool HplusRun::finalize() {
  string mOutName = mConfig->GetCommonSetting("OutName")[0];
  fFile = TFile::Open(mOutName.c_str(), "RECREATE");

  // finalize working classes
  mMH->finalize(fFile);
  mPI->finalize();

  // Close File!
  fFile->Close();
  printf("HplusRun:: finalize:: HplusRun Has Finished Running!\n");
  return true;
}

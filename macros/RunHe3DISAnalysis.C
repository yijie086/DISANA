#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "./../DreamAN/Cuts/QADBCuts.h"
#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/EventProcessor.h"
#include "./../DreamAN/core/He3DISAnalysis.h"

void RunHe3DISAnalysis(const std::string& inputDir, int nfile,
                       int nthreads = 0) {
  if (nthreads <= 0) {
    const auto hc = std::thread::hardware_concurrency();
    nthreads = hc ? static_cast<int>(hc) : 2;
  }

  if (nthreads > 1) {
    ROOT::EnableImplicitMT(nthreads);
    std::cout << "[RunHe3DISAnalysis] IMT enabled with " << nthreads
              << " thread(s)\n";
  } else {
    std::cout << "[RunHe3DISAnalysis] IMT disabled (single thread mode).\n";
  }

  bool IsMC = true;
  bool IsreprocRootFile = false;
  bool IsMinimalBook = true;
  bool FiducalSystematics = false;
  bool PIDSystematics = false;
  float SFthreeSigmaFactor = 5.0f / 6.0f;
  std::cout << "Running He3 DIS Analysis: epp OR epn OR ed" << std::endl;

  std::string inputFileDir = inputDir;
  std::string inputRootFileName = " ";
  std::string inputRootTreeName = " ";
  const std::string inputHipoDirTest =
      "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/"
      "NewAnalysisFrameWork/temp_hipofile/";

  if (inputDir.empty()) {
    inputFileDir = inputHipoDirTest;
    std::cout << "No custom input directory provided; using "
              << inputFileDir << std::endl;
  }
  if (IsreprocRootFile) {
    inputRootFileName = "dfSelected.root";
    inputRootTreeName = "dfSelected";
  }

  AnalysisTaskManager mgr;
  const std::string outputFileDir = "./";
  mgr.SetOututDir(outputFileDir);

  auto trackCuts = std::make_shared<TrackCut>();
  auto edgeRegionsElectron = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edgeRegionsProton = std::vector<float>{3.0f, 3.0f, 5.0f};
  auto edgeRegionsDeuteron = std::vector<float>{3.0f, 3.0f, 5.0f};
  auto cvtEdgesProton = std::vector<float>{0, 0, 0, 0, 0};
  auto cvtEdgesDeuteron = std::vector<float>{0, 0, 0, 0, 0};
  if (FiducalSystematics) {
    edgeRegionsElectron = {4.0f, 4.0f, 11.0f};
    edgeRegionsProton = {4.0f, 4.0f, 6.0f};
    edgeRegionsDeuteron = {4.0f, 4.0f, 6.0f};
    cvtEdgesProton = {0.1f, 0.1f, 0.1f, 0.1f, 0.1f};
    cvtEdgesDeuteron = {0.1f, 0.1f, 0.1f, 0.1f, 0.1f};
  }
  trackCuts->SetDCEdgeCuts(11, edgeRegionsElectron);
  trackCuts->SetDCEdgeCuts(2212, edgeRegionsProton);
  trackCuts->SetCVTEdgeCuts(2212, cvtEdgesProton);
  trackCuts->SetDCEdgeCuts(45, edgeRegionsDeuteron);
  trackCuts->SetCVTEdgeCuts(45, cvtEdgesDeuteron);

  const float pcalEdge = FiducalSystematics ? 14.5f : 13.5f;
  for (int sector = 1; sector <= 6; ++sector) {
    trackCuts->AddPCalFiducialRange(
        11, sector, "lw", 0.0f, pcalEdge);
    trackCuts->AddPCalFiducialRange(
        11, sector, "lv", 0.0f, pcalEdge);
  }

  trackCuts->SetMinECALEnergyCut(11, 1, 0.06);
  trackCuts->SetSFCut(true, 11, PIDSystematics ? 0.21 : 0.19, 4.9);
  trackCuts->SetSamplingFractionSystematics(
      PIDSystematics, SFthreeSigmaFactor);

  auto* eventCuts = new EventCut();

  ParticleCut electron;
  electron.pid = 11;
  electron.charge = -1;
  electron.minCount = 1;
  electron.maxCount = 1;
  electron.minFDMomentum = 2.0f;

  // The hadron cuts mark acceptable candidates. Their required multiplicities
  // are imposed by the alternative topologies below, not by the common AND.
  ParticleCut proton;
  proton.pid = 2212;
  proton.charge = 1;
  proton.minCount = 0;
  proton.maxCount = 100;
  proton.minCDMomentum = 0.0f;
  proton.minFDMomentum = 0.0f;

  ParticleCut neutron;
  neutron.pid = 2112;
  neutron.charge = 0;
  neutron.minCount = 0;
  neutron.maxCount = 100;
  neutron.minFTMomentum = 0.0f;
  neutron.minFDMomentum = 0.0f;
  neutron.minCDMomentum = 0.0f;
  neutron.minVz = -8.0f;
  neutron.maxVz = 2.0f;

  ParticleCut deuteron;
  deuteron.pid = 45;
  deuteron.charge = 1;
  deuteron.minCount = 0;
  deuteron.maxCount = 100;
  deuteron.minCDMomentum = 0.0f;
  deuteron.minFDMomentum = 0.0f;
  deuteron.minVz = -8.0f;
  deuteron.maxVz = 2.0f;

  eventCuts->AddParticleCut("electron", electron);
  eventCuts->AddParticleCut("proton", proton);
  eventCuts->AddParticleCut("neutron", neutron);
  eventCuts->AddParticleCut("deuteron", deuteron);

  //eventCuts->AddAlternativeTopology("epp", {{"proton", 2}});
  //eventCuts->AddAlternativeTopology(
  //    "epn", {{"proton", 1}, {"neutron", 1}});
  //eventCuts->AddAlternativeTopology("ed", {{"deuteron", 1}});

  eventCuts->AddAlternativeTopology("epp", {{"proton", 0}});
  eventCuts->AddAlternativeTopology(
      "epn", {{"proton", 0}, {"neutron", 0}});
  eventCuts->AddAlternativeTopology("ed", {{"deuteron", 0}});

  auto qadbCuts = std::make_shared<QADBCuts>();
  const char* qaDefects[] = {
      "TotalOutlier", "TerminalOutlier", "MarginalOutlier", "SectorLoss",
      "LowLiveTime", "Misc", "TotalOutlierFT", "TerminalOutlierFT",
      "MarginalOutlierFT", "LossFT", "BSAWrong", "BSAUnknown", "TSAWrong",
      "TSAUnknown", "DSAWrong", "DSAUnknown", "ChargeHigh",
      "ChargeNegative", "ChargeUnknown", "PossiblyNoBeam"};
  qadbCuts->SetDefects(qaDefects);

  auto disTask =
      std::make_unique<He3DISAnalysis>(
          IsMC, IsreprocRootFile, IsMinimalBook);
  disTask->SetTrackCuts(trackCuts);
  disTask->SetEventCuts(eventCuts);
  disTask->SetBeamEnergy(10.6);
  disTask->SetFTonConfig(false);
  disTask->SetDoFiducialCut(true);
  disTask->SetDoInvMassCut(false);
  disTask->SetDoMomentumCorrection(false);
  disTask->SetMaxEvents(0);
  disTask->SetAcceptEverything(false);
  disTask->SetQADBCuts(qadbCuts);
  if (IsMC) {
    disTask->SetDoQADBCuts(false);
    disTask->SetOptimizeColumns(true);
  } else {
    disTask->SetDoQADBCuts(true);
    disTask->SetOptimizeColumns(false);
  }

  mgr.AddTask(std::move(disTask));
  EventProcessor processor(
      mgr, inputFileDir, outputFileDir, IsreprocRootFile,
      inputRootTreeName, inputRootFileName, nfile, nthreads);
  processor.ProcessEvents();
}

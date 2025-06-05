#include <iostream>

#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/DVCSAnalysis.h"
#include "./../DreamAN/core/EventProcessor.h"

void RunDVCSAnalysis(const std::string& inputDir) {
  std::string inputHipoDir = inputDir;
  const std::string inputHipoDirTest = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/temp_hipofile/";

  if (inputDir.empty()) {
    std::cout << "Input directory for hipo files is empty." << std::endl;
    std::cout << "For testing purpose: I have picked up one file from " << inputHipoDirTest << std::endl;
    inputHipoDir = inputHipoDirTest;
  }

  AnalysisTaskManager mgr;
  mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork");

  // Track Cuts
  auto* trackCuts = new TrackCut();
  trackCuts->SetECALEdgeCut(9, 100000000);
  trackCuts->SetDCEdgeCut(1, 100000000);

  // Photon Cuts
  auto* photonCuts = new EventCut();
  photonCuts->SetChargeCut(0);
  photonCuts->SetPIDCountCut(22, 1, 1);

  // Electron Cuts
  auto* electronCuts = new EventCut();
  electronCuts->SetChargeCut(-1);
  electronCuts->SetPIDCountCut(11, 1, 1);

  // Proton Cuts
  auto* protonCuts = new EventCut();
  protonCuts->SetChargeCut(1);
  protonCuts->SetPIDCountCut(2212, 1, 1);

  // Task
  auto dvcsTask = std::make_unique<DVCSAnalysis>();
  dvcsTask->SetTrackCuts(trackCuts);
  dvcsTask->SetPhotonCuts(photonCuts);
  dvcsTask->SetElectronCuts(electronCuts);
  dvcsTask->SetProtonCuts(protonCuts);

  mgr.AddTask(std::move(dvcsTask));

  // Processor
  EventProcessor processor(inputHipoDir, mgr);
  processor.ProcessEvents();
}

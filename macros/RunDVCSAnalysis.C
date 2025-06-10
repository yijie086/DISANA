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
  mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/afterCalorimeterCuts/");

  // fiducial cuts///
  std::shared_ptr<TrackCut> trackCuts = std::make_shared<TrackCut>();
  auto theta_bins = std::vector<std::pair<float, float>>({{5.0 * M_PI / 180, 10.0 * M_PI / 180},
                                                          {10.0 * M_PI / 180, 15.0 * M_PI / 180},
                                                          {15.0 * M_PI / 180, 20.0 * M_PI / 180},
                                                          {20.0 * M_PI / 180, 25.0 * M_PI / 180},
                                                          {25.0 * M_PI / 180, 30.0 * M_PI / 180}});
  auto edge_regions = std::vector<float>{3.0f, 5.0f, 10.0f};

  trackCuts->SetThetaBins(theta_bins);
  trackCuts->SetEdgeCuts(edge_regions);
  trackCuts->SetSectorCut_Bhawani({1, 2, 3, 4, 5, 6}, 11, 6, true);

  // Sector 1, PCal
  trackCuts->AddPCalFiducialRange(1, "lw", 72.0, 94.5);
  trackCuts->AddPCalFiducialRange(1, "lw", 220.5, 234.0);
  // Sector 2, PCal
  trackCuts->AddPCalFiducialRange(2, "lv", 99.0, 117.5);
  // Sector 3, PCal
  trackCuts->AddPCalFiducialRange(3, "lv", 346.5, 378.0);
  // Sector 4, PCal
  trackCuts->AddPCalFiducialRange(4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(4, "lv", 229.5, 243.0);
  // Sector 6, PCal
  trackCuts->AddPCalFiducialRange(6, "lw", 166.5, 193.5);

  // Sector 1, ECin only
  trackCuts->AddECinFiducialRange(1, "lv", 67.5, 94.5);
  trackCuts->AddECinFiducialRange(4, "lw", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(5, "lv", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(6, "lw", 0.0, 23.5);
  // Sctor 5, ECout only
  trackCuts->AddECoutFiducialRange(1, "lv", 0, 40.5);
  trackCuts->AddECoutFiducialRange(5, "lv", 193.5, 216.0);

  // DVCS particle cuts
  auto* photonCuts = EventCut::PhotonCuts();
  auto* electronCuts = EventCut::ElectronCuts();
  auto* protonCuts = EventCut::ProtonCuts();

  // Task
  auto dvcsTask = std::make_unique<DVCSAnalysis>();
  dvcsTask->SetTrackCuts(trackCuts);
  dvcsTask->SetPhotonCuts(photonCuts);
  dvcsTask->SetElectronCuts(electronCuts);
  dvcsTask->SetProtonCuts(protonCuts);
  dvcsTask->SetBeamEnergy(6.535);
  dvcsTask->SetDoFiducialCut(true);

  mgr.AddTask(std::move(dvcsTask));

  // Processor
  EventProcessor processor(inputHipoDir, mgr);
  processor.ProcessEvents();
}

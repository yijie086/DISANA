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
  mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DC_fiducialcuts/");
  //mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/CheckWithInclusiveData_electron_photon/");

  // fiducial cuts///
  std::shared_ptr<TrackCut> trackCuts = std::make_shared<TrackCut>();
  auto theta_bins = std::vector<std::pair<float, float>>({{5.0 * M_PI / 180, 10.0 * M_PI / 180},
                                                          {10.0 * M_PI / 180, 15.0 * M_PI / 180},
                                                          {15.0 * M_PI / 180, 20.0 * M_PI / 180},
                                                          {20.0 * M_PI / 180, 25.0 * M_PI / 180},
                                                          {25.0 * M_PI / 180, 30.0 * M_PI / 180}});
  // defined here the DC edge cuts for electron and proton                                                        
  auto edge_regions_e = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edge_regions_p = std::vector<float>{3.0f, 3.0f, 5.0f};

  trackCuts->SetDCEdgeCuts(11, edge_regions_e); // DC edge cuts for electrons
  trackCuts->SetDCEdgeCuts(2212, edge_regions_p); // DC edge cuts for protons

  // Cal fiducial cuts for eletron,
  // Sector 1, PCal args PID, sector, side, min, max
  trackCuts->AddPCalFiducialRange(11,1, "lw", 72.0, 94.5);
  trackCuts->AddPCalFiducialRange(11,1, "lw", 220.5, 234.0);
  // Sector 2, PCal
  trackCuts->AddPCalFiducialRange(11,2, "lv", 99.0, 117.5);
  // Sector 3, PCal,
  trackCuts->AddPCalFiducialRange(11,3, "lv", 346.5, 378.0);
  // Sector 4, PCal,
  trackCuts->AddPCalFiducialRange(11,4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11,4, "lv", 229.5, 243.0);
  // Sector 6, PCal,
  trackCuts->AddPCalFiducialRange(11,6, "lw", 166.5, 193.5);

  // Sector 1, ECin only, sector, side, min, max
  trackCuts->AddECinFiducialRange(11,1, "lv", 67.5, 94.5);
  trackCuts->AddECinFiducialRange(11,4, "lw", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(11,5, "lv", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(11,6, "lw", 0.0, 23.5);

  // Sctor 5, ECout only
  trackCuts->AddECoutFiducialRange(11,1, "lv", 0.0, 40.5);
  trackCuts->AddECoutFiducialRange(11,5, "lv", 193.5, 216.0);

  // Cal fiducial cuts for photon, sector, side, min, max
  trackCuts->AddPCalFiducialRange(22,1, "lw", 72.0, 94.5);
  trackCuts->AddPCalFiducialRange(22,1, "lw", 220.5, 234.0);
  // Sector 2, PCal
  trackCuts->AddPCalFiducialRange(22,2, "lv", 99.0, 117.5);
  // Sector 3, PCal,
  trackCuts->AddPCalFiducialRange(22,3, "lv", 346.5, 378.0);
  // Sector 4, PCal,
  trackCuts->AddPCalFiducialRange(22,4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22,4, "lv", 229.5, 243.0);
  // Sector 6, PCal,
  trackCuts->AddPCalFiducialRange(22,6, "lw", 166.5, 193.5);

  // Sector 1, ECin only, sector, side, min, max
  trackCuts->AddECinFiducialRange(22,1, "lv", 67.5, 94.5);
  trackCuts->AddECinFiducialRange(22,4, "lw", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(22,5, "lv", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(22,6, "lw", 0.0, 23.5);

  // Sctor 5, ECout only, sector, side, min, max
  trackCuts->AddECoutFiducialRange(22,1, "lv", 0, 40.5);
  trackCuts->AddECoutFiducialRange(22,5, "lv", 193.5, 216.0);

  // particles for the reaction DVCS: e, p, and gamma
  auto* photonCuts = EventCut::PhotonCuts();
  auto* electronCuts = EventCut::ElectronCuts();
  auto* protonCuts = EventCut::ProtonCuts();

  // Task
  auto dvcsTask = std::make_unique<DVCSAnalysis>();
  dvcsTask->SetTrackCuts(trackCuts);
  dvcsTask->SetPhotonCuts(photonCuts);
  dvcsTask->SetElectronCuts(electronCuts);
  dvcsTask->SetProtonCuts(protonCuts);
  dvcsTask->SetBeamEnergy(10.6);
  dvcsTask->SetDoFiducialCut(true);

  mgr.AddTask(std::move(dvcsTask));
  
  // Processor
  EventProcessor processor(inputHipoDir, mgr);
  processor.ProcessEvents();
}

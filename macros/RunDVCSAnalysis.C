#include <iostream>

#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/DVCSAnalysis.h"
#include "./../DreamAN/core/EventProcessor.h"

void RunDVCSAnalysis(const std::string& inputDir) {
  bool IsMC = false;  // Set to true if you want to run on MC data
  bool IsreprocRootFile = false;  // Set to true if you want to reprocess ROOT files
  std::string inputFileDir = inputDir;
  std::string inputRootFileName = " ";
  std::string inputRootTreeName = " ";
  std::string inputHipoDirTest = " ";

  if (inputDir.empty()) {
    std::cout << "No input custom input dir provided: I have picked up one file from " << inputHipoDirTest << std::endl;
    inputHipoDirTest = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/temp_hipofile/";
    inputFileDir = inputHipoDirTest;
  }
  // If you are reprocossing the existing output you may want to change the path and ttree name if you are using a different ROOT file
  if (IsreprocRootFile) {
    //inputFileDir = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DC_fiducialcuts/";
    inputFileDir = "./";
    inputRootFileName = "dfSelected_before_fiducialCuts.root";
    inputRootTreeName = "dfSelected_before";
  }

  AnalysisTaskManager mgr;
  //mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DC_fiducialcuts/");
  mgr.SetOututDir("./");
  // mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/CheckWithInclusiveData_electron_photon/");

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
  //auto CVT_edge_layers_p = std::vector<float>{-100.0f, -100.0f, -100.0f, -100.0f, -100.0f};  // CVT edge cuts for test
  auto CVT_edge_layers_p = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for protons

  trackCuts->SetDCEdgeCuts(11, edge_regions_e);    // DC edge cuts for electrons
  trackCuts->SetDCEdgeCuts(2212, edge_regions_p);  // DC edge cuts for protons
  trackCuts->SetCVTEdgeCuts(2212, CVT_edge_layers_p);  // CVT edge cuts for protons

  trackCuts->AddCVTFiducialRange(2212, 1, "phi", -110.0, -100.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 1, "phi", 10.0, 20.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 1, "phi", 140.0, 150.0);  // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 3, "phi", -108.0, -98.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 3, "phi", 12.0, 22.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 3, "phi", 140.0, 150.0);  // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 5, "phi", -105.0, -95.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 5, "phi", 15.0, 25.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 5, "phi", 142.0, 152.0);  // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 7, "phi", -102.0, -92.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 7, "phi", 18.0, 28.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 7, "phi", 145.0, 155.0);  // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 12, "phi", -99.0, -89.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 12, "phi", 21.0, 31.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 12, "phi", 148.0, 158.0);  // Proton CVT fiducial cuts

  trackCuts->AddFTCalFiducialRange(22, 1, 0, 0, 0.0, 8.5);  // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, 0, 0, 15.5, 100.0);  // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, -8.42, 9.89, 0.0, 1.6);  // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, -9.89, -5.33, 0.0, 1.6);  // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, -6.15, -13.00, 0.0, 2.3);  // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, 3.7, -6.5, 0.0, 2.0);  // Photon FTCal fiducial cuts

  trackCuts->AddFTCalFiducialRange(11, 1, 0, 0, 0.0, 8.5);  // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, 0, 0, 15.5, 100.0);  // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, -8.42, 9.89, 0.0, 1.6);  // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, -9.89, -5.33, 0.0, 1.6);  // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, -6.15, -13.00, 0.0, 2.3);  // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, 3.7, -6.5, 0.0, 2.0);  // Electron FTCal fiducial cuts


  // Cal fiducial edge cuts for electron and photon,
  trackCuts->AddPCalFiducialRange(11, 1, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 1, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 2, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 2, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 3, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 3, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 4, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 5, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 5, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 6, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 6, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 1, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 1, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 2, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 2, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 3, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 3, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 4, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 5, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 5, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 6, "lw", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 6, "lv", 0.0, 13.5);
  // End of edge fiducial cuts



  // Cal fiducial cuts for eletron,
  // Sector 1, PCal args PID, sector, side, min, max
  trackCuts->AddPCalFiducialRange(11, 1, "lw", 72.0, 94.5);
  trackCuts->AddPCalFiducialRange(11, 1, "lw", 211.5, 234.0);
  // Sector 2, PCal
  trackCuts->AddPCalFiducialRange(11, 2, "lv", 99.0, 117.5);
  // Sector 3, PCal,
  trackCuts->AddPCalFiducialRange(11, 3, "lv", 346.5, 378.0);
  // Sector 4, PCal,
  trackCuts->AddPCalFiducialRange(11, 4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 4, "lv", 229.5, 243.0);
  // Sector 6, PCal,
  trackCuts->AddPCalFiducialRange(11, 6, "lw", 166.5, 193.5);

  // Sector 1, ECin only, sector, side, min, max
  trackCuts->AddECinFiducialRange(11, 1, "lv", 67.5, 94.5);
  trackCuts->AddECinFiducialRange(11, 4, "lw", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(11, 5, "lv", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(11, 6, "lw", 0.0, 23.5);

  // Sctor 5, ECout only
  trackCuts->AddECoutFiducialRange(11, 1, "lv", 0.0, 40.5);
  trackCuts->AddECoutFiducialRange(11, 5, "lv", 193.5, 216.0);

  // Cal fiducial cuts for photon, sector, side, min, max
  trackCuts->AddPCalFiducialRange(22, 1, "lw", 72.0, 94.5);
  trackCuts->AddPCalFiducialRange(22, 1, "lw", 220.5, 234.0);
  // Sector 2, PCal
  trackCuts->AddPCalFiducialRange(22, 2, "lv", 99.0, 117.5);
  // Sector 3, PCal,
  trackCuts->AddPCalFiducialRange(22, 3, "lv", 346.5, 378.0);
  // Sector 4, PCal,
  trackCuts->AddPCalFiducialRange(22, 4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(22, 4, "lv", 229.5, 243.0);
  // Sector 6, PCal,
  trackCuts->AddPCalFiducialRange(22, 6, "lw", 166.5, 193.5);

  // Sector 1, ECin only, sector, side, min, max
  trackCuts->AddECinFiducialRange(22, 1, "lv", 67.5, 94.5);
  trackCuts->AddECinFiducialRange(22, 4, "lw", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(22, 5, "lv", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(22, 6, "lw", 0.0, 23.5);

  // Sctor 5, ECout only, sector, side, min, max
  trackCuts->AddECoutFiducialRange(22, 1, "lv", 0, 40.5);
  trackCuts->AddECoutFiducialRange(22, 5, "lv", 193.5, 216.0);

  // particles for the reaction DVCS: e, p, and gamma
  auto* photonCuts = EventCut::PhotonCuts();
  auto* electronCuts = EventCut::ElectronCuts();
  auto* protonCuts = EventCut::ProtonCuts();


  auto corr = std::make_shared<MomentumCorrection>();
  // AddPiecewiseCorrection pid=11(electron),p∈[0.5,2]GeV, θ∈[10,60]°, φ∈[0,360]°
  corr->AddPiecewiseCorrection(
    11,
    {0.0, 10.0, 0.0*M_PI/180, 180.0*M_PI/180, 0.0*M_PI/180, 360.0*M_PI/180},
    [](double p, double theta, double phi) {
        return 10.0;  // example correction
    }
  );

  corr->AddPiecewiseCorrection(
    2212,
    {0.0, 10.0, 0.0*M_PI/180, 180.0*M_PI/180, 0.0*M_PI/180, 360.0*M_PI/180},
    [](double p, double theta, double phi) {
        return 11.0;  // example correction
    }
  );

  // Task
  auto dvcsTask = std::make_unique<DVCSAnalysis>(IsMC, IsreprocRootFile);
  dvcsTask->SetTrackCuts(trackCuts);
  dvcsTask->SetPhotonCuts(photonCuts);
  dvcsTask->SetElectronCuts(electronCuts);
  dvcsTask->SetProtonCuts(protonCuts);
  dvcsTask->SetBeamEnergy(7.546);
  dvcsTask->SetFTonConfig(true);  // Set to true if you have FT (eq. RGK Fall2018 Pass2 6.535GeV is FT-off)
  dvcsTask->SetDoFiducialCut(true);
  dvcsTask->SetDoMomentumCorrection(true);  // Set to true if you want to apply momentum correction
  dvcsTask->SetMomentumCorrection(corr);  // Set the momentum correction object

  mgr.AddTask(std::move(dvcsTask));

  // Processor
  EventProcessor processor(mgr, inputFileDir, IsreprocRootFile, inputRootTreeName, inputRootFileName);
  processor.ProcessEvents();
}

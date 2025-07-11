#include <iostream>

#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/DVCSAnalysis.h"
#include "./../DreamAN/core/EventProcessor.h"

void RunDVCSAnalysis(const std::string& inputDir, int nfile) {
  bool IsMC = false;              // Set to true if you want to run on MC data
  bool IsreprocRootFile = false;  // Set to true if you want to reprocess ROOT files
  bool IsInbending = false;        // Set to true if you want to run on inbending data
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
    inputFileDir = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Inb/";
    // inputFileDir = "./";
    inputRootFileName = "dfSelected.root";
    inputRootTreeName = "dfSelected";
  }

  AnalysisTaskManager mgr;
  // mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/CheckWithInclusiveData_electron_photon/Inb/");
  // mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/CheckWithInclusiveData_electron_photon/Outb/");
  //mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DVCS_wagon/inb/");
  // mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_sims/test/");
  mgr.SetOututDir("./");
  // mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/test/");
  // mgr.SetOututDir("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/test/");


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
  // auto CVT_edge_layers_p = std::vector<float>{-100.0f, -100.0f, -100.0f, -100.0f, -100.0f};  // CVT edge cuts for test
  auto CVT_edge_layers_p = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for protons

  trackCuts->SetDCEdgeCuts(11, edge_regions_e);        // DC edge cuts for electrons
  trackCuts->SetDCEdgeCuts(2212, edge_regions_p);      // DC edge cuts for protons
  trackCuts->SetCVTEdgeCuts(2212, CVT_edge_layers_p);  // CVT edge cuts for protons

  trackCuts->AddCVTFiducialRange(2212, 1, "phi", -110.0, -100.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 1, "phi", 10.0, 20.0);      // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 1, "phi", 140.0, 150.0);    // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 3, "phi", -108.0, -98.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 3, "phi", 12.0, 22.0);     // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 3, "phi", 140.0, 150.0);   // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 5, "phi", -105.0, -95.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 5, "phi", 15.0, 25.0);     // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 5, "phi", 142.0, 152.0);   // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 7, "phi", -102.0, -92.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 7, "phi", 18.0, 28.0);     // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 7, "phi", 145.0, 155.0);   // Proton CVT fiducial cuts

  trackCuts->AddCVTFiducialRange(2212, 12, "phi", -99.0, -89.0);  // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 12, "phi", 21.0, 31.0);    // Proton CVT fiducial cuts
  trackCuts->AddCVTFiducialRange(2212, 12, "phi", 148.0, 158.0);  // Proton CVT fiducial cuts

  trackCuts->AddFTCalFiducialRange(22, 1, 0, 0, 0.0, 8.5);           // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, 0, 0, 15.5, 100.0);        // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, -8.42, 9.89, 0.0, 1.6);    // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, -9.89, -5.33, 0.0, 1.6);   // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, -6.15, -13.00, 0.0, 2.3);  // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, 3.7, -6.5, 0.0, 2.0);      // Photon FTCal fiducial cuts

  trackCuts->AddFTCalFiducialRange(11, 1, 0, 0, 0.0, 8.5);           // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, 0, 0, 15.5, 100.0);        // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, -8.42, 9.89, 0.0, 1.6);    // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, -9.89, -5.33, 0.0, 1.6);   // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, -6.15, -13.00, 0.0, 2.3);  // Electron FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(11, 1, 3.7, -6.5, 0.0, 2.0);      // Electron FTCal fiducial cuts

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
  /// check the edge cut for photon is different its loose in RGA analysis note!
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
  // trackCuts->AddPCalFiducialRange(11, 2, "lv", 99.0, 117.5); /// RGA spring2018 does not have dead zone
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

  /// set sampling fraction for the particle in detector
  trackCuts->SetMinECALEnergyCut(11, 1, 0.06);  // Electron PCal layer 1
  // apply sampling fraction and diagolal cuts tbd!!
  trackCuts->SetSFCut(true,11,0.19,4.9);  // Set to true if you want to apply sampling fraction cut
  trackCuts->AddSamplingFractionMinCut(11, 1, 0.167070, 0.010930, -0.001085);  // Electronsector 1
  trackCuts->AddSamplingFractionMaxCut(11, 1, 0.302129, -0.000735, -0.000409);  // Electro sector 1
  trackCuts->AddSamplingFractionMinCut(11, 2, 0.146344, 0.023547, -0.002660);  // Electronsector 2
  trackCuts->AddSamplingFractionMaxCut(11, 2, 0.313502, -0.004386, 0.000036);  // Electronsector 2
  trackCuts->AddSamplingFractionMinCut(11, 3, 0.136665, 0.031226, -0.003705);  // Electronsector 3
  trackCuts->AddSamplingFractionMaxCut(11, 3, 0.310201, -0.001985, -0.000242);  // Electro sector 3
  trackCuts->AddSamplingFractionMinCut(11, 4, 0.152493, 0.018234, -0.001604);  // Electronsector 4
  trackCuts->AddSamplingFractionMaxCut(11, 4, 0.315967, -0.004037, 0.000058);  // Electronsector 4
  trackCuts->AddSamplingFractionMinCut(11, 5, 0.143086, 0.023846, -0.002510);  // Electronsector 5
  trackCuts->AddSamplingFractionMaxCut(11, 5, 0.321858, -0.009117, 0.000657);  // Electronsector 5
  trackCuts->AddSamplingFractionMinCut(11, 6, 0.152500, 0.021085, -0.002384);  // Electronsector 6
  trackCuts->AddSamplingFractionMaxCut(11, 6, 0.313124, -0.003053, -0.000193);  // Electro sector 6

  // particles for the reaction DVCS: e, p, and gamma
  EventCut* eventCuts = new EventCut();


  ParticleCut proton;
  proton.pid = 2212;                              // Proton PID
  proton.charge = 1;                              // Proton charge
  proton.minCount = 1;                            // Minimum count of protons
  proton.maxCount = 1;                            // Maximum count of protons
  proton.minCDMomentum = 0.3f;                    // Minimum momentum for protons
  if (IsInbending) proton.minFDMomentum = 0.42f;  // Minimum momentum for protons in FD Inbending
  if (!IsInbending) proton.minFDMomentum = 0.5f;  // Minimum momentum for protons in FD Outbending
  proton.maxCDMomentum = 1.2f;                    // Maximum momentum for protons
  proton.maxFDMomentum = 1.2f;                    // Maximum momentum for protons
  proton.minTheta = 0.0f;                         // Minimum theta for protons
  proton.maxTheta = 64.23 * M_PI / 180.0;         // Maximum theta for protons (approximately 64.23 degrees)
  proton.minPhi = 0.0f;                           // Minimum phi for protons
  proton.maxPhi = 2.0f * M_PI;

/*
  ParticleCut proton;
  proton.pid = 2212;                              // Proton PID
  proton.charge = 1;                              // Proton charge
  proton.minCount = 1;                            // Minimum count of protons
  proton.maxCount = 1;                            // Maximum count of protons
  proton.minCDMomentum = 0.3f;                    // Minimum momentum for protons
  if (IsInbending) proton.minFDMomentum = 0.3f;  // Minimum momentum for protons in FD Inbending
  if (!IsInbending) proton.minFDMomentum = 0.4f;  // Minimum momentum for protons in FD Outbending
*/
  ParticleCut electron;
  electron.pid = 11;              // Electron PID
  electron.charge = -1;           // Electron charge
  electron.minCount = 1;          // Minimum count of electrons
  electron.maxCount = 1;          // Maximum count of electrons
  electron.minFDMomentum = 2.0f;  // Minimum momentum for electrons

  ParticleCut photon;
  photon.pid = 22;      // Photon PID
  photon.charge = 0;    // Photon charge
  photon.minCount = 1;  // Minimum count of photons
  photon.minFDMomentum = 2.0f;
  photon.minFTMomentum = 2.0f;
  photon.minBeta = 0.9f;  // Minimum beta for photons
  photon.maxBeta = 1.1f;  // Maximum beta for photons

  // pi0 background studies, similar is valid for exclusive meson production
  TwoBodyMotherCut pi0;
  pi0.expectedMotherMass = 0.132f;  // determined from the fit inv mass of two photon
  pi0.massSigma = 0.0129f;          // determined from the fit inv mass of two photon
  pi0.nSigmaMass = 3.0;             // choice for nsigma is user dependent

  eventCuts->AddParticleCut("proton", proton);      // Applies defaults automatically
  eventCuts->AddParticleCut("electron", electron);  // Applies defaults automatically
  eventCuts->AddParticleCut("photon", photon);      // Applies defaults automatically
  eventCuts->AddParticleMotherCut("pi0", pi0);      // Applies defaults automatically

  auto corr = std::make_shared<MomentumCorrection>();
  corr->AddPiecewiseCorrection(  // Momentum correction for proton RGA Fa18 in
      2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD}, [](double p, double theta, double phi) {
        theta = theta * 180.0 / M_PI;  // Convert theta to degrees
        float A_p = -0.0285153 + 0.000678814 * theta;
        float B_p = 0.0468792 - 0.00089179 * theta;
        float C_p = -0.026879 + 0.000431868 * theta;
        return p + (A_p + B_p * p + C_p * p * p);
      });

  corr->AddPiecewiseCorrection(  // Momentum correction for proton RGA Fa18 in
      2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD}, [](double p, double theta, double phi) {
        theta = theta * 180.0 / M_PI;  // Convert theta to degrees
        float A_p = 0.0120379 - 0.00173099 * theta + 0.000158133 * theta * theta - 6.23728e-06 * theta * theta * theta + 9.24082e-08 * theta * theta * theta * theta;
        float B_p = -0.0128831 + 0.00384179 * theta - 0.00044874 * theta * theta + 1.96184e-05 * theta * theta * theta - 2.81305e-07 * theta * theta * theta * theta;
        float C_p = -0.0102071 + 0.00370448 * theta - 0.000250471 * theta * theta + 6.16952e-06 * theta * theta * theta - 3.72598e-08 * theta * theta * theta * theta;

        return p + (A_p + B_p / p + C_p / (p * p));
      });

  // Task
  auto dvcsTask = std::make_unique<DVCSAnalysis>(IsMC, IsreprocRootFile);
  dvcsTask->SetTrackCuts(trackCuts);
  dvcsTask->SetEventCuts(eventCuts);
  dvcsTask->SetBeamEnergy(10.6);
  dvcsTask->SetFTonConfig(true);  // Set to true if you have FT (eq. RGK Fall2018 Pass2 6.535GeV is FT-off)
  dvcsTask->SetDoFiducialCut(true);

  dvcsTask->SetDoInvMassCut(true); // in this case pi0 background for two-photon pairs in the event
  dvcsTask->SetDoMomentumCorrection(true);  // Set to true if you want to apply momentum correction
  dvcsTask->SetMomentumCorrection(corr);  // Set the momentum correction object
  dvcsTask->SetMaxEvents(0);  // Set the maximum number of events to process, 0 means no limit


  mgr.AddTask(std::move(dvcsTask));

  // Processor
  EventProcessor processor(mgr, inputFileDir, IsreprocRootFile, inputRootTreeName, inputRootFileName, nfile);
  processor.ProcessEvents();
}

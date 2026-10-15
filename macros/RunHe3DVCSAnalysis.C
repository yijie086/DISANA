#include <iostream>

#include "./../DreamAN/Cuts/QADBCuts.h"
#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/He3DVCSAnalysis.h"
#include "./../DreamAN/core/EventProcessor.h"

void RunHe3DVCSAnalysis(const std::string& inputDir, int nfile, int nthreads = 0) {
  if (nthreads <= 0) {
    const auto hc = std::thread::hardware_concurrency();
    nthreads = hc ? static_cast<int>(hc) : 2;  // sensible default
  }

  // Enable implicit multi-threading for ROOT.
  if (nthreads > 1) {
    ROOT::EnableImplicitMT();
    std::cout << "[RunHe3DVCSAnalysis] IMT enabled with " << nthreads << " thread(s)\n";
  } else {
    std::cout << "[RunHe3DVCSAnalysis] IMT disabled (single thread mode).\n";
  }
  bool IsMC = true;              // Set to true if you want to run on MC data
  bool IsreprocRootFile = false;  // Set to true if you want to reprocess ROOT files
  bool IsMinimalBook = true;
  bool FiducalSystematics = false;  // Set to true if you want to run systematic variations for fiducial cuts
  bool PIDSystematics = false;  // Set to true if you want to run systematic variations for PID cuts
  float SFthreeSigmaFactor = 5.0/6.0;
  std::cout << "Running He3 DVCS Analysis" << std::endl;

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
    inputRootFileName = "dfSelected.root";
    inputRootTreeName = "dfSelected";
  }

  AnalysisTaskManager mgr;
  std::string outputFileDir = "./";
  mgr.SetOututDir(outputFileDir);
  // fiducial cuts///
  std::shared_ptr<TrackCut> trackCuts = std::make_shared<TrackCut>();
  /*auto theta_bins = std::vector<std::pair<float, float>>({{5.0 * M_PI / 180, 10.0 * M_PI / 180},
                                                          {10.0 * M_PI / 180, 15.0 * M_PI / 180},
                                                          {15.0 * M_PI / 180, 20.0 * M_PI / 180},
                                                          {20.0 * M_PI / 180, 25.0 * M_PI / 180},
                                                          {25.0 * M_PI / 180, 30.0 * M_PI / 180}});*/
  // DC/CVT edge cuts for the electron and recoil He3.
  auto edge_regions_e = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edge_regions_he3 = std::vector<float>{3.0f, 3.0f, 5.0f};
  auto CVT_edge_layers_he3 = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
  if (FiducalSystematics){
    edge_regions_e = std::vector<float>{4.0f, 4.0f, 11.0f};  // DC edge cuts for electrons for systematic variation
    edge_regions_he3 = std::vector<float>{4.0f, 4.0f, 6.0f};
    CVT_edge_layers_he3 = std::vector<float>{0.1f, 0.1f, 0.1f, 0.1f, 0.1f};
  }

  trackCuts->SetDCEdgeCuts(11, edge_regions_e);        // DC edge cuts for electrons
  trackCuts->SetDCEdgeCuts(49, edge_regions_he3);
  trackCuts->SetCVTEdgeCuts(49, CVT_edge_layers_he3);

  trackCuts->AddFTCalFiducialRange(22, 1, 0, 0, 0.0, 8.5);           // Photon FTCal fiducial cuts
  trackCuts->AddFTCalFiducialRange(22, 1, 0, 0, 15.5, 100.0);        // Photon FTCal fiducial cuts

  if (!FiducalSystematics){
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
    /// Photon edge cut is intentionally loose for this He3 selection.
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
    std::cout << "Applied standard PCal fiducial cuts (0.0 to 13.5) for electrons and photons in all sectors.\n";
    // End of edge fiducial cuts
  } else {
    trackCuts->AddPCalFiducialRange(11, 1, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 1, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 2, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 2, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 3, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 3, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 4, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 4, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 5, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 5, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 6, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(11, 6, "lv", 0.0, 14.5);
    
    trackCuts->AddPCalFiducialRange(22, 1, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 1, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 2, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 2, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 3, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 3, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 4, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 4, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 5, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 5, "lv", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 6, "lw", 0.0, 14.5);
    trackCuts->AddPCalFiducialRange(22, 6, "lv", 0.0, 14.5);
    std::cout << "Applied systematic variation PCal fiducial cuts (0.0 to 14.5) for electrons and photons in all sectors.\n";
  }


  if (!PIDSystematics){
    /// set sampling fraction for the particle in detector
    trackCuts->SetMinECALEnergyCut(11, 1, 0.06);  // Electron PCal layer 1
    // apply sampling fraction and diagolal cuts tbd!!
    trackCuts->SetSFCut(true, 11, 0.19, 4.9);  // Set to true if you want to apply sampling fraction cut
    std::cout << "Applied standard sampling fraction cut for electrons in PCal layer 1 with min energy cut of 0.06 and SF cut of 0.19.\n";
  } else {
    trackCuts->SetMinECALEnergyCut(11, 1, 0.06);  // Electron PCal layer 1
    trackCuts->SetSFCut(true, 11, 0.21, 4.9);  // Set to true if you want to apply sampling fraction cut
    std::cout << "Applied systematic variation sampling fraction cut for electrons in PCal layer 1 with min energy cut of 0.06 and SF cut of 0.21.\n";
  }
  trackCuts->SetSamplingFractionSystematics(PIDSystematics, SFthreeSigmaFactor);
  if (PIDSystematics) {
    std::cout << "Sampling-fraction min/max curves rescaled around their midpoint "
              << "with SFthreeSigmaFactor = " << SFthreeSigmaFactor << ".\n";
  }
  // Final state: e + He3 + gamma.
  EventCut* eventCuts = new EventCut();

  ParticleCut he3;
  he3.pid = 49;
  he3.charge = 1;  // REC::Particle stores the reconstructed positive-track charge sign.
  he3.minCount = 0;
  he3.maxCount = 100;
  he3.minCDMomentum = 0.0f;
  he3.minFDMomentum = 0.0f;

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
  photon.minFDMomentum = 0.15f;
  photon.minFTMomentum = 0.15f;
  photon.minBeta = 0.9f;  // Minimum beta for photons
  photon.maxBeta = 1.1f;  // Maximum beta for photons

  // pi0 background studies, similar is valid for exclusive meson production
  //TwoBodyMotherCut pi0;
  //pi0.expectedMotherMass = 0.132f;  // determined from the fit inv mass of two photon
  //pi0.massSigma = 0.0129f;          // determined from the fit inv mass of two photon
  //pi0.nSigmaMass = 3.0;             // choice for nsigma is user dependent

  eventCuts->AddParticleCut("he3", he3);
  eventCuts->AddParticleCut("electron", electron);  // Applies defaults automatically
  eventCuts->AddParticleCut("photon", photon);      // Applies defaults automatically
  //eventCuts->AddParticleMotherCut("pi0", pi0);      // Applies defaults automatically

  /// QADB cut parameters for the He3 data selection

  auto qadbCuts = std::make_shared<QADBCuts>();
  const char* qaDefects[] = {"TotalOutlier",      "TerminalOutlier",   "MarginalOutlier", "SectorLoss",     "LowLiveTime",   "Misc",          "TotalOutlierFT",
                             "TerminalOutlierFT", "MarginalOutlierFT", "LossFT",          "BSAWrong",       "BSAUnknown",    "TSAWrong",      "TSAUnknown",
                             "DSAWrong",          "DSAUnknown",        "ChargeHigh",      "ChargeNegative", "ChargeUnknown", "PossiblyNoBeam"};

  qadbCuts->SetDefects(qaDefects);

                                   // qadbCuts->AddDefect("SomeExtraDefect");
                                   // qadbCuts->SetAccumulateCharge(true); // default true
  // Task
  auto dvcsTask = std::make_unique<He3DVCSAnalysis>(IsMC, IsreprocRootFile, IsMinimalBook);
  dvcsTask->SetTrackCuts(trackCuts);
  dvcsTask->SetEventCuts(eventCuts);
  dvcsTask->SetBeamEnergy(10.6);
  dvcsTask->SetFTonConfig(true);
  dvcsTask->SetDoFiducialCut(true);

  dvcsTask->SetDoInvMassCut(true);          // in this case pi0 background for two-photon pairs in the event
  dvcsTask->SetMaxEvents(0);                // Set the maximum number of events to process, 0 means no limit
  dvcsTask->SetAcceptEverything(false);     // Set to true to accept all events, false to apply cuts
  dvcsTask->SetQADBCuts(qadbCuts);          // <-- this now matters
  if(IsMC) {
    dvcsTask->SetDoQADBCuts(false);  // for MC we usually do not apply QADB false rejection
    dvcsTask->SetOptimizeColumns(true);       // Set to true to optimize the columns in the output file based on the cuts applied
  } else {
    dvcsTask->SetDoQADBCuts(true);   // for data we usually apply QADB false rejection
    dvcsTask->SetOptimizeColumns(false);       // Set to true to optimize the columns in the output file based on the cuts applied
  }

  mgr.AddTask(std::move(dvcsTask));

  // Processor
  EventProcessor processor(mgr, inputFileDir, outputFileDir, IsreprocRootFile, inputRootTreeName, inputRootFileName, nfile, nthreads);
  processor.ProcessEvents();
}

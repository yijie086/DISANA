#include <TROOT.h>
#include <TSystem.h>

#include <ROOT/RDataFrame.hxx>  // NEW
#include <iostream>
#include <thread>  // NEW

#include "./../DreamAN/Cuts/QADBCuts.h"
#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/EventProcessor.h"
#include "./../DreamAN/core/PhiAnalysis.h"

void RunPhiAnalysis(const std::string& inputDir, int nfile, int nthreads,
                    const std::string& outputDir,
                    const std::string dataconfig,
                    bool IsMC, bool IsreprocRootFile,
                    bool IsInbending /* you can ignore and derive from dataconfig if you like */,
                    bool IsMinimalBook,
                    bool IsMissingKm,
                    const std::string& reprocRootFile = "",
                    const std::string& reprocTreeName = "") {
  // Determine the number of threads to use.
  if (nthreads <= 0) {
    const auto hc = std::thread::hardware_concurrency();
    nthreads = hc ? static_cast<int>(hc) : 2;  // sensible default
  }

  // Enable implicit multi-threading for ROOT.
  if (nthreads > 1) {
    ROOT::EnableImplicitMT();
    std::cout << "[RunPhiAnalysis] IMT enabled with " << nthreads << " thread(s)\n";
  } else {
    std::cout << "[RunPhiAnalysis] IMT disabled (single thread mode).\n";
  }

  if (dataconfig == "rgkfa18_7546") {
    IsInbending = false;  // Set to false for outbending data
  }
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgafall18_inb" || dataconfig == "rgkfa18_7546" || dataconfig == "rgasp19_inb") {
    IsInbending = true;  // Set to true for inbending data
  }

  std::cout << "Running Phi Analysis with configuration: " << dataconfig << std::endl;

    // ======================================================================
  // NEW QADB Configuration: Allow 'Misc' bit for specific runs (e.g., empty target)
  // This is the C++ equivalent of the Groovy script's global 'allowMiscBit' call.
  // It configures the static QA::QADB instance wrapped by QADBCuts.
  std::vector<int> miscAllowedRuns = {
    5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
    5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
    5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
    5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
    6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
    6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
    6757,                         // RGA runs
    16194, 16089, 16185, 16308, 16184, 16307, 16309, // RGC Su22 He/ET
    16872, 16975,                                  // RGC Fa22 He/ET
    17763, 17764, 17765, 17766, 17767, 17768,      // RGC Sp23 He/ET
    17179, 17180, 17181, 17182, 17183, 17188, 17189, // RICH off/partially down
    17252
  };
  
  std::string inputFileDir = inputDir;
  std::string outputFileDir = outputDir;
  std::string inputRootFileName = " ";
  std::string inputRootTreeName = " ";
  std::string inputHipoDirTest = " ";

  if (inputDir.empty()) {
    std::cout << "No input custom input dir provided: I have picked up one file from " << inputHipoDirTest << std::endl;
    inputHipoDirTest = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/temp_hipofile/";
    inputFileDir = inputHipoDirTest;
  }

  if (IsreprocRootFile) {
    // Use custom paths if provided, otherwise use defaults
    if (!reprocRootFile.empty()) {
      inputRootFileName = reprocRootFile;
      std::cout << "[RunPhiAnalysis] Using custom reproc file: " << inputRootFileName << std::endl;
    } else {
      // Default filename if not specified
      inputRootFileName = "DVKpKm.root";
      std::cout << "[RunPhiAnalysis] Using default reproc file: " << inputRootFileName << std::endl;
    }

    if (!reprocTreeName.empty()) {
      inputRootTreeName = reprocTreeName;
      std::cout << "[RunPhiAnalysis] Using custom tree name: " << inputRootTreeName << std::endl;
    } else {
      // Default tree name if not specified
      inputRootTreeName = "hipo";
      std::cout << "[RunPhiAnalysis] Using default tree name: " << inputRootTreeName << std::endl;
    }

    // If inputDir was not changed from default, use a default reprocessed data path
    if (inputFileDir == ".") {
      inputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/HIPO_test/newhipo2root/hipo-utils/test_outputs/";
      std::cout << "[RunPhiAnalysis] No input dir specified for reproc, using default: " << inputFileDir << std::endl;
    }
  }

  AnalysisTaskManager mgr;
  mgr.SetOututDir(outputFileDir);

  // fiducial cuts///
  std::shared_ptr<TrackCut> trackCuts = std::make_shared<TrackCut>();
  // defined here the DC edge cuts for electron and proton
  auto edge_regions_e = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edge_regions_p = std::vector<float>{3.0f, 3.0f, 7.0f};
  auto edge_regions_kM = std::vector<float>{3.0f, 3.0f, 7.0f};
  auto edge_regions_kP = std::vector<float>{3.0f, 3.0f, 7.0f};

  // auto CVT_edge_layers_p = std::vector<float>{-100.0f, -100.0f, -100.0f, -100.0f, -100.0f};  // CVT edge cuts for test
  auto CVT_edge_layers_p = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};   // CVT edge cuts for protons
  auto CVT_edge_layers_kM = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for Kaon M
  auto CVT_edge_layers_kP = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for KaonP

  trackCuts->SetDCEdgeCuts(11, edge_regions_e);     // DC edge cuts for electrons
  trackCuts->SetDCEdgeCuts(2212, edge_regions_p);   // DC edge cuts for protons
  trackCuts->SetDCEdgeCuts(-321, edge_regions_kM);  // DC edge cuts for kM
  trackCuts->SetDCEdgeCuts(321, edge_regions_kP);   // DC edge cuts for kP

  trackCuts->SetCVTEdgeCuts(-321, CVT_edge_layers_p);   // CVT edge cuts for proton
  trackCuts->SetCVTEdgeCuts(-321, CVT_edge_layers_kM);  // CVT edge cuts for Kaon M
  trackCuts->SetCVTEdgeCuts(321, CVT_edge_layers_kP);   // CVT edge cuts for KaonP

  // CVT fiducial cuts for proton from RGA common analysis note
  const std::array<int, 3> pids = {2212, 321, -321};
  // For each layer, keep the exact phi ranges you had (min,max) triplets
  const std::map<int, std::array<std::pair<double, double>, 3>> phiRangesByLayer = {
      {1, {std::pair{-110.0, -100.0}, std::pair{10.0, 20.0}, std::pair{140.0, 150.0}}}, {3, {std::pair{-108.0, -98.0}, std::pair{12.0, 22.0}, std::pair{140.0, 150.0}}},
      {5, {std::pair{-105.0, -95.0}, std::pair{15.0, 25.0}, std::pair{142.0, 152.0}}},  {7, {std::pair{-102.0, -92.0}, std::pair{18.0, 28.0}, std::pair{145.0, 155.0}}},
      {12, {std::pair{-99.0, -89.0}, std::pair{21.0, 31.0}, std::pair{148.0, 158.0}}},
  };

  // Now loop over species, layers, and ranges
  for (const int pid : pids) {
    for (const auto& [layer, ranges] : phiRangesByLayer) {
      for (const auto& [phiMin, phiMax] : ranges) {
        trackCuts->AddCVTFiducialRange(pid, layer, "phi", phiMin, phiMax);
      }
    }
  }

  // ===== FT cal fiducial cuts for electron ===== checked with RGA common analysis note
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgafall18_inb" || dataconfig == "rgasp19_inb" || dataconfig == "rgasp18_outb" || dataconfig == "rgafall18_outb") {
    trackCuts->AddFTCalFiducialRange(11, 1, 0, 0, 0.0, 8.5);           // Electron FTCal fiducial cuts
    trackCuts->AddFTCalFiducialRange(11, 1, 0, 0, 15.5, 100.0);        // Electron FTCal fiducial cuts
    trackCuts->AddFTCalFiducialRange(11, 1, -8.42, 9.89, 0.0, 1.6);    // Electron FTCal fiducial cuts
    trackCuts->AddFTCalFiducialRange(11, 1, -9.89, -5.33, 0.0, 1.6);   // Electron FTCal fiducial cuts
    trackCuts->AddFTCalFiducialRange(11, 1, -6.15, -13.00, 0.0, 2.3);  // Electron FTCal fiducial cuts
    trackCuts->AddFTCalFiducialRange(11, 1, 3.7, -6.5, 0.0, 2.0);      // Electron FTCal fiducial cuts
  }

  // Cal fiducial edge cuts for electron medium cuts for the lower edge
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

  // momentum correction
  std::shared_ptr<MomentumCorrection> corr = std::make_shared<MomentumCorrection>();

  // Spring 2018
  if (dataconfig == "rgasp18_inb") {
    // CD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = -0.2399996 + 0.0124935 * theta - 0.0001621 * theta * theta;
                                   double B_p = 0.62191181 - 0.03248928 * theta + 0.00042729 * theta * theta;
                                   double C_p = -0.45094157 + 0.02272091 * theta - 0.00028791 * theta * theta;
                                   return p + (A_p + B_p * p + C_p * p * p);
                                 });

    // FD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = 0.0089759 - 0.0001757 * theta - 0.0000031 * theta * theta;
                                   double B_p = -0.01186058 + 0.00023068 * theta + 0.00001345 * theta * theta;
                                   double C_p = 0.01162939 - 0.00047203 * theta + 0.00000630 * theta * theta;
                                   return p + (A_p + B_p / p + C_p / (p * p));
                                 });
  }

  // Spring 2018 Outbending
  if (dataconfig == "rgasp18_outb") {
    // CD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = -0.1998669 + 0.0102981 * theta - 0.0001355 * theta * theta;
                                   double B_p = 0.47190846 - 0.02474115 * theta + 0.00033280 * theta * theta;
                                   double C_p = -0.35149154 + 0.01764687 * theta - 0.00022812 * theta * theta;
                                   return p + (A_p + B_p * p + C_p * p * p);
                                 });

    // FD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = 0.0130398 - 0.0005187 * theta;
                                   double B_p = -0.01972407 + 0.00113010 * theta;
                                   double C_p = 0.0;
                                   return p + (A_p + B_p / p + C_p / (p * p));
                                 });
  }

  // ---------- Fall 2018 (inb) ----------
  if (dataconfig == "rgafall18_inb") {
    // CD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = -0.1904655 + 0.0099100 * theta - 0.0001281 * theta * theta;
                                   double B_p = 0.46119797 - 0.02395829 * theta + 0.00031211 * theta * theta;
                                   double C_p = -0.34158626 + 0.01713694 * theta - 0.00021835 * theta * theta;
                                   return p + (A_p + B_p * p + C_p * p * p);
                                 });

    // FD branch (forward detector)
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = 0.0099626 - 0.0002414 * theta - 0.0000020 * theta * theta;
                                   double B_p = -0.01428267 + 0.00042833 * theta + 0.00001081 * theta * theta;
                                   double C_p = 0.01197102 - 0.00055673 * theta + 0.00000785 * theta * theta;
                                   return p + (A_p + B_p / p + C_p / (p * p));
                                 });
  }

  // ---------- Fall 2018 (outb) ----------
  if (dataconfig == "rgafall18_outb") {
    // CD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = -0.1927861 + 0.0099546 * theta - 0.0001299 * theta * theta;
                                   double B_p = 0.44307822 - 0.02309469 * theta + 0.00030784 * theta * theta;
                                   double C_p = -0.32938000 + 0.01648659 * theta - 0.00021181 * theta * theta;
                                   return p + (A_p + B_p * p + C_p * p * p);
                                 });

    // FD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = 0.0135790 - 0.0005303 * theta;
                                   double B_p = -0.02165929 + 0.00121123 * theta;
                                   double C_p = 0.0;
                                   return p + (A_p + B_p / p + C_p / (p * p));
                                 });
  }
  // ---------- Spring 2019 (inb) ----------
  if (dataconfig == "rgasp19_inb") {
    // CD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = -0.2716918 + 0.0142491 * theta - 0.0001862 * theta * theta;
                                   double B_p = 0.65945101 - 0.03431360 * theta + 0.00045036 * theta * theta;
                                   double C_p = -0.46602726 + 0.02335623 * theta - 0.00029720 * theta * theta;
                                   return p + (A_p + B_p * p + C_p * p * p);
                                 });

    // FD branch
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = 0.0095205 - 0.0001914 * theta - 0.0000031 * theta * theta;
                                   double B_p = -0.01365658 + 0.00036322 * theta + 0.00001217 * theta * theta;
                                   double C_p = 0.01175256 - 0.00053407 * theta + 0.00000742 * theta * theta;
                                   return p + (A_p + B_p / p + C_p / (p * p));
                                 });
  }

  // kaon energy loss corrections to be added here: TBD

  // QADB cuts
  auto qadbCuts = std::make_shared<QADBCuts>();
  const char* qaDefects[] = {"TotalOutlier", "TerminalOutlier", "MarginalOutlier", "SectorLoss",   "Misc", "ChargeHigh", "ChargeNegative", "ChargeUnknown", "PossiblyNoBeam"};
  qadbCuts->SetDefects(qaDefects);
  qadbCuts->AddExcludedRun(3262);  // outbendng RGA spring 2018
  qadbCuts->SetAllowedMiscRuns(miscAllowedRuns); 
  std::cout << "[RunPhiAnalysis] Configured QADB to allow Misc bit for " 
            << miscAllowedRuns.size() << " runs (Dilution Factor Runs).\n";
  // particles for the reaction DVCS: e, p, and gamma
  EventCut* eventCuts = new EventCut();
  ParticleCut proton;
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgasp19_inb" || dataconfig == "rgafall18_inb" || dataconfig == "rgasp18_outb" || dataconfig == "rgafall18_outb" ||
      dataconfig == "rgasp19_outb") {
    proton.pid = 2212;    // Proton PID
    proton.charge = 1;    // Proton charge
    proton.minCount = 1;  // Minimum count of protons
  }

  ParticleCut electron;
  electron.pid = 11;              // Electron PID
  electron.charge = -1;           // Electron charge
  electron.minCount = 1;          // Minimum count of electrons
  electron.minFDMomentum = 1.0f;  // Minimum momentum for electrons

  ParticleCut kMinus;
  kMinus.pid = -321;    // KaonM PID
  kMinus.charge = -1;   // charge
  kMinus.minCount = 1;  // Minimum count of Negative kaons

  ParticleCut kPos;
  kPos.pid = 321;  // KaonPos PID
  kPos.charge = 1;
  kPos.minCount = 1;  // Minimum count of possible kaons

  // pi0 background studies, similar is valid for exclusive meson production
  TwoBodyMotherCut phi;
  phi.expectedMotherMass = 1.019f;  // determined from the fit inv mass of two photon
  phi.massSigma = 0.5f;             // determined from the fit inv mass of two photon
  phi.nSigmaMass = 100;             // choice for nsigma is user dependent

  eventCuts->AddParticleCut("proton", proton);      // Applies defaults automatically
  eventCuts->AddParticleCut("electron", electron);  // Applies defaults automatically
  
  if (IsMissingKm && IsInbending) {
    std::cout << "Adding missing K+ particle cut.\n";
    eventCuts->AddParticleCut("Pos Kaon", kPos);
  } else if (IsMissingKm && !IsInbending) {
    std::cout << "Adding missing K- particle cut.\n";
    eventCuts->AddParticleCut("Neg Kaon", kMinus);
  } else {
    std::cout << "Adding missing K- and K+ particle cut.\n";
    eventCuts->AddParticleCut("Neg Kaon", kMinus);  // Applies defaults automatically
    eventCuts->AddParticleCut("Pos Kaon", kPos);    // Applies defaults automatically
  }
  // Task
  auto PhiTask = std::make_unique<PhiAnalysis>(IsMC, IsreprocRootFile, IsMinimalBook);
  PhiTask->SetTrackCuts(trackCuts);
  PhiTask->SetEventCuts(eventCuts);
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgasp18_outb" || dataconfig == "rgafall18_inb" || dataconfig == "rgafall18_outb") {
    PhiTask->SetBeamEnergy(10.6);
  }else if (dataconfig == "rgasp19_outb") {
     PhiTask->SetBeamEnergy(10.2);
  }

  PhiTask->SetFTonConfig(true);  // Set to true if you have FT (eq. RGK Fall2018 Pass2 6.535GeV is FT-off)
  PhiTask->SetDoFiducialCut(true);

  PhiTask->SetDoInvMassCut(true);          // in this case pi0 background for two-photon pairs in the event
  PhiTask->SetDoMomentumCorrection(true);  // Set to true if you want to apply momentum correction
  PhiTask->SetMomentumCorrection(corr);    // Set the momentum correction object
  PhiTask->SetMaxEvents(0);                // Set the maximum number of events to process, 0 means no
  PhiTask->SetQADBCuts(qadbCuts);          // <-- this now matters

  if (IsMC) {
    PhiTask->SetDoQADBCuts(false);  // for MC we usually do not apply QADB false rejection
  } else {
    PhiTask->SetDoQADBCuts(true);  // for RGA sp18 inb data we skip QADB cuts due to missing QA info
  }

  mgr.AddTask(std::move(PhiTask));
  // Processor
    // --------------------------------------------------
  // Print summary of input arguments
  // --------------------------------------------------
  std::cout << "\n================ RunPhiAnalysis INPUT ARGUMENTS ================\n";
  std::cout << "inputDir            : " << (inputDir.empty() ? "<empty>" : inputDir) << "\n";
  std::cout << "nfile               : " << nfile << "\n";
  std::cout << "nthreads (requested): " << nthreads << "\n";
  std::cout << "outputDir           : " << outputDir << "\n";
  std::cout << "dataconfig          : " << dataconfig << "\n";
  std::cout << "IsMC                : " << (IsMC ? "true" : "false") << "\n";
  std::cout << "IsreprocRootFile    : " << (IsreprocRootFile ? "true" : "false") << "\n";
  std::cout << "IsInbending (arg)   : " << (IsInbending ? "true" : "false") << "\n";
  std::cout << "IsMinimalBook       : " << (IsMinimalBook ? "true" : "false") << "\n";
  std::cout << "IsMissingKm         : " << (IsMissingKm ? "true" : "false") << "\n";
  std::cout << "reprocRootFile      : " 
            << (reprocRootFile.empty() ? "<default>" : reprocRootFile) << "\n";
  std::cout << "reprocTreeName      : " 
            << (reprocTreeName.empty() ? "<default>" : reprocTreeName) << "\n";
  std::cout << "================================================================\n\n";
  EventProcessor processor(mgr, inputFileDir, outputFileDir, IsreprocRootFile, inputRootTreeName, inputRootFileName, nfile, nthreads);
  processor.ProcessEvents();
}
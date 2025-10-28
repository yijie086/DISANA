#include <TROOT.h>
#include <TSystem.h>

#include <ROOT/RDataFrame.hxx>  // NEW
#include <iostream>
#include <thread>  // NEW

#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/EventProcessor.h"
#include "./../DreamAN/core/PhiAnalysis.h"

void RunPhiAnalysis(const std::string& inputDir, int nfile, int nthreads = 0) {
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

  bool IsMC = false;              // Set to true if you want to run on MC data
  bool IsreprocRootFile = false;  // Set to true if you want to reprocess ROOT files
  bool IsInbending = false;       // Set to true if you want to run on inbending data
  bool IsMinimalBook = false;     // Set to true if you want to run minimal booking (only essential histograms)
  // std::string dataconfig = "rgasp18_inb";
  // std::string dataconfig = "rgasp18_outb";
  // std::string dataconfig = "rgafall18_inb";
  // std::string dataconfig = "rgafall18_outb";
  std::string dataconfig = "rgasp19_inb";
  // std::string dataconfig = "rgkfa18_7546";

  if (dataconfig == "rgkfa18_7546") {
    IsInbending = false;  // Set to false for outbending data
  }
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgafall18_inb" || dataconfig == "rgkfa18_7546" || dataconfig == "rgasp19_inb") {
    IsInbending = true;  // Set to true for inbending data
  }

  std::cout << "Running Phi Analysis with configuration: " << dataconfig << std::endl;

  std::string inputFileDir = inputDir;
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA/data_processed/spring2019/inb/DVKpKm_wagon/";  // Default output directory
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/nsidis_wagon/missing_Km_output/";  // Default output
  // directory std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/inb/nsidis/missing_Km_output/";  // Default output
  // directory std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/DVKpKm_wagon/";  // Default output directory
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/nsidis_wagon/missing_Km_output/";  // Default output
  // directory std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/fall2018/sims/DVCS/inb/";  // Default output directory
  //  DVCS wagons
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/inb/DVKpKm_wagon/";  // Default output directory
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/outb/DVKpKm_wagon/";  // Default output directory
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/fall2018/inb/DVKpKm_wagon/";  // Default output directory
  // std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/fall2018/outb/DVKpKm_wagon/";  // Default output directory
  /*std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/DVKpKm_wagon/";  // Default output directory*/
  std::string outputFileDir = "./";  // local output directory
  std::string inputRootFileName = " ";
  std::string inputRootTreeName = " ";
  std::string inputHipoDirTest = " ";

  if (inputDir.empty()) {
    std::cout << "No input custom input dir provided: I have picked up one file from " << inputHipoDirTest << std::endl;
    inputHipoDirTest = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/temp_hipofile/";
    inputFileDir = inputHipoDirTest;
  }

  if (IsreprocRootFile) {
    // inputFileDir = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Inb/";
    // inputFileDir = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis316/DISANA/build/rgk7546dataSFCorr/";
    // inputFileDir = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim/fall2018/nsidis_wagon/inb/";
    inputFileDir = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim/sp2018/nsidis_wagon/inb/";
    // inputFileDir = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim/fall2018/nsidis_wagon/outb/";
    inputRootFileName = "dfSelected.root";
    inputRootTreeName = "dfSelected";
  }

  AnalysisTaskManager mgr;
  mgr.SetOututDir(outputFileDir);

  // fiducial cuts///
  std::shared_ptr<TrackCut> trackCuts = std::make_shared<TrackCut>();
  // defined here the DC edge cuts for electron and proton
  auto edge_regions_e = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edge_regions_p = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edge_regions_kM = std::vector<float>{3.0f, 3.0f, 10.0f};
  auto edge_regions_kP = std::vector<float>{3.0f, 3.0f, 10.0f};

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
  // End of edge fiducial cuts
  // Sector 1, PCal args PID, sector, side, min, max
  trackCuts->AddPCalFiducialRange(11, 1, "lw", 72.0, 94.5);
  trackCuts->AddPCalFiducialRange(11, 1, "lw", 211.5, 234.0);
  if (dataconfig == "rgasp19_inb" || dataconfig == "rgafall18_inb" || dataconfig == "rgafall18_outb") {
    trackCuts->AddPCalFiducialRange(11, 2, "lv", 99.0, 117.5);
  }
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgasp18_outb" || dataconfig == "rgasp19_inb") {
    trackCuts->AddPCalFiducialRange(11, 2, "lv", 31.5, 49.5);
  }
  trackCuts->AddPCalFiducialRange(11, 3, "lv", 346.5, 378.0);
  trackCuts->AddPCalFiducialRange(11, 4, "lv", 0.0, 13.5);
  trackCuts->AddPCalFiducialRange(11, 4, "lv", 229.5, 243.0);
  trackCuts->AddPCalFiducialRange(11, 6, "lw", 166.5, 193.5);
  // Sector 1, ECin only, sector, side, min, max
  trackCuts->AddECinFiducialRange(11, 1, "lv", 67.5, 94.5);
  trackCuts->AddECinFiducialRange(11, 4, "lw", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(11, 5, "lv", 0.0, 23.5);
  trackCuts->AddECinFiducialRange(11, 6, "lw", 0.0, 23.5);
  trackCuts->AddECoutFiducialRange(11, 1, "lv", 0.0, 40.5);
  trackCuts->AddECoutFiducialRange(11, 5, "lv", 193.5, 216.0);

  /// set sampling fraction for the particle in detector
  trackCuts->SetMinECALEnergyCut(11, 1, 0.06);  // Electron PCal layer 1
  // apply sampling fraction and diagolal cuts tbd!!
  /*trackCuts->SetSFCut(true,11,0.19,4.9);  // Set to true if you want to apply sampling fraction cut
  //3 sigma cuts for electron in rga sp18 in
  if (dataconfig == "rgasp18_inb") {
      trackCuts->AddSamplingFractionMinCut(11, 1, 0.160145, 0.0121428, -0.00130927);  // Electronsector 1
      trackCuts->AddSamplingFractionMaxCut(11, 1, 0.288592, 0.00348667, -6.33249e-05);  // Electro sector 1
      trackCuts->AddSamplingFractionMinCut(11, 2, 0.135106, 0.0249842, -0.00270237);  // Electronsector 2
      trackCuts->AddSamplingFractionMaxCut(11, 2, 0.331777, -0.0134885, 0.00144937);  // Electronsector 2
      trackCuts->AddSamplingFractionMinCut(11, 3, 0.135934, 0.0294809, -0.00307425);  // Electronsector 3
      trackCuts->AddSamplingFractionMaxCut(11, 3, 0.324863, -0.0116574, 0.00104596);  // Electro sector 3
      trackCuts->AddSamplingFractionMinCut(11, 4, 0.145492, 0.0211262, -0.00191841);  // Electronsector 4
      trackCuts->AddSamplingFractionMaxCut(11, 4, 0.325025, -0.0139604, 0.00143415);  // Electronsector 4
      trackCuts->AddSamplingFractionMinCut(11, 5, 0.158792, 0.0159181, -0.00139627);  // Electronsector 5
      trackCuts->AddSamplingFractionMaxCut(11, 5, 0.310616, -0.00565024, 0.000488421);  // Electronsector 5
      trackCuts->AddSamplingFractionMinCut(11, 6, 0.154587, 0.0202241, -0.00200259);  // Electronsector 6
      trackCuts->AddSamplingFractionMaxCut(11, 6, 0.315009, -0.00725113, 0.000452379);  // Electro sector 6
  }
  if (dataconfig == "rgasp18_outb") {
      trackCuts->AddSamplingFractionMinCut(11, 1, 0.153109, 0.018688, -0.00156815);  // Electronsector 1
      trackCuts->AddSamplingFractionMaxCut(11, 1, 0.27668, 0.00393487, -0.000495404);  // Electro sector 1
      trackCuts->AddSamplingFractionMinCut(11, 2, 0.156875, 0.0164895, -0.00139281);  // Electronsector 2
      trackCuts->AddSamplingFractionMaxCut(11, 2, 0.285319, 0.000615206, 7.59911e-06);  // Electronsector 2
      trackCuts->AddSamplingFractionMinCut(11, 3, 0.159718, 0.0176749, -0.00147327);  // Electronsector 3
      trackCuts->AddSamplingFractionMaxCut(11, 3, 0.280747, 0.00494729, -0.00071406);  // Electro sector 3
      trackCuts->AddSamplingFractionMinCut(11, 4, 0.158643, 0.0173334, -0.00137467);  // Electronsector 4
      trackCuts->AddSamplingFractionMaxCut(11, 4, 0.279673, 0.00510175, -0.00068517);  // Electronsector 4
      trackCuts->AddSamplingFractionMinCut(11, 5, 0.155656, 0.0166533, -0.0012774);  // Electronsector 5
      trackCuts->AddSamplingFractionMaxCut(11, 5, 0.285419, 1.77522e-05, -8.01094e-05);  // Electronsector 5
      trackCuts->AddSamplingFractionMinCut(11, 6, 0.151147, 0.0212032, -0.00181402);   // Electronsector6
      trackCuts->AddSamplingFractionMaxCut(11, 6, 0.281681, 0.00402848, -0.000666406);  // Electro sector 6
  }
  // rga fall18 in
  if (dataconfig == "rgafall18_inb") {
      trackCuts->AddSamplingFractionMinCut(11, 1, 0.160145, 0.0121428, -0.00130927);  // Electronsector 1
      trackCuts->AddSamplingFractionMaxCut(11, 1, 0.288592, 0.00348667, -6.33249e-05);  // Electro sector 1
      trackCuts->AddSamplingFractionMinCut(11, 2, 0.135106, 0.0249842, -0.00270237);  // Electronsector 2
      trackCuts->AddSamplingFractionMaxCut(11, 2, 0.331777, -0.0134885, 0.00144937);  // Electronsector 2
      trackCuts->AddSamplingFractionMinCut(11, 3, 0.135934, 0.0294809, -0.00307425);  // Electronsector 3
      trackCuts->AddSamplingFractionMaxCut(11, 3, 0.324863, -0.0116574, 0.00104596);  // Electro sector 3
      trackCuts->AddSamplingFractionMinCut(11, 4, 0.145492, 0.0211262, -0.00191841);  // Electronsector 4
      trackCuts->AddSamplingFractionMaxCut(11, 4, 0.325025, -0.0139604, 0.00143415);  // Electronsector 4
      trackCuts->AddSamplingFractionMinCut(11, 5, 0.158792, 0.0159181, -0.00139627);  // Electronsector 5
      trackCuts->AddSamplingFractionMaxCut(11, 5, 0.310616, -0.00565024, 0.000488421);  // Electronsector 5
      trackCuts->AddSamplingFractionMinCut(11, 6, 0.154587, 0.0202241, -0.00200259);  // Electronsector 6
      trackCuts->AddSamplingFractionMaxCut(11, 6, 0.315009, -0.00725113, 0.000452379);  // Electro sector 6
  }

  //fall18 out
  if (dataconfig == "rgafall18_outb") {
      trackCuts->AddSamplingFractionMinCut(11, 1, 0.153109, 0.018688, -0.00156815);  // Electronsector 1
      trackCuts->AddSamplingFractionMaxCut(11, 1, 0.27668, 0.00393487, -0.000495404);  // Electro sector 1
      trackCuts->AddSamplingFractionMinCut(11, 2, 0.156875, 0.0164895, -0.00139281);  // Electronsector 2
      trackCuts->AddSamplingFractionMaxCut(11, 2, 0.285319, 0.000615206, 7.59911e-06);  // Electronsector 2
      trackCuts->AddSamplingFractionMinCut(11, 3, 0.159718, 0.0176749, -0.00147327);  // Electronsector 3
      trackCuts->AddSamplingFractionMaxCut(11, 3, 0.280747, 0.00494729, -0.00071406);  // Electro sector 3
      trackCuts->AddSamplingFractionMinCut(11, 4, 0.158643, 0.0173334, -0.00137467);  // Electronsector 4
      trackCuts->AddSamplingFractionMaxCut(11, 4, 0.279673, 0.00510175, -0.00068517);  // Electronsector 4
      trackCuts->AddSamplingFractionMinCut(11, 5, 0.155656, 0.0166533, -0.0012774);  // Electronsector 5
      trackCuts->AddSamplingFractionMaxCut(11, 5, 0.285419, 1.77522e-05, -8.01094e-05);  // Electronsector 5
      trackCuts->AddSamplingFractionMinCut(11, 6, 0.151147, 0.0212032, -0.00181402);   // Electronsector6
      trackCuts->AddSamplingFractionMaxCut(11, 6, 0.281681, 0.00402848, -0.000666406);  // Electro sector 6
  }
      //// //rga sp19 in
  if (dataconfig == "rgasp19_inb") {
      trackCuts->AddSamplingFractionMinCut(11, 1, 0.160145, 0.0121428, -0.00130927);  // Electronsector 1
      trackCuts->AddSamplingFractionMaxCut(11, 1, 0.288592, 0.00348667, -6.33249e-05);  // Electro sector 1
      trackCuts->AddSamplingFractionMinCut(11, 2, 0.135106, 0.0249842, -0.00270237);  // Electronsector 2
      trackCuts->AddSamplingFractionMaxCut(11, 2, 0.331777, -0.0134885, 0.00144937);  // Electronsector 2
      trackCuts->AddSamplingFractionMinCut(11, 3, 0.135934, 0.0294809, -0.00307425);  // Electronsector 3
      trackCuts->AddSamplingFractionMaxCut(11, 3, 0.324863, -0.0116574, 0.00104596);  // Electro sector 3
      trackCuts->AddSamplingFractionMinCut(11, 4, 0.145492, 0.0211262, -0.00191841);  // Electronsector 4
      trackCuts->AddSamplingFractionMaxCut(11, 4, 0.325025, -0.0139604, 0.00143415);  // Electronsector 4
      trackCuts->AddSamplingFractionMinCut(11, 5, 0.158792, 0.0159181, -0.00139627);  // Electronsector 5
      trackCuts->AddSamplingFractionMaxCut(11, 5, 0.310616, -0.00565024, 0.000488421);  // Electronsector 5
      trackCuts->AddSamplingFractionMinCut(11, 6, 0.154587, 0.0202241, -0.00200259);  // Electronsector 6
      trackCuts->AddSamplingFractionMaxCut(11, 6, 0.315009, -0.00725113, 0.000452379);  // Electro sector 6
  }

*/
  auto corr = std::make_shared<MomentumCorrection>();

  /// Momentum correction for the proton
  if (dataconfig == "rgasp18_inb") {
    corr->AddPiecewiseCorrection(  // Momentum correction for proton RGA sp18 inb
        2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD}, [](double p, double theta, double phi) {
          theta = theta * 180.0 / M_PI;  // Convert theta to degrees
          float A_p = -0.229055 + 0.00924571 * theta - 9.09927e-05 * theta * theta;
          float B_p = 0.371002 - 0.0146818 * theta + 0.000146548 * theta * theta;
          float C_p = -0.174565 + 0.00680452 * theta - 6.9e-05 * theta * theta;
          return p + (A_p + B_p * p + C_p * p * p);
        });
    corr->AddPiecewiseCorrection(  // Momentum correction for proton RGA sp18 inb
        2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD}, [](double p, double theta, double phi) {
          theta = theta * 180.0 / M_PI;  // Convert theta to degrees
          float A_p = 0.0146275 - 0.00124929 * theta + 3.64154e-05 * theta * theta;
          float B_p = -0.00743169 + 0.000458648 * theta - 6.45703e-06 * theta * theta;
          float C_p = 0.0175282 - 0.00128554 * theta + 3.5249e-05 * theta * theta;

          return p + (A_p + B_p / p + C_p / (p * p));
        });
  }
  if (dataconfig == "rgasp18_outb") {
    corr->AddPiecewiseCorrection(  // Momentum correction for proton RGK sp18 out
        2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD}, [](double p, double theta, double phi) {
          theta = theta * 180.0 / M_PI;  // Convert theta to degrees
          float A_p = -0.204359 + 0.00857339 * theta - 8.79867e-05 * theta * theta;
          float B_p = 0.402543 - 0.0168624 * theta + 0.000178539 * theta * theta;
          float C_p = -0.217865 + 0.00908787 * theta - 9.77617e-05 * theta * theta;
          return p + (A_p + B_p * p + C_p * p * p);
        });
    corr->AddPiecewiseCorrection(  // Momentum correction for proton RGA sp18 out
        2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::FD}, [](double p, double theta, double phi) {
          theta = theta * 180.0 / M_PI;  // Convert theta to degrees
          float A_p = 0.00523188 - 9.43614e-05 * theta;
          float B_p = -0.00887291 + 0.000759277 * theta;
          float C_p = 0;

          return p + (A_p + B_p / p + C_p / (p * p));
        });
  }
  // ---------- Fall 2018 (inb) ----------
  if (dataconfig == "rgafall18_inb") {
    // CD branch (central detector)
    corr->AddPiecewiseCorrection(2212, {0.0, 10.0, 0.0 * M_PI / 180, 180.0 * M_PI / 180, 0.0 * M_PI / 180, 360.0 * M_PI / 180, MomentumCorrection::CD},
                                 [](double p, double theta, double /*phi*/) {
                                   theta = theta * 180.0 / M_PI;  // degrees
                                   double A_p = -0.2383991 + 0.0124992 * theta - 0.0001646 * theta * theta;
                                   double B_p = 0.60123885 - 0.03128464 * theta + 0.00041314 * theta * theta;
                                   double C_p = -0.44080146 + 0.02209857 * theta - 0.00028224 * theta * theta;
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
  eventCuts->AddParticleCut("Neg Kaon", kMinus);    // Applies defaults automatically
  eventCuts->AddParticleCut("Pos Kaon", kPos);      // Applies defaults automatically
  // eventCuts->AddParticleMotherCut("phi", phi);      // Applies defaults automatically

  // Task
  auto PhiTask = std::make_unique<PhiAnalysis>(IsMC, IsreprocRootFile, IsMinimalBook);
  PhiTask->SetTrackCuts(trackCuts);
  PhiTask->SetEventCuts(eventCuts);
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgasp18_outb" || dataconfig == "rgafall18_inb" || dataconfig == "rgafall18_outb" || dataconfig == "rgasp19_inb" ||
      dataconfig == "rgasp19_outb") {
    PhiTask->SetBeamEnergy(10.6);
  }

  PhiTask->SetFTonConfig(true);  // Set to true if you have FT (eq. RGK Fall2018 Pass2 6.535GeV is FT-off)
  PhiTask->SetDoFiducialCut(true);

  PhiTask->SetDoInvMassCut(true);  // in this case pi0 background for two-photon pairs in the event
  // PhiTask->SetDoMomentumCorrection(true);  // Set to true if you want to apply momentum correction
  // PhiTask->SetMomentumCorrection(corr);  // Set the momentum correction object
  PhiTask->SetMaxEvents(0);  // Set the maximum number of events to process, 0 means no limit

  mgr.AddTask(std::move(PhiTask));
  // Processor
  EventProcessor processor(mgr, inputFileDir, outputFileDir, IsreprocRootFile, inputRootTreeName, inputRootFileName, nfile, nthreads);
  processor.ProcessEvents();
}

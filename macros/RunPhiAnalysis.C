#include <iostream>

#include "./../DreamAN/core/AnalysisTaskManager.h"
#include "./../DreamAN/core/PhiAnalysis.h"
#include "./../DreamAN/core/EventProcessor.h"

void RunPhiAnalysis(const std::string& inputDir, int nfile) {
  bool IsMC = false;              // Set to true if you want to run on MC data
  bool IsreprocRootFile = false;  // Set to true if you want to reprocess ROOT files
  bool IsInbending = true;        // Set to true if you want to run on inbending data

  //std::string dataconfig = "rgasp18_inb";
  //std::string dataconfig = "rgasp18_outb";
  //std::string dataconfig = "rgafall18_inb";
  //std::string dataconfig = "rgafall18_outb";
  std::string dataconfig = "rgasp19_inb";
 
  

  if (dataconfig == "rgkfa18_7546") {
    IsInbending = false;       // Set to false for outbending data
  }
  if (dataconfig == "rgasp18_inb"|| dataconfig == "rgafall18_inb"|| dataconfig == "rgkfa18_7546"||dataconfig == "rgasp19_inb") {
    IsInbending = true;        // Set to true for inbending data
  }

  std::cout << "Running DVCS Analysis with configuration: " << dataconfig << std::endl;

  std::string inputFileDir = inputDir;
  std::string outputFileDir = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA/data_processed/spring2019/inb/DVKpKm_wagon/";  // Default output directory
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
    //inputFileDir = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Inb/";
    inputFileDir = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis316/DISANA/build/rgk7546dataSFCorr/";
    // inputFileDir = "./";
    inputRootFileName = "dfSelected.root";
    inputRootTreeName = "dfSelected";
  }

  AnalysisTaskManager mgr;
  mgr.SetOututDir(outputFileDir);
  
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
  auto edge_regions_kM = std::vector<float>{3.0f, 3.0f, 5.0f};
  auto edge_regions_kP = std::vector<float>{3.0f, 3.0f, 5.0f};
  // auto CVT_edge_layers_p = std::vector<float>{-100.0f, -100.0f, -100.0f, -100.0f, -100.0f};  // CVT edge cuts for test
  auto CVT_edge_layers_p = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for protons
  auto CVT_edge_layers_kM = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for protons
  auto CVT_edge_layers_kP = std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};  // CVT edge cuts for protons

  trackCuts->SetDCEdgeCuts(11, edge_regions_e);        // DC edge cuts for electrons
  trackCuts->SetDCEdgeCuts(2212, edge_regions_p);      // DC edge cuts for protons
  trackCuts->SetDCEdgeCuts(-321, edge_regions_p);      // DC edge cuts for kM
  trackCuts->SetDCEdgeCuts(321, edge_regions_p);      // DC edge cuts for kP
  trackCuts->SetCVTEdgeCuts(-321, CVT_edge_layers_p);  // CVT edge cuts for proton
  trackCuts->SetCVTEdgeCuts(-321, CVT_edge_layers_kM);  // CVT edge cuts for Kaon M
  trackCuts->SetCVTEdgeCuts(321, CVT_edge_layers_kP);  // CVT edge cuts for KaonP

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
  /*trackCuts->AddPCalFiducialRange(321, 1, "lw", 0.0, 13.5);
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
  trackCuts->AddPCalFiducialRange(22, 6, "lv", 0.0, 13.5);*/
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
  /*trackCuts->AddPCalFiducialRange(22, 1, "lw", 72.0, 94.5);
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
*/
  /// set sampling fraction for the particle in detector
  trackCuts->SetMinECALEnergyCut(11, 1, 0.06);  // Electron PCal layer 1
  // apply sampling fraction and diagolal cuts tbd!!
  trackCuts->SetSFCut(true,11,0.19,4.9);  // Set to true if you want to apply sampling fraction cut
  
  //rga sp18 in
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
   ////   //rga sp19 out
  if (dataconfig == "rgasp19_outb") {
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


  // particles for the reaction DVCS: e, p, and gamma
  EventCut* eventCuts = new EventCut();


  ParticleCut proton;
  if (dataconfig == "rgasp18_inb"|| dataconfig == "rgasp19_inb"|| dataconfig == "rgafall18_inb"|| dataconfig == "rgasp18_outb"|| dataconfig == "rgafall18_outb"|| dataconfig == "rgasp19_outb") {
    proton.pid = 2212;                              // Proton PID
    proton.charge = 1;                              // Proton charge
    proton.minCount = 1;                            // Minimum count of protons
    proton.maxCount = 1;                            // Maximum count of protons
    proton.minCDMomentum = 0.3f;                    // Minimum momentum for protons
    if (IsInbending) proton.minFDMomentum = 0.42f;  // Minimum momentum for protons in FD Inbending
    if (!IsInbending) proton.minFDMomentum = 0.5f;  // Minimum momentum for protons in FD Outbending
    //proton.maxCDMomentum = 1.2f;                    // Maximum momentum for protons
    //proton.maxFDMomentum = 1.2f;                    // Maximum momentum for protons
    //proton.minTheta = 0.0f;                         // Minimum theta for protons
    //proton.maxTheta = 64.23 * M_PI / 180.0;         // Maximum theta for protons (approximately 64.23 degrees)
    proton.minPhi = 0.0f;                           // Minimum phi for protons
    proton.maxPhi = 2.0f * M_PI;
  }

  ParticleCut electron;
  electron.pid = 11;              // Electron PID
  electron.charge = -1;           // Electron charge
  electron.minCount = 1;          // Minimum count of electrons
  electron.maxCount = 1;          // Maximum count of electrons
  electron.minFDMomentum = 2.0f;  // Minimum momentum for electrons

  ParticleCut kMinus;
  kMinus.pid = -321;      // KaonM PID
  kMinus.charge = -1;    // charge
  kMinus.minCount = 1;  // Minimum count of Negative kaons

  ParticleCut kPos;
  kPos.pid = 321;      // KaonPos PID
  kPos.charge = 1;    
  kPos.minCount = 1;  // Minimum count of possible kaons
  
  // pi0 background studies, similar is valid for exclusive meson production
  TwoBodyMotherCut phi;
  phi.expectedMotherMass = 1.019f;  // determined from the fit inv mass of two photon
  phi.massSigma = 0.5f;          // determined from the fit inv mass of two photon
  phi.nSigmaMass = 100;             // choice for nsigma is user dependent

  eventCuts->AddParticleCut("proton", proton);      // Applies defaults automatically
  eventCuts->AddParticleCut("electron", electron);  // Applies defaults automatically
  eventCuts->AddParticleCut("Neg Kaon", kMinus);      // Applies defaults automatically
  eventCuts->AddParticleCut("Pos Kaon", kPos);      // Applies defaults automatically
  eventCuts->AddParticleMotherCut("phi", phi);      // Applies defaults automatically

  auto corr = std::make_shared<MomentumCorrection>();

  // Task
  auto PhiTask = std::make_unique<PhiAnalysis>(IsMC, IsreprocRootFile);
  PhiTask->SetTrackCuts(trackCuts);
  PhiTask->SetEventCuts(eventCuts);
  if (dataconfig == "rgasp18_inb" || dataconfig == "rgasp18_outb"|| dataconfig == "rgafall18_inb"|| dataconfig == "rgafall18_outb"||dataconfig == "rgasp19_inb"|| dataconfig == "rgasp19_outb") {
    PhiTask->SetBeamEnergy(10.6);
  }

  PhiTask->SetFTonConfig(true);  // Set to true if you have FT (eq. RGK Fall2018 Pass2 6.535GeV is FT-off)
  PhiTask->SetDoFiducialCut(true);

  PhiTask->SetDoInvMassCut(true); // in this case pi0 background for two-photon pairs in the event
 // PhiTask->SetDoMomentumCorrection(true);  // Set to true if you want to apply momentum correction
  //PhiTask->SetMomentumCorrection(corr);  // Set the momentum correction object
  PhiTask->SetMaxEvents(0);  // Set the maximum number of events to process, 0 means no limit


  mgr.AddTask(std::move(PhiTask));

  // Processor
  EventProcessor processor(mgr, inputFileDir, IsreprocRootFile, inputRootTreeName, inputRootFileName, nfile);
  processor.ProcessEvents();
}

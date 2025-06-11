#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>
#include <map>

// Function to draw sector lines (for visual aid)
void DrawSectorLines() {
  const float rMax = 400.0;

  TLine *vline = new TLine(0, -rMax, 0, rMax);
  vline->SetLineColor(kRed);
  vline->Draw("same");

  const float angles[] = {30, 90, 150, 210, 270, 330};
  for (float angleDeg : angles) {
    float angleRad = angleDeg * TMath::Pi() / 180.0;
    float x = rMax * std::cos(angleRad);
    float y = rMax * std::sin(angleRad);
    TLine *l = new TLine(0, 0, x, y);
    l->SetLineColor(kRed);
    l->Draw("same");
  }
}

void DrawPCALFromFile(
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // 只处理R1
  auto dfPCALS1 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS2 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});
  
  auto dfPCALS3 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS4 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS5 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS6 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});


  // 只画R1
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  auto hs1 = dfPCALS1.Histo2D({"hs1", "ECAL PCAL sector 1;x [cm];y [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs2 = dfPCALS2.Histo2D({"hs2", "ECAL PCAL sector 2;x [cm];y [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs3 = dfPCALS3.Histo2D({"hs3", "ECAL PCAL sector 3;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs4 = dfPCALS4.Histo2D({"hs4", "ECAL PCAL sector 4;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs5 = dfPCALS5.Histo2D({"hs5", "ECAL PCAL sector 5;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs6 = dfPCALS6.Histo2D({"hs6", "ECAL PCAL sector 6;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");

  TCanvas *c1 = new TCanvas("c1", "ECAL PCAL", 6000, 4000);
  c1->Divide(3, 2);
  gStyle->SetOptStat(0);

  c1->cd(1);
  gPad->SetLogz();
  hs1->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  hs2->Draw("colz");
  c1->cd(3);
  gPad->SetLogz();
  hs3->Draw("colz");
  c1->cd(4);
  gPad->SetLogz();
  hs4->Draw("colz");
  c1->cd(5);
  gPad->SetLogz();
  hs5->Draw("colz");
  c1->cd(6);
  gPad->SetLogz();
  hs6->Draw("colz");

  std::string outname = "Calorimeter_PCAL_" + treename + ".png";
  c1->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}

void DrawECINFromFile(
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // 只处理ECIN（layer==4）
  auto dfECINS1 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS2 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});
  
  auto dfECINS3 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS4 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS5 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS6 =
      df.Define("xlw",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});


  // 只画ECIN
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  auto hs1 = dfECINS1.Histo2D({"hs1", "ECAL ECIN sector 1;x [cm];y [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs2 = dfECINS2.Histo2D({"hs2", "ECAL ECIN sector 2;x [cm];y [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs3 = dfECINS3.Histo2D({"hs3", "ECAL ECIN sector 3;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs4 = dfECINS4.Histo2D({"hs4", "ECAL ECIN sector 4;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs5 = dfECINS5.Histo2D({"hs5", "ECAL ECIN sector 5;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs6 = dfECINS6.Histo2D({"hs6", "ECAL ECIN sector 6;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");

  TCanvas *c1 = new TCanvas("c1", "ECAL ECIN", 6000, 4000);
  c1->Divide(3, 2);
  gStyle->SetOptStat(0);

  c1->cd(1);
  gPad->SetLogz();
  hs1->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  hs2->Draw("colz");
  c1->cd(3);
  gPad->SetLogz();
  hs3->Draw("colz");
  c1->cd(4);
  gPad->SetLogz();
  hs4->Draw("colz");
  c1->cd(5);
  gPad->SetLogz();
  hs5->Draw("colz");
  c1->cd(6);
  gPad->SetLogz();
  hs6->Draw("colz");

  std::string outname = "Calorimeter_ECIN_" + treename + ".png";
  c1->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}

void DrawECOUTFromFile(
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // 只处理ECOUT（layer==7），画xlu/ylv
  auto dfECOUTS1 =
      df.Define("xlu",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lu,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 7 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(lu[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 7 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  // 2~6 sector
  auto dfECOUTS2 = df.Define("xlu", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS3 = df.Define("xlu", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS4 = df.Define("xlu", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS5 = df.Define("xlu", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS6 = df.Define("xlu", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  // 只画ECOUT
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  auto hs1 = dfECOUTS1.Histo2D({"hs1", "ECAL ECOUT sector 1;x [cm];y [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs2 = dfECOUTS2.Histo2D({"hs2", "ECAL ECOUT sector 2;x [cm];y [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs3 = dfECOUTS3.Histo2D({"hs3", "ECAL ECOUT sector 3;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs4 = dfECOUTS4.Histo2D({"hs4", "ECAL ECOUT sector 4;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs5 = dfECOUTS5.Histo2D({"hs5", "ECAL ECOUT sector 5;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs6 = dfECOUTS6.Histo2D({"hs6", "ECAL ECOUT sector 6;x [cm];y [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");

  TCanvas *c1 = new TCanvas("c1", "ECAL ECOUT", 6000, 4000);
  c1->Divide(3, 2);
  gStyle->SetOptStat(0);

  c1->cd(1);
  gPad->SetLogz();
  hs1->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  hs2->Draw("colz");
  c1->cd(3);
  gPad->SetLogz();
  hs3->Draw("colz");
  c1->cd(4);
  gPad->SetLogz();
  hs4->Draw("colz");
  c1->cd(5);
  gPad->SetLogz();
  hs5->Draw("colz");
  c1->cd(6);
  gPad->SetLogz();
  hs6->Draw("colz");

  std::string outname = "Calorimeter_ECOUT_" + treename + ".png";
  c1->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}


void plotECALFids() {

  //std::string input_path_from_analysisRun = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/CheckWithInclusiveData_electron_photon/";
  std::string input_path_from_analysisRun = "./../Build/";
  std::string filename_after = Form("%s/dfSelected_after_fiducialCuts.root", input_path_from_analysisRun.c_str());
  std::string filename_before = Form("%s/dfSelected_before_fiducialCuts.root", input_path_from_analysisRun.c_str());

  DrawPCALFromFile(filename_after, "dfSelected_after");
  DrawPCALFromFile(filename_before, "dfSelected_before");
  DrawECINFromFile(filename_after, "dfSelected_after");
  DrawECINFromFile(filename_before, "dfSelected_before");
  DrawECOUTFromFile(filename_after, "dfSelected_after");
  DrawECOUTFromFile(filename_before, "dfSelected_before");
  gApplication->Terminate(0);
}
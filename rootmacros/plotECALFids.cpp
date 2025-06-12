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
    const int &selectedPid,
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // 只处理R1
  auto dfPCALS1 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS2 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});
  
  auto dfPCALS3 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS4 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS5 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfPCALS6 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 1 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 1 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});


  // 只画R1
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  auto hs1 = dfPCALS1.Histo2D({"hs1", ";w [cm];v [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs2 = dfPCALS2.Histo2D({"hs2", ";w [cm];v [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs3 = dfPCALS3.Histo2D({"hs3", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs4 = dfPCALS4.Histo2D({"hs4", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs5 = dfPCALS5.Histo2D({"hs5", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs6 = dfPCALS6.Histo2D({"hs6", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");

  TCanvas *c1 = new TCanvas("c1", "ECAL PCAL", 6000, 4000);
  c1->Divide(3, 2);
  gStyle->SetOptStat(0);

  c1->cd(1);
  gPad->SetLogz();
  hs1->Draw("colz");
  if (selectedPid == 11) {
    hs1->SetTitle("electron PCAL Sector 1");
  }
  else if (selectedPid == 22) {
    hs1->SetTitle("photon PCAL Sector 1");
  } else {
    hs1->SetTitle("PCAL Sector 1");
  }
  c1->cd(2);
  gPad->SetLogz();
  hs2->Draw("colz");
  if (selectedPid == 11) {
    hs2->SetTitle("electron PCAL Sector 2");
  }
  else if (selectedPid == 22) {
    hs2->SetTitle("photon PCAL Sector 2");
  } else {
    hs2->SetTitle("PCAL Sector 2");
  }
  c1->cd(3);
  gPad->SetLogz();
  hs3->Draw("colz");
  if (selectedPid == 11) {
    hs3->SetTitle("electron PCAL Sector 3");
  }
  else if (selectedPid == 22) {
    hs3->SetTitle("photon PCAL Sector 3");
  } else {
    hs3->SetTitle("PCAL Sector 3");
  }
  c1->cd(4);
  gPad->SetLogz();
  hs4->Draw("colz");
  if (selectedPid == 11) {
    hs4->SetTitle("electron PCAL Sector 4");
  }
  else if (selectedPid == 22) {
    hs4->SetTitle("photon PCAL Sector 4");
  } else {
    hs4->SetTitle("PCAL Sector 4");
  }
  c1->cd(5);
  gPad->SetLogz();
  hs5->Draw("colz");
  if (selectedPid == 11) {
    hs5->SetTitle("electron PCAL Sector 5");
  }
  else if (selectedPid == 22) {
    hs5->SetTitle("photon PCAL Sector 5");
  } else {
    hs5->SetTitle("PCAL Sector 5");
  }
  c1->cd(6);
  gPad->SetLogz();
  hs6->Draw("colz");
  if (selectedPid == 11) {
    hs6->SetTitle("electron PCAL Sector 6");
  }
  else if (selectedPid == 22) {
    hs6->SetTitle("photon PCAL Sector 6");
  } else {
    hs6->SetTitle("PCAL Sector 6");
  }

  std::string outname = "Calorimeter_PCAL_" + treename + ".png";
  if (selectedPid == 11) {
    outname = "Calorimeter_PCAL_electron_" + treename + ".png";
  } else if (selectedPid == 22) {
    outname = "Calorimeter_PCAL_photon_" + treename + ".png";
  }
  c1->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}

void DrawECINFromFile(
    const int &selectedPid,
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // 只处理ECIN（layer==4）
  auto dfECINS1 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS2 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});
  
  auto dfECINS3 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS4 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS5 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECINS6 =
      df.Define("xlw",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lw,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 4 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lw[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lw", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 4 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});


  // 只画ECIN
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  auto hs1 = dfECINS1.Histo2D({"hs1", ";w [cm];v [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs2 = dfECINS2.Histo2D({"hs2", ";w [cm];v [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs3 = dfECINS3.Histo2D({"hs3", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs4 = dfECINS4.Histo2D({"hs4", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs5 = dfECINS5.Histo2D({"hs5", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");
  auto hs6 = dfECINS6.Histo2D({"hs6", ";w [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlw", "ylv");

  TCanvas *c1 = new TCanvas("c1", "ECAL ECIN", 6000, 4000);
  c1->Divide(3, 2);
  gStyle->SetOptStat(0);

  c1->cd(1);
  gPad->SetLogz();
  hs1->Draw("colz");
  if (selectedPid == 11) {
    hs1->SetTitle("electron ECIN Sector 1");
  }
  else if (selectedPid == 22) {
    hs1->SetTitle("photon ECIN Sector 1");
  } else {
    hs1->SetTitle("ECIN Sector 1");
  }
  c1->cd(2);
  gPad->SetLogz();
  hs2->Draw("colz");
  if (selectedPid == 11) {
    hs2->SetTitle("electron ECIN Sector 2");
  }
  else if (selectedPid == 22) {
    hs2->SetTitle("photon ECIN Sector 2");
  } else {
    hs2->SetTitle("ECIN Sector 2");
  }
  c1->cd(3);
  gPad->SetLogz();
  hs3->Draw("colz");
  if (selectedPid == 11) {
    hs3->SetTitle("electron ECIN Sector 3");
  }
  else if (selectedPid == 22) {
    hs3->SetTitle("photon ECIN Sector 3");
  } else {
    hs3->SetTitle("ECIN Sector 3");
  }
  c1->cd(4);
  gPad->SetLogz();
  hs4->Draw("colz");
  if (selectedPid == 11) {
    hs4->SetTitle("electron ECIN Sector 4");
  }
  else if (selectedPid == 22) {
    hs4->SetTitle("photon ECIN Sector 4");
  } else {
    hs4->SetTitle("ECIN Sector 4");
  }
  c1->cd(5);
  gPad->SetLogz();
  hs5->Draw("colz");
  if (selectedPid == 11) {
    hs5->SetTitle("electron ECIN Sector 5");
  }
  else if (selectedPid == 22) {
    hs5->SetTitle("photon ECIN Sector 5");
  } else {
    hs5->SetTitle("ECIN Sector 5");
  }
  c1->cd(6);
  gPad->SetLogz();
  hs6->Draw("colz");
  if (selectedPid == 11) {
    hs6->SetTitle("electron ECIN Sector 6");
  }
  else if (selectedPid == 22) {
    hs6->SetTitle("photon ECIN Sector 6");
  } else {
    hs6->SetTitle("ECIN Sector 6");
  }

  std::string outname = "Calorimeter_ECIN_" + treename + ".png";
  if (selectedPid == 11) {
    outname = "Calorimeter_ECIN_electron_" + treename + ".png";
  } else if (selectedPid == 22) {
    outname = "Calorimeter_ECIN_photon_" + treename + ".png";
  }
  c1->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}

void DrawECOUTFromFile(
    const int &selectedPid,
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // 只处理ECOUT（layer==7），画xlu/ylv
  auto dfECOUTS1 =
      df.Define("xlu",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<short> &sector,
                   const ROOT::VecOps::RVec<float> &lu,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 7 && layer[i] == 7 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(lu[i]);
                  }
                  return out;
                },
                {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
          .Define("ylv",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<short> &sector,
                     const ROOT::VecOps::RVec<float> &lv,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      if (det[i] == 7 && layer[i] == 7 && sector[i] == 1 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(lv[i]);
                    return out;
                  },
                  {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  // 2~6 sector
  auto dfECOUTS2 = df.Define("xlu", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 2 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS3 = df.Define("xlu", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 3 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS4 = df.Define("xlu", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 4 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS5 = df.Define("xlu", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 5 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  auto dfECOUTS6 = df.Define("xlu", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lu, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) { if (det[i] == 7 && layer[i] == 7 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lu[i]); } return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lu", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"})
    .Define("ylv", [selectedPid](const ROOT::VecOps::RVec<short> &det, const ROOT::VecOps::RVec<short> &layer, const ROOT::VecOps::RVec<short> &sector, const ROOT::VecOps::RVec<float> &lv, const ROOT::VecOps::RVec<float> &p, const ROOT::VecOps::RVec<int16_t> &pindex, const ROOT::VecOps::RVec<int> &pid) {
    ROOT::VecOps::RVec<float> out; for (size_t i = 0; i < layer.size(); ++i) if (det[i] == 7 && layer[i] == 7 && sector[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid) out.push_back(lv[i]); return out;
  }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector", "REC_Calorimeter_lv", "REC_Particle_p", "REC_Calorimeter_pindex","REC_Particle_pid"});

  // 只画ECOUT
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  auto hs1 = dfECOUTS1.Histo2D({"hs1", ";u [cm];v [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs2 = dfECOUTS2.Histo2D({"hs2", ";u [cm];v [cm]", 
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs3 = dfECOUTS3.Histo2D({"hs3", ";u [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs4 = dfECOUTS4.Histo2D({"hs4", ";u [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs5 = dfECOUTS5.Histo2D({"hs5", ";u [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");
  auto hs6 = dfECOUTS6.Histo2D({"hs6", ";u [cm];v [cm]",
                          1000, 0, 500, 1000, 0, 500},
                          "xlu", "ylv");

  TCanvas *c1 = new TCanvas("c1", "ECAL ECOUT", 6000, 4000);
  c1->Divide(3, 2);
  gStyle->SetOptStat(0);

  c1->cd(1);
  gPad->SetLogz();
  hs1->Draw("colz");
  if (selectedPid == 11) {
    hs1->SetTitle("electron ECOUT Sector 1");
  }
  else if (selectedPid == 22) {
    hs1->SetTitle("photon ECOUT Sector 1");
  } else {
    hs1->SetTitle("ECOUT Sector 1");
  }
  c1->cd(2);
  gPad->SetLogz();
  hs2->Draw("colz");
  if (selectedPid == 11) {
    hs2->SetTitle("electron ECOUT Sector 2");
  }
  else if (selectedPid == 22) {
    hs2->SetTitle("photon ECOUT Sector 2");
  } else {
    hs2->SetTitle("ECOUT Sector 2");
  }
  c1->cd(3);
  gPad->SetLogz();
  hs3->Draw("colz");
  if (selectedPid == 11) {
    hs3->SetTitle("electron ECOUT Sector 3");
  }
  else if (selectedPid == 22) {
    hs3->SetTitle("photon ECOUT Sector 3");
  } else {
    hs3->SetTitle("ECOUT Sector 3");
  }
  c1->cd(4);
  gPad->SetLogz();
  hs4->Draw("colz");
  if (selectedPid == 11) {
    hs4->SetTitle("electron ECOUT Sector 4");
  }
  else if (selectedPid == 22) {
    hs4->SetTitle("photon ECOUT Sector 4");
  } else {
    hs4->SetTitle("ECOUT Sector 4");
  }
  c1->cd(5);
  gPad->SetLogz();
  hs5->Draw("colz");
  if (selectedPid == 11) {
    hs5->SetTitle("electron ECOUT Sector 5");
  }
  else if (selectedPid == 22) {
    hs5->SetTitle("photon ECOUT Sector 5");
  } else {
    hs5->SetTitle("ECOUT Sector 5");
  }
  c1->cd(6);
  gPad->SetLogz();
  hs6->Draw("colz");
  if (selectedPid == 11) {
    hs6->SetTitle("electron ECOUT Sector 6");
  }
  else if (selectedPid == 22) {
    hs6->SetTitle("photon ECOUT Sector 6");
  } else {
    hs6->SetTitle("ECOUT Sector 6");
  }

  std::string outname = "Calorimeter_ECOUT_" + treename + ".png";
  if (selectedPid == 11) {
    outname = "Calorimeter_ECOUT_electron_" + treename + ".png";
  } else if (selectedPid == 22) {
    outname = "Calorimeter_ECOUT_photon_" + treename + ".png";
  }
  c1->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}


void plotECALFids() {

  //std::string input_path_from_analysisRun = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/CheckWithInclusiveData_electron_photon/";
  std::string input_path_from_analysisRun = "./../build/";
  std::string filename_after = Form("%s/dfSelected_after_fiducialCuts.root", input_path_from_analysisRun.c_str());
  std::string filename_before = Form("%s/dfSelected_before_fiducialCuts.root", input_path_from_analysisRun.c_str());

  DrawPCALFromFile(11,filename_after, "dfSelected_after");
  DrawPCALFromFile(11,filename_before, "dfSelected_before");
  DrawECINFromFile(11,filename_after, "dfSelected_after");
  DrawECINFromFile(11,filename_before, "dfSelected_before");
  DrawECOUTFromFile(11,filename_after, "dfSelected_after");
  DrawECOUTFromFile(11,filename_before, "dfSelected_before");
  DrawPCALFromFile(22,filename_after, "dfSelected_after");
  DrawPCALFromFile(22,filename_before, "dfSelected_before");
  DrawECINFromFile(22,filename_after, "dfSelected_after");
  DrawECINFromFile(22,filename_before, "dfSelected_before");
  DrawECOUTFromFile(22,filename_after, "dfSelected_after");
  DrawECOUTFromFile(22,filename_before, "dfSelected_before");
  gApplication->Terminate(0);
}
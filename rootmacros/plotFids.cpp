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

void DrawDriftChamberRegionsFromFile(
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // Flatten DC hits per particle for each layer
  auto dfR1 =
      df.Define("xR1",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<float> &x,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid,
                   const ROOT::VecOps::RVec<short> &charge) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    //if (det[i] == 6 && layer[i] == 6 && p[pindex[i]] > 0.02 && pid[pindex[i]] == 11 && charge[pindex[i]]==-1)
                    if (det[i] == 6 && layer[i] == 6 && pid[pindex[i]] == 11)
                      out.push_back(x[i]);
                  }
                  
                  return out;
                },
                {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Particle_p", "REC_Traj_pindex","REC_Particle_pid","REC_Particle_charge"})
          .Define("yR1",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<float> &y,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid,
                     const ROOT::VecOps::RVec<short> &charge) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      //if (det[i] == 6 && layer[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11 && charge[pindex[i]]==-1)
                      if (det[i] == 6 && layer[i] == 6 && pid[pindex[i]] == 11)  
                        out.push_back(y[i]);
                    return out;
                  },
                  {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_y", "REC_Particle_p", "REC_Traj_pindex","REC_Particle_pid","REC_Particle_charge"});

  auto dfR2 =
      df.Define("xR2",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<float> &x,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 6 && layer[i] == 18 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == 11)
                      out.push_back(x[i]);
                  }
                  return out;
                },
                {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"})
          .Define("yR2",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<float> &y,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i) {
                      if (det[i] == 6 && layer[i] == 18 && p[pindex[i]] > 0.02 && pid[pindex[i]] == 11)
                        out.push_back(y[i]);
                    }
                    return out;
                  },
                  {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_y", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"});

  auto dfR3 =
      df.Define("xR3",
                [](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<float> &x,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 6 && layer[i] == 36 && p[pindex[i]] > 0.02 && pid[pindex[i]] == 11)
                      out.push_back(x[i]);
                  }
                  return out;
                },
                {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"})
          .Define("yR3",
                  [](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<float> &y,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i) {
                      if (det[i] == 6 && layer[i] == 36 && p[pindex[i]] > 0.02 && pid[pindex[i]] == 11)
                        out.push_back(y[i]);
                    }
                    return out;
                  },
                  {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_y", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"});

  // Draw histograms
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);

  TCanvas *c = new TCanvas("c", "DC Regions", 4000, 1000);
  c->Divide(3, 1);

  auto hR1 = dfR1.Histo2D({"hR1", "DC Region R1 (Layer 6);x [cm];y [cm]", 500,
                           -350, 350, 500, -400, 400},
                          "xR1", "yR1");
  auto hR2 = dfR2.Histo2D({"hR2", "DC Region R2 (Layer 18);x [cm];y [cm]", 500,
                           -350, 350, 500, -400, 400},
                          "xR2", "yR2");
  auto hR3 = dfR3.Histo2D({"hR3", "DC Region R3 (Layer 36);x [cm];y [cm]", 500,
                           -350, 350, 500, -400, 400},
                          "xR3", "yR3");

  c->cd(1);
  gPad->SetLogz();
  hR1->Draw("colz");
  DrawSectorLines();
  c->cd(2);
  gPad->SetLogz();
  hR2->Draw("colz");
  DrawSectorLines();
  c->cd(3);
  gPad->SetLogz();
  hR3->Draw("colz");
  DrawSectorLines();

  std::string outname = "DriftChamberRegions_" + treename + ".png";
  c->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}

void plotFids() {
  std::string
      filename_after =
          "../build/"
          "dfSelected_after_fiducialCuts.root"; // Update with your file path
  std::string
      filename_before =
          "../build/"
          "dfSelected_before_fiducialCuts.root"; // Update with your file path
  DrawDriftChamberRegionsFromFile(filename_after, "dfSelected_after");
  DrawDriftChamberRegionsFromFile(filename_before, "dfSelected_before");
  gApplication->Terminate(0);
}
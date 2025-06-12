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
    const int &selectedPid,
    const std::string &filename,
    const std::string &treename = "dfSelected_after") {
  ROOT::RDataFrame df(treename, filename);

  using namespace ROOT::VecOps;

  // Flatten DC hits per particle for each layer
  auto dfR1 =
      df.Define("xR1",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<float> &x,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid,
                   const ROOT::VecOps::RVec<short> &charge,
                   const ROOT::VecOps::RVec<int> &passFid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    //if (det[i] == 6 && layer[i] == 6 && p[pindex[i]] > 0.02 && pid[pindex[i]] == selectedPid && charge[pindex[i]]==-1)
                    if (det[i] == 6 && layer[i] == 6 && pid[pindex[i]] == selectedPid && passFid[pindex[i]])
                      out.push_back(x[i]);
                  }
                  
                  return out;
                },
                {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Particle_p", "REC_Traj_pindex","REC_Particle_pid","REC_Particle_charge", "REC_Track_pass_nofid"})
          .Define("yR1",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<float> &y,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid,
                     const ROOT::VecOps::RVec<short> &charge,
                    const ROOT::VecOps::RVec<int> &passFid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i)
                      //if (det[i] == 6 && layer[i] == 6 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid && charge[pindex[i]]==-1)
                      if (det[i] == 6 && layer[i] == 6 && pid[pindex[i]] == selectedPid && passFid[pindex[i]])  
                        out.push_back(y[i]);
                    return out;
                  },
                  {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_y", "REC_Particle_p", "REC_Traj_pindex","REC_Particle_pid","REC_Particle_charge", "REC_Track_pass_nofid"});

  auto dfR2 =
      df.Define("xR2",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<float> &x,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 6 && layer[i] == 18 && p[pindex[i]] >= 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(x[i]);
                  }
                  return out;
                },
                {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"})
          .Define("yR2",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<float> &y,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i) {
                      if (det[i] == 6 && layer[i] == 18 && p[pindex[i]] > 0.02 && pid[pindex[i]] == selectedPid)
                        out.push_back(y[i]);
                    }
                    return out;
                  },
                  {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_y", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"});

  auto dfR3 =
      df.Define("xR3",
                [selectedPid](const ROOT::VecOps::RVec<short> &det,
                   const ROOT::VecOps::RVec<short> &layer,
                   const ROOT::VecOps::RVec<float> &x,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<int16_t> &pindex,
                   const ROOT::VecOps::RVec<int> &pid) {
                  ROOT::VecOps::RVec<float> out;
                  for (size_t i = 0; i < layer.size(); ++i) {
                    if (det[i] == 6 && layer[i] == 36 && p[pindex[i]] > 0.02 && pid[pindex[i]] == selectedPid)
                      out.push_back(x[i]);
                  }
                  return out;
                },
                {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_x", "REC_Particle_p", "REC_Traj_pindex", "REC_Particle_pid"})
          .Define("yR3",
                  [selectedPid](const ROOT::VecOps::RVec<short> &det,
                     const ROOT::VecOps::RVec<short> &layer,
                     const ROOT::VecOps::RVec<float> &y,
                     const ROOT::VecOps::RVec<float> &p,
                     const ROOT::VecOps::RVec<int16_t> &pindex,
                     const ROOT::VecOps::RVec<int> &pid) {
                    ROOT::VecOps::RVec<float> out;
                    for (size_t i = 0; i < layer.size(); ++i) {
                      if (det[i] == 6 && layer[i] == 36 && p[pindex[i]] > 0.02 && pid[pindex[i]] == selectedPid)
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
  
  auto hR1 = dfR1.Histo2D({"hR1", ";x [cm];y [cm]", 500,
                           -350, 350, 500, -400, 400},
                          "xR1", "yR1");
  auto hR2 = dfR2.Histo2D({"hR2", ";x [cm];y [cm]", 500,
                           -350, 350, 500, -400, 400},
                          "xR2", "yR2");
  auto hR3 = dfR3.Histo2D({"hR3", ";x [cm];y [cm]", 500,
                           -350, 350, 500, -400, 400},
                          "xR3", "yR3");


  c->cd(1);
  gPad->SetLogz();
  hR1->Draw("colz");
  if (selectedPid == 11) {
    hR1->SetTitle("electron Drift Chamber R1");
  } else if (selectedPid == 2212) {
    hR1->SetTitle("proton Drift Chamber R1");
  } else {
    hR1->SetTitle("Drift Chamber R1");
  }
  DrawSectorLines();
  c->cd(2);
  gPad->SetLogz();
  hR2->Draw("colz");
  if (selectedPid == 11) {
    hR2->SetTitle("electron Drift Chamber R2");
  } else if (selectedPid == 2212) {
    hR2->SetTitle("proton Drift Chamber R2");
  } else {
    hR2->SetTitle("Drift Chamber R2");
  }
  DrawSectorLines();
  c->cd(3);
  gPad->SetLogz();
  hR3->Draw("colz");
  if (selectedPid == 11) {
    hR3->SetTitle("electron Drift Chamber R3");
  } else if (selectedPid == 2212) {
    hR3->SetTitle("proton Drift Chamber R3");
  } else {
    hR3->SetTitle("Drift Chamber R3");
  }
  DrawSectorLines();

  std::string outname = "DriftChamberRegions_" + treename + ".png";

  if (selectedPid == 11) {
    outname = "DriftChamberRegions_electron_" + treename + ".png";
  }
  else if (selectedPid == 2212) {
    outname = "DriftChamberRegions_proton_" + treename + ".png";
  }
  
  c->SaveAs(outname.c_str());
  std::cout << "Saved plot to: " << outname << std::endl;
}

void plotFids() {
  std::string input_path_from_analysisRun = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DC_fiducialcuts/";
  //std::string input_path_from_analysisRun = "./../build";
  std::string filename_after = Form("%s/dfSelected_after_fiducialCuts.root", input_path_from_analysisRun.c_str());
  std::string filename_before = Form("%s/dfSelected_before_fiducialCuts.root", input_path_from_analysisRun.c_str());
  DrawDriftChamberRegionsFromFile(11,filename_after, "dfSelected_after");
  DrawDriftChamberRegionsFromFile(11,filename_before, "dfSelected_before");
  DrawDriftChamberRegionsFromFile(2212,filename_after, "dfSelected_after");
  DrawDriftChamberRegionsFromFile(2212,filename_before, "dfSelected_before");
  gApplication->Terminate(0);
}
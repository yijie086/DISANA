#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

using namespace ROOT;
using namespace ROOT::VecOps;

// ------------- helpers -------------
static inline std::string MakeSafeName(std::string s) {
  for (auto &ch : s) if (!std::isalnum((unsigned char)ch)) ch = '_';
  return s;
}
static inline void MakeDir(const std::string& dir) {
  gSystem->Exec(std::string("mkdir -p \"" + dir + "\"").c_str());
}

// ------------------------
// Style helper functions
void StyleCanvas(TCanvas* c) {
  gStyle->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");

  c->SetTicks(1, 1);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.03);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.12);

  // white backgrounds
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  c->SetFillColor(0);
  if (c->GetPad(0)) c->GetPad(0)->SetFillColor(0);
}
void ApplyHistStyle(TH1* hist, double fit_range_min = 0.9, double fit_range_max = 1.25, TString xtitle = "M(K^{+}K^{-}) [GeV]", TString ytitle = "Counts") {
  hist->SetTitle("");
  hist->SetLineColor(kBlue + 1);
  hist->SetLineWidth(2);
  hist->SetFillColorAlpha(kBlue - 9, 0.3);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetMarkerColor(kBlue + 2);

  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.0);

  hist->GetYaxis()->SetTitle(ytitle);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetRangeUser(fit_range_min, fit_range_max);
}
void ApplyHistStyleTH2(TH2* hist, double fit_range_min_x = 0.9, double fit_range_max_x = 1.25, TString xtitle = "M(K^{+}K^{-}) [GeV]", TString ytitle = "Counts") {
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetYaxis()->SetTitle(ytitle);
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetRangeUser(fit_range_min_x, fit_range_max_x);
}

void Draw2DHistogram(TH2D* hist, const std::string& title, const std::string& outname, TString xtitle, TString ytitle, double xmin, double xmax, double ymin, double ymax) {
  TCanvas* c = new TCanvas(("c_" + outname).c_str(), title.c_str(), 1200, 1000);
  StyleCanvas(c);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.13);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.13);

  ApplyHistStyleTH2(hist, xmin, xmax, xtitle, ytitle);
  hist->GetXaxis()->SetRangeUser(xmin, xmax);
  hist->GetYaxis()->SetRangeUser(ymin, ymax);
  hist->Draw("COLZ");
  c->SaveAs((outname + ".png").c_str());
  c->SaveAs((outname + ".pdf").c_str());
  delete c;
}

// ------------------------
// Histogram Fitting & Drawing
// ------------------------
void FitAndDrawHistogram(TH1D* hist, const std::string& title, const std::string& outname, TString xtitle = "M(K^{+}K^{-}) [GeV]", TString ytitle = "Counts",
                         double fit_range_min = 0.9, double fit_range_max = 1.25, bool fit = true, const std::string& tag = "") {
  const std::string safe = MakeSafeName(tag);
  TCanvas* c = new TCanvas(("c_" + outname).c_str(), title.c_str(), 1200, 1000);
  StyleCanvas(c);
  ApplyHistStyle(hist, fit_range_min - .02, fit_range_max, xtitle, ytitle);
  hist->Draw("PE");
  // Legend
  TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist, "K^{+}K^{-} Inv. Mass", "lep");

  if (fit) {
    // unique TF1 names per dataset
    TF1* fitTotal = new TF1(("fitTotal_" + safe).c_str(),
      "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874)) + [3]*TMath::Gaus(x,[4],[5])", fit_range_min, fit_range_max);

    fitTotal->SetParameters(10, 0.9, 2, 100, 1.0, 0.04);
    fitTotal->SetParLimits(4, .990, 1.04);
    fitTotal->SetParLimits(5, 0.001, 0.03);
    fitTotal->SetLineColor(kRed + 1);
    fitTotal->SetLineWidth(3);
    hist->Fit(fitTotal, "R0");
    fitTotal->Draw("SAME C");

    double A = fitTotal->GetParameter(0);
    double alpha = fitTotal->GetParameter(1);
    double lambda = fitTotal->GetParameter(2);
    double N = fitTotal->GetParameter(3);
    double mu = fitTotal->GetParameter(4);
    double sigma = fitTotal->GetParameter(5);
    double chi2 = fitTotal->GetChisquare();
    double ndf = fitTotal->GetNDF();

    TF1* fitSignal = new TF1(("fitSignal_" + safe).c_str(), "[0]*TMath::Gaus(x,[1],[2])", fit_range_min, fit_range_max);
    fitSignal->SetParameters(N, mu, sigma);
    fitSignal->SetLineColor(kOrange + 1);
    fitSignal->SetLineStyle(2);
    fitSignal->SetLineWidth(2);
    fitSignal->SetFillColorAlpha(kOrange - 3, 0.3);
    fitSignal->SetFillStyle(1001);
    fitSignal->Draw("SAME FC");

    TF1* fitBkg = new TF1(("fitBkg_" + safe).c_str(), "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))", fit_range_min, fit_range_max);
    fitBkg->SetParameters(A, alpha, lambda);
    fitBkg->SetLineColor(kGreen + 2);
    fitBkg->SetLineStyle(3);
    fitBkg->SetLineWidth(2);
    fitBkg->SetFillColorAlpha(kGreen - 7, 0.3);
    fitBkg->SetFillStyle(1001);
    fitBkg->Draw("SAME FC");

    TLine* line1 = new TLine(mu - 3 * sigma, 0, mu - 3 * sigma, hist->GetMaximum() * 0.5);
    TLine* line2 = new TLine(mu + 3 * sigma, 0, mu + 3 * sigma, hist->GetMaximum() * 0.5);
    line1->SetLineColor(kMagenta + 2);
    line2->SetLineColor(kMagenta + 2);
    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line1->SetLineWidth(2);
    line2->SetLineWidth(2);
    line1->Draw("SAME");
    line2->Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.92, Form("%s", tag.empty() ? "" : tag.c_str()));
    latex.DrawLatex(0.15, 0.87, Form("#mu = %.3f GeV", mu));
    latex.DrawLatex(0.15, 0.82, Form("#sigma = %.4f GeV", sigma));
    latex.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf = %.2f", chi2 / ndf));
    latex.DrawLatex(0.15, 0.72, Form("Fit range: %.2f-%.2f GeV", fitTotal->GetXmin(), fitTotal->GetXmax()));

    leg->AddEntry(fitTotal, "Total Fit: Exp + Gauss", "l");
    leg->AddEntry(fitSignal, "Signal (Gauss)", "f");
    leg->AddEntry(fitBkg, "Background (Exp)", "f");
    leg->AddEntry(line1, "#pm 3#sigma", "l");
  }

  leg->Draw();
  c->SaveAs((outname + ".png").c_str());
  c->SaveAs((outname + ".pdf").c_str());
  delete c;
}

// ------------------------
// Main Analysis Routine (single dataset)
// ------------------------
void DrawPhiMass(const std::string& filename, const std::string& treename,
                 const std::string& outputDir = "PhiMassPlots", const std::string& tag = "") {
  const std::string safe = MakeSafeName(tag);
  TStopwatch timer;
  timer.Start();
  ROOT::EnableImplicitMT();

  RDataFrame df(treename, filename);
  auto df_filtered =
      df.Filter("REC_MotherMass.size() > 0")
        .Define("REC_DaughterParticle_pass_int",
                [](const std::vector<bool>& passVec) { return RVec<int>(passVec.begin(), passVec.end()); },
                {"REC_DaughterParticle_pass"});

  auto df_beforeCut = df_filtered.Define("MotherMass_all",
      [](const std::vector<float>& masses) { return RVec<float>(masses.begin(), masses.end()); }, {"REC_MotherMass"});

  auto df_afterCut = df_filtered.Define("MotherMass_passed",
      [](const RVec<float>& masses, const RVec<int>& pass) {
        RVec<float> out;
        for (size_t i = 0; i < masses.size(); ++i)
          if (pass[i] == 0 && masses[i] > 0) out.push_back(masses[i]);
        return out;
      },
      {"REC_MotherMass", "REC_DaughterParticle_pass_int"});

  double fit_range_min = 0.9874;
  double fit_range_max = 1.120;

  auto df_selectedKaons = df_filtered
    .Define("KaonIndices",
      [](const RVec<int>& pid, const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz) {
        RVec<size_t> indices;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (std::abs(pid[i]) == 321) {
            float p = std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            if (p < 3.5) indices.push_back(i);
          }
        }
        return indices;
      }, {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz"})
    .Define("KK_mass",
      [](const RVec<int>& pid, const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, const RVec<size_t>& indices) {
        RVec<float> masses;
        const float mK = 0.493677;
        for (size_t i = 0; i < indices.size(); ++i) {
          for (size_t j = i + 1; j < indices.size(); ++j) {
            int pid1 = pid[indices[i]];
            int pid2 = pid[indices[j]];
            if (pid1 * pid2 == -321 * 321) {
              float px1 = px[indices[i]], py1 = py[indices[i]], pz1 = pz[indices[i]];
              float px2 = px[indices[j]], py2 = py[indices[j]], pz2 = pz[indices[j]];
              float E1 = std::sqrt(px1*px1 + py1*py1 + pz1*pz1 + mK*mK);
              float E2 = std::sqrt(px2*px2 + py2*py2 + pz2*pz2 + mK*mK);
              float px_tot = px1 + px2;
              float py_tot = py1 + py2;
              float pz_tot = pz1 + pz2;
              float E_tot = E1 + E2;
              float mass2 = E_tot*E_tot - (px_tot*px_tot + py_tot*py_tot + pz_tot*pz_tot);
              if (mass2 > 0) masses.push_back(std::sqrt(mass2));
            }
          }
        }
        return masses;
      }, {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz","KaonIndices"});

  auto df_withPKmPKp = df_selectedKaons
    .Define("PKm_mass_vs_KK",
      [](const RVec<int>& pid, const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, const RVec<size_t>& indices, const RVec<float>& KK_mass_vec) {
        RVec<std::pair<float,float>> out;
        const float mK = 0.493677;
        const float mP = 0.938272;
        for (size_t i = 0; i < indices.size(); ++i) {
          int kaon_pid = pid[indices[i]];
          if (kaon_pid != -321) continue;
          float k_px = px[indices[i]], k_py = py[indices[i]], k_pz = pz[indices[i]];
          float E_K = std::sqrt(k_px*k_px + k_py*k_py + k_pz*k_pz + mK*mK);
          for (size_t j = 0; j < pid.size(); ++j) {
            if (pid[j] != 2212) continue;
            float p_px = px[j], p_py = py[j], p_pz = pz[j];
            float E_P = std::sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz + mP*mP);
            float E_tot = E_K + E_P;
            float px_tot = k_px + p_px;
            float py_tot = k_py + p_py;
            float pz_tot = k_pz + p_pz;
            float mass2 = E_tot*E_tot - (px_tot*px_tot + py_tot*py_tot + pz_tot*pz_tot);
            if (mass2 > 0 && KK_mass_vec.size() > 0) out.emplace_back(std::sqrt(mass2), KK_mass_vec[0]);
          }
        }
        return out;
      }, {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz","KaonIndices","KK_mass"})
    .Define("PKp_mass_vs_KK",
      [](const RVec<int>& pid, const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, const RVec<size_t>& indices, const RVec<float>& KK_mass_vec) {
        RVec<std::pair<float,float>> out;
        const float mK = 0.493677;
        const float mP = 0.938272;
        for (size_t i = 0; i < indices.size(); ++i) {
          int kaon_pid = pid[indices[i]];
          if (kaon_pid != 321) continue;
          float k_px = px[indices[i]], k_py = py[indices[i]], k_pz = pz[indices[i]];
          float E_K = std::sqrt(k_px*k_px + k_py*k_py + k_pz*k_pz + mK*mK);
          for (size_t j = 0; j < pid.size(); ++j) {
            if (pid[j] != 2212) continue;
            float p_px = px[j], p_py = py[j], p_pz = pz[j];
            float E_P = std::sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz + mP*mP);
            float E_tot = E_K + E_P;
            float px_tot = k_px + p_px;
            float py_tot = k_py + p_py;
            float pz_tot = k_pz + p_pz;
            float mass2 = E_tot*E_tot - (px_tot*px_tot + py_tot*py_tot + pz_tot*pz_tot);
            if (mass2 > 0 && KK_mass_vec.size() > 0) out.emplace_back(std::sqrt(mass2), KK_mass_vec[0]);
          }
        }
        return out;
      }, {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz","KaonIndices","KK_mass"});

  auto df_withPK_only = df_withPKmPKp
    .Define("PKm_mass", [](const RVec<std::pair<float,float>>& pairs){ RVec<float> r; for (auto &p: pairs) r.push_back(p.first); return r; }, {"PKm_mass_vs_KK"})
    .Define("PKp_mass", [](const RVec<std::pair<float,float>>& pairs){ RVec<float> r; for (auto &p: pairs) r.push_back(p.first); return r; }, {"PKp_mass_vs_KK"});

  auto df_split2D = df_withPK_only
    .Define("PKm_mass_X", [](const RVec<std::pair<float,float>>& v){ RVec<float> x; for (auto &p: v) x.push_back(p.first); return x; }, {"PKm_mass_vs_KK"})
    .Define("PKm_mass_Y", [](const RVec<std::pair<float,float>>& v){ RVec<float> y; for (auto &p: v) y.push_back(p.second); return y; }, {"PKm_mass_vs_KK"})
    .Define("PKp_mass_X", [](const RVec<std::pair<float,float>>& v){ RVec<float> x; for (auto &p: v) x.push_back(p.first); return x; }, {"PKp_mass_vs_KK"})
    .Define("PKp_mass_Y", [](const RVec<std::pair<float,float>>& v){ RVec<float> y; for (auto &p: v) y.push_back(p.second); return y; }, {"PKp_mass_vs_KK"});

  // unique histogram names per dataset
  auto h_KKmass = df_selectedKaons.Histo1D(
      {Form("hKKMass_%s", safe.c_str()), "After Mass Cut;M(K^{+}K^{-}) [GeV];Counts", 200, 0.80, 1.6}, "KK_mass");

  auto h2_PKm_vs_KK = df_split2D.Histo2D(
      {Form("h2PKmKK_%s", safe.c_str()), "M(pK^{-}) vs M(K^{+}K^{-});M(K^{+}K^{-}) [GeV];M(pK^{-}) [GeV]", 200, .980, 2.5, 100, 1.4, 4.5},
      "PKm_mass_Y", "PKm_mass_X");

  auto h2_PKp_vs_KK = df_split2D.Histo2D(
      {Form("h2PKpKK_%s", safe.c_str()), "M(pK^{+}) vs M(K^{+}K^{-});M(K^{+}K^{-}) [GeV];M(pK^{+}) [GeV]", 200, 0.980, 2.5, 100, 1.4, 4.5},
      "PKp_mass_Y", "PKp_mass_X");

  auto h_PKm = df_withPK_only.Histo1D(
      {Form("hPKm_%s", safe.c_str()), "Invariant Mass of pK^{-};M(pK^{-}) [GeV];Counts", 200, 1.2, 2.6}, "PKm_mass");

  auto h_PKp = df_withPK_only.Histo1D(
      {Form("hPKp_%s", safe.c_str()), "Invariant Mass of pK^{+};M(pK^{+}) [GeV];Counts", 200, 1.2, 2.6}, "PKp_mass");

  MakeDir(outputDir);

  FitAndDrawHistogram(h_KKmass.GetPtr(), "After InvMass Cut", outputDir + "/KpKm_mass_after",
                      "M(K^{+}K^{-}) [GeV]", "Counts", fit_range_min, fit_range_max, true, tag);

  FitAndDrawHistogram(h_PKm.GetPtr(), "After InvMass Cut", outputDir + "/pKm_mass_after",
                      "M(pK^{-}) [GeV]", "Counts", 1.0, 4.0, false, tag);

  FitAndDrawHistogram(h_PKp.GetPtr(), "After InvMass Cut", outputDir + "/pKp_mass_after",
                      "M(pK^{+}) [GeV]", "Counts", 1.0, 4.0, false, tag);

  Draw2DHistogram(h2_PKm_vs_KK.GetPtr(), "M(K^{+}K^{-}) vs M(pK^{-})",
                  outputDir + "/pKm_vs_KK", "M(K^{+}K^{-}) [GeV]","M(pK^{-}) [GeV]", 0.98, 2.2, 1.3, 2.8);

  Draw2DHistogram(h2_PKp_vs_KK.GetPtr(), "M(K^{+}K^{-}) vs M(pK^{+})",
                  outputDir + "/pKp_vs_KK", "M(K^{+}K^{-}) [GeV]","M(pK^{+}) [GeV]", 0.98, 2.2, 1.3, 2.8);

  std::cout << "Plots saved in: " << outputDir << std::endl;
  timer.Stop();
  std::cout << "Time for DrawPhiMass: ";
  timer.Print();
}

// ------------------------
// Batch driver over all periods / polarities
// ------------------------
void analysisPhiMass() {
  // ---------- your provided paths ----------
  std::string input_path_from_analysisRun_SP18inb_data  = "./../data_processed/spring2018/inb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_SP18outb_data = "./../data_processed/spring2018/outb/DVKpKm_wagon/";

  std::string input_path_from_analysisRun_Fall18inb_data  = "./../data_processed/fall2018/inb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_Fall18outb_data = "./../data_processed/fall2018/outb/DVKpKm_wagon/";

  std::string input_path_from_analysisRun_SP19inb_data  = "./../data_processed/spring2019/inb/DVKpKm_wagon/";

  std::string filename_afterFid_SP18inb_data  = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18inb_data.c_str());
  std::string filename_afterFid_SP18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18outb_data.c_str());

  std::string filename_afterFid_Fall18inb_data  = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18inb_data.c_str());
  std::string filename_afterFid_Fall18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18outb_data.c_str());

  std::string filename_afterFid_SP19inb_data  = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP19inb_data.c_str());

  // ---------- run config ----------
  struct Job {
    std::string filename;
    std::string treename;
    std::string outdir;
    std::string tag;     // shown on canvas & used to uniquify object names
  };

  const std::string baseOut = "PhiMassPlots";
  std::vector<Job> jobs = {
    { filename_afterFid_SP18inb_data,  "dfSelected_afterFid", baseOut + "/spring2018/inb",  "Spring 2018 INB"  },
    { filename_afterFid_SP18outb_data, "dfSelected_afterFid", baseOut + "/spring2018/outb", "Spring 2018 OUTB" },
    { filename_afterFid_Fall18inb_data,"dfSelected_afterFid", baseOut + "/fall2018/inb",    "Fall 2018 INB"    },
    { filename_afterFid_Fall18outb_data,"dfSelected_afterFid",baseOut + "/fall2018/outb",   "Fall 2018 OUTB"   },
    { filename_afterFid_SP19inb_data,  "dfSelected_afterFid", baseOut + "/spring2019/inb",  "Spring 2019 INB"  }
  };

  // Create top-level dirs first (and then per job)
  MakeDir(baseOut);
  MakeDir(baseOut + "/spring2018/inb");
  MakeDir(baseOut + "/spring2018/outb");
  MakeDir(baseOut + "/fall2018/inb");
  MakeDir(baseOut + "/fall2018/outb");
  MakeDir(baseOut + "/spring2019/inb");

  // ---------- execute ----------
  for (const auto& j : jobs) {
    std::cout << "\n>>> Processing " << j.tag << " from file: " << j.filename << std::endl;
    DrawPhiMass(j.filename, j.treename, j.outdir, j.tag);
  }

  gApplication->Terminate(0);
}

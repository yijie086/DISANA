#pragma once

// ROOT
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

// STL
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

// Your helpers
#include "DISANAMathFitUtils.h"   // FitWindowFromGaussianCore
#include "DISANAMath.h"           // columns like Mx2_epKp, Mx2_epKm (+ others)

using namespace ROOT;
using namespace ROOT::VecOps;

namespace DISANA {
namespace PhiMass {

// -------------------------------------------------------------------------------------
// constants (kept for completeness; we no longer re-calc MM from 4-vectors here)
// -------------------------------------------------------------------------------------
static constexpr double kMe = 0.000511;
static constexpr double kMp = 0.938272;
static constexpr double kMK = 0.493677;

// -------------------------------------------------------------------------------------
// Local fit model for missing-mass: gaus(0) + pol3(3)  (unchanged style)
// -------------------------------------------------------------------------------------
inline TF1* DISANA_FitGausPlusPoly3(TH1* h,
                                    double mu_guess = 0.4937,
                                    double sigma_guess = 0.030,
                                    double xmin = -1,
                                    double xmax = -1) {
  if (!h) throw std::runtime_error("DISANA_FitGausPlusPoly3: null histogram");
  if (xmin < 0 || xmax < 0) { xmin = h->GetXaxis()->GetXmin(); xmax = h->GetXaxis()->GetXmax(); }
  std::string fname = Form("fKminusMx_gP3_%p", (void*)h);
  TF1* f = new TF1(fname.c_str(), "gaus(0)+pol3(3)", xmin, xmax);
  f->SetParNames("A", "mu", "sigma", "p0", "p1", "p2", "p3");
  f->SetParameters(h->GetMaximum(), mu_guess, sigma_guess, 0, 0, 0, 0);
  f->SetParLimits(2, 0.003, 0.120);  // 3–120 MeV width
  h->Fit(f, "RQ0");                  // quiet, bounded, don’t draw
  return f;
}

// -------------------------------------------------------------------------------------
// Styling & Dir helpers (unchanged)
// -------------------------------------------------------------------------------------
inline std::string MakeSafeName(std::string s) {
  for (auto& ch : s) if (!std::isalnum((unsigned char)ch)) ch = '_';
  return s;
}
inline void MakeDir(const std::string& dir) { gSystem->Exec(std::string("mkdir -p \"" + dir + "\"").c_str()); }

inline void StyleCanvas(TCanvas* c) {
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
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  c->SetFillColor(0);
  if (c->GetPad(0)) c->GetPad(0)->SetFillColor(0);
}

inline void ApplyHistStyle(TH1* hist,
                           double fit_range_min = 0.9,
                           double fit_range_max = 1.25,
                           TString xtitle = "M(K^{+}K^{-}) [GeV]",
                           TString ytitle = "Counts") {
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

inline void ApplyHistStyleTH2(TH2* hist,
                              double fit_range_min_x = 0.9,
                              double fit_range_max_x = 1.25,
                              TString xtitle = "M(K^{+}K^{-}) [GeV]",
                              TString ytitle = "Counts") {
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

inline void Draw2DHistogram(TH2D* hist,
                            const std::string& title,
                            const std::string& outname,
                            TString xtitle, TString ytitle,
                            double xmin, double xmax, double ymin, double ymax) {
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

// -------------------------------------------------------------------------------------
// φ-mass fitter/drawer (unchanged style)
// -------------------------------------------------------------------------------------
inline void FitAndDrawHistogram(TH1D* hist,
                                const std::string& title,
                                const std::string& outname,
                                TString xtitle = "M(K^{+}K^{-}) [GeV]",
                                TString ytitle = "Counts",
                                double fit_range_min = 0.9,
                                double fit_range_max = 1.25,
                                bool fit = true,
                                const std::string& tag = "",
                                double annotate_sigma_for_count = -1.0) {
  const std::string safe = MakeSafeName(tag);
  TCanvas* c = new TCanvas(("c_" + outname).c_str(), title.c_str(), 1200, 1000);
  StyleCanvas(c);
  ApplyHistStyle(hist, fit_range_min - .02, fit_range_max, xtitle, ytitle);
  hist->Draw("PE");

  TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist, "K^{+}K^{-} Inv. Mass", "lep");

  double mu = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();

  if (fit) {
    TF1* fitTotal = new TF1(("fitTotal_" + safe).c_str(),
                            "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874)) + [3]*TMath::Gaus(x,[4],[5])",
                            fit_range_min, fit_range_max);

    fitTotal->SetParameters(10, 0.9, 2, 100, 1.02, 0.010);
    fitTotal->SetParLimits(4, .990, 1.050);
    fitTotal->SetParLimits(5, 0.001, 0.03);
    fitTotal->SetLineColor(kRed + 1);
    fitTotal->SetLineWidth(3);
    hist->Fit(fitTotal, "R0Q");
    fitTotal->Draw("SAME C");

    double A = fitTotal->GetParameter(0);
    double alpha = fitTotal->GetParameter(1);
    double lambda = fitTotal->GetParameter(2);
    double N = fitTotal->GetParameter(3);
    mu = fitTotal->GetParameter(4);
    sigma = fitTotal->GetParameter(5);
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

    TF1* fitBkg = new TF1(("fitBkg_" + safe).c_str(),
                           "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))",
                           fit_range_min, fit_range_max);
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
    latex.DrawLatex(0.15, 0.87, Form("#mu = %.4f GeV", mu));
    latex.DrawLatex(0.15, 0.82, Form("#sigma = %.4f GeV", sigma));
    latex.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf = %.2f", chi2 / std::max(1.0, ndf)));
    latex.DrawLatex(0.15, 0.72, Form("Fit range: %.3f-%.3f GeV", fit_range_min, fit_range_max));

    if (annotate_sigma_for_count > 0 && std::isfinite(mu) && std::isfinite(sigma) && sigma > 0) {
      const double lo = mu - annotate_sigma_for_count * sigma;
      const double hi = mu + annotate_sigma_for_count * sigma;
      int binLo = hist->GetXaxis()->FindBin(lo);
      int binHi = hist->GetXaxis()->FindBin(hi);
      long long nCand = 0;
      for (int b = binLo; b <= binHi; ++b) nCand += (long long)std::llround(hist->GetBinContent(b));
      latex.DrawLatex(0.15, 0.67, Form("N_{|M-#mu|<%.0f#sigma} = %lld", annotate_sigma_for_count, nCand));
    }

    leg->AddEntry(fitTotal, "Total Fit: Thr.Exp + Gauss", "l");
    leg->AddEntry(fitSignal, "Signal (Gauss)", "f");
    leg->AddEntry(fitBkg, "Background (Thr.Exp)", "f");
    leg->AddEntry(line1, "#pm 3#sigma", "l");
  }

  leg->Draw();
  c->SaveAs((outname + ".png").c_str());
  c->SaveAs((outname + ".pdf").c_str());
  delete c;
}

// -------------------------------------------------------------------------------------
// Nice drawers for missing mass with identical style (separate K− / K+ labels)
// -------------------------------------------------------------------------------------
struct KmMMWindow { double mu{0}, sigma{0}, low{0}, high{0}; };

inline KmMMWindow DrawFit_KaonMM(TH1D* hMx,
                                 const std::string& outBase,
                                 const std::string& tag,
                                 double xMin, double xMax, double nSigma,
                                 const char* axisTitle) {
  const std::string safe = MakeSafeName(tag);
  TCanvas* c = new TCanvas(Form("c_KaonMx_%s", safe.c_str()), "Kaon missing mass", 1200, 1000);
  StyleCanvas(c);
  ApplyHistStyle(hMx, xMin, xMax, axisTitle, "Counts");
  hMx->Draw("PE");

  TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hMx, "Kaon missing mass", "lep");

  KmMMWindow win;
  try {
    TF1* fTot = DISANA_FitGausPlusPoly3(hMx, 0.4937, 0.030, xMin, xMax);
    fTot->SetLineColor(kRed + 1);
    fTot->SetLineWidth(3);
    fTot->Draw("SAME C");

    const double A = fTot->GetParameter(0);
    win.mu      = fTot->GetParameter(1);
    win.sigma   = fTot->GetParameter(2);
    const double p0 = fTot->GetParameter(3);
    const double p1 = fTot->GetParameter(4);
    const double p2 = fTot->GetParameter(5);
    const double p3 = fTot->GetParameter(6);

    TF1* fSig = new TF1(Form("fSig_%s", safe.c_str()), "gaus(0)", xMin, xMax);
    fSig->SetParameters(A, win.mu, win.sigma);
    fSig->SetLineColor(kOrange + 1);
    fSig->SetLineStyle(2);
    fSig->SetLineWidth(2);
    fSig->SetFillColorAlpha(kOrange - 3, 0.3);
    fSig->SetFillStyle(1001);
    fSig->Draw("SAME FC");

    TF1* fBkg = new TF1(Form("fBkg_%s", safe.c_str()), "pol3(0)", xMin, xMax);
    fBkg->SetParameters(p0, p1, p2, p3);
    fBkg->SetLineColor(kGreen + 2);
    fBkg->SetLineStyle(3);
    fBkg->SetLineWidth(2);
    fBkg->SetFillColorAlpha(kGreen - 7, 0.3);
    fBkg->SetFillStyle(1001);
    fBkg->Draw("SAME FC");

    const auto lr = FitWindowFromGaussianCore(fTot, nSigma);
    win.low  = lr.first;
    win.high = lr.second;

    TLine L1(win.low, 0, win.low, hMx->GetMaximum() * 0.6);
    TLine L2(win.high, 0, win.high, hMx->GetMaximum() * 0.6);
    L1.SetLineColor(kMagenta + 2);
    L2.SetLineColor(kMagenta + 2);
    L1.SetLineStyle(2);
    L2.SetLineStyle(2);
    L1.SetLineWidth(2);
    L2.SetLineWidth(2);
    L1.Draw("SAME");
    L2.Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.92, Form("%s", tag.c_str()));
    latex.DrawLatex(0.15, 0.87, Form("#mu = %.4f GeV", win.mu));
    latex.DrawLatex(0.15, 0.82, Form("#sigma = %.4f GeV", win.sigma));
    latex.DrawLatex(0.15, 0.77, Form("Fit range: %.3f-%.3f GeV", xMin, xMax));

    leg->AddEntry(fTot, "Total Fit: Gauss + poly3", "l");
    leg->AddEntry(fSig, "Signal (Gauss)", "f");
    leg->AddEntry(fBkg, "Background (poly3)", "f");
    leg->AddEntry(&L1, Form("#pm %.0f#sigma", nSigma), "l");
  } catch (...) {
    win.mu = 0.4937;
    win.sigma = 0.030;
    win.low = 0.40;
    win.high = 0.58;
  }

  leg->Draw();
  gPad->RedrawAxis();
  gPad->Update();

  TCanvas* pad = (TCanvas*)gPad;
  pad->SaveAs((outBase + ".png").c_str());
  pad->SaveAs((outBase + ".pdf").c_str());
  delete pad;
  return win;
}

// Keep the legacy function name for K−-missing style compatibility
inline KmMMWindow DrawFit_KminusMM(TH1D* hMx, const std::string& outBase, const std::string& tag,
                                   double xMin, double xMax, double nSigma) {
  return DrawFit_KaonMM(hMx, outBase, tag, xMin, xMax, nSigma, "M_{X}(epK^{+}) [GeV]");
}

// Matching drawer for K⁺-missing (epK⁻)
inline KmMMWindow DrawFit_KplusMM(TH1D* hMx, const std::string& outBase, const std::string& tag,
                                  double xMin, double xMax, double nSigma) {
  return DrawFit_KaonMM(hMx, outBase, tag, xMin, xMax, nSigma, "M_{X}(epK^{-}) [GeV]");
}

// -------------------------------------------------------------------------------------
// API 1a: K− missing-mass → ±nσ window cut (uses DISANAMath::Mx2_epKp)
// -------------------------------------------------------------------------------------
inline std::tuple<ROOT::RDF::RNode, KmMMWindow>
WireKminusMMCut(ROOT::RDF::RNode df_in,
                float /*beam_energy_unused_now*/,
                const std::string& outputDir,
                const std::string& tag,
                int nBins = 220, double xMin = 0.30, double xMax = 0.70, double nSigma = 3.0) {
  const std::string safe = MakeSafeName(tag);

  // Use the column produced by InitKinematics/define_DISCAT
  auto df = df_in
              .Define("Mx_epKp_forCut", [](double mx2) { return (mx2 > 0.0) ? std::sqrt(mx2) : -999.0; }, {"Mx2_epKp"})
              .Filter("Mx_epKp_forCut>0", "valid K^{-} missing mass");

  auto hMx = df.Histo1D({Form("hKminusMx_%s", safe.c_str()),
                         "K^{-} missing mass;M_{X}(epK^{+}) [GeV];Counts",
                         nBins, xMin, xMax},
                        "Mx_epKp_forCut");

  MakeDir(outputDir);
  const std::string outBase = outputDir + "/KminusMM_" + safe;
  auto win = DrawFit_KminusMM(hMx.GetPtr(), outBase, tag, xMin, xMax, nSigma);

  auto df_cut = df.Filter([lo = win.low, hi = win.high](double x) { return x > lo && x < hi; },
                          {"Mx_epKp_forCut"},
                          Form("Cut: M_{X}(epK^{+}) in [%.3f, %.3f]", win.low, win.high));
  return {df_cut, win};
}

// -------------------------------------------------------------------------------------
// API 1b: K⁺ missing-mass → ±nσ window cut (uses DISANAMath::Mx2_epKm)
// -------------------------------------------------------------------------------------
inline std::tuple<ROOT::RDF::RNode, KmMMWindow>
WireKplusMMCut(ROOT::RDF::RNode df_in,
               float /*beam_energy_unused_now*/,
               const std::string& outputDir,
               const std::string& tag,
               int nBins = 220, double xMin = 0.30, double xMax = 0.70, double nSigma = 3.0) {
  const std::string safe = MakeSafeName(tag);

  auto df = df_in
              .Define("Mx_epKm_forCut", [](double mx2) { return (mx2 > 0.0) ? std::sqrt(mx2) : -999.0; }, {"Mx2_epKm"})
              .Filter("Mx_epKm_forCut>0", "valid K^{+} missing mass");

  auto hMx = df.Histo1D({Form("hKplusMx_%s", safe.c_str()),
                         "K^{+} missing mass;M_{X}(epK^{-}) [GeV];Counts",
                         nBins, xMin, xMax},
                        "Mx_epKm_forCut");

  MakeDir(outputDir);
  const std::string outBase = outputDir + "/KplusMM_" + safe;
  auto win = DrawFit_KplusMM(hMx.GetPtr(), outBase, tag, xMin, xMax, nSigma);

  auto df_cut = df.Filter([lo = win.low, hi = win.high](double x) { return x > lo && x < hi; },
                          {"Mx_epKm_forCut"},
                          Form("Cut: M_{X}(epK^{-}) in [%.3f, %.3f]", win.low, win.high));
  return {df_cut, win};
}

// -------------------------------------------------------------------------------------
// API 2: φ→K⁺K⁻ mass plotter from node
// Prefers InitKinematics’ "invMass_KpKm" column, falls back to a local Define.
// -------------------------------------------------------------------------------------
inline void DrawPhiMass_FromNode(ROOT::RDF::RNode df_in,
                                 const std::string& outputDir,
                                 const std::string& tag,
                                 double fit_range_min = 0.9874,
                                 double fit_range_max = 1.120) {
  const std::string safe = MakeSafeName(tag);

  // Try to use existing column "invMass_KpKm" if already defined by InitKinematics
  // If the column doesn't exist, define it from REC_* (graceful fallback).
  // The try/catch here is compile-time-unrelated; RDataFrame throws at run-time if column is missing.
  ROOT::RDF::RNode df = df_in;
  bool need_define = false;
  try {
    // cheap probe: create a trivial pass-through Define that depends on invMass_KpKm;
    // if column is missing, an exception will trigger upon book/Run.
    df = df.Define("invMass_KpKm_passthru", [](double v){ return v; }, {"invMass_KpKm"})
           .Filter("invMass_KpKm_passthru>0", "has K+K- pair");
  } catch (...) {
    need_define = true;
  }

  if (need_define) {
    df = df_in
           .Define("invMass_KpKm",
                   [](const RVec<int>& pid, const RVec<bool>& pass,
                      const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz) {
                     int iKp = -1, iKm = -1;
                     for (size_t i = 0; i < pid.size(); ++i) {
                       if (!pass[i]) continue;
                       if (pid[i] == 321  && iKp < 0)  iKp = (int)i;
                       if (pid[i] == -321 && iKm < 0)  iKm = (int)i;
                     }
                     if (iKp < 0 || iKm < 0) return -1.0f;
                     const float mK = 0.493677f;
                     auto E = [&](int i) { return std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] + mK*mK); };
                     float Ex = E(iKp) + E(iKm);
                     float pxs = px[iKp] + px[iKm], pys = py[iKp] + py[iKm], pzs = pz[iKp] + pz[iKm];
                     float m2 = Ex*Ex - (pxs*pxs + pys*pys + pzs*pzs);
                     return (m2 > 0) ? std::sqrt(m2) : -1.0f;
                   },
                   {"REC_Particle_pid", "REC_Particle_pass", "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz"})
           .Filter("invMass_KpKm>0", "has K+K- pair");
  }

  auto h = df.Histo1D({Form("hKK_fromNode_%s", safe.c_str()),
                       "K^{+}K^{-} invariant mass;M(K^{+}K^{-}) [GeV];Counts",
                       200, 0.80, 1.60},
                      need_define ? "invMass_KpKm" : "invMass_KpKm_passthru");

  MakeDir(outputDir);
  FitAndDrawHistogram(h.GetPtr(),
                      "After kaon missing-mass window",
                      outputDir + "/KpKm_mass_afterCut_" + safe,
                      "M(K^{+}K^{-}) [GeV]", "Counts",
                      fit_range_min, fit_range_max,
                      /*fit=*/true, tag,
                      /*annotate_sigma_for_count=*/3.0);
}

}  // namespace PhiMass
}  // namespace DISANA

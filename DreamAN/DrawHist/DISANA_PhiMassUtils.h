#pragma once

// ROOT
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TString.h>

// RDataFrame
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

// STL
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Your helper (for FitWindowFromGaussianCore)
#include "DISANAMathFitUtils.h"

namespace DISANA {
namespace PhiMass {

using namespace ROOT;
using namespace ROOT::VecOps;

// -------------------------------------------------------------------------------------
// constants
// -------------------------------------------------------------------------------------
static constexpr double kMe = 0.000511;
static constexpr double kMp = 0.938272;
static constexpr double kMK = 0.493677; // GeV (kaon mass)

// -------------------------------------------------------------------------------------
// Small styling helpers
// -------------------------------------------------------------------------------------
inline std::string MakeSafeName(std::string s) {
  for (auto &ch : s) if (!std::isalnum((unsigned char)ch)) ch = '_';
  return s;
}
inline void MakeDir(const std::string &dir) { gSystem->Exec(std::string("mkdir -p \"" + dir + "\"").c_str()); }

inline void StyleCanvas(TCanvas *c) {
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

inline void StyleH1(TH1 *h, TString xtitle, TString ytitle, double xmin, double xmax) {
  h->SetTitle("");
  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(2);
  h->SetFillColorAlpha(kBlue - 9, 0.30);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(kBlue + 2);

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleOffset(1.35);
  h->GetXaxis()->SetRangeUser(xmin, xmax);
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
// Model: total = gaus(0) + pol3(3)  (for missing-mass peak on smooth background)
// -------------------------------------------------------------------------------------
inline TF1 *FitGausPlusPoly3(TH1 *h, double mu_guess, double sigma_guess,
                             double xmin, double xmax, const std::string &tag) {
  if (!h) throw std::runtime_error("FitGausPlusPoly3: histogram is null");
  TF1 *f = new TF1(("fGMp3_"+tag).c_str(), "gaus(0)+pol3(3)", xmin, xmax);
  f->SetParNames("A","mu","sigma","p0","p1","p2","p3");
  f->SetParameters(h->GetMaximum(), mu_guess, sigma_guess, 0, 0, 0, 0);
  f->SetParLimits(2, 0.003, 0.120);
  h->Fit(f, "RQ0");
  return f;
}

// -------------------------------------------------------------------------------------
// Window result
// -------------------------------------------------------------------------------------
struct Window { double mu{0}, sigma{0}, low{0}, high{0}; };
struct KmMMWindow { double mu{0}, sigma{0}, low{0}, high{0}; };
struct PhiMassWindow { double mu{0}, sigma{0}, low{0}, high{0}; };

// -------------------------------------------------------------------------------------
// Unified missing mass drawer
// -------------------------------------------------------------------------------------
inline KmMMWindow DrawFit_KaonMM(TH1D* hMx,
                                 const std::string& outBase,
                                 const std::string& tag,
                                 double xMin, double xMax, double nSigma,
                                 const char* axisTitle,
                                 const char* legendEntry) {
  const std::string safe = MakeSafeName(tag);
  TCanvas* c = new TCanvas(Form("c_KaonMx_%s", safe.c_str()), "Kaon missing mass", 1200, 1000);
  StyleCanvas(c);
  StyleH1(hMx, axisTitle, "Counts", xMin, xMax);
  hMx->Draw("PE");

  TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hMx, legendEntry, "lep");

  KmMMWindow win;
  try {
    TF1* fTot = FitGausPlusPoly3(hMx, kMK, 0.030, xMin, xMax, safe);
    fTot->SetLineColor(kRed + 1);
    fTot->SetLineWidth(3);
    fTot->Draw("SAME C");

    const double A = fTot->GetParameter(0);
    win.mu = fTot->GetParameter(1);
    win.sigma = fTot->GetParameter(2);
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
    win.low = lr.first;
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
    win.mu = kMK;
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

// -------------------------------------------------------------------------------------
// GENERIC: build |Mx| from an Mx² column that already exists (from DISANAMath)
//          then fit, draw, and return a ±Nσ window + filtered node.
// - mx2_col:  "Mx2_epKp" (K− missing)  or  "Mx2_epKm" (K+ missing)
// - kaon_label_for_axes: "K^{-}" or "K^{+}" just for pretty axis titles
// -------------------------------------------------------------------------------------
inline std::tuple<ROOT::RDF::RNode, KmMMWindow>
WireKaonMMCut_FromMx2(ROOT::RDF::RNode df_in,
                      const std::string &mx2_col,
                      const std::string &outputDir,
                      const std::string &tag,
                      const std::string &kaon_label_for_axes,
                      int nBins = 220, double xMin = 0.30, double xMax = 0.70,
                      double nSigma = 3.0,
                      double mu_seed = kMK, double sigma_seed = 0.030) {

  const std::string safe = MakeSafeName(tag);
  const std::string mx_col = "Mx_abs_forCut__" + safe + "__" + mx2_col; // unique

  // Build axis and histogram titles safely (std::string on the left side)
  const bool kmiss = (kaon_label_for_axes == "K^{-}");
  const std::string parent_close = kmiss ? "K^{+})" : "K^{-})";
  const std::string xTitle = std::string("M_{X}(ep") + parent_close + " [GeV]";
  const std::string legendEntry = kaon_label_for_axes + " missing mass";

  // Always derive |Mx| from the provided Mx² column
  auto df = df_in
              .Define(mx_col, [](double mx2){ return (mx2>0) ? std::sqrt(mx2) : -999.0; }, {mx2_col})
              .Filter(mx_col + ">0", ("valid " + mx2_col).c_str());

  // Book hist
  auto h = df.Histo1D({("hMx_"+safe).c_str(), legendEntry.c_str(), nBins, xMin, xMax}, mx_col.c_str());

  MakeDir(outputDir);

  // Draw & fit
  const std::string outBase = outputDir + "/MM_"+safe;
  auto win = DrawFit_KaonMM(h.GetPtr(), outBase, tag, xMin, xMax, nSigma, xTitle.c_str(), legendEntry.c_str());

  auto df_cut = df.Filter([lo=win.low, hi=win.high](double x){ return x>lo && x<hi; }, {mx_col},
                          Form("Cut: %s in [%.3f, %.3f]", mx_col.c_str(), win.low, win.high));
  return {df_cut, win};
}


// -------------------------------------------------------------------------------------
// Generic invariant-mass from two 3-momenta columns (both already defined by InitKin)
// (labels purely cosmetic for plot titles / file names)
// -------------------------------------------------------------------------------------
inline ROOT::RDF::RNode DefineInvMassFromVecs(ROOT::RDF::RNode df_in,
                                              const std::string &name_out,
                                              const std::string &px1, const std::string &py1, const std::string &pz1,
                                              const std::string &px2, const std::string &py2, const std::string &pz2) {
  return df_in.Define(name_out,
    [](float px1, float py1, float pz1, float px2, float py2, float pz2) {
      const float mK = static_cast<float>(kMK);
      const float E1 = std::sqrt(px1*px1 + py1*py1 + pz1*pz1 + mK*mK);
      const float E2 = std::sqrt(px2*px2 + py2*py2 + pz2*pz2 + mK*mK);
      const float Ex = E1 + E2;
      const float px = px1 + px2, py = py1 + py2, pz = pz1 + pz2;
      const float m2 = Ex*Ex - (px*px + py*py + pz*pz);
      return (m2>0) ? std::sqrt(m2) : -1.0f;
    }, {px1,py1,pz1, px2,py2,pz2});
}

// -------------------------------------------------------------------------------------
// Phi mass fitter/drawer with signal and background components
// -------------------------------------------------------------------------------------
inline PhiMassWindow FitAndDrawPhiMassHistogram(TH1D* hist,
                                const std::string& title,
                                const std::string& outname,
                                TString xtitle = "M(K^{+}K^{-}) [GeV]",
                                TString ytitle = "Counts",
                                double fit_range_min = 0.9,
                                double fit_range_max = 1.25,
                                bool fit = true,
                                const std::string& tag = "",
                                double annotate_sigma_for_count = 3.0) {
  const std::string safe = MakeSafeName(tag);
  TCanvas* c = new TCanvas(("c_" + outname).c_str(), title.c_str(), 1200, 1000);
  StyleCanvas(c);
  StyleH1(hist, xtitle, ytitle, fit_range_min - .02, fit_range_max);
  hist->Draw("PE");

  TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist, "K^{+}K^{-} Inv. Mass", "lep");

  PhiMassWindow win;

  if (fit) {
    TF1* fitTotal = new TF1(("fitTotal_" + safe).c_str(),
                            "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874)) + [3]*TMath::Gaus(x,[4],[5])",
                            fit_range_min, fit_range_max);

    fitTotal->SetParameters(10, 0.9, 2, 100, 1.02, 0.010);
    fitTotal->SetParLimits(4, .90, 1.070);
    fitTotal->SetParLimits(5, 0.001, 0.3);
    fitTotal->SetLineColor(kRed + 1);
    fitTotal->SetLineWidth(3);
    hist->Fit(fitTotal, "R0Q");
    fitTotal->Draw("SAME C");

    double A = fitTotal->GetParameter(0);
    double alpha = fitTotal->GetParameter(1);
    double lambda = fitTotal->GetParameter(2);
    double N = fitTotal->GetParameter(3);
    win.mu = fitTotal->GetParameter(4);
    win.sigma = fitTotal->GetParameter(5);
    double chi2 = fitTotal->GetChisquare();
    double ndf = fitTotal->GetNDF();

    TF1* fitSignal = new TF1(("fitSignal_" + safe).c_str(), "[0]*TMath::Gaus(x,[1],[2])", fit_range_min, fit_range_max);
    fitSignal->SetParameters(N, win.mu, win.sigma);
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

    const auto lr = FitWindowFromGaussianCore(fitTotal, annotate_sigma_for_count);
    win.low = lr.first;
    win.high = lr.second;

    TLine* line1 = new TLine(win.low, 0, win.low, hist->GetMaximum() * 0.5);
    TLine* line2 = new TLine(win.high, 0, win.high, hist->GetMaximum() * 0.5);
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
    latex.DrawLatex(0.15, 0.87, Form("#mu = %.4f GeV", win.mu));
    latex.DrawLatex(0.15, 0.82, Form("#sigma = %.4f GeV", win.sigma));
    latex.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf = %.2f", chi2 / std::max(1.0, ndf)));
    latex.DrawLatex(0.15, 0.72, Form("Fit range: %.3f-%.3f GeV", fit_range_min, fit_range_max));

    if (annotate_sigma_for_count > 0 && std::isfinite(win.mu) && std::isfinite(win.sigma) && win.sigma > 0) {
      const double lo = win.mu - annotate_sigma_for_count * win.sigma;
      const double hi = win.mu + annotate_sigma_for_count * win.sigma;
      int binLo = hist->GetXaxis()->FindBin(lo);
      int binHi = hist->GetXaxis()->FindBin(hi);
      long long nCand = 0;
      for (int b = binLo; b <= binHi; ++b) nCand += (long long)std::llround(hist->GetBinContent(b));
      latex.DrawLatex(0.15, 0.67, Form("N_{|M-#mu|<%.0f#sigma} = %lld", annotate_sigma_for_count, nCand));
    }

    leg->AddEntry(fitTotal, "Total Fit: Thr.Exp + Gauss", "l");
    leg->AddEntry(fitSignal, "Signal (Gauss)", "f");
    leg->AddEntry(fitBkg, "Background (Thr.Exp)", "f");
    leg->AddEntry(line1, Form("#pm %.0f#sigma", annotate_sigma_for_count), "l");
  }

  leg->Draw();
  c->SaveAs((outname + ".png").c_str());
  c->SaveAs((outname + ".pdf").c_str());
  delete c;
  return win;
}

// -------------------------------------------------------------------------------------
// φ mass drawers (you pick the one that matches the InitKin path you used)
// -------------------------------------------------------------------------------------
inline std::tuple<ROOT::RDF::RNode, PhiMassWindow> DrawPhiMass_Measured(ROOT::RDF::RNode df_in,
                                 const std::string &outputDir, const std::string &tag,
                                 int nBins = 200, double xMin = 0.8, double xMax = 1.6, double fit_range_min = 0.9874, double fit_range_max = 1.12, double nSigma = 3.0) {
  const std::string safe = MakeSafeName(tag);
  auto df = DefineInvMassFromVecs(df_in, "Mphi_meas", "kPlus_px","kPlus_py","kPlus_pz",
                                                 "kMinus_px","kMinus_py","kMinus_pz")
              .Filter("Mphi_meas>0", "valid K+K- mass");
  auto h = df.Histo1D({("hPhi_meas_"+safe).c_str(), "K^{+}K^{-} invariant mass;M(K^{+}K^{-}) [GeV];Counts", nBins, xMin, xMax}, "Mphi_meas");

  MakeDir(outputDir);
  PhiMassWindow win = FitAndDrawPhiMassHistogram(h.GetPtr(), "Measured K+K- Inv. Mass", outputDir + "/Phi_meas_"+safe, "M(K^{+}K^{-}) [GeV]", "Counts", fit_range_min, fit_range_max, true, tag, nSigma);

  auto df_cut = df.Filter([lo=win.low, hi=win.high](double x){ return x>lo && x<hi; }, {"Mphi_meas"},
                          Form("Cut: Mphi_meas in [%.3f, %.3f]", win.low, win.high));

  return {df_cut, win};
}

inline std::tuple<ROOT::RDF::RNode, PhiMassWindow> DrawPhiMass_Kp_plus_KmMissing(ROOT::RDF::RNode df_in,
                                          const std::string &outputDir, const std::string &tag,
                                          int nBins = 200, double xMin = 0.8, double xMax = 1.6, double fit_range_min = 0.9874, double fit_range_max = 1.12, double nSigma = 3.0) {
  const std::string safe = MakeSafeName(tag);
  auto df = DefineInvMassFromVecs(df_in, "Mphi_Kp_KmMiss",
                                  "kPlus_px","kPlus_py","kPlus_pz",
                                  "kMinus_miss_px","kMinus_miss_py","kMinus_miss_pz")
              .Filter("Mphi_Kp_KmMiss>0", "valid M(K+ + K-_{miss})");
  auto h = df.Histo1D({("hPhi_KpKmMiss_"+safe).c_str(), "M(K^{+} + K^{-}_{miss}) ;M [GeV];Counts", nBins, xMin, xMax}, "Mphi_Kp_KmMiss");

  MakeDir(outputDir);
  PhiMassWindow win = FitAndDrawPhiMassHistogram(h.GetPtr(), "M(K+ + K-_{miss})", outputDir + "/Phi_KpKmMiss_"+safe, "M(K^{+}K^{-}_{miss}) [GeV]", "Counts", fit_range_min, fit_range_max, true, tag, nSigma);

  auto df_cut = df.Filter([lo=win.low, hi=win.high](double x){ return x>lo && x<hi; }, {"Mphi_Kp_KmMiss"},
                          Form("Cut: Mphi_Kp_KmMiss in [%.3f, %.3f]", win.low, win.high));

  return {df_cut, win};
}

inline std::tuple<ROOT::RDF::RNode, PhiMassWindow> DrawPhiMass_KpMissing_plus_Km(ROOT::RDF::RNode df_in,
                                          const std::string &outputDir, const std::string &tag,
                                          int nBins = 200, double xMin = 0.8, double xMax = 1.6, double fit_range_min = 0.9874, double fit_range_max = 1.12, double nSigma = 3.0) {
  const std::string safe = MakeSafeName(tag);
  auto df = DefineInvMassFromVecs(df_in, "Mphi_KpMiss_Km",
                                  "kPlus_miss_px","kPlus_miss_py","kPlus_miss_pz",
                                  "kMinus_px","kMinus_py","kMinus_pz")
              .Filter("Mphi_KpMiss_Km>0", "valid M(K+_{miss} + K-)");
  auto h = df.Histo1D({("hPhi_KpMissKm_"+safe).c_str(), "M(K^{+}_{miss} + K^{-}) ;M [GeV];Counts", nBins, xMin, xMax}, "Mphi_KpMiss_Km");

  MakeDir(outputDir);
  PhiMassWindow win = FitAndDrawPhiMassHistogram(h.GetPtr(), "M(K+_{miss} + K-)", outputDir + "/Phi_KpMissKm_"+safe, "M(K^{+}_{miss}K^{-}) [GeV]", "Counts", fit_range_min, fit_range_max, true, tag, nSigma);

  auto df_cut = df.Filter([lo=win.low, hi=win.high](double x){ return x>lo && x<hi; }, {"Mphi_KpMiss_Km"},
                          Form("Cut: Mphi_KpMiss_Km in [%.3f, %.3f]", win.low, win.high));

  return {df_cut, win};
}

// -------------------------------------------------------------------------------------
// API 2: φ→K⁺K⁻ mass plotter from node
// Prefers InitKinematics’ "invMass_KpKm" column, falls back to a local Define.
// -------------------------------------------------------------------------------------
inline std::tuple<ROOT::RDF::RNode, PhiMassWindow> DrawPhiMass_FromNode(ROOT::RDF::RNode df_in,
                                 const std::string& outputDir,
                                 const std::string& tag,
                                 double fit_range_min = 0.9874,
                                 double fit_range_max = 1.120,
                                 double nSigma = 3.0) {
  const std::string safe = MakeSafeName(tag);

  ROOT::RDF::RNode df = df_in;
  bool need_define = false;
  try {
    df = df.Define("invMass_KpKm_passthru", [](double v){ return v; }, {"invMass_KpKm"})
           .Filter("invMass_KpKm_passthru>0", "has K+K- pair");
  } catch (...) {
    need_define = true;
  }

  const std::string mass_col_name = need_define ? "invMass_KpKm" : "invMass_KpKm_passthru";
  
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
                      mass_col_name.c_str());

  MakeDir(outputDir);
  PhiMassWindow win = FitAndDrawPhiMassHistogram(h.GetPtr(),
                      "After kaon missing-mass window",
                      outputDir + "/KpKm_mass_afterCut_" + safe,
                      "M(K^{+}K^{-}) [GeV]", "Counts",
                      fit_range_min, fit_range_max,
                      /*fit=*/true, tag,
                      nSigma);

  auto df_cut = df.Filter([lo=win.low, hi=win.high](double x){ return x>lo && x<hi; }, {mass_col_name},
                          Form("Cut: %s in [%.3f, %.3f]", mass_col_name.c_str(), win.low, win.high));

  return {df_cut, win};
}

} // namespace PhiMass
} // namespace DISANA
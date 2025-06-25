#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TLegend.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>

using namespace ROOT;
using namespace ROOT::VecOps;

// ------------------------
// Style helper functions
// ------------------------
void ApplyHistStyle(TH1* hist) {
  hist->SetTitle("");
  hist->SetLineColor(kBlue + 1);
  hist->SetLineWidth(2);
  hist->SetFillColorAlpha(kBlue - 9, 0.3);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetMarkerColor(kBlue + 2);

  hist->GetXaxis()->SetTitle("M(#gamma#gamma) [GeV]");
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.0);

  hist->GetYaxis()->SetTitle("Counts");
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetRangeUser(0.05, 0.23);
}

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
}

// ------------------------
// Histogram Fitting & Drawing
// ------------------------
void FitAndDrawHistogram(TH1D* hist, const std::string& title, const std::string& outname) {
  TCanvas* c = new TCanvas(("c_" + outname).c_str(), title.c_str(), 1200, 1000);
  StyleCanvas(c);
  ApplyHistStyle(hist);
  hist->Draw("PE");

  // Fit: Total (Exp + Gauss)
  TF1* fitTotal = new TF1("fitTotal", "[0]*exp([1]*x) + [2]*TMath::Gaus(x,[3],[4])", 0.06, 0.23);
  fitTotal->SetParameters(100, -5, 300, 0.135, 0.01);
  fitTotal->SetLineColor(kRed + 1);
  fitTotal->SetLineWidth(3);
  hist->Fit(fitTotal, "R0");
  fitTotal->Draw("SAME C");

  // Extract parameters
  double A = fitTotal->GetParameter(0);
  double lambda = fitTotal->GetParameter(1);
  double N = fitTotal->GetParameter(2);
  double mu = fitTotal->GetParameter(3);
  double sigma = fitTotal->GetParameter(4);
  double chi2 = fitTotal->GetChisquare();
  double ndf  = fitTotal->GetNDF();

  // Signal (Gaussian)
  TF1* fitSignal = new TF1("fitSignal", "[0]*TMath::Gaus(x,[1],[2])", 0.05, 0.23);
  fitSignal->SetParameters(N, mu, sigma);
  fitSignal->SetLineColor(kOrange + 1);
  fitSignal->SetLineStyle(2);
  fitSignal->SetLineWidth(2);
  fitSignal->SetFillColorAlpha(kOrange - 3, 0.3);
  fitSignal->SetFillStyle(1001);
  fitSignal->Draw("SAME FC");

  // Background (Exp)
  TF1* fitBkg = new TF1("fitBkg", "[0]*exp([1]*x)", 0.05, 0.25);
  fitBkg->SetParameters(A, lambda);
  fitBkg->SetLineColor(kGreen + 2);
  fitBkg->SetLineStyle(3);
  fitBkg->SetLineWidth(2);
  fitBkg->SetFillColorAlpha(kGreen - 7, 0.3);
  fitBkg->SetFillStyle(1001);
  fitBkg->Draw("SAME FC");

  // Vertical lines ±3σ
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

  // Labels
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.15, 0.87, Form("#mu = %.3f GeV", mu));
  latex.DrawLatex(0.15, 0.82, Form("#sigma = %.4f GeV", sigma));
  latex.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf = %.2f", chi2 / ndf));
  latex.DrawLatex(0.15, 0.72, Form("Fit range: %.2f-%.2f GeV", fitTotal->GetXmin(), fitTotal->GetXmax()));



  // Legend
  TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist, "#gamma#gamma Inv. Mass", "lep");
  leg->AddEntry(fitTotal, "Total Fit: Exp + Gauss", "l");
  leg->AddEntry(fitSignal, "Signal (Gauss)", "f");
  leg->AddEntry(fitBkg, "Background (Exp)", "f");
  leg->AddEntry(line1, "#pm 3#sigma", "l");
  leg->Draw();

  // Save
  c->SaveAs((outname + ".png").c_str());
  c->SaveAs((outname + ".pdf").c_str());
  delete c;
}

// ------------------------
// Main Analysis Routine
// ------------------------
void DrawPi0Mass(const std::string& filename, const std::string& treename, const std::string& outputDir = "Pi0MassPlots") {
  TStopwatch timer;
  timer.Start();
  ROOT::EnableImplicitMT();

  RDataFrame df(treename, filename);
  auto df_filtered = df.Filter("REC_MotherMass.size() > 0")
    .Define("REC_DaughterParticle_pass_int", [](const std::vector<bool>& passVec) {
      return RVec<int>(passVec.begin(), passVec.end());
    }, {"REC_DaughterParticle_pass"});

  auto df_beforeCut = df_filtered.Define("MotherMass_all", [](const std::vector<float>& masses) {
    return RVec<float>(masses.begin(), masses.end());
  }, {"REC_MotherMass"});

  auto df_afterCut = df_filtered.Define("MotherMass_passed", [](const RVec<float>& masses, const RVec<int>& pass) {
    RVec<float> out;
    for (size_t i = 0; i < masses.size(); ++i)
      if (pass[i] == 1 && masses[i] > 0) out.push_back(masses[i]);
    return out;
  }, {"REC_MotherMass", "REC_DaughterParticle_pass_int"});

  auto h_before = df_beforeCut.Histo1D({"hBefore", "Before Mass Cut;M(#gamma#gamma) [GeV];Counts", 100, 0.05, 0.25}, "MotherMass_all");
  auto h_after  = df_afterCut.Histo1D({"hAfter",  "After Mass Cut;M(#gamma#gamma) [GeV];Counts", 100, 0.05, 0.25}, "MotherMass_passed");

  gSystem->Exec(("mkdir -p " + outputDir).c_str());
  FitAndDrawHistogram(h_before.GetPtr(), "Before InvMass Cut", outputDir + "/pi0_mass_before");
  FitAndDrawHistogram(h_after.GetPtr(),  "After InvMass Cut",  outputDir + "/pi0_mass_after");

  std::cout << "Plots saved in: " << outputDir << std::endl;
  timer.Stop();
  std::cout << "Time for DrawPi0Mass: ";
  timer.Print();
}

void analysisPi0Mass() {
  std::string path = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/test/";
  DrawPi0Mass(path + "dfSelected_afterFid.root", "dfSelected_afterFid");
  gApplication->Terminate(0);
}

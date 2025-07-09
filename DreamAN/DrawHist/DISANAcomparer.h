#ifndef DISANA_COMPARER_H
#define DISANA_COMPARER_H

// ROOT headers
#include <TCanvas.h>
#include <TLegend.h>

// STL headers
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project-specific headers
#include "DISANAplotter.h"
#include "DrawStyle.h"

namespace fs = std::filesystem;
// color palettes for different models
std::vector<std::tuple<double, double, double>> modelShades = {
    {0.20, 0.30, 0.85},  // Blue
    {0.90, 0.45, 0.10},  // Orange
    {0.00, 0.60, 0.60},  // Teal green
    {0.00, 0.70, 0.00},  // Green
    {0.60, 0.30, 0.80},  // Purple
    {0.85, 0.10, 0.25},  // Red
    {0.40, 0.40, 0.40}   // Gray (fallback)
};

class DISANAcomparer {
 public:
  // Set the bin ranges used for cross-section calculations and plotting
  void SetXBinsRanges(BinManager bins) { fXbins = bins; }

  void NormalizeHistogram(TH1* hist) {
    if (!hist) return;
    double integral = hist->Integral();
    if (integral > 0) hist->Scale(1.0 / integral);
  }
  // Add a new model with its DataFrame, label, and beam energy
  void AddModelwithPi0Corr(ROOT::RDF::RNode df_dvcs_data,
                ROOT::RDF::RNode df_pi0_data,
                ROOT::RDF::RNode df_dvcs_mc, 
                ROOT::RDF::RNode df_pi0_mc, 
                const std::string& label, double beamEnergy, bool fCorrection = false) {
    auto plotter = std::make_unique<DISANAplotter>(df_dvcs_data,
                                                   beamEnergy,
                                                   df_pi0_data,
                                                   df_dvcs_mc,
                                                   df_pi0_mc
                                                   );
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV with Pi0 Correction: " << fCorrection << std::endl;
    plotter->SetPlotApplyCorrection(fCorrection);
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("pho");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  void AddModel(ROOT::RDF::RNode df, const std::string& label, double beamEnergy) {
    auto plotter = std::make_unique<DISANAplotter>(df, beamEnergy);
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV without Pi0 Correction" << std::endl;
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("pho");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  // Set the output directory for saving plots
  void SetOutputDir(const std::string& outdir) {
    outputDir = outdir;
    if (!fs::exists(outputDir)) {
      fs::create_directories(outputDir);
    }
  }

  // Enable or disable individual variable plotting
  void PlotIndividual(bool plotInd) { plotIndividual = plotInd; }

  // Set plot styles for various plot types
  void SetKinStyle(const DrawStyle& style) { styleKin_ = style; }
  void SetDVCSStyle(const DrawStyle& style) { styleDVCS_ = style; }
  void SetCrossSectionStyle(const DrawStyle& style) { styleCrossSection_ = style; }
  void SetBSAStyle(const DrawStyle& style) { styleBSA_ = style; }

  // Enable or disable correctio
  void SetApplyCorrection(bool apply) { applyCorrection = apply; }

  // Load correction histogram from ROOT file
  void LoadCorrectionHistogram(const std::string& filename, const std::string& histoname = "h_correction") {
    correctionHist = nullptr;  // Reset after applying to avoid reusing the same histogram
    TFile* f = TFile::Open(filename.c_str(), "READ");

    if (!f || f->IsZombie()) {
      std::cerr << "Error: Cannot open correction file: " << filename << "\n";
      return;
    }

    correctionHist = dynamic_cast<THnSparseD*>(f->Get(histoname.c_str()));
    if (!correctionHist) {
      std::cerr << "Error: Correction histogram '" << histoname << "' not found in file: " << filename << "\n";
      return;
    }
    
    // correctionHist->SetDirectory(0);  // Detach from file
    f->Close();
    delete f; 
    std::cout << "✅ Correction histogram loaded: " << histoname << "\n";
  }

  // Plot all basic kinematic distributions (p, theta, phi) for all particle types
  void PlotKinematicComparison() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 3);

    std::vector<std::string> types = {"el", "pro", "pho"};
    std::vector<std::string> vars = {"p", "theta", "phi"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotVariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "KinematicComparison.pdf").c_str());

    // Optionally save individual plots
    if (plotIndividual) {
      for (const auto& type : types) {
        for (const auto& var : vars) {
          PlotSingleVariableComparison(type, var);
        }
      }
    }

    std::cout << "Saved kinematic comparison plots to: " << outputDir + "/KinematicComparison.pdf" << std::endl;
    delete canvas;
  }

  // Plot one specific variable (e.g., p) for a given particle type (e.g., "el")
  void PlotVariableComparison(const std::string& type, const std::string& var, int padIndex, TCanvas* canvas) {
    canvas->cd(padIndex);
    std::string hname_target = "rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;
    styleKin_.StylePad((TPad*)gPad);

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }
      NormalizeHistogram(target);
      styleKin_.StyleTH1(target);
      target->SetLineColor(i + 2);
      target->SetTitle(Form("%s;%s;Count", typeToParticle[type].c_str(), VarName[var].c_str()));

      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw();
  }

  // Save an individual variable comparison plot as PNG
  void PlotSingleVariableComparison(const std::string& type, const std::string& var) {
    TCanvas* canvas = new TCanvas(("c_" + type + "_" + var).c_str(), ("Comparison " + type + " " + var).c_str(), 800, 600);
    gPad->SetGrid();

    std::string hname_target = "rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotSingleVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }

      target->SetLineColor(i + 1);

      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw();
    canvas->Update();
    canvas->SaveAs((outputDir + "/compare_" + type + "_" + var + ".pdf").c_str());
    delete canvas;
  }

  // Plot DVCS-specific kinematic distributions and a 2D Q² vs xB histogram
  void PlotDVCSKinematicsComparison() {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleKin_.StylePad((TPad*)gPad);

      bool first = true;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto rdf = plotters[i]->GetRDF();

        if (!rdf.HasColumn(var)) {
          std::cerr << "[ERROR] Column " << var << " not found in RDF for model " << labels[i] << "\n";
          continue;
        }

        double min = *(rdf.Min(var));
        double max = *(rdf.Max(var));
        if (min == max) {
          min -= 0.1;
          max += 0.1;
        }

        double margin = std::max(1e-3, 0.05 * (max - min));
        auto h = rdf.Histo1D({Form("h_%s_%zu", var.c_str(), i), titles[var].c_str(), 100, min - margin, max + margin}, var);

        styleKin_.StylePad((TPad*)gPad);
        styleKin_.StyleTH1(h.GetPtr());
        NormalizeHistogram(h.GetPtr());
        h->SetLineColor(i + 2);

        // Set axis titles
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetXaxis()->SetNoExponent(true);  // ✅ Explicitly prevent scientific notation

        h->GetYaxis()->SetTitle("Counts");
        TGaxis::SetMaxDigits(1);  // ✅ Apply 10^n style for Y axis
        h->GetYaxis()->SetNoExponent(false);
        h->GetYaxis()->SetLabelFont(42);
        h->GetYaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetTitleOffset(1.2);
        if (first) {
          h->DrawCopy("HIST");
          first = false;
        } else {
          h->DrawCopy("HIST SAME");
        }
      }
    }

    // 🔳 Final pad: 2D Q² vs x_B plot
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 60, 0, 1.0, 60, 0, 10.0}, "xB", "Q2");

    styleKin_.StyleTH2(h2d.GetPtr());
    h2d->GetYaxis()->SetNoExponent(true);
    TGaxis::SetMaxDigits(1);
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.05);
    h2d->GetYaxis()->SetTitleOffset(1.2);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.05);
    h2d->GetXaxis()->SetTitleOffset(0.2);

    // ✅ Z-axis scientific formatting (palette)
    TGaxis::SetMaxDigits(1);
    h2d->DrawCopy("COLZ");
    gPad->Update();  // Needed to initialize palette

    TPaletteAxis* palette = (TPaletteAxis*)h2d->GetListOfFunctions()->FindObject("palette");
    if (palette) {
      palette->SetLabelFont(42);
      palette->SetLabelSize(0.035);
    }

    // Final save and cleanup
    canvas->SaveAs((outputDir + "/DVCS_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved DVCS kinematics comparison to: " << outputDir + "/DVCS_Kinematics_Comparison.pdf" << std::endl;
    delete canvas;

    // 🔁 Restore TGaxis global state
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  /// For exclusivity cuts, you can use the following function to select one triplet
  void PlotExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", -2.0, 2.0},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -2, 3.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -1.0, 1.0},
        {"Theta_gamma_gamma", "#theta(#gamma, #vec{q})", "#theta_{#gamma#gamma'} [deg]", -10.0, 30},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 90},
        {"Mx2_epg", "Missing Mass Squared (ep#gamma)", "MM^{2}(ep#gamma) [GeV^{2}]", -1.0, 1.0},
        {"Mx2_eg", "Invariant Mass (e#gamma)", "M^{2}(e#gamma) [GeV^{2}]", -5.5, 5.5},
        {"Theta_e_gamma", "Angle: e-#gamma", "#theta(e, #gamma) [deg]", 0.0, 180.0},
        {"DeltaE", "Energy Balance", "#DeltaE [GeV]", -2.0, 4.0},
    };

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 1800, 1200);
      int cols = 3;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_cut = plotters[m]->GetRDF().Filter(cutExpr, cutLabel);
          if (!rdf_cut.HasColumn(var)) continue;

          auto h = rdf_cut.Histo1D({Form("h_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
          h.GetValue();

          TH1D* h_clone = (TH1D*)h.GetPtr()->Clone();
          h_clone->SetDirectory(0);
          NormalizeHistogram(h_clone);

          styleKin_.StyleTH1(h_clone);
          h_clone->SetLineColor(m + 2);
          h_clone->SetLineWidth(2);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColor(m + 2);
          line2->SetLineColor(m + 2);
          line1->SetLineStyle(2);  // Dashed
          line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(2) << mean << ", #sigma = " << std::fixed << std::setprecision(2) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };

  // cross-section comparison plots
  void PlotDISCrossSectionComparison(double luminosity, double pol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }

    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allCSHists;
    for (auto& plotter : plotters) {
      auto hists = plotter->ComputeDVCS_CrossSection(fXbins, luminosity);
      allCSHists.push_back(std::move(hists));
    }

    const auto& q2_edges = fXbins.GetQ2Bins();
    const auto& t_edges = fXbins.GetTBins();
    const auto& xb_edges = fXbins.GetXBBins();

    const size_t n_q2 = q2_edges.size() - 1;
    const size_t n_t = t_edges.size() - 1;
    const size_t n_xb = xb_edges.size() - 1;

    const int rows = n_q2;
    const int cols = n_xb;

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      TString cname = Form("DIS_CrossSection_t[%zu]", t_bin);
      TCanvas* c = new TCanvas(cname, cname, 2200, 1600);
      c->Divide(cols, rows, 0.00, 0.00);  // No gaps
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          int visualRow = rows - 1 - q2_bin;
          int pad = visualRow * cols + xb_bin + 1;
          TPad* thisPad = (TPad*)c->cd(pad);

          // Pad margins
          double l = (xb_bin == 0) ? 0.22 : 0.00;
          double r = (xb_bin == cols - 1) ? 0.05 : 0.00;
          double b = (visualRow == rows - 1) ? 0.16 : 0.00;
          double t = (visualRow == 0) ? 0.001 : 0.00;

          styleBSA_.StylePad(thisPad, l, r, b, t);
          thisPad->SetLogy();
          thisPad->SetTicks(1, 0);

          const int idx = q2_bin * (n_t * n_xb) + t_bin * n_xb + xb_bin;

          TLegend* leg = new TLegend(0.35, 0.85, 0.85, 0.95);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.06);

          bool first = true;
          gStyle->SetCanvasPreferGL(true);

          for (size_t m = 0; m < allCSHists.size(); ++m) {
            TH1D* h = allCSHists[m][xb_bin][q2_bin][t_bin];
            styleBSA_.StyleTH1(h);

            auto [r, g, b] = modelShades[m % modelShades.size()];

            int colorIdx = 3000 + m * 20;  // Avoid low TColor indices

            if (!gROOT->GetColor(colorIdx)) {
              new TColor(colorIdx, r, g, b);
            }

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);
            h->SetLineWidth(1);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.6);
            h->SetStats(0);

            h->GetXaxis()->SetTitle((visualRow == rows - 1) ? "#phi [deg]" : "");
            h->GetYaxis()->SetTitle((xb_bin == 0) ? "d#sigma/d#phi [nb/deg]" : "");
            h->GetXaxis()->SetLabelSize((visualRow == rows - 1) ? 0.075 : 0.0);
            h->GetXaxis()->SetTitleSize((visualRow == rows - 1) ? 0.09 : 0.0);
            h->GetYaxis()->SetLabelSize((xb_bin == 0) ? 0.075 : 0.0);
            h->GetYaxis()->SetTitleSize((xb_bin == 0) ? 0.09 : 0.0);

            h->GetXaxis()->SetTitleOffset((visualRow == rows - 1) ? 0.8 : 0.0);
            h->GetYaxis()->SetTitleOffset((xb_bin == 0) ? 1.1 : 0.0);

            h->GetXaxis()->SetNdivisions(4, false);
            h->GetYaxis()->SetNdivisions(6, false);
           // h->GetYaxis()->SetRangeUser(-0.6, 0.6);
            h->GetXaxis()->SetRangeUser(0.0, 360.0);
            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);

            h->Draw(first ? "E1X0" : "E1X0 SAME");
            first = false;

            leg->AddEntry(h, labels[m].c_str(), "p");
          }

          // Annotate bin ranges
          double xB_low = xb_edges[xb_bin], xB_high = xb_edges[xb_bin + 1];
          double Q2_low = q2_edges[q2_bin], Q2_high = q2_edges[q2_bin + 1];
          TString labelText = Form("x_{B} #in [%.2f, %.2f], Q^{2} #in [%.1f, %.1f]", xB_low, xB_high, Q2_low, Q2_high);
          TLatex* latex = new TLatex(0.25, 0.82, labelText.Data());
          latex->SetTextSize(0.055);
          latex->SetNDC();
          latex->SetTextFont(42);
          latex->Draw();

          leg->Draw();
        }
      }

      TString outfile = Form("%s/DIS_CrossSection_t_%.2f-%.2f.pdf", outputDir.c_str(), t_edges[t_bin], t_edges[t_bin + 1]);
      c->SaveAs(outfile);
      std::cout << "Saved: " << outfile << '\n';
      delete c;
    }
  }

  /// BSA a plots
  void PlotDIS_BSA_Comparison(double luminosity, double pol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }

    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allBSA;
    for (auto& p : plotters) {
      auto h = p->ComputeBSA(fXbins, luminosity, pol);
      allBSA.push_back(std::move(h));
    }

    const auto& q2_edges = fXbins.GetQ2Bins();
    const auto& t_edges = fXbins.GetTBins();
    const auto& xb_edges = fXbins.GetXBBins();

    const size_t n_q2 = q2_edges.size() - 1;
    const size_t n_t = t_edges.size() - 1;
    const size_t n_xb = xb_edges.size() - 1;

    const int rows = n_q2;
    const int cols = n_xb;

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      TString cname = Form("DIS_BSA_t[%zu]", t_bin);
      TCanvas* c = new TCanvas(cname, cname, 2200, 1600);
      c->Divide(cols, rows, 0.00, 0.00);  // No gaps

      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          int visualRow = rows - 1 - q2_bin;
          int pad = visualRow * cols + xb_bin + 1;
          TPad* thisPad = (TPad*)c->cd(pad);

          // Pad margins
          double l = (xb_bin == 0) ? 0.2 : 0.00;
          double r = (xb_bin == cols - 1) ? 0.05 : 0.00;
          double b = (visualRow == rows - 1) ? 0.16 : 0.00;
          double t = (visualRow == 0) ? 0.001 : 0.00;

          styleBSA_.StylePad(thisPad, l, r, b, t);
          thisPad->SetTicks(1, 0);

          const int idx = q2_bin * (n_t * n_xb) + t_bin * n_xb + xb_bin;

          TLegend* leg = new TLegend(0.35, 0.85, 0.85, 0.95);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.06);

          TLegend* legParams = new TLegend(0.35, 0.16, 0.85, 0.32);  // Bottom legend for a₁
          legParams->SetBorderSize(0);
          legParams->SetFillStyle(0);
          legParams->SetTextSize(0.06);

          bool first = true;
          gStyle->SetCanvasPreferGL(true);

          for (size_t m = 0; m < allBSA.size(); ++m) {
            TH1D* h = allBSA[m][xb_bin][q2_bin][t_bin];
            styleBSA_.StyleTH1(h);

            auto [r, g, b] = modelShades[m % modelShades.size()];

            int colorIdx = 3000 + m * 20;  // Avoid low TColor indices

            if (!gROOT->GetColor(colorIdx)) {
              new TColor(colorIdx, r, g, b);
            }

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);
            h->SetLineWidth(1);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.5);
            h->SetStats(0);

            h->GetXaxis()->SetTitle((visualRow == rows - 1) ? "#phi [deg]" : "");
            h->GetYaxis()->SetTitle((xb_bin == 0) ? "A_{LU}" : "");
            h->GetXaxis()->SetLabelSize((visualRow == rows - 1) ? 0.075 : 0.0);
            h->GetXaxis()->SetTitleSize((visualRow == rows - 1) ? 0.09 : 0.0);
            h->GetYaxis()->SetLabelSize((xb_bin == 0) ? 0.075 : 0.0);
            h->GetYaxis()->SetTitleSize((xb_bin == 0) ? 0.09 : 0.0);

            h->GetXaxis()->SetTitleOffset((visualRow == rows - 1) ? 0.8 : 0.0);
            h->GetYaxis()->SetTitleOffset((xb_bin == 0) ? 1.1 : 0.0);

            h->GetXaxis()->SetNdivisions(4, false);
            h->GetYaxis()->SetNdivisions(6, false);
            h->GetYaxis()->SetRangeUser(-0.6, 0.6);
            h->GetXaxis()->SetRangeUser(0.01, 355);
            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);

            h->Draw(first ? "E1X0" : "E1X0 SAME");
            first = false;

            // Fit function and extract a₁
            TF1* fitFunc = new TF1(Form("fit_%zu_%zu_%zu_%zu", m, t_bin, q2_bin, xb_bin), "[0] + ([1]*sin(x*TMath::DegToRad())) / (1 + [2]*cos(x*TMath::DegToRad()))", 0, 360);
            fitFunc->SetParameters(0.0, 0.2, 0.1);
            fitFunc->SetFillColorAlpha(colorIdx, 0.5);
            fitFunc->SetLineColorAlpha(colorIdx, 0.5);
            fitFunc->SetLineStyle(2);
            fitFunc->SetLineWidth(1);
            h->Fit(fitFunc, "Q0");
            fitFunc->Draw("SAME");

            double a1 = fitFunc->GetParameter(1);
            double a1e = fitFunc->GetParError(1);
            TString a1label = Form("a_{1} = %.2f #pm %.2f", a1, a1e);
            legParams->AddEntry(fitFunc, a1label, "l");

            leg->AddEntry(h, labels[m].c_str(), "p");
          }

          // Annotate bin ranges
          double xB_low = xb_edges[xb_bin], xB_high = xb_edges[xb_bin + 1];
          double Q2_low = q2_edges[q2_bin], Q2_high = q2_edges[q2_bin + 1];
          TString labelText = Form("x_{B} #in [%.2f, %.2f], Q^{2} #in [%.1f, %.1f]", xB_low, xB_high, Q2_low, Q2_high);
          TLatex* latex = new TLatex(0.25, 0.82, labelText.Data());
          latex->SetTextSize(0.055);
          latex->SetNDC();
          latex->SetTextFont(42);
          latex->Draw();

          leg->Draw();
          legParams->Draw();
        }
      }

      TString outfile = Form("%s/DIS_BSA_t_%.2f-%.2f.pdf", outputDir.c_str(), t_edges[t_bin], t_edges[t_bin + 1]);
      c->SaveAs(outfile);
      std::cout << "Saved: " << outfile << '\n';

      delete c;
    }
  }
// pi0 correction plots for the BSA and Cross
  /// BSA a plots
  void PlotDIS_Pi0CorrComparison() {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }

    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0Corr;
    for (auto& p : plotters) {
      if (!p->getDoPi0Corr()) { allPi0Corr.emplace_back(); continue; }
      auto h = p->ComputePi0Corr(fXbins);
      allPi0Corr.push_back(std::move(h));
    }

    const auto& q2_edges = fXbins.GetQ2Bins();
    const auto& t_edges = fXbins.GetTBins();
    const auto& xb_edges = fXbins.GetXBBins();

    const size_t n_q2 = q2_edges.size() - 1;
    const size_t n_t = t_edges.size() - 1;
    const size_t n_xb = xb_edges.size() - 1;

    const int rows = n_q2;
    const int cols = n_xb;

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      TString cname = Form("DIS_pi0Corr_t[%zu]", t_bin);
      TCanvas* c = new TCanvas(cname, cname, 2200, 1600);
      c->Divide(cols, rows, 0.00, 0.00);  // No gaps

      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          int visualRow = rows - 1 - q2_bin;
          int pad = visualRow * cols + xb_bin + 1;
          TPad* thisPad = (TPad*)c->cd(pad);

          // Pad margins
          double l = (xb_bin == 0) ? 0.2 : 0.00;
          double r = (xb_bin == cols - 1) ? 0.05 : 0.00;
          double b = (visualRow == rows - 1) ? 0.16 : 0.00;
          double t = (visualRow == 0) ? 0.001 : 0.00;

          styleBSA_.StylePad(thisPad, l, r, b, t);
          thisPad->SetTicks(1, 0);

          const int idx = q2_bin * (n_t * n_xb) + t_bin * n_xb + xb_bin;

          TLegend* leg = new TLegend(0.35, 0.85, 0.85, 0.95);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.06);

          bool first = true;
          gStyle->SetCanvasPreferGL(true);

          for (size_t m = 0; m < allPi0Corr.size(); ++m) {
            if (allPi0Corr[m].empty())  continue;

            TH1D* h = allPi0Corr[m][xb_bin][q2_bin][t_bin];
            if (!h) continue;
            styleBSA_.StyleTH1(h);

            auto [r, g, b] = modelShades[m % modelShades.size()];

            int colorIdx = 3000 + m * 20;  // Avoid low TColor indices

            if (!gROOT->GetColor(colorIdx)) {
              new TColor(colorIdx, r, g, b);
            }

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);
            h->SetLineWidth(1);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.5);
            h->SetStats(0);

            h->GetXaxis()->SetTitle((visualRow == rows - 1) ? "#phi [deg]" : "");
            h->GetYaxis()->SetTitle((xb_bin == 0) ? "#eta^{#pi^{0}}" : "");
            h->GetXaxis()->SetLabelSize((visualRow == rows - 1) ? 0.075 : 0.0);
            h->GetXaxis()->SetTitleSize((visualRow == rows - 1) ? 0.09 : 0.0);
            h->GetYaxis()->SetLabelSize((xb_bin == 0) ? 0.075 : 0.0);
            h->GetYaxis()->SetTitleSize((xb_bin == 0) ? 0.09 : 0.0);

            h->GetXaxis()->SetTitleOffset((visualRow == rows - 1) ? 0.8 : 0.0);
            h->GetYaxis()->SetTitleOffset((xb_bin == 0) ? 1.1 : 0.0);

            h->GetXaxis()->SetNdivisions(4, false);
            h->GetYaxis()->SetNdivisions(6, false);
            h->GetYaxis()->SetRangeUser(-0.05, 0.55);
            h->GetXaxis()->SetRangeUser(0.0, 360);
            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);

            h->Draw(first ? "E1X0" : "E1X0 SAME");
            first = false;

            leg->AddEntry(h, labels[m].c_str(), "p");
          }

          // Annotate bin ranges
          double xB_low = xb_edges[xb_bin], xB_high = xb_edges[xb_bin + 1];
          double Q2_low = q2_edges[q2_bin], Q2_high = q2_edges[q2_bin + 1];
          TString labelText = Form("x_{B} #in [%.2f, %.2f], Q^{2} #in [%.1f, %.1f]", xB_low, xB_high, Q2_low, Q2_high);
          TLatex* latex = new TLatex(0.25, 0.82, labelText.Data());
          latex->SetTextSize(0.055);
          latex->SetNDC();
          latex->SetTextFont(42);
          latex->Draw();

          leg->Draw();
        }
      }

      TString outfile = Form("%s/DIS_pi0Corr_t_%.2f-%.2f.pdf", outputDir.c_str(), t_edges[t_bin], t_edges[t_bin + 1]);
      c->SaveAs(outfile);
      std::cout << "Saved: " << outfile << '\n';

      delete c;
    }
  }

  

 private:
  BinManager fXbins;
  bool plotIndividual = false;

  DrawStyle style_;              // Default style
  DrawStyle styleKin_;           // Kin plot style
  DrawStyle styleDVCS_;          // DVCS plot style
  DrawStyle styleCrossSection_;  // Cross-section plot style
  DrawStyle styleBSA_;           // BSA plot style

  bool applyCorrection = false;
  THnSparseD* correctionHist = nullptr;

  std::unique_ptr<ROOT::RDF::RNode> rdf;
  std::string outputDir = ".";

  std::vector<std::unique_ptr<DISANAplotter>> plotters;
  std::vector<std::string> labels;

  std::vector<std::string> particleName = {"e", "p", "#gamma"};
  std::map<std::string, std::string> typeToParticle = {{"el", "electron"}, {"pro", "proton"}, {"pho", "#gamma"}};
  std::map<std::string, std::string> VarName = {{"p", "p (GeV/#it{c})"}, {"theta", "#theta (rad)"}, {"phi", "#phi(rad)"}};
};
#endif  // DISANA_COMPARER_H

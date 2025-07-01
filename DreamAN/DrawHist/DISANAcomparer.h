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
  void AddModel(ROOT::RDF::RNode df, const std::string& label, double beamEnergy, bool applyCorrection = false, const std::string& correctionFile = "",
                const std::string& histoname = "h_correction") {
    auto plotter = std::make_unique<DISANAplotter>(df, beamEnergy);
    if (applyCorrection) {
      LoadCorrectionHistogram(correctionFile, histoname);
      if (!correctionHist) {
        std::cerr << "Error: Correction histogram not loaded. Skipping background correction.\n";
      } else {
        std::cout << "Applying background correction for model: " << label << "\n";
        plotter->ApplyPi0BkgCorr(correctionHist);
      }
    }
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

  // Compare DIS cross sections across models in Q^2-xB-t bins
  void PlotDISCrossSectionComparison(double luminosity) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }

    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allCSHists;
    for (auto& plotter : plotters) {
      auto hists = plotter->ComputeDVCS_CrossSection(fXbins, luminosity);
      allCSHists.push_back(std::move(hists));
    }

    const auto& q2_bins = fXbins.GetQ2Bins();
    const auto& xb_bins = fXbins.GetXBBins();
    const auto& t_bins = fXbins.GetTBins();

    size_t n_q2 = q2_bins.size() - 1;
    size_t n_xb = xb_bins.size() - 1;
    size_t n_t = t_bins.size() - 1;

    for (size_t it = 0; it < n_t; ++it) {
      TString canvasName = Form("DIS_CS_Comparison_tbin_%zu", it);
      TString canvasTitle = Form("DIS Cross-Section Comparison (t-bin %.1f-%.1f)", t_bins[it], t_bins[it + 1]);
      TCanvas* canvas = new TCanvas(canvasName, canvasTitle, 1800, 1200);
      canvas->Divide(n_xb, n_q2);

      for (size_t iq = 0; iq < n_q2; ++iq) {
        for (size_t ix = 0; ix < n_xb; ++ix) {
          int pad_idx = iq * n_xb + ix + 1;
          canvas->cd(pad_idx);
          styleCrossSection_.StylePad((TPad*)gPad);

          TLegend* legend = new TLegend(0.6, 0.7, 0.7, 0.8);
          legend->SetBorderSize(0);
          legend->SetFillStyle(0);
          legend->SetTextSize(0.035);

          bool first = true;

          for (size_t m = 0; m < allCSHists.size(); ++m) {
            TH1D* hist = allCSHists[m][ix][iq][it];
            if (!hist) continue;
            styleCrossSection_.StyleTH1(hist);
            hist->SetLineColor(static_cast<int>(m + 2));
            hist->SetLineWidth(2);
            hist->SetMarkerStyle(20);
            hist->SetStats(0);
            hist->SetMarkerColor(static_cast<int>(m + 2));
            hist->SetMarkerSize(2);
            hist->SetStats(0);
            hist->GetXaxis()->SetTitle("#phi [deg]");
            hist->GetYaxis()->SetTitle("d#sigma/d#phi [nb/deg]");

            if (first) {
              hist->Draw("E1X0");
              first = false;
            } else {
              hist->Draw("E1X0 SAME");
            }

            legend->AddEntry(hist, labels[m].c_str(), "l");
          }

          double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
          double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];
          double tmin = t_bins[it], tmax = t_bins[it + 1];

          TString labelText = Form("Q^{2} #in [%.1f-%.1f], x_{B} #in [%.2f-%.2f]", qmin, qmax, xbmin, xbmax);
          TLatex* label = new TLatex(0.45, 0.85, labelText.Data());
          label->SetNDC();
          label->SetTextSize(0.04);
          label->Draw();

          legend->Draw();
        }
      }

      TString filename = outputDir + Form("/DISCrossSectionComparison_t%.0f.pdf", 1000 * t_bins[it]);
      canvas->SaveAs(filename);
      std::cout << "Saved canvas to " << filename << std::endl;
      delete canvas;
    }

    allCSHists.clear();
  }

  void PlotDIS_BSA_Comparison(double luminosity, double pol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }
    std::vector<std::vector<TH1D*>> allBSA;
    for (auto& p : plotters) {
      auto h = p->ComputeBSA(fXbins, luminosity, pol);  
      if (h.empty()) std::cerr << "Warning: empty BSA histos.\n";
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
      TCanvas* c = new TCanvas(cname, cname, 2000, 1500);
      c->Divide(cols, rows);

      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          const int pad = q2_bin * cols + xb_bin + 1;
          c->cd(pad);
          styleBSA_.StylePad((TPad*)gPad);

          const int idx = q2_bin * (n_t * n_xb) + t_bin * n_xb + xb_bin;

          TLegend* leg =  new TLegend(0.6, 0.7, 0.7, 0.8);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.04);

          bool first = true;
          for (size_t m = 0; m < allBSA.size(); ++m) {
            TH1D* h = allBSA[m][idx];
            styleBSA_.StyleTH1(h);
            h->SetLineColor(static_cast<int>(m + 2));
            h->SetLineWidth(2);
            h->SetMarkerStyle(20);
            h->SetStats(0);
            h->SetMarkerColor(static_cast<int>(m + 2));
            h->SetMarkerSize(2);
            h->SetStats(0);

            h->GetXaxis()->SetTitle("#phi [deg]");
            h->GetYaxis()->SetTitle("A_{LU}");

            h->Draw(first ? "E1X0" : "E1X0 SAME");
            first = false;
            leg->AddEntry(h, labels[m].c_str(), "l");
          }

          double xB_low = xb_edges[xb_bin], xB_high = xb_edges[xb_bin + 1];
          double Q2_low = q2_edges[q2_bin], Q2_high = q2_edges[q2_bin + 1];
          double t_low = t_edges[t_bin], t_high = t_edges[t_bin + 1];

       
          leg->Draw();
          TString labelText = Form("Q^{2} #in [%.1f-%.1f], x_{B} #in [%.2f-%.2f]", Q2_low, Q2_high, xB_low, xB_high);
          TLatex* label = new TLatex(0.45, 0.85, labelText.Data());
          label->SetNDC();
          label->SetTextSize(0.04);
          label->Draw();
        }
      }

      c->SetLogy(false);
      TString outfile = Form("%s/DIS_BSA_t_%.2f-%.2f.pdf", outputDir.c_str(), t_edges[t_bin], t_edges[t_bin + 1]);
      c->SaveAs(outfile);
      std::cout << "Saved: " << outfile << '\n';

      delete c;
    }
  }

  ///
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

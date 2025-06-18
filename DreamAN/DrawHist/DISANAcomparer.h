#ifndef DISANA_COMPARER_H
#define DISANA_COMPARER_H

// ROOT headers
#include <TCanvas.h>
#include <TLegend.h>

// STL headers
#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <iostream>

// Project-specific headers
#include "DISANAplotter.h"
#include "DrawStyle.h"

namespace fs = std::filesystem;

class DISANAcomparer {
 public:
  // Set the bin ranges used for cross-section calculations and plotting
  void SetXBinsRanges(BinManager bins) { fXbins = bins; }

  // Add a new model with its DataFrame, label, and beam energy
  void AddModel(ROOT::RDF::RNode df, const std::string& label, double beamEnergy) {
    auto plotter = std::make_unique<DISANAplotter>(df, beamEnergy);
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
    canvas->SaveAs((outputDir + "KinematicComparison.png").c_str());

    // Optionally save individual plots
    if (plotIndividual) {
      for (const auto& type : types) {
        for (const auto& var : vars) {
          PlotSingleVariableComparison(type, var);
        }
      }
    }

    std::cout << "Saved kinematic comparison plots to: " << outputDir + "/KinematicComparison.png" << std::endl;
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
    canvas->SaveAs((outputDir + "/compare_" + type + "_" + var + ".png").c_str());
    delete canvas;
  }

  // Plot DVCS-specific kinematic distributions and a 2D Q² vs xB histogram
  void PlotDVCSKinematicsComparison() {
    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {
      {"Q2", "Q^{2} [GeV^{2}]"},
      {"xB", "x_{B}"},
      {"t", "-t [GeV^{2}]"},
      {"W", "W [GeV]"},
      {"phi", "#phi [deg]"}
    };

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);

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

        styleDVCS_.StylePad((TPad*)gPad);
        styleDVCS_.StyleTH1(h.GetPtr());
        h->SetLineColor(i + 2);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Counts");

        if (first) {
          h->DrawCopy("HIST");
          first = false;
        } else {
          h->DrawCopy("HIST SAME");
        }
      }
    }

    // Final pad: 2D Q² vs x_B for the first model
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 60, 0, 1.0, 60, 0, 10.0}, "xB", "Q2");
    styleDVCS_.StyleTH2(h2d.GetPtr());
    h2d->DrawCopy("COLZ");

    canvas->SaveAs((outputDir + "/DVCS_Kinematics_Comparison.png").c_str());
    std::cout << "Saved DVCS kinematics comparison to: " << outputDir + "/DVCS_Kinematics_Comparison.png" << std::endl;
    delete canvas;
  }

  // Compare DIS cross sections across models in Q^2-xB-t bins
  void PlotDISCrossSectionComparison(double luminosity) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }

    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allCSHists;
    for (auto& plotter : plotters) {
      auto hists = plotter->ComputeDVCS_CrossSection2(fXbins, luminosity);
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
      TCanvas* canvas = new TCanvas(canvasName, canvasTitle, 1800, 1600);
      canvas->Divide(n_xb, n_q2);

      for (size_t iq = 0; iq < n_q2; ++iq) {
        for (size_t ix = 0; ix < n_xb; ++ix) {
          int pad_idx = iq * n_xb + ix + 1;
          canvas->cd(pad_idx);
          styleCrossSection_.StylePad((TPad*)gPad);

          TLegend* legend = new TLegend(0.45, 0.7, 0.85, 0.88);
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

          TString labelText = Form("Q^{2} #in [%.1f-%.1f], x_{B} #in [%.1f-%.1f]", qmin, qmax, xbmin, xbmax);
          TLatex* label = new TLatex(0.45, 0.85, labelText.Data());
          label->SetNDC();
          label->SetTextSize(0.04);
          label->Draw();

          legend->Draw();
        }
      }

      TString filename = outputDir + Form("/DISCrossSectionComparison_t%.0f.png", 1000 * t_bins[it]);
      canvas->SaveAs(filename);
      std::cout << "Saved canvas to " << filename << std::endl;
      delete canvas;
    }

    allCSHists.clear();
  }

  // Compare BSA (Beam Spin Asymmetry) distributions across models
  void PlotDIS_BSA_Comparison(double luminosity) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }

    std::vector<std::vector<TH1D*>> allCSHists;
    for (int i = 0; i < plotters.size(); ++i) {
      auto& plotter = plotters[i];
      auto histos = plotter->ComputeBSA(fXbins, luminosity);
      if (histos.empty()) {
        std::cerr << "Warning: Empty BSA histograms for one of the models.\n";
      }
      allCSHists.push_back(histos);
    }

    size_t n_q2 = fXbins.GetQ2Bins().size() - 1;
    size_t n_t = fXbins.GetTBins().size() - 1;
    size_t n_xb = fXbins.GetXBBins().size() - 1;

    int cols = n_t * n_xb;
    int rows = n_q2;
    int totalBins = n_q2 * n_t * n_xb;

    TCanvas* canvas = new TCanvas("DIS_BSA_Comparison", "DIS BSA Comparison", 2000, 1500);
    canvas->Divide(cols, rows);

    for (int idx = 0; idx < totalBins; ++idx) {
      canvas->cd(idx + 1);
      styleBSA_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.45, 0.7, 0.8, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.04);

      int xb_bin = idx % n_xb;
      int t_bin = (idx / n_xb) % n_t;
      int q2_bin = idx / (n_xb * n_t);

      double xB_low = fXbins.GetXBBins()[xb_bin];
      double xB_high = fXbins.GetXBBins()[xb_bin + 1];
      double Q2_low = fXbins.GetQ2Bins()[q2_bin];
      double Q2_high = fXbins.GetQ2Bins()[q2_bin + 1];

      bool first = true;
      for (size_t m = 0; m < allCSHists.size(); ++m) {
        TH1D* hist = allCSHists[m][idx];
        styleBSA_.StyleTH1(hist);
        hist->SetLineColor(static_cast<int>(m + 2));
        hist->SetLineWidth(1);
        hist->GetXaxis()->SetTitle("#phi [deg]");
        hist->GetYaxis()->SetTitle("A_{LU}");

        if (first) {
          hist->Draw("E1X0");
          first = false;
        } else {
          hist->Draw("E1X0 SAME");
        }

        legend->AddEntry(hist, labels[m].c_str(), "l");
      }

      TString labelText = Form("Q^{2} #in [%.1f-%.1f] GeV^{2},  x_{B} #in [%.1f-%.1f]", Q2_low, Q2_high, xB_low, xB_high);
      TLatex* label = new TLatex(0.45, 0.85, labelText.Data());
      label->SetNDC();
      label->SetTextSize(0.04);
      label->Draw();
      legend->Draw();
    }

    canvas->SetLogy(false);
    canvas->SaveAs((outputDir + "/DIS_BSA_Comparison.png").c_str());
    std::cout << "Saved DIS BSA comparison to: " << outputDir + "/DIS_BSA_Comparison.png" << std::endl;

    allCSHists.clear();
    delete canvas;
  }


 private:
  BinManager fXbins;
  bool plotIndividual = false;

  DrawStyle style_;               // Default style
  DrawStyle styleKin_;           // Kin plot style
  DrawStyle styleDVCS_;           // DVCS plot style
  DrawStyle styleCrossSection_;   // Cross-section plot style
  DrawStyle styleBSA_;            // BSA plot style

  std::unique_ptr<ROOT::RDF::RNode> rdf;
  std::string outputDir = ".";

  std::vector<std::unique_ptr<DISANAplotter>> plotters;
  std::vector<std::string> labels;

  std::vector<std::string> particleName = {"e", "p", "#gamma"};
  std::map<std::string, std::string> typeToParticle = {
    {"el", "electron"},
    {"pro", "proton"},
    {"pho", "#gamma"}
  };
  std::map<std::string, std::string> VarName = {
    {"p", "p (GeV/#it{c})"},
    {"theta", "#theta (rad)"},
    {"phi", "#phi(rad)"}
  };
};

#endif  // DISANA_COMPARER_H

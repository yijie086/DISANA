#ifndef DISANA_COMPARER_H
#define DISANA_COMPARER_H

// ROOT headers
#include <TCanvas.h>
#include <TBox.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TLine.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TVirtualFitter.h>
#include <ROOT/RDFHelpers.hxx>
#include <sys/stat.h>
// STL headers
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project-specific headers
#include <chrono>

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
  // Shared display binning for He3DIS spectator distributions and acceptance.
  // These ranges are NOT event-selection cuts. pT and alpha use the q axis.
  struct He3DISSpectatorPlotSpec {
    std::string variable, axis;
    int bins;
    double lo, hi;
  };
  static std::vector<He3DISSpectatorPlotSpec> He3DISSpectatorPlotSpecs() {
    return {{"p_sum", "|#Sigma #vec{p}_{spectator}| [GeV]", 50, 0., 0.5},
            {"alpha", "#alpha_{S} = 3(E_{S}-p_{L,S})/M_{^{3}He}", 50, 1.5, 2.5},
            {"pT", "p_{T,S} relative to #vec{q} [GeV]", 50, 0., 0.2}};
  }
  struct He3DISModel {
    std::unique_ptr<ROOT::RDF::RNode> all;
    std::unique_ptr<ROOT::RDF::RNode> epp;
    std::unique_ptr<ROOT::RDF::RNode> epn;
    std::unique_ptr<ROOT::RDF::RNode> ed;
    std::string label;
  };

 public:
  void SetHe3DISXActiveBinning(int bins,double xMin,double xMax,
                               double modelXMax) {
    if(bins<1 || !std::isfinite(xMin) || !std::isfinite(xMax) ||
       !std::isfinite(modelXMax) || !(xMin<xMax) ||
       !(modelXMax>xMin && modelXMax<=xMax))
      throw std::invalid_argument(
          "He3DIS x_active binning requires bins>=1, xMin<xMax, and "
          "xMin<modelXMax<=xMax");
    he3DISXActiveBins=bins;
    he3DISXActiveMin=xMin;
    he3DISXActiveMax=xMax;
    he3DISModelXMax=modelXMax;
  }

  // Set the bin ranges used for cross-section calculations and plotting
  void SetXBinsRanges(BinManager bins) { fXbins = bins; }
  void SetUseSDHEP(bool enabled = true) { fXbins.SetUseSDHEP(enabled); }
  bool UsesSDHEP() const { return fXbins.UsesSDHEP(); }

  void SetDVCSWeightFunction(DVCSWeightFunction weightFunc) {
    dvcs_weight_function_ = std::move(weightFunc);
    for (auto& plotter : plotters) {
      plotter->SetDVCSWeightFunction(dvcs_weight_function_);
    }
  }

  void NormalizeHistogram(TH1* hist) {
    if (!hist) return;
    double integral = hist->Integral();
    if (integral > 0) hist->Scale(1.0 / integral);
  }
  // New interface: radiative and P1 corrections may use different
  // no-radiation reference samples.
  void AddModelwithPi0Corr(ROOT::RDF::RNode df_dvcs_data, ROOT::RDF::RNode df_pi0_data, ROOT::RDF::RNode df_dvcs_pi0mc, ROOT::RDF::RNode df_pi0_pi0mc,
                           ROOT::RDF::RNode df_gen_dvcsmc, ROOT::RDF::RNode df_accept_dvcsmc, ROOT::RDF::RNode df_dvcsmc_bkg, ROOT::RDF::RNode df_dvcsmc_nobkg,
                           ROOT::RDF::RNode df_dvcsmc_rad, ROOT::RDF::RNode df_dvcsmc_norad_rad, ROOT::RDF::RNode df_dvcsmc_norad_p1,
                           ROOT::RDF::RNode df_dvcsmc_p1cut, const std::string& label, double beamEnergy,
                           bool fPi0Correction = false, bool fAcceptanceCorrection = false, bool fEfficiencyCorrection = false, bool fRadiativeCorrection = false,
                           bool fP1cut = false, double luminosity = 1.0, double I_avg = 60.0, double I_mc = 60.0, double eff_corr = 1.0) {
    auto plotter = std::make_unique<DISANAplotter>(DVCSModeTag{}, df_dvcs_data, beamEnergy, luminosity, I_avg, I_mc, eff_corr, df_pi0_data, df_dvcs_pi0mc, df_pi0_pi0mc, df_gen_dvcsmc, df_accept_dvcsmc,
                                                   df_dvcsmc_bkg, df_dvcsmc_nobkg, df_dvcsmc_rad, df_dvcsmc_norad_rad, df_dvcsmc_p1cut, df_dvcsmc_norad_p1);
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV with Pi0 Correction: " << fPi0Correction
              << ", Acceptance Correction: " << fAcceptanceCorrection << ", Background Merging efficiency: " << fEfficiencyCorrection
              << ", Radiative Correction: " << fRadiativeCorrection << ", P1 cut: " << fP1cut << std::endl;
    plotter->SetPlotApplyCorrection(fPi0Correction);
    plotter->SetPlotApplyAcceptanceCorrection(fAcceptanceCorrection);
    plotter->SetPlotApplyEfficiencyCorrection(fEfficiencyCorrection);
    plotter->SetPlotApplyRadiativeCorrection(fRadiativeCorrection);
    plotter->SetPlotApplyP1Cut(fP1cut);
    plotter->SetDVCSWeightFunction(dvcs_weight_function_);
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("pho");
    plotter->GeneratePi0KinematicHistos("el");
    plotter->GeneratePi0KinematicHistos("pro");
    plotter->GeneratePi0KinematicHistos("pho");
    plotter->GeneratePi0KinematicHistos("pho2");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
    WriteDVCSPhysicalVolumeCSV(beamEnergy);
  }

  // Backward-compatible interface used by existing macros.  Its behavior is
  // unchanged: the same no-radiation sample is used for both corrections.
  void AddModelwithPi0Corr(ROOT::RDF::RNode df_dvcs_data, ROOT::RDF::RNode df_pi0_data, ROOT::RDF::RNode df_dvcs_pi0mc, ROOT::RDF::RNode df_pi0_pi0mc,
                           ROOT::RDF::RNode df_gen_dvcsmc, ROOT::RDF::RNode df_accept_dvcsmc, ROOT::RDF::RNode df_dvcsmc_bkg, ROOT::RDF::RNode df_dvcsmc_nobkg,
                           ROOT::RDF::RNode df_dvcsmc_rad, ROOT::RDF::RNode df_dvcsmc_norad, ROOT::RDF::RNode df_dvcsmc_p1cut, const std::string& label, double beamEnergy,
                           bool fPi0Correction = false, bool fAcceptanceCorrection = false, bool fEfficiencyCorrection = false, bool fRadiativeCorrection = false,
                           bool fP1cut = false, double luminosity = 1.0, double I_avg = 60.0, double I_mc = 60.0, double eff_corr = 1.0) {
    AddModelwithPi0Corr(
        df_dvcs_data, df_pi0_data, df_dvcs_pi0mc, df_pi0_pi0mc,
        df_gen_dvcsmc, df_accept_dvcsmc, df_dvcsmc_bkg, df_dvcsmc_nobkg,
        df_dvcsmc_rad, df_dvcsmc_norad, df_dvcsmc_norad,
        df_dvcsmc_p1cut, label, beamEnergy,
        fPi0Correction, fAcceptanceCorrection, fEfficiencyCorrection,
        fRadiativeCorrection, fP1cut, luminosity, I_avg, I_mc, eff_corr);
  }

  void AddModel(ROOT::RDF::RNode df, const std::string& label, double beamEnergy, double luminosity = 1.0) {
    auto plotter = std::make_unique<DISANAplotter>(DVCSModeTag{}, df, beamEnergy, luminosity);
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV without Pi0 Correction" << std::endl;
    plotter->SetDVCSWeightFunction(dvcs_weight_function_);
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("pho");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
    WriteDVCSPhysicalVolumeCSV(beamEnergy);
  }

  // Register a generator-level He3 MC dataframe.  This intentionally does
  // not construct a DVCS/Phi DISANAplotter or book any histograms: He3DIS
  // currently provides only the inclusive electron plus nuclear spectators.
  // He3-specific comparison and plotting behavior can be added here later
  // without changing the calling macros.
  void AddModelHe3(ROOT::RDF::RNode df_mc) {
    he3Models.push_back(
        std::make_unique<ROOT::RDF::RNode>(std::move(df_mc)));
    std::cout << "Adding He3 MC model (no histograms booked). Total He3 models: "
              << he3Models.size() << std::endl;
  }

  // Register coherent e' gamma He3' generator-level events separately from
  // both inclusive He3DIS and the proton-DVCS DISANAplotter.
  void AddModelHe3DVCS(ROOT::RDF::RNode df_mc) {
    const auto eventCount = df_mc.Count().GetValue();
    he3DVCSModels.push_back(
        std::make_unique<ROOT::RDF::RNode>(std::move(df_mc)));
    he3DVCSModelEventCounts.push_back(eventCount);
    std::cout << "Adding coherent He3 DVCS MC model. Total He3DVCS models: "
              << he3DVCSModels.size() << ", selected events: "
              << FormatMillions(eventCount) << std::endl;
  }

  // Overlay GEN and REC He3DVCS distributions with the same binning and
  // independent unit-area normalization. Model 0 is GEN and model 1 is REC.
  // Variables that require a reconstructed He3 are still booked for REC; if
  // no He3 exists their sentinel values stay outside the displayed range.
  void PlotNormalizedComparisonHe3DVCS(bool drawFakeHe3Rec = false) {
    if (he3DVCSModels.size() < 2) {
      std::cerr << "[PlotNormalizedComparisonHe3DVCS] requires GEN and REC models.\n";
      return;
    }
    struct Spec {
      std::string var, title, xTitle;
      int bins;
      double lo, hi;
      bool requiresHe3 = false;
    };
    const std::vector<Spec> particles = {
      {"recel_p", "Electron momentum", "p_{e'} [GeV]", 120, 0, 11},
      {"recel_theta", "Electron polar angle", "#theta_{e'} [deg]", 60, 0, 30},
      {"recel_phi", "Electron azimuth", "#phi_{e'} [deg]", 60, 0, 360},
      {"recpho_p", "Photon momentum", "p_{#gamma} [GeV]", 120, 0, 11},
      {"recpho_theta", "Photon polar angle", "#theta_{#gamma} [deg]", 60, 0, 20},
      {"recpho_phi", "Photon azimuth", "#phi_{#gamma} [deg]", 60, 0, 360},
      {"reche3_p", "He3 momentum", "p_{He3} [GeV]", 120, 0, 0.7, true},
      {"reche3_theta", "He3 polar angle", "#theta_{He3} [deg]", 60, 0, 100, true},
      {"reche3_phi", "He3 azimuth", "#phi_{He3} [deg]", 60, 0, 360, true}};
    const std::vector<Spec> dvcs = {
      {"Q2", "Q^{2}", "Q^{2} [GeV^{2}]", 120, 0, 6},
      {"xB", "x_{B}", "x_{B}", 120, 0, 0.4},
      {"WN", "Nucleon invariant mass", "W_{N} [GeV]", 120, 1, 5},
      {"WHe3", "He3 invariant mass", "W_{^{3}He} [GeV]", 120, 2.5, 9},
      {"y", "y", "y", 120, 0, 1},
      {"nu", "Energy transfer", "#nu [GeV]", 120, 0, 11},
      {"minus_t", "-t from q-#gamma", "-t [GeV^{2}]", 120, 0, 1},
      {"phi", "Trento #phi", "#phi [deg]", 120, 0, 360}};
    const std::vector<Spec> exclusivity = {
      {"mm2_ehe3", "MM^{2}(eHe3)", "MM^{2}(eHe3) [GeV^{2}]", 160, -0.3, 0.3, true},
      {"mm2_ehe3gamma", "MM^{2}(eHe3#gamma)", "MM^{2}(eHe3#gamma) [GeV^{2}]", 160, -0.2, 0.2, true},
      {"mm2_egamma", "MM^{2}(e#gamma)", "MM^{2}(e#gamma) [GeV^{2}]", 160, 6.5, 9.5},
      {"missing_energy", "Missing energy", "E_{miss} [GeV]", 160, -0.7, 0.7, true},
      {"missing_pt", "Missing transverse momentum", "P_{T,miss} [GeV]", 160, 0, 0.4, true},
      {"theta_gamma_gamma_pred", "Measured/predicted photon angle", "#theta_{#gamma#gamma'} [deg]", 160, 0, 5, true},
      {"theta_e_gamma", "Electron-photon angle", "#theta_{e#gamma} [deg]", 180, 0, 60},
      {"coplanarity", "Coplanarity", "#Delta#phi [deg]", 180, 0, 40, true}};

    fs::create_directories(outputDir);
    auto drawGroup = [&](const std::vector<Spec>& specs,
                         const std::string& tag, int columns,
                         bool appendRecQ2VsXB = false) {
      const int padCount = static_cast<int>(specs.size()) +
                           (appendRecQ2VsXB ? 1 : 0);
      const int rows = (padCount + columns - 1) / columns;
      TCanvas canvas(("c_he3dvcs_compare_" + tag).c_str(),
                     ("He3DVCS GEN/REC " + tag).c_str(),
                     500 * columns, 400 * rows);
      canvas.Divide(columns, rows);
      std::vector<std::unique_ptr<TH1D>> owned;
      owned.reserve(specs.size() * 2);
      for (size_t i = 0; i < specs.size(); ++i) {
        canvas.cd(static_cast<int>(i + 1));
        gPad->SetTicks();
        TLegend legend(0.63, 0.76, 0.88, 0.89);
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        double maximum = 0.0;
        TH1D* first = nullptr;
        bool haveRecMoments = false;
        double recMean = 0.0, recSigma = 0.0;
        const size_t modelCount = (!specs[i].requiresHe3 || drawFakeHe3Rec) ? 2 : 1;
        for (size_t model = 0; model < modelCount; ++model) {
          const auto& s = specs[i];
          const std::string name = "h_he3dvcs_compare_" + tag + "_" +
                                   std::to_string(i) + "_m" +
                                   std::to_string(model);
          // For fake-He3 REC, He3-dependent observables must use only events
          // where the simulated ALERT acceptance/PID decision succeeded.
          ROOT::RDF::RNode plotNode = *he3DVCSModels[model];
          if (model == 1 && s.requiresHe3 && drawFakeHe3Rec)
            plotNode = plotNode.Filter("fake_he3_pass");
          auto action = plotNode.Histo1D(
              {name.c_str(), (s.title + ";" + s.xTitle + ";Normalized events").c_str(),
               s.bins, s.lo, s.hi}, s.var);
          auto hist = std::unique_ptr<TH1D>(
              static_cast<TH1D*>(action.GetValue().Clone((name + "_draw").c_str())));
          hist->SetDirectory(nullptr);
          // Report the REC exclusivity peak directly from histogram moments.
          // This is not a Gaussian fit: mu=<x> and sigma=sqrt(<x^2>-<x>^2),
          // evaluated over the displayed histogram range.
          if (tag == "exclusivity" && model == 1) {
            recMean = hist->GetMean();
            recSigma = hist->GetStdDev();
            haveRecMoments = true;
            std::cout << "[He3DVCS exclusivity REC] " << s.var
                      << ": mu = " << recMean
                      << ", sigma = " << recSigma << '\n';
          }
          const double integral = hist->Integral(1, hist->GetNbinsX());
          if (integral > 0.0) hist->Scale(1.0 / integral);
          const bool scaleRecByTen = tag == "exclusivity" && model == 1 &&
                                     s.var != "theta_e_gamma";
          if (scaleRecByTen) hist->Scale(10.0);
          hist->SetLineColor(model == 0 ? kBlue + 1 : kRed + 1);
          hist->SetLineWidth(2);
          hist->SetStats(0);
          maximum = std::max(maximum, hist->GetMaximum());
          if (!first) {
            first = hist.get();
            hist->Draw("hist");
          } else {
            hist->Draw("hist same");
          }
          const std::string countLabel =
              model < he3DVCSModelEventCounts.size()
                  ? FormatMillions(he3DVCSModelEventCounts[model]) : "? M";
          const std::string modelLabel = model == 0
              ? "GEN (" + countLabel + ")"
              : (scaleRecByTen ? "10#times REC (" + countLabel + ")"
                               : "REC (" + countLabel + ")");
          legend.AddEntry(hist.get(),modelLabel.c_str(),"l");
          owned.push_back(std::move(hist));
        }
        if (first) first->SetMaximum(maximum > 0.0 ? 1.22 * maximum : 1.0);
        legend.DrawClone();
        if (haveRecMoments) {
          std::ostringstream moments;
          moments << std::setprecision(4)
                  << "#mu = " << recMean << ",  #sigma = " << recSigma;
          TLatex label;
          label.SetNDC();
          label.SetTextSize(0.036);
          label.SetTextColor(kRed + 1);
          label.DrawLatex(0.14,0.70,moments.str().c_str());
        }
      }
      std::unique_ptr<TH2D> q2VsXB;
      if (appendRecQ2VsXB) {
        canvas.cd(static_cast<int>(specs.size()+1));
        gPad->SetTicks();
        auto action=he3DVCSModels[1]->Histo2D(
            {("h_he3dvcs_"+tag+"_rec_q2_vs_xb").c_str(),
             "REC Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}];Events",
             120,0.0,0.4,120,0.0,6.0},"xB","Q2");
        q2VsXB=std::unique_ptr<TH2D>(static_cast<TH2D*>(
            action.GetValue().Clone(("h_he3dvcs_"+tag+
                                     "_rec_q2_vs_xb_draw").c_str())));
        q2VsXB->SetDirectory(nullptr);
        q2VsXB->SetStats(0);
        q2VsXB->Draw("COLZ");
        gPad->SetLogz();
        TLatex count; count.SetNDC(); count.SetTextSize(0.035);
        count.DrawLatex(0.13,0.92,
            ("GEN ("+FormatMillions(he3DVCSModelEventCounts[0])+")  REC ("+
             FormatMillions(he3DVCSModelEventCounts[1])+")").c_str());
      }
      const std::string output = outputDir + "/He3DVCS_GEN_REC_" + tag + ".pdf";
      canvas.SaveAs(output.c_str());
      std::cout << "[PlotNormalizedComparisonHe3DVCS] wrote " << output << '\n';
    };
    // One combined particle canvas: each row is e', gamma, or He3; columns
    // are normalized p/theta/phi overlays followed by REC-only p-vs-theta.
    TCanvas particleCanvas("c_he3dvcs_particles_combined",
                           "He3DVCS particle kinematics",2000,1350);
    particleCanvas.Divide(4,3);
    std::vector<std::unique_ptr<TH1D>> particleOwned1D;
    std::vector<std::unique_ptr<TH2D>> particleOwned2D;
    for (size_t row=0; row<3; ++row) {
      const bool requiresHe3 = row==2;
      for (size_t column=0; column<3; ++column) {
        const Spec& s = particles[3*row+column];
        particleCanvas.cd(static_cast<int>(4*row+column+1));
        gPad->SetTicks();
        TLegend legend(0.63,0.76,0.88,0.89);
        legend.SetBorderSize(0); legend.SetFillStyle(0);
        TH1D* first=nullptr; double maximum=0.0;
        const size_t modelCount = (!requiresHe3 || drawFakeHe3Rec) ? 2 : 1;
        for (size_t model=0; model<modelCount; ++model) {
          ROOT::RDF::RNode node=*he3DVCSModels[model];
          if (model==1 && requiresHe3) node=node.Filter("fake_he3_pass");
          const std::string name="h_he3dvcs_particle_combined_"+
              std::to_string(row)+"_"+std::to_string(column)+"_m"+
              std::to_string(model);
          auto action=node.Histo1D(
              {name.c_str(),(s.title+";"+s.xTitle+";Normalized events").c_str(),
               s.bins,s.lo,s.hi},s.var);
          auto hist=std::unique_ptr<TH1D>(static_cast<TH1D*>(
              action.GetValue().Clone((name+"_draw").c_str())));
          hist->SetDirectory(nullptr); hist->SetStats(0);
          const double integral=hist->Integral(1,hist->GetNbinsX());
          if (integral>0.0) hist->Scale(1.0/integral);
          hist->SetLineColor(model==0?kBlue+1:kRed+1);
          hist->SetLineWidth(2);
          maximum=std::max(maximum,hist->GetMaximum());
          if (!first) { first=hist.get(); hist->Draw("hist"); }
          else hist->Draw("hist same");
          const std::string modelLabel=(model==0?"GEN (":"REC (")+
              FormatMillions(he3DVCSModelEventCounts[model])+")";
          legend.AddEntry(hist.get(),modelLabel.c_str(),"l");
          particleOwned1D.push_back(std::move(hist));
        }
        if (first) first->SetMaximum(maximum>0.0?1.22*maximum:1.0);
        legend.DrawClone();
      }

      particleCanvas.cd(static_cast<int>(4*row+4));
      gPad->SetTicks();
      if (requiresHe3 && !drawFakeHe3Rec) {
        TLatex label; label.SetNDC(); label.SetTextAlign(22); label.SetTextSize(0.045);
        label.DrawLatex(0.5,0.55,"REC He3 not drawn");
        label.DrawLatex(0.5,0.47,"fake He3 is OFF");
        continue;
      }
      ROOT::RDF::RNode recNode=*he3DVCSModels[1];
      if (requiresHe3) recNode=recNode.Filter("fake_he3_pass");
      const Spec& pSpec=particles[3*row];
      const Spec& thetaSpec=particles[3*row+1];
      const std::string name="h_he3dvcs_rec_p_vs_theta_"+std::to_string(row);
      auto action=recNode.Histo2D(
          {name.c_str(),(std::string("REC ")+pSpec.title+
           ";#theta [deg];p [GeV];Events").c_str(),
           thetaSpec.bins,thetaSpec.lo,thetaSpec.hi,
           pSpec.bins,pSpec.lo,pSpec.hi},thetaSpec.var,pSpec.var);
      auto hist=std::unique_ptr<TH2D>(static_cast<TH2D*>(
          action.GetValue().Clone((name+"_draw").c_str())));
      hist->SetDirectory(nullptr); hist->SetStats(0); hist->Draw("COLZ");
      gPad->SetLogz();
      TLatex count; count.SetNDC(); count.SetTextSize(0.035);
      count.DrawLatex(0.13,0.92,
          ("GEN ("+FormatMillions(he3DVCSModelEventCounts[0])+")  REC ("+
           FormatMillions(he3DVCSModelEventCounts[1])+")").c_str());
      particleOwned2D.push_back(std::move(hist));
    }
    const std::string particleOutput=outputDir+"/He3DVCS_GEN_REC_particles.pdf";
    particleCanvas.SaveAs(particleOutput.c_str());
    std::cout << "[PlotNormalizedComparisonHe3DVCS] wrote "
              << particleOutput << '\n';

    drawGroup(dvcs, "kinematics", 4, true);
    drawGroup(exclusivity, "exclusivity", 4);
  }

  // He3DVCS acceptance: probability-weighted He3 response, direct REC/GEN
  // event acceptance, and electron acceptance inferred by factorization.
  void PlotHe3DVCSAcc() {
    if (he3DVCSModels.size()<2) {
      std::cerr<<"[PlotHe3DVCSAcc] requires GEN and REC models.\n"; return;
    }
    fs::create_directories(outputDir);
    ROOT::RDF::RNode gen=*he3DVCSModels[0];
    ROOT::RDF::RNode rec=he3DVCSModels[1]->Filter("fake_he3_pass");
    const std::string selectedCounts =
        "GEN ("+FormatMillions(he3DVCSModelEventCounts[0])+")  REC ("+
        FormatMillions(he3DVCSModelEventCounts[1])+")";
    auto drawSelectedCounts=[&]() {
      TLatex label; label.SetNDC(); label.SetTextSize(0.032);
      label.DrawLatex(0.13,0.92,selectedCounts.c_str());
    };
    auto r1=[](const TH1D& n,const TH1D& d,const std::string& name,
               const std::string& title) {
      auto h=std::unique_ptr<TH1D>(static_cast<TH1D*>(n.Clone(name.c_str())));
      h->SetDirectory(nullptr); h->Divide(&d); h->SetTitle(title.c_str());
      h->SetStats(0); h->SetMinimum(0); h->SetMaximum(1.05); h->SetLineWidth(2);
      return h;
    };
    auto r2=[](const TH2D& n,const TH2D& d,const std::string& name,
               const std::string& title) {
      auto h=std::unique_ptr<TH2D>(static_cast<TH2D*>(n.Clone(name.c_str())));
      h->SetDirectory(nullptr); h->Divide(&d); h->SetTitle(title.c_str());
      h->SetStats(0); h->SetMinimum(0); h->SetMaximum(1.0);
      h->GetZaxis()->SetTitle("Acceptance"); return h;
    };
    auto b1=[](ROOT::RDF::RNode n,const std::string& name,const std::string& v,
               int bins,double lo,double hi,const std::string& w="") {
      return w.empty()?n.Histo1D({name.c_str(),"",bins,lo,hi},v)
                      :n.Histo1D({name.c_str(),"",bins,lo,hi},v,w);
    };
    auto b2=[](ROOT::RDF::RNode n,const std::string& name,
               const std::string& x,const std::string& y,
               int nx,double xl,double xh,int ny,double yl,double yh,
               const std::string& w="") {
      return w.empty()?n.Histo2D({name.c_str(),"",nx,xl,xh,ny,yl,yh},x,y)
                      :n.Histo2D({name.c_str(),"",nx,xl,xh,ny,yl,yh},x,y,w);
    };

    auto hpD=b1(gen,"h_h3pa_pd","reche3_p",60,0,.7);
    auto htD=b1(gen,"h_h3pa_td","reche3_theta",30,0,100);
    auto hfD=b1(gen,"h_h3pa_fd","reche3_phi",36,0,360);
    auto hptD=b2(gen,"h_h3pa_ptd","reche3_theta","reche3_p",30,0,100,60,0,.7);
    auto hpW=b1(gen,"h_h3pa_pw","reche3_p",60,0,.7,"fake_he3_acceptance");
    auto htW=b1(gen,"h_h3pa_tw","reche3_theta",30,0,100,"fake_he3_acceptance");
    auto hfW=b1(gen,"h_h3pa_fw","reche3_phi",36,0,360,"fake_he3_acceptance");
    auto hptW=b2(gen,"h_h3pa_ptw","reche3_theta","reche3_p",30,0,100,60,0,.7,"fake_he3_acceptance");
    auto hp=r1(hpW.GetValue(),hpD.GetValue(),"h_h3pa_p","He3 particle acceptance;p(He3) [GeV];Acceptance");
    auto ht=r1(htW.GetValue(),htD.GetValue(),"h_h3pa_t","He3 particle acceptance;#theta(He3) [deg];Acceptance");
    auto hf=r1(hfW.GetValue(),hfD.GetValue(),"h_h3pa_f","He3 particle acceptance;#phi(He3) [deg];Acceptance");
    auto hpt=r2(hptW.GetValue(),hptD.GetValue(),"h_h3pa_pt","He3 particle acceptance;#theta(He3) [deg];p(He3) [GeV];Acceptance");
    TCanvas pc("c_h3pa","He3 particle acceptance",2000,450); pc.Divide(4,1);
    pc.cd(1);hp->Draw("HIST");drawSelectedCounts();
    pc.cd(2);ht->Draw("HIST");drawSelectedCounts();
    pc.cd(3);hf->Draw("HIST");drawSelectedCounts();
    pc.cd(4);hpt->Draw("COLZ");drawSelectedCounts();
    pc.SaveAs((outputDir+"/He3DVCS_acceptance_particle_He3.pdf").c_str());

    auto hpR=b1(rec,"h_h3ea_pr","reche3_p",60,0,.7);
    auto htR=b1(rec,"h_h3ea_tr","reche3_theta",30,0,100);
    auto hfR=b1(rec,"h_h3ea_fr","reche3_phi",36,0,360);
    auto hptR=b2(rec,"h_h3ea_ptr","reche3_theta","reche3_p",30,0,100,60,0,.7);
    auto evp=r1(hpR.GetValue(),hpD.GetValue(),"h_h3ea_p","He3DVCS event acceptance;p(He3) [GeV];Event acceptance");
    auto evt=r1(htR.GetValue(),htD.GetValue(),"h_h3ea_t","He3DVCS event acceptance;#theta(He3) [deg];Event acceptance");
    auto evf=r1(hfR.GetValue(),hfD.GetValue(),"h_h3ea_f","He3DVCS event acceptance;#phi(He3) [deg];Event acceptance");
    auto evpt=r2(hptR.GetValue(),hptD.GetValue(),"h_h3ea_pt","He3DVCS event acceptance;#theta(He3) [deg];p(He3) [GeV];Event acceptance");

    auto epD=b1(gen,"h_h3el_pd","recel_p",60,0,11);
    auto etD=b1(gen,"h_h3el_td","recel_theta",30,0,30);
    auto efD=b1(gen,"h_h3el_fd","recel_phi",36,0,360);
    auto eptD=b2(gen,"h_h3el_ptd","recel_theta","recel_p",30,0,30,60,0,11);
    auto epR=b1(rec,"h_h3el_pr","recel_p",60,0,11);
    auto etR=b1(rec,"h_h3el_tr","recel_theta",30,0,30);
    auto efR=b1(rec,"h_h3el_fr","recel_phi",36,0,360);
    auto eptR=b2(rec,"h_h3el_ptr","recel_theta","recel_p",30,0,30,60,0,11);
    auto epS=b1(gen,"h_h3el_ps","recel_p",60,0,11,"fake_he3_acceptance");
    auto etS=b1(gen,"h_h3el_ts","recel_theta",30,0,30,"fake_he3_acceptance");
    auto efS=b1(gen,"h_h3el_fs","recel_phi",36,0,360,"fake_he3_acceptance");
    auto eptS=b2(gen,"h_h3el_pts","recel_theta","recel_p",30,0,30,60,0,11,"fake_he3_acceptance");
    auto tep=r1(epR.GetValue(),epD.GetValue(),"h_h3el_tp","");
    auto tet=r1(etR.GetValue(),etD.GetValue(),"h_h3el_tt","");
    auto tef=r1(efR.GetValue(),efD.GetValue(),"h_h3el_tf","");
    auto tept=r2(eptR.GetValue(),eptD.GetValue(),"h_h3el_tpt","");
    auto sep=r1(epS.GetValue(),epD.GetValue(),"h_h3el_sp","");
    auto set=r1(etS.GetValue(),etD.GetValue(),"h_h3el_st","");
    auto sef=r1(efS.GetValue(),efD.GetValue(),"h_h3el_sf","");
    auto sept=r2(eptS.GetValue(),eptD.GetValue(),"h_h3el_spt","");
    auto eap=r1(*tep,*sep,"h_h3el_p","Inferred electron acceptance;p(e') [GeV];Acceptance");
    auto eat=r1(*tet,*set,"h_h3el_t","Inferred electron acceptance;#theta(e') [deg];Acceptance");
    auto eaf=r1(*tef,*sef,"h_h3el_f","Inferred electron acceptance;#phi(e') [deg];Acceptance");
    auto eapt=r2(*tept,*sept,"h_h3el_pt","Inferred electron acceptance;#theta(e') [deg];p(e') [GeV];Acceptance");
    TCanvas ec("c_h3ea","He3DVCS event and electron acceptance",2000,900); ec.Divide(4,2);
    ec.cd(1);evp->Draw("HIST");drawSelectedCounts();
    ec.cd(2);evt->Draw("HIST");drawSelectedCounts();
    ec.cd(3);evf->Draw("HIST");drawSelectedCounts();
    ec.cd(4);evpt->Draw("COLZ");drawSelectedCounts();
    ec.cd(5);eap->Draw("HIST");drawSelectedCounts();
    ec.cd(6);eat->Draw("HIST");drawSelectedCounts();
    ec.cd(7);eaf->Draw("HIST");drawSelectedCounts();
    ec.cd(8);eapt->Draw("COLZ");drawSelectedCounts();
    ec.SaveAs((outputDir+"/He3DVCS_acceptance_event.pdf").c_str());

    struct KV {std::string v,label;int bins;double lo,hi;};
    const std::vector<KV> kv={{"xB","x_{B}",60,0,.4},{"Q2","Q^{2} [GeV^{2}]",60,0,6},
      {"W","W [GeV]",60,1,5},{"y","y",60,0,1},{"nu","#nu [GeV]",60,0,11}};
    std::vector<std::unique_ptr<TH1D>> ka;
    for(size_t i=0;i<kv.size();++i){
      auto d=b1(gen,"h_h3ka_d"+std::to_string(i),kv[i].v,kv[i].bins,kv[i].lo,kv[i].hi);
      auto n=b1(rec,"h_h3ka_n"+std::to_string(i),kv[i].v,kv[i].bins,kv[i].lo,kv[i].hi);
      ka.push_back(r1(n.GetValue(),d.GetValue(),"h_h3ka_a"+std::to_string(i),
          "He3DVCS event acceptance;"+kv[i].label+";Event acceptance"));
    }
    auto xqd=b2(gen,"h_h3ka_xqd","xB","Q2",60,0,.4,60,0,6);
    auto xqn=b2(rec,"h_h3ka_xqn","xB","Q2",60,0,.4,60,0,6);
    auto xqa=r2(xqn.GetValue(),xqd.GetValue(),"h_h3ka_xqa",
        "He3DVCS event acceptance;x_{B};Q^{2} [GeV^{2}];Event acceptance");
    xqa->GetZaxis()->SetTitle("Event acceptance");
    TCanvas kc("c_h3ka","He3DVCS event acceptance",1800,900);kc.Divide(3,2);
    for(size_t i=0;i<ka.size();++i){
      kc.cd(int(i+1));ka[i]->Draw("HIST");drawSelectedCounts();
    }
    kc.cd(6);xqa->Draw("COLZ");drawSelectedCounts();
    kc.SaveAs((outputDir+"/He3DVCS_acceptance_event_DIS.pdf").c_str());
    std::cout<<"[PlotHe3DVCSAcc] wrote He3DVCS acceptance PDFs\n";
  }

  // Plot coherent-He3 longitudinal spin asymmetries versus Trento phi.
  // The four yields are N_{beam,target}; generator events are already
  // unweighted, so event_weight is one for both GEN and REC.  Both models use
  // identical spin definitions and binning and are overlaid in every panel.
  void PlotAsymmetryHe3(double beamPolarization=0.80,
                        double targetPolarization=0.50,
                        bool drawCorrected=true) {
    if (he3DVCSModels.size() < 2) {
      std::cerr << "[PlotAsymmetryHe3] requires GEN and REC models.\n";
      return;
    }
    if (!std::isfinite(beamPolarization) || !std::isfinite(targetPolarization) ||
        beamPolarization <= 0.0 || beamPolarization > 1.0 ||
        targetPolarization <= 0.0 || targetPolarization > 1.0)
      throw std::runtime_error("He3 beam/target polarizations must be in (0,1]");

    struct SpinYields {
      std::unique_ptr<TH1D> pp, pm, mp, mm;
    };
    struct AsymmetrySpec {
      const char* key;
      const char* yTitle;
      // First pair minus second pair, divided by their sum.
      std::array<int, 4> signs;
    };
    const std::array<AsymmetrySpec, 3> specs = {{
        {"ALU", "A_{LU}", {{+1, +1, -1, -1}}},
        {"AUL", "A_{UL}", {{+1, -1, +1, -1}}},
        {"ALL", "A_{LL}", {{+1, -1, -1, +1}}}
    }};

    fs::create_directories(outputDir);
    std::vector<SpinYields> yields;
    yields.reserve(2);
    for (size_t model = 0; model < 2; ++model) {
      SpinYields spin;
      const std::array<std::pair<const char*, std::unique_ptr<TH1D>*>, 4> channels = {{
          {"real_beam_spin > 0 && real_target_spin > 0", &spin.pp},
          {"real_beam_spin > 0 && real_target_spin < 0", &spin.pm},
          {"real_beam_spin < 0 && real_target_spin > 0", &spin.mp},
          {"real_beam_spin < 0 && real_target_spin < 0", &spin.mm}
      }};
      for (size_t channel = 0; channel < channels.size(); ++channel) {
        const std::string name = "h_he3_asym_yield_m" + std::to_string(model) +
                                 "_s" + std::to_string(channel);
        auto action = he3DVCSModels[model]->Filter(channels[channel].first)
            .Histo1D({name.c_str(), ";#phi [deg];Weighted yield",
                      24, 0.0, 360.0}, "phi", "event_weight");
        *channels[channel].second = std::unique_ptr<TH1D>(
            static_cast<TH1D*>(action.GetValue().Clone((name + "_draw").c_str())));
        (*channels[channel].second)->SetDirectory(nullptr);
      }
      yields.push_back(std::move(spin));
    }

    auto makeAsymmetry = [drawCorrected](const std::string& name,
                            const AsymmetrySpec& spec,const SpinYields& spin,
                            double polarizationFactor) {
      auto result = std::make_unique<TH1D>(
          name.c_str(), (std::string(";#phi [deg];") + spec.yTitle).c_str(),
          24, 0.0, 360.0);
      result->SetDirectory(nullptr);
      result->SetStats(0);
      const std::array<const TH1D*, 4> h = {{
          spin.pp.get(), spin.pm.get(), spin.mp.get(), spin.mm.get()}};
      for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
        double numerator = 0.0, denominator = 0.0, varianceNumerator = 0.0;
        for (size_t i = 0; i < h.size(); ++i) {
          const double n = h[i]->GetBinContent(bin);
          numerator += spec.signs[i] * n;
          denominator += n;
        }
        if (denominator <= 0.0) continue;
        const double asymmetry = numerator / denominator;
        // General propagation for A=(sum s_i N_i)/(sum N_i), including Sumw2.
        for (size_t i = 0; i < h.size(); ++i) {
          const double derivative = (spec.signs[i] - asymmetry) / denominator;
          varianceNumerator += derivative * derivative *
                               std::pow(h[i]->GetBinError(bin), 2);
        }
        const double correction=drawCorrected?polarizationFactor:1.0;
        result->SetBinContent(bin,asymmetry/correction);
        result->SetBinError(bin,
            std::sqrt(std::max(0.0,varianceNumerator))/correction);
      }
      return result;
    };

    TCanvas canvas("c_he3_asymmetry", "He3DVCS spin asymmetries", 1800, 600);
    canvas.Divide(3, 1);
    std::vector<std::unique_ptr<TH1D>> owned;
    owned.reserve(6);
    for (size_t observable = 0; observable < specs.size(); ++observable) {
      canvas.cd(static_cast<int>(observable + 1));
      gPad->SetTicks();
      TLegend legend(0.68, 0.78, 0.88, 0.90);
      legend.SetBorderSize(0);
      legend.SetFillStyle(0);
      const double polarizationFactor=observable==0?beamPolarization:
          (observable==1?targetPolarization:beamPolarization*targetPolarization);
      for (size_t model = 0; model < 2; ++model) {
        auto hist = makeAsymmetry(
            "h_he3_" + std::string(specs[observable].key) + "_m" +
                std::to_string(model), specs[observable], yields[model],
                polarizationFactor);
        const double yRange = observable == 0 ? 1.0 : 0.1;
        hist->SetMinimum(-yRange);
        hist->SetMaximum(+yRange);
        hist->SetLineColor(model == 0 ? kBlue + 1 : kRed + 1);
        hist->SetMarkerColor(model == 0 ? kBlue + 1 : kRed + 1);
        hist->SetMarkerStyle(model == 0 ? 20 : 24);
        hist->SetMarkerSize(0.85);
        hist->SetLineWidth(2);
        hist->Draw(model == 0 ? "E1" : "E1 SAME");
        const std::string modelLabel=(model==0?"GEN (":"REC (")+
            FormatMillions(he3DVCSModelEventCounts[model])+")";
        legend.AddEntry(hist.get(),modelLabel.c_str(),"lep");
        owned.push_back(std::move(hist));
      }
      TLine zero(0.0, 0.0, 360.0, 0.0);
      zero.SetLineStyle(2);
      zero.SetLineColor(kGray + 1);
      zero.DrawClone();
      legend.DrawClone();
      std::ostringstream polarizationLabel;
      polarizationLabel<<std::fixed<<std::setprecision(2)
          <<(drawCorrected?"Reconstructed A":"Raw A")
          <<", P_{b}="<<beamPolarization<<", P_{t}="<<targetPolarization;
      TLatex label;label.SetNDC();label.SetTextSize(0.035);
      label.DrawLatex(0.13,0.92,polarizationLabel.str().c_str());
    }
    const std::string output = outputDir + "/He3DVCS_GEN_REC_asymmetries.pdf";
    canvas.SaveAs(output.c_str());
    std::cout << "[PlotAsymmetryHe3] Pb=" << beamPolarization
              << ", Pt=" << targetPolarization
              << ", output=" << (drawCorrected?"reconstructed A":"raw A")
              << "\n[PlotAsymmetryHe3] wrote " << output << '\n';
  }

  void AddModelHe3DIS(ROOT::RDF::RNode all, ROOT::RDF::RNode epp,
                      ROOT::RDF::RNode epn, ROOT::RDF::RNode ed,
                      const std::string& label) {
    he3DISModels.push_back({
        std::make_unique<ROOT::RDF::RNode>(std::move(all)),
        std::make_unique<ROOT::RDF::RNode>(std::move(epp)),
        std::make_unique<ROOT::RDF::RNode>(std::move(epn)),
        std::make_unique<ROOT::RDF::RNode>(std::move(ed)), label});
    std::cout << "Adding He3DIS model: " << label
              << ". Total He3DIS models: " << he3DISModels.size() << '\n';
  }

  // Overlay unit-area GEN and REC distributions for inclusive He3DIS and
  // the epp, epn, and ed topology views.
  void PlotNormalizedComparisonHe3DIS() {
    if (he3DISModels.size() < 2) {
      std::cerr << "[PlotNormalizedComparisonHe3DIS] requires GEN and REC.\n";
      return;
    }
    struct Spec {
      std::string var, title, xTitle;
      int bins;
      double lo, hi;
    };
    struct PairSpec {
      std::string p, theta, title;
      double pHi, thetaHi;
    };
    const std::vector<Spec> electron = {
        {"recel_p", "Electron momentum", "p_{e'} [GeV]", 120, 0, 11},
        {"recel_theta", "Electron polar angle", "#theta_{e'} [deg]", 60, 0, 35},
        {"recel_phi", "Electron azimuth", "#phi_{e'} [deg]", 60, 0, 360}};
    const std::vector<Spec> dis = {
        {"Q2", "Q^{2}", "Q^{2} [GeV^{2}]", 120, 0, 10},
        {"xB", "x_{B}", "x_{B}", 120, 0, 1},
        {"WN", "Nucleon invariant mass", "W_{N} [GeV]", 120, 0, 5},
        {"WHe3", "He3 invariant mass", "W_{^{3}He} [GeV]", 120, 2.5, 9},
        {"nu", "nu", "#nu [GeV]", 120, 0, 11},
        {"y", "y", "y", 120, 0, 1}};
    const std::vector<Spec> protons = {
        {"pro1_p", "Leading proton momentum", "p_{p1} [GeV]", 120, 0, 0.5},
        {"pro1_theta", "Leading proton angle", "#theta_{p1} [deg]", 60, 0, 180},
        {"pro1_phi", "Leading proton azimuth", "#phi_{p1} [deg]", 60, 0, 360},
        {"pro2_p", "Subleading proton momentum", "p_{p2} [GeV]", 120, 0, 0.5},
        {"pro2_theta", "Subleading proton angle", "#theta_{p2} [deg]", 60, 0, 180},
        {"pro2_phi", "Subleading proton azimuth", "#phi_{p2} [deg]", 60, 0, 360}};
    const std::vector<Spec> neutron = {
        {"neu_p", "Neutron momentum", "p_{n} [GeV]", 120, 0, 3},
        {"neu_theta", "Neutron angle", "#theta_{n} [deg]", 60, 0, 180},
        {"neu_phi", "Neutron azimuth", "#phi_{n} [deg]", 60, 0, 360}};
    const std::vector<Spec> deuteron = {
        {"deu_p", "Deuteron momentum", "p_{d} [GeV]", 120, 0, 0.5},
        {"deu_theta", "Deuteron angle", "#theta_{d} [deg]", 60, 0, 180},
        {"deu_phi", "Deuteron azimuth", "#phi_{d} [deg]", 60, 0, 360}};

    auto combine = [](std::initializer_list<std::vector<Spec>> groups) {
      std::vector<Spec> result;
      for (const auto& group : groups)
        result.insert(result.end(), group.begin(), group.end());
      return result;
    };
    fs::create_directories(outputDir);
    auto nodeFor = [](He3DISModel& model, const std::string& topology)
        -> ROOT::RDF::RNode& {
      if (topology == "epp") return *model.epp;
      if (topology == "epn") return *model.epn;
      if (topology == "ed") return *model.ed;
      return *model.all;
    };
    auto draw = [&](const std::string& category,
                    const std::string& topology,
                    const std::vector<Spec>& specs,
                    const std::vector<PairSpec>& pairs) {
      const std::string tag = category + "_" + topology;
      const int pads = static_cast<int>(specs.size() + pairs.size());
      const int columns = 4;
      const int rows = (pads + columns - 1)/columns;
      // Particle specifications are stored as repeating p/theta/phi triples.
      // When each triple has a matching 2D entry, reserve the fourth pad in
      // that same row for p-vs-theta.
      const bool particleRows = !pairs.empty() && specs.size() == 3*pairs.size();
      TCanvas canvas(("c_he3dis_" + tag).c_str(),
                     ("He3DIS GEN/REC " + tag).c_str(),
                     500*columns, 400*rows);
      canvas.Divide(columns, rows);
      std::vector<std::unique_ptr<TH1D>> owned1D;
      std::vector<std::unique_ptr<TH2D>> owned2D;
      for (size_t i = 0; i < specs.size(); ++i) {
        const size_t pad = particleRows ? 4*(i/3)+(i%3)+1 : i+1;
        canvas.cd(static_cast<int>(pad));
        gPad->SetTicks();
        TLegend legend(0.66, 0.77, 0.88, 0.90);
        legend.SetBorderSize(0); legend.SetFillStyle(0);
        TH1D* first = nullptr;
        double maximum = 0.0;
        for (size_t model = 0; model < 2; ++model) {
          const auto& spec = specs[i];
          const std::string name = "h_he3dis_" + tag + "_" +
              std::to_string(i) + "_m" + std::to_string(model);
          auto action = nodeFor(he3DISModels[model], topology).Histo1D(
              {name.c_str(), (spec.title + ";" + spec.xTitle +
               ";Normalized events").c_str(), spec.bins, spec.lo, spec.hi},
              spec.var);
          auto hist = std::unique_ptr<TH1D>(static_cast<TH1D*>(
              action.GetValue().Clone((name + "_draw").c_str())));
          hist->SetDirectory(nullptr); hist->SetStats(0);
          const double integral = hist->Integral(1, hist->GetNbinsX());
          if (integral > 0.0) hist->Scale(1.0/integral);
          hist->SetLineColor(model == 0 ? kBlue+1 : kRed+1);
          hist->SetLineWidth(2);
          maximum = std::max(maximum, hist->GetMaximum());
          if (!first) { first = hist.get(); hist->Draw("HIST"); }
          else hist->Draw("HIST SAME");
          const std::string label = he3DISModels[model].label +
              (hist->GetEntries() == 0 ? " (no entries)" : "");
          legend.AddEntry(hist.get(), label.c_str(), "l");
          owned1D.push_back(std::move(hist));
        }
        if (first) first->SetMaximum(maximum > 0.0 ? 1.22*maximum : 1.0);
        legend.DrawClone();
      }
      for (size_t i = 0; i < pairs.size(); ++i) {
        const size_t pad = particleRows ? 4*i+4 : specs.size()+i+1;
        canvas.cd(static_cast<int>(pad));
        gPad->SetTicks();
        const auto& pair = pairs[i];
        const std::string name = "h_he3dis_" + tag + "_rec_2d_" +
                                 std::to_string(i);
        auto action = nodeFor(he3DISModels[1], topology).Histo2D(
            {name.c_str(), ("REC " + pair.title +
             ";#theta [deg];p [GeV];Events").c_str(),
             50, 0, pair.thetaHi, 100, 0, pair.pHi},
            pair.theta, pair.p);
        auto hist = std::unique_ptr<TH2D>(static_cast<TH2D*>(
            action.GetValue().Clone((name + "_draw").c_str())));
        hist->SetDirectory(nullptr); hist->SetStats(0);
        hist->Draw("COLZ"); gPad->SetLogz();
        owned2D.push_back(std::move(hist));
      }
      const std::string output =
          outputDir + "/He3DIS_GEN_REC_" + tag + ".pdf";
      canvas.SaveAs(output.c_str());
      std::cout << "[PlotNormalizedComparisonHe3DIS] wrote "
                << output << '\n';
    };

    // Every general plot class is produced for every topology.  Keeping the
    // category and topology as separate arguments is important: the latter
    // selects the correct RDF, while both appear in the output file name.
    for (const auto& topology : {std::string("all"), std::string("epp"),
                                 std::string("epn"), std::string("ed")}) {
      auto specs = dis;
      if (topology != "all") {
        for (const auto& s : He3DISSpectatorPlotSpecs())
          specs.push_back({topology+"_"+s.variable,
              topology+" spectator "+s.variable, s.axis, s.bins, s.lo, s.hi});
      }
      draw("kinematics", topology, specs, {});
    }
    const auto allSpectators = combine({protons, neutron, deuteron});
    const std::vector<PairSpec> allParticlePairs = {
        {"recel_p", "recel_theta", "electron p vs theta", 11, 35},
        {"pro1_p", "pro1_theta", "leading proton p vs theta", 0.5, 180},
        {"pro2_p", "pro2_theta", "subleading proton p vs theta", 0.5, 180},
        {"neu_p", "neu_theta", "neutron p vs theta", 3, 180},
        {"deu_p", "deu_theta", "deuteron p vs theta", 0.5, 180}};
    draw("particles", "all", combine({electron,allSpectators}),
         allParticlePairs);
    draw("particles", "epp", combine({electron,protons}),
         {{"recel_p", "recel_theta", "electron p vs theta", 11, 35},
          {"pro1_p", "pro1_theta", "leading proton p vs theta", 0.5, 180},
          {"pro2_p", "pro2_theta", "subleading proton p vs theta", 0.5, 180}});
    draw("particles", "epn",
         combine({electron,
                  std::vector<Spec>(protons.begin(),protons.begin()+3),neutron}),
         {{"recel_p", "recel_theta", "electron p vs theta", 11, 35},
          {"pro1_p", "pro1_theta", "proton p vs theta", 0.5, 180},
          {"neu_p", "neu_theta", "neutron p vs theta", 3, 180}});
    draw("particles", "ed", combine({electron,deuteron}),
         {{"recel_p", "recel_theta", "electron p vs theta", 11, 35},
          {"deu_p", "deu_theta", "deuteron p vs theta", 0.5, 180}});
  }

  // Particle acceptance uses the probability-weighted fake ALERT response.
  // Event acceptance is counted directly as REC/GEN.  Their ratio is shown as
  // an inferred electron acceptance under the factorization Aevent=Ae*Aspect.
  void PlotHe3DISAcc(ROOT::RDF::RNode gen, ROOT::RDF::RNode rec) {
    fs::create_directories(outputDir);
    struct AccSpec {
      std::string tag, title, p, theta, phi, weight, selection;
    };
    const std::vector<AccSpec> particleSpecs = {
        {"particle_p", "Per-particle proton", "acc_p_p", "acc_p_theta",
         "acc_p_phi", "acc_p_weight", "GEN_pass_epp || GEN_pass_epn"},
        {"particle_d", "Per-particle deuteron", "acc_d_p", "acc_d_theta",
         "acc_d_phi", "acc_d_weight", "GEN_pass_ed"}};

    auto ratio1D=[](const TH1D& numerator,const TH1D& denominator,
                    const std::string& name,const std::string& title) {
      auto result=std::unique_ptr<TH1D>(static_cast<TH1D*>(
          numerator.Clone(name.c_str())));
      result->SetDirectory(nullptr); result->Divide(&denominator);
      result->SetTitle(title.c_str()); result->SetStats(0);
      result->SetMinimum(0.0); result->SetMaximum(1.05);
      result->SetLineWidth(2); return result;
    };
    auto ratio2D=[](const TH2D& numerator,const TH2D& denominator,
                    const std::string& name,const std::string& title) {
      auto result=std::unique_ptr<TH2D>(static_cast<TH2D*>(
          numerator.Clone(name.c_str())));
      result->SetDirectory(nullptr); result->Divide(&denominator);
      result->SetTitle(title.c_str()); result->SetStats(0);
      result->SetMinimum(0.0); result->SetMaximum(1.0);
      result->GetZaxis()->SetTitle("Acceptance"); return result;
    };

    auto drawOne = [&](const AccSpec& spec) {
      auto node=gen.Filter(spec.selection);
      const double pHi=0.5;
      const std::string prefix="h_he3dis_acc_"+spec.tag;
      ROOT::RDF::TH1DModel pModel((prefix+"_p_den").c_str(),
          (spec.title+";p [GeV];Acceptance").c_str(),50,0,pHi);
      ROOT::RDF::TH1DModel thetaModel((prefix+"_theta_den").c_str(),
          (spec.title+";#theta [deg];Acceptance").c_str(),30,0,180);
      ROOT::RDF::TH1DModel phiModel((prefix+"_phi_den").c_str(),
          (spec.title+";#phi [deg];Acceptance").c_str(),36,0,360);
      ROOT::RDF::TH2DModel ptModel((prefix+"_pt_den").c_str(),
          (spec.title+";#theta [deg];p [GeV]").c_str(),30,0,180,50,0,pHi);

      auto pDen=node.Histo1D(pModel,spec.p);
      auto pNum=node.Histo1D(
          {(prefix+"_p_num").c_str(),"",50,0,pHi},spec.p,spec.weight);
      auto thetaDen=node.Histo1D(thetaModel,spec.theta);
      auto thetaNum=node.Histo1D(
          {(prefix+"_theta_num").c_str(),"",30,0,180},spec.theta,spec.weight);
      auto phiDen=node.Histo1D(phiModel,spec.phi);
      auto phiNum=node.Histo1D(
          {(prefix+"_phi_num").c_str(),"",36,0,360},spec.phi,spec.weight);
      auto ptDen=node.Histo2D(ptModel,spec.theta,spec.p);
      auto ptNum=node.Histo2D(
          {(prefix+"_pt_num").c_str(),"",30,0,180,50,0,pHi},
          spec.theta,spec.p,spec.weight);

      auto rp=ratio1D(pNum.GetValue(),pDen.GetValue(),prefix+"_p_acc",
          spec.title+" acceptance;p [GeV];Acceptance");
      auto rt=ratio1D(thetaNum.GetValue(),thetaDen.GetValue(),prefix+"_theta_acc",
          spec.title+" acceptance;#theta [deg];Acceptance");
      auto rf=ratio1D(phiNum.GetValue(),phiDen.GetValue(),prefix+"_phi_acc",
          spec.title+" acceptance;#phi [deg];Acceptance");
      auto r2=ratio2D(ptNum.GetValue(),ptDen.GetValue(),prefix+"_pt_acc",
          spec.title+" acceptance;#theta [deg];p [GeV];Acceptance");

      TCanvas canvas(("c_"+prefix).c_str(),spec.title.c_str(),2000,450);
      canvas.Divide(4,1);
      canvas.cd(1); rp->Draw("HIST");
      canvas.cd(2); rt->Draw("HIST");
      canvas.cd(3); rf->Draw("HIST");
      canvas.cd(4); r2->Draw("COLZ");
      const std::string output=outputDir+"/He3DIS_acceptance_"+spec.tag+".pdf";
      canvas.SaveAs(output.c_str());
      std::cout<<"[PlotHe3DISAcc] wrote "<<output<<'\n';
    };
    for (const auto& spec:particleSpecs) drawOne(spec);

    struct EventSpec {
      std::string tag,title,genP,genTheta,genPhi,recP,recTheta,recPhi;
      std::string spectatorWeight,genSelection,recSelection,axisParticle;
    };
    const std::vector<EventSpec> eventSpecs = {
      {"event_epp","epp","acc_epp_p","acc_epp_theta","acc_epp_phi",
       "pro1_p","pro1_theta","pro1_phi","acc_event_epp","GEN_pass_epp",
       "REC_Event_pass_epp","leading proton"},
      {"event_ed","ed","acc_ed_p","acc_ed_theta","acc_ed_phi",
       "deu_p","deu_theta","deu_phi","acc_event_ed","GEN_pass_ed",
       "REC_Event_pass_ed","deuteron"}};
    for (const auto& spec:eventSpecs) {
      auto gd=gen.Filter(spec.genSelection);
      auto rn=rec.Filter(spec.recSelection);
      const double pHi=0.5;
      const std::string prefix="h_he3dis_acc_"+spec.tag;
      auto pDen=gd.Histo1D({(prefix+"_p_den").c_str(),"",50,0,pHi},spec.genP);
      auto tDen=gd.Histo1D({(prefix+"_t_den").c_str(),"",30,0,180},spec.genTheta);
      auto fDen=gd.Histo1D({(prefix+"_f_den").c_str(),"",36,0,360},spec.genPhi);
      auto ptDen=gd.Histo2D({(prefix+"_pt_den").c_str(),"",30,0,180,50,0,pHi},spec.genTheta,spec.genP);
      auto pRec=rn.Histo1D({(prefix+"_p_rec").c_str(),"",50,0,pHi},spec.recP);
      auto tRec=rn.Histo1D({(prefix+"_t_rec").c_str(),"",30,0,180},spec.recTheta);
      auto fRec=rn.Histo1D({(prefix+"_f_rec").c_str(),"",36,0,360},spec.recPhi);
      auto ptRec=rn.Histo2D({(prefix+"_pt_rec").c_str(),"",30,0,180,50,0,pHi},spec.recTheta,spec.recP);
      // The inferred electron acceptance must be binned in electron variables,
      // not in spectator variables.  GEN recel_* columns are MC electron truth.
      auto epDen=gd.Histo1D({(prefix+"_ep_den").c_str(),"",60,0,11},"recel_p");
      auto etDen=gd.Histo1D({(prefix+"_et_den").c_str(),"",30,0,35},"recel_theta");
      auto efDen=gd.Histo1D({(prefix+"_ef_den").c_str(),"",36,0,360},"recel_phi");
      auto eptDen=gd.Histo2D({(prefix+"_ept_den").c_str(),"",30,0,35,60,0,11},"recel_theta","recel_p");
      auto epRec=rn.Histo1D({(prefix+"_ep_rec").c_str(),"",60,0,11},"recel_p");
      auto etRec=rn.Histo1D({(prefix+"_et_rec").c_str(),"",30,0,35},"recel_theta");
      auto efRec=rn.Histo1D({(prefix+"_ef_rec").c_str(),"",36,0,360},"recel_phi");
      auto eptRec=rn.Histo2D({(prefix+"_ept_rec").c_str(),"",30,0,35,60,0,11},"recel_theta","recel_p");
      auto epSpec=gd.Histo1D({(prefix+"_ep_spec").c_str(),"",60,0,11},"recel_p",spec.spectatorWeight);
      auto etSpec=gd.Histo1D({(prefix+"_et_spec").c_str(),"",30,0,35},"recel_theta",spec.spectatorWeight);
      auto efSpec=gd.Histo1D({(prefix+"_ef_spec").c_str(),"",36,0,360},"recel_phi",spec.spectatorWeight);
      auto eptSpec=gd.Histo2D({(prefix+"_ept_spec").c_str(),"",30,0,35,60,0,11},"recel_theta","recel_p",spec.spectatorWeight);

      const std::string base="Total "+spec.title+" event acceptance";
      auto totalP=ratio1D(pRec.GetValue(),pDen.GetValue(),prefix+"_total_p",
          base+";p("+spec.axisParticle+") [GeV];Acceptance");
      auto totalT=ratio1D(tRec.GetValue(),tDen.GetValue(),prefix+"_total_t",
          base+";#theta("+spec.axisParticle+") [deg];Acceptance");
      auto totalF=ratio1D(fRec.GetValue(),fDen.GetValue(),prefix+"_total_f",
          base+";#phi("+spec.axisParticle+") [deg];Acceptance");
      auto totalPT=ratio2D(ptRec.GetValue(),ptDen.GetValue(),prefix+"_total_pt",
          base+";#theta("+spec.axisParticle+") [deg];p("+spec.axisParticle+") [GeV];Acceptance");
      auto totalEP=ratio1D(epRec.GetValue(),epDen.GetValue(),prefix+"_total_ep","");
      auto totalET=ratio1D(etRec.GetValue(),etDen.GetValue(),prefix+"_total_et","");
      auto totalEF=ratio1D(efRec.GetValue(),efDen.GetValue(),prefix+"_total_ef","");
      auto totalEPT=ratio2D(eptRec.GetValue(),eptDen.GetValue(),prefix+"_total_ept","");
      auto specEP=ratio1D(epSpec.GetValue(),epDen.GetValue(),prefix+"_spec_ep","");
      auto specET=ratio1D(etSpec.GetValue(),etDen.GetValue(),prefix+"_spec_et","");
      auto specEF=ratio1D(efSpec.GetValue(),efDen.GetValue(),prefix+"_spec_ef","");
      auto specEPT=ratio2D(eptSpec.GetValue(),eptDen.GetValue(),prefix+"_spec_ept","");
      auto eleP=ratio1D(*totalEP,*specEP,prefix+"_electron_p",
          "Inferred electron acceptance from "+spec.title+";p(e') [GeV];Acceptance");
      auto eleT=ratio1D(*totalET,*specET,prefix+"_electron_t",
          "Inferred electron acceptance from "+spec.title+";#theta(e') [deg];Acceptance");
      auto eleF=ratio1D(*totalEF,*specEF,prefix+"_electron_f",
          "Inferred electron acceptance from "+spec.title+";#phi(e') [deg];Acceptance");
      auto elePT=ratio2D(*totalEPT,*specEPT,prefix+"_electron_pt",
          "Inferred electron acceptance from "+spec.title+";#theta(e') [deg];p(e') [GeV];Acceptance");

      TCanvas canvas(("c_"+prefix).c_str(),spec.title.c_str(),2000,900);
      canvas.Divide(4,2);
      canvas.cd(1); totalP->Draw("HIST"); canvas.cd(2); totalT->Draw("HIST");
      canvas.cd(3); totalF->Draw("HIST"); canvas.cd(4); totalPT->Draw("COLZ");
      canvas.cd(5); eleP->Draw("HIST"); canvas.cd(6); eleT->Draw("HIST");
      canvas.cd(7); eleF->Draw("HIST"); canvas.cd(8); elePT->Draw("COLZ");
      const std::string output=outputDir+"/He3DIS_acceptance_"+spec.tag+".pdf";
      canvas.SaveAs(output.c_str());
      std::cout<<"[PlotHe3DISAcc] wrote "<<output<<'\n';

    }
    // DIS and spectator-system acceptance for every topology, including epn
    // (zero REC by construction while neutron detection is disabled).
    for (const std::string topology : {"epp", "epn", "ed"}) {
      auto gd=gen.Filter("GEN_pass_"+topology);
      auto rn=rec.Filter("REC_Event_pass_"+topology);
      const std::string prefix="h_he3dis_acc_event_"+topology;
      // Direct REC/GEN event acceptance in the inclusive DIS variables.
      auto xbDen=gd.Histo1D({(prefix+"_xb_den").c_str(),"",60,0,1},"xB");
      auto q2Den=gd.Histo1D({(prefix+"_q2_den").c_str(),"",60,0,10},"Q2");
      auto wDen=gd.Histo1D({(prefix+"_w_den").c_str(),"",60,0,5},"W");
      auto yDen=gd.Histo1D({(prefix+"_y_den").c_str(),"",60,0,1},"y");
      auto nuDen=gd.Histo1D({(prefix+"_nu_den").c_str(),"",60,0,11},"nu");
      auto xbq2Den=gd.Histo2D({(prefix+"_xbq2_den").c_str(),"",60,0,1,60,0,10},"xB","Q2");
      auto xbRec=rn.Histo1D({(prefix+"_xb_rec").c_str(),"",60,0,1},"xB");
      auto q2Rec=rn.Histo1D({(prefix+"_q2_rec").c_str(),"",60,0,10},"Q2");
      auto wRec=rn.Histo1D({(prefix+"_w_rec").c_str(),"",60,0,5},"W");
      auto yRec=rn.Histo1D({(prefix+"_y_rec").c_str(),"",60,0,1},"y");
      auto nuRec=rn.Histo1D({(prefix+"_nu_rec").c_str(),"",60,0,11},"nu");
      auto xbq2Rec=rn.Histo2D({(prefix+"_xbq2_rec").c_str(),"",60,0,1,60,0,10},"xB","Q2");
      const std::string disTitle=topology+" event acceptance";
      auto accXB=ratio1D(xbRec.GetValue(),xbDen.GetValue(),prefix+"_acc_xb",
          disTitle+";x_{B};Event acceptance");
      auto accQ2=ratio1D(q2Rec.GetValue(),q2Den.GetValue(),prefix+"_acc_q2",
          disTitle+";Q^{2} [GeV^{2}];Event acceptance");
      auto accW=ratio1D(wRec.GetValue(),wDen.GetValue(),prefix+"_acc_w",
          disTitle+";W [GeV];Event acceptance");
      auto accY=ratio1D(yRec.GetValue(),yDen.GetValue(),prefix+"_acc_y",
          disTitle+";y;Event acceptance");
      auto accNu=ratio1D(nuRec.GetValue(),nuDen.GetValue(),prefix+"_acc_nu",
          disTitle+";#nu [GeV];Event acceptance");
      auto accXBQ2=ratio2D(xbq2Rec.GetValue(),xbq2Den.GetValue(),
          prefix+"_acc_xbq2",
          disTitle+";x_{B};Q^{2} [GeV^{2}];Event acceptance");
      accXBQ2->GetZaxis()->SetTitle("Event acceptance");
      TCanvas disCanvas(("c_"+prefix+"_dis").c_str(),
                        (topology+" DIS event acceptance").c_str(),1800,1350);
      disCanvas.Divide(3,3);
      disCanvas.cd(1); accXB->Draw("HIST");
      disCanvas.cd(2); accQ2->Draw("HIST");
      disCanvas.cd(3); accW->Draw("HIST");
      disCanvas.cd(4); accY->Draw("HIST");
      disCanvas.cd(5); accNu->Draw("HIST");
      disCanvas.cd(6); accXBQ2->Draw("COLZ");
      std::vector<std::unique_ptr<TH1D>> spectatorFrames;
      int pad = 7;
      for (const auto& s : He3DISSpectatorPlotSpecs()) {
        const std::string column=topology+"_"+s.variable;
        const std::string name=prefix+"_"+s.variable;
        auto denominator=gd.Histo1D(
            {(name+"_den").c_str(),"",s.bins,s.lo,s.hi},column);
        auto numerator=rn.Histo1D(
            {(name+"_rec").c_str(),"",s.bins,s.lo,s.hi},column);
        auto& den=denominator.GetValue();
        auto& num=numerator.GetValue();
        auto frame=ratio1D(num,den,name+"_acc",
            disTitle+";"+s.axis+";Event acceptance (REC/GEN)");
        int empty=0;
        for (int bin=1; bin<=s.bins; ++bin) {
          if (den.GetBinContent(bin)<=0) ++empty;
        }
        disCanvas.cd(pad++);
        // Match the other acceptance panels: histogram line, no error bars.
        // Preserve ratio contents/errors; only the display range is limited.
        frame->SetMinimum(0.0); frame->SetMaximum(1.0);
        frame->Draw("HIST");
        TLatex note; note.SetNDC(); note.SetTextSize(0.029);
        if (empty>0)
          note.DrawLatex(0.13,0.86,"GEN-empty bins are undefined (shown as 0)");
        if (num.GetEntries()==0)
          note.DrawLatex(0.13,0.81,"REC: no accepted entries");
        std::cout<<"[He3DIS spectator acceptance] "<<column
                 <<": "<<empty<<" bins without GEN entries\n";
        spectatorFrames.push_back(std::move(frame));
      }
      const std::string disOutput=
          outputDir+"/He3DIS_acceptance_event_"+topology+"_DIS.pdf";
      disCanvas.SaveAs(disOutput.c_str());
      std::cout<<"[PlotHe3DISAcc] wrote "<<disOutput<<'\n';
    }
  }

  // He3DIS sample-level projection. Kinematic yields may be scaled to a
  // projected sample size. GEN stays at 100% polarization; only REC spin
  // counts model partial polarization and are corrected by 1/(Pb*Pt).
  void PlotProjection(double projectionToCountRatio=1.0,
                      double beamPolarization=0.80,
                      double targetPolarization=0.60,
                      bool applyPNToLegacyCounts=true) {
    const double polarizationProduct=beamPolarization*targetPolarization;
    if (he3DISModels.size()!=2 || !(projectionToCountRatio>0.0) ||
        !std::isfinite(projectionToCountRatio) ||
        !(beamPolarization>0.0 && beamPolarization<=1.0) ||
        !(targetPolarization>0.0 && targetPolarization<=1.0)) {
      std::cerr << "[PlotProjection] requires GEN, REC and a positive finite "
                   "projection/count ratio, with Pb and Pt in (0,1].\n";
      return;
    }
    const std::string dir=outputDir+"/projection";
    fs::create_directories(dir);
    struct Spec { std::string col,axis; int bins; double lo,hi; };
    struct Pair { Spec x,y; };
    using H1=ROOT::RDF::RResultPtr<TH1D>;
    using H2=ROOT::RDF::RResultPtr<TH2D>;
    using Profile=ROOT::RDF::RResultPtr<TProfile>;
    auto nodeFor=[](He3DISModel& m,const std::string& t) -> ROOT::RDF::RNode {
      if(t=="epp") return *m.epp;
      if(t=="epn") return *m.epn;
      return *m.ed;
    };
    for(const std::string topo:{"epp","epn","ed"}) {
      std::vector<Spec> specs={{"xB","x_{B}",50,0,1},
        {"Q2","Q^{2} [GeV^{2}]",50,0,10},{"W","W [GeV]",50,1,5},
        {"y","y",50,0,1},{"nu","#nu [GeV]",50,0,11}};
      for(const auto& s:He3DISSpectatorPlotSpecs())
        specs.push_back({topo+"_"+s.variable,s.axis,s.bins,s.lo,s.hi});
      // Yields use twice the configured x_active bin count because the
      // spin-count projection below deliberately rebins by a factor of two.
      specs.push_back({topo+"_x_active","x_{N} = Q^{2}/(2p_{active}#upointq)",
                       2*he3DISXActiveBins,he3DISXActiveMin,
                       he3DISXActiveMax});
      std::vector<Pair> pairs={{specs[0],specs[1]},
                               {specs[6],specs[7]},{specs[0],specs[6]}};
      // Correlate each detected charged spectator's laboratory polar angle
      // with the topology's existing spectator momentum sum.  Keep the two
      // protons separate in epp; for ed the single spectator is the deuteron.
      // No corresponding epn panel is added while neutron detection is off.
      if(topo=="epp") {
        pairs.push_back({specs[5],
                         {"pro1_theta","#theta_{p1} [deg]",60,0,180}});
        pairs.push_back({specs[5],
                         {"pro2_theta","#theta_{p2} [deg]",60,0,180}});
      }
      else if(topo=="ed") {
        pairs.push_back({specs[5],
                         {"deu_theta","#theta_{d} [deg]",60,0,180}});
      }
      std::array<std::vector<H1>,2> yields,anti,parallel;
      std::array<std::vector<H2>,2> maps;
      std::array<std::vector<Profile>,2> meanPN;
      std::vector<ROOT::RDF::RResultPtr<ULong64_t>> totals,invalid,na,np;
      std::vector<ROOT::RDF::RResultPtr<double>> averagePN;
      // Book every action for this topology before any GetValue: GEN and REC
      // each scan once for this page set, not once per histogram/spin state.
      for(int m=0;m<2;++m) {
        auto node=nodeFor(he3DISModels[m],topo);
        auto a=node.Filter("projection_spin == -1");
        auto p=node.Filter("projection_spin == 1");
        totals.push_back(node.Count()); na.push_back(a.Count()); np.push_back(p.Count());
        invalid.push_back(node.Filter("projection_spin == 0").Count());
        averagePN.push_back(node.Mean(topo+"_PN"));
        for(const auto& s:specs) {
          const std::string name="proj_"+topo+"_"+std::to_string(m)+"_"+s.col;
          ROOT::RDF::TH1DModel model(name.c_str(),(";"+s.axis+";MC events / bin").c_str(),s.bins,s.lo,s.hi);
          yields[m].push_back(node.Histo1D(model,s.col));
          // Coarser spin-count bins only; keep yields and coverage unchanged.
          model.fNbinsX=std::max(1,s.bins/2);
          model.fName=name+"_anti"; anti[m].push_back(a.Histo1D(model,s.col));
          model.fName=name+"_parallel"; parallel[m].push_back(p.Histo1D(model,s.col));
          ROOT::RDF::TProfile1DModel pnModel(
              (name+"_mean_pn").c_str(),"",model.fNbinsX,s.lo,s.hi);
          meanPN[m].push_back(node.Profile1D(pnModel,s.col,topo+"_PN"));
        }
        for(size_t i=0;i<pairs.size();++i) {
          const auto& v=pairs[i];
          const std::string name="proj_map_"+topo+"_"+std::to_string(m)+"_"+std::to_string(i);
          maps[m].push_back(node.Histo2D({name.c_str(),(";"+v.x.axis+";"+v.y.axis+";MC events / bin").c_str(),
            v.x.bins,v.x.lo,v.x.hi,v.y.bins,v.y.lo,v.y.hi},v.x.col,v.y.col));
        }
      }
      const double ng=totals[0].GetValue(),nr=totals[1].GetValue();
      std::ostringstream summary;
      summary<<topo<<" | GEN "<<ng/1.e6<<" M | REC tagged "<<nr/1.e6<<" M | REC/GEN ";
      if(ng>0) summary<<nr/ng; else summary<<"N/A";
      std::cout<<"[PlotProjection] "<<summary.str()<<'\n';
      for(int m=0;m<2;++m)
        std::cout<<"  "<<(m==0?"GEN":"REC")<<" anti/parallel/invalid spin: "
          <<na[m].GetValue()<<" / "<<np[m].GetValue()<<" / "<<invalid[m].GetValue()
          <<"; <PN>="<<averagePN[m].GetValue()<<'\n';
      const std::string pdf=dir+"/He3DIS_projection_"+topo+".pdf";
      auto note=[](const std::string& s) {
        TLatex text; text.SetNDC(); text.SetTextSize(0.030);
        text.DrawLatex(0.12,0.92,s.c_str());
      };
      TCanvas counts(("proj_counts_"+topo).c_str(),"Sample yields",2000,1800);
      counts.Divide(4,3);
      // Scale the displayed kinematic yields to the requested projection.
      // The existing REC x10 is retained only for visibility. RDF counts,
      // spin histograms and acceptance inputs remain unmodified.
      constexpr double recYieldDisplayScale=10.0;
      std::vector<std::unique_ptr<TH1D>> owned;
      for(size_t i=0;i<specs.size();++i) {
        counts.cd(i+1);
        double maximum=0;
        for(int m=0;m<2;++m)
          maximum=std::max(maximum,projectionToCountRatio*
              (m==0?1.0:recYieldDisplayScale)*
              yields[m][i].GetValue().GetMaximum());
        TLegend legend(0.65,0.75,0.89,0.88); legend.SetBorderSize(0);
        for(int m=0;m<2;++m) {
          auto h=std::unique_ptr<TH1D>(static_cast<TH1D*>(yields[m][i].GetValue().Clone()));
          h->SetDirectory(nullptr); h->SetStats(0); h->SetLineWidth(2);
          h->Scale(projectionToCountRatio*
                   (m==0?1.0:recYieldDisplayScale));
          h->GetYaxis()->SetTitle("Projected events / bin");
          h->SetLineColor(m==0?kBlue+1:kRed+1); h->SetMinimum(0);
          h->SetMaximum(maximum>0?maximum*1.3:1);
          h->Draw(m==0?"HIST":"HIST SAME");
          legend.AddEntry(h.get(),m==0?"R_{N} #times GEN":
              "10 R_{N} #times REC tagged","l");
          owned.push_back(std::move(h));
        }
        legend.DrawClone();
      }
      counts.cd(10); note(summary.str());
      TLatex text; text.SetNDC(); text.SetTextSize(0.04);
      std::ostringstream yieldScaleNote;
      yieldScaleNote<<"Yield projection: R_{N}=N_{projection}/N_{current}="
                    <<projectionToCountRatio;
      text.DrawLatex(0.08,0.70,yieldScaleNote.str().c_str());
      text.DrawLatex(0.08,0.62,"REC yield curves additionally #times 10; counts above raw");
      if(nr==0) text.DrawLatex(0.08,0.55,"REC: no accepted tagged events");
      counts.Print((pdf+"(").c_str());
      constexpr size_t standardCoverageRows=3;
      TCanvas coverage(("proj_coverage_"+topo).c_str(),"Coverage",1800,1350);
      coverage.Divide(3,static_cast<int>(standardCoverageRows));
      std::vector<std::unique_ptr<TH2D>> owned2;
      std::vector<std::unique_ptr<TGraphErrors>> masks;
      auto drawCoverageRows=[&](TCanvas& canvas,size_t first,size_t rowCount) {
        for(size_t row=0;row<rowCount;++row) {
        const size_t i=first+row;
        auto& den=maps[0][i].GetValue(); auto& num=maps[1][i].GetValue();
        const double maximum=projectionToCountRatio*
            std::max(1.,std::max(den.GetMaximum(),num.GetMaximum()));
        for(int m=0;m<3;++m) {
          canvas.cd(3*row+m+1);
          auto h=std::unique_ptr<TH2D>(static_cast<TH2D*>((m==0?den:num).Clone()));
          h->SetDirectory(nullptr); h->SetStats(0);
          h->SetTitle(m==0?"GEN":m==1?"REC tagged":"REC/GEN; ; ");
          h->GetXaxis()->SetTitle(pairs[i].x.axis.c_str());
          h->GetYaxis()->SetTitle(pairs[i].y.axis.c_str());
          if(m<2) {
            h->Scale(projectionToCountRatio);
            h->SetMinimum(0.5);h->SetMaximum(maximum);
            h->GetZaxis()->SetTitle("Projected events / bin");gPad->SetLogz();
          }
          else {h->Divide(&den);h->SetMinimum(0);h->SetMaximum(0.5);h->GetZaxis()->SetRangeUser(0.,0.5);h->GetZaxis()->SetTitle("REC/GEN (includes migration)");}
          h->Draw("COLZ");
          if(m==2) {
            auto mask=std::make_unique<TGraphErrors>();
            for(int x=1;x<=den.GetNbinsX();++x) for(int y=1;y<=den.GetNbinsY();++y)
              if(den.GetBinContent(x,y)==0) mask->SetPoint(mask->GetN(),den.GetXaxis()->GetBinCenter(x),den.GetYaxis()->GetBinCenter(y));
            mask->SetMarkerStyle(5);mask->SetMarkerSize(0.35);mask->SetMarkerColor(kGray+1);mask->Draw("P SAME");
            note("Gray crosses: no GEN statistics"); masks.push_back(std::move(mask));
          }
          owned2.push_back(std::move(h));
        }
        }
      };
      drawCoverageRows(coverage,0,standardCoverageRows);
      coverage.Print(pdf.c_str());
      std::unique_ptr<TCanvas> angleCoverage;
      if(pairs.size()>standardCoverageRows) {
        const size_t angleRows=pairs.size()-standardCoverageRows;
        angleCoverage=std::make_unique<TCanvas>(
            ("proj_theta_sum_p_"+topo).c_str(),"Theta versus momentum sum",
            1800,600*static_cast<int>(angleRows));
        angleCoverage->Divide(3,static_cast<int>(angleRows));
        drawCoverageRows(*angleCoverage,standardCoverageRows,angleRows);
        angleCoverage->Print(pdf.c_str());
      }
      TCanvas asymRaw(("proj_asym_raw_"+topo).c_str(),
                      "Tagged longitudinal double spin",2000,1800);
      TCanvas asymOverPN(("proj_asym_over_pn_"+topo).c_str(),
                         "Nucleon longitudinal double spin",2000,1800);
      asymRaw.Divide(4,3); asymOverPN.Divide(4,3);
      std::array<std::vector<std::unique_ptr<TH1D>>,2> frames;
      std::array<std::vector<std::unique_ptr<TGraphErrors>>,2> graphs;
      for(size_t i=0;i<specs.size();++i) {
        const std::string projectionKey=topo+"/"+specs[i].col;
        auto& savedRawErrors=he3DISProjectionErrors[projectionKey];
        auto& savedOverPNErrors=he3DISProjectionErrorsOverPN[projectionKey];
        auto& savedOverPNValues=he3DISProjectionValuesOverPN[projectionKey];
        for(auto& errors:savedRawErrors)
          errors.assign(anti[0][i].GetValue().GetNbinsX()+1,
                        std::numeric_limits<double>::quiet_NaN());
        for(auto& errors:savedOverPNErrors)
          errors.assign(anti[0][i].GetValue().GetNbinsX()+1,
                        std::numeric_limits<double>::quiet_NaN());
        for(auto& values:savedOverPNValues)
          values.assign(anti[0][i].GetValue().GetNbinsX()+1,
                        std::numeric_limits<double>::quiet_NaN());
        std::array<double,2> extent{{0.02,0.02}};
        std::array<std::array<TGraphErrors*,2>,2> drawn{};
        for(int m=0;m<2;++m) {
          auto rawGraph=std::make_unique<TGraphErrors>();
          auto overPNGraph=std::make_unique<TGraphErrors>();
          // Both configurations must exist in this model. Invalid labels are
          // excluded only here, never silently guessed or removed from yields.
          if(na[m].GetValue()>0 && np[m].GetValue()>0) {
            for(int bin=1;bin<=anti[0][i].GetValue().GetNbinsX();++bin) {
              const double a=anti[m][i].GetValue().GetBinContent(bin);
              const double p=parallel[m][i].GetValue().GetBinContent(bin);
              const double n=a+p;if(n<=0) continue;
              const double exposurePolarization=m==0 ? 1.0 : polarizationProduct;
              const double countAsymmetry=(a-p)/n;
              const double baseValue=countAsymmetry/exposurePolarization;
              const double pn=meanPN[m][i].GetValue().GetBinContent(bin);
              if(!std::isfinite(pn)) continue;
              // ON repairs legacy samples generated with PN=1. OFF assumes
              // the spin counts already contain the physical PN dependence.
              const double taggedValue=applyPNToLegacyCounts
                                           ? pn*baseValue : baseValue;
              const double observedValue=exposurePolarization*taggedValue;
              const double taggedError=std::sqrt(
                  std::max(0.,1-observedValue*observedValue)/n)/
                  exposurePolarization;
              savedRawErrors[m][bin]=taggedError;
              rawGraph->SetPoint(rawGraph->GetN(),
                                 anti[0][i].GetValue().GetBinCenter(bin),
                                 taggedValue);
              rawGraph->SetPointError(rawGraph->GetN()-1,0,taggedError);
              extent[0]=std::max(extent[0],
                  1.15*std::abs(taggedValue));
              if(std::abs(pn)>1.e-6) {
                const double nucleonValue=taggedValue/pn;
                const double nucleonError=taggedError/std::abs(pn);
                savedOverPNValues[m][bin]=nucleonValue;
                savedOverPNErrors[m][bin]=nucleonError;
                overPNGraph->SetPoint(overPNGraph->GetN(),
                    anti[0][i].GetValue().GetBinCenter(bin),nucleonValue);
                overPNGraph->SetPointError(overPNGraph->GetN()-1,0,nucleonError);
                extent[1]=std::max(extent[1],
                    1.15*std::abs(nucleonValue));
              }
            }
          }
          for(auto* graph:{rawGraph.get(),overPNGraph.get()}) {
            graph->SetMarkerStyle(m==0?20:24);graph->SetMarkerSize(0.6);
            graph->SetLineColor(m==0?kBlue+1:kRed+1);
            graph->SetMarkerColor(m==0?kBlue+1:kRed+1);
          }
          drawn[0][m]=rawGraph.get();drawn[1][m]=overPNGraph.get();
          graphs[0].push_back(std::move(rawGraph));
          graphs[1].push_back(std::move(overPNGraph));
        }
        he3DISProjectionYExtents[projectionKey]=extent[0];
        he3DISProjectionYExtentsOverPN[projectionKey]=extent[1];
        for(int view=0;view<2;++view) {
          (view==0?asymRaw:asymOverPN).cd(i+1);
          auto frame=std::unique_ptr<TH1D>(static_cast<TH1D*>(
              anti[0][i].GetValue().Clone()));
          frame->SetDirectory(nullptr);frame->Reset();frame->SetStats(0);
          frame->GetYaxis()->SetTitle(view==0
              ? "A_{#parallel}^{raw} (P_{N} retained)"
              : "A_{#parallel}^{raw}/P_{N}");
          frame->SetMinimum(-extent[view]);frame->SetMaximum(extent[view]);
          frame->Draw("AXIS");
          TLine zero(specs[i].lo,0,specs[i].hi,0);
          zero.SetLineStyle(2);zero.SetLineColor(kGray+2);zero.DrawClone();
          TLegend legend(0.64,0.75,0.89,0.88);legend.SetBorderSize(0);
          for(int m=0;m<2;++m) {
            if(drawn[view][m]->GetN()>0) drawn[view][m]->Draw("P SAME");
            std::string label=m==0?"GEN":"REC tagged";
            if(drawn[view][m]->GetN()==0) label+=" (unavailable)";
            legend.AddEntry(drawn[view][m],label.c_str(),"lp");
          }
          legend.DrawClone();frames[view].push_back(std::move(frame));
        }
      }
      for(int view=0;view<2;++view) {
        auto& canvas=view==0?asymRaw:asymOverPN;
        canvas.cd(10);note(summary.str());
        text.SetTextSize(0.035);
        std::ostringstream polarizationNote;
        polarizationNote<<"P_{b}="<<beamPolarization<<", P_{t}="
                        <<targetPolarization<<", P_{b}P_{t}="
                        <<polarizationProduct;
        text.DrawLatex(0.07,0.79,polarizationNote.str().c_str());
        text.DrawLatex(0.07,0.67,view==0
            ? "Tagged A_{#parallel}: P_{N} retained"
            : "Extracted nucleon A_{#parallel}: divided by P_{N}");
        text.DrawLatex(0.07,0.55,applyPNToLegacyCounts
            ? "PN projection ON: legacy count centers repaired"
            : "PN projection OFF: input counts assumed to include PN");
        text.DrawLatex(0.07,0.43,"GEN: 100% polarization; REC corrected for P_{b}P_{t}");
        text.DrawLatex(0.07,0.31,"Generator-copy correlations are not included");
        text.DrawLatex(0.07,0.19,"Near-zero P_{N} bins excluded after PN division");
      }
      asymRaw.Print(pdf.c_str());
      asymOverPN.Print((pdf+")").c_str());
      std::cout<<"[PlotProjection] wrote "<<pdf<<'\n';
      // Pads keep non-owning references to our unique_ptr objects. Detach all
      // primitives before those objects are destroyed at the end of the loop.
      asymRaw.Clear(); asymOverPN.Clear();
      if(angleCoverage) angleCoverage->Clear();
      coverage.Clear(); counts.Clear();
    }
  }

  // Smooth generator-model centers evaluated at the mean kinematics of each
  // GEN or tagged-REC bin. This is intentionally separate from PlotProjection,
  // whose asymmetry points continue to come from finite spin-state counts.
  void PlotGeneratorModelProjection(
      const std::function<double(bool,double,double,double,double,double,double)>&
          evaluateAparallel,
      double projectionToCountRatio=1.0,
      double beamPolarization=0.80,
      double targetPolarization=0.60,
      bool applyPNToLegacyCounts=true) {
    const double polarizationProduct=beamPolarization*targetPolarization;
    // This flag changes legacy count centers in PlotProjection. The smooth
    // model always draws both PN*A_N and A_N by construction.
    (void)applyPNToLegacyCounts;
    if (he3DISModels.size()!=2 || !evaluateAparallel ||
        !(projectionToCountRatio>0.0) ||
        !std::isfinite(projectionToCountRatio) ||
        !(beamPolarization>0.0 && beamPolarization<=1.0) ||
        !(targetPolarization>0.0 && targetPolarization<=1.0)) {
      std::cerr << "[PlotGeneratorModelProjection] requires GEN, REC, a model "
                   "evaluator, and Pb/Pt in (0,1].\n";
      return;
    }
    const double errorScale=1.0/std::sqrt(projectionToCountRatio);
    const double modelXMax=he3DISModelXMax;
    const std::string dir=outputDir+"/projection_model";
    fs::create_directories(dir);
    struct Spec {std::string col,axis;int bins;double lo,hi;};
    using Profile=ROOT::RDF::RResultPtr<TProfile>;
    auto nodeFor=[](He3DISModel& m,const std::string& t)->ROOT::RDF::RNode {
      if(t=="epp") return *m.epp;
      if(t=="epn") return *m.epn;
      return *m.ed;
    };
    for(const std::string topo:{"epp","epn","ed"}) {
      std::vector<Spec> specs={{"xB","x_{B}",25,0,1},
        {"Q2","Q^{2} [GeV^{2}]",25,0,10},{"W","W [GeV]",25,1,5},
        {"y","y",25,0,1},{"nu","#nu [GeV]",25,0,11}};
      for(const auto& s:He3DISSpectatorPlotSpecs())
        specs.push_back({topo+"_"+s.variable,s.axis,std::max(1,s.bins/2),s.lo,s.hi});
      specs.push_back({topo+"_x_active","x_{N} = Q^{2}/(2p_{active}#upointq)",
                       he3DISXActiveBins,he3DISXActiveMin,
                       he3DISXActiveMax});
      std::array<std::vector<Profile>,2> meanXB,meanQ2,meanY,meanAlpha,
          meanActiveX,meanCosQ,meanPN;
      std::vector<ROOT::RDF::RResultPtr<ULong64_t>> totals;
      for(int m=0;m<2;++m) {
        auto node=nodeFor(he3DISModels[m],topo);
        totals.push_back(node.Count());
        for(size_t i=0;i<specs.size();++i) {
          const auto& s=specs[i];
          const std::string base="model_center_"+topo+"_"+std::to_string(m)+"_"+std::to_string(i);
          auto model=[&](const std::string& suffix) {
            return ROOT::RDF::TProfile1DModel((base+suffix).c_str(),"",s.bins,s.lo,s.hi);
          };
          meanXB[m].push_back(node.Profile1D(model("_xb"),s.col,"xB"));
          meanQ2[m].push_back(node.Profile1D(model("_q2"),s.col,"Q2"));
          meanY[m].push_back(node.Profile1D(model("_y"),s.col,"y"));
          meanAlpha[m].push_back(node.Profile1D(model("_alpha"),s.col,topo+"_alpha_active"));
          meanActiveX[m].push_back(node.Profile1D(
              model("_xactive"),s.col,topo+"_x_active"));
          meanCosQ[m].push_back(node.Profile1D(model("_cosq"),s.col,"model_cos_theta_q"));
          meanPN[m].push_back(node.Profile1D(model("_pn"),s.col,topo+"_PN"));
        }
      }
      const double ngen=totals[0].GetValue(),nrec=totals[1].GetValue();
      for(int view=0;view<2;++view) {
      const std::string viewName=view==0?"raw":"raw_over_pn";
      TCanvas canvas(("c_model_center_"+viewName+"_"+topo).c_str(),
                     view==0?"Tagged generator model centers":
                             "Nucleon generator model centers",2000,1800);
      canvas.Divide(4,3);
      std::vector<std::unique_ptr<TH1D>> frames;
      std::vector<std::unique_ptr<TGraphErrors>> graphs;
      for(size_t i=0;i<specs.size();++i) {
        canvas.cd(i+1);gPad->SetTicks();
        auto frame=std::make_unique<TH1D>(
            ("frame_model_"+viewName+"_"+topo+"_"+std::to_string(i)).c_str(),
            (";"+specs[i].axis+
             (view==0?";Generator-model A_{#parallel}^{raw} (P_{N} retained)":
                      ";Generator-model A_{#parallel}^{raw}/P_{N}")).c_str(),
            specs[i].bins,specs[i].lo,specs[i].hi);
        frame->SetDirectory(nullptr);frame->SetStats(0);
        // Determine the vertical range from the values that are actually
        // returned by the generator model. A fixed +/-0.05 range made small
        // asymmetries look indistinguishable from zero.
        double minValue=std::numeric_limits<double>::infinity();
        double maxValue=-std::numeric_limits<double>::infinity();
        double modelExtent=0.02;
        std::array<int,2> validPoints{{0,0}};
        const std::string projectionKey=topo+"/"+specs[i].col;
        const auto& errorMap=view==0 ? he3DISProjectionErrors
                                    : he3DISProjectionErrorsOverPN;
        const auto savedErrorIt=errorMap.find(projectionKey);
        std::array<TGraphErrors*,2> drawn{};
        TLegend legend(0.52,0.73,0.89,0.89);legend.SetBorderSize(0);legend.SetFillStyle(0);
        for(int m=0;m<2;++m) {
          auto graph=std::make_unique<TGraphErrors>();
          graph->SetName(("g_model_"+viewName+"_"+topo+"_"+
              std::to_string(i)+"_"+std::to_string(m)).c_str());
          auto& xb=meanXB[m][i].GetValue();auto& q2=meanQ2[m][i].GetValue();
          auto& yy=meanY[m][i].GetValue();auto& aa=meanAlpha[m][i].GetValue();
          auto& xx=meanActiveX[m][i].GetValue();
          auto& cq=meanCosQ[m][i].GetValue();
          auto& pnProfile=meanPN[m][i].GetValue();
          for(int bin=1;bin<=specs[i].bins;++bin) {
            if(xb.GetBinEntries(bin)<=0) continue;
            const double alphaMean=aa.GetBinContent(bin);
            if(!(alphaMean>0.0)) continue;
            const double modelX=xx.GetBinContent(bin);
            // JAMpol25 is constrained only through approximately x=0.85.
            // Do not display its endpoint extrapolation as a prediction.
            if(!(modelX>0.0 && modelX<modelXMax)) continue;
            const double cosq=std::clamp(cq.GetBinContent(bin),-1.0,1.0);
            // The external API forms x_N=xB/alpha internally and also uses
            // xB itself in finite-Q2 factors.  Preserve the measured xB and
            // adapt only its alpha argument so the internal quotient is the
            // exact invariant x_N.
            const double measuredXB=xb.GetBinContent(bin);
            const double modelAlpha=measuredXB/modelX;
            const double nucleonValue=evaluateAparallel(topo=="epp",
                measuredXB,modelAlpha,q2.GetBinContent(bin),
                yy.GetBinContent(bin),-cosq,std::sqrt(std::max(0.0,1.0-cosq*cosq)));
            const double pn=pnProfile.GetBinContent(bin);
            const double value=view==0 ? pn*nucleonValue : nucleonValue;
            if(!std::isfinite(value)) continue;
            if(savedErrorIt==errorMap.end() ||
               bin>=static_cast<int>(savedErrorIt->second[m].size())) continue;
            const double countError=savedErrorIt->second[m][bin];
            if(!std::isfinite(countError)) continue;
            graph->SetPoint(graph->GetN(),xb.GetBinCenter(bin),value);
            graph->SetPointError(graph->GetN()-1,0.0,countError*errorScale);
            modelExtent=std::max(modelExtent,
                1.15*std::abs(value));
            ++validPoints[m];
            minValue=std::min(minValue,value);
            maxValue=std::max(maxValue,value);
          }
          graph->SetLineColor(m==0?kBlue+1:kRed+1);
          graph->SetMarkerColor(m==0?kBlue+1:kRed+1);
          graph->SetMarkerStyle(m==0?20:24);graph->SetMarkerSize(0.6);
          std::string label=m==0?"Model at GEN-bin means":"Model at REC-tagged-bin means";
          if(graph->GetN()==0) label+=" (unavailable)";
          legend.AddEntry(graph.get(),label.c_str(),"p");
          drawn[m]=graph.get();graphs.push_back(std::move(graph));
        }
        // Reuse the corresponding count-projection range so the two PDF sets
        // remain directly comparable. Fall back to the model points if the
        // count plot has not been called first.
        const auto& rangeMap=view==0 ? he3DISProjectionYExtents
                                    : he3DISProjectionYExtentsOverPN;
        const auto savedRange=rangeMap.find(projectionKey);
        const double extent=savedRange!=rangeMap.end()
                                ? savedRange->second : modelExtent;
        frame->SetMinimum(-extent);frame->SetMaximum(extent);frame->Draw("AXIS");
        TLine zero(specs[i].lo,0,specs[i].hi,0);
        zero.SetLineStyle(2);zero.SetLineColor(kGray+2);zero.DrawClone();
        for(auto* graph:drawn) if(graph && graph->GetN()>0) graph->Draw("P SAME");
        if(validPoints[0]+validPoints[1]==0) {
          TLatex unavailable;unavailable.SetNDC();unavailable.SetTextColor(kRed+1);
          unavailable.SetTextSize(0.055);
          unavailable.DrawLatex(0.20,0.52,"No valid model points");
        }
        legend.DrawClone();frames.push_back(std::move(frame));
        std::cout<<"[PlotGeneratorModelProjection] "<<topo<<" "<<specs[i].col
                 <<": valid GEN/REC points = "<<validPoints[0]<<"/"
                 <<validPoints[1];
        if(validPoints[0]+validPoints[1]>0)
          std::cout<<", A_parallel range = ["<<minValue<<", "<<maxValue<<"]";
        std::cout<<'\n';
      }
      canvas.cd(10);
      TLatex text;text.SetNDC();text.SetTextSize(0.036);
      text.DrawLatex(0.07,0.82,(topo+(view==0
          ? " tagged generator-model centers"
          : " extracted nucleon generator-model centers")).c_str());
      std::ostringstream counts;counts<<"GEN "<<ngen/1.e6<<" M, REC tagged "<<nrec/1.e6<<" M";
      text.DrawLatex(0.07,0.70,counts.str().c_str());
      text.DrawLatex(0.07,0.57,"Evaluated at bin-mean reconstructed kinematics");
      std::ostringstream errorNote;
      errorNote<<"Count-projection errors / sqrt(R_{N}), R_{N}="
               <<projectionToCountRatio;
      text.DrawLatex(0.07,0.44,errorNote.str().c_str());
      std::ostringstream polarizationNote;
      polarizationNote<<"REC correction: P_{b}="<<beamPolarization
                      <<", P_{t}="<<targetPolarization
                      <<", P_{b}P_{t}="<<polarizationProduct;
      text.DrawLatex(0.07,0.31,polarizationNote.str().c_str());
      text.DrawLatex(0.07,0.19,view==0
          ? "Tagged result: generator-model A_{#parallel} #times P_{N}"
          : "Nucleon result: tagged A_{#parallel} divided by P_{N}");
      const std::string output=dir+"/He3DIS_generator_model_asymmetry_"+
          viewName+"_"+topo+".pdf";
      canvas.SaveAs(output.c_str());
      std::cout<<"[PlotGeneratorModelProjection] wrote "<<output<<'\n';
      canvas.Clear();
      }
    }
  }

  // Extract nucleon A1 from the exact finite-count A_parallel/PN values saved
  // by PlotProjection.  The plotted center always uses A_perp=0.  The absolute
  // shift obtained when JAM25 supplies A_perp/d=A2-xi*A1 is drawn as a
  // symmetric, semi-transparent systematic-error box around that center.
  // Each topology PDF contains a current-GEMC page with finite-count centers
  // and a projected-statistics page with smooth JAM25 A_parallel centers.
  // The latter scales count-derived errors by
  // 1/sqrt(projectionToCountRatio).
  void PlotA1ProjectionHe3DIS(
      const std::function<double(bool,double,double,double,double,double,double)>&
          evaluateAparallel,
      const std::function<std::array<double,7>(bool,double,double,double)>&
          evaluateNucleon,
      double projectionToCountRatio=1.0) {
    if(he3DISModels.size()!=2 || !evaluateAparallel || !evaluateNucleon ||
       !(projectionToCountRatio>0.0) ||
       !std::isfinite(projectionToCountRatio)) {
      std::cerr<<"[PlotA1ProjectionHe3DIS] requires GEN, REC, the nucleon "
                   "model, and a positive projection/count ratio.\n";
      return;
    }
    if(he3DISProjectionValuesOverPN.empty() ||
       he3DISProjectionErrorsOverPN.empty()) {
      std::cerr<<"[PlotA1ProjectionHe3DIS] call PlotProjection first so the "
                   "same finite-count A_parallel values can be reused.\n";
      return;
    }
    const std::string dir=outputDir+"/A1_projection";
    fs::create_directories(dir);
    const double projectedErrorScale=1.0/std::sqrt(projectionToCountRatio);
    const double modelXMax=he3DISModelXMax;
    constexpr double nucleonMass=0.93892;
    struct Spec {std::string col,axis;int bins;double lo,hi;};
    using Profile=ROOT::RDF::RResultPtr<TProfile>;
    auto nodeFor=[](He3DISModel& m,const std::string& t)->ROOT::RDF::RNode {
      return t=="epp" ? *m.epp : *m.ed;
    };
    auto finite=[](const std::array<double,7>& v) {
      return std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);});
    };
    // Return D, eta and xi in the convention of Sec. 8 of arXiv:2106.08805.
    auto factors=[&](double xb,double q2,double y,double r) {
      std::array<double,3> out{{NAN,NAN,NAN}};
      if(!(xb>0.0) || !(q2>0.0) || !(y>0.0 && y<1.0) ||
         !std::isfinite(r)) return out;
      const double gamma=2.0*nucleonMass*xb/std::sqrt(q2);
      const double gamma2=gamma*gamma;
      const double common=4.0*(1.0-y)-gamma2*y*y;
      const double twoPlus=2.0+gamma2*y;
      const double dDen=2.0*(1.0+gamma2)*y*y+common*(1.0+r);
      if(!(dDen>0.0) || !(twoPlus>0.0) || !(2.0-y>0.0)) return out;
      out[0]=y*(2.0-y)*twoPlus/dDen;
      out[1]=gamma*common/((2.0-y)*twoPlus);
      out[2]=gamma*(2.0-y)/twoPlus;
      return out;
    };
    for(const std::string topo:{"epp","ed"}) {
      std::vector<Spec> specs={{"xB","x_{B}",25,0,1},
        {"Q2","Q^{2} [GeV^{2}]",25,0,10},{"W","W [GeV]",25,1,5},
        {"y","y",25,0,1},{"nu","#nu [GeV]",25,0,11}};
      for(const auto& s:He3DISSpectatorPlotSpecs())
        specs.push_back({topo+"_"+s.variable,s.axis,std::max(1,s.bins/2),s.lo,s.hi});
      specs.push_back({topo+"_x_active","x_{N} = Q^{2}/(2p_{active}#upointq)",
                       he3DISXActiveBins,he3DISXActiveMin,
                       he3DISXActiveMax});
      std::array<std::vector<Profile>,2> meanXB,meanQ2,meanY,meanAlpha,
          meanActiveX,meanCosQ;
      for(int m=0;m<2;++m) {
        auto node=nodeFor(he3DISModels[m],topo);
        for(size_t i=0;i<specs.size();++i) {
          const auto& s=specs[i];
          const std::string base="a1_mean_"+topo+"_"+std::to_string(m)+"_"+
                                 std::to_string(i);
          auto model=[&](const std::string& suffix) {
            return ROOT::RDF::TProfile1DModel((base+suffix).c_str(),"",
                                              s.bins,s.lo,s.hi);
          };
          meanXB[m].push_back(node.Profile1D(model("_xb"),s.col,"xB"));
          meanQ2[m].push_back(node.Profile1D(model("_q2"),s.col,"Q2"));
          meanY[m].push_back(node.Profile1D(model("_y"),s.col,"y"));
          meanAlpha[m].push_back(node.Profile1D(
              model("_alpha"),s.col,topo+"_alpha_active"));
          meanActiveX[m].push_back(node.Profile1D(
              model("_xactive"),s.col,topo+"_x_active"));
          meanCosQ[m].push_back(node.Profile1D(
              model("_cosq"),s.col,"model_cos_theta_q"));
        }
      }
      TCanvas current(("a1_current_"+topo).c_str(),
                      "Current GEMC A1",2200,1800);
      TCanvas projected(("a1_projected_"+topo).c_str(),
                        "Projected A1",2200,1800);
      current.Divide(4,3); projected.Divide(4,3);
      std::vector<std::unique_ptr<TH1D>> frames;
      std::vector<std::unique_ptr<TGraphErrors>> points;
      std::vector<std::unique_ptr<TGraph>> curves;
      int skippedProjectedLargeError=0;
      for(size_t i=0;i<specs.size();++i) {
        const std::string key=topo+"/"+specs[i].col;
        const auto valueIt=he3DISProjectionValuesOverPN.find(key);
        const auto errorIt=he3DISProjectionErrorsOverPN.find(key);
        if(valueIt==he3DISProjectionValuesOverPN.end() ||
           errorIt==he3DISProjectionErrorsOverPN.end()) continue;
        auto genCurrent=std::make_unique<TGraphErrors>();
        auto recCurrent=std::make_unique<TGraphErrors>();
        auto recProjected=std::make_unique<TGraphErrors>();
        auto jamCurrent=std::make_unique<TGraph>();
        auto jamProjected=std::make_unique<TGraph>();
        struct BoxSpec {double x1,y1,x2,y2;int color;};
        std::vector<BoxSpec> currentBoxes,projectedBoxes;
        double extent=0.02;
        const double boxHalfWidth=0.35*(specs[i].hi-specs[i].lo)/specs[i].bins;
        for(int m=0;m<2;++m) {
          auto& xb=meanXB[m][i].GetValue();
          auto& q2=meanQ2[m][i].GetValue();
          auto& yy=meanY[m][i].GetValue();
          auto& aa=meanAlpha[m][i].GetValue();
          auto& xx=meanActiveX[m][i].GetValue();
          auto& cq=meanCosQ[m][i].GetValue();
          const auto& values=valueIt->second[m];
          const auto& errors=errorIt->second[m];
          for(int bin=1;bin<=specs[i].bins;++bin) {
            if(xb.GetBinEntries(bin)<=0 ||
               bin>=static_cast<int>(values.size()) ||
               bin>=static_cast<int>(errors.size()) ||
               !std::isfinite(values[bin]) ||
               !std::isfinite(errors[bin])) continue;
            const double xbMean=xb.GetBinContent(bin);
            const double alphaMean=aa.GetBinContent(bin);
            if(!(alphaMean>0.0)) continue;
            const double modelX=xx.GetBinContent(bin);
            if(!(modelX>0.0 && modelX<modelXMax)) continue;
            const double modelAlpha=xbMean/modelX;
            const auto model=evaluateNucleon(topo=="epp",
                xbMean,modelAlpha,
                q2.GetBinContent(bin));
            if(!finite(model)) continue;
            const auto k=factors(xb.GetBinContent(bin),q2.GetBinContent(bin),
                                 yy.GetBinContent(bin),model[4]);
            if(!std::isfinite(k[0]) || !std::isfinite(k[1]) ||
               !std::isfinite(k[2])) continue;
            const double mixing=1.0+k[1]*k[2];
            if(std::abs(k[0]*mixing)<1.e-12) continue;

            // Central value: set A_T (A_perp) to zero.  The full difference
            // to the JAM25 A_T correction is treated as a symmetric systematic
            // uncertainty; it is not divided by two because the central value
            // deliberately remains the A_T=0 result.
            const double a1Zero=values[bin]/(k[0]*mixing);
            const double transverseShift=
                -k[1]*(model[6]-k[2]*model[5])/mixing;
            const double systematic=std::abs(transverseShift);
            const double error=errors[bin]/std::abs(k[0]*mixing);
            const double x=xb.GetBinCenter(bin);
            auto* graph=m==0 ? genCurrent.get() : recCurrent.get();
            graph->SetPoint(graph->GetN(),x,a1Zero);
            graph->SetPointError(graph->GetN()-1,0,error);
            currentBoxes.push_back({x-boxHalfWidth,a1Zero-systematic,
                                    x+boxHalfWidth,a1Zero+systematic,
                                    m==0 ? kBlue+1 : kRed+1});

            if(m==1) {
              // Target-exposure center uses the same smooth generator-model
              // A_parallel as PlotGeneratorModelProjection, evaluated at the
              // REC-bin mean kinematics.  Its displayed center still assumes
              // A_T=0; the possible JAM25 correction is the shaded box.
              const double cosq=std::clamp(cq.GetBinContent(bin),-1.0,1.0);
              const double modelAparallel=evaluateAparallel(topo=="epp",
                  xbMean,modelAlpha,
                  q2.GetBinContent(bin),yy.GetBinContent(bin),-cosq,
                  std::sqrt(std::max(0.0,1.0-cosq*cosq)));
              if(std::isfinite(modelAparallel)) {
                const double projectedA1Zero=modelAparallel/(k[0]*mixing);
                const double projectedError=error*projectedErrorScale;
                // Very large A1 uncertainties are non-informative and obscure
                // the projection.  Omit those points and their boxes together.
                if(std::isfinite(projectedError) && projectedError<=1.0) {
                  recProjected->SetPoint(recProjected->GetN(),x,projectedA1Zero);
                  recProjected->SetPointError(recProjected->GetN()-1,0,
                                               projectedError);
                  projectedBoxes.push_back({x-boxHalfWidth,
                      projectedA1Zero-systematic,x+boxHalfWidth,
                      projectedA1Zero+systematic,kRed+1});
                  if(std::abs(projectedA1Zero)<=1.0)
                    extent=std::max(extent,1.15*std::abs(projectedA1Zero));
                } else {
                  ++skippedProjectedLargeError;
                }
              }
              jamProjected->SetPoint(jamProjected->GetN(),x,model[5]);
            } else {
              jamCurrent->SetPoint(jamCurrent->GetN(),x,model[5]);
            }
            // Statistical errors, systematic boxes, the JAM25 curve, and
            // fluctuation-driven |A1|>1 points do not set the display zoom.
            if(std::abs(a1Zero)<=1.0)
              extent=std::max(extent,1.15*std::abs(a1Zero));
          }
        }
        for(auto* graph:{genCurrent.get(),recCurrent.get(),recProjected.get()}) {
          graph->SetLineStyle(1);graph->SetLineWidth(1);
        }
        genCurrent->SetMarkerStyle(20);genCurrent->SetMarkerSize(0.55);
        genCurrent->SetMarkerColor(kBlue+1);genCurrent->SetLineColor(kBlue+1);
        recCurrent->SetMarkerStyle(24);recCurrent->SetMarkerSize(0.55);
        recCurrent->SetMarkerColor(kRed+1);recCurrent->SetLineColor(kRed+1);
        recProjected->SetMarkerStyle(20);recProjected->SetMarkerSize(0.60);
        recProjected->SetMarkerColor(kRed+1);recProjected->SetLineColor(kRed+1);
        for(auto* curve:{jamCurrent.get(),jamProjected.get()}) {
          curve->SetLineColor(kBlack);curve->SetLineWidth(2);
        }
        // Keep every finite A1 point in the display graph, including
        // fluctuation-driven |A1|>1 values, without allowing an overflow
        // point or its statistical error bar to paint outside the frame.
        // The un-clipped centers above still determine all physics results;
        // this is only a display transform applied after the y range has been
        // fixed from the |A1|<=1 centers.  Leave a small inset so an overflow
        // marker remains visible rather than being cut in half at the frame.
        const double markerLimit=0.985*extent;
        auto clipPointsToFrame=[&](TGraphErrors& graph) {
          for(int point=0;point<graph.GetN();++point) {
            double x=0.0,y=0.0;
            graph.GetPoint(point,x,y);
            const double displayY=std::clamp(y,-markerLimit,markerLimit);
            const double room=std::max(0.0,extent-std::abs(displayY));
            const double displayEY=std::min(graph.GetErrorY(point),room);
            graph.SetPoint(point,x,displayY);
            graph.SetPointError(point,graph.GetErrorX(point),displayEY);
          }
        };
        clipPointsToFrame(*genCurrent);
        clipPointsToFrame(*recCurrent);
        clipPointsToFrame(*recProjected);
        for(int page=0;page<2;++page) {
          (page==0?current:projected).cd(static_cast<int>(i)+1);
          gPad->SetTicks();
          auto frame=std::make_unique<TH1D>(
              ("frame_a1_"+topo+"_"+std::to_string(i)+"_"+
               std::to_string(page)).c_str(),
              (";"+specs[i].axis+";A_{1}^{"+
               std::string(topo=="epp"?"n":"p")+"}").c_str(),
              specs[i].bins,specs[i].lo,specs[i].hi);
          frame->SetDirectory(nullptr);frame->SetStats(0);
          frame->SetMinimum(-extent);frame->SetMaximum(extent);frame->Draw("AXIS");
          TLine zero(specs[i].lo,0,specs[i].hi,0);
          zero.SetLineStyle(2);zero.SetLineColor(kGray+2);zero.DrawClone();
          const auto& boxes=page==0 ? currentBoxes : projectedBoxes;
          for(const auto& b:boxes) {
            // TBox itself is not reliably clipped by every ROOT PDF backend.
            // Center it on the same display-clipped point and symmetrically
            // shorten its half-height when the full systematic interval would
            // cross the frame.  The point therefore remains at the visible
            // box center even for an overflow value.
            const double boxCenter=0.5*(b.y1+b.y2);
            const double displayCenter=
                std::clamp(boxCenter,-markerLimit,markerLimit);
            const double originalHalfHeight=0.5*std::abs(b.y2-b.y1);
            const double displayHalfHeight=std::min(
                originalHalfHeight,std::max(0.0,extent-std::abs(displayCenter)));
            if(!(displayHalfHeight>0.0)) continue;
            TBox box(b.x1,displayCenter-displayHalfHeight,
                     b.x2,displayCenter+displayHalfHeight);
            box.SetFillColorAlpha(b.color,0.22);
            box.SetLineColor(b.color);box.SetLineWidth(1);box.DrawClone("F SAME");
          }
          if(page==0) {
            if(jamCurrent->GetN()>1) jamCurrent->Draw("L SAME");
            if(genCurrent->GetN()>0) genCurrent->Draw("P SAME");
            if(recCurrent->GetN()>0) recCurrent->Draw("P SAME");
          } else {
            if(jamProjected->GetN()>1) jamProjected->Draw("L SAME");
            if(recProjected->GetN()>0) recProjected->Draw("P SAME");
          }
          frames.push_back(std::move(frame));
        }
        points.push_back(std::move(genCurrent));
        points.push_back(std::move(recCurrent));
        points.push_back(std::move(recProjected));
        curves.push_back(std::move(jamCurrent));
        curves.push_back(std::move(jamProjected));
      }

      // Put one shared legend after the nine distributions.  A_T avoids
      // ROOT's unreliable #perp PDF glyph.
      TGraphErrors dummyGen,dummyRec,dummyProjected;
      dummyGen.SetMarkerStyle(20);dummyGen.SetMarkerColor(kBlue+1);
      dummyRec.SetMarkerStyle(24);dummyRec.SetMarkerColor(kRed+1);
      dummyProjected.SetMarkerStyle(20);dummyProjected.SetMarkerColor(kRed+1);
      TGraph dummyJam;dummyJam.SetLineColor(kBlack);dummyJam.SetLineWidth(2);
      TBox dummyGenBox,dummyRecBox;
      dummyGenBox.SetFillColorAlpha(kBlue+1,0.22);dummyGenBox.SetLineColor(kBlue+1);
      dummyRecBox.SetFillColorAlpha(kRed+1,0.22);dummyRecBox.SetLineColor(kRed+1);
      current.cd(10);
      {
        TLegend legend(0.08,0.48,0.92,0.90);legend.SetBorderSize(0);
        legend.SetFillStyle(0);legend.SetTextSize(0.042);
        legend.AddEntry(&dummyGen,"GEMC GEN count; center: A_{T}=0","p");
        legend.AddEntry(&dummyGenBox,"GEN |#Delta A_{1}(JAM25 A_{T})| syst.","f");
        legend.AddEntry(&dummyRec,"GEMC REC count; center: A_{T}=0","p");
        legend.AddEntry(&dummyRecBox,"REC |#Delta A_{1}(JAM25 A_{T})| syst.","f");
        legend.AddEntry(&dummyJam,"JAM25 A_{1}","l");
        legend.DrawClone();
        TLatex note;note.SetNDC();note.SetTextSize(0.040);
        note.DrawLatex(0.08,0.34,(topo+" current GEMC statistics").c_str());
        note.DrawLatex(0.08,0.25,"Boxes show #pm full A_{T}-model shift");
        note.DrawLatex(0.08,0.16,"Overflow points/boxes are clipped at frame edge");
      }
      projected.cd(10);
      {
        TLegend legend(0.08,0.59,0.92,0.90);legend.SetBorderSize(0);
        legend.SetFillStyle(0);legend.SetTextSize(0.044);
        std::ostringstream projectedLabel;
        projectedLabel<<"Projected REC (R_{N}="<<std::setprecision(4)
                      <<projectionToCountRatio<<"); center: A_{T}=0";
        legend.AddEntry(&dummyProjected,projectedLabel.str().c_str(),"p");
        legend.AddEntry(&dummyRecBox,"|#Delta A_{1}(JAM25 A_{T})| syst.","f");
        legend.AddEntry(&dummyJam,"JAM25 A_{1}","l");
        legend.DrawClone();
        TLatex note;note.SetNDC();note.SetTextSize(0.040);
        note.DrawLatex(0.08,0.43,(topo+" target-statistics projection").c_str());
        note.DrawLatex(0.08,0.34,"Boxes show #pm full A_{T}-model shift");
        note.DrawLatex(0.08,0.25,"Points with #deltaA_{1}^{stat}>1 omitted");
        note.DrawLatex(0.08,0.16,"Overflow points/boxes are clipped at frame edge");
      }
      const std::string pdf=dir+"/He3DIS_A1_"+topo+".pdf";
      current.Print((pdf+"(").c_str());
      projected.Print((pdf+")").c_str());
      std::cout<<"[PlotA1ProjectionHe3DIS] wrote "<<pdf
               <<" (current GEMC + projected REC; R_N="
               <<projectionToCountRatio<<"; omitted projected points with "
               <<"deltaA1>1: "<<skippedProjectedLargeError<<")\n";
      current.Clear();projected.Clear();
    }
  }

  // Extract g1 directly from the measured longitudinal asymmetry with three
  // explicitly separated treatments of the unmeasured transverse response:
  //   1) SMC: A2=0;
  //   2) CLAS: use the JAM25 g2/g1 correction in the central value;
  //   3) current: A_perp=0, with the JAM25 transverse shift shown as a box.
  // All methods reuse the exact same finite-count A_parallel values, bin-mean
  // kinematics and projected statistical errors.  Thus differences between
  // the PDFs are caused only by the g1 extraction prescription.
  void PlotG1ProjectionHe3DIS(
      const std::function<double(bool,double,double,double,double,double,double)>&
          evaluateAparallel,
      const std::function<std::array<double,7>(bool,double,double,double)>&
          evaluateNucleon,
      double projectionToCountRatio=1.0) {
    if(he3DISModels.size()!=2 || !evaluateAparallel || !evaluateNucleon ||
       !(projectionToCountRatio>0.0) ||
       !std::isfinite(projectionToCountRatio)) {
      std::cerr<<"[PlotG1ProjectionHe3DIS] requires GEN, REC, the nucleon "
                   "model, and a positive projection/count ratio.\n";
      return;
    }
    if(he3DISProjectionValuesOverPN.empty() ||
       he3DISProjectionErrorsOverPN.empty()) {
      std::cerr<<"[PlotG1ProjectionHe3DIS] call PlotProjection first so the "
                   "same finite-count A_parallel values can be reused.\n";
      return;
    }
    const std::string dir=outputDir+"/g1_projection";
    fs::create_directories(dir);
    const std::string xDir=outputDir+"/xg1_projection";
    fs::create_directories(xDir);
    const double projectedErrorScale=1.0/std::sqrt(projectionToCountRatio);
    const double modelXMax=he3DISModelXMax;
    constexpr double nucleonMass=0.93892;
    struct Spec {std::string col,axis;int bins;double lo,hi;};
    using Profile=ROOT::RDF::RResultPtr<TProfile>;
    auto nodeFor=[](He3DISModel& m,const std::string& t)->ROOT::RDF::RNode {
      return t=="epp" ? *m.epp : *m.ed;
    };
    auto finite=[](const std::array<double,7>& v) {
      return std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);});
    };
    // Return D, eta and xi in the convention of Sec. 8 of arXiv:2106.08805.
    auto factors=[&](double xb,double q2,double y,double r) {
      std::array<double,3> out{{NAN,NAN,NAN}};
      if(!(xb>0.0) || !(q2>0.0) || !(y>0.0 && y<1.0) ||
         !std::isfinite(r)) return out;
      const double gamma=2.0*nucleonMass*xb/std::sqrt(q2);
      const double gamma2=gamma*gamma;
      const double common=4.0*(1.0-y)-gamma2*y*y;
      const double twoPlus=2.0+gamma2*y;
      const double dDen=2.0*(1.0+gamma2)*y*y+common*(1.0+r);
      if(!(dDen>0.0) || !(twoPlus>0.0) || !(2.0-y>0.0)) return out;
      out[0]=y*(2.0-y)*twoPlus/dDen;
      out[1]=gamma*common/((2.0-y)*twoPlus);
      out[2]=gamma*(2.0-y)/twoPlus;
      return out;
    };
    struct MethodSpec {
      int number;
      std::string tag;
      std::string title;
      std::string centerLabel;
      bool drawModelShift;
    };
    const std::array<MethodSpec,3> methods{{
      {1,"SMC_A2_zero","SMC: A_{2}=0","A_{2}=0",true},
      {2,"CLAS_g2_model","CLAS: JAM25 g_{2}/g_{1}",
       "JAM25 g_{2}/g_{1} correction",false},
      {3,"Aperp_zero","Current: A_{#perp}=0","A_{#perp}=0",true}
    }};
    for(const std::string topo:{"epp","ed"}) {
      std::vector<Spec> specs={{"xB","x_{B}",25,0,1},
        {"Q2","Q^{2} [GeV^{2}]",25,0,10},{"W","W [GeV]",25,1,5},
        {"y","y",25,0,1},{"nu","#nu [GeV]",25,0,11}};
      for(const auto& s:He3DISSpectatorPlotSpecs())
        specs.push_back({topo+"_"+s.variable,s.axis,std::max(1,s.bins/2),s.lo,s.hi});
      specs.push_back({topo+"_x_active","x_{N} = Q^{2}/(2p_{active}#upointq)",
                       he3DISXActiveBins,he3DISXActiveMin,
                       he3DISXActiveMax});
      std::array<std::vector<Profile>,2> meanXB,meanQ2,meanY,meanAlpha,
          meanActiveX,meanCosQ;
      for(int m=0;m<2;++m) {
        auto node=nodeFor(he3DISModels[m],topo);
        for(size_t i=0;i<specs.size();++i) {
          const auto& s=specs[i];
          const std::string base="g1_mean_"+topo+"_"+std::to_string(m)+"_"+
                                 std::to_string(i);
          auto profileModel=[&](const std::string& suffix) {
            return ROOT::RDF::TProfile1DModel((base+suffix).c_str(),"",
                                              s.bins,s.lo,s.hi);
          };
          meanXB[m].push_back(node.Profile1D(profileModel("_xb"),s.col,"xB"));
          meanQ2[m].push_back(node.Profile1D(profileModel("_q2"),s.col,"Q2"));
          meanY[m].push_back(node.Profile1D(profileModel("_y"),s.col,"y"));
          meanAlpha[m].push_back(node.Profile1D(
              profileModel("_alpha"),s.col,topo+"_alpha_active"));
          meanActiveX[m].push_back(node.Profile1D(
              profileModel("_xactive"),s.col,topo+"_x_active"));
          meanCosQ[m].push_back(node.Profile1D(
              profileModel("_cosq"),s.col,"model_cos_theta_q"));
        }
      }
      // The same cached bin-mean profiles produce both g1 and x*g1.  This
      // avoids a second RDataFrame event scan for the x-weighted PDF.
      for(int quantity=0;quantity<2;++quantity) {
      const bool plotXG1=quantity==1;
      const std::string quantityName=plotXG1 ? "xg1" : "g1";
      for(const auto& method:methods) {
      const std::string objectPrefix=quantityName+"_method"+
          std::to_string(method.number)+"_"+topo;
      TCanvas current((objectPrefix+"_current").c_str(),
                      ((plotXG1?"Current GEMC xg1 - ":"Current GEMC g1 - ")+
                       method.title).c_str(),
                      2200,1800);
      TCanvas projected((objectPrefix+"_projected").c_str(),
                        ((plotXG1?"Projected xg1 - ":"Projected g1 - ")+
                         method.title).c_str(),
                        2200,1800);
      current.Divide(4,3);projected.Divide(4,3);
      std::vector<std::unique_ptr<TH1D>> frames;
      std::vector<std::unique_ptr<TGraphErrors>> points;
      std::vector<std::unique_ptr<TGraph>> curves;
      int skippedProjectedLargeA1Error=0;
      for(size_t i=0;i<specs.size();++i) {
        const std::string key=topo+"/"+specs[i].col;
        const auto valueIt=he3DISProjectionValuesOverPN.find(key);
        const auto errorIt=he3DISProjectionErrorsOverPN.find(key);
        if(valueIt==he3DISProjectionValuesOverPN.end() ||
           errorIt==he3DISProjectionErrorsOverPN.end()) continue;
        auto genCurrent=std::make_unique<TGraphErrors>();
        auto recCurrent=std::make_unique<TGraphErrors>();
        auto recProjected=std::make_unique<TGraphErrors>();
        auto jamCurrent=std::make_unique<TGraph>();
        auto jamProjected=std::make_unique<TGraph>();
        struct BoxSpec {double x1,y1,x2,y2;int color;};
        std::vector<BoxSpec> currentBoxes,projectedBoxes;
        double extent=0.02;
        const double boxHalfWidth=0.35*(specs[i].hi-specs[i].lo)/specs[i].bins;
        for(int m=0;m<2;++m) {
          auto& xb=meanXB[m][i].GetValue();
          auto& q2=meanQ2[m][i].GetValue();
          auto& yy=meanY[m][i].GetValue();
          auto& aa=meanAlpha[m][i].GetValue();
          auto& xx=meanActiveX[m][i].GetValue();
          auto& cq=meanCosQ[m][i].GetValue();
          const auto& values=valueIt->second[m];
          const auto& errors=errorIt->second[m];
          for(int bin=1;bin<=specs[i].bins;++bin) {
            if(xb.GetBinEntries(bin)<=0 ||
               bin>=static_cast<int>(values.size()) ||
               bin>=static_cast<int>(errors.size()) ||
               !std::isfinite(values[bin]) ||
               !std::isfinite(errors[bin])) continue;
            const double xbMean=xb.GetBinContent(bin);
            const double q2Mean=q2.GetBinContent(bin);
            const double alphaMean=aa.GetBinContent(bin);
            if(!(alphaMean>0.0)) continue;
            const double activeX=xx.GetBinContent(bin);
            if(!(activeX>0.0 && activeX<modelXMax)) continue;
            const double modelAlpha=xbMean/activeX;
            const auto nucleon=evaluateNucleon(topo=="epp",xbMean,
                                               modelAlpha,q2Mean);
            if(!finite(nucleon)) continue;
            const auto k=factors(xbMean,q2Mean,yy.GetBinContent(bin),nucleon[4]);
            if(!std::isfinite(k[0]) || !std::isfinite(k[1]) ||
               !std::isfinite(k[2])) continue;
            const double mixing=1.0+k[1]*k[2];
            if(std::abs(k[0]*mixing)<1.e-12) continue;
            const double gamma=2.0*nucleonMass*xbMean/std::sqrt(q2Mean);
            const double gamma2=gamma*gamma;
            const double f1Factor=nucleon[0]/(1.0+gamma*gamma);
            const double quantityScale=plotXG1 ? activeX : 1.0;

            // Coefficient converting the measured nucleon A_parallel into g1.
            // Method 1 (SMC): set A2=0.
            // Method 2 (CLAS): retain a model g2/g1 in
            //   A_parallel = D/F1 * [(1+eta*gamma)g1
            //               +(eta*gamma-gamma^2)g2].
            // Method 3 (current): set A_perp=d(A2-xi*A1)=0.
            double longitudinalCoefficient=NAN;
            if(method.number==1) {
              longitudinalCoefficient=f1Factor/k[0];
            } else if(method.number==2) {
              // Algebraically this is F1/{D[(1+eta*gamma)+
              // (eta*gamma-gamma^2)g2/g1]}.  Multiplying numerator and
              // denominator by g1 avoids a numerical singularity where the
              // JAM25 g1 curve crosses zero.
              const double modelResponse=
                  (1.0+k[1]*gamma)*nucleon[2]+
                  (k[1]*gamma-gamma2)*nucleon[3];
              const double denominator=k[0]*modelResponse;
              if(!std::isfinite(denominator) ||
                 std::abs(denominator)<1.e-14) continue;
              longitudinalCoefficient=nucleon[0]*nucleon[2]/denominator;
            } else {
              longitudinalCoefficient=
                  f1Factor*(1.0+gamma*k[2])/(k[0]*mixing);
            }
            if(!std::isfinite(longitudinalCoefficient)) continue;

            // values[] already contains the nucleon A_parallel after the P_N
            // correction made in PlotProjection.  Do not divide by P_N again.
            const double g1Zero=
                quantityScale*longitudinalCoefficient*values[bin];
            const double g1Error=quantityScale*
                std::abs(longitudinalCoefficient)*errors[bin];
            // For the two zero-transverse-input methods, show the full JAM25
            // correction as a symmetric model-dependence box.  The CLAS
            // method already puts the JAM25 g2/g1 correction in its center.
            double modelShift=0.0;
            if(method.number==1) {
              modelShift=quantityScale*f1Factor*(gamma-k[1])*nucleon[6];
            } else if(method.number==3) {
              modelShift=quantityScale*f1Factor*(gamma-k[1])*
                  (nucleon[6]-k[2]*nucleon[5])/mixing;
            }
            const double systematic=std::abs(modelShift);
            const double x=xb.GetBinCenter(bin);
            auto* graph=m==0 ? genCurrent.get() : recCurrent.get();
            graph->SetPoint(graph->GetN(),x,g1Zero);
            graph->SetPointError(graph->GetN()-1,0,g1Error);
            if(method.drawModelShift)
              currentBoxes.push_back({x-boxHalfWidth,g1Zero-systematic,
                                      x+boxHalfWidth,g1Zero+systematic,
                                      m==0 ? kBlue+1 : kRed+1});

            if(m==1) {
              const double cosq=std::clamp(cq.GetBinContent(bin),-1.0,1.0);
              const double modelAparallel=evaluateAparallel(topo=="epp",
                  xbMean,modelAlpha,q2Mean,yy.GetBinContent(bin),
                  -cosq,std::sqrt(std::max(0.0,1.0-cosq*cosq)));
              if(std::isfinite(modelAparallel)) {
                const double projectedG1Zero=
                    quantityScale*longitudinalCoefficient*modelAparallel;
                const double projectedG1Error=g1Error*projectedErrorScale;
                // Use exactly the A1 plot's validity criterion so the two
                // projections contain the same bins.
                const double projectedA1Error=errors[bin]/
                    std::abs(k[0]*mixing)*projectedErrorScale;
                if(std::isfinite(projectedA1Error) &&
                   projectedA1Error<=1.0) {
                  recProjected->SetPoint(recProjected->GetN(),x,projectedG1Zero);
                  recProjected->SetPointError(recProjected->GetN()-1,0,
                                               projectedG1Error);
                  if(method.drawModelShift)
                    projectedBoxes.push_back({x-boxHalfWidth,
                        projectedG1Zero-systematic,x+boxHalfWidth,
                        projectedG1Zero+systematic,kRed+1});
                  if(std::abs(projectedG1Zero)<=1.0)
                    extent=std::max(extent,1.15*std::abs(projectedG1Zero));
                } else {
                  ++skippedProjectedLargeA1Error;
                }
              }
              jamProjected->SetPoint(jamProjected->GetN(),x,
                                     quantityScale*nucleon[2]);
            } else {
              jamCurrent->SetPoint(jamCurrent->GetN(),x,
                                   quantityScale*nucleon[2]);
            }
            if(std::abs(g1Zero)<=1.0)
              extent=std::max(extent,1.15*std::abs(g1Zero));
          }
        }
        for(auto* graph:{genCurrent.get(),recCurrent.get(),recProjected.get()}) {
          graph->SetLineStyle(1);graph->SetLineWidth(1);
        }
        genCurrent->SetMarkerStyle(20);genCurrent->SetMarkerSize(0.55);
        genCurrent->SetMarkerColor(kBlue+1);genCurrent->SetLineColor(kBlue+1);
        recCurrent->SetMarkerStyle(24);recCurrent->SetMarkerSize(0.55);
        recCurrent->SetMarkerColor(kRed+1);recCurrent->SetLineColor(kRed+1);
        recProjected->SetMarkerStyle(20);recProjected->SetMarkerSize(0.60);
        recProjected->SetMarkerColor(kRed+1);recProjected->SetLineColor(kRed+1);
        for(auto* curve:{jamCurrent.get(),jamProjected.get()}) {
          curve->SetLineColor(kBlack);curve->SetLineWidth(2);
        }
        for(int page=0;page<2;++page) {
          (page==0?current:projected).cd(static_cast<int>(i)+1);
          gPad->SetTicks();
          auto frame=std::make_unique<TH1D>(
              ("frame_"+objectPrefix+"_"+std::to_string(i)+"_"+
               std::to_string(page)).c_str(),
              (";"+specs[i].axis+";"+
               std::string(plotXG1?"x g_{1}^{":"g_{1}^{")+
               std::string(topo=="epp"?"n":"p")+"}").c_str(),
              specs[i].bins,specs[i].lo,specs[i].hi);
          frame->SetDirectory(nullptr);frame->SetStats(0);
          frame->SetMinimum(-extent);frame->SetMaximum(extent);frame->Draw("AXIS");
          TLine zero(specs[i].lo,0,specs[i].hi,0);
          zero.SetLineStyle(2);zero.SetLineColor(kGray+2);zero.DrawClone();
          const auto& boxes=page==0 ? currentBoxes : projectedBoxes;
          for(const auto& b:boxes) {
            TBox box(b.x1,b.y1,b.x2,b.y2);
            box.SetFillColorAlpha(b.color,0.22);
            box.SetLineColor(b.color);box.SetLineWidth(1);box.DrawClone("F SAME");
          }
          if(page==0) {
            if(jamCurrent->GetN()>1) jamCurrent->Draw("L SAME");
            if(genCurrent->GetN()>0) genCurrent->Draw("P SAME");
            if(recCurrent->GetN()>0) recCurrent->Draw("P SAME");
          } else {
            if(jamProjected->GetN()>1) jamProjected->Draw("L SAME");
            if(recProjected->GetN()>0) recProjected->Draw("P SAME");
          }
          frames.push_back(std::move(frame));
        }
        points.push_back(std::move(genCurrent));
        points.push_back(std::move(recCurrent));
        points.push_back(std::move(recProjected));
        curves.push_back(std::move(jamCurrent));
        curves.push_back(std::move(jamProjected));
      }

      TGraphErrors dummyGen,dummyRec,dummyProjected;
      dummyGen.SetMarkerStyle(20);dummyGen.SetMarkerColor(kBlue+1);
      dummyRec.SetMarkerStyle(24);dummyRec.SetMarkerColor(kRed+1);
      dummyProjected.SetMarkerStyle(20);dummyProjected.SetMarkerColor(kRed+1);
      TGraph dummyJam;dummyJam.SetLineColor(kBlack);dummyJam.SetLineWidth(2);
      TBox dummyGenBox,dummyRecBox;
      dummyGenBox.SetFillColorAlpha(kBlue+1,0.22);dummyGenBox.SetLineColor(kBlue+1);
      dummyRecBox.SetFillColorAlpha(kRed+1,0.22);dummyRecBox.SetLineColor(kRed+1);
      const std::string shiftSource=method.number==1 ? "JAM25 A_{2}" :
                                                    "JAM25 A_{#perp}";
      current.cd(10);
      {
        TLegend legend(0.08,0.48,0.92,0.90);legend.SetBorderSize(0);
        legend.SetFillStyle(0);legend.SetTextSize(0.042);
        legend.AddEntry(&dummyGen,("GEMC GEN count; "+
            method.centerLabel).c_str(),"p");
        if(method.drawModelShift)
          legend.AddEntry(&dummyGenBox,("GEN |#Delta "+
              std::string(plotXG1?"xg_{1}":"g_{1}")+"("+shiftSource+")|").c_str(),"f");
        legend.AddEntry(&dummyRec,("GEMC REC count; "+
            method.centerLabel).c_str(),"p");
        if(method.drawModelShift)
          legend.AddEntry(&dummyRecBox,("REC |#Delta "+
              std::string(plotXG1?"xg_{1}":"g_{1}")+"("+shiftSource+")|").c_str(),"f");
        legend.AddEntry(&dummyJam,plotXG1?"JAM25 xg_{1}":"JAM25 g_{1}","l");
        legend.DrawClone();
        TLatex note;note.SetNDC();note.SetTextSize(0.040);
        note.DrawLatex(0.08,0.34,(topo+" current GEMC statistics").c_str());
        note.DrawLatex(0.08,0.25,method.title.c_str());
      }
      projected.cd(10);
      {
        TLegend legend(0.08,0.59,0.92,0.90);legend.SetBorderSize(0);
        legend.SetFillStyle(0);legend.SetTextSize(0.044);
        std::ostringstream projectedLabel;
        projectedLabel<<"Projected REC (R_{N}="<<std::setprecision(4)
                      <<projectionToCountRatio<<"); "<<method.centerLabel;
        legend.AddEntry(&dummyProjected,projectedLabel.str().c_str(),"p");
        if(method.drawModelShift)
          legend.AddEntry(&dummyRecBox,("|#Delta "+
              std::string(plotXG1?"xg_{1}":"g_{1}")+"("+shiftSource+")|").c_str(),"f");
        legend.AddEntry(&dummyJam,plotXG1?"JAM25 xg_{1}":"JAM25 g_{1}","l");
        legend.DrawClone();
        TLatex note;note.SetNDC();note.SetTextSize(0.040);
        note.DrawLatex(0.08,0.43,(topo+" target-statistics projection").c_str());
        note.DrawLatex(0.08,0.34,method.title.c_str());
        note.DrawLatex(0.08,0.25,"Bins with #deltaA_{1}^{stat}>1 omitted");
      }
      const std::string pdf=(plotXG1?xDir:dir)+"/He3DIS_"+
          quantityName+"_method"+std::to_string(method.number)+"_"+
          method.tag+"_"+topo+".pdf";
      current.Print((pdf+"(").c_str());
      projected.Print((pdf+")").c_str());
      std::cout<<"[PlotG1ProjectionHe3DIS] wrote "<<pdf
               <<" using method "<<method.number<<" ("<<method.tag<<")"
               <<" (current GEMC + projected REC; R_N="
               <<projectionToCountRatio<<"; omitted projected bins with "
               <<"deltaA1>1: "<<skippedProjectedLargeA1Error<<")\n";
      current.Clear();projected.Clear();
      }
      }
    }
  }

  // Polarized-EMC projection for the ed two-body-breakup tag.  The active
  // nucleon is a proton.  JAMpol25 supplies the trusted free-proton g1 and its
  // replica uncertainty, while the external table supplies only the bound /
  // free modification ratio.  Statistical inputs are reused from the exact
  // REC ed A_parallel/PN count projection; Pb, Pt and PN must not be applied a
  // second time here.
  void PlotPolarizedEMCEd(
      const std::function<std::array<double,7>(bool,double,double,double)>&
          evaluateNucleon,
      const std::function<std::array<double,4>(bool,double,double,double)>&
          evaluateG1Uncertainty,
      const std::string& ratioTablePath,
      double projectionToCountRatio=1.0) {
    if(he3DISModels.size()!=2 || !evaluateNucleon || !evaluateG1Uncertainty ||
       !(projectionToCountRatio>0.0) ||
       !std::isfinite(projectionToCountRatio)) {
      std::cerr<<"[PlotPolarizedEMCEd] requires GEN/REC, JAM25 central and "
                   "replica evaluators, and a positive projection/count ratio.\n";
      return;
    }
    const std::string key="ed/ed_x_active";
    const auto errorIt=he3DISProjectionErrorsOverPN.find(key);
    if(errorIt==he3DISProjectionErrorsOverPN.end()) {
      std::cerr<<"[PlotPolarizedEMCEd] call PlotProjection first so the REC "
                   "ed A_parallel/PN statistical errors are available.\n";
      return;
    }

    struct RatioPoint {double x,ratio;};
    std::vector<RatioPoint> ratioTable;
    {
      std::ifstream input(ratioTablePath);
      if(!input) {
        std::cerr<<"[PlotPolarizedEMCEd] cannot open ratio table: "
                 <<ratioTablePath<<'\n';
        return;
      }
      std::string line;
      while(std::getline(input,line)) {
        const auto first=line.find_first_not_of(" \t");
        if(first==std::string::npos || line[first]=='#') continue;
        std::istringstream values(line);
        double x=0.,g1pFree=0.,g1pBound=0.,protonRatio=0.;
        double g1nFree=0.,g1nBound=0.,neutronRatio=0.;
        if(!(values>>x>>g1pFree>>g1pBound>>protonRatio
                   >>g1nFree>>g1nBound>>neutronRatio)) continue;
        if(std::isfinite(x) && std::isfinite(protonRatio) && x>0.0)
          ratioTable.push_back({x,protonRatio});
      }
    }
    std::sort(ratioTable.begin(),ratioTable.end(),
              [](const RatioPoint& a,const RatioPoint& b){return a.x<b.x;});
    ratioTable.erase(std::unique(ratioTable.begin(),ratioTable.end(),
        [](const RatioPoint& a,const RatioPoint& b){return a.x==b.x;}),
        ratioTable.end());
    if(ratioTable.size()<2) {
      std::cerr<<"[PlotPolarizedEMCEd] ratio table has fewer than two valid "
                   "proton points: "<<ratioTablePath<<'\n';
      return;
    }
    auto interpolateRatio=[&](double x) {
      if(!std::isfinite(x) || x<ratioTable.front().x ||
         x>ratioTable.back().x)
        return std::numeric_limits<double>::quiet_NaN();
      const auto upper=std::lower_bound(ratioTable.begin(),ratioTable.end(),x,
          [](const RatioPoint& point,double value){return point.x<value;});
      if(upper==ratioTable.begin()) return upper->ratio;
      if(upper==ratioTable.end()) return ratioTable.back().ratio;
      if(upper->x==x) return upper->ratio;
      const auto lower=std::prev(upper);
      const double fraction=(x-lower->x)/(upper->x-lower->x);
      return lower->ratio+fraction*(upper->ratio-lower->ratio);
    };

    const int bins=he3DISXActiveBins;
    const double xMin=he3DISXActiveMin;
    const double xMax=he3DISXActiveMax;
    const double modelXMax=he3DISModelXMax;
    constexpr double nucleonMass=0.93892;
    const double projectedErrorScale=1.0/std::sqrt(projectionToCountRatio);
    auto rec=*he3DISModels[1].ed;
    auto profileModel=[=](const char* name) {
      return ROOT::RDF::TProfile1DModel(name,"",bins,xMin,xMax);
    };
    // Book before GetValue so this is one additional REC ed event loop.
    auto meanXB=rec.Profile1D(profileModel("emc_ed_mean_xb"),
                              "ed_x_active","xB");
    auto meanQ2=rec.Profile1D(profileModel("emc_ed_mean_q2"),
                              "ed_x_active","Q2");
    auto meanY=rec.Profile1D(profileModel("emc_ed_mean_y"),
                             "ed_x_active","y");
    auto meanAlpha=rec.Profile1D(profileModel("emc_ed_mean_alpha"),
                                 "ed_x_active","ed_alpha_active");
    auto meanXN=rec.Profile1D(profileModel("emc_ed_mean_xn"),
                              "ed_x_active","ed_x_active");

    auto& xbProfile=meanXB.GetValue();
    auto& q2Profile=meanQ2.GetValue();
    auto& yProfile=meanY.GetValue();
    auto& alphaProfile=meanAlpha.GetValue();
    auto& xNProfile=meanXN.GetValue();
    const auto& recAparallelErrors=errorIt->second[1];
    auto finite7=[](const std::array<double,7>& values) {
      return std::all_of(values.begin(),values.end(),
                         [](double value){return std::isfinite(value);});
    };
    auto finite4=[](const std::array<double,4>& values) {
      return std::all_of(values.begin(),values.end(),
                         [](double value){return std::isfinite(value);});
    };
    auto factors=[=](double xb,double q2,double y,double r) {
      std::array<double,3> out{{NAN,NAN,NAN}};
      if(!(xb>0.0) || !(q2>0.0) || !(y>0.0 && y<1.0) ||
         !std::isfinite(r)) return out;
      const double gamma=2.0*nucleonMass*xb/std::sqrt(q2);
      const double gamma2=gamma*gamma;
      const double common=4.0*(1.0-y)-gamma2*y*y;
      const double twoPlus=2.0+gamma2*y;
      const double denominator=2.0*(1.0+gamma2)*y*y+common*(1.0+r);
      if(!(denominator>0.0) || !(twoPlus>0.0) || !(2.0-y>0.0))
        return out;
      out[0]=y*(2.0-y)*twoPlus/denominator;       // D
      out[1]=gamma*common/((2.0-y)*twoPlus);      // eta
      out[2]=gamma*(2.0-y)/twoPlus;               // xi
      return out;
    };

    struct BoxSpec {double x1,y1,x2,y2;int color;double alpha;};
    TGraph freeCurve,boundCurve,ratioCurve;
    TGraph alphaControl,xbControl,q2Control;
    TGraphAsymmErrors freeJamBand;
    std::array<TGraphErrors,3> boundProjection,totalBoundProjection;
    std::array<TGraphErrors,3> ratioProjection,totalRatioProjection;
    TGraphErrors statSensitivity,totalSensitivity;
    std::vector<BoxSpec> ratioJamBoxes;
    std::array<std::vector<BoxSpec>,3> boundMethodBoxes;
    std::array<std::vector<BoxSpec>,3> ratioMethodBoxes;
    std::array<double,3> g1Min{{0.0,0.0,0.0}},g1Max{{0.0,0.0,0.0}};
    double sensitivityMax=0.0,chi2Stat=0.0,chi2Total=0.0;
    int validBins=0;
    std::cout<<"[PlotPolarizedEMCEd] xN Q2 R_EMC stat Aperp JAM "
                 "S_stat S_total\n";
    for(int bin=1;bin<=bins;++bin) {
      if(xbProfile.GetBinEntries(bin)<=0 ||
         bin>=static_cast<int>(recAparallelErrors.size()) ||
         !std::isfinite(recAparallelErrors[bin])) continue;
      const double xb=xbProfile.GetBinContent(bin);
      const double q2=q2Profile.GetBinContent(bin);
      const double y=yProfile.GetBinContent(bin);
      const double alpha=alphaProfile.GetBinContent(bin);
      const double binLow=xbProfile.GetXaxis()->GetBinLowEdge(bin);
      const double binHigh=xbProfile.GetXaxis()->GetBinUpEdge(bin);
      if(!(alpha>0.0)) continue;
      const double xN=xNProfile.GetBinContent(bin);
      if(!(xN>0.0 && xN<modelXMax)) continue;
      const double ratio=interpolateRatio(xN);
      if(!std::isfinite(ratio)) continue;
      // The external model API forms x_N=xB/alpha internally but also uses xB
      // separately in finite-Q2 factors.  Keep the measured electron xB and
      // adapt only the API alpha argument so xB/modelAlpha=xN_exact.
      const double modelAlpha=xb/xN;
      const auto nucleon=evaluateNucleon(false,xb,modelAlpha,q2);
      const auto jam=evaluateG1Uncertainty(false,xb,modelAlpha,q2);
      if(!finite7(nucleon) || !finite4(jam)) continue;
      const double g1Free=jam[0];
      if(!std::isfinite(g1Free) || std::abs(g1Free)<1.e-12) continue;
      const auto k=factors(xb,q2,y,nucleon[4]);
      if(!std::isfinite(k[0]) || !std::isfinite(k[1]) ||
         !std::isfinite(k[2])) continue;
      const double mixing=1.0+k[1]*k[2];
      if(std::abs(k[0]*mixing)<1.e-12) continue;
      const double gamma=2.0*nucleonMass*xb/std::sqrt(q2);
      const double gamma2=gamma*gamma;
      const double f1Factor=nucleon[0]/(1.0+gamma*gamma);
      const double modelResponse=(1.0+k[1]*gamma)*nucleon[2]+
          (k[1]*gamma-gamma2)*nucleon[3];
      const double clasDenominator=k[0]*modelResponse;
      if(!std::isfinite(clasDenominator) ||
         std::abs(clasDenominator)<1.e-14) continue;
      // The three g1 extraction prescriptions used by
      // PlotG1ProjectionHe3DIS: SMC A2=0, CLAS JAM25 g2/g1, and
      // the current A_perp=0 method.
      const std::array<double,3> coefficients{{
          f1Factor/k[0],
          nucleon[0]*nucleon[2]/clasDenominator,
          f1Factor*(1.0+gamma*k[2])/(k[0]*mixing)}};
      if(!std::all_of(coefficients.begin(),coefficients.end(),
                      [](double value){return std::isfinite(value);})) continue;
      const double g1Bound=ratio*g1Free;
      std::array<double,3> statG1,transverseG1{{
          std::abs(f1Factor*(gamma-k[1])*nucleon[6]),
          0.0,
          std::abs(f1Factor*(gamma-k[1])*
              (nucleon[6]-k[2]*nucleon[5])/mixing)}};
      for(size_t method=0;method<coefficients.size();++method)
        statG1[method]=std::abs(coefficients[method])*
            recAparallelErrors[bin]*projectedErrorScale;
      // Keep the ratio/sensitivity summary on the current A_perp=0 method.
      const double currentStatG1=statG1[2];
      const double aperpG1=transverseG1[2];
      const double statRatio=currentStatG1/std::abs(g1Free);
      const double aperpRatio=aperpG1/std::abs(g1Free);
      const double jamRatioPlus=std::abs(ratio)*jam[1]/std::abs(g1Free);
      const double jamRatioMinus=std::abs(ratio)*jam[2]/std::abs(g1Free);
      const double jamRatioSymm=std::abs(ratio)*jam[3]/std::abs(g1Free);
      if(!std::isfinite(statRatio) || !std::isfinite(aperpRatio) ||
         !std::isfinite(jamRatioSymm)) continue;

      freeCurve.SetPoint(freeCurve.GetN(),xN,xN*g1Free);
      boundCurve.SetPoint(boundCurve.GetN(),xN,xN*g1Bound);
      ratioCurve.SetPoint(ratioCurve.GetN(),xN,ratio);
      alphaControl.SetPoint(alphaControl.GetN(),xN,alpha);
      xbControl.SetPoint(xbControl.GetN(),xN,xb);
      q2Control.SetPoint(q2Control.GetN(),xN,q2);
      for(size_t method=0;method<boundProjection.size();++method) {
        auto& graph=boundProjection[method];
        graph.SetPoint(graph.GetN(),xN,xN*g1Bound);
        graph.SetPointError(graph.GetN()-1,0.0,xN*statG1[method]);
        const double boundJamG1=std::abs(ratio)*jam[3];
        const double totalG1=std::sqrt(statG1[method]*statG1[method]+
            transverseG1[method]*transverseG1[method]+
            boundJamG1*boundJamG1);
        auto& totalGraph=totalBoundProjection[method];
        totalGraph.SetPoint(totalGraph.GetN(),xN,xN*g1Bound);
        totalGraph.SetPointError(totalGraph.GetN()-1,0.0,xN*totalG1);
        const double methodStatRatio=statG1[method]/std::abs(g1Free);
        const double methodTransverseRatio=
            transverseG1[method]/std::abs(g1Free);
        const double totalRatioError=std::sqrt(
            methodStatRatio*methodStatRatio+
            methodTransverseRatio*methodTransverseRatio+
            jamRatioSymm*jamRatioSymm);
        auto& ratioGraph=ratioProjection[method];
        ratioGraph.SetPoint(ratioGraph.GetN(),xN,ratio);
        ratioGraph.SetPointError(ratioGraph.GetN()-1,0.0,methodStatRatio);
        auto& totalRatioGraph=totalRatioProjection[method];
        totalRatioGraph.SetPoint(totalRatioGraph.GetN(),xN,ratio);
        totalRatioGraph.SetPointError(
            totalRatioGraph.GetN()-1,0.0,totalRatioError);
        if(method!=1)
          boundMethodBoxes[method].push_back({
              binLow,xN*(g1Bound-transverseG1[method]),
              binHigh,xN*(g1Bound+transverseG1[method]),
              kAzure+1,0.24});
        if(method!=1)
          ratioMethodBoxes[method].push_back({
              binLow,ratio-methodTransverseRatio,
              binHigh,ratio+methodTransverseRatio,
              kAzure+1,0.24});
        g1Min[method]=std::min({g1Min[method],xN*(g1Free-jam[2]),
                               xN*(g1Bound-totalG1)});
        g1Max[method]=std::max({g1Max[method],xN*(g1Free+jam[1]),
                               xN*(g1Bound+totalG1)});
      }
      freeJamBand.SetPoint(freeJamBand.GetN(),xN,xN*g1Free);
      freeJamBand.SetPointError(freeJamBand.GetN()-1,0.0,0.0,
                                xN*jam[2],xN*jam[1]);
      // The ratio JAM25 baseline uncertainty belongs to each projected bin;
      // keep it as a per-bin box rather than a continuous model band.
      ratioJamBoxes.push_back({binLow,ratio-jamRatioPlus,
                               binHigh,ratio+jamRatioMinus,
                               kGray+1,0.28});
      const double difference=std::abs(ratio-1.0);
      const double totalError=std::sqrt(statRatio*statRatio+
          aperpRatio*aperpRatio+jamRatioSymm*jamRatioSymm);
      const double statS=statRatio>0.0 ? difference/statRatio : 0.0;
      const double totalS=totalError>0.0 ? difference/totalError : 0.0;
      statSensitivity.SetPoint(statSensitivity.GetN(),xN,statS);
      totalSensitivity.SetPoint(totalSensitivity.GetN(),xN,totalS);
      sensitivityMax=std::max({sensitivityMax,statS,totalS});
      chi2Stat+=statS*statS;chi2Total+=totalS*totalS;
      ++validBins;
      std::cout<<std::setprecision(6)<<"  "<<xN<<' '<<q2<<' '<<ratio<<' '
               <<statRatio<<' '<<aperpRatio<<' '<<jamRatioSymm<<' '
               <<statS<<' '<<totalS<<'\n';
    }
    if(validBins==0) {
      std::cerr<<"[PlotPolarizedEMCEd] no valid REC ed bins.\n";
      return;
    }

    // The projected data retain their coarse statistical bins, but the
    // JAM25 and bound/free model objects are evaluated on an independent
    // dense x_N grid.  Interpolate only the bin-mean alpha and Q2 needed to
    // map x_N back to the model inputs; g1 and its replica uncertainty are
    // reevaluated at every dense point.
    if(alphaControl.GetN()>=2 && xbControl.GetN()>=2 &&
       q2Control.GetN()>=2) {
      alphaControl.Sort();xbControl.Sort();q2Control.Sort();
      double denseXMin=0.0,denseXMax=0.0,unused=0.0;
      alphaControl.GetPoint(0,denseXMin,unused);
      alphaControl.GetPoint(alphaControl.GetN()-1,denseXMax,unused);
      constexpr int denseModelPoints=201;
      TGraph denseFreeCurve,denseBoundCurve,denseRatioCurve;
      TGraphAsymmErrors denseFreeJamBand;
      for(int point=0;point<denseModelPoints;++point) {
        const double fraction=static_cast<double>(point)/
                              static_cast<double>(denseModelPoints-1);
        const double xN=denseXMin+fraction*(denseXMax-denseXMin);
        const double alpha=alphaControl.Eval(xN);
        const double xb=xbControl.Eval(xN);
        const double q2=q2Control.Eval(xN);
        if(!(alpha>0.0) || !(xb>0.0) || !(q2>0.0)) continue;
        const double ratio=interpolateRatio(xN);
        const double modelAlpha=xb/xN;
        const auto jam=evaluateG1Uncertainty(false,xb,modelAlpha,q2);
        if(!std::isfinite(ratio) || !finite4(jam) ||
           std::abs(jam[0])<1.e-12) continue;
        const double g1Free=jam[0];
        const double g1Bound=ratio*g1Free;
        denseFreeCurve.SetPoint(denseFreeCurve.GetN(),xN,xN*g1Free);
        denseBoundCurve.SetPoint(denseBoundCurve.GetN(),xN,xN*g1Bound);
        denseRatioCurve.SetPoint(denseRatioCurve.GetN(),xN,ratio);
        denseFreeJamBand.SetPoint(
            denseFreeJamBand.GetN(),xN,xN*g1Free);
        denseFreeJamBand.SetPointError(
            denseFreeJamBand.GetN()-1,0.0,0.0,xN*jam[2],xN*jam[1]);
      }
      if(denseFreeCurve.GetN()>=2) {
        freeCurve=denseFreeCurve;
        boundCurve=denseBoundCurve;
        ratioCurve=denseRatioCurve;
        freeJamBand=denseFreeJamBand;
      }
    }

    auto style=[](TGraph& graph,int color,int style,int width) {
      graph.SetLineColor(color);graph.SetLineStyle(style);graph.SetLineWidth(width);
    };
    // Filled TGraph bands require monotonically ordered x coordinates.
    freeJamBand.Sort();
    style(freeCurve,kBlack,1,3);style(boundCurve,kBlue+1,2,3);
    style(ratioCurve,kBlue+1,2,3);
    freeJamBand.SetFillColorAlpha(kGray+1,0.28);
    freeJamBand.SetLineColor(kGray+1);
    for(auto& graph:boundProjection) {
      graph.SetMarkerStyle(20);graph.SetMarkerSize(0.8);
      graph.SetMarkerColor(kRed+1);graph.SetLineColor(kRed+1);
    }
    for(auto& graph:totalBoundProjection) {
      graph.SetMarkerStyle(20);graph.SetMarkerSize(0.8);
      graph.SetMarkerColor(kMagenta+2);graph.SetLineColor(kMagenta+2);
      graph.SetLineWidth(2);
    }
    for(auto& graph:ratioProjection) {
      graph.SetMarkerStyle(20);graph.SetMarkerSize(0.8);
      graph.SetMarkerColor(kRed+1);graph.SetLineColor(kRed+1);
    }
    for(auto& graph:totalRatioProjection) {
      graph.SetMarkerStyle(20);graph.SetMarkerSize(0.8);
      graph.SetMarkerColor(kMagenta+2);graph.SetLineColor(kMagenta+2);
      graph.SetLineWidth(2);
    }
    statSensitivity.SetMarkerStyle(20);statSensitivity.SetMarkerColor(kBlue+1);
    statSensitivity.SetLineColor(kBlue+1);
    totalSensitivity.SetMarkerStyle(24);totalSensitivity.SetMarkerColor(kRed+1);
    totalSensitivity.SetLineColor(kRed+1);
    const std::string dir=outputDir+"/polarized_EMC";
    fs::create_directories(dir);
    const std::string pdf=dir+"/He3DIS_polarized_EMC_ed.pdf";
    auto drawBoxes=[](const std::vector<BoxSpec>& boxes) {
      for(const auto& b:boxes) {
        TBox box(b.x1,b.y1,b.x2,b.y2);
        box.SetFillColorAlpha(b.color,b.alpha);box.SetLineColor(b.color);
        box.DrawClone("F SAME");
      }
    };

    TBox dummyGray,dummyBlue;
    dummyGray.SetFillColorAlpha(kGray+1,0.28);
    dummyGray.SetLineColor(kGray+1);
    dummyBlue.SetFillColorAlpha(kAzure+1,0.24);
    const std::array<std::string,3> methodTitles{{
        "Method 1: SMC, A_{2}=0",
        "Method 2: CLAS, JAM25 g_{2}/g_{1}",
        "Method 3: Current, A_{#perp}=0"}};
    const std::array<std::string,3> systematicLabels{{
        "JAM25 A_{2} extraction systematic","",
        "JAM25 A_{#perp} extraction systematic"}};
    for(size_t method=0;method<boundProjection.size();++method) {
      const std::string suffix=std::to_string(method+1);
      TCanvas g1Canvas(("c_polarized_emc_ed_xg1_method"+suffix).c_str(),
                       ("ed polarized EMC xg1 "+methodTitles[method]).c_str(),
                       1100,850);
      const double padding=0.12*std::max(1.e-6,g1Max[method]-g1Min[method]);
      TH1D g1Frame(("h_polarized_emc_ed_xg1_method"+suffix).c_str(),
          (methodTitles[method]+";x_{N}=x_{B}/#alpha_{active};x_{N}g_{1}^{p}").c_str(),
          bins,xMin,modelXMax);
      g1Frame.SetStats(0);g1Frame.SetMinimum(g1Min[method]-padding);
      g1Frame.SetMaximum(g1Max[method]+padding);g1Frame.Draw("AXIS");
      TLine zeroG1(xMin,0.0,modelXMax,0.0);zeroG1.SetLineStyle(2);
      zeroG1.SetLineColor(kGray+2);zeroG1.DrawClone();
      freeJamBand.Draw("3 SAME");drawBoxes(boundMethodBoxes[method]);
      freeCurve.Draw("L SAME");boundCurve.Draw("L SAME");
      boundProjection[method].Draw("P SAME");
      TLegend g1Legend(0.13,0.13,0.56,0.40);g1Legend.SetBorderSize(0);
      g1Legend.SetFillStyle(0);g1Legend.AddEntry(&freeCurve,"JAM25 free p","l");
      g1Legend.AddEntry(&freeJamBand,"JAM25 replica uncertainty","f");
      g1Legend.AddEntry(&boundCurve,"^{3}He bound/free ratio #times JAM25","l");
      g1Legend.AddEntry(&boundProjection[method],"Projected REC ed (stat.)","p");
      if(method!=1)
        g1Legend.AddEntry(&dummyBlue,systematicLabels[method].c_str(),"f");
      g1Legend.DrawClone();
      g1Canvas.Print((pdf+(method==0?"(":"")).c_str());

      TCanvas totalCanvas(
          ("c_polarized_emc_ed_xg1_method"+suffix+"_total").c_str(),
          ("ed polarized EMC xg1 total uncertainty "+methodTitles[method]).c_str(),
          1100,850);
      TH1D totalFrame(
          ("h_polarized_emc_ed_xg1_method"+suffix+"_total").c_str(),
          (methodTitles[method]+" (combined uncertainty);x_{N}=x_{B}/#alpha_{active};x_{N}g_{1}^{p}").c_str(),
          bins,xMin,modelXMax);
      totalFrame.SetStats(0);totalFrame.SetMinimum(g1Min[method]-padding);
      totalFrame.SetMaximum(g1Max[method]+padding);totalFrame.Draw("AXIS");
      zeroG1.DrawClone();
      freeJamBand.Draw("3 SAME");
      freeCurve.Draw("L SAME");boundCurve.Draw("L SAME");
      totalBoundProjection[method].Draw("P SAME");
      TLegend totalLegend(0.13,0.13,0.56,0.34);
      totalLegend.SetBorderSize(0);totalLegend.SetFillStyle(0);
      totalLegend.AddEntry(&freeCurve,"JAM25 free p","l");
      totalLegend.AddEntry(&freeJamBand,"JAM25 replica uncertainty","f");
      totalLegend.AddEntry(&boundCurve,"^{3}He bound/free ratio #times JAM25","l");
      totalLegend.AddEntry(&totalBoundProjection[method],
          "Projected REC ed (total uncertainty)","p");
      totalLegend.DrawClone();
      TLatex totalNote;totalNote.SetNDC();totalNote.SetTextSize(0.034);
      totalNote.DrawLatex(0.14,0.92,
          "Single bar: stat. #oplus transverse extraction #oplus JAM25");
      totalCanvas.Print(pdf.c_str());
    }

    for(size_t method=0;method<ratioProjection.size();++method) {
      const std::string suffix=std::to_string(method+1);
      TCanvas ratioCanvas(
          ("c_polarized_emc_ed_ratio_method"+suffix).c_str(),
          ("ed polarized EMC ratio "+methodTitles[method]).c_str(),1100,850);
      TH1D ratioFrame(
          ("h_polarized_emc_ed_ratio_method"+suffix).c_str(),
          (methodTitles[method]+";x_{N}=x_{B}/#alpha_{active};R_{EMC}^{pol,p}=g_{1}^{p,*}/g_{1,JAM25}^{p}").c_str(),
          bins,xMin,modelXMax);
      ratioFrame.SetStats(0);ratioFrame.SetMinimum(0.8);
      ratioFrame.SetMaximum(1.2);ratioFrame.Draw("AXIS");
      TLine unity(xMin,1.0,modelXMax,1.0);unity.SetLineStyle(2);
      unity.SetLineColor(kBlack);unity.DrawClone();
      drawBoxes(ratioJamBoxes);drawBoxes(ratioMethodBoxes[method]);
      ratioCurve.Draw("L SAME");ratioProjection[method].Draw("P SAME");
      TLegend ratioLegend(0.13,0.13,0.56,0.40);
      ratioLegend.SetBorderSize(0);ratioLegend.SetFillStyle(0);
      ratioLegend.AddEntry(&unity,"No polarized EMC: R=1","l");
      ratioLegend.AddEntry(&ratioCurve,"^{3}He bound/free model","l");
      ratioLegend.AddEntry(&ratioProjection[method],"Projected REC ed (stat.)","p");
      if(method!=1)
        ratioLegend.AddEntry(&dummyBlue,systematicLabels[method].c_str(),"f");
      ratioLegend.AddEntry(&dummyGray,"JAM25 free-baseline uncertainty (per bin)","f");
      ratioLegend.DrawClone();
      ratioCanvas.Print(pdf.c_str());

      TCanvas totalRatioCanvas(
          ("c_polarized_emc_ed_ratio_method"+suffix+"_total").c_str(),
          ("ed polarized EMC ratio total uncertainty "+methodTitles[method]).c_str(),
          1100,850);
      TH1D totalRatioFrame(
          ("h_polarized_emc_ed_ratio_method"+suffix+"_total").c_str(),
          (methodTitles[method]+" (combined uncertainty);x_{N}=x_{B}/#alpha_{active};R_{EMC}^{pol,p}").c_str(),
          bins,xMin,modelXMax);
      totalRatioFrame.SetStats(0);
      totalRatioFrame.SetMinimum(0.8);
      totalRatioFrame.SetMaximum(1.2);
      totalRatioFrame.Draw("AXIS");unity.DrawClone();
      ratioCurve.Draw("L SAME");totalRatioProjection[method].Draw("P SAME");
      TLegend totalRatioLegend(0.13,0.13,0.56,0.34);
      totalRatioLegend.SetBorderSize(0);totalRatioLegend.SetFillStyle(0);
      totalRatioLegend.AddEntry(&unity,"No polarized EMC: R=1","l");
      totalRatioLegend.AddEntry(&ratioCurve,"^{3}He bound/free model","l");
      totalRatioLegend.AddEntry(&totalRatioProjection[method],
          "Projected REC ed (total uncertainty)","p");
      totalRatioLegend.DrawClone();
      TLatex totalRatioNote;totalRatioNote.SetNDC();
      totalRatioNote.SetTextSize(0.034);
      totalRatioNote.DrawLatex(0.14,0.92,
          "Single bar: stat. #oplus transverse extraction #oplus JAM25");
      totalRatioNote.DrawLatex(0.14,0.87,
          ("g_{1} extraction: "+methodTitles[method]).c_str());
      totalRatioCanvas.Print(pdf.c_str());
    }

    TCanvas sensitivityCanvas("c_polarized_emc_ed_sensitivity",
                              "ed polarized EMC sensitivity",1100,850);
    TH1D sensitivityFrame("h_polarized_emc_ed_sensitivity",
        ";x_{N}=x_{B}/#alpha_{active};Expected |R_{EMC}^{pol,p}-1|/#delta R",
        bins,xMin,modelXMax);
    sensitivityFrame.SetStats(0);sensitivityFrame.SetMinimum(0.0);
    sensitivityFrame.SetMaximum(std::max(1.0,1.18*sensitivityMax));
    sensitivityFrame.Draw("AXIS");
    statSensitivity.Draw("PL SAME");totalSensitivity.Draw("PL SAME");
    TLegend sensitivityLegend(0.49,0.73,0.89,0.89);
    sensitivityLegend.SetBorderSize(0);sensitivityLegend.SetFillStyle(0);
    sensitivityLegend.AddEntry(&statSensitivity,"Statistical only","lp");
    sensitivityLegend.AddEntry(&totalSensitivity,
        "Stat. #oplus A_{T} syst. #oplus JAM25","lp");
    sensitivityLegend.DrawClone();
    TLatex sensitivityNote;sensitivityNote.SetNDC();sensitivityNote.SetTextSize(0.035);
    std::ostringstream note;
    note<<"Independent-bin estimate: Z_{stat}="<<std::setprecision(3)
        <<std::sqrt(chi2Stat)<<", Z_{total}="<<std::sqrt(chi2Total);
    sensitivityNote.DrawLatex(0.14,0.92,note.str().c_str());
    sensitivityCanvas.Print((pdf+")").c_str());
    std::cout<<"[PlotPolarizedEMCEd] wrote "<<pdf<<" with "<<validBins
             <<" valid bins; independent-bin Z_stat="<<std::sqrt(chi2Stat)
             <<", Z_total="<<std::sqrt(chi2Total)<<'\n';
  }

  void PlotTopologyMigrationHe3DIS(ROOT::RDF::RNode rec) {
    fs::create_directories(outputDir);
    const std::array<std::string,3> gen = {
        "GEN_pass_epp", "GEN_pass_epn", "GEN_pass_ed"};
    const std::array<std::string,3> reco = {
        "REC_Event_pass_epp", "REC_Event_pass_epn", "REC_Event_pass_ed"};
    const std::array<const char*,3> labels = {"epp", "epn", "ed"};
    TH2D matrix("h_he3dis_topology_migration",
                "He3DIS topology flags;GEN topology;REC topology;Events",
                3, 0, 3, 3, 0, 3);
    matrix.SetDirectory(nullptr); matrix.SetStats(0);
    for (int i = 0; i < 3; ++i) {
      matrix.GetXaxis()->SetBinLabel(i+1, labels[i]);
      matrix.GetYaxis()->SetBinLabel(i+1, labels[i]);
      for (int j = 0; j < 3; ++j) {
        auto count = rec.Filter(gen[i] + " && " + reco[j]).Count();
        matrix.SetBinContent(i+1, j+1, count.GetValue());
      }
    }
    TCanvas canvas("c_he3dis_topology_migration",
                   "He3DIS topology migration", 800, 700);
    matrix.Draw("COLZ TEXT");
    const std::string output =
        outputDir + "/He3DIS_topology_migration.pdf";
    canvas.SaveAs(output.c_str());
    std::cout << "[PlotTopologyMigrationHe3DIS] wrote " << output << '\n';
  }

  void PlotParticleKinematicHe3DVCS() {
    if (he3DVCSModels.empty()) {
      std::cerr << "[PlotParticleKinematicHe3DVCS] no He3DVCS model was added.\n";
      return;
    }

    struct ParticleSpec {
      std::string key;
      std::string title;
      std::string p;
      std::string theta;
      std::string phi;
      double pMin;
      double pMax;
      double thetaMax;
    };
    const std::array<ParticleSpec, 3> particles = {{
        {"electron", "e'", "recel_p", "recel_theta", "recel_phi",
         0.0, 6.0, 180.0},
        {"photon", "#gamma", "recpho_p", "recpho_theta", "recpho_phi",
         4.0, 11.0, 10.0},
        {"he3", "^{3}He'", "reche3_p", "reche3_theta", "reche3_phi",
         0.0, 0.6, 100.0}
    }};

    fs::create_directories(outputDir);
    for (size_t model = 0; model < he3DVCSModels.size(); ++model) {
      struct HistSet {
        ROOT::RDF::RResultPtr<TH1D> p;
        ROOT::RDF::RResultPtr<TH1D> theta;
        ROOT::RDF::RResultPtr<TH1D> phi;
        ROOT::RDF::RResultPtr<TH2D> pVsTheta;
      };
      std::vector<HistSet> histograms;
      for (const auto& particle : particles) {
        const std::string prefix = "h_he3dvcs_m" + std::to_string(model) +
                                   "_" + particle.key;
        histograms.push_back({
            he3DVCSModels[model]->Histo1D(
                {(prefix + "_p").c_str(),
                 (particle.title + ";p [GeV];Events").c_str(),
                 120, particle.pMin, particle.pMax}, particle.p),
            he3DVCSModels[model]->Histo1D(
                {(prefix + "_theta").c_str(),
                 (particle.title + ";#theta [deg];Events").c_str(),
                 120, 0.0, particle.thetaMax}, particle.theta),
            he3DVCSModels[model]->Histo1D(
                {(prefix + "_phi").c_str(),
                 (particle.title + ";#phi [deg];Events").c_str(),
                 120, 0.0, 360.0}, particle.phi),
            he3DVCSModels[model]->Histo2D(
                {(prefix + "_p_vs_theta").c_str(),
                 (particle.title +
                  ";#theta [deg];p [GeV]").c_str(),
                 100, 0.0, particle.thetaMax,
                 100, particle.pMin, particle.pMax},
                particle.theta, particle.p)});
      }

      TCanvas canvas(("c_he3dvcs_particles_m" + std::to_string(model)).c_str(),
                     "He3 DVCS particle kinematics", 2000, 1200);
      canvas.Divide(4, 3);
      for (size_t row = 0; row < histograms.size(); ++row) {
        canvas.cd(static_cast<int>(4 * row + 1));
        histograms[row].p->Draw("hist");
        canvas.cd(static_cast<int>(4 * row + 2));
        histograms[row].theta->Draw("hist");
        canvas.cd(static_cast<int>(4 * row + 3));
        histograms[row].phi->Draw("hist");
        canvas.cd(static_cast<int>(4 * row + 4));
        gPad->SetLogz();
        histograms[row].pVsTheta->Draw("COLZ");
      }
      std::string output = outputDir + "/He3DVCS_particle_kinematics";
      if (he3DVCSModels.size() > 1) output += "_model" + std::to_string(model);
      output += ".pdf";
      canvas.SaveAs(output.c_str());
      std::cout << "[PlotParticleKinematicHe3DVCS] wrote " << output << std::endl;
    }
  }

  void PlotDVCSKinematicHe3DVCS() {
    if (he3DVCSModels.empty()) {
      std::cerr << "[PlotDVCSKinematicHe3DVCS] no He3DVCS model was added.\n";
      return;
    }
    fs::create_directories(outputDir);
    for (size_t model = 0; model < he3DVCSModels.size(); ++model) {
      const std::string prefix = "h_he3dvcs_kin_m" + std::to_string(model);
      auto h_xb = he3DVCSModels[model]->Histo1D(
          {(prefix + "_xb").c_str(), ";x_{B};Events", 120, 0.0, 0.6}, "xB");
      auto h_q2 = he3DVCSModels[model]->Histo1D(
          {(prefix + "_q2").c_str(), ";Q^{2} [GeV^{2}];Events", 120, 0.0, 8.0}, "Q2");
      auto h_w = he3DVCSModels[model]->Histo1D(
          {(prefix + "_w").c_str(), ";W [GeV];Events", 120, 0.0, 5.0}, "W");
      auto h_y = he3DVCSModels[model]->Histo1D(
          {(prefix + "_y").c_str(), ";y;Events", 120, 0.0, 1.0}, "y");
      auto h_nu = he3DVCSModels[model]->Histo1D(
          {(prefix + "_nu").c_str(), ";#nu [GeV];Events", 120, 0.0, 11.0}, "nu");
      auto h_t = he3DVCSModels[model]->Histo1D(
          {(prefix + "_t").c_str(), ";-t [GeV^{2}];Events", 120, 0.0, 1.0}, "minus_t");
      auto h_phi = he3DVCSModels[model]->Histo1D(
          {(prefix + "_phi").c_str(), ";#phi [deg];Events", 120, 0.0, 360.0}, "phi");
      auto h_xbq2 = he3DVCSModels[model]->Histo2D(
          {(prefix + "_xbq2").c_str(), ";x_{B};Q^{2} [GeV^{2}]",
           100, 0.0, 0.6, 100, 0.0, 8.0}, "xB", "Q2");

      TCanvas canvas(("c_he3dvcs_kin_m" + std::to_string(model)).c_str(),
                     "He3 DVCS kinematics", 1600, 900);
      canvas.Divide(4, 2);
      canvas.cd(1); gPad->SetLogy(); h_xb->Draw("hist");
      canvas.cd(2); gPad->SetLogy(); h_q2->Draw("hist");
      canvas.cd(3); gPad->SetLogy(); h_w->Draw("hist");
      canvas.cd(4); gPad->SetLogy(); h_y->Draw("hist");
      canvas.cd(5); gPad->SetLogy(); h_nu->Draw("hist");
      canvas.cd(6); gPad->SetLogy(); h_t->Draw("hist");
      canvas.cd(7); h_phi->Draw("hist");
      canvas.cd(8); gPad->SetLogz(); h_xbq2->Draw("colz");
      std::string output = outputDir + "/He3DVCS_DVCS_kinematics";
      if (he3DVCSModels.size() > 1) output += "_model" + std::to_string(model);
      output += ".pdf";
      canvas.SaveAs(output.c_str());
      std::cout << "[PlotDVCSKinematicHe3DVCS] wrote " << output << std::endl;
    }
  }

  void PlotExclusivityHe3DVCS() {
    if (he3DVCSModels.empty()) {
      std::cerr << "[PlotExclusivityHe3DVCS] no He3DVCS model was added.\n";
      return;
    }
    fs::create_directories(outputDir);
    for (size_t model = 0; model < he3DVCSModels.size(); ++model) {
      const std::string prefix = "h_he3dvcs_excl_m" + std::to_string(model);
      auto h_mm2_ehe3 = he3DVCSModels[model]->Histo1D(
          {(prefix + "_mm2_ehe3").c_str(),
           ";MM^{2}_{e^{3}He} [GeV^{2}];Events",
           200, -0.2, 0.2}, "mm2_ehe3");
      auto h_mm2_ehe3gamma = he3DVCSModels[model]->Histo1D(
          {(prefix + "_mm2_ehe3gamma").c_str(),
           ";MM^{2}_{e^{3}He#gamma} [GeV^{2}];Events",
           200, -0.1, 0.1}, "mm2_ehe3gamma");
      auto h_mm2_egamma = he3DVCSModels[model]->Histo1D(
          {(prefix + "_mm2_egamma").c_str(),
           ";MM^{2}_{e#gamma} [GeV^{2}];Events",
           200, 7.0, 9.0}, "mm2_egamma");
      auto h_emiss = he3DVCSModels[model]->Histo1D(
          {(prefix + "_emiss").c_str(),
           ";E_{miss}(e^{3}He#gamma) [GeV];Events",
           200, -0.2, 0.2}, "missing_energy");
      auto h_ptmiss = he3DVCSModels[model]->Histo1D(
          {(prefix + "_ptmiss").c_str(),
           ";P^{miss}_{T}(e^{3}He#gamma) [GeV];Events",
           200, 0.0, 0.2}, "missing_pt");
      auto h_theta_gg = he3DVCSModels[model]->Histo1D(
          {(prefix + "_theta_gg").c_str(),
           ";#theta_{#gamma#gamma'} [deg];Events",
           200, 0.0, 10.0}, "theta_gamma_gamma_pred");
      auto h_theta_eg = he3DVCSModels[model]->Histo1D(
          {(prefix + "_theta_eg").c_str(),
           ";#theta_{e#gamma} [deg];Events",
           180, 0.0, 180.0}, "theta_e_gamma");
      auto h_dphi = he3DVCSModels[model]->Histo1D(
          {(prefix + "_dphi").c_str(),
           ";#Delta#phi [deg];Events",
           180, 0.0, 30.0}, "coplanarity");

      TCanvas canvas(("c_he3dvcs_excl_m" + std::to_string(model)).c_str(),
                     "He3 DVCS exclusivity", 1800, 900);
      canvas.Divide(4, 2);
      canvas.cd(1); h_mm2_ehe3->Draw("hist");
      canvas.cd(2); h_mm2_ehe3gamma->Draw("hist");
      canvas.cd(3); h_mm2_egamma->Draw("hist");
      canvas.cd(4); h_emiss->Draw("hist");
      canvas.cd(5); h_ptmiss->Draw("hist");
      canvas.cd(6); h_theta_gg->Draw("hist");
      canvas.cd(7); h_theta_eg->Draw("hist");
      canvas.cd(8); h_dphi->Draw("hist");
      std::string output = outputDir + "/He3DVCS_exclusivity";
      if (he3DVCSModels.size() > 1) output += "_model" + std::to_string(model);
      output += ".pdf";
      canvas.SaveAs(output.c_str());
      std::cout << "[PlotExclusivityHe3DVCS] wrote " << output << std::endl;
    }
  }

  void PlotPolarizationHe3DVCS(
      const std::vector<double>& xBbins = {},
      const std::vector<double>& Q2bins = {},
      const std::vector<double>& minusTbins = {}) {
    if (he3DVCSModels.empty()) {
      std::cerr << "[PlotPolarizationHe3DVCS] no He3DVCS model was added.\n";
      return;
    }
    const bool makeBinnedPlots =
        !xBbins.empty() || !Q2bins.empty() || !minusTbins.empty();
    auto validEdges = [](const std::vector<double>& edges) {
      return edges.size() >= 2 &&
             std::adjacent_find(edges.begin(), edges.end(),
                                std::greater_equal<double>()) == edges.end();
    };
    if (makeBinnedPlots &&
        (!validEdges(xBbins) || !validEdges(Q2bins) ||
         !validEdges(minusTbins))) {
      throw std::invalid_argument(
          "PlotPolarizationHe3DVCS requires strictly increasing xB, Q2, "
          "and -t bin-edge vectors, each with at least two entries");
    }
    fs::create_directories(outputDir);
    for (size_t model = 0; model < he3DVCSModels.size(); ++model) {
      auto& df = *he3DVCSModels[model];
      const std::string prefix = "h_he3dvcs_pol_m" + std::to_string(model);
      auto h_beam = df.Histo1D({(prefix + "_beam").c_str(),
          ";Beam helicity;Events", 3, -1.5, 1.5}, "MC_Event_helicity");
      auto h_target = df.Histo1D({(prefix + "_target").c_str(),
          ";Target polarization;Events", 3, -1.5, 1.5}, "MC_Event_targetPolarization");
      auto h_spin = df.Histo1D({(prefix + "_spin").c_str(),
          ";Beam #times target spin;Events", 3, -1.5, 1.5}, "beam_target_spin");
      auto h_beam_pos = df.Filter("MC_Event_helicity == 1").Histo1D(
          {(prefix + "_phi_beam_pos").c_str(), ";#phi [deg];Events", 24, 0.0, 360.0}, "phi");
      auto h_beam_neg = df.Filter("MC_Event_helicity == -1").Histo1D(
          {(prefix + "_phi_beam_neg").c_str(), ";#phi [deg];Events", 24, 0.0, 360.0}, "phi");
      auto h_target_pos = df.Filter("MC_Event_targetPolarization == 1").Histo1D(
          {(prefix + "_phi_target_pos").c_str(), ";#phi [deg];Events", 24, 0.0, 360.0}, "phi");
      auto h_target_neg = df.Filter("MC_Event_targetPolarization == -1").Histo1D(
          {(prefix + "_phi_target_neg").c_str(), ";#phi [deg];Events", 24, 0.0, 360.0}, "phi");
      auto h_parallel = df.Filter("beam_target_spin == 1").Histo1D(
          {(prefix + "_phi_parallel").c_str(), ";#phi [deg];Events", 24, 0.0, 360.0}, "phi");
      auto h_antiparallel = df.Filter("beam_target_spin == -1").Histo1D(
          {(prefix + "_phi_antiparallel").c_str(), ";#phi [deg];Events", 24, 0.0, 360.0}, "phi");

      // For A=(N_1-N_2)/(N_1+N_2), N_1 and N_2 are disjoint Poisson
      // counts.  The correlated numerator/denominator propagation gives
      // sigma_A=sqrt((1-A^2)/(N_1+N_2)); TH1::Divide is not appropriate.
      auto makeAsymmetry = [](const TH1D& first, const TH1D& second,
                              const std::string& name,
                              const std::string& yTitle,
                              const std::string& label) {
        auto result = std::unique_ptr<TH1D>(
            static_cast<TH1D*>(first.Clone(name.c_str())));
        result->Reset("ICES");
        result->SetDirectory(nullptr);
        result->SetTitle((";#phi [deg];" + yTitle).c_str());
        result->SetMinimum(-1.0);
        result->SetMaximum(1.0);

        std::cout << "[He3DVCS " << label
                  << "] bin  phi_center  N1  N2  asymmetry  stat_error\n";
        for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
          const double n1 = first.GetBinContent(bin);
          const double n2 = second.GetBinContent(bin);
          const double total = n1 + n2;
          double asymmetry = 0.0;
          double error = 0.0;
          if (total > 0.0) {
            asymmetry = (n1 - n2) / total;
            error = std::sqrt(std::max(0.0, 1.0 - asymmetry * asymmetry) /
                              total);
          }
          result->SetBinContent(bin, asymmetry);
          result->SetBinError(bin, error);
          std::cout << "  " << bin << "  " << result->GetBinCenter(bin)
                    << "  " << n1 << "  " << n2 << "  " << asymmetry
                    << "  " << error << '\n';
        }
        return result;
      };

      struct BinnedActions {
        size_t ix;
        size_t iq;
        size_t it;
        ROOT::RDF::RResultPtr<TH1D> beamPos;
        ROOT::RDF::RResultPtr<TH1D> beamNeg;
        ROOT::RDF::RResultPtr<TH1D> targetPos;
        ROOT::RDF::RResultPtr<TH1D> targetNeg;
        ROOT::RDF::RResultPtr<TH1D> parallel;
        ROOT::RDF::RResultPtr<TH1D> antiparallel;
      };
      std::vector<BinnedActions> binnedActions;

      // Book the inclusive and every binned action before the first
      // GetValue().  ROOT then fills everything in one event loop.
      if (makeBinnedPlots) {
        for (size_t ix = 0; ix + 1 < xBbins.size(); ++ix) {
          for (size_t iq = 0; iq + 1 < Q2bins.size(); ++iq) {
            for (size_t it = 0; it + 1 < minusTbins.size(); ++it) {
              const double xLo = xBbins[ix], xHi = xBbins[ix + 1];
              const double qLo = Q2bins[iq], qHi = Q2bins[iq + 1];
              const double tLo = minusTbins[it], tHi = minusTbins[it + 1];
              const bool lastX = ix + 2 == xBbins.size();
              const bool lastQ = iq + 2 == Q2bins.size();
              const bool lastT = it + 2 == minusTbins.size();
              const std::string binLabel =
                  "He3DVCS xB[" + std::to_string(xLo) + "," +
                  std::to_string(xHi) + "] Q2[" + std::to_string(qLo) +
                  "," + std::to_string(qHi) + "] -t[" +
                  std::to_string(tLo) + "," + std::to_string(tHi) + "]";
              auto binDf = df.Filter(
                  [=](double xb, double q2, double minusT) {
                    const bool passX =
                        xb >= xLo && (xb < xHi || (lastX && xb <= xHi));
                    const bool passQ =
                        q2 >= qLo && (q2 < qHi || (lastQ && q2 <= qHi));
                    const bool passT =
                        minusT >= tLo &&
                        (minusT < tHi || (lastT && minusT <= tHi));
                    return passX && passQ && passT;
                  },
                  {"xB", "Q2", "minus_t"}, binLabel);
              const std::string binPrefix =
                  prefix + "_xb" + std::to_string(ix) + "_q2" +
                  std::to_string(iq) + "_t" + std::to_string(it);
              auto phiModel = [&](const std::string& suffix) {
                return ROOT::RDF::TH1DModel(
                    (binPrefix + suffix).c_str(), ";#phi [deg];Events",
                    24, 0.0, 360.0);
              };
              binnedActions.push_back({
                  ix, iq, it,
                  binDf.Filter("MC_Event_helicity == 1").Histo1D(
                      phiModel("_beam_pos"), "phi"),
                  binDf.Filter("MC_Event_helicity == -1").Histo1D(
                      phiModel("_beam_neg"), "phi"),
                  binDf.Filter("MC_Event_targetPolarization == 1").Histo1D(
                      phiModel("_target_pos"), "phi"),
                  binDf.Filter("MC_Event_targetPolarization == -1").Histo1D(
                      phiModel("_target_neg"), "phi"),
                  binDf.Filter("beam_target_spin == 1").Histo1D(
                      phiModel("_parallel"), "phi"),
                  binDf.Filter("beam_target_spin == -1").Histo1D(
                      phiModel("_antiparallel"), "phi")});
            }
          }
        }
      }

      // This first dereference triggers the single common RDataFrame loop.
      const TH1D& beam_pos = h_beam_pos.GetValue();
      const TH1D& beam_neg = h_beam_neg.GetValue();
      const TH1D& target_pos = h_target_pos.GetValue();
      const TH1D& target_neg = h_target_neg.GetValue();
      const TH1D& parallel = h_parallel.GetValue();
      const TH1D& antiparallel = h_antiparallel.GetValue();

      auto alu = makeAsymmetry(beam_pos, beam_neg, prefix + "_alu",
                               "A_{LU}", "A_LU: beam + versus beam -");
      auto aul = makeAsymmetry(target_pos, target_neg, prefix + "_aul",
                               "A_{UL}", "A_UL: target + versus target -");
      auto all = makeAsymmetry(parallel, antiparallel, prefix + "_all",
                               "A_{LL}", "A_LL: parallel versus antiparallel");

      TCanvas canvas(("c_he3dvcs_pol_m" + std::to_string(model)).c_str(),
                     "He3 DVCS polarization", 1800, 900);
      canvas.Divide(3, 2);
      canvas.cd(1); h_beam->Draw("hist");
      canvas.cd(2); h_target->Draw("hist");
      canvas.cd(3); h_spin->Draw("hist");
      canvas.cd(4); alu->Draw("E1");
      canvas.cd(5); aul->Draw("E1");
      canvas.cd(6); all->Draw("E1");
      std::string output = outputDir + "/He3DVCS_polarization";
      if (he3DVCSModels.size() > 1) output += "_model" + std::to_string(model);
      output += ".pdf";
      canvas.SaveAs(output.c_str());
      std::cout << "[PlotPolarizationHe3DVCS] wrote " << output << std::endl;

      if (!makeBinnedPlots) continue;

      for (auto& actions : binnedActions) {
        const std::string binKey =
            "xb" + std::to_string(actions.ix) + "_q2" +
            std::to_string(actions.iq) + "_t" +
            std::to_string(actions.it);
        const std::string ranges =
            "x_{B} #in [" + std::to_string(xBbins[actions.ix]) + ", " +
            std::to_string(xBbins[actions.ix + 1]) + "], Q^{2} #in [" +
            std::to_string(Q2bins[actions.iq]) + ", " +
            std::to_string(Q2bins[actions.iq + 1]) + "], -t #in [" +
            std::to_string(minusTbins[actions.it]) + ", " +
            std::to_string(minusTbins[actions.it + 1]) + "]";
        auto binAlu = makeAsymmetry(
            actions.beamPos.GetValue(), actions.beamNeg.GetValue(),
            prefix + "_" + binKey + "_alu", "A_{LU}",
            binKey + " A_LU");
        auto binAul = makeAsymmetry(
            actions.targetPos.GetValue(), actions.targetNeg.GetValue(),
            prefix + "_" + binKey + "_aul", "A_{UL}",
            binKey + " A_UL");
        auto binAll = makeAsymmetry(
            actions.parallel.GetValue(), actions.antiparallel.GetValue(),
            prefix + "_" + binKey + "_all", "A_{LL}",
            binKey + " A_LL");
        binAlu->SetTitle((ranges + ";#phi [deg];A_{LU}").c_str());
        binAul->SetTitle((ranges + ";#phi [deg];A_{UL}").c_str());
        binAll->SetTitle((ranges + ";#phi [deg];A_{LL}").c_str());

        TCanvas binCanvas(("c_he3dvcs_pol_" + binKey + "_m" +
                           std::to_string(model)).c_str(),
                          "Binned He3 DVCS polarization", 1800, 550);
        binCanvas.Divide(3, 1);
        binCanvas.cd(1); binAlu->Draw("E1");
        binCanvas.cd(2); binAul->Draw("E1");
        binCanvas.cd(3); binAll->Draw("E1");
        std::string binOutput = outputDir + "/He3DVCS_polarization_" + binKey;
        if (he3DVCSModels.size() > 1)
          binOutput += "_model" + std::to_string(model);
        binOutput += ".pdf";
        binCanvas.SaveAs(binOutput.c_str());
        std::cout << "[PlotPolarizationHe3DVCS] wrote " << binOutput
                  << std::endl;
      }
    }
  }

  // Draw p, theta, phi, and p-vs-theta for every final-state particle in each of the
  // three He3DIS topologies. With one registered He3 model this produces
  // exactly three PDF files: epp, ed, and epn.
  void PlotParticleKinematicHe3() {
    if (he3Models.empty()) {
      std::cerr << "[PlotParticleKinematicHe3] no He3 MC model was added.\n";
      return;
    }

    struct ParticleSpec {
      int pid;
      int occurrence;
      std::string key;
      std::string title;
      double momentumMax;
      double thetaMax;
    };
    struct TopologySpec {
      int finalState;
      std::string key;
      std::string title;
      std::vector<ParticleSpec> particles;
    };

    const std::vector<TopologySpec> topologies = {
        {1, "epp", "e' + p + p",
         {{11, 0, "electron", "e'", 11.0, 60.0},
          {2212, 0, "proton1", "p_{1}", 0.8, 180.0},
          {2212, 1, "proton2", "p_{2}", 0.8, 180.0}}},
        {2, "ed", "e' + d",
         {{11, 0, "electron", "e'", 11.0, 60.0},
          {1000010020, 0, "deuteron", "d", 0.8, 180.0}}},
        {3, "epn", "e' + p + n",
         {{11, 0, "electron", "e'", 11.0, 60.0},
          {2212, 0, "proton", "p", 0.8, 180.0},
          {2112, 0, "neutron", "n", 0.8, 180.0}}}
    };

    auto particleKinematic =
        [](int targetPid, int occurrence, int observable) {
          return [targetPid, occurrence, observable](
                     const ROOT::VecOps::RVec<int>& pid,
                     const ROOT::VecOps::RVec<float>& px,
                     const ROOT::VecOps::RVec<float>& py,
                     const ROOT::VecOps::RVec<float>& pz) {
            const size_t size = std::min(
                {pid.size(), px.size(), py.size(), pz.size()});
            int found = 0;
            for (size_t i = 0; i < size; ++i) {
              if (pid[i] != targetPid) continue;
              if (found++ != occurrence) continue;

              const double momentum = std::sqrt(
                  static_cast<double>(px[i]) * px[i] +
                  static_cast<double>(py[i]) * py[i] +
                  static_cast<double>(pz[i]) * pz[i]);
              if (observable == 0) return momentum;
              if (observable == 1) {
                if (momentum == 0.0) return 0.0;
                const double cosine = std::clamp(
                    static_cast<double>(pz[i]) / momentum, -1.0, 1.0);
                return std::acos(cosine) * 180.0 / M_PI;
              }
              double phi = std::atan2(py[i], px[i]);
              if (phi < 0.0) phi += 2.0 * M_PI;
              return phi * 180.0 / M_PI;
            }
            return -999.0;
          };
        };

    fs::create_directories(outputDir);

    for (size_t modelIndex = 0; modelIndex < he3Models.size(); ++modelIndex) {
      for (const auto& topology : topologies) {
        ROOT::RDF::RNode df = he3Models[modelIndex]->Filter(
            [state = topology.finalState](int eventState) {
              return eventState == state;
            },
            {"he3_final_state"},
            "He3 topology " + topology.key);

        struct ParticleHistograms {
          ROOT::RDF::RResultPtr<TH1D> momentum;
          ROOT::RDF::RResultPtr<TH1D> theta;
          ROOT::RDF::RResultPtr<TH1D> phi;
          ROOT::RDF::RResultPtr<TH2D> pVsTheta;
        };
        std::vector<ParticleHistograms> histograms;
        histograms.reserve(topology.particles.size());

        for (const auto& particle : topology.particles) {
          const std::string prefix =
              "he3_m" + std::to_string(modelIndex) + "_" + topology.key +
              "_" + particle.key;
          const std::string pColumn = prefix + "_p";
          const std::string thetaColumn = prefix + "_theta";
          const std::string phiColumn = prefix + "_phi";

          df = df
              .Define(pColumn,
                      particleKinematic(particle.pid, particle.occurrence, 0),
                      {"MC_Particle_pid", "MC_Particle_px",
                       "MC_Particle_py", "MC_Particle_pz"})
              .Define(thetaColumn,
                      particleKinematic(particle.pid, particle.occurrence, 1),
                      {"MC_Particle_pid", "MC_Particle_px",
                       "MC_Particle_py", "MC_Particle_pz"})
              .Define(phiColumn,
                      particleKinematic(particle.pid, particle.occurrence, 2),
                      {"MC_Particle_pid", "MC_Particle_px",
                       "MC_Particle_py", "MC_Particle_pz"});

          const std::string pName = "h_" + prefix + "_p";
          const std::string thetaName = "h_" + prefix + "_theta";
          const std::string phiName = "h_" + prefix + "_phi";
          const std::string pTitle = topology.title + ";p(" +
                                     particle.title + ") [GeV];Events";
          const std::string thetaTitle = topology.title + ";#theta(" +
                                         particle.title + ") [deg];Events";
          const std::string phiTitle = topology.title + ";#phi(" +
                                       particle.title + ") [deg];Events";
          const std::string pVsThetaName =
              "h_" + prefix + "_p_vs_theta";
          const std::string pVsThetaTitle =
              topology.title + ";#theta(" + particle.title +
              ") [deg];p(" + particle.title + ") [GeV]";

          histograms.push_back({
              df.Histo1D({pName.c_str(), pTitle.c_str(), 120, 0.0,
                          particle.momentumMax},
                         pColumn),
              df.Histo1D({thetaName.c_str(), thetaTitle.c_str(), 120, 0.0,
                          particle.thetaMax},
                         thetaColumn),
              df.Histo1D({phiName.c_str(), phiTitle.c_str(), 120, 0.0,
                          360.0},
                         phiColumn),
              df.Histo2D({pVsThetaName.c_str(), pVsThetaTitle.c_str(),
                          100, 0.0, particle.thetaMax,
                          100, 0.0, particle.momentumMax},
                         thetaColumn, pColumn)});
        }

        const int rows = static_cast<int>(topology.particles.size());
        const std::string canvasName =
            "c_he3_m" + std::to_string(modelIndex) + "_" + topology.key;
        TCanvas canvas(canvasName.c_str(), topology.title.c_str(), 2000,
                       400 * rows);
        canvas.Divide(4, rows);

        for (size_t row = 0; row < histograms.size(); ++row) {
          canvas.cd(static_cast<int>(4 * row + 1));
          histograms[row].momentum->Draw("hist");
          canvas.cd(static_cast<int>(4 * row + 2));
          histograms[row].theta->Draw("hist");
          canvas.cd(static_cast<int>(4 * row + 3));
          histograms[row].phi->Draw("hist");
          canvas.cd(static_cast<int>(4 * row + 4));
          gPad->SetLogz();
          histograms[row].pVsTheta->Draw("COLZ");
        }

        std::string pdfName = outputDir + "/He3_" + topology.key +
                              "_particle_kinematics";
        if (he3Models.size() > 1) {
          pdfName += "_model" + std::to_string(modelIndex);
        }
        pdfName += ".pdf";
        canvas.SaveAs(pdfName.c_str());
        std::cout << "[PlotParticleKinematicHe3] wrote " << pdfName
                  << std::endl;
      }
    }
  }

  // Compare inclusive DIS kinematics for the three visible He3 final states.
  // Rows are e'pp, e'd, and e'pn; columns are xB, Q2, W, and xB vs Q2.
  void PlotDISKinematicHe3() {
    if (he3Models.empty()) {
      std::cerr << "[PlotDISKinematicHe3] no He3 MC model was added.\n";
      return;
    }

    struct TopologySpec {
      int finalState;
      std::string key;
      std::string title;
    };
    const std::array<TopologySpec, 3> topologies = {{
        {1, "epp", "e' + p + p"},
        {2, "ed", "e' + d"},
        {3, "epn", "e' + p + n"}
    }};

    fs::create_directories(outputDir);

    for (size_t modelIndex = 0; modelIndex < he3Models.size(); ++modelIndex) {
      struct DISHistograms {
        ROOT::RDF::RResultPtr<TH1D> xB;
        ROOT::RDF::RResultPtr<TH1D> Q2;
        ROOT::RDF::RResultPtr<TH1D> W;
        ROOT::RDF::RResultPtr<TH2D> xBvsQ2;
      };
      std::vector<DISHistograms> histograms;
      histograms.reserve(topologies.size());

      for (const auto& topology : topologies) {
        auto df = he3Models[modelIndex]->Filter(
            [state = topology.finalState](int eventState) {
              return eventState == state;
            },
            {"he3_final_state"},
            "He3 DIS topology " + topology.key);

        const std::string prefix =
            "he3_dis_m" + std::to_string(modelIndex) + "_" + topology.key;
        const std::string xBName = "h_" + prefix + "_xB";
        const std::string Q2Name = "h_" + prefix + "_Q2";
        const std::string WName = "h_" + prefix + "_W";
        const std::string xBvsQ2Name = "h_" + prefix + "_xB_vs_Q2";
        const std::string xBTitle =
            topology.title + ";x_{B};Events";
        const std::string Q2Title =
            topology.title + ";Q^{2} [GeV^{2}];Events";
        const std::string WTitle =
            topology.title + ";W [GeV];Events";
        const std::string xBvsQ2Title =
            topology.title + ";x_{B};Q^{2} [GeV^{2}]";

        histograms.push_back({
            df.Histo1D({xBName.c_str(), xBTitle.c_str(), 120, 0.0, 1.0},
                       "xB"),
            df.Histo1D({Q2Name.c_str(), Q2Title.c_str(), 120, 0.0, 12.0},
                       "Q2"),
            df.Histo1D({WName.c_str(), WTitle.c_str(), 120, 0.0, 5.0},
                       "W"),
            df.Histo2D({xBvsQ2Name.c_str(), xBvsQ2Title.c_str(), 100,
                        0.0, 1.0, 100, 0.0, 12.0},
                       "xB", "Q2")});
      }

      const std::string canvasName =
          "c_he3_dis_m" + std::to_string(modelIndex);
      TCanvas canvas(canvasName.c_str(), "He3 DIS kinematics", 1900, 1200);
      canvas.Divide(4, 3);

      for (size_t row = 0; row < histograms.size(); ++row) {
        canvas.cd(static_cast<int>(4 * row + 1));
        histograms[row].xB->Draw("hist");
        canvas.cd(static_cast<int>(4 * row + 2));
        histograms[row].Q2->Draw("hist");
        canvas.cd(static_cast<int>(4 * row + 3));
        histograms[row].W->Draw("hist");
        canvas.cd(static_cast<int>(4 * row + 4));
        histograms[row].xBvsQ2->Draw("colz");
      }

      std::string pdfName = outputDir + "/He3_DIS_kinematics";
      if (he3Models.size() > 1) {
        pdfName += "_model" + std::to_string(modelIndex);
      }
      pdfName += ".pdf";
      canvas.SaveAs(pdfName.c_str());
      std::cout << "[PlotDISKinematicHe3] wrote " << pdfName << std::endl;
    }
  }

  void AddModelPhi(ROOT::RDF::RNode df_data, const std::string& label, double beamEnergy, double luminosity) {
    auto plotter = std::make_unique<DISANAplotter>(PhiModeTag{}, df_data, beamEnergy, luminosity);
    labels.push_back(label);
    plotter->GeneratePhiKinematicHistos("el");
    plotter->GeneratePhiKinematicHistos("pro");
    plotter->GeneratePhiKinematicHistos("kMinus");
    plotter->GeneratePhiKinematicHistos("kPlus");
    plotters.push_back(std::move(plotter));
  }

  void AddModelPhi(ROOT::RDF::RNode df_data, const std::string& label, double beamEnergy, double luminosity, ROOT::RDF::RNode df_gen, ROOT::RDF::RNode df_rec,
                   ROOT::RDF::RNode df_radRatio, bool fAcc = false, bool fEff = false, bool fRad = false) {
    auto plotter = std::make_unique<DISANAplotter>(PhiModeTag{}, df_data, beamEnergy, luminosity, df_gen, df_rec, std::nullopt, std::nullopt, df_radRatio);
    plotter->GeneratePhiKinematicHistos("el");
    plotter->GeneratePhiKinematicHistos("pro");
    plotter->GeneratePhiKinematicHistos("kMinus");
    plotter->GeneratePhiKinematicHistos("kPlus");
    plotter->SetPlotApplyAcceptanceCorrection(fAcc);
    plotter->SetPlotApplyEfficiencyCorrection(fEff);
    plotter->SetPlotApplyRadiativeCorrection(fRad);

    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  void AddModelPhi(ROOT::RDF::RNode df_data, const std::string& label, double beamEnergy, double luminosity, ROOT::RDF::RNode df_gen, ROOT::RDF::RNode df_rec,
                   ROOT::RDF::RNode df_bkg, ROOT::RDF::RNode df_nobkg, ROOT::RDF::RNode df_radRatio, bool fAcc = false, bool fEff = false, bool fRad = false) {
    auto plotter = std::make_unique<DISANAplotter>(PhiModeTag{}, df_data, beamEnergy, luminosity, df_gen, df_rec, df_bkg, df_nobkg, df_radRatio);
    plotter->GeneratePhiKinematicHistos("el");
    plotter->GeneratePhiKinematicHistos("pro");
    plotter->GeneratePhiKinematicHistos("kMinus");
    plotter->GeneratePhiKinematicHistos("kPlus");
    plotter->SetPlotApplyAcceptanceCorrection(fAcc);
    plotter->SetPlotApplyEfficiencyCorrection(fEff);
    plotter->SetPlotApplyRadiativeCorrection(fRad);

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

  // Write one row per configured (xB,Q2,|t|) bin.  Called automatically when
  // a DVCS model is added, so every DISANA_Xplotter2csv*.cpp entry point that
  // uses this comparer gets the diagnostic without per-macro changes.
  void WriteDVCSPhysicalVolumeCSV(double beamEnergy,
                                  int nXBSteps = 240,
                                  int nQ2Steps = 240) {
    if (fXbins.UsesSDHEP()) {
      std::cout << "[Veff] SDHEP bins use (xi,cos(theta),t), not (xB,Q2,t); "
                   "skipping the xB-Q2-t physical-volume CSV.\n";
      return;
    }

    const long long energyKey = std::llround(beamEnergy * 1.0e6);
    if (physicalVolumeCSVWritten_[energyKey]) return;

    const BinManager::GeneratorPhaseSpaceCuts cuts{};
    std::cout << "[Veff] Generator phase-space cuts used in the integral: "
              << cuts.minXB << " < xB < " << cuts.maxXB << ", "
              << cuts.minQ2 << " < Q2 < " << cuts.maxQ2 << " GeV^2, "
              << cuts.minY << " < y < " << cuts.maxY << ", "
              << "W2 > " << cuts.minW2 << " GeV^2, "
              << cuts.minAbsT << " < |t| < " << cuts.maxAbsT << " GeV^2\n"
              << "[Veff] Analysis cuts p_e' > 2 GeV and p_gamma > 2 GeV are NOT "
                 "included in Veff; they are left to the acceptance correction.\n";

    std::ostringstream energyTag;
    energyTag << std::fixed << std::setprecision(3) << beamEnergy;
    std::string tag = energyTag.str();
    std::replace(tag.begin(), tag.end(), '.', 'p');
    const std::string filename = outputDir + "/effective_bin_volume_Eb" + tag + ".csv";

    std::ofstream csv(filename);
    if (!csv.is_open()) {
      std::cerr << "[Veff] Cannot open " << filename << " for writing\n";
      return;
    }

    csv << "mode,bin_index,ixB,iQ2,it,xB_min,xB_max,Q2_min,Q2_max,"
           "t_min,t_max,beam_energy,nominal_volume,effective_volume,"
           "Veff_over_Vnominal,gen_xB_min,gen_xB_max,gen_Q2_min,gen_Q2_max,"
           "gen_y_min,gen_y_max,gen_W2_min,gen_t_min,gen_t_max,"
           "n_xB_steps,n_Q2_steps\n";
    csv << std::setprecision(12);

    size_t binIndex = 0;
    auto writeRow = [&](const char *mode, long long ix, long long iq, long long it,
                        const BinManager::KinematicBin &bin) {
      const auto result = BinManager::ComputeDVCSPhysicalVolume(
          bin, beamEnergy, cuts, nXBSteps, nQ2Steps);
      csv << mode << ',' << binIndex++ << ',' << ix << ',' << iq << ',' << it << ','
          << bin.xBMin << ',' << bin.xBMax << ','
          << bin.Q2Min << ',' << bin.Q2Max << ','
          << bin.tMin << ',' << bin.tMax << ',' << beamEnergy << ','
          << result.nominalVolume << ',' << result.effectiveVolume << ','
          << result.fraction << ',' << cuts.minXB << ',' << cuts.maxXB << ','
          << cuts.minQ2 << ',' << cuts.maxQ2 << ','
          << cuts.minY << ',' << cuts.maxY << ',' << cuts.minW2 << ','
          << cuts.minAbsT << ',' << cuts.maxAbsT << ','
          << nXBSteps << ',' << nQ2Steps << '\n';
    };

    if (fXbins.HasCustomBins()) {
      const auto &bins = fXbins.GetCustomBins();
      for (size_t ib = 0; ib < bins.size(); ++ib) {
        writeRow("custom", static_cast<long long>(ib), -1, -1, bins[ib]);
      }
    } else {
      const auto &q2Edges = fXbins.GetQ2Bins();
      const auto &tEdges = fXbins.GetTBins();
      for (size_t iq = 0; iq + 1 < q2Edges.size(); ++iq) {
        const auto &xbEdges = fXbins.GetXBBins(iq);
        for (size_t ix = 0; ix + 1 < xbEdges.size(); ++ix) {
          for (size_t it = 0; it + 1 < tEdges.size(); ++it) {
            writeRow("conventional", static_cast<long long>(ix),
                     static_cast<long long>(iq), static_cast<long long>(it),
                     {xbEdges[ix], xbEdges[ix + 1],
                      q2Edges[iq], q2Edges[iq + 1],
                      tEdges[it], tEdges[it + 1]});
          }
        }
      }
    }

    csv.close();
    physicalVolumeCSVWritten_[energyKey] = true;
    std::cout << "[Veff] Wrote " << binIndex << " bins to " << filename << '\n';
  }

  // Enable or disable individual variable plotting
  void PlotIndividual(bool plotInd) { plotIndividual = plotInd; }

  // Set plot styles for various plot types
  void SetKinStyle(const DrawStyle& style) { styleKin_ = style; }
  void SetDVCSStyle(const DrawStyle& style) { styleDVCS_ = style; }
  void SetCrossSectionStyle(const DrawStyle& style) { styleCrossSection_ = style; }
  void SetBSAStyle(const DrawStyle& style) { styleBSA_ = style; }
  void UseFittedPhiYields(bool on = true) { useFittedYields_ = on; }

  // Enable or disable correctio
  void SetApplyCorrection(bool apply) { applyCorrection = apply; }

  void PlotMomentumCorrection(const std::string& subdir = "MomentumCorrection") {
    const std::string outdir = outputDir + "/" + subdir;
    const std::string ratioHistDir = outdir + "/trueOverExpHist2D";
    fs::create_directories(outdir);
    fs::create_directories(ratioHistDir);

    constexpr double kElectronMass = 0.00051099895;
    constexpr double kProtonMass = 0.9382720813;
    constexpr double kRadToDeg = 180.0 / M_PI;

    auto cleanLabel = [](std::string s) {
      for (auto& ch : s) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '-') ch = '_';
      }
      return s;
    };

    auto unitFromAngles = [](double theta, double phi) {
      return TVector3(std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi),
                      std::cos(theta));
    };

    struct PlotSpec {
      std::string tag;
      std::string ref;
      std::string calc;
      std::string title;
      std::string detectorColumn;
      int bins;
      double min;
      double max;
    };
    const std::vector<PlotSpec> specs = {
        {"Q2", "Q2", "momcorr_Q2", "Q^{2} [GeV^{2}]", "", 180, 0.0, 10.0},
        {"xB", "xB", "momcorr_xB", "x_{B}", "", 160, 0.0, 1.0},
        {"t", "momcorr_ref_abs_t", "momcorr_t", "|t| [GeV^{2}]", "", 180, 0.0, 2.0},
        {"phi", "phi", "momcorr_phi", "#phi_{BMK} [deg]", "", 180, 0.0, 360.0},
        {"photonE", "recpho_p", "momcorr_pho_E", "E_{#gamma} [GeV]", "pho_det_region", 180, 1.5, 7.0},
        {"protonP", "recpro_p", "momcorr_pro_p", "p_{p} [GeV]", "pro_det_region", 180, 0.0, 1.5},
    };
    const std::vector<std::pair<int, std::string>> detectorCases = {
        {0, "FT"},
        {1, "FD"},
        {2, "CD"},
    };
    auto detectorAllowed = [](const std::string& tag, int detector) {
      if (tag == "photonE") return detector == 0 || detector == 1;  // photon: FT, FD
      if (tag == "protonP") return detector == 1 || detector == 2;  // proton: FD, CD
      return true;
    };
    auto thetaRangeFor = [](const std::string& tag, int detector) {
      if (tag == "photonE" && detector == 0) return std::pair<double, double>{2.0, 5.0};
      if (tag == "photonE" && detector == 1) return std::pair<double, double>{0.0, 40.0};
      if (tag == "protonP" && detector == 1) return std::pair<double, double>{10.0, 50.0};
      if (tag == "protonP" && detector == 2) return std::pair<double, double>{30.0, 80.0};
      return std::pair<double, double>{0.0, 70.0};
    };
    auto ratioXBinsFor = [](const std::string& tag, int detector, const std::string& xTag) {
      if (tag == "photonE" && detector == 0 && xTag == "theta") return 20;  // photon FT theta
      if (tag == "photonE" && detector == 0 && xTag == "phi") return 36;    // photon FT phi
      if (tag == "photonE" && detector == 1 && xTag == "theta") return 10;  // photon FD theta
      if (tag == "photonE" && detector == 1 && xTag == "phi") return 18;    // photon FD phi
      if (tag == "protonP" && detector == 1 && xTag == "theta") return 12;  // proton FD theta
      if (tag == "protonP" && detector == 1 && xTag == "phi") return 12;    // proton FD phi
      if (tag == "protonP" && detector == 2 && xTag == "theta") return 20;  // proton CD theta
      if (tag == "protonP" && detector == 2 && xTag == "phi") return 18;    // proton CD phi
      return 50;
    };
    auto ratioHistXBinsFor = [](const std::string& tag, int detector, const std::string& xTag) {
      if (tag == "photonE" && detector == 0 && xTag == "theta") return 120;  // photon FT theta 2D
      if (tag == "photonE" && detector == 0 && xTag == "phi") return 180;    // photon FT phi 2D
      if (tag == "photonE" && detector == 1 && xTag == "theta") return 160;  // photon FD theta 2D
      if (tag == "photonE" && detector == 1 && xTag == "phi") return 180;    // photon FD phi 2D
      if (tag == "protonP" && detector == 1 && xTag == "theta") return 160;  // proton FD theta 2D
      if (tag == "protonP" && detector == 1 && xTag == "phi") return 180;    // proton FD phi 2D
      if (tag == "protonP" && detector == 2 && xTag == "theta") return 160;  // proton CD theta 2D
      if (tag == "protonP" && detector == 2 && xTag == "phi") return 180;    // proton CD phi 2D
      return 180;
    };
    auto minFitEntriesFor = [](const std::string& tag, int detector, const std::string& xTag) {
      if (tag == "photonE" && detector == 0 && xTag == "theta") return 25;
      if (tag == "photonE" && detector == 0 && xTag == "phi") return 25;
      if (tag == "photonE" && detector == 1 && xTag == "theta") return 50;
      if (tag == "photonE" && detector == 1 && xTag == "phi") return 50;
      if (tag == "protonP" && detector == 1 && xTag == "theta") return 50;
      if (tag == "protonP" && detector == 1 && xTag == "phi") return 50;
      if (tag == "protonP" && detector == 2 && xTag == "theta") return 40;
      if (tag == "protonP" && detector == 2 && xTag == "phi") return 40;
      return 50;
    };
    const std::vector<std::pair<double, double>> photonEnergyRanges = {
        {2.0, 2.5},
        {2.5, 3.0},
        {3.0, 3.5},
        {3.5, 4.0},
        {4.0, 4.5},
        {4.5, 5.0},
        {5.0, 5.5},
    };
    const std::vector<std::pair<double, double>> protonMomentumRanges = {
        {0.2, 0.4},
        {0.4, 0.6},
        {0.6, 0.8},
        {0.8, 1.0},
        {1.0, 1.2},
    };
    auto rangeTag = [](double min, double max) {
      std::string s = Form("_%.2f_%.2f", min, max);
      for (auto& ch : s) {
        if (ch == '.') ch = 'p';
        if (ch == '-') ch = 'm';
      }
      return s;
    };
    struct RatioPlotSpec {
      std::string tag;
      std::string ratio;
      std::string particle;
      std::string detectorColumn;
      std::string thetaColumn;
      std::string phiColumn;
      std::string rangeColumn;
      std::string rangeTitle;
      std::vector<std::pair<double, double>> ranges;
    };
    const std::vector<RatioPlotSpec> ratioSpecs = {
        {"photonE", "momcorr_pho_E_over_recpho_p", "#gamma", "pho_det_region", "recpho_theta_deg", "recpho_phi_deg", "recpho_p", "E_{#gamma}", photonEnergyRanges},
        {"protonP", "momcorr_pro_p_over_recpro_p", "p", "pro_det_region", "recpro_theta_deg", "recpro_phi_deg", "recpro_p", "p_{p}", protonMomentumRanges},
    };

    for (size_t imodel = 0; imodel < plotters.size(); ++imodel) {
      const std::string label = imodel < labels.size() ? labels[imodel] : Form("model_%zu", imodel);
      const std::string tag = cleanLabel(label);
      const double modelBeamEnergy = plotters[imodel]->GetBeamEnergy();

      auto q2FromElectron = [=](double p, double theta) {
        if (p <= 0.0 || theta < 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double eprime = std::sqrt(p * p + kElectronMass * kElectronMass);
        return 2.0 * modelBeamEnergy * (eprime - p * std::cos(theta)) - kElectronMass * kElectronMass;
      };

      auto nuFromElectron = [=](double p) {
        if (p <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double eprime = std::sqrt(p * p + kElectronMass * kElectronMass);
        return modelBeamEnergy - eprime;
      };

      auto xBFromElectron = [=](double p, double theta) {
        const double q2 = q2FromElectron(p, theta);
        const double nu = nuFromElectron(p);
        if (!std::isfinite(q2) || !std::isfinite(nu) || nu <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return q2 / (2.0 * kProtonMass * nu);
      };

      auto qVector = [=](double p, double theta, double phi) {
        const TVector3 kin(0.0, 0.0, modelBeamEnergy);
        const TVector3 kout = p * unitFromAngles(theta, phi);
        return kin - kout;
      };

      auto tFromElectronPhotonAngles = [=](double p, double theta, double phi, double phoTheta, double phoPhi) {
        const double q2 = q2FromElectron(p, theta);
        const double nu = nuFromElectron(p);
        const TVector3 qvec = qVector(p, theta, phi);
        const TVector3 ghat = unitFromAngles(phoTheta, phoPhi);
        const double qmag = qvec.Mag();
        if (!std::isfinite(q2) || !std::isfinite(nu) || nu <= 0.0 || qmag <= 0.0) return std::numeric_limits<double>::quiet_NaN();

        double cosThetaGG = qvec.Dot(ghat) / qmag;
        cosThetaGG = std::max(-1.0, std::min(1.0, cosThetaGG));
        const double sqrtNuQ = std::sqrt(nu * nu + q2);
        const double numerator = q2 * kProtonMass + 2.0 * nu * kProtonMass * (nu - sqrtNuQ * cosThetaGG);
        const double denominator = sqrtNuQ * cosThetaGG - nu - kProtonMass;
        if (std::abs(denominator) < 1e-12) return std::numeric_limits<double>::quiet_NaN();
        return std::abs(numerator / denominator);
      };

      auto phiBMKFromElectronPhotonAngles = [=](double p, double theta, double phi, double phoTheta, double phoPhi) {
        const TVector3 kin(0.0, 0.0, modelBeamEnergy);
        const TVector3 kout = p * unitFromAngles(theta, phi);
        const TVector3 qvec = kin - kout;
        const TVector3 gvec = unitFromAngles(phoTheta, phoPhi);
        if (qvec.Mag() <= 0.0) return std::numeric_limits<double>::quiet_NaN();

        const TVector3 nL = kin.Cross(kout).Unit();
        const TVector3 nH = gvec.Cross(qvec).Unit();
        if (nL.Mag() <= 0.0 || nH.Mag() <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double cosPhi = std::max(-1.0, std::min(1.0, nL.Dot(nH)));
        const double sinPhi = (nL.Cross(nH)).Dot(qvec.Unit());
        return (std::atan2(sinPhi, cosPhi) + M_PI) * kRadToDeg;
      };

      auto photonEnergyFromElectronPhotonAngles = [=](double p, double theta, double phi, double phoTheta, double phoPhi) {
        const double q2 = q2FromElectron(p, theta);
        const double absT = tFromElectronPhotonAngles(p, theta, phi, phoTheta, phoPhi);
        const double nu = nuFromElectron(p);
        const TVector3 qvec = qVector(p, theta, phi);
        const TVector3 ghat = unitFromAngles(phoTheta, phoPhi);
        const double qmag = qvec.Mag();
        if (!std::isfinite(q2) || !std::isfinite(absT) || !std::isfinite(nu) || qmag <= 0.0) return std::numeric_limits<double>::quiet_NaN();

        double cosThetaGG = qvec.Dot(ghat) / qmag;
        cosThetaGG = std::max(-1.0, std::min(1.0, cosThetaGG));
        const double denominator = 2.0 * (qmag * cosThetaGG - nu);
        if (std::abs(denominator) < 1e-12) return std::numeric_limits<double>::quiet_NaN();
        const double egamma = (q2 - absT) / denominator;
        return egamma > 0.0 ? egamma : std::numeric_limits<double>::quiet_NaN();
      };

      auto protonMomentumFromElectronPhotonAngles = [=](double p, double theta, double phi, double phoTheta, double phoPhi) {
        const double egamma = photonEnergyFromElectronPhotonAngles(p, theta, phi, phoTheta, phoPhi);
        if (!std::isfinite(egamma)) return std::numeric_limits<double>::quiet_NaN();
        const TVector3 qvec = qVector(p, theta, phi);
        const TVector3 gvec = egamma * unitFromAngles(phoTheta, phoPhi);
        return (qvec - gvec).Mag();
      };

      auto df = plotters[imodel]->GetRDF()
          .Define("momcorr_Q2", q2FromElectron, {"recel_p", "recel_theta"})
          .Define("momcorr_xB", xBFromElectron, {"recel_p", "recel_theta"})
          .Define("momcorr_t", tFromElectronPhotonAngles, {"recel_p", "recel_theta", "recel_phi", "recpho_theta", "recpho_phi"})
          .Define("momcorr_phi", phiBMKFromElectronPhotonAngles, {"recel_p", "recel_theta", "recel_phi", "recpho_theta", "recpho_phi"})
          .Define("momcorr_pho_E", photonEnergyFromElectronPhotonAngles, {"recel_p", "recel_theta", "recel_phi", "recpho_theta", "recpho_phi"})
          .Define("momcorr_pro_p", protonMomentumFromElectronPhotonAngles, {"recel_p", "recel_theta", "recel_phi", "recpho_theta", "recpho_phi"})
          .Define("momcorr_ref_abs_t", [](double t) { return std::abs(t); }, {"t"})
          .Define("momcorr_pho_E_over_recpho_p",
                  [](double calc, double rec) {
                    return (std::isfinite(calc) && std::isfinite(rec) && rec > 0.0)
                               ? calc / rec
                               : std::numeric_limits<double>::quiet_NaN();
                  },
                  {"momcorr_pho_E", "recpho_p"})
          .Define("momcorr_pro_p_over_recpro_p",
                  [](double calc, double rec) {
                    return (std::isfinite(calc) && std::isfinite(rec) && rec > 0.0)
                               ? calc / rec
                               : std::numeric_limits<double>::quiet_NaN();
                  },
                  {"momcorr_pro_p", "recpro_p"})
          .Define("recpho_theta_deg", [](double theta) { return theta * 180.0 / M_PI; }, {"recpho_theta"})
          .Define("recpro_theta_deg", [](double theta) { return theta * 180.0 / M_PI; }, {"recpro_theta"})
          .Define("recpho_phi_deg",
                  [](double phi) {
                    double phiDeg = phi * 180.0 / M_PI;
                    while (phiDeg < 0.0) phiDeg += 360.0;
                    while (phiDeg >= 360.0) phiDeg -= 360.0;
                    return phiDeg;
                  },
                  {"recpho_phi"})
          .Define("recpro_phi_deg",
                  [](double phi) {
                    double phiDeg = phi * 180.0 / M_PI;
                    while (phiDeg < 0.0) phiDeg += 360.0;
                    while (phiDeg >= 360.0) phiDeg -= 360.0;
                    return phiDeg;
                  },
                  {"recpro_phi"});

      std::map<std::string, std::unique_ptr<TGraphErrors>> a0VsMomentumGraphs;
      std::map<std::string, std::unique_ptr<TGraphErrors>> a1VsMomentumGraphs;
      std::map<std::string, std::string> paramGraphTitles;
      std::map<std::string, std::string> paramGraphXTitles;
      std::array<std::array<double, 2>, 3> cdPhiA0MomentumCoefficients{};
      std::array<std::array<double, 2>, 3> cdPhiA1MomentumCoefficients{};
      bool cdPhiCorrectionModelValid = false;

      auto shouldCollectParam = [](const std::string& ratioTag, int detector,
                                   const std::string& xTag, const std::string& rangeName) {
        if (rangeName == "_all") return false;
        if (xTag == "theta") {
          if (ratioTag == "photonE") return detector == 0 || detector == 1;  // photon FT/FD
          if (ratioTag == "protonP") return detector == 1 || detector == 2;   // proton FD/CD
        }
        return xTag == "phi" && ratioTag == "protonP" && detector == 2;
      };

      auto fillParamGraphs = [&](const std::string& ratioTag, const std::string& detectorName,
                                 const std::string& keySuffix,
                                 const std::string& particleLabel, const std::string& rangeAxisTitle,
                                 double rangeCenter, double rangeError, TF1* fit) {
        if (!fit || !std::isfinite(rangeCenter) || !std::isfinite(rangeError) || rangeError <= 0.0) return;
        const double a0Err = fit->GetParError(0);
        if (!std::isfinite(a0Err) || a0Err > 0.1) return;

        const std::string key = ratioTag + "_" + detectorName + keySuffix;
        if (!a0VsMomentumGraphs.count(key)) {
          a0VsMomentumGraphs[key] = std::make_unique<TGraphErrors>();
          a1VsMomentumGraphs[key] = std::make_unique<TGraphErrors>();
          a0VsMomentumGraphs[key]->SetName(Form("g_%s_%s_a0_vs_momentum", tag.c_str(), key.c_str()));
          a1VsMomentumGraphs[key]->SetName(Form("g_%s_%s_a1_vs_momentum", tag.c_str(), key.c_str()));
          const std::string segmentTitle = keySuffix.empty() ? "" : " " + keySuffix.substr(1);
          paramGraphTitles[key] = Form("%s %s %s%s", label.c_str(), particleLabel.c_str(),
                                       detectorName.c_str(), segmentTitle.c_str());
          paramGraphXTitles[key] = rangeAxisTitle + " [GeV]";
        }

        auto addPoint = [&](TGraphErrors* graph, double y, double ey) {
          if (!graph || !std::isfinite(y) || !std::isfinite(ey)) return;
          const int ip = graph->GetN();
          graph->SetPoint(ip, rangeCenter, y);
          graph->SetPointError(ip, rangeError, ey);
        };

        addPoint(a0VsMomentumGraphs[key].get(), fit->GetParameter(0), fit->GetParError(0));
        addPoint(a1VsMomentumGraphs[key].get(), fit->GetParameter(1), fit->GetParError(1));
      };

      for (const auto& spec : specs) {
        auto drawOne = [&](ROOT::RDF::RNode dfIn, const std::string& detTag, const std::string& detTitle) {
          auto dfGood = dfIn.Filter([](double ref, double calc) {
            return std::isfinite(ref) && std::isfinite(calc);
          }, {spec.ref, spec.calc});

          auto h = dfGood.Histo2D({Form("h_momcorr_%s_%s%s", tag.c_str(), spec.tag.c_str(), detTag.c_str()),
                                   Form("%s%s;%s;%s from e momentum,angle and #gamma angle;Events",
                                        label.c_str(), detTitle.c_str(), spec.title.c_str(), spec.title.c_str()),
                                   spec.bins, spec.min, spec.max,
                                   spec.bins, spec.min, spec.max},
                                  spec.ref, spec.calc);

          TCanvas c(Form("c_momcorr_%s_%s%s", tag.c_str(), spec.tag.c_str(), detTag.c_str()), "", 950, 850);
          c.SetRightMargin(0.16);
          c.SetLeftMargin(0.13);
          c.SetBottomMargin(0.13);
          h->SetStats(0);
          h->Draw("COLZ");
          TLine diag(spec.min, spec.min, spec.max, spec.max);
          diag.SetLineColor(kRed + 1);
          diag.SetLineStyle(2);
          diag.SetLineWidth(2);
          diag.Draw("SAME");
          c.SaveAs(Form("%s/%s_%s%s_momentumCorrection2D.png", outdir.c_str(), tag.c_str(), spec.tag.c_str(), detTag.c_str()));
        };

        if (spec.detectorColumn.empty()) {
          drawOne(df, "", "");
          continue;
        }

        if (!df.HasColumn(spec.detectorColumn)) {
          std::cerr << "[PlotMomentumCorrection] Column " << spec.detectorColumn
                    << " not found for model " << label << "; drawing " << spec.tag
                    << " without detector split.\n";
          drawOne(df, "", "");
          continue;
        }

        for (const auto& detCase : detectorCases) {
          if (!detectorAllowed(spec.tag, detCase.first)) continue;
          auto dfDet = df.Filter([region = detCase.first](int detRegion) {
            return detRegion == region;
          }, {spec.detectorColumn});
          drawOne(dfDet, "_" + detCase.second, " " + detCase.second);
        }
      }

      for (const auto& ratioSpec : ratioSpecs) {
        if (!df.HasColumn(ratioSpec.detectorColumn)) {
          std::cerr << "[PlotMomentumCorrection] Column " << ratioSpec.detectorColumn
                    << " not found for model " << label << "; skipping " << ratioSpec.tag
                    << " true/exp detector plots.\n";
          continue;
        }

        for (const auto& detCase : detectorCases) {
          if (!detectorAllowed(ratioSpec.tag, detCase.first)) continue;
          auto dfDet = df.Filter([region = detCase.first](int detRegion) {
            return detRegion == region;
          }, {ratioSpec.detectorColumn});

          auto drawRatio = [&](ROOT::RDF::RNode dfRange, const std::string& rangeName, const std::string& rangeTitle,
                               double rangeCenter, double rangeError,
                               const std::string& xColumn, const std::string& xTag,
                               const std::string& xTitle, int histXBins, int fitXBins,
                               double xMin, double xMax, const std::string& ratioColumn,
                               const std::string& stageTag, const std::string& ratioTitle,
                               bool collectFitParams) {
            const int minFitEntries = minFitEntriesFor(ratioSpec.tag, detCase.first, xTag);
            const std::string outputRatioTag = stageTag.empty()
                                                   ? "trueOverExp"
                                                   : "trueOverPhiCorrectedExp";
            const std::string outputPeakRatioTag = stageTag.empty()
                                                       ? "TrueOverExp"
                                                       : "TrueOverPhiCorrectedExp";
            auto dfGood = dfRange.Filter([](double x, double ratio) {
              return std::isfinite(x) && std::isfinite(ratio);
            }, {xColumn, ratioColumn});

            auto h = dfGood.Histo2D({Form("h_momcorr_%s_%s_%s%s%s_vs_%s", tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(), stageTag.c_str(), rangeName.c_str(), xTag.c_str()),
                                     Form("%s %s %s %s;%s;%s %s;Events",
                                          label.c_str(), ratioSpec.particle.c_str(), detCase.second.c_str(), rangeTitle.c_str(), xTitle.c_str(), ratioSpec.particle.c_str(), ratioTitle.c_str()),
                                     histXBins, xMin, xMax,
                                     160, 0.5, 1.5},
                                    xColumn, ratioColumn);

            TCanvas c(Form("c_momcorr_%s_%s_%s%s%s_vs_%s", tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(), stageTag.c_str(), rangeName.c_str(), xTag.c_str()), "", 950, 850);
            c.SetRightMargin(0.16);
            c.SetLeftMargin(0.13);
            c.SetBottomMargin(0.13);
            h->SetStats(0);
            h->Draw("COLZ");
            TLine unity(xMin, 1.0, xMax, 1.0);
            unity.SetLineColor(kRed + 1);
            unity.SetLineStyle(2);
            unity.SetLineWidth(2);
            unity.Draw("SAME");
            c.SaveAs(Form("%s/%s_%s_%s%s%s_%s_vs_%s_momentumCorrection2D.png",
                          ratioHistDir.c_str(), tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                          stageTag.c_str(), rangeName.c_str(), outputRatioTag.c_str(), xTag.c_str()));

            auto hFit = dfGood.Histo2D({Form("hfit_momcorr_%s_%s_%s%s%s_vs_%s", tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(), stageTag.c_str(), rangeName.c_str(), xTag.c_str()),
                                        Form("%s %s %s %s;%s;%s %s;Events",
                                             label.c_str(), ratioSpec.particle.c_str(), detCase.second.c_str(), rangeTitle.c_str(), xTitle.c_str(), ratioSpec.particle.c_str(), ratioTitle.c_str()),
                                        fitXBins, xMin, xMax,
                                        160, 0.5, 1.5},
                                       xColumn, ratioColumn);
            TH2D* h2 = hFit.GetPtr();
            TGraphErrors graph;
            graph.SetName(Form("g_peak_%s_%s_%s%s%s_vs_%s", tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(), stageTag.c_str(), rangeName.c_str(), xTag.c_str()));
            graph.SetTitle(Form("%s %s %s %s;%s;%s peak %s",
                                label.c_str(), ratioSpec.particle.c_str(), detCase.second.c_str(), rangeTitle.c_str(), xTitle.c_str(), ratioSpec.particle.c_str(), ratioTitle.c_str()));

            int point = 0;
            for (int ibin = 1; ibin <= h2->GetNbinsX(); ++ibin) {
              std::unique_ptr<TH1D> proj(h2->ProjectionY(Form("py_%s_%s_%s%s%s_%s_%d",
                                                              tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                                                              stageTag.c_str(), rangeName.c_str(), xTag.c_str(), ibin),
                                                        ibin, ibin));
              if (!proj || proj->GetEntries() < minFitEntries || proj->Integral() <= 0.0) continue;

              const double mean = proj->GetMean();
              const double rms = proj->GetRMS();
              if (!std::isfinite(mean) || !std::isfinite(rms) || rms <= 0.0) continue;

              TF1 gaus(Form("fit_%s_%s_%s%s%s_%s_%d",
                            tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                            stageTag.c_str(), rangeName.c_str(), xTag.c_str(), ibin),
                       "gaus", 0.5, 1.5);
              gaus.SetParameters(proj->GetMaximum(), mean, rms);
              const int fitStatus = proj->Fit(&gaus, "QNR");
              const double peak = gaus.GetParameter(1);
              const double peakErr = gaus.GetParError(1);
              if (fitStatus != 0 || !std::isfinite(peak) || !std::isfinite(peakErr) || peakErr <= 0.0) continue;
              if (peak < 0.5 || peak > 1.5) continue;

              const double xCenter = h2->GetXaxis()->GetBinCenter(ibin);
              const double xErr = 0.5 * h2->GetXaxis()->GetBinWidth(ibin);
              graph.SetPoint(point, xCenter, peak);
              graph.SetPointError(point, xErr, peakErr);
              ++point;
            }

            TCanvas cFit(Form("c_peak_%s_%s_%s%s%s_vs_%s", tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(), stageTag.c_str(), rangeName.c_str(), xTag.c_str()), "", 950, 850);
            cFit.SetRightMargin(0.05);
            cFit.SetLeftMargin(0.13);
            cFit.SetBottomMargin(0.13);
            TH1D frame(Form("frame_peak_%s_%s_%s%s%s_vs_%s", tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(), stageTag.c_str(), rangeName.c_str(), xTag.c_str()),
                       graph.GetTitle(), fitXBins, xMin, xMax);
            frame.SetMinimum(0.5);
            frame.SetMaximum(1.5);
            frame.SetStats(0);
            frame.Draw();
            const bool splitProtonCDPhi =
                ratioSpec.tag == "protonP" && detCase.first == 2 && xTag == "phi";
            if (splitProtonCDPhi) {
              struct PhiSegment {
                const char* name;
                double phiMin;
                double phiMax;
                double phiOrigin;
                bool wrapsZero;
              };
              const std::array<PhiSegment, 3> phiSegments = {{
                  {"phiSeg1", 20.0, 140.0, 20.0, false},
                  {"phiSeg2", 150.0, 250.0, 150.0, false},
                  {"phiSeg3", 260.0, 370.0, 260.0, true},
              }};
              const std::array<int, 3> segmentColors = {kRed + 1, kBlue + 1, kGreen + 2};

              for (size_t iseg = 0; iseg < phiSegments.size(); ++iseg) {
                const auto& segment = phiSegments[iseg];
                TGraphErrors segmentGraph;
                segmentGraph.SetName(Form("gfit_%s_%s_%s%s%s_%s",
                                          tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                                          stageTag.c_str(), rangeName.c_str(), segment.name));
                for (int ip = 0; ip < graph.GetN(); ++ip) {
                  double phiDeg = 0.0;
                  double peak = 0.0;
                  graph.GetPoint(ip, phiDeg, peak);
                  double phiUnwrapped = phiDeg;
                  if (segment.wrapsZero && phiDeg < 20.0) phiUnwrapped += 360.0;
                  const bool inSegment = segment.wrapsZero
                                             ? phiUnwrapped >= segment.phiMin && phiUnwrapped <= segment.phiMax
                                             : phiUnwrapped >= segment.phiMin && phiUnwrapped < segment.phiMax;
                  if (!inSegment) continue;
                  const int jp = segmentGraph.GetN();
                  segmentGraph.SetPoint(jp, phiUnwrapped - segment.phiOrigin, peak);
                  segmentGraph.SetPointError(jp, graph.GetErrorX(ip), graph.GetErrorY(ip));
                }
                if (segmentGraph.GetN() <= 3) continue;

                auto* segmentFit = new TF1(
                    Form("linfit_%s_%s_%s%s%s_%s",
                         tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                         stageTag.c_str(), rangeName.c_str(), segment.name),
                    "[0]+[1]*x", 0.0, segment.phiMax - segment.phiOrigin);
                const int fitStatus = segmentGraph.Fit(segmentFit, "QNEX0");
                if (fitStatus != 0) {
                  delete segmentFit;
                  continue;
                }

                if (collectFitParams && shouldCollectParam(ratioSpec.tag, detCase.first, xTag, rangeName)) {
                  fillParamGraphs(ratioSpec.tag, detCase.second,
                                  std::string("_") + segment.name,
                                  ratioSpec.particle, ratioSpec.rangeTitle,
                                  rangeCenter, rangeError, segmentFit);
                }

                auto drawSegmentLine = [&](const char* suffix, const char* formula,
                                           double drawMin, double drawMax, double offset) {
                  auto* line = new TF1(
                      Form("displayfit_%s_%s_%s%s%s_%s_%s",
                           tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                           stageTag.c_str(), rangeName.c_str(), segment.name, suffix),
                      formula, drawMin, drawMax);
                  line->SetParameter(0, segmentFit->GetParameter(0));
                  line->SetParameter(1, segmentFit->GetParameter(1));
                  line->SetParameter(2, offset);
                  line->SetLineColor(segmentColors[iseg]);
                  line->SetLineWidth(3);
                  line->SetBit(TObject::kCanDelete);
                  line->Draw("SAME");
                };

                if (!segment.wrapsZero) {
                  drawSegmentLine("main", "[0]+[1]*(x-[2])",
                                  segment.phiMin, segment.phiMax, segment.phiOrigin);
                } else {
                  drawSegmentLine("high", "[0]+[1]*(x-[2])", 260.0, 360.0, 260.0);
                  drawSegmentLine("low", "[0]+[1]*(x+[2])", 0.0, 10.0, 100.0);
                }
                delete segmentFit;
              }
            } else if (graph.GetN() > 3) {
              auto* linearFit = new TF1(Form("linfit_%s_%s_%s%s%s_vs_%s",
                                             tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                                             stageTag.c_str(), rangeName.c_str(), xTag.c_str()),
                                        "[0]+[1]*x", xMin, xMax);
              const int linearFitStatus = graph.Fit(linearFit, "QEX0");
              if (linearFitStatus == 0) {
                auto* fitBand = new TGraphErrors(200);
                fitBand->SetName(Form("band_%s_%s_%s%s%s_vs_%s",
                                      tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                                      stageTag.c_str(), rangeName.c_str(), xTag.c_str()));
                for (int ip = 0; ip < 200; ++ip) {
                  const double x = xMin + (xMax - xMin) * static_cast<double>(ip) / 199.0;
                  fitBand->SetPoint(ip, x, linearFit->Eval(x));
                  fitBand->SetPointError(ip, 0.0, 0.0);
                }
                TVirtualFitter* fitter = TVirtualFitter::GetFitter();
                if (fitter) {
                  fitter->GetConfidenceIntervals(fitBand, 0.68);
                  fitBand->SetFillColorAlpha(kRed, 0.22);
                  fitBand->SetLineColor(kRed);
                  fitBand->SetBit(TObject::kCanDelete);
                  fitBand->Draw("E3 SAME");
                }
                linearFit->SetLineColor(kRed + 1);
                linearFit->SetLineWidth(3);
                linearFit->SetBit(TObject::kCanDelete);
                linearFit->Draw("SAME");
                if (collectFitParams && shouldCollectParam(ratioSpec.tag, detCase.first, xTag, rangeName)) {
                  const std::string particleLabel = stageTag.empty()
                                                        ? ratioSpec.particle
                                                        : ratioSpec.particle + " after #phi correction";
                  fillParamGraphs(ratioSpec.tag, detCase.second, "", particleLabel, ratioSpec.rangeTitle,
                                  rangeCenter, rangeError, linearFit);
                }
              }
            }
            if (graph.GetN() > 0) {
              graph.SetMarkerStyle(20);
              graph.SetMarkerSize(0.8);
              graph.SetMarkerColor(kBlack);
              graph.SetLineColor(kBlack);
              graph.Draw("P SAME");
            }
            TLine unityPeak(xMin, 1.0, xMax, 1.0);
            unityPeak.SetLineColor(kGray + 2);
            unityPeak.SetLineStyle(2);
            unityPeak.SetLineWidth(2);
            unityPeak.Draw("SAME");
            cFit.SaveAs(Form("%s/%s_%s_%s%s%s_peak%s_vs_%s_momentumCorrection.png",
                             outdir.c_str(), tag.c_str(), ratioSpec.tag.c_str(), detCase.second.c_str(),
                             stageTag.c_str(), rangeName.c_str(), outputPeakRatioTag.c_str(), xTag.c_str()));

          };

          const auto thetaRange = thetaRangeFor(ratioSpec.tag, detCase.first);
          const bool sequentialProtonCD = ratioSpec.tag == "protonP" && detCase.first == 2;
          if (!sequentialProtonCD) {
            drawRatio(dfDet, "_all", "all", std::numeric_limits<double>::quiet_NaN(), 0.0,
                      ratioSpec.thetaColumn, "theta", "#theta [deg]",
                      ratioHistXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                      ratioXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                      thetaRange.first, thetaRange.second,
                      ratioSpec.ratio, "", "true/exp", true);
          }
          drawRatio(dfDet, "_all", "all", std::numeric_limits<double>::quiet_NaN(), 0.0,
                    ratioSpec.phiColumn, "phi", "#phi [deg]",
                    ratioHistXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                    ratioXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                    0.0, 360.0,
                    ratioSpec.ratio, "", "true/exp", true);
          for (size_t irange = 0; irange < ratioSpec.ranges.size(); ++irange) {
            const auto& range = ratioSpec.ranges[irange];
            const bool lastRange = irange + 1 == ratioSpec.ranges.size();
            auto dfRange = dfDet.Filter([min = range.first, max = range.second, lastRange](double value) {
              return std::isfinite(value) && value >= min && (lastRange ? value <= max : value < max);
            }, {ratioSpec.rangeColumn});
            const std::string thisRangeTag = rangeTag(range.first, range.second);
            const std::string thisRangeTitle = Form("%s [%.2f, %.2f%s GeV",
                                                    ratioSpec.rangeTitle.c_str(), range.first, range.second, lastRange ? "]" : ")");
            const double rangeCenter = 0.5 * (range.first + range.second);
            const double rangeError = 0.5 * (range.second - range.first);
            if (!sequentialProtonCD) {
              drawRatio(dfRange, thisRangeTag, thisRangeTitle, rangeCenter, rangeError,
                        ratioSpec.thetaColumn, "theta", "#theta [deg]",
                        ratioHistXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                        ratioXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                        thetaRange.first, thetaRange.second,
                        ratioSpec.ratio, "", "true/exp", true);
            }
            drawRatio(dfRange, thisRangeTag, thisRangeTitle, rangeCenter, rangeError,
                      ratioSpec.phiColumn, "phi", "#phi [deg]",
                      ratioHistXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                      ratioXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                      0.0, 360.0,
                      ratioSpec.ratio, "", "true/exp", true);
          }

          if (sequentialProtonCD) {
            auto fitMomentumDependence = [&](TGraphErrors* source, const std::string& fitName,
                                             double errorFloor, std::array<double, 2>& coefficients) {
              if (!source || source->GetN() < 2) return false;
              TGraphErrors fitGraph;
              double pMin = std::numeric_limits<double>::infinity();
              double pMax = -std::numeric_limits<double>::infinity();
              for (int ip = 0; ip < source->GetN(); ++ip) {
                double p = 0.0;
                double value = 0.0;
                source->GetPoint(ip, p, value);
                const int jp = fitGraph.GetN();
                fitGraph.SetPoint(jp, p, value);
                fitGraph.SetPointError(jp, source->GetErrorX(ip),
                                       std::max(source->GetErrorY(ip), errorFloor));
                pMin = std::min(pMin, p);
                pMax = std::max(pMax, p);
              }
              if (!std::isfinite(pMin) || !std::isfinite(pMax) || pMax <= pMin) return false;
              TF1 fit(fitName.c_str(), "[0]+[1]*x", pMin, pMax);
              if (fitGraph.Fit(&fit, "QNEX0") != 0) return false;
              coefficients = {fit.GetParameter(0), fit.GetParameter(1)};
              return std::isfinite(coefficients[0]) && std::isfinite(coefficients[1]);
            };

            bool allPhiFitsValid = true;
            for (size_t iseg = 0; iseg < 3; ++iseg) {
              const std::string segmentKey = "protonP_CD_phiSeg" + std::to_string(iseg + 1);
              const auto a0It = a0VsMomentumGraphs.find(segmentKey);
              const auto a1It = a1VsMomentumGraphs.find(segmentKey);
              const bool a0Valid = a0It != a0VsMomentumGraphs.end() &&
                                   fitMomentumDependence(a0It->second.get(),
                                       Form("fit_%s_%s_a0_forSequentialCorrection", tag.c_str(), segmentKey.c_str()),
                                       0.006, cdPhiA0MomentumCoefficients[iseg]);
              const bool a1Valid = a1It != a1VsMomentumGraphs.end() &&
                                   fitMomentumDependence(a1It->second.get(),
                                       Form("fit_%s_%s_a1_forSequentialCorrection", tag.c_str(), segmentKey.c_str()),
                                       0.00035, cdPhiA1MomentumCoefficients[iseg]);
              allPhiFitsValid = allPhiFitsValid && a0Valid && a1Valid;
            }
            cdPhiCorrectionModelValid = allPhiFitsValid;

            if (!cdPhiCorrectionModelValid) {
              std::cerr << "[PlotMomentumCorrection] Could not build all three CD proton phi models for "
                        << label << "; skipping the phi-corrected residual fits.\n";
            } else {
              auto phiCorrectionScale = [a0Coefficients = cdPhiA0MomentumCoefficients,
                                         a1Coefficients = cdPhiA1MomentumCoefficients]
                                        (double p, double phiDeg) {
                if (!std::isfinite(p) || !std::isfinite(phiDeg) || p <= 0.0)
                  return std::numeric_limits<double>::quiet_NaN();
                while (phiDeg < 0.0) phiDeg += 360.0;
                while (phiDeg >= 360.0) phiDeg -= 360.0;

                size_t segment = 0;
                double phiLocal = 0.0;
                if (phiDeg >= 20.0 && phiDeg < 140.0) {
                  segment = 0;
                  phiLocal = phiDeg - 20.0;
                } else if (phiDeg >= 150.0 && phiDeg < 250.0) {
                  segment = 1;
                  phiLocal = phiDeg - 150.0;
                } else if (phiDeg >= 260.0 || phiDeg < 10.0) {
                  segment = 2;
                  const double phiUnwrapped = phiDeg < 10.0 ? phiDeg + 360.0 : phiDeg;
                  phiLocal = phiUnwrapped - 260.0;
                } else {
                  return std::numeric_limits<double>::quiet_NaN();
                }

                const double a0 = a0Coefficients[segment][0] + a0Coefficients[segment][1] * p;
                const double a1 = a1Coefficients[segment][0] + a1Coefficients[segment][1] * p;
                const double scale = a0 + a1 * phiLocal;
                return std::isfinite(scale) && scale > 0.0
                           ? scale
                           : std::numeric_limits<double>::quiet_NaN();
              };

              auto dfAfterPhi = dfDet
                  .Define("momcorr_pro_phi_scale", phiCorrectionScale,
                          {"recpro_p", "recpro_phi_deg"})
                  .Define("momcorr_pro_p_after_phi",
                          [](double p, double scale) {
                            return std::isfinite(p) && std::isfinite(scale)
                                       ? p * scale
                                       : std::numeric_limits<double>::quiet_NaN();
                          },
                          {"recpro_p", "momcorr_pro_phi_scale"})
                  .Define("momcorr_pro_ratio_after_phi",
                          [](double expected, double corrected) {
                            return std::isfinite(expected) && std::isfinite(corrected) && corrected > 0.0
                                       ? expected / corrected
                                       : std::numeric_limits<double>::quiet_NaN();
                          },
                          {"momcorr_pro_p", "momcorr_pro_p_after_phi"});

              drawRatio(dfAfterPhi, "_all", "all after #phi correction",
                        std::numeric_limits<double>::quiet_NaN(), 0.0,
                        ratioSpec.phiColumn, "phi", "#phi [deg]",
                        ratioHistXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                        ratioXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                        0.0, 360.0, "momcorr_pro_ratio_after_phi",
                        "_afterPhiCorrection", "true/(#phi-corrected exp)", false);
              drawRatio(dfAfterPhi, "_all", "all after #phi correction",
                        std::numeric_limits<double>::quiet_NaN(), 0.0,
                        ratioSpec.thetaColumn, "theta", "#theta [deg]",
                        ratioHistXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                        ratioXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                        thetaRange.first, thetaRange.second, "momcorr_pro_ratio_after_phi",
                        "_afterPhiCorrection", "true/(#phi-corrected exp)", false);

              for (size_t irange = 0; irange < ratioSpec.ranges.size(); ++irange) {
                const auto& range = ratioSpec.ranges[irange];
                const bool lastRange = irange + 1 == ratioSpec.ranges.size();
                auto dfRangeAfterPhi = dfAfterPhi.Filter(
                    [min = range.first, max = range.second, lastRange](double value) {
                      return std::isfinite(value) && value >= min &&
                             (lastRange ? value <= max : value < max);
                    },
                    {"momcorr_pro_p_after_phi"});
                const std::string thisRangeTag = rangeTag(range.first, range.second);
                const std::string thisRangeTitle = Form("p_{p}^{#phi corr} [%.2f, %.2f%s GeV",
                                                        range.first, range.second, lastRange ? "]" : ")");
                const double rangeCenter = 0.5 * (range.first + range.second);
                const double rangeError = 0.5 * (range.second - range.first);
                drawRatio(dfRangeAfterPhi, thisRangeTag, thisRangeTitle, rangeCenter, rangeError,
                          ratioSpec.phiColumn, "phi", "#phi [deg]",
                          ratioHistXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                          ratioXBinsFor(ratioSpec.tag, detCase.first, "phi"),
                          0.0, 360.0, "momcorr_pro_ratio_after_phi",
                          "_afterPhiCorrection", "true/(#phi-corrected exp)", false);
                drawRatio(dfRangeAfterPhi, thisRangeTag, thisRangeTitle, rangeCenter, rangeError,
                          ratioSpec.thetaColumn, "theta", "#theta [deg]",
                          ratioHistXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                          ratioXBinsFor(ratioSpec.tag, detCase.first, "theta"),
                          thetaRange.first, thetaRange.second, "momcorr_pro_ratio_after_phi",
                          "_afterPhiCorrection", "true/(#phi-corrected exp)", true);
              }
            }
          }
        }
      }

      std::map<std::string, std::array<double, 2>> momentumParameterFitCoefficients;
      auto drawParamGraph = [&](TGraphErrors* graph, const std::string& key,
                                const std::string& paramName, const std::string& yTitle) {
        if (!graph || graph->GetN() <= 0) return;

        double xMin = std::numeric_limits<double>::infinity();
        double xMax = -std::numeric_limits<double>::infinity();
        double yMin = std::numeric_limits<double>::infinity();
        double yMax = -std::numeric_limits<double>::infinity();
        for (int ip = 0; ip < graph->GetN(); ++ip) {
          double x = 0.0;
          double y = 0.0;
          graph->GetPoint(ip, x, y);
          const double ex = graph->GetErrorX(ip);
          const double ey = graph->GetErrorY(ip);
          xMin = std::min(xMin, x - ex);
          xMax = std::max(xMax, x + ex);
          yMin = std::min(yMin, y - ey);
          yMax = std::max(yMax, y + ey);
        }
        if (!std::isfinite(xMin) || !std::isfinite(xMax) || xMax <= xMin) return;
        if (!std::isfinite(yMin) || !std::isfinite(yMax)) return;
        if (yMax <= yMin) {
          const double center = 0.5 * (yMin + yMax);
          const double halfRange = std::max(1e-6, 0.05 * std::abs(center));
          yMin = center - halfRange;
          yMax = center + halfRange;
        }
        const double yPad = std::max(1e-6, 0.18 * (yMax - yMin));

        TCanvas cParam(Form("c_%s_%s_%s_vs_momentum", tag.c_str(), key.c_str(), paramName.c_str()), "", 950, 850);
        cParam.SetRightMargin(0.05);
        cParam.SetLeftMargin(0.14);
        cParam.SetBottomMargin(0.13);
        TH1D frame(Form("frame_%s_%s_%s_vs_momentum", tag.c_str(), key.c_str(), paramName.c_str()),
                   Form("%s;%s;%s", paramGraphTitles[key].c_str(), paramGraphXTitles[key].c_str(), yTitle.c_str()),
                   100, xMin, xMax);
        frame.SetMinimum(yMin - yPad);
        frame.SetMaximum(yMax + yPad);
        frame.SetStats(0);
        frame.GetXaxis()->SetTitleSize(0.052);
        frame.GetYaxis()->SetTitleSize(0.052);
        frame.GetXaxis()->SetLabelSize(0.044);
        frame.GetYaxis()->SetLabelSize(0.044);
        frame.GetYaxis()->SetTitleOffset(paramName == "a1" ? 1.45 : 1.22);
        if (paramName == "a1") frame.GetYaxis()->SetNoExponent(false);
        frame.Draw();

        TF1* fit = nullptr;
        if (graph->GetN() >= 2) {
          const double yErrorFloor = paramName == "a1" ? 0.00035 : 0.006;
          auto* fitGraph = new TGraphErrors();
          fitGraph->SetName(Form("fitGraph_%s_%s_%s_vs_momentum", tag.c_str(), key.c_str(), paramName.c_str()));
          for (int ip = 0; ip < graph->GetN(); ++ip) {
            double x = 0.0;
            double y = 0.0;
            graph->GetPoint(ip, x, y);
            fitGraph->SetPoint(ip, x, y);
            fitGraph->SetPointError(ip, graph->GetErrorX(ip), std::max(graph->GetErrorY(ip), yErrorFloor));
          }
          fit = new TF1(Form("fit_%s_%s_%s_vs_momentum", tag.c_str(), key.c_str(), paramName.c_str()),
                        "[0]+[1]*x", xMin, xMax);
          const int fitStatus = fitGraph->Fit(fit, "QEX0");
          if (fitStatus == 0) {
            momentumParameterFitCoefficients[key + "_" + paramName] = {
                fit->GetParameter(0), fit->GetParameter(1)};
            auto* fitBand = new TGraphErrors(200);
            fitBand->SetName(Form("band_%s_%s_%s_vs_momentum", tag.c_str(), key.c_str(), paramName.c_str()));
            for (int ip = 0; ip < 200; ++ip) {
              const double x = xMin + (xMax - xMin) * static_cast<double>(ip) / 199.0;
              fitBand->SetPoint(ip, x, fit->Eval(x));
              fitBand->SetPointError(ip, 0.0, 0.0);
            }
            TVirtualFitter* fitter = TVirtualFitter::GetFitter();
            if (fitter) {
              fitter->GetConfidenceIntervals(fitBand, 0.68);
              fitBand->SetFillColorAlpha(kRed, 0.22);
              fitBand->SetLineColor(kRed);
              fitBand->SetBit(TObject::kCanDelete);
              fitBand->Draw("E3 SAME");
            }
            fit->SetLineColor(kRed + 1);
            fit->SetLineWidth(3);
            fit->SetBit(TObject::kCanDelete);
            fit->Draw("SAME");
          } else {
            delete fit;
            fit = nullptr;
          }
          delete fitGraph;
        }

        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.9);
        graph->SetMarkerColor(kBlack);
        graph->SetLineColor(kBlack);
        graph->Draw("PE SAME");

        cParam.SaveAs(Form("%s/%s_%s_%s_vs_momentum.png", outdir.c_str(), tag.c_str(), key.c_str(), paramName.c_str()));
        std::ofstream fitOut(Form("%s/%s_%s_%s_vs_momentum_fit.txt", outdir.c_str(), tag.c_str(), key.c_str(), paramName.c_str()));
        fitOut << "# Fit expression: [0]+[1]*p\n";
        fitOut << "# Graph: " << key << " " << paramName << "\n";
        fitOut << "# x axis: " << paramGraphXTitles[key] << "\n";
        fitOut << std::setprecision(12);
        if (fit) {
          fitOut << "p0 = " << fit->GetParameter(0) << " +/- " << fit->GetParError(0) << "\n";
          fitOut << "p1 = " << fit->GetParameter(1) << " +/- " << fit->GetParError(1) << "\n";
          fitOut << "chi2 = " << fit->GetChisquare() << "\n";
          fitOut << "ndf = " << fit->GetNDF() << "\n";
        } else {
          fitOut << "# Fit was not performed or failed.\n";
        }
      };

      for (const auto& item : a0VsMomentumGraphs) {
        const std::string& key = item.first;
        drawParamGraph(item.second.get(), key, "a0", "a_{0}");
        if (a1VsMomentumGraphs.count(key)) {
          drawParamGraph(a1VsMomentumGraphs[key].get(), key, "a1", "a_{1} [1/deg]");
        }
      }

      const auto thetaA0It = momentumParameterFitCoefficients.find("protonP_CD_a0");
      const auto thetaA1It = momentumParameterFitCoefficients.find("protonP_CD_a1");
      if (cdPhiCorrectionModelValid &&
          thetaA0It != momentumParameterFitCoefficients.end() &&
          thetaA1It != momentumParameterFitCoefficients.end()) {
        std::ofstream totalFitOut(
            Form("%s/%s_protonP_CD_totalCorrection_fit.txt", outdir.c_str(), tag.c_str()));
        totalFitOut << std::setprecision(12);
        totalFitOut << "# Sequential CD proton momentum correction\n";
        totalFitOut << "# Step 1: Cphi = a0_phi(p) + a1_phi(p) * phiLocal\n";
        totalFitOut << "# Step 2: pPhi = p * Cphi\n";
        totalFitOut << "# Step 3: Ctheta = a0_theta(pPhi) + a1_theta(pPhi) * thetaDeg\n";
        totalFitOut << "# Final: pCorrected = pPhi * Ctheta = p * Cphi * Ctheta\n\n";
        const std::array<std::string, 3> segmentRanges = {
            "20 <= phiDeg < 140; phiLocal = phiDeg - 20",
            "150 <= phiDeg < 250; phiLocal = phiDeg - 150",
            "260 <= phiDeg < 360 or 0 <= phiDeg < 10; phiLocal = (phiDeg < 10 ? phiDeg + 360 : phiDeg) - 260"};
        for (size_t iseg = 0; iseg < 3; ++iseg) {
          totalFitOut << "phiSeg" << iseg + 1 << ": " << segmentRanges[iseg] << "\n";
          totalFitOut << "  a0_phi(p) = " << cdPhiA0MomentumCoefficients[iseg][0]
                      << " + (" << cdPhiA0MomentumCoefficients[iseg][1] << ") * p\n";
          totalFitOut << "  a1_phi(p) = " << cdPhiA1MomentumCoefficients[iseg][0]
                      << " + (" << cdPhiA1MomentumCoefficients[iseg][1] << ") * p\n";
        }
        totalFitOut << "\na0_theta(pPhi) = " << thetaA0It->second[0]
                    << " + (" << thetaA0It->second[1] << ") * pPhi\n";
        totalFitOut << "a1_theta(pPhi) = " << thetaA1It->second[0]
                    << " + (" << thetaA1It->second[1] << ") * pPhi\n";
      }
    }
  }

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
  /// get mean values of Q^2 and x_B
  std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> getMeanQ2xBt_old(const BinManager& bins, std::unique_ptr<DISANAplotter>& plotter) {
    const auto& q2_bins = bins.GetQ2Bins();
    const auto& t_bins = bins.GetTBins();

    size_t n_xb = bins.GetMaxXBBinCount();
    size_t n_q2 = q2_bins.size() - 1;
    size_t n_t = t_bins.size() - 1;

    auto rdf = plotter->GetRDF();

    std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> result(
        n_xb, std::vector<std::vector<std::tuple<double, double, double>>>(n_q2, std::vector<std::tuple<double, double, double>>(n_t)));

    for (size_t ix = 0; ix < n_xb; ++ix) {
      for (size_t iq = 0; iq < n_q2; ++iq) {
        const auto& xb_bins = bins.GetXBBins(iq);
        if (ix + 1 >= xb_bins.size()) continue;
        for (size_t it = 0; it < n_t; ++it) {
          double xb_lo = xb_bins[ix], xb_hi = xb_bins[ix + 1];
          double q2_lo = q2_bins[iq], q2_hi = q2_bins[iq + 1];
          double t_lo = t_bins[it], t_hi = t_bins[it + 1];

          // Apply filter
          auto rdf_cut = rdf.Filter(Form("xB >= %f && xB < %f", xb_lo, xb_hi)).Filter(Form("Q2 >= %f && Q2 < %f", q2_lo, q2_hi)).Filter(Form("t >= %f && t < %f", t_lo, t_hi));

          // Compute means
          double mean_xB = rdf_cut.Mean("xB").GetValue();
          double mean_Q2 = rdf_cut.Mean("Q2").GetValue();
          double mean_t = rdf_cut.Mean("t").GetValue();

          result[ix][iq][it] = std::make_tuple(mean_xB, mean_Q2, mean_t);
        }
      }
    }

    return result;
  }

  std::vector<std::vector<std::vector<std::tuple<double, double, double>>>>
  getMeanQ2xBt(const BinManager& bins, std::unique_ptr<DISANAplotter>& plotter) {
    if (bins.HasCustomBins()) {
      struct Accumulator {
        double sum_xB = 0.0;
        double sum_Q2 = 0.0;
        double sum_t = 0.0;
        unsigned long long count = 0;
      };

      const auto& customBins = bins.GetCustomBins();
      auto rdf = plotter->GetRDF();
      const unsigned int nSlots = rdf.GetNSlots();
      std::vector<std::vector<Accumulator>> slotAccumulators(
          nSlots, std::vector<Accumulator>(customBins.size()));

      rdf.ForeachSlot(
          [&slotAccumulators, customBins](unsigned int slot, double xB, double Q2, double t) {
            for (size_t ib = 0; ib < customBins.size(); ++ib) {
              const auto& bin = customBins[ib];
              if (xB >= bin.xBMin && xB < bin.xBMax &&
                  Q2 >= bin.Q2Min && Q2 < bin.Q2Max &&
                  t >= bin.tMin && t < bin.tMax) {
                auto& acc = slotAccumulators[slot][ib];
                acc.sum_xB += xB;
                acc.sum_Q2 += Q2;
                acc.sum_t += t;
                ++acc.count;
              }
            }
          },
          {bins.GetFirstColumn(), bins.GetSecondColumn(), "t"});

      std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> result(
          1,
          std::vector<std::vector<std::tuple<double, double, double>>>(
              customBins.size(),
              std::vector<std::tuple<double, double, double>>(
                  1, std::make_tuple(0.0, 0.0, 0.0))));

      for (size_t ib = 0; ib < customBins.size(); ++ib) {
        Accumulator total;
        for (unsigned int slot = 0; slot < nSlots; ++slot) {
          const auto& acc = slotAccumulators[slot][ib];
          total.sum_xB += acc.sum_xB;
          total.sum_Q2 += acc.sum_Q2;
          total.sum_t += acc.sum_t;
          total.count += acc.count;
        }
        if (total.count > 0) {
          result[0][ib][0] = std::make_tuple(
              total.sum_xB / total.count,
              total.sum_Q2 / total.count,
              total.sum_t / total.count);
        }
      }
      return result;
    }

    const auto& q2_bins = bins.GetQ2Bins();
    const auto& t_bins  = bins.GetTBins();
  
    const size_t n_xb = bins.GetMaxXBBinCount();
    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t  = t_bins.size() - 1;
  
    auto rdf = plotter->GetRDF();
  
    // 保持原来的返回结构
    std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> result(
        n_xb,
        std::vector<std::vector<std::tuple<double, double, double>>>(
            n_q2,
            std::vector<std::tuple<double, double, double>>(n_t, std::make_tuple(0.0, 0.0, 0.0))
        )
    );
  
    struct Accumulator {
      double sum_xB = 0.0;
      double sum_Q2 = 0.0;
      double sum_t  = 0.0;
      unsigned long long count = 0;
    };
  
    auto findBin = [](double x, const std::vector<double>& edges) -> int {
      // 与你原来的条件一致：lo <= x < hi
      auto it = std::upper_bound(edges.begin(), edges.end(), x);
      if (it == edges.begin() || it == edges.end()) return -1;
      return static_cast<int>(std::distance(edges.begin(), it) - 1);
    };
  
    const size_t nCells = n_xb * n_q2 * n_t;
    const unsigned int nSlots = rdf.GetNSlots();
  
    // 每个 slot 一份局部缓存，避免多线程竞争
    std::vector<std::vector<Accumulator>> slotAccs(
        nSlots, std::vector<Accumulator>(nCells)
    );
  
    auto rdf_binned = rdf
      .Define("iq_bin", [q2_bins, findBin](double Q2) {
        return findBin(Q2, q2_bins);
      }, {bins.GetSecondColumn()})
      .Define("ix_bin", [bins, findBin](double xB, int iq) {
        if (iq < 0) return -1;
        return findBin(xB, bins.GetXBBins(static_cast<size_t>(iq)));
      }, {bins.GetFirstColumn(), "iq_bin"})
      .Define("it_bin", [t_bins, findBin](double t) {
        return findBin(t, t_bins);
      }, {"t"})
      .Filter([](int ix, int iq, int it) {
        return ix >= 0 && iq >= 0 && it >= 0;
      }, {"ix_bin", "iq_bin", "it_bin"});
  
    rdf_binned.ForeachSlot(
        [&slotAccs, n_q2, n_t](unsigned int slot, double xB, double Q2, double t,
                               int ix, int iq, int it) {
          const size_t idx = static_cast<size_t>(ix) * n_q2 * n_t
                           + static_cast<size_t>(iq) * n_t
                           + static_cast<size_t>(it);
  
          auto& acc = slotAccs[slot][idx];
          acc.sum_xB += xB;
          acc.sum_Q2 += Q2;
          acc.sum_t  += t;
          ++acc.count;
        },
        {bins.GetFirstColumn(), bins.GetSecondColumn(), "t", "ix_bin", "iq_bin", "it_bin"}
    );
  
    // 合并各个 slot 的结果
    std::vector<Accumulator> totalAccs(nCells);
    for (unsigned int slot = 0; slot < nSlots; ++slot) {
      for (size_t idx = 0; idx < nCells; ++idx) {
        totalAccs[idx].sum_xB += slotAccs[slot][idx].sum_xB;
        totalAccs[idx].sum_Q2 += slotAccs[slot][idx].sum_Q2;
        totalAccs[idx].sum_t  += slotAccs[slot][idx].sum_t;
        totalAccs[idx].count  += slotAccs[slot][idx].count;
      }
    }
  
    // 写回原来的 3D tuple 结构
    for (size_t ix = 0; ix < n_xb; ++ix) {
      for (size_t iq = 0; iq < n_q2; ++iq) {
        for (size_t it = 0; it < n_t; ++it) {
          const size_t idx = ix * n_q2 * n_t + iq * n_t + it;
          const auto& acc = totalAccs[idx];
  
          if (acc.count > 0) {
            result[ix][iq][it] = std::make_tuple(
                acc.sum_xB / acc.count,
                acc.sum_Q2 / acc.count,
                acc.sum_t  / acc.count
            );
          } else {
            // 空 bin 时保持默认值
            result[ix][iq][it] = std::make_tuple(0.0, 0.0, 0.0);
            // 如果你更想看 NaN，可以改成：
            // const double nan = std::numeric_limits<double>::quiet_NaN();
            // result[ix][iq][it] = std::make_tuple(nan, nan, nan);
          }
        }
      }
    }
  
    return result;
  }

  // Plot all basic kinematic distributions (p, theta, phi) for all particle types
  void PlotKinematicComparison_phiAna() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(4, 4);

    std::vector<std::string> types = {"el", "pro", "kPlus", "kMinus"};
    std::vector<std::string> vars = {"p", "theta", "phi", "vz"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotVariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "KinematicComparison_phiAna.pdf").c_str());

    // Optionally save individual plots
    if (plotIndividual) {
      for (const auto& type : types) {
        for (const auto& var : vars) {
          PlotSingleVariableComparison(type, var);
        }
      }
    }

    std::cout << "Saved kinematic comparison plots to: " << outputDir + "/KinematicComparison_phiAna.pdf" << std::endl;
    delete canvas;
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
      // styleKin_.StyleTH1(target, typeToParticle[type].c_str());
      auto [cr, cg, cb] = modelShades[i % modelShades.size()];
      const int colorIdx = 4000 + int(i) * 20;
      if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
      target->SetMarkerColor(colorIdx);
      target->SetLineColorAlpha(colorIdx, 0.8);
      target->SetLineWidth(1);
      // target->SetTitleSize(0.09);
      target->SetTitle(typeToParticle[type].c_str());
      target->GetXaxis()->SetTitle(VarName[var].c_str());
      target->GetYaxis()->SetTitle("Count");
      gStyle->SetTitleFillColor(0);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetOptTitle(1);
      // After the first draw, remove the title box
      gStyle->SetTitleX(0.5);
      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw("SAME");
  }

  void PlotPi0KinematicComparison() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 4);

    std::vector<std::string> types = {"el", "pro", "pho", "pho2"};
    std::vector<std::string> vars = {"p", "theta", "phi"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotPi0VariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "Pi0KinematicComparison.pdf").c_str());

    std::cout << "Saved pi0 kinematic comparison plots to: " << outputDir + "/Pi0KinematicComparison.pdf" << std::endl;
    delete canvas;
  }

  // Plot one specific variable (e.g., p) for a given particle type (e.g., "el")
  void PlotPi0VariableComparison(const std::string& type, const std::string& var, int padIndex, TCanvas* canvas) {
    canvas->cd(padIndex);
    std::string hname_target = "pi0_rec" + type + "_" + var;

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
      auto [cr, cg, cb] = modelShades[i % modelShades.size()];
      const int colorIdx = 4000 + int(i) * 20;
      if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
      target->SetMarkerColor(colorIdx);
      target->SetLineColorAlpha(colorIdx, 0.8);
      target->SetLineWidth(1);
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
      // target->SetTitle(Form("%s;%s;Count", typeToParticle[type].c_str()), VarName[var].c_str());

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

  void PlotKinematic1D() {
    if (plotters.empty()) {
      std::cerr << "PlotKinematic1D requires at least one loaded model.\n";
      return;
    }

    struct KinSpec {
      std::string var;
      std::string title;
      std::string xlabel;
      int bins;
      double xmin;
      double xmax;
    };

    const std::vector<KinSpec> vars = {
        {"Q2", "Q^{2}", "Q^{2} [GeV^{2}]", 120, 1.0, 5.0},
        {"xB", "x_{B}", "x_{B}", 120, 0.1, 0.5},
        {"t", "-t", "-t [GeV^{2}]", 120, 0.1, 0.8},
        {"phi", "#phi", "#phi [deg]", 120, 0.0, 360.0}};

    struct BookedKinHistograms {
      size_t plotterIndex;
      size_t variableIndex;
      ROOT::RDF::RResultPtr<TH1D> data;
      ROOT::RDF::RResultPtr<TH1D> mc;
    };

    std::vector<BookedKinHistograms> bookedHistograms;
    std::vector<ROOT::RDF::RResultHandle> resultHandles;

    for (size_t m = 0; m < plotters.size(); ++m) {
      auto rdfData = plotters[m]->GetRDF();
      auto rdfMC = plotters[m]->GetRDF_DVCSMC();

      for (size_t i = 0; i < vars.size(); ++i) {
        const auto& spec = vars[i];
        if (!rdfData.HasColumn(spec.var) || !rdfMC.HasColumn(spec.var)) {
          std::cerr << "[PlotKinematic1D] skipping " << spec.var
                    << " for " << labels[m] << ": missing column.\n";
          continue;
        }

        BookedKinHistograms booked{
            m,
            i,
            rdfData.Histo1D(
                {Form("h_kin1d_data_%s_%zu", spec.var.c_str(), m),
                 (spec.title + ";" + spec.xlabel + ";Normalized counts").c_str(),
                 spec.bins, spec.xmin, spec.xmax},
                spec.var),
            rdfMC.Histo1D(
                {Form("h_kin1d_dvcsmc_%s_%zu", spec.var.c_str(), m),
                 (spec.title + ";" + spec.xlabel + ";Normalized counts").c_str(),
                 spec.bins, spec.xmin, spec.xmax},
                spec.var)};

        resultHandles.emplace_back(booked.data);
        resultHandles.emplace_back(booked.mc);
        bookedHistograms.emplace_back(std::move(booked));
      }
    }

    if (resultHandles.empty()) {
      std::cerr << "[PlotKinematic1D] no histograms were booked.\n";
      return;
    }

    const unsigned int eventLoopCount = ROOT::RDF::RunGraphs(resultHandles);
    std::cout << "[PlotKinematic1D] completed " << eventLoopCount
              << " independent RDataFrame event loop(s).\n";

    auto normalizeVisibleRange = [](TH1* hist) {
      if (!hist) return;
      const double integral = hist->Integral(1, hist->GetNbinsX());
      if (integral > 0.0) hist->Scale(1.0 / integral);
    };

    for (size_t i = 0; i < vars.size(); ++i) {
      const auto& spec = vars[i];
      TCanvas* canvas = new TCanvas(
          Form("c_Kinematic1D_%s", spec.var.c_str()),
          Form("DVCS data vs MC %s", spec.var.c_str()), 900, 700);
      styleKin_.StylePad((TPad*)gPad);
      gPad->SetTicks();

      TLegend* legend = new TLegend(0.58, 0.70, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.035);

      bool first = true;
      TH1D* frameHistogram = nullptr;
      double maxY = 0.0;

      for (auto& booked : bookedHistograms) {
        if (booked.variableIndex != i) continue;

        TH1D* hData = (TH1D*)booked.data.GetPtr()->Clone();
        hData->SetDirectory(0);
        TH1D* hMC = (TH1D*)booked.mc.GetPtr()->Clone();
        hMC->SetDirectory(0);

        normalizeVisibleRange(hData);
        normalizeVisibleRange(hMC);

        styleKin_.StyleTH1(hData);
        hData->SetLineColor(kRed + 1);
        hData->SetLineWidth(1);
        hData->SetLineStyle(1);

        styleKin_.StyleTH1(hMC);
        hMC->SetLineColor(kGreen + 1);
        hMC->SetLineWidth(1);
        hMC->SetLineStyle(1);

        if (first) {
          hData->Draw("HIST");
          frameHistogram = hData;
          first = false;
        } else {
          hData->Draw("HIST SAME");
        }
        hMC->Draw("HIST SAME");

        maxY = std::max(maxY, hData->GetMaximum());
        maxY = std::max(maxY, hMC->GetMaximum());
        const std::string& label = labels[booked.plotterIndex];
        legend->AddEntry(hData, (label + " Data").c_str(), "l");
        legend->AddEntry(hMC, (label + " MC").c_str(), "l");
      }

      if (frameHistogram) {
        frameHistogram->SetTitle("");
        frameHistogram->GetXaxis()->SetTitle(spec.xlabel.c_str());
        frameHistogram->GetYaxis()->SetTitle("Normalized counts");
        frameHistogram->SetMaximum(maxY * 1.2);
      }
      legend->Draw();
      gPad->Modified();
      gPad->Update();

      const std::string outpath = outputDir + "/Kinematic1D_" + spec.var + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved DVCS kinematic 1D comparison to: " << outpath << "\n";
      delete legend;
      delete canvas;
    }
  }

  void PlotKinematic2D() {
    if (plotters.empty()) {
      std::cerr << "PlotKinematic2D requires at least one loaded model.\n";
      return;
    }

    struct Kin2DSpec {
      std::string tag;
      std::string xvar;
      std::string yvar;
      std::string xlabel;
      std::string ylabel;
      int xbins;
      double xmin;
      double xmax;
      int ybins;
      double ymin;
      double ymax;
    };

    const std::vector<Kin2DSpec> vars = {
        {"xB_vs_Q2", "xB", "Q2", "x_{B}", "Q^{2} [GeV^{2}]",
         160, 0.1, 0.7, 160, 1.0, 6.0},
        {"phi_vs_t", "phi", "t", "#phi [deg]", "-t [GeV^{2}]",
         180, 0.0, 360.0, 160, 0.0, 2.0}};

    struct BookedKin2DHistogram {
      size_t plotterIndex;
      size_t variableIndex;
      ROOT::RDF::RResultPtr<TH2D> data;
    };

    std::vector<BookedKin2DHistogram> bookedHistograms;
    std::vector<ROOT::RDF::RResultHandle> resultHandles;

    for (size_t m = 0; m < plotters.size(); ++m) {
      auto rdfData = plotters[m]->GetRDF();

      for (size_t i = 0; i < vars.size(); ++i) {
        const auto& spec = vars[i];
        if (!rdfData.HasColumn(spec.xvar) || !rdfData.HasColumn(spec.yvar)) {
          std::cerr << "[PlotKinematic2D] skipping " << spec.tag
                    << " for " << labels[m] << ": missing column.\n";
          continue;
        }

        BookedKin2DHistogram booked{
            m,
            i,
            rdfData.Histo2D(
                {Form("h_kin2d_data_%s_%zu", spec.tag.c_str(), m),
                 (";" + spec.xlabel + ";" + spec.ylabel + ";Counts").c_str(),
                 spec.xbins, spec.xmin, spec.xmax,
                 spec.ybins, spec.ymin, spec.ymax},
                spec.xvar, spec.yvar)};

        resultHandles.emplace_back(booked.data);
        bookedHistograms.emplace_back(std::move(booked));
      }
    }

    if (resultHandles.empty()) {
      std::cerr << "[PlotKinematic2D] no histograms were booked.\n";
      return;
    }

    const unsigned int eventLoopCount = ROOT::RDF::RunGraphs(resultHandles);
    std::cout << "[PlotKinematic2D] completed " << eventLoopCount
              << " independent RDataFrame event loop(s).\n";

    for (size_t i = 0; i < vars.size(); ++i) {
      const auto& spec = vars[i];
      for (auto& booked : bookedHistograms) {
        if (booked.variableIndex != i) continue;

        const std::string& label = labels[booked.plotterIndex];
        TCanvas* canvas = new TCanvas(
            Form("c_Kinematic2D_%s_%zu", spec.tag.c_str(), booked.plotterIndex),
            Form("%s %s", label.c_str(), spec.tag.c_str()), 900, 750);
        styleKin_.StylePad((TPad*)gPad);
        gPad->SetTicks();
        gPad->SetRightMargin(0.16);

        TH2D* hist = (TH2D*)booked.data.GetPtr()->Clone();
        hist->SetDirectory(0);
        hist->SetStats(0);
        hist->SetTitle(label.c_str());
        hist->GetXaxis()->SetTitle(spec.xlabel.c_str());
        hist->GetYaxis()->SetTitle(spec.ylabel.c_str());
        hist->GetZaxis()->SetTitle("Counts");
        hist->Draw("COLZ");
        gPad->Modified();
        gPad->Update();

        std::string outpath = outputDir + "/Kinematic2D_" + spec.tag;
        if (plotters.size() > 1) outpath += "_" + label;
        std::replace(outpath.begin(), outpath.end(), ' ', '_');
        std::replace(outpath.begin(), outpath.end(), ',', '_');
        outpath += ".pdf";

        canvas->SaveAs(outpath.c_str());
        std::cout << "Saved DVCS kinematic 2D plot to: " << outpath << "\n";
        delete hist;
        delete canvas;
      }
    }
  }

  void Plot2DParticle() {
    if (plotters.empty()) {
      std::cerr << "Plot2DParticle requires at least one loaded model.\n";
      return;
    }

    struct ParticleSpec {
      std::string tag;
      std::string title;
      std::string pColumn;
      std::string thetaColumn;
      double thetaMin;
      double thetaMax;
      double pMin;
      double pMax;
    };

    const std::vector<ParticleSpec> particles = {
        {"e", "e", "recel_p", "recel_theta", 5.0, 40.0, 1.0, 7.0},
        {"p", "p", "recpro_p", "recpro_theta", 0.0, 80.0, 0.0, 2.0},
        {"gamma", "#gamma", "recpho_p", "recpho_theta", 0.0, 40.0, 1.0, 7.0}};

    struct BookedParticle2DHistogram {
      size_t plotterIndex;
      size_t particleIndex;
      ROOT::RDF::RResultPtr<TH2D> data;
    };

    std::vector<BookedParticle2DHistogram> bookedHistograms;
    std::vector<ROOT::RDF::RResultHandle> resultHandles;

    for (size_t m = 0; m < plotters.size(); ++m) {
      auto rdfData = plotters[m]->GetRDF();

      for (size_t i = 0; i < particles.size(); ++i) {
        const auto& spec = particles[i];
        if (!rdfData.HasColumn(spec.pColumn) || !rdfData.HasColumn(spec.thetaColumn)) {
          std::cerr << "[Plot2DParticle] skipping " << spec.tag
                    << " for " << labels[m] << ": missing column.\n";
          continue;
        }

        const std::string thetaDegColumn =
            Form("plot2dparticle_%s_theta_deg_%zu", spec.tag.c_str(), m);
        auto rdfParticle = rdfData.Define(
            thetaDegColumn,
            [](double theta) { return theta * 180.0 / M_PI; },
            {spec.thetaColumn});
        if (spec.tag == "e") {
          rdfParticle = rdfParticle.Filter(
              [](double p) { return p >= 2.0; },
              {spec.pColumn},
              "electron p >= 2 GeV/c for Plot2DParticle");
        }

        BookedParticle2DHistogram booked{
            m,
            i,
            rdfParticle.Histo2D(
                {Form("h_particle2d_%s_%zu", spec.tag.c_str(), m),
                 (spec.title + ";#theta [deg];p [GeV/c];Counts").c_str(),
                 160, spec.thetaMin, spec.thetaMax,
                 160, spec.pMin, spec.pMax},
                thetaDegColumn, spec.pColumn)};

        resultHandles.emplace_back(booked.data);
        bookedHistograms.emplace_back(std::move(booked));
      }
    }

    if (resultHandles.empty()) {
      std::cerr << "[Plot2DParticle] no histograms were booked.\n";
      return;
    }

    const unsigned int eventLoopCount = ROOT::RDF::RunGraphs(resultHandles);
    std::cout << "[Plot2DParticle] completed " << eventLoopCount
              << " independent RDataFrame event loop(s).\n";

    for (size_t i = 0; i < particles.size(); ++i) {
      const auto& spec = particles[i];
      for (auto& booked : bookedHistograms) {
        if (booked.particleIndex != i) continue;

        const std::string& label = labels[booked.plotterIndex];
        TCanvas* canvas = new TCanvas(
            Form("c_Particle2D_%s_%zu", spec.tag.c_str(), booked.plotterIndex),
            Form("%s %s theta vs p", label.c_str(), spec.tag.c_str()), 900, 750);
        styleKin_.StylePad((TPad*)gPad);
        gPad->SetTicks();
        gPad->SetRightMargin(0.16);

        TH2D* hist = (TH2D*)booked.data.GetPtr()->Clone();
        hist->SetDirectory(0);
        hist->SetStats(0);
        hist->SetTitle(label.c_str());
        hist->GetXaxis()->SetTitle("#theta [deg]");
        hist->GetYaxis()->SetTitle("p [GeV/c]");
        hist->GetZaxis()->SetTitle("Counts");
        hist->Draw("COLZ");
        gPad->Modified();
        gPad->Update();

        std::string cleanLabel = label;
        std::replace(cleanLabel.begin(), cleanLabel.end(), ' ', '_');
        std::replace(cleanLabel.begin(), cleanLabel.end(), ',', '_');

        std::string outpath = outputDir + "/Particle2D_" + spec.tag;
        if (plotters.size() > 1) outpath += "_" + cleanLabel;
        outpath += ".pdf";

        canvas->SaveAs(outpath.c_str());
        std::cout << "Saved DVCS particle 2D plot to: " << outpath << "\n";
        delete hist;
        delete canvas;
      }
    }
  }

  void PlotPhiDVEPKinematicsPlots(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();
    // variables to draw (9 one-dimensional plots)
    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi", "mtprime", "tmin", "cos_thetaKK", "cos_phiKK"};

    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"},
                                                 {"xB", "x_{B}"},
                                                 {"t", "-t [GeV^{2}]"},
                                                 {"mtprime", "-t' \\equiv (t_{min}-t) [GeV^{2}]"},
                                                 {"tmin", "t_{min} [GeV^{2}]"},
                                                 {"W", "W [GeV]"},
                                                 {"phi", "#phi [deg]"},
                                                 {"cos_thetaKK", "cos(#theta_{K^{+}K^{-}})"},
                                                 {"cos_phiKK", "cos(#phi_{K^{+}K^{-}})"}};

    // First pass: compute global [min,max] per variable across all plotters for consistent binning
    struct Range {
      double lo{+1e300}, hi{-1e300};
    };
    std::map<std::string, Range> globalRange;

    for (const auto& var : variables) {
      Range R;
      for (size_t i = 0; i < plotters.size(); ++i) {
        auto r = plotters[i]->GetRDF();
        const std::string col = var;
        if (!r.HasColumn(col)) continue;

        double lo = *(r.Min(col));
        double hi = *(r.Max(col));
        if (std::isfinite(lo) && std::isfinite(hi)) {
          R.lo = std::min(R.lo, lo);
          R.hi = std::max(R.hi, hi);
        }
      }
      if (R.lo > R.hi) {  // fallback if empty
        R.lo = 0.0;
        R.hi = 1.0;
      }
      // pad a margin
      double margin = std::max(1e-3, 0.05 * (R.hi - R.lo));
      globalRange[var].lo = R.lo - margin;
      globalRange[var].hi = R.hi + margin;
    }

    // Canvas: 3x4 = 12 pads → 9×1D + 3×2D
    TCanvas* canvas = new TCanvas("PhiEPVars", "Phi electroproduction — kinematic comparison", 1800, 1400);
    canvas->Divide(3, 4);

    int pad = 1;

    // ---------------------------------------------------------------------------
    //  9 × 1D plots
    // ---------------------------------------------------------------------------
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.60, 0.55, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.04);

      std::vector<TH1D*> histos_to_draw;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto r = plotters[i]->GetRDF();
        if (!r.HasColumn(var)) {
          std::cerr << "[WARN] Column " << var << " not present after derivation for model " << labels[i] << "\n";
          continue;
        }

        const auto& R = globalRange[var];
        const int nbins = 100;
        const std::string hname = Form("h_%s_%zu_%u", var.c_str(), i, gRandom->Integer(1u << 30));

        auto htmp = r.Histo1D({hname.c_str(), titles[var].c_str(), nbins, R.lo, R.hi}, var);
        auto* h = (TH1D*)htmp->Clone((hname + "_clone").c_str());
        if (!h) continue;

        h->SetDirectory(0);
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);

        auto [cr, cg, cb] = modelShades[i % modelShades.size()];
        const int colorIdx = 4000 + int(i) * 20;
        if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
        h->SetMarkerColor(colorIdx);
        h->SetLineColorAlpha(colorIdx, 0.95);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Normalized counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }
      if (!histos_to_draw.empty()) legend->Draw();

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi" || var == "mtprime" || var == "tmin")) {  // <- fixed name
        PlotSingleVariableComparison("el", var);
      }
    }

    // For 2D plots we just use the first plotter (change if you need overlays)
    auto r2d = plotters.front()->GetRDF();

    // Common style lambda
    auto style_h2d = [&]() {
      styleDVCS_.StylePad((TPad*)gPad);
      gPad->SetRightMargin(0.16);
      TH2* h = (TH2*)gPad->GetListOfPrimitives()->At(0);
      if (!h) return;
      h->SetStats(0);
      h->SetTitle("");
      h->GetYaxis()->SetNoExponent(true);
      h->GetYaxis()->SetLabelFont(42);
      h->GetYaxis()->SetLabelSize(0.06);
      h->GetYaxis()->SetTitleOffset(1.0);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetNdivisions(410);
      h->GetXaxis()->SetTitleSize(0.065);
      h->GetXaxis()->SetLabelFont(42);
      h->GetXaxis()->SetLabelSize(0.06);
      h->GetXaxis()->SetTitleOffset(0.9);
      h->GetXaxis()->SetNdivisions(205);
      h->GetZaxis()->SetNdivisions(410);
      h->GetZaxis()->SetLabelSize(0.06);
      h->GetZaxis()->SetTitleOffset(1.5);
      h->GetZaxis()->SetTitleSize(0.06);
      TGaxis::SetMaxDigits(3);
    };

    // ---------------------------------------------------------------------------
    //  3 × 2D plots (pads 10, 11, 12)
    // ---------------------------------------------------------------------------

    // 1) Q2 vs -t' (mtprime)
    if (pad <= 12) {
      canvas->cd(pad++);
      auto h2d = r2d.Histo2D({"h_Q2_vs_tprime", "Q^{2} vs -t';-t' [GeV^{2}];Q^{2} [GeV^{2}]", 60, std::max(0.0, globalRange["mtprime"].lo), globalRange["mtprime"].hi, 60,
                              globalRange["Q2"].lo, globalRange["Q2"].hi},
                             "mtprime", "Q2");
      h2d->DrawCopy("COLZ");
      style_h2d();
    }

    // 2) Q2 vs W
    if (pad <= 12) {
      canvas->cd(pad++);
      auto h2d2 =
          r2d.Histo2D({"h_Q2_vs_W", "Q^{2} vs W;W [GeV];Q^{2} [GeV^{2}]", 60, globalRange["W"].lo, globalRange["W"].hi, 60, globalRange["Q2"].lo, globalRange["Q2"].hi}, "W", "Q2");
      h2d2->DrawCopy("COLZ");
      style_h2d();
    }

    // 3) Q2 vs tmin  (requested last plot)
    if (pad <= 12) {
      canvas->cd(pad++);
      auto h2d3 = r2d.Histo2D({"h_Q2_vs_tmin", "Q^{2} vs t_{min};t_{min} [GeV^{2}];Q^{2} [GeV^{2}]", 60, globalRange["tmin"].lo, globalRange["tmin"].hi, 60, globalRange["Q2"].lo,
                               globalRange["Q2"].hi},
                              "tmin", "Q2");
      h2d3->DrawCopy("COLZ");
      style_h2d();
    }

    // Final save and cleanup
    canvas->SaveAs((outputDir + "/PhiAna_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved phi electroproduction kinematics comparison to: " << (outputDir + "/PhiAna_Kinematics_Comparison.pdf") << std::endl;

    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void DrawCustomGrid(const std::vector<double>& xB_lines,
                    const std::vector<double>& Q2_lines,
                    double xBmin, double xBmax,
                    double Q2min, double Q2max,
                    int lineStyle=2, int lineWidth=1, int lineColor=kRed)
  {
    for (double x : xB_lines){
      if (x < xBmin || x > xBmax) continue;
      TLine* lv = new TLine(x, Q2min, x, Q2max);
      lv->SetLineStyle(lineStyle);
      lv->SetLineWidth(lineWidth);
      lv->SetLineColor(lineColor);
      lv->Draw("SAME");
    }
    for (double q2 : Q2_lines){
      if (q2 < Q2min || q2 > Q2max) continue;
      TLine* lh = new TLine(xBmin, q2, xBmax, q2);
      lh->SetLineStyle(lineStyle);
      lh->SetLineWidth(lineWidth);
      lh->SetLineColor(lineColor);
      lh->Draw("SAME");
    }
  }

  std::vector<double> GetAllConfiguredXBEdges() const {
    std::vector<double> edges;
    if (fXbins.HasCustomBins()) {
      for (const auto& bin : fXbins.GetCustomBins()) {
        edges.push_back(bin.xBMin);
        edges.push_back(bin.xBMax);
      }
      std::sort(edges.begin(), edges.end());
      edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
      return edges;
    }
    const auto& q2_edges = fXbins.GetQ2Bins();
    for (size_t iq = 0; iq + 1 < q2_edges.size(); ++iq) {
      const auto& row_edges = fXbins.GetXBBins(iq);
      edges.insert(edges.end(), row_edges.begin(), row_edges.end());
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
  }

  std::vector<double> GetAllConfiguredQ2Edges() const {
    if (!fXbins.HasCustomBins()) return fXbins.GetQ2Bins();
    std::vector<double> edges;
    for (const auto& bin : fXbins.GetCustomBins()) {
      edges.push_back(bin.Q2Min);
      edges.push_back(bin.Q2Max);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
  }

  std::vector<double> GetAllConfiguredTEdges() const {
    if (!fXbins.HasCustomBins()) return fXbins.GetTBins();
    std::vector<double> edges;
    for (const auto& bin : fXbins.GetCustomBins()) {
      edges.push_back(bin.tMin);
      edges.push_back(bin.tMax);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
  }

  void DrawConfiguredXBQ2Grid(double xBmin, double xBmax,
                              double Q2min, double Q2max,
                              int lineStyle=2, int lineWidth=1,
                              int lineColor=kRed) {
    if (fXbins.HasCustomBins()) {
      for (const auto& bin : fXbins.GetCustomBins()) {
        const double xlo = std::max(bin.xBMin, xBmin);
        const double xhi = std::min(bin.xBMax, xBmax);
        const double ylo = std::max(bin.Q2Min, Q2min);
        const double yhi = std::min(bin.Q2Max, Q2max);
        if (xlo >= xhi || ylo >= yhi) continue;
        for (const auto& segment : {
                 std::array<double, 4>{xlo, ylo, xhi, ylo},
                 std::array<double, 4>{xhi, ylo, xhi, yhi},
                 std::array<double, 4>{xhi, yhi, xlo, yhi},
                 std::array<double, 4>{xlo, yhi, xlo, ylo}}) {
          TLine* line = new TLine(segment[0], segment[1], segment[2], segment[3]);
          line->SetLineStyle(lineStyle);
          line->SetLineWidth(lineWidth);
          line->SetLineColor(lineColor);
          line->Draw("SAME");
        }
      }
      return;
    }

    const auto& q2_edges = fXbins.GetQ2Bins();
    for (size_t iq = 0; iq + 1 < q2_edges.size(); ++iq) {
      const double qlo = q2_edges[iq];
      const double qhi = q2_edges[iq + 1];
      const auto& xb_edges = fXbins.GetXBBins(iq);
      if (qhi < Q2min || qlo > Q2max || xb_edges.size() < 2) continue;

      const double ylo = std::max(qlo, Q2min);
      const double yhi = std::min(qhi, Q2max);
      for (double xb : xb_edges) {
        if (xb < xBmin || xb > xBmax) continue;
        TLine* line = new TLine(xb, ylo, xb, yhi);
        line->SetLineStyle(lineStyle);
        line->SetLineWidth(lineWidth);
        line->SetLineColor(lineColor);
        line->Draw("SAME");
      }

      const double xlo = std::max(xb_edges.front(), xBmin);
      const double xhi = std::min(xb_edges.back(), xBmax);
      if (xlo > xhi) continue;
      const double q_edges_to_draw[2] = {qlo, qhi};
      for (double q : q_edges_to_draw) {
        if (q < Q2min || q > Q2max) continue;
        TLine* line = new TLine(xlo, q, xhi, q);
        line->SetLineStyle(lineStyle);
        line->SetLineWidth(lineWidth);
        line->SetLineColor(lineColor);
        line->Draw("SAME");
      }
    }
  }

  void DrawVerticalBinLines(const std::vector<double>& edges,
                            double xmin, double xmax,
                            double ymin, double ymax,
                            int lineStyle=1, int lineWidth=1,
                            int lineColor=kRed) {
    for (double x : edges) {
      if (x < xmin || x > xmax) continue;
      TLine* line = new TLine(x, ymin, x, ymax);
      line->SetLineStyle(lineStyle);
      line->SetLineWidth(lineWidth);
      line->SetLineColor(lineColor);
      line->Draw("SAME");
    }
  }

  void WriteXBQ2tEventsCSV(const std::string& filename = "xBQ2t_events.csv") const {
    const std::string outPath = outputDir + "/" + filename;
    std::ofstream csv(outPath);
    if (!csv) {
      std::cerr << "[WriteXBQ2tEventsCSV] ERROR: cannot open output file: "
                << outPath << std::endl;
      return;
    }

    auto csvEscape = [](const std::string& value) {
      if (value.find_first_of(",\"\n") == std::string::npos) return value;
      std::string escaped = "\"";
      for (char c : value) {
        if (c == '"') escaped += "\"\"";
        else escaped += c;
      }
      escaped += "\"";
      return escaped;
    };

    csv << std::setprecision(12);
    csv << "model_index,model_label,event_index,xB,Q2,t\n";

    unsigned long long totalEvents = 0;
    for (size_t imodel = 0; imodel < plotters.size(); ++imodel) {
      auto rdf = plotters[imodel]->GetRDF();
      if (!rdf.HasColumn("xB") || !rdf.HasColumn("Q2") || !rdf.HasColumn("t")) {
        std::cerr << "[WriteXBQ2tEventsCSV] skipping model " << imodel
                  << ": missing xB, Q2, or t column." << std::endl;
        continue;
      }

      auto xBValues = rdf.Take<double>("xB");
      auto Q2Values = rdf.Take<double>("Q2");
      auto tValues = rdf.Take<double>("t");
      std::vector<ROOT::RDF::RResultHandle> resultHandles = {
          xBValues, Q2Values, tValues};
      ROOT::RDF::RunGraphs(resultHandles);

      const auto& xB = *xBValues;
      const auto& Q2 = *Q2Values;
      const auto& t = *tValues;
      const size_t nEvents = std::min({xB.size(), Q2.size(), t.size()});
      const std::string label =
          imodel < labels.size() ? csvEscape(labels[imodel]) : Form("model_%zu", imodel);

      for (size_t iev = 0; iev < nEvents; ++iev) {
        csv << imodel << "," << label << "," << iev << ","
            << xB[iev] << "," << Q2[iev] << "," << t[iev] << "\n";
      }
      totalEvents += nEvents;
    }

    std::cout << "Saved xB Q2 t event CSV to: " << outPath
              << " (" << totalEvents << " events)" << std::endl;
  }


  void PlotDVCSKinematicsComparison(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1400);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      bool first = true;
      std::vector<TH1D*> histos_to_draw;

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

        // Get histogram (RResultPtr) and clone it
        auto htmp = rdf.Histo1D({Form("h_%s_%zu", var.c_str(), i), titles[var].c_str(), 100, min - margin, max + margin}, var);
        auto h = (TH1D*)htmp->Clone(Form("h_%s_%zu_clone", var.c_str(), i));

        if (!h) continue;  // guard against failed clone

        h->SetDirectory(0);  // prevent ROOT from managing ownership
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);
        h->SetLineColor(i + 2);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }

      if (!histos_to_draw.empty()) {
        legend->Draw();
      }

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi")) {
        PlotSingleVariableComparison("el", var);
      }
    }
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 60, 0, 1.0, 60, 0, 10.0}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/DVCS_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved DVCS kinematics comparison to: " << outputDir + "/DVCS_Kinematics_Comparison.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotDVPi0KinematicsComparison(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1400);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      bool first = true;
      std::vector<TH1D*> histos_to_draw;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto rdf = plotters[i]->GetRDF_Pi0Data();
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

        // Get histogram (RResultPtr) and clone it
        auto htmp = rdf.Histo1D({Form("h_%s_%zu", var.c_str(), i), titles[var].c_str(), 100, min - margin, max + margin}, var);
        auto h = (TH1D*)htmp->Clone(Form("h_%s_%zu_clone", var.c_str(), i));

        if (!h) continue;  // guard against failed clone

        h->SetDirectory(0);  // prevent ROOT from managing ownership
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);
        h->SetLineColor(i + 2);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }

      if (!histos_to_draw.empty()) {
        legend->Draw();
      }

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi")) {
        PlotSingleVariableComparison("el", var);
      }
    }
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF_Pi0Data();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 60, 0, 1.0, 60, 0, 10.0}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/DVPi0_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved DVPi0 kinematics comparison to: " << outputDir + "/DVPi0_Kinematics_Comparison.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotxBQ2tBin(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF();
    const auto xB_lines = GetAllConfiguredXBEdges();
    const auto Q2_lines = GetAllConfiguredQ2Edges();
    const auto t_lines = GetAllConfiguredTEdges();
    double xBmin = 0.05, xBmax = 0.7;
    double Q2min = 0.5, Q2max = 7.0;
    double tmin = 0.0, tmax = 2.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawConfiguredXBQ2Grid(xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBin.pdf").c_str());
    std::cout << "Saved xBQ2tBin kinematics to: " << outputDir + "/xBQ2tBin.pdf" << std::endl;
    WriteXBQ2tEventsCSV();
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  // SDHEP analogue of PlotxBQ2tBin.  The first two pads show the variables
  // used for binning; qT is included as a hard-scale diagnostic.
  void PlotSDHEPPhaseSpace() {
    if (plotters.empty()) {
      std::cerr << "[PlotSDHEPPhaseSpace] No models loaded.\n";
      return;
    }
    auto rdf = plotters.front()->GetRDF();
    for (const char* col : {"xi_SDHEP", "costheta_SDHEP", "t", "qT_SDHEP"}) {
      if (!rdf.HasColumn(col)) {
        std::cerr << "[PlotSDHEPPhaseSpace] Missing column: " << col << "\n";
        return;
      }
    }

    const auto xiEdges = GetAllConfiguredXBEdges();
    const auto cosEdges = GetAllConfiguredQ2Edges();
    const auto tEdges = GetAllConfiguredTEdges();
    const double xiMin = xiEdges.empty() ? 0.0 : xiEdges.front();
    const double xiMax = xiEdges.empty() ? 0.7 : xiEdges.back();
    const double cosMin = -1.0;
    const double cosMax = 1.0;
    const double tMin = tEdges.empty() ? 0.0 : tEdges.front();
    const double tMax = tEdges.empty() ? 2.0 : tEdges.back();

    TCanvas* canvas = new TCanvas("SDHEPPhaseSpace", "SDHEP phase space", 5400, 1800);
    canvas->Divide(3, 1);
    auto style2D = [&](TH2* h) {
      styleDVCS_.StylePad((TPad*)gPad);
      gPad->SetRightMargin(0.16);
      h->SetStats(0);
      h->SetTitle("");
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetZaxis()->SetLabelSize(0.05);
      h->DrawCopy("COLZ");
    };

    canvas->cd(1);
    auto hXiCos = rdf.Histo2D(
        {"h_SDHEP_costheta_vs_xi", ";#xi;cos#theta_{SDHEP}",
         500, xiMin, xiMax, 500, cosMin, cosMax},
        "xi_SDHEP", "costheta_SDHEP");
    style2D(hXiCos.GetPtr());
    DrawConfiguredXBQ2Grid(xiMin, xiMax, cosMin, cosMax, 1, 2, kRed);

    canvas->cd(2);
    auto hTCos = rdf.Histo2D(
        {"h_SDHEP_costheta_vs_t", ";-t [GeV^{2}];cos#theta_{SDHEP}",
         500, tMin, tMax, 500, cosMin, cosMax},
        "t", "costheta_SDHEP");
    style2D(hTCos.GetPtr());
    DrawVerticalBinLines(tEdges, tMin, tMax, cosMin, cosMax, 1, 2, kRed);

    canvas->cd(3);
    auto hQTXi = rdf.Histo2D(
        {"h_SDHEP_qT_vs_xi", ";#xi;q_{T}^{SDHEP} [GeV]",
         500, xiMin, xiMax, 500, 0.0, 0.5 * std::sqrt(2.0 * 0.9382720813 * plotters.front()->GetBeamEnergy())},
        "xi_SDHEP", "qT_SDHEP");
    style2D(hQTXi.GetPtr());

    const std::string path = outputDir + "/SDHEP_xi_costheta_t_PhaseSpace.pdf";
    canvas->SaveAs(path.c_str());
    std::cout << "Saved SDHEP phase space to: " << path << std::endl;
    delete canvas;
  }

  void PlotxBQ2tBinMC(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF_DVCSMC();
    const auto xB_lines = GetAllConfiguredXBEdges();
    const auto Q2_lines = GetAllConfiguredQ2Edges();
    const auto t_lines = GetAllConfiguredTEdges();
    double xBmin = 0.05, xBmax = 0.7;
    double Q2min = 0.5, Q2max = 7.0;
    double tmin = 0.0, tmax = 2.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawConfiguredXBQ2Grid(xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBinMC.pdf").c_str());
    std::cout << "Saved xBQ2tBinMC kinematics to: " << outputDir + "/xBQ2tBinMC.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotxBQ2tBinPi0(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF_Pi0Data();
    const auto xB_lines = GetAllConfiguredXBEdges();
    const auto Q2_lines = GetAllConfiguredQ2Edges();
    const auto t_lines = GetAllConfiguredTEdges();
    double xBmin = 0.05, xBmax = 0.7;
    double Q2min = 0.5, Q2max = 7.0;
    double tmin = 0.0, tmax = 2.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawConfiguredXBQ2Grid(xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBinPi0.pdf").c_str());
    std::cout << "Saved xBQ2tBinPi0 kinematics to: " << outputDir + "/xBQ2tBinPi0.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }


  void PlotxBQ2tBinPi0MC(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF_Pi0MC();
    const auto xB_lines = GetAllConfiguredXBEdges();
    const auto Q2_lines = GetAllConfiguredQ2Edges();
    const auto t_lines = GetAllConfiguredTEdges();
    double xBmin = 0.05, xBmax = 0.7;
    double Q2min = 0.5, Q2max = 7.0;
    double tmin = 0.0, tmax = 2.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawConfiguredXBQ2Grid(xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawVerticalBinLines(t_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBinPi0MC.pdf").c_str());
    std::cout << "Saved xBQ2tBinPi0MC kinematics to: " << outputDir + "/xBQ2tBinPi0MC.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }
  /// For exclusivity cuts, you can use the following function to select one triplet
  void PlotExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", -0.6, 0.6},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.4},
        {"Theta_gamma_gamma", "#theta(#gamma, #vec{q})", "#theta_{#gamma#gamma'} [deg]", -2.0, 4.0},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
        {"Mx2_epg", "Missing Mass Squared (ep#gamma)", "MM^{2}(ep#gamma) [GeV^{2}]", -0.05, 0.05},
        {"Mx2_eg", "Invariant Mass (e#gamma)", "M^{2}(e#gamma) [GeV^{2}]", -0.5, 3},
        {"Theta_e_gamma", "Angle: e-#gamma", "#theta(e, #gamma) [deg]", 0.0, 60.0}};
        //{"DeltaE", "Energy Balance", "#DeltaE [GeV]", -1.0, 2.0}};

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
          // double x1 = mean - 3 * sigma;
          // double x2 = mean + 3 * sigma;

          // TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          // TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          // line1->SetLineColor(m + 2);
          // line2->SetLineColor(m + 2);
          // line1->SetLineStyle(2);  // Dashed
          // line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(3) << mean << ", #sigma = " << std::fixed << std::setprecision(3) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          // line1->Draw("SAME");
          // line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };

  void PlotExclusivityComparisonByDetectorCaseswithPi0(
      const std::vector<std::pair<std::string, std::string>>& detectorCuts,
      const std::vector<double>& tbin,
      bool drawPi0Subtracted = false) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", -0.6, 0.6},
        {"Mx2_eg", "Missing Mass Squared (e#gamma)", "MM^{2}(e#gamma) [GeV^{2}]", -0.5, 3.0},
        {"Mx2_epg", "Missing Mass Squared (ep#gamma)", "MM^{2}(ep#gamma) [GeV^{2}]", -0.05, 0.05}};

    if (tbin.size() < 2) {
      std::cerr << "PlotExclusivityComparisonByDetectorCaseswithPi0 requires at least two t-bin edges.\n";
      return;
    }

    struct BookedHistograms {
      ROOT::RDF::RResultPtr<TH1D> dvcs;
      ROOT::RDF::RResultPtr<TH1D> pi0Data;
      ROOT::RDF::RResultPtr<TH1D> dvcsMC;
      ROOT::RDF::RResultPtr<TH1D> dvcsPi0MC;
      ROOT::RDF::RResultPtr<TH1D> pi0Pi0MC;
    };
    using HistogramKey = std::tuple<size_t, size_t, size_t, size_t>;
    std::map<HistogramKey, BookedHistograms> bookedHistograms;
    std::vector<ROOT::RDF::RResultHandle> resultHandles;

    // Book all detector/t-bin branches first so each underlying data source is
    // scanned once when RunGraphs executes the complete graph.
    for (size_t detectorIndex = 0; detectorIndex < detectorCuts.size(); ++detectorIndex) {
      const auto& [cutExpr, cutLabel] = detectorCuts[detectorIndex];
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      for (size_t tbinIndex = 0; tbinIndex + 1 < tbin.size(); ++tbinIndex) {
        const double tLow = tbin[tbinIndex];
        const double tHigh = tbin[tbinIndex + 1];
        if (tHigh <= tLow) continue;

        const std::string tCut = Form("t >= %.12g && t < %.12g", tLow, tHigh);
        const std::string tLabel = Form("%.3f <= t < %.3f GeV^{2}", tLow, tHigh);

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdfDVCS = plotters[m]->GetRDF().Filter(cutExpr, cutLabel).Filter(tCut, tLabel);
          auto rdfPi0 = plotters[m]->GetRDF_Pi0Data().Filter(cutExpr, cutLabel).Filter(tCut, tLabel);
          auto rdfDVCSMC = plotters[m]->GetRDF_DVCSMC().Filter(cutExpr, cutLabel).Filter(tCut, tLabel);
          auto rdfDVCSPi0MC = plotters[m]->GetRDF_DVCSSelectedPi0MC().Filter(cutExpr, cutLabel).Filter(tCut, tLabel);
          auto rdfPi0Pi0MC = plotters[m]->GetRDF_Pi0SelectedPi0MC().Filter(cutExpr, cutLabel).Filter(tCut, tLabel);

          for (size_t variableIndex = 0; variableIndex < vars.size(); ++variableIndex) {
            const auto& [var, title, xlabel, xmin, xmax] = vars[variableIndex];
            if (!rdfDVCS.HasColumn(var) || !rdfPi0.HasColumn(var) ||
                !rdfDVCSMC.HasColumn(var) || !rdfDVCSPi0MC.HasColumn(var) ||
                !rdfPi0Pi0MC.HasColumn(var)) {
              continue;
            }

            BookedHistograms booked{
                rdfDVCS.Histo1D(
                    {Form("h_dvcs_%s_%s_tbin%zu_%zu", var.c_str(), cleanName.c_str(), tbinIndex, m),
                     (title + ";" + xlabel + ";Counts / N_{DVCS}").c_str(), 100, xmin, xmax},
                    var),
                rdfPi0.Histo1D(
                    {Form("h_pi0data_%s_%s_tbin%zu_%zu", var.c_str(), cleanName.c_str(), tbinIndex, m),
                     (title + ";" + xlabel + ";Counts / N_{DVCS}").c_str(), 100, xmin, xmax},
                    var),
                rdfDVCSMC.Histo1D(
                    {Form("h_dvcsmc_%s_%s_tbin%zu_%zu", var.c_str(), cleanName.c_str(), tbinIndex, m),
                     "", 100, xmin, xmax},
                    var),
                rdfDVCSPi0MC.Histo1D(
                    {Form("h_dvcs_pi0mc_%s_%s_tbin%zu_%zu", var.c_str(), cleanName.c_str(), tbinIndex, m),
                     "", 100, xmin, xmax},
                    var),
                rdfPi0Pi0MC.Histo1D(
                    {Form("h_pi0_pi0mc_%s_%s_tbin%zu_%zu", var.c_str(), cleanName.c_str(), tbinIndex, m),
                     "", 100, xmin, xmax},
                    var)};

            resultHandles.emplace_back(booked.dvcs);
            resultHandles.emplace_back(booked.pi0Data);
            resultHandles.emplace_back(booked.dvcsMC);
            resultHandles.emplace_back(booked.dvcsPi0MC);
            resultHandles.emplace_back(booked.pi0Pi0MC);
            bookedHistograms.emplace(
                HistogramKey{detectorIndex, tbinIndex, variableIndex, m},
                std::move(booked));
          }
        }
      }
    }

    unsigned int eventLoopCount = 0;
    if (!resultHandles.empty()) eventLoopCount = ROOT::RDF::RunGraphs(resultHandles);
    std::cout << "[PlotExclusivityComparisonByDetectorCaseswithPi0] completed "
              << eventLoopCount << " independent RDataFrame event loop(s).\n";

    for (size_t detectorIndex = 0; detectorIndex < detectorCuts.size(); ++detectorIndex) {
      const std::string& cutLabel = detectorCuts[detectorIndex].second;
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      for (size_t tbinIndex = 0; tbinIndex + 1 < tbin.size(); ++tbinIndex) {
        const double tLow = tbin[tbinIndex];
        const double tHigh = tbin[tbinIndex + 1];
        if (tHigh <= tLow) {
          std::cerr << "Skipping invalid t bin [" << tLow << ", " << tHigh << ").\n";
          continue;
        }

        const std::string tLabel = Form("%.3f <= t < %.3f GeV^{2}", tLow, tHigh);
        const std::string canvasName = Form("c_DVCS_vs_Pi0Data_%s_tbin%zu", cleanName.c_str(), tbinIndex);
        TCanvas* canvas = new TCanvas(canvasName.c_str(), (cutLabel + ", " + tLabel).c_str(), 1800, 600);
        const int cols = 3;
        const int rows = 1;
        canvas->Divide(cols, rows);

        for (size_t i = 0; i < vars.size(); ++i) {
          canvas->cd(i + 1);
          const std::string& var = std::get<0>(vars[i]);
          gPad->SetTicks();
          styleKin_.StylePad((TPad*)gPad);
          if (var == "Mx2_epg") {
            gPad->SetBottomMargin(0.16);
          }

        TLegend* legend = new TLegend(0.54, 0.68, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.036);

        TH1D* frameHistogram = nullptr;
        double maxY = 0.0;
        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          const auto bookedIt = bookedHistograms.find(
              HistogramKey{detectorIndex, tbinIndex, i, m});
          if (bookedIt == bookedHistograms.end()) continue;
          auto& booked = bookedIt->second;

          TH1D* hDVCS = (TH1D*)booked.dvcs.GetPtr()->Clone();
          hDVCS->SetDirectory(0);
          TH1D* hDVCSMC = (TH1D*)booked.dvcsMC.GetPtr()->Clone();
          hDVCSMC->SetDirectory(0);

          const double dvcsIntegral = hDVCS->Integral();
          if (dvcsIntegral <= 0.0) {
            delete hDVCS;
            delete hDVCSMC;
            continue;
          }
          const double dvcsMCIntegral = hDVCSMC->Integral();

          TH1D* hSubtracted = nullptr;
          if (drawPi0Subtracted) {
            const double nDVCSPi0MC = booked.dvcsPi0MC.GetPtr()->Integral();
            const double nPi0Pi0MC = booked.pi0Pi0MC.GetPtr()->Integral();
            const double pi0TransferFactor =
                nPi0Pi0MC > 0.0 ? nDVCSPi0MC / nPi0Pi0MC : 0.0;

            TH1D* hPi0Contamination = (TH1D*)booked.pi0Data.GetPtr()->Clone();
            hPi0Contamination->SetDirectory(0);
            hPi0Contamination->Scale(pi0TransferFactor);

            hSubtracted = (TH1D*)hDVCS->Clone(
                Form("h_pi0_subtracted_%s_%s_tbin%zu_%zu",
                     var.c_str(), cleanName.c_str(), tbinIndex, m));
            hSubtracted->SetDirectory(0);
            hSubtracted->Add(hPi0Contamination, -1.0);
            for (int bin = 0; bin <= hSubtracted->GetNbinsX() + 1; ++bin) {
              if (hSubtracted->GetBinContent(bin) < 0.0) {
                hSubtracted->SetBinContent(bin, 0.0);
              }
            }
            delete hPi0Contamination;
          }

          hDVCS->Scale(1.0 / dvcsIntegral);
          styleKin_.StyleTH1(hDVCS);
          hDVCS->SetLineColor(m + 2);
          hDVCS->SetLineWidth(1);
          hDVCS->SetLineStyle(drawPi0Subtracted ? 2 : 1);
          if (var == "Mx2_epg") {
            hDVCS->GetXaxis()->SetNdivisions(404);
          }

          if (hSubtracted) {
            hSubtracted->Scale(1.0 / dvcsIntegral);
            styleKin_.StyleTH1(hSubtracted);
            hSubtracted->SetLineColor(m + 2);
            hSubtracted->SetLineWidth(1);
            hSubtracted->SetLineStyle(1);
            if (var == "Mx2_epg") {
              hSubtracted->GetXaxis()->SetNdivisions(404);
            }
          }

          if (dvcsMCIntegral > 0.0) {
            hDVCSMC->Scale(1.0 / dvcsMCIntegral);
            styleKin_.StyleTH1(hDVCSMC);
            hDVCSMC->SetLineColor(kGreen + 1);
            hDVCSMC->SetLineWidth(1);
            hDVCSMC->SetLineStyle(1);
            if (var == "Mx2_epg") {
              hDVCSMC->GetXaxis()->SetNdivisions(404);
            }
          }

          if (first) {
            hDVCS->Draw("HIST");
            frameHistogram = hDVCS;
            first = false;
          } else {
            hDVCS->Draw("HIST SAME");
          }
          if (hSubtracted) hSubtracted->Draw("HIST SAME");
          if (dvcsMCIntegral > 0.0) hDVCSMC->Draw("HIST SAME");

          maxY = std::max(maxY, hDVCS->GetMaximum());
          if (hSubtracted) maxY = std::max(maxY, hSubtracted->GetMaximum());
          if (dvcsMCIntegral > 0.0) maxY = std::max(maxY, hDVCSMC->GetMaximum());
          legend->AddEntry(hDVCS, (labels[m] + " DVCS Data").c_str(), "l");
          if (hSubtracted) legend->AddEntry(hSubtracted, (labels[m] + " Pi0 Subtracted").c_str(), "l");
          if (dvcsMCIntegral > 0.0) legend->AddEntry(hDVCSMC, (labels[m] + " DVCS MC").c_str(), "l");
        }

          if (frameHistogram) frameHistogram->SetMaximum(maxY * 1.2);
          legend->Draw();
          gPad->Modified();
          gPad->Update();
        }

        std::string outpath = outputDir + Form("/Exclusivity_DVCS_vs_Pi0Data_%s_tbin%zu_%.3f_%.3f.pdf",
                                               cleanName.c_str(), tbinIndex, tLow, tHigh);
        canvas->SaveAs(outpath.c_str());
        std::cout << "Saved DVCS vs Pi0 data exclusivity comparison to: " << outpath << "\n";
        delete canvas;
      }
    }
  }

  void PlotPi0ExclusivityComparisonByDetectorCases(
      const std::vector<std::pair<std::string, std::string>>& detectorCuts,
      bool drawDVPi0MC = true,
      bool outputWideMpi0 = false) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mass_pi0", "Pi0 Mass", "M_{#pi0} [GeV]", 0.07, 0.2},
        {"Emiss_pi0", "Missing Energy", "E_{miss} [GeV]", -1.0, 1.0},
        {"PTmiss_pi0", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", 0.0, 0.3},
        {"Theta_pi0pi0", "#theta(#pi0, #pi')", "#theta_{#pi#pi'} [deg]", 0.0, 3.0},
        {"DeltaPhi_pi0", "Coplanarity Angle", "#Delta#phi [deg]", 0.0, 15.0},
        {"Mx2_eppi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep#pi) [GeV^{2}]", -0.03, 0.03},
        {"Mx2_ep_pi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep) [GeV^{2}]", -0.4, 0.6},
        {"Mx2_epi0", "Missing Mass Squared (e#pi)", "MM^{2}(e#pi) [GeV^{2}]", -0.2, 2.0},
        {"Theta_epho1", "Angle: e-#gamma_{1}", "#theta(e, #gamma_{1}) [deg]", 0.0, 60.0},
        {"Theta_epho2", "Angle: e-#gamma_{2}", "#theta(e, #gamma_{2}) [deg]", 0.0, 60.0},
        {"Theta_pho1pho2", "Angle: #gamma_{1}-#gamma_{2}", "#theta(#gamma_{1}, #gamma_{2}) [deg]", 0.0, 15.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 4800, 3200);
      int cols = 4;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];

        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.50, 0.70, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.040);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_data = plotters[m]->GetRDF_Pi0Data().Filter(cutExpr, cutLabel);
          if (!rdf_data.HasColumn(var)) continue;

          double xlo = xmin, xhi = xmax;

          auto h_data = rdf_data.Histo1D({Form("h_data_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xlo, xhi}, var);
          h_data.GetValue();  // materialize
          TH1D* hD = (TH1D*)h_data.GetPtr()->Clone();
          hD->SetDirectory(0);
          NormalizeHistogram(hD);
          styleKin_.StyleTH1(hD);
          hD->SetLineColor(m + 2);
          hD->SetLineWidth(2);
          hD->SetLineStyle(1);

          TH1D* hM = nullptr;
          if (drawDVPi0MC) {
            auto rdf_mc = plotters[m]->GetRDF_DVCSPi0MC().Filter(cutExpr, cutLabel);
            if (rdf_mc.HasColumn(var)) {
              auto h_mc = rdf_mc.Histo1D({Form("h_mc_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
              h_mc.GetValue();
              hM = (TH1D*)h_mc.GetPtr()->Clone();
              hM->SetDirectory(0);
              NormalizeHistogram(hM);
              styleKin_.StyleTH1(hM);
              hM->SetLineColor(m + 2);
              hM->SetLineWidth(2);
              hM->SetLineStyle(2);
            }
          }

          if (first) {
            hD->Draw("HIST");
            if (hM) hM->Draw("HIST SAME");
            first = false;
          } else {
            hD->Draw("HIST SAME");
            if (hM) hM->Draw("HIST SAME");
          }
          double maxY = hD->GetMaximum();
          if (hM) maxY = std::max(maxY, hM->GetMaximum());
          hD->SetMaximum(maxY * 1.2);
          gPad->Modified();
          gPad->Update();

          legend->AddEntry(hD, (labels[m] + " Data").c_str(), "l");
          if (hM) legend->AddEntry(hM, (labels[m] + " MC").c_str(), "l");

          double mean = hD->GetMean();
          double sigma = hD->GetStdDev();
          // double x1 = mean - 3.0 * sigma;
          // double x2 = mean + 3.0 * sigma;

          double meanM = hM ? hM->GetMean() : 0.0;
          double sigmaM = hM ? hM->GetStdDev() : 0.0;
          // double x1M = meanM - 3.0 * sigmaM;
          // double x2M = meanM + 3.0 * sigmaM;

          // double ymax = std::max(hD->GetMaximum(), hM ? hM->GetMaximum() : 0.0);
          // TLine* line1 = new TLine(x1, 0.0, x1, ymax * 0.5);
          // TLine* line2 = new TLine(x2, 0.0, x2, ymax * 0.5);
          // line1->SetLineColor(m + 2);
          // line2->SetLineColor(m + 2);
          // line1->SetLineStyle(3);
          // line2->SetLineStyle(3);
          // line1->Draw("SAME");
          // line2->Draw("SAME");

          // TLine* line1M = new TLine(x1M, 0.0, x1M, ymax * 0.5);
          // TLine* line2M = new TLine(x2M, 0.0, x2M, ymax * 0.5);
          // line1M->SetLineColor(m + 2);
          // line2M->SetLineColor(m + 2);
          // line1M->SetLineStyle(3);
          // line2M->SetLineStyle(3);
          // line1M->Draw("SAME");
          // line2M->Draw("SAME");

          std::ostringstream stats;
          stats << "Data #mu = " << std::fixed << std::setprecision(3) << mean << ", #sigma = " << std::fixed << std::setprecision(3) << sigma;
          legend->AddEntry((TObject*)nullptr, stats.str().c_str(), "");

          if (hM) {
            std::ostringstream statsM;
            statsM << "MC #mu = " << std::fixed << std::setprecision(3) << meanM << ", #sigma = " << std::fixed << std::setprecision(3) << sigmaM;
            legend->AddEntry((TObject*)nullptr, statsM.str().c_str(), "");
          }
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Pi0Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;

      if (outputWideMpi0) {
        TCanvas* mpi0Canvas = new TCanvas(("c_Mass_pi0_wide_" + cleanName).c_str(), cutLabel.c_str(), 1600, 1200);
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* mpi0Legend = new TLegend(0.50, 0.62, 0.88, 0.88);
        mpi0Legend->SetBorderSize(0);
        mpi0Legend->SetFillStyle(0);
        mpi0Legend->SetTextSize(0.032);

        TH1D* frameHistogram = nullptr;
        double maxY = 0.0;
        bool firstMpi0 = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_data = plotters[m]->GetRDF_Pi0Data().Filter(cutExpr, cutLabel);
          if (!rdf_data.HasColumn("Mass_pi0")) continue;

          auto h_data = rdf_data.Histo1D(
              {Form("h_data_Mass_pi0_wide_%s_%zu", cleanName.c_str(), m),
               "Wide Pi0 Mass;M_{#gamma#gamma} [GeV];Counts", 200, 0.0, 1.0},
              "Mass_pi0");
          h_data.GetValue();
          TH1D* hD = (TH1D*)h_data.GetPtr()->Clone();
          hD->SetDirectory(0);
          auto countInRange = [](const TH1D* hist, double xmin, double xmax) {
            const int firstBin = hist->GetXaxis()->FindFixBin(xmin + 1.0e-9);
            const int lastBin = hist->GetXaxis()->FindFixBin(xmax - 1.0e-9);
            return hist->Integral(firstBin, lastBin);
          };
          const double dataPi0Count = countInRange(hD, 0.05, 0.265);
          const double dataEtaCount = countInRange(hD, 0.4, 0.7);
          const double dataEtaPi0Ratio = dataPi0Count > 0.0 ? dataEtaCount / dataPi0Count : 0.0;
          NormalizeHistogram(hD);
          styleKin_.StyleTH1(hD);
          hD->GetXaxis()->SetTitleSize(0.040);
          hD->GetYaxis()->SetTitleSize(0.040);
          hD->GetXaxis()->SetLabelSize(0.034);
          hD->GetYaxis()->SetLabelSize(0.034);
          hD->SetLineColor(m + 2);
          hD->SetLineWidth(2);
          hD->SetLineStyle(1);

          TH1D* hM = nullptr;
          double mcEtaPi0Ratio = 0.0;
          if (drawDVPi0MC) {
            auto rdf_mc = plotters[m]->GetRDF_DVCSPi0MC().Filter(cutExpr, cutLabel);
            if (rdf_mc.HasColumn("Mass_pi0")) {
              auto h_mc = rdf_mc.Histo1D(
                  {Form("h_mc_Mass_pi0_wide_%s_%zu", cleanName.c_str(), m),
                   "Wide Pi0 Mass;M_{#gamma#gamma} [GeV];Counts", 200, 0.0, 1.0},
                  "Mass_pi0");
              h_mc.GetValue();
              hM = (TH1D*)h_mc.GetPtr()->Clone();
              hM->SetDirectory(0);
              const double mcPi0Count = countInRange(hM, 0.05, 0.265);
              const double mcEtaCount = countInRange(hM, 0.4, 0.7);
              mcEtaPi0Ratio = mcPi0Count > 0.0 ? mcEtaCount / mcPi0Count : 0.0;
              NormalizeHistogram(hM);
              styleKin_.StyleTH1(hM);
              hM->SetLineColor(m + 2);
              hM->SetLineWidth(2);
              hM->SetLineStyle(2);
            }
          }

          if (firstMpi0) {
            hD->Draw("HIST");
            frameHistogram = hD;
            firstMpi0 = false;
          } else {
            hD->Draw("HIST SAME");
          }
          if (hM) hM->Draw("HIST SAME");

          maxY = std::max(maxY, hD->GetMaximum());
          if (hM) maxY = std::max(maxY, hM->GetMaximum());
          mpi0Legend->AddEntry(hD, (labels[m] + " Data").c_str(), "l");
          std::ostringstream dataRatioText;
          dataRatioText << "N_{#eta}/N_{#pi^{0}} = " << std::fixed << std::setprecision(3) << dataEtaPi0Ratio;
          mpi0Legend->AddEntry((TObject*)nullptr, dataRatioText.str().c_str(), "");
          if (hM) mpi0Legend->AddEntry(hM, (labels[m] + " MC").c_str(), "l");
          if (hM) {
            std::ostringstream mcRatioText;
            mcRatioText << "N_{#eta}/N_{#pi^{0}} = " << std::fixed << std::setprecision(3) << mcEtaPi0Ratio;
            mpi0Legend->AddEntry((TObject*)nullptr, mcRatioText.str().c_str(), "");
          }
        }

        if (frameHistogram) {
          const double displayMax = maxY * 1.2;
          frameHistogram->SetMaximum(displayMax);
          const double lineHeight = displayMax / 3.0;
          for (double boundary : {0.05, 0.265, 0.4, 0.7}) {
            TLine* rangeLine = new TLine(boundary, 0.0, boundary, lineHeight);
            rangeLine->SetLineColor(kGray + 2);
            rangeLine->SetLineStyle(2);
            rangeLine->SetLineWidth(2);
            rangeLine->Draw("SAME");
          }
          mpi0Legend->Draw();
          gPad->Modified();
          gPad->Update();
        }

        std::string mpi0Outpath = outputDir + "/Pi0Exclusivity_Mass_pi0_wide_" + cleanName + ".pdf";
        mpi0Canvas->SaveAs(mpi0Outpath.c_str());
        std::cout << "Saved wide-range Mpi0 comparison to: " << mpi0Outpath << "\n";
        delete mpi0Canvas;
      }
    }
  }

  void PlotPi0ExclusivityComparisonByDetectorCases_old(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {{"Mass_pi0", "Pi0 Mass", "M_{#pi0} [GeV]", -0.6, 0.6},
                                                                                           {"Emiss_pi0", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
                                                                                           {"PTmiss_pi0", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.4},
                                                                                           {"Theta_pi0pi0", "#theta(#pi0, #pi')", "#theta_{#pi#pi'} [deg]", -2.0, 4.0},
                                                                                           {"DeltaPhi_pi0", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
                                                                                           {"Mx2_eppi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep#pi) [GeV^{2}]", -0.05, 0.05},
                                                                                           {"Mx2_ep_pi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep) [GeV^{2}]", -0.5, 1},
                                                                                           {"Mx2_epi0", "Missing Mass Squared (e#pi)", "MM^{2}(e#pi) [GeV^{2}]", 0.0, 2.0},
                                                                                           {"Theta_epho1", "Angle: e-#gamma_{1}", "#theta(e, #gamma_{1}) [deg]", 0.0, 60.0},
                                                                                           {"Theta_epho2", "Angle: e-#gamma_{2}", "#theta(e, #gamma_{2}) [deg]", 0.0, 60.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 1800, 1200);
      int cols = 4;
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
          auto rdf_cutdata = plotters[m]->GetRDF_Pi0Data().Filter(cutExpr, cutLabel);
          auto rdf_cut = plotters[m]->GetRDF_DVCSPi0MC().Filter(cutExpr, cutLabel);
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
          stats << "#mu = " << std::fixed << std::setprecision(3) << mean << ", #sigma = " << std::fixed << std::setprecision(3) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Pi0Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };

  // phi analysis
  /// For exclusivity cuts, you can use the following function to select one triplet
  /*void PlotPhiAnaExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", 0.8, 1.3},
        {"Mx2_epKpKm", "Missing Mass Squared (epK^{+}K^{-})", "MM^{2}(epK^{+}K^{-}) [GeV^{2}]", -0.07, 0.07},
        {"Mx2_eKpKm", "Invariant Mass (eK^{+}K^{-})", "M^{2}(eK^{+}K^{-}) [GeV^{2}]", -0.5, 3},
        {"Mx2_epKp", "Missing Mass Squared (epK^{+})", "MM^{2}(epK^{+}) [GeV^{2}]", -0.5, 1.5},
        {"Mx2_epKm", "Missing Mass Squared (epK^{-})", "MM^{2}(epK^{-}) [GeV^{2}]", -0.5, 1.5},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.5},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
        {"Theta_e_phimeson", "Angle: e-#phi", "#theta(e, #phi) [deg]", 0.0, 60.0}};

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

        TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.88);
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
          h_clone->SetLineColorAlpha(m + 4,0.8);
          auto [cr, cg, cb] = modelShades[m % modelShades.size()];
          const int colorIdx = 4000 + int(m) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h_clone->SetMarkerColor(colorIdx);
          h_clone->SetLineColorAlpha(colorIdx, 0.8);
          h_clone->SetLineWidth(1);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColorAlpha(colorIdx, 0.8);
          line2->SetLineColorAlpha(colorIdx, 0.8);
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
          stats << "#mu = " << std::fixed << std::setprecision(3) << mean << ", #sigma = " << std::fixed << std::setprecision(3) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_Phi_Ana" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };*/

  void PlotPhiAnaExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", 0.8, 1.3},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.5},
        {"Mx2_epKpKm", "Missing Mass Squared (epK^{+}K^{-})", "MM^{2}(epK^{+}K^{-}) [GeV^{2}]", -0.07, 0.07},
        {"Mx2_eKpKm", "Invariant Mass Squared (eK^{+}K^{-})", "M^{2}(eK^{+}K^{-}) [GeV^{2}]", -0.5, 3.0},
        {"Mx2_epKm", "Missing Mass Squared (epK^{-})", "MM^{2}(epK^{-}) [GeV^{2}]", -0.5, 1.5},
        {"Mx2_epKp", "Missing Mass Squared (epK^{+})", "MM^{2}(epK^{+}) [GeV^{2}]", -0.5, 1.5},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0.0, 20.0},
        {"Theta_g_phimeson", "Angle: #gamma - #phi", "#theta(#gamma, #phi) [deg]", 0.0, 10.0},
        {"Theta_e_phimeson", "Angle: e - #phi", "#theta(e, #phi) [deg]", 0.0, 60.0},
        {"DeltaE", "Energy Difference", "#DeltaE [GeV]", -1.0, 1.0},
        {"Cone_p", "Cone Angle (p)", "Cone(p) [deg]", 0.0, 20.0},
        {"Cone_Kp", "Cone Angle (K^{+})", "Cone(K^{+}) [deg]", 0.0, 20.0},
        {"Cone_Km", "Cone Angle (K^{-})", "Cone(K^{-}) [deg]", 0.0, 20.0},
        {"Coplanarity_had_normals_deg", "Coplanarity of Hadronic Normals", "Coplanarity_{had} [deg]", 0.0, 20.0},

        // New: Missing mass of eK+
        {"Mx_eKp", "Missing Mass (eK^{+})", "MM(eK^{+}) [GeV]", 0.0, 4.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      // 4 columns now
      int cols = 4;
      int rows = (vars.size() + cols - 1) / cols;
      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 2000, 1200);
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.88);
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
          h_clone->SetLineColorAlpha(m + 4, 0.8);
          auto [cr, cg, cb] = modelShades[m % modelShades.size()];
          const int colorIdx = 4000 + int(m) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h_clone->SetMarkerColor(colorIdx);
          h_clone->SetLineColorAlpha(colorIdx, 0.8);
          h_clone->SetLineWidth(1);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColorAlpha(colorIdx, 0.8);
          line2->SetLineColorAlpha(colorIdx, 0.8);
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
          stats << "#mu = " << std::fixed << std::setprecision(3) << mean << ", #sigma = " << std::fixed << std::setprecision(3) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_Phi_Ana" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  }

  void PlotDIS_BSA_Cross_Section_AndCorr_Comparison(double pol = 1.0, bool plotBSA = true, bool plotDVCSCross = false, bool plotPi0Corr = false, bool plotAccCorr = false,
                                                    bool plotEffCorr = false, bool plotRadCorr = false, bool plotP1Cut = false, bool meanKinVar = false) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allBSA;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allBSAraw;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allDVCSCross;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allDVCSPolCross_postive;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allDVCSPolCross_negative;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0Corr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0DVCSdiffmc;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0DVCSdiffexp;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0BSA;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allAccCorr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allEffCorr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allRadCorr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allP1Cut;
    // job for chatgpt
    std::vector<std::vector<std::vector<std::vector<std::tuple<double, double, double>>>>> allBSAmeans;
    bool hasPi0CorrectedBSA = false;

    for (auto& p : plotters) {
      if (plotBSA) {
        auto [hraw, hcorr] = p->ComputeBSAWithRaw(fXbins, pol);
        if (p->getDoPi0Corr() && !hcorr.empty()) {
          hasPi0CorrectedBSA = true;
          allBSA.push_back(std::move(hcorr));
        } else {
          allBSA.push_back(hraw);
        }
        allBSAraw.push_back(std::move(hraw));
      }
      if (plotDVCSCross) {
        auto hists = p->ComputeDVCS_CrossSection(fXbins);
        allDVCSCross.push_back(std::move(hists));
        auto hists_pol_postive = p->ComputePolDVCS_CrossSection(1,fXbins);
        auto hists_pol_negative = p->ComputePolDVCS_CrossSection(-1,fXbins);
        allDVCSPolCross_postive.push_back(std::move(hists_pol_postive));
        allDVCSPolCross_negative.push_back(std::move(hists_pol_negative));
      }
      if (plotPi0Corr) {
        auto hcorr = p->ComputePi0Corr(fXbins);
        auto hpi0dvcsdiffmc = p->ComputePi0DVCSdiffmc(fXbins);
        auto hpi0dvcsdiffexp = p->ComputePi0DVCSdiffexp(fXbins);
        auto hpi0bsa = p->ComputePi0BSA(fXbins, pol);
        allPi0Corr.push_back(std::move(hcorr));
        allPi0DVCSdiffmc.push_back(std::move(hpi0dvcsdiffmc));
        allPi0DVCSdiffexp.push_back(std::move(hpi0dvcsdiffexp));
        allPi0BSA.push_back(std::move(hpi0bsa));
      }
      if (plotAccCorr) {
        auto hacc = p->ComputeAccCorr(fXbins);
        allAccCorr.push_back(std::move(hacc));
      }
      if (plotEffCorr) {
        auto heff = p->ComputeEffCorr(fXbins);
        allEffCorr.push_back(std::move(heff));
      }
      if (plotRadCorr) {
        auto hrad = p->ComputeRadCorr(fXbins);
        allRadCorr.push_back(std::move(hrad));
      }
      if (plotP1Cut) {
        auto hp1 = p->ComputeP1CutEffect(fXbins);
        allP1Cut.push_back(std::move(hp1));
      }
      if (meanKinVar) {
        allBSAmeans.push_back(getMeanQ2xBt(fXbins, p));
      }
    }

    const std::string modePrefix = fXbins.UsesSDHEP() ? "SDHEP_" : "";
    if (plotBSA) MakeTiledGridComparison(modePrefix + "DIS_BSA", "A_{LU}", allBSA, &allBSAmeans, -0.65, 0.65, "pdf", true, true, false, false, meanKinVar);
    if (plotBSA && hasPi0CorrectedBSA) MakeTiledGridComparison(modePrefix + "DIS_BSAraw", "A_{LU}^{raw}", allBSAraw, &allBSAmeans, -0.65, 0.65, "pdf", true, true, false, false, meanKinVar);
    if (plotDVCSCross){
        MakeTiledGridComparison(modePrefix + "DIS_Cross_Section", "d#sigma/d#phi [nb/GeV^4]", allDVCSCross, &allBSAmeans, 0.0001, 1, "pdf", false, false, true, true, meanKinVar);
        MakeTiledGridComparison(modePrefix + "DIS_PolCross_Section_Positive", "d#sigma^{+}/d#phi [nb/GeV^4]", allDVCSPolCross_postive, &allBSAmeans, 0.0001, 1, "pdf", false, false, true, true, meanKinVar);
        MakeTiledGridComparison(modePrefix + "DIS_PolCross_Section_Negative", "d#sigma^{-}/d#phi [nb/GeV^4]", allDVCSPolCross_negative, &allBSAmeans, 0.0001, 1, "pdf", false, false, true, true, meanKinVar);
      }
    if (plotPi0Corr) MakeTiledGridComparison(modePrefix + "DIS_pi0Corr", "#eta^{#pi^{0}}", allPi0Corr, &allBSAmeans, 0.0, 1, "pdf", false, true, false, false, meanKinVar);
    if (plotPi0Corr) MakeTiledGridComparison(modePrefix + "DIS_pi0DVCSdiffmc", "d_{mc}", allPi0DVCSdiffmc, &allBSAmeans, 0.0, 2, "pdf", false, true, false, false, meanKinVar);
    if (plotPi0Corr) MakeTiledGridComparison(modePrefix + "DIS_pi0DVCSdiffexp", "d_{exp}", allPi0DVCSdiffexp, &allBSAmeans, 0.0, 2, "pdf", false, true, false, false, meanKinVar);
    if (plotPi0Corr) MakeTiledGridComparison(modePrefix + "DIS_pi0BSA", "A_{LU}^{#pi^{0}}", allPi0BSA, &allBSAmeans, -0.65, 0.65, "pdf", true, true, false, false, meanKinVar);
    if (plotAccCorr) MakeTiledGridComparison("DIS_accCorr", "A_{acc}", allAccCorr, &allBSAmeans, 0.01, 1.0, "pdf", false, true, true, false, meanKinVar);
    if (plotEffCorr) MakeTiledGridComparison("DIS_effCorr", "A_{eff}", allEffCorr, &allBSAmeans, 0.1, 1.1, "pdf", false, true, false, false, meanKinVar);
    if (plotRadCorr) MakeTiledGridComparison("DIS_radCorr", "C_{rad}", allRadCorr, &allBSAmeans, 0.5, 1.5, "pdf", false, true, false, false, meanKinVar);
    if (plotP1Cut) MakeTiledGridComparison("DIS_P1Cut", "C_{P1}", allP1Cut, &allBSAmeans, 0.0, 1.2, "pdf", false, true, false, false, meanKinVar);
  }

  bool file_exists(const char* name) {
    struct stat buffer;
    return (stat(name, &buffer) == 0);
  }

  void dumpHistogram(TH1D* h,
                     double xB, double Q2, double t,
                     double xBmin, double xBmax,
                     double Q2min, double Q2max,
                     double tmin, double tmax,
                     const std::string& label,
                     const char* filename = "h_data.txt") {

    bool exists = file_exists(filename);

    std::ofstream fout(filename, std::ios::out | std::ios::app);
    if (!fout.is_open()) {
      std::cerr << "cannot open " << filename << " to write!\n";
      return;
    }

    if (!exists) {
      if (fXbins.UsesSDHEP())
        fout << "# xi_SDHEP\tcostheta_SDHEP\t-t\tphi_SDHEP\tvalue\terror\t"
             << "xiMin\txiMax\tcosthetaMin\tcosthetaMax\ttmin\ttmax\tlabel\n";
      else
        fout << "# xB\tQ2\t-t\tphi\tvalue\terror\t"
             << "xBmin\txBmax\tQ2min\tQ2max\ttmin\ttmax\tlabel\n";
    }

    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      double phi   = h->GetBinCenter(ibin);
      double value = h->GetBinContent(ibin);
      double err   = h->GetBinError(ibin);
      if (!std::isfinite(value) || !std::isfinite(err)) continue;

      fout << xB << "\t"
           << Q2 << "\t"
           << t << "\t"
           << phi << "\t"
           << value << "\t"
           << err << "\t"
           << xBmin << "\t"
           << xBmax << "\t"
           << Q2min << "\t"
           << Q2max << "\t"
           << tmin << "\t"
           << tmax << "\t"
           << label << "\n";
    }

    fout.close();
    std::cout << "Data " << h->GetName()
              << " written into " << filename
              << (exists ? " (appended)" : "") << "\n";
  }

  void MakeTiledGridComparison(const std::string& observableName, const std::string& yAxisTitle, const std::vector<std::vector<std::vector<std::vector<TH1D*>>>>& histograms,
                               const std::vector<std::vector<std::vector<std::vector<std::tuple<double, double, double>>>>>* meanValues, double yMin, double yMax,
                               const std::string& suffix = "png", bool fitSinusoid = false, bool setManualYrange = false, bool setLogY = false, bool fixlineYrang = true,
                               bool showMeanKin = false) {
    if (histograms.empty() || histograms[0].empty() || histograms[0][0].empty() || histograms[0][0][0].empty()) {
      std::cerr << "No histograms to compare.\n";
      return;
    }
    if (fXbins.HasCustomBins()) {
      const auto& customBins = fXbins.GetCustomBins();
      const std::string customOutputDir = outputDir + "/" + observableName;
      fs::create_directories(customOutputDir);

      for (size_t binIndex = 0; binIndex < customBins.size(); ++binIndex) {
        const auto& bin = customBins[binIndex];
        TCanvas* canvas = new TCanvas(
            Form("%s_bin%zu", observableName.c_str(), binIndex),
            Form("%s bin %zu", observableName.c_str(), binIndex),
            1200, 900);
        styleBSA_.StylePad((TPad*)canvas);
        if (setLogY) canvas->SetLogy();

        TLegend* legend = new TLegend(0.56, 0.72, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);

        bool first = true;
        std::vector<std::unique_ptr<TGraphErrors>> finiteGraphs;
        std::vector<std::unique_ptr<TF1>> fits;
        for (size_t model = 0; model < histograms.size(); ++model) {
          if (histograms[model].empty() ||
              binIndex >= histograms[model][0].size() ||
              histograms[model][0][binIndex].empty()) {
            continue;
          }
          TH1D* h = histograms[model][0][binIndex][0];
          if (!h) continue;

          styleBSA_.StyleTH1(h);
          auto [r, g, b] = modelShades[model % modelShades.size()];
          const int colorIndex = 3000 + model * 20;
          if (!gROOT->GetColor(colorIndex)) new TColor(colorIndex, r, g, b);
          h->SetLineColor(colorIndex);
          h->SetMarkerColor(colorIndex);
          h->SetLineWidth(1);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetStats(0);
          h->GetXaxis()->SetTitle(Form("%s [deg]", fXbins.GetPhiLabel()));
          h->GetYaxis()->SetTitle(yAxisTitle.c_str());
          h->GetXaxis()->SetRangeUser(0.0, 360.0);
          if (setManualYrange) h->GetYaxis()->SetRangeUser(yMin, yMax);

          auto graph = std::make_unique<TGraphErrors>();
          for (int histBin = 1; histBin <= h->GetNbinsX(); ++histBin) {
            const double x = h->GetBinCenter(histBin);
            const double y = h->GetBinContent(histBin);
            const double ey = h->GetBinError(histBin);
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(ey)) continue;
            const int point = graph->GetN();
            graph->SetPoint(point, x, y);
            graph->SetPointError(point, 0.0, ey);
          }
          if (graph->GetN() == 0) {
            std::cout << "No finite points for " << observableName
                      << " custom bin " << binIndex
                      << ", model " << labels[model] << "; skipping draw.\n";
            continue;
          }

          graph->SetName(Form("graph_%s_bin%zu_model%zu", observableName.c_str(), binIndex, model));
          graph->SetTitle("");
          graph->SetLineColor(colorIndex);
          graph->SetMarkerColor(colorIndex);
          graph->SetLineWidth(1);
          graph->SetMarkerStyle(20);
          graph->SetMarkerSize(1.0);
          graph->Draw(first ? "AP" : "P SAME");
          if (first) {
            graph->GetXaxis()->SetTitle(Form("%s [deg]", fXbins.GetPhiLabel()));
            graph->GetYaxis()->SetTitle(yAxisTitle.c_str());
            graph->GetXaxis()->SetLimits(0.0, 360.0);
            if (setManualYrange) {
              graph->SetMinimum(yMin);
              graph->SetMaximum(yMax);
            }
          }
          first = false;
          legend->AddEntry(graph.get(), labels[model].c_str(), "p");

          if (fitSinusoid && graph->GetN() >= 4) {
            auto fit = std::make_unique<TF1>(
                Form("fit_%s_bin%zu_model%zu", observableName.c_str(), binIndex, model),
                "[0] + ([1]*sin(x*TMath::DegToRad())) / (1 + [2]*cos(x*TMath::DegToRad()))",
                0.0, 360.0);
            fit->SetParameters(0.0, 0.2, 0.1);
            fit->SetLineColor(colorIndex);
            fit->SetLineStyle(2);
            fit->SetLineWidth(1);
            graph->Fit(fit.get(), "Q0");
            fit->Draw("SAME");
            fits.push_back(std::move(fit));
          }
          finiteGraphs.push_back(std::move(graph));

          if (showMeanKin && meanValues &&
              model < meanValues->size() &&
              !(*meanValues)[model].empty() &&
              binIndex < (*meanValues)[model][0].size() &&
              !(*meanValues)[model][0][binIndex].empty()) {
            auto [mean_xB, mean_Q2, mean_t] = (*meanValues)[model][0][binIndex][0];
            dumpHistogram(
                h, mean_xB, mean_Q2, mean_t,
                bin.xBMin, bin.xBMax,
                bin.Q2Min, bin.Q2Max,
                bin.tMin, bin.tMax,
                labels[model],
                Form("datapoint_%s.txt", observableName.c_str()));
          }
        }

        TLatex binLabel;
        binLabel.SetNDC();
        binLabel.SetTextFont(42);
        binLabel.SetTextSize(0.032);
        binLabel.DrawLatex(
            0.14, 0.92,
            Form("%s:[%.3f,%.3f), %s:[%.3f,%.3f), t:[%.3f,%.3f)",
                 fXbins.GetFirstLabel(), bin.xBMin, bin.xBMax,
                 fXbins.GetSecondLabel(), bin.Q2Min, bin.Q2Max,
                 bin.tMin, bin.tMax));
        legend->Draw();
        canvas->SaveAs(
            Form("%s/%s_bin%03zu.%s",
                 customOutputDir.c_str(), observableName.c_str(), binIndex, suffix.c_str()));
        delete canvas;
      }
      return;
    }

    const auto& q2_edges = fXbins.GetQ2Bins();
    const auto& t_edges = fXbins.GetTBins();

    const size_t n_q2 = q2_edges.size() - 1;
    const size_t n_t = t_edges.size() - 1;
    const size_t n_xb = fXbins.GetMaxXBBinCount();

    const int rows = n_q2;
    const int cols = n_xb;

    bool Doplot = true;
    int first_perbin_xb = 0;
    bool first_perbin_q2 = true;
    bool first_first_perbin_q2 = true;
    double fixlineYMin = 0.00001;
    double fixlineYMax = 1.0;

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      TString cname = Form("DIS_BSA_t[%zu]", t_bin);
      TCanvas* c = new TCanvas(cname, cname, 2200, 1600);
      first_perbin_xb = 0;
      first_first_perbin_q2 = true;
      double canvasBorderX = 0.06;
      double canvasBorderY = 0.08;
      double gpad_margin_ratio = 0.2;

      double cellW = (1 - 2 * canvasBorderX) / cols, cellH = (1 - 2 * canvasBorderY) / rows;

      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        const auto& xb_edges = fXbins.GetXBBins(q2_bin);
        first_perbin_q2 = true;

        /// xbin loop
        for (size_t xb_bin = 0; xb_bin + 1 < xb_edges.size(); ++xb_bin) {
          int visualRow = rows - 1 - q2_bin;
          int pad = visualRow * cols + xb_bin + 1;
          c->cd();

          bool first = true;
          gStyle->SetCanvasPreferGL(true);

          TLegend* leg = new TLegend(0.35, 0.85, 0.85, 0.95);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.08);

          TLegend* legParams = new TLegend(0.35, 0.16, 0.85, 0.32);  // Bottom legend for a₁
          legParams->SetBorderSize(0);
          legParams->SetFillStyle(0);
          legParams->SetTextSize(0.08);

          TPad* thisPad = new TPad(Form("%zu_%zu", xb_bin, q2_bin), Form("%zu_%zu", xb_bin, q2_bin), cellW * xb_bin + canvasBorderX, cellH * (q2_bin) + canvasBorderY,
                                   cellW * (xb_bin + 1) + canvasBorderX, cellH * (q2_bin + 1) + canvasBorderY);
          double l = 0.00, r = 0.00, b = 0.00, t = 0.00;
          Doplot = false;

          for (size_t m = 0; m < histograms.size(); ++m) {
            // if(q2_bin == 2 && xb_bin == 2) continue; // save first per bin xb
            //  Pad margins
            /*
            double l = (first_perbin_q2) ? 0.2 : 0.00;
            double r = (!first_perbin_q2) ? 0.00 : 0.00;
            double b = (xb_bin == first_perbin_xb) ? 0.16 : 0.00;
            double t = (visualRow == 0) ? 0.000 : 0.00;
            */

            // const int idx = q2_bin * (n_t * n_xb) + t_bin * n_xb + xb_bin;

            TH1D* h = histograms[m][xb_bin][q2_bin][t_bin];
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
            h->SetMarkerSize(1.0);
            h->SetStats(0);

            if (first) {
              l = (first_perbin_q2) ? (gpad_margin_ratio) / (1 + gpad_margin_ratio) : 0.00;
              r = (xb_bin == first_perbin_xb) ? (gpad_margin_ratio) / (1 + gpad_margin_ratio) : 0.00;
              b = (xb_bin == first_perbin_xb) ? (gpad_margin_ratio) / (1 + gpad_margin_ratio) : 0.00;
              t = (visualRow == 0) ? 0.000 : 0.00;

              l = (first_perbin_q2 && xb_bin == first_perbin_xb) ? (gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio) : l;
              r = (xb_bin == first_perbin_xb && first_perbin_q2) ? (gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio) : r;

              styleBSA_.StylePad(thisPad, l, r, b, t);
              thisPad->SetTicks(1, 0);
              thisPad->SetFillStyle(4000);
              // std::cout << "xb_bin: " << xb_bin << ", q2_bin: " << q2_bin << ", first_perbin_xb: " << first_perbin_xb << ", first_perbin_q2: " << first_perbin_q2 << ",
              // first_first_perbin_q2: " <<first_first_perbin_q2<< std::endl;

              h->GetXaxis()->SetTitle((xb_bin == first_perbin_xb) ? Form("%s [deg]", fXbins.GetPhiLabel()) : "");
              h->GetYaxis()->SetTitle((first_perbin_q2) ? yAxisTitle.c_str() : "");
              h->GetXaxis()->SetLabelSize((xb_bin == first_perbin_xb) ? 0.085 : 0.0);
              h->GetXaxis()->SetTitleSize((xb_bin == first_perbin_xb) ? 0.095 : 0.0);
              h->GetYaxis()->SetLabelSize((first_perbin_q2) ? 0.085 : 0.0);
              h->GetYaxis()->SetTitleSize((first_perbin_q2) ? 0.1 : 0.0);
              if (xb_bin == first_perbin_xb && first_perbin_q2) {
                h->GetYaxis()->SetLabelSize(0.085 * (1 + gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio));
                h->GetYaxis()->SetTitleSize(0.1 * (1 + gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio));
              }
            }

            h->GetXaxis()->SetTitleOffset((xb_bin == first_perbin_xb) ? 0.82 : 0.0);
            h->GetYaxis()->SetTitleOffset((first_perbin_q2) ? 0.82 : 0.0);

            h->GetXaxis()->SetNdivisions(4, false);
            h->GetYaxis()->SetNdivisions(6, true);
            if (setManualYrange) h->GetYaxis()->SetRangeUser(yMin, yMax);
            if (fixlineYrang && first_perbin_q2 && h->GetMinimum() > 0) {
              fixlineYMin = h->GetMinimum() * 0.3;
              fixlineYMax = h->GetMaximum() * 3;
              h->GetYaxis()->SetRangeUser(fixlineYMin, fixlineYMax);
            }
            if (fixlineYrang && !first_perbin_q2) {
              h->GetYaxis()->SetRangeUser(fixlineYMin, fixlineYMax);
            }
            h->GetXaxis()->SetRangeUser(0, 360);

            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);
            Doplot = !(!h || h->GetBinContent(5) == 0) || Doplot;
            if (!h || h->GetBinContent(5) == 0) {
              continue;
            }

            if (first) {
              if (xb_bin == first_perbin_xb && first_perbin_q2) {
                thisPad->SetPad(cellW * (xb_bin - gpad_margin_ratio) + canvasBorderX, cellH * (q2_bin - gpad_margin_ratio) + canvasBorderY,
                                cellW * (xb_bin + 1 + gpad_margin_ratio) + canvasBorderX, cellH * (q2_bin + 1) + canvasBorderY);
              } else if (xb_bin == first_perbin_xb && !first_perbin_q2) {
                h->GetXaxis()->ChangeLabel(1, -1, 0, -1, -1, -1, "");  // blank it out
                thisPad->SetPad(cellW * (xb_bin) + canvasBorderX, cellH * (q2_bin - gpad_margin_ratio) + canvasBorderY, cellW * (xb_bin + 1 + gpad_margin_ratio) + canvasBorderX,
                                cellH * (q2_bin + 1) + canvasBorderY);
              } else if (first_perbin_q2 && xb_bin != first_perbin_xb) {
                thisPad->SetPad(cellW * (xb_bin - gpad_margin_ratio) + canvasBorderX, cellH * (q2_bin) + canvasBorderY, cellW * (xb_bin + 1) + canvasBorderX,
                                cellH * (q2_bin + 1) + canvasBorderY);
              } else if (!first_perbin_q2 && xb_bin != first_perbin_xb) {
                thisPad->SetPad(cellW * (xb_bin) + canvasBorderX, cellH * (q2_bin) + canvasBorderY, cellW * (xb_bin + 1) + canvasBorderX, cellH * (q2_bin + 1) + canvasBorderY);
              }
            }

            if (first) thisPad->Draw();
            thisPad->cd();
            if (setLogY) thisPad->SetLogy();
            h->Draw(first ? "E1X0" : "E1X0 SAME");
            first = false;
            first_perbin_q2 = false;

            // Fit function and extract a₁
            if (fitSinusoid) {
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
            }
            leg->AddEntry(h, labels[m].c_str(), "p");
            // auto [mean_xB, mean_Q2, mean_t] = meanValues[m][xb_bin][q2_bin][t_bin];
            if (showMeanKin) {
              auto [mean_xB, mean_Q2, mean_t] = (*meanValues)[m][xb_bin][q2_bin][t_bin];
              TString meanText = Form("<%s> = %.2f, <%s> = %.2f, <|t|> = %.2f",
                                      fXbins.GetFirstLabel(), mean_xB,
                                      fXbins.GetSecondLabel(), mean_Q2, mean_t);
              TLatex* meanLatex = new TLatex(0.25, 0.78 - m * 0.10, meanText.Data());
              meanLatex->SetTextSize(0.05);
              meanLatex->SetNDC();
              meanLatex->SetTextFont(42);
              meanLatex->Draw();
              dumpHistogram(h, mean_xB, mean_Q2, mean_t, xb_edges[xb_bin], xb_edges[xb_bin + 1], q2_edges[q2_bin], q2_edges[q2_bin + 1], t_edges[t_bin], t_edges[t_bin + 1], labels[m],
                            Form("datapoint_%s.txt", observableName.c_str()));
            }
          }
          if (!Doplot) {
            if (first_first_perbin_q2) first_perbin_xb++;
            std::cout << "No data for this bin combination, skipping...\n";
            continue;
          }
          /*
                    // Annotate bin ranges
                    double xB_low = xb_edges[xb_bin], xB_high = xb_edges[xb_bin + 1];
                    double Q2_low = q2_edges[q2_bin], Q2_high = q2_edges[q2_bin + 1];
                    TString labelText = Form("x_{B} #in [%.2f, %.2f], Q^{2} #in [%.1f, %.1f]", xB_low, xB_high, Q2_low, Q2_high);
                    TLatex* latex = new TLatex(0.25, 0.82, labelText.Data());
                    latex->SetTextSize(0.055);
                    latex->SetNDC();
                    latex->SetTextFont(42);
                    latex->Draw();
          */
          leg->Draw();
          if (fitSinusoid) legParams->Draw();

          thisPad->Modified();
          thisPad->Update();
          c->Modified();
          c->Update();

          if (first_first_perbin_q2 && xb_bin == first_perbin_xb) {
            first_perbin_xb++;
            first_first_perbin_q2 = false;
            // std::cout<< "1"<<std::endl;
          } else if (!first_first_perbin_q2 && xb_bin == first_perbin_xb) {
            first_perbin_xb++;
            // std::cout<< "2"<<std::endl;
          } else if (first_first_perbin_q2 && xb_bin != first_perbin_xb) {
            first_perbin_xb++;
            // std::cout<< "3"<<std::endl;
          } else if (!first_first_perbin_q2 && xb_bin != first_perbin_xb) {
            // std::cout<< "4"<<std::endl;
          }
        }
      }
      TString outfile = Form("%s/%s_t_%.2f-%.2f.%s", outputDir.c_str(), observableName.c_str(), t_edges[t_bin], t_edges[t_bin + 1], suffix.c_str());
      c->SaveAs(outfile);

      // std::cout << "Saved: " << outfile << '\n';

      delete c;
    }
  }
  // === Add to DISANAcomparer (public): =========================
  void PlotPhiInvMassPerBin_AllModels(const std::string& baseOutDir = "PhiInvMassFits", int nBins = 120, double mMin = 0.98, double mMax = 1.08, bool constrainSigma = true,
                                      double sigmaRef = 0.004, double sigmaFrac = 0.25, double branching = 1.0, bool doAcceptanceCorr = false, bool doEfficiencyCorr = false,
                                      bool doRadCorr = false) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiInvMassPerBin_AllModels] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      DISANAplotter* p = plotters[i].get();
      if (!p) continue;

      const std::string subdir = baseOutDir + "/" + labels[i];
      gSystem->Exec(Form("mkdir -p %s", subdir.c_str()));

      std::cout << "→ Fitting/drawing per-bin K^{+}K^{-} mass for model: " << labels[i] << " → " << subdir << "  (AccCorr=" << (doAcceptanceCorr ? "ON" : "OFF")
                << ", RadCorr=" << (doRadCorr ? "ON" : "OFF") << ")\n";

      if (doAcceptanceCorr) {
        p->MakePhiAcceptanceCorrection3D(fXbins, subdir + "/acc/", "mtprime", "mtprime", /*minGenEntries=*/10, constrainSigma, sigmaRef, sigmaFrac);
      }

      if (doRadCorr) {
        p->MakePhiRadiativeCorrection3D_FromRatio(fXbins, subdir + "/rad/", "rad_corr", "mtprime");
      }
      if (doEfficiencyCorr) {
        p->MakePhiEfficiencyCorrection3D(fXbins, subdir + "/eff/");
      }

      p->SetPlotApplyAcceptanceCorrection(doAcceptanceCorr);
      p->SetPlotApplyRadiativeCorrection(doRadCorr);

      (void)p->MakePhiMassFitCanvases3D(fXbins, subdir, nBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, branching);
    }
  }

  // =========================================================================
  // PlotPhiInvMassPerBin_XBBins_AllModels
  // =========================================================================
  // Drives MakePhiMassFitCanvases_XBBins() for every loaded model.
  // Bins dσ/dt' in (x_B, t') — the correct binning for gluon radius
  // extraction because at fixed x_B the value of W (and hence |t_min|)
  // is kinematically determined event-by-event instead of smeared across
  // a wide W range.
  // =========================================================================
  void PlotPhiInvMassPerBin_XBBins_AllModels(const std::string& baseOutDir = "PhiInvMassFits_xBBins", int nBins = 120, double mMin = 0.98, double mMax = 1.08,
                                             bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.25, double branching = 1.0, bool doAcceptanceCorr = false,
                                             bool doEfficiencyCorr = false, bool doRadCorr = false) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiInvMassPerBin_XBBins_AllModels] no models.\n";
      return;
    }

    // Require x_B binning to be set
    const auto& xbEdges = fXbins.GetXBBins();
    if (xbEdges.size() < 2) {
      std::cerr << "[PlotPhiInvMassPerBin_XBBins_AllModels] "
                   "No x_B bins found in BinManager. "
                   "Call SetXBinsRanges() with a BinManager that has GetXBBins() set.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      DISANAplotter* p = plotters[i].get();
      if (!p) continue;

      const std::string subdir = baseOutDir + "/" + labels[i];
      gSystem->Exec(Form("mkdir -p %s", subdir.c_str()));

      std::cout << "→ [xB-binned] Fitting K^{+}K^{-} mass for model: " << labels[i] << " → " << subdir << "  (AccCorr=" << (doAcceptanceCorr ? "ON" : "OFF")
                << ", RadCorr=" << (doRadCorr ? "ON" : "OFF") << ")\n";

      p->SetPlotApplyAcceptanceCorrection(doAcceptanceCorr);
      p->SetPlotApplyRadiativeCorrection(doRadCorr);

      p->MakePhiMassFitCanvases_XBBins(fXbins, subdir, nBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, branching);
    }
  }

  // =========================================================================
  // PlotPhiDSigmaDt_XBBins_FromCache
  // =========================================================================
  // Plots and exports dσ/dt' vs t' in each x_B bin for all models.
  // The CSV written per (model, x_B bin) has columns designed so that the
  // Python script can compute a physically correct |t_min| from the
  // data-weighted mean W (stored as W_mean) rather than from bin-boundary
  // arithmetic. The key columns written are:
  //
  //   tprime_center, tprime_lo, tprime_hi
  //   xB_lo, xB_hi
  //   W_mean        ← data-weighted mean W  (use this for tmin)
  //   Q2_mean       ← data-weighted mean Q²
  //   xB_mean       ← data-weighted mean x_B
  //   tmin_mean     ← mean |t_min|(Q²_i, W_i) averaged event-by-event
  //   GammaV_mean
  //   CrossSection, CrossSection_Err
  //   RawCounts, RawCounts_Err
  // =========================================================================
  void PlotPhiDSigmaDt_XBBins_FromCache(bool writeCSV = true, bool plotData = true) {
    if (plotters.empty()) return;

    const auto& xbEdges = fXbins.GetXBBins();
    const auto& tpEdges = fXbins.GetTprimeBins();
    if (xbEdges.size() < 2 || tpEdges.size() < 2) {
      std::cerr << "[PlotPhiDSigmaDt_XBBins_FromCache] "
                   "Missing x_B or t' bin edges in BinManager.\n";
      return;
    }

    const size_t nXB = xbEdges.size() - 1;
    const size_t nT = tpEdges.size() - 1;

    // ---- CSV writer --------------------------------------------------------
    auto writeCSVForXBBin = [&](size_t im, size_t ixB) {
      TH1D* hXS = nullptr;
      TH1D* hAcc = nullptr;
      TH1D* hEff = nullptr;
      TH1D* hRad = nullptr;
      TH1D* hNsig = nullptr;

      {
        auto& H = plotters[im]->GetPhiDSigmaDt_XBT();
        if (ixB < H.size()) hXS = H[ixB];
      }
      {
        auto& A = plotters[im]->GetPhiAcceptance_XBT();
        if (ixB < A.size()) hAcc = A[ixB];
      }
      {
        auto& E = plotters[im]->GetPhiEfficiency_XBT();
        if (ixB < E.size()) hEff = E[ixB];
      }
      {
        auto& R = plotters[im]->GetPhiRadCorr_XBT();
        if (ixB < R.size()) hRad = R[ixB];
      }
      {
        auto& N = plotters[im]->GetPhiRawCounts_XBT();
        if (ixB < N.size()) hNsig = N[ixB];
      }
      if (!hXS) return;

      // Sanitise label for directory name
      std::string safeLabel = labels[im];
      for (char& ch : safeLabel)
        if (ch == ' ' || ch == '/' || ch == '\\') ch = '_';

      const std::string csvDir = outputDir + "/CSVs_xBBins/" + safeLabel;
      gSystem->Exec(Form("mkdir -p \"%s\"", csvDir.c_str()));

      const std::string fname = Form("%s/dsdt_xB%zu.csv", csvDir.c_str(), ixB);

      std::ofstream csv(fname);
      if (!csv.is_open()) {
        std::cerr << "[writeCSV_xB] Cannot open " << fname << "\n";
        return;
      }

      // Header — all quantities needed for tmin and xB extraction downstream
      csv << "tprime_center,tprime_lo,tprime_hi"
          << ",xB_lo,xB_hi"
          << ",xB_mean"    // data-weighted mean x_B in this (xB, t') cell
          << ",W_mean"     // data-weighted mean W  — USE THIS for tmin/xB
          << ",Q2_mean"    // data-weighted mean Q²
          << ",tmin_mean"  // mean |t_min|(Q²_i,W_i) averaged event-by-event
          << ",GammaV_mean"
          << ",CrossSection,CrossSection_Err"
          << ",Acceptance,Efficiency,RadCorr"
          << ",RawCounts,RawCounts_Err\n";
      csv << std::setprecision(6) << std::scientific;

      const double xbLo = xbEdges[ixB];
      const double xbHi = xbEdges[ixB + 1];

      for (int ib = 1; ib <= hXS->GetNbinsX(); ++ib) {
        const size_t it = static_cast<size_t>(ib - 1);
        const double tp_cen = hXS->GetBinCenter(ib);
        const double tp_lo = hXS->GetBinLowEdge(ib);
        const double tp_hi = hXS->GetBinLowEdge(ib + 1);

        const double xs = hXS->GetBinContent(ib);
        const double xs_err = hXS->GetBinError(ib);
        const double acc = hAcc ? hAcc->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double eff = hEff ? hEff->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double rad = hRad ? hRad->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig = hNsig ? hNsig->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig_err = hNsig ? hNsig->GetBinError(ib) : std::numeric_limits<double>::quiet_NaN();

        // Per-(xB,t') data-weighted means from new getters
        const double xb_mean = plotters[im]->GetPhiMeanXB_XBT(ixB, it);
        const double w_mean = plotters[im]->GetPhiMeanW_XBT(ixB, it);
        const double q2_mean = plotters[im]->GetPhiMeanQ2_XBT(ixB, it);
        const double tmin_mean = plotters[im]->GetPhiMeanTmin_XBT(ixB, it);
        const double gv_mean = plotters[im]->GetPhiMeanGammaV_XBT(ixB, it);

        csv << tp_cen << "," << tp_lo << "," << tp_hi << "," << xbLo << "," << xbHi << "," << xb_mean << "," << w_mean << "," << q2_mean << "," << tmin_mean << "," << gv_mean
            << "," << xs << "," << xs_err << "," << acc << "," << eff << "," << rad << "," << nsig << "," << nsig_err << "\n";
      }
      csv.close();
      std::cout << "[CSV-xB] Written → " << fname << "\n";
    };

    // ---- Plotting ----------------------------------------------------------
    auto doOnePlot = [&](size_t ixB) {
      const double xbLo = xbEdges[ixB];
      const double xbHi = xbEdges[ixB + 1];
      const TString head = Form("x_{B} #in [%.3f, %.3f]", xbLo, xbHi);

      TCanvas* c = new TCanvas(Form("c_phi_dsdt_xB%zu", ixB), "", 1200, 900);
      styleCrossSection_.StylePad((TPad*)gPad);
      gPad->SetTicks(1, 1);
      gPad->SetLogy();

      TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.035);

      bool first = true;
      for (size_t im = 0; im < plotters.size(); ++im) {
        auto& H = plotters[im]->GetPhiDSigmaDt_XBT();
        if (ixB >= H.size() || !H[ixB]) continue;
        TH1D* h = H[ixB];

        styleCrossSection_.StyleTH1(h);
        auto [cr, cg, cb] = modelShades[im % modelShades.size()];
        const int colorIdx = 5500 + (int)im * 20;
        if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
        h->SetLineColor(colorIdx);
        h->SetMarkerColor(colorIdx);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->SetTitle("");
        h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
        h->GetYaxis()->SetTitle("d#sigma/dt' [nb/GeV^{2}]");

        if (first) {
          h->Draw("E1X0");
          TLatex latex;
          latex.SetNDC();
          latex.SetTextFont(42);
          latex.SetTextSize(0.040);
          latex.DrawLatex(0.14, 0.93, head);

          // Annotate with mean kinematics from the lowest-t' bin (bin 0)
          const double wmean = plotters[im]->GetPhiMeanW_XBT(ixB, 0);
          const double q2mean = plotters[im]->GetPhiMeanQ2_XBT(ixB, 0);
          if (std::isfinite(wmean) && std::isfinite(q2mean)) latex.DrawLatex(0.14, 0.88, Form("<W> = %.2f GeV,  <Q^{2}> = %.2f GeV^{2}", wmean, q2mean));
          first = false;
        } else {
          h->Draw("E1X0 SAME");
        }
        leg->AddEntry(h, labels[im].c_str(), "lep");
      }
      leg->Draw();

      const TString out = Form("%s/phi_dsdt_xBBin_%zu.pdf", outputDir.c_str(), ixB);
      c->SaveAs(out);
      delete leg;
      delete c;
    };

    // ---- Main loop ---------------------------------------------------------
    for (size_t ixB = 0; ixB < nXB; ++ixB) {
      if (plotData) doOnePlot(ixB);
      if (writeCSV) {
        for (size_t im = 0; im < plotters.size(); ++im) writeCSVForXBBin(im, ixB);
      }
    }
  }

  void PlotPhiDSigmaDt_FromCache(bool plotData = true, bool plotAcc = false, bool plotRadCorr = false, bool plotEff = false, bool writeCSV = true) {
    if (plotters.empty()) return;

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const auto& t_edges = fXbins.GetTBins();  // |t| bin edges (may be empty)
    const bool hasW = !w.empty();
    const bool hasT = (t_edges.size() > 1);

    const size_t nQ = (q2.size() > 1) ? (q2.size() - 1) : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    // -----------------------------------------------------------------------
    //  CSV helper – writes one CSV per (model, Q2-bin, W-bin)
    //  Columns: tprime_center, tprime_lo, tprime_hi,
    //           Q2_center, Q2_lo, Q2_hi,
    //           t_center, t_lo, t_hi,       (only when t-bins exist)
    //           W_center, W_lo, W_hi,        (only when W-bins exist)
    //           CrossSection, Mean, Err,
    //           Acceptance, RadCorr
    // -----------------------------------------------------------------------
    auto writeCSVForBin = [&](size_t im, size_t iq, size_t iw) {
      // --- gather the three histograms (xs mandatory, acc/rad optional) ---
      TH1D* hXS = nullptr;
      TH1D* hAcc = nullptr;
      TH1D* hEff = nullptr;
      TH1D* hRad = nullptr;
      TH1D* hNsig = nullptr;

      {
        auto& H = plotters[im]->GetPhiDSigmaDt3D();
        if (iq < H.size() && iw < H[iq].size()) hXS = H[iq][iw];
      }
      {
        auto& A = plotters[im]->GetPhiAcceptance3D();
        if (iq < A.size() && iw < A[iq].size()) hAcc = A[iq][iw];
      }
      {
        auto& E = plotters[im]->GetPhiEfficiency3D();
        if (iq < E.size() && iw < E[iq].size()) hEff = E[iq][iw];
      }
      {
        auto& R = plotters[im]->GetPhiRadCorr3D();
        if (iq < R.size() && iw < R[iq].size()) hRad = R[iq][iw];
      }
      {
        auto& N = plotters[im]->GetPhiRawCounts3D();
        if (iq < N.size() && iw < N[iq].size()) hNsig = N[iq][iw];
      }

      if (!hXS) return;  // nothing to write

      // --- build output path: CSVs/<ModelName>/ ---
      // Sanitise label for use as directory name (replace spaces / slashes)
      std::string safeLabel = labels[im];
      for (char& ch : safeLabel)
        if (ch == ' ' || ch == '/' || ch == '\\') ch = '_';

      const std::string csvDir = outputDir + "/CSVs/" + safeLabel;
      gSystem->Exec(Form("mkdir -p \"%s\"", csvDir.c_str()));

      // One CSV per (Q2-bin, W-bin)
      std::string fname;
      if (hasW)
        fname = Form("%s/dsdt_Q%zu_W%zu.csv", csvDir.c_str(), iq, iw);
      else
        fname = Form("%s/dsdt_Q%zu.csv", csvDir.c_str(), iq);

      std::ofstream csv(fname);
      if (!csv.is_open()) {
        std::cerr << "[writeCSV] Cannot open " << fname << "\n";
        return;
      }

      // --- convenience bin-centre/edge values for outer axes ---
      const double q2_lo = q2[iq];
      const double q2_hi = q2[iq + 1];
      const double q2_cen = 0.5 * (q2_lo + q2_hi);

      const double w_lo = hasW ? w[iw] : std::numeric_limits<double>::quiet_NaN();
      const double w_hi = hasW ? w[iw + 1] : std::numeric_limits<double>::quiet_NaN();
      // NOTE: w_cen is the arithmetic midpoint of the W bin boundaries (e.g. 5.9 GeV
      // for a 1.8–10 GeV bin).  It is only used for W2 columns; the per-row W_mean
      // below comes from the actual data-weighted mean and is the physically meaningful
      // quantity for tmin and xB computation.
      const double w_cen = hasW ? 0.5 * (w_lo + w_hi) : std::numeric_limits<double>::quiet_NaN();

      // --- write header ---
      csv << "tprime_center,tprime_lo,tprime_hi"
          << ",t_center,t_lo,t_hi"
          << ",Q2_center,Q2_lo,Q2_hi";
      if (hasW)
        csv << ",W2_center,W2_lo,W2_hi"
            << ",W_lo,W_hi"  // bin boundaries kept for reference
            << ",W_center";  // data-weighted mean W (correct for tmin/xB)

      // NEW mean-kin columns (per (Q2,W,t') bin)
      csv << ",xB_mean,W_mean,GammaV_mean";

      csv << ",CrossSection,CrossSection_Err,tprime_mean"
          << ",Acceptance,Efficiency,RadCorr"
          << ",RawCounts,RawCounts_Err\n";
      csv << std::setprecision(6) << std::scientific;

      const int nBins = hXS->GetNbinsX();
      for (int ib = 1; ib <= nBins; ++ib) {
        const double tp_cen = hXS->GetBinCenter(ib);
        const double tp_lo = hXS->GetBinLowEdge(ib);
        const double tp_hi = hXS->GetBinLowEdge(ib + 1);

        // |t| from t-edges array if available, else fall back to tprime bin
        double t_cen, t_lo_val, t_hi_val;
        if (hasT && (ib - 1) < (int)(t_edges.size() - 1)) {
          t_lo_val = t_edges[ib - 1];
          t_hi_val = t_edges[ib];
          t_cen = 0.5 * (t_lo_val + t_hi_val);
        } else {
          t_cen = tp_cen;
          t_lo_val = tp_lo;
          t_hi_val = tp_hi;
        }

        const double xs = hXS->GetBinContent(ib);
        const double err = hXS->GetBinError(ib);
        const double mean_val = tp_cen;

        const double acc = hAcc ? hAcc->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double eff = hEff ? hEff->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double rad = hRad ? hRad->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig = hNsig ? hNsig->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig_err = hNsig ? hNsig->GetBinError(ib) : std::numeric_limits<double>::quiet_NaN();
        const int ib0 = ib - 1;

        const double xBmean = plotters[im]->GetPhiMeanXB(iq, iw, ib0);
        const double Wmean = plotters[im]->GetPhiMeanW(iq, iw, ib0);
        const double Gvmean = plotters[im]->GetPhiMeanGammaV(iq, iw, ib0);

        // W_center: use data-weighted mean W when available (finite and physical),
        // fall back to bin-boundary midpoint only as last resort.
        const double W_center_out = (std::isfinite(Wmean) && Wmean > 0.5) ? Wmean : w_cen;

        csv << tp_cen << "," << tp_lo << "," << tp_hi << "," << t_cen << "," << t_lo_val << "," << t_hi_val << "," << q2_cen << "," << q2_lo << "," << q2_hi;

        if (hasW) csv << "," << w_cen * w_cen << "," << w_lo * w_lo << "," << w_hi * w_hi << "," << w_lo << "," << w_hi << "," << W_center_out;

        // NEW means
        csv << "," << xBmean << "," << Wmean << "," << Gvmean;

        csv << "," << xs << "," << err << "," << mean_val << "," << acc << "," << eff << "," << rad << "," << nsig << "," << nsig_err << "\n";
      }

      csv.close();
      std::cout << "[CSV] Written → " << fname << "\n";
    };

    // -----------------------------------------------------------------------
    //  Main plotting lambda (unchanged from original)
    // -----------------------------------------------------------------------
    auto doOne = [&](const TString& tag, const TString& yTitle, double yMin, double yMax, std::function<TH1D*(DISANAplotter*, size_t, size_t)> getHist) {
      for (size_t iq = 0; iq < nQ; ++iq)
        for (size_t iw = 0; iw < nW; ++iw) {
          TCanvas* c = new TCanvas(Form("c_%s_Q%zu_W%zu", tag.Data(), iq, iw), "", 1200, 900);
          styleCrossSection_.StylePad((TPad*)gPad);
          gPad->SetTicks(1, 1);

          TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.035);

          TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

          bool first = true;

          for (size_t im = 0; im < plotters.size(); ++im) {
            TH1D* h = getHist(plotters[im].get(), iq, iw);
            if (!h) continue;

            styleCrossSection_.StyleTH1(h);

            auto [cr, cg, cb] = modelShades[im % modelShades.size()];
            const int colorIdx = 5000 + int(im) * 20;
            if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.0);

            h->SetTitle("");
            h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
            h->GetYaxis()->SetTitle(yTitle);
            if (yTitle == "d#sigma/dt' [nb/GeV^{2}]") {
              gPad->SetLogy();
            }

            if (first) {
              h->Draw("E1X0");
              TLatex latex;
              latex.SetNDC();
              latex.SetTextFont(42);
              latex.SetTextSize(0.040);
              latex.DrawLatex(0.14, 0.93, head);
              first = false;
            } else {
              h->Draw("E1X0 SAME");
            }

            leg->AddEntry(h, labels[im].c_str(), "lep");
          }

          leg->Draw();

          TString out = hasW ? Form("%s/%s_Q%zu_W%zu.pdf", outputDir.c_str(), tag.Data(), iq, iw) : Form("%s/%s_Q%zu.pdf", outputDir.c_str(), tag.Data(), iq);

          c->SaveAs(out);

          delete leg;
          delete c;
        }
    };

    // ---- Run enabled plot types ----
    if (plotData) {
      doOne("phi_dsdt", "d#sigma/dt' [nb/GeV^{2}]", 0.0, 10.0, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& H = P->GetPhiDSigmaDt3D();
        if (iq >= H.size() || iw >= H[iq].size()) return nullptr;
        return H[iq][iw];
      });
    }

    if (plotAcc) {
      doOne("phi_accept", "Acceptance", 0.0, .05, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& A = P->GetPhiAcceptance3D();
        if (iq >= A.size() || iw >= A[iq].size()) return nullptr;
        return A[iq][iw];
      });
    }
    if (plotEff) {
      doOne("phi_efficiency", "Efficiency", 0.0, .05, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& E = P->GetPhiEfficiency3D();
        if (iq >= E.size() || iw >= E[iq].size()) return nullptr;
        return E[iq][iw];
      });
    }

    if (plotRadCorr) {
      doOne("phi_radCorr", "C_{rad}", 0.5, 1.5, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& R = P->GetPhiRadCorr3D();
        if (iq >= R.size() || iw >= R[iq].size()) return nullptr;
        return R[iq][iw];
      });
    }

    // ---- CSV output: one CSV per (model × Q2-bin × W-bin) ----
    // Requires #include <fstream> and <iomanip> (add at top if not already present)
    if (writeCSV) {
      for (size_t im = 0; im < plotters.size(); ++im)
        for (size_t iq = 0; iq < nQ; ++iq)
          for (size_t iw = 0; iw < nW; ++iw) writeCSVForBin(im, iq, iw);
    }
  }

  // Generic per-(Q2,W) grid plotter of TH1 vs t'
  void PlotPhiPerBin_FromCache(const std::string& tag,                                       // used in output filename + canvas name
                               const std::string& yTitle,                                    // axis title
                               double yMin, double yMax,                                     // y range
                               std::function<TH1D*(DISANAplotter*, size_t, size_t)> getHist  // how to fetch TH1
  ) {
    if (plotters.empty()) return;

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = (q2.size() > 1) ? (q2.size() - 1) : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_%s_Q%zu_W%zu", tag.c_str(), iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;

        for (size_t im = 0; im < plotters.size(); ++im) {
          DISANAplotter* P = plotters[im].get();
          TH1D* h = getHist(P, iq, iw);
          if (!h) continue;

          // Style (same logic as your cross-section plot)
          styleCrossSection_.StyleTH1(h);

          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 5000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);

          h->SetTitle("");
          h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
          h->GetYaxis()->SetTitle(yTitle.c_str());
          h->GetYaxis()->SetRangeUser(yMin, yMax);

          if (first) {
            h->Draw("E1X0");
            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
            first = false;
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
        }

        leg->Draw();

        TString out = hasW ? Form("%s/%s_Q%zu_W%zu.pdf", outputDir.c_str(), tag.c_str(), iq, iw) : Form("%s/%s_Q%zu.pdf", outputDir.c_str(), tag.c_str(), iq);

        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  // === NEW: A_LU(cos(theta_KK)) workflow (mirrors inv-mass fits + cached dσ/dt pattern) ===
  void PlotPhiALUCosThetaPerBin_AllModels(const std::string& baseOutDir = "PhiALUCosThetaFits", int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                          bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUCosThetaPerBin] no models.\n";
      return;
    }
    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing helicity-separated K^{+}K^{-} mass per cos(theta_KK) bin for model: " << labels[i] << " → " << subdir << std::endl;
      plotters[i]->MakePhiBSAMassFitCanvases3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamPol);
    }
  }
  void PlotPhiALUCosTheta_FromCache(const std::string& outDir = "PhiALUCosTheta") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_alu_costh_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        double yAbsMax = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiALUCosTheta3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMaximum()));
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMinimum()));
        }
        yAbsMax = std::max(0.20, 1.2 * yAbsMax);

        const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiALUCosTheta3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;

          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);
          h->SetTitle("");
          h->GetXaxis()->SetTitle("cos#theta_{K^{+}K^{-}}");
          h->GetYaxis()->SetTitle("A_{LU}");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          h->GetYaxis()->SetRangeUser(-yAbsMax * 2.0, yAbsMax * 2.0);

          if (first) {
            h->Draw("E1X0");
            TLine z(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
            z.SetLineStyle(2);
            z.SetLineWidth(2);
            z.Draw("SAME");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        TString out = hasW ? Form("%s/BSA_vs_CosKK_Q%zu_W%zu.pdf", outBase.c_str(), iq, iw) : Form("%s/BSA_vs_CosKK_Q%zu.pdf", outBase.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  void PlotPhiBSATrentoPhiPerBin_AllModels(const std::string& baseOutDir = "PhiBSATrentoPhiFits", int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                           bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiBSATrentoPhiPerBin] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing helicity-separated K^{+}K^{-} mass per Trento-phi bin for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiBSATrentoPhiMassFitCanvases3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamPol);
    }
  }

  void PlotPhiBSATrentoPhi_FromCache(const std::string& outDir = "PhiBSATrentoPhi") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_bsa_trentophi_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        double yAbsMax = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiBSATrentoPhi3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMaximum()));
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMinimum()));
        }
        yAbsMax = std::max(0.20, 1.2 * yAbsMax);

        const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiBSATrentoPhi3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;

          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);

          h->SetTitle("");
          h->GetXaxis()->SetTitle("#phi_{Trento} [deg]");
          h->GetYaxis()->SetTitle("A_{LU}");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          h->GetYaxis()->SetRangeUser(-0.750, 0.750);

          if (first) {
            h->Draw("E1X0");
            TLine z(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
            z.SetLineStyle(2);
            z.SetLineWidth(2);
            z.Draw("SAME");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        TString out = hasW ? Form("%s/BSA_vs_trentoPhi_Q%zu_W%zu.pdf", outBase.c_str(), iq, iw) : Form("%s/BSA_vs_trentoPhi_Q%zu.pdf", outBase.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  // === NEW: A_LU(z_phi) workflow, analogous to cos(theta_KK) ===
  void PlotPhiALUZPhiPerBin_AllModels(const std::string& baseOutDir = "PhiALUZPhiFits", int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120, bool constrainSigma = true,
                                      double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUZPhiPerBin_AllModels] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing helicity-separated K^{+}K^{-} mass per z_phi bin for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiALUZPhiMassFitCanvases3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamPol);
    }
  }

  void PlotPhiALUZPhi_FromCache(const std::string& outDir = "PhiALUZPhi") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_alu_zphi_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        double yAbsMax = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& alu3D = plotters[im]->GetPhiALUZPhi3D();
          if (iq >= alu3D.size() || iw >= alu3D[iq].size()) continue;
          TH1D* h = alu3D[iq][iw];
          if (!h) continue;
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMaximum()));
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMinimum()));
        }
        yAbsMax = std::max(0.20, 1.2 * yAbsMax);

        const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;

        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& alu3D = plotters[im]->GetPhiALUZPhi3D();
          if (iq >= alu3D.size() || iw >= alu3D[iq].size()) continue;
          TH1D* h = alu3D[iq][iw];
          if (!h) continue;

          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);

          h->SetTitle("");
          h->GetXaxis()->SetTitle("z_{#phi}");
          h->GetYaxis()->SetTitle("A_{LU}");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          h->GetYaxis()->SetRangeUser(-yAbsMax, yAbsMax);

          if (first) {
            h->Draw("E1X0");
            TLine z(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
            z.SetLineStyle(2);
            z.SetLineWidth(2);
            z.Draw("SAME");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        TString out = hasW ? Form("%s/ALU_vs_zphi_Q%zu_W%zu.pdf", outBase.c_str(), iq, iw) : Form("%s/ALU_vs_zphi_Q%zu.pdf", outBase.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  // ================================================================
  //  R = sigma_L / sigma_T  extraction from cos(theta_H) angular fit
  // ================================================================

  /// Step 1 — run mass fits per (Q2, W, t', cos θ) bin and fill the
  ///          r04_00 / R cache in each plotter.
  void PlotPhiRLTPerBin_AllModels(const std::string& baseOutDir = "PhiRLTFits", int nMassBins = 40, double mMin = 0.988, double mMax = 1.15, bool constrainSigma = true,
                                  double sigmaRef = 0.004, double sigmaFrac = 0.25, double beamEnergyGeV = 10.6) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiRLTPerBin_AllModels] no models.\n";
      return;
    }
    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Extracting R=sigma_L/sigma_T from cos(theta_KK) angular fit for model: " << labels[i] << " → " << subdir << std::endl;
      plotters[i]->MakePhiRLTFromCosTheta3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamEnergyGeV);
    }
  }

  /// Step 2 — plot R and r04_00 vs t' (one canvas per Q2 bin) and
  ///          write a summary CSV.
  void PlotPhiRLT_FromCache(const std::string& outDir = "PhiRLT") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();
    const size_t nQ = q2.size() > 1 ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    // ---- Helper lambda: draw one observable (R or r04_00) ----
    auto drawObs = [&](const std::string& tag, const std::string& yTitle, double yLo, double yHi, std::function<TH1D*(DISANAplotter*, size_t, size_t)> getter, bool drawZeroLine) {
      for (size_t iq = 0; iq < nQ; ++iq) {
        for (size_t iw = 0; iw < nW; ++iw) {
          TCanvas* c = new TCanvas(Form("c_%s_Q%zu_W%zu", tag.c_str(), iq, iw), "", 1200, 900);
          styleCrossSection_.StylePad((TPad*)gPad);
          gPad->SetTicks(1, 1);

          TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.035);

          const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

          bool first = true;
          for (size_t im = 0; im < plotters.size(); ++im) {
            TH1D* h = getter(plotters[im].get(), iq, iw);
            if (!h) continue;
            styleCrossSection_.StyleTH1(h);
            auto [cr, cg, cb] = modelShades[im % modelShades.size()];
            const int colorIdx = 6000 + (int)im * 20;
            if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.0);
            h->SetTitle("");
            h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
            h->GetYaxis()->SetTitle(yTitle.c_str());
            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);
            h->GetYaxis()->SetRangeUser(yLo, yHi);
            if (first) {
              h->Draw("E1X0");
              if (drawZeroLine) {
                TLine zl(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
                zl.SetLineStyle(2);
                zl.SetLineWidth(2);
                zl.Draw("SAME");
              }
              TLatex latex;
              latex.SetNDC();
              latex.SetTextFont(42);
              latex.SetTextSize(0.040);
              latex.DrawLatex(0.14, 0.93, head);
              first = false;
            } else {
              h->Draw("E1X0 SAME");
            }
            leg->AddEntry(h, labels[im].c_str(), "lep");
          }
          leg->Draw();
          TString out = hasW ? Form("%s/%s_Q%zu_W%zu.pdf", outBase.c_str(), tag.c_str(), iq, iw) : Form("%s/%s_Q%zu.pdf", outBase.c_str(), tag.c_str(), iq);
          c->SaveAs(out);
          delete leg;
          delete c;
        }
      }
    };

    // ---- Draw R vs t' ----
    drawObs(
        "phi_RLT", "R = #sigma_{L}/#sigma_{T}", 0.0, 5.0,
        [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
          auto& H = P->GetPhiRLT3D();
          if (iq >= H.size() || iw >= H[iq].size()) return nullptr;
          return H[iq][iw];
        },
        false);

    // ---- Draw r04_00 vs t' ----
    drawObs(
        "phi_r04_00", "r^{04}_{00}", 0.0, 1.0,
        [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
          auto& H = P->GetPhiR04_003D();
          if (iq >= H.size() || iw >= H[iq].size()) return nullptr;
          return H[iq][iw];
        },
        false);

    // ---- Write summary CSV ----
    const std::string csvPath = outBase + "/RLT_summary.csv";
    std::ofstream csv(csvPath);
    if (csv.is_open()) {
      csv << "model,Q2_lo,Q2_hi,W_lo,W_hi,tprime_center,r04_00,r04_00_err,R,R_err\n";
      csv << std::setprecision(6) << std::scientific;
      for (size_t im = 0; im < plotters.size(); ++im) {
        auto& Hrlt = plotters[im]->GetPhiRLT3D();
        auto& Hr04 = plotters[im]->GetPhiR04_003D();
        for (size_t iq = 0; iq < nQ; ++iq) {
          for (size_t iw = 0; iw < nW; ++iw) {
            TH1D* hR = (iq < Hrlt.size() && iw < Hrlt[iq].size()) ? Hrlt[iq][iw] : nullptr;
            TH1D* hr04 = (iq < Hr04.size() && iw < Hr04[iq].size()) ? Hr04[iq][iw] : nullptr;
            if (!hR || !hr04) continue;
            for (int b = 1; b <= hR->GetNbinsX(); ++b) {
              csv << labels[im] << "," << q2[iq] << "," << q2[iq + 1] << "," << (hasW ? w[iw] : 0.0) << "," << (hasW ? w[iw + 1] : 0.0) << "," << hR->GetBinCenter(b) << ","
                  << hr04->GetBinContent(b) << "," << hr04->GetBinError(b) << "," << hR->GetBinContent(b) << "," << hR->GetBinError(b) << "\n";
            }
          }
        }
      }
      csv.close();
      std::cout << "[CSV] R=sigma_L/sigma_T summary → " << csvPath << "\n";
    }
  }

  // === A_LU^{sin(phi_Trento)}(cos(theta_KK)) using the moment method ===
  // Uses the invariant-mass window [mMin, mMax] instead of mass fits.
  void PlotPhiALUCosThetaPerBin_AllModels_SinPhiMoment(const std::string& baseOutDir = "PhiALUCosTheta_SinPhiMoment", double mMin = 0.9874, double mMax = 1.120,
                                                       double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUCosThetaPerBin_AllModels_SinPhiMoment] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];

      std::cout << "→ Computing A_{LU}^{sin#phi}(cos#theta_{KK}) via sin(phi_Trento) moment "
                << "in m(K^{+}K^{-}) ∈ [" << mMin << ", " << mMax << "]"
                << " for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiALUCosThetaSinPhiMoment3D(fXbins, subdir, mMin, mMax, beamPol);
    }
  }

  // A_LU(cos(theta_KK)) using a fit A*sin(phi)/(1 + b cos(phi))
  // in each cos(theta_KK) bin. Uses invMass_KpKm window [mMin, mMax].
  void PlotPhiALUCosThetaPerBin_AllModels_SinOver1PlusbCosFit(const std::string& baseOutDir = "PhiALUCosTheta_SinOver1PlusbCosFit", double mMin = 1.01, double mMax = 1.03,
                                                              double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUCosThetaPerBin_AllModels_SinOver1PlusbCosFit] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];

      std::cout << "→ Computing A_{LU}(cos#theta_{KK}) via fit "
                << "A*sin(#phi)/(1 + b cos(#phi)) "
                << "in m(K^{+}K^{-}) ∈ [" << mMin << ", " << mMax << "]"
                << " for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiALUCosTheta_SinOver1PlusbCosFit3D(fXbins, subdir, mMin, mMax, beamPol);
    }
  }

 private:
  static std::string FormatMillions(ULong64_t count) {
    std::ostringstream text;
    text << std::fixed << std::setprecision(3)
         << static_cast<double>(count)/1.0e6 << " M";
    return text.str();
  }

  BinManager fXbins;
  bool plotIndividual = false;
  bool useFittedYields_ = true;
  bool applyCorrection = false;

  DrawStyle style_;              // Default style
  DrawStyle styleKin_;           // Kin plot style
  DrawStyle styleDVCS_;          // DVCS plot style
  DrawStyle styleCrossSection_;  // Cross-section plot style
  DrawStyle styleBSA_;           // BSA plot style

  THnSparseD* correctionHist = nullptr;
  DVCSWeightFunction dvcs_weight_function_;

  std::unique_ptr<ROOT::RDF::RNode> rdf;
  std::string outputDir = ".";
  std::map<long long, bool> physicalVolumeCSVWritten_;

  std::vector<std::unique_ptr<DISANAplotter>> plotters;
  std::vector<std::string> labels;

  // Kept separate from the exclusive-DVCS and Phi plotters because the
  // current He3 model has only one generator-level MC dataframe.
  std::vector<std::unique_ptr<ROOT::RDF::RNode>> he3Models;

  // Coherent e' gamma He3' generator-level models. These deliberately avoid
  // the proton-mass assumptions in the standard exclusive-DVCS plotter.
  std::vector<std::unique_ptr<ROOT::RDF::RNode>> he3DVCSModels;
  std::vector<ULong64_t> he3DVCSModelEventCounts;

  // Independent He3DIS GEN/REC topology models.
  std::vector<He3DISModel> he3DISModels;
  int he3DISXActiveBins=6;
  double he3DISXActiveMin=0.0;
  double he3DISXActiveMax=0.96;
  double he3DISModelXMax=0.85;
  // Adaptive symmetric y half-ranges from the count projection, reused by
  // the corresponding generator-model raw and raw/PN panels.
  std::map<std::string,double> he3DISProjectionYExtents;
  std::map<std::string,double> he3DISProjectionYExtentsOverPN;
  // Per-bin GEN/REC statistical errors from the finite-count projection.
  std::map<std::string,std::array<std::vector<double>,2>>
      he3DISProjectionErrors;
  std::map<std::string,std::array<std::vector<double>,2>>
      he3DISProjectionErrorsOverPN;
  // Per-bin nucleon A_parallel centers extracted from the finite GEN/REC
  // spin counts.  PlotA1ProjectionHe3DIS reuses these exact values rather
  // than rebuilding a second asymmetry pipeline.
  std::map<std::string,std::array<std::vector<double>,2>>
      he3DISProjectionValuesOverPN;

  std::vector<std::string> particleName = {"e", "p", "#gamma"};
  std::map<std::string, std::string> typeToParticle = {{"el", "electron"},     {"pro", "proton"},   {"pho", "#gamma_{1}"},
                                                       {"pho2", "#gamma_{2}"}, {"kMinus", "K^{-}"}, {"kPlus", "K^{+}"}};
  std::map<std::string, std::string> VarName = {{"p", "p (GeV/#it{c})"}, {"theta", "#theta (rad)"}, {"phi", "#phi(rad)"}, {"vz", "v_{z}(cm)"}};
};
#endif  // DISANA_COMPARER_H

#ifndef DISANA_PLOTTER_H
#define DISANA_PLOTTER_H

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include <string>
#include <vector>

#include "DISANAMath.h"

class DISANAplotter {
 public:
DISANAplotter(
        ROOT::RDF::RNode                 df_dvcs_data,
        double                           beamEnergy,
        std::optional<ROOT::RDF::RNode>  df_pi0_data   = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcs_pi0mc    = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_pi0_pi0mc     = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_gen_dvcsmc = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_accept_dvcsmc = std::nullopt
        )
        :rdf(std::move(df_dvcs_data)),
        rdf_pi0_data(std::move(df_pi0_data)),
        rdf_dvcs_pi0mc(std::move(df_dvcs_pi0mc)),
        rdf_pi0_pi0mc(std::move(df_pi0_pi0mc)),
        rdf_gen_dvcsmc(std::move(df_gen_dvcsmc)),
        rdf_accept_dvcsmc(std::move(df_accept_dvcsmc)),
        beam_energy(beamEnergy)
    {}

  ~DISANAplotter() {
    disHistos.clear();
    kinematicHistos.clear();
  }
   
ROOT::RDF::RNode GetRDF() { return rdf; }

void SetPlotApplyCorrection(bool apply) {dopi0corr = apply;}
void SetPlotApplyAcceptanceCorrection(bool apply) {doacceptcorr = apply;}
bool getDoPi0Corr() const { return dopi0corr; }



  void GenerateKinematicHistos(const std::string& type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi"};
    for (const auto& v : vars) {
      std::string base = "rec" + type + "_" + v;
      //std::cout << " base is " << base << std::endl;
      double histMin, histMax;
      auto it = kinematicAxisRanges.find(base);
      if (it != kinematicAxisRanges.end()) {
        histMin = it->second.first;
        histMax = it->second.second;
      } else {
        auto minVal = rdf.Min(base);
        auto maxVal = rdf.Max(base);
        double lo = *minVal;
        double hi = *maxVal;
        double margin = std::max(1e-3, 0.05 * (hi - lo));
        histMin = lo - margin;
        histMax = hi + margin;
        kinematicAxisRanges[base] = {histMin, histMax};  // store it
      }
    
       // auto h_all = rdf.Histo1D({(base + "_all").c_str(), "", 100, histMin, histMax}, base);
    auto h_acc = rdf.Histo1D({(base).c_str(), "", 100, histMin, histMax}, base);
    kinematicHistos.push_back(h_acc);
    }

  std::map<std::string, std::pair<double, double>> cachedRanges;

  for (const auto& var : disvars) {
    double histMin = axisRanges[var].first;
    double histMax = axisRanges[var].second;
    // auto h_all = rdf.Histo1D({(var + "_all").c_str(), (var + " All").c_str(), 100, histMin, histMax}, var);
    auto h_acc = rdf.Histo1D({var.c_str(), var.c_str(), 100, histMin, histMax}, var);
    disHistos.push_back(h_acc);

    //auto h_acceptance = dynamic_cast<TH1*>(h_acc->Clone((var + "_acceptance").c_str()));
    //h_acceptance->SetTitle(("Acceptance for " + var).c_str());
    //acceptHistos.push_back(std::shared_ptr<TH1>(h_acceptance));
  }
  }

  std::vector<TH1*> GetDISHistograms() {
    std::vector<TH1*> allDIShisto;
    for (auto& h : disHistos) allDIShisto.push_back(h.GetPtr());
    return allDIShisto;
  }

  std::vector<TH1*> GetAllHistograms() {
    std::vector<TH1*> all;
    for (auto& h : kinematicHistos) all.push_back(h.GetPtr());
    for (auto& h : disHistos) all.push_back(h.GetPtr());
    return all;
  }

  std::vector<std::vector<std::vector<TH1D*>>> ComputeDVCS_CrossSection(const BinManager& bins, double luminosity) {
    auto result = kinCalc.ComputeDVCS_CrossSection(rdf, bins, luminosity);
    if (dopi0corr) {
      auto sigma_pi0_3d = kinCalc.ComputeDVCS_CrossSection(*rdf_pi0_data, bins, luminosity);
      result = UsePi0Correction(result,sigma_pi0_3d,ComputePi0Corr(bins));
    }
    if (doacceptcorr) {
      auto acc3D = ComputeAccCorr(bins);
      result = UseAccCorrection(result, acc3D);
    }
    return result;
  }

  /// BSA computations // we need to have refined version of this codes
  std::vector<std::vector<std::vector<TH1D*>>> ComputeBSA(const BinManager& bins, double luminosity, double pol = 1.0) {
    // 1) Helicity selection
    auto rdf_pos = rdf.Filter("REC_Event_helicity ==  1");
    auto rdf_neg = rdf.Filter("REC_Event_helicity == -1");
    // 2) One-pass cross-section per helicity

    auto sigma_pos_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pos, bins, luminosity);
    auto sigma_neg_3d = kinCalc.ComputeDVCS_CrossSection(rdf_neg, bins, luminosity);

    auto result = kinCalc.ComputeBeamSpinAsymmetry(sigma_pos_3d, sigma_neg_3d, pol);

    if (dopi0corr) {
      auto rdf_pi0_data_pos = rdf_pi0_data->Filter("REC_Event_helicity ==  1");
      auto rdf_pi0_data_neg = rdf_pi0_data->Filter("REC_Event_helicity == -1");
    // 2) One-pass cross-section per helicity

      auto sigma_pi0_pos_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pi0_data_pos, bins, luminosity);
      auto sigma_pi0_neg_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pi0_data_neg, bins, luminosity);
      
      auto corr3D = ComputePi0Corr(bins);
      auto Api0 = kinCalc.ComputeBeamSpinAsymmetry(sigma_pi0_pos_3d, sigma_pi0_neg_3d, pol);
      result = UsePi0CorrectionForBSA(result, Api0, corr3D);
    }
    return result;

  }

   /// BSA computations // we need to have refined version of this codes
  std::vector<std::vector<std::vector<TH1D*>>> ComputePi0Corr(const BinManager& bins) {
    if (!rdf_dvcs_pi0mc || !rdf_pi0_pi0mc || !rdf_pi0_data) {
      std::cerr << "[ComputePi0Corr] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcPi0Corr(*rdf_dvcs_pi0mc, *rdf_pi0_pi0mc, rdf, *rdf_pi0_data, bins);

  }

  std::vector<std::vector<std::vector<TH1D*>>> ComputeAccCorr(const BinManager& bins) {
    if (!rdf_gen_dvcsmc || !rdf_accept_dvcsmc) {
      std::cerr << "[ComputeAccCorr] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcAcceptanceCorr(*rdf_gen_dvcsmc, *rdf_accept_dvcsmc, bins);
  }

  std::vector<std::vector<std::vector<TH1D*>>> UseAccCorrection(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D,
                                                                  const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
    if (xs3D.size() != corr3D.size()) {
        std::cerr << "[UseAccCorrection] Dimension mismatch (level-0)\n";
        return {};
    }
    std::vector<std::vector<std::vector<TH1D*>>> xsCorr3D(xs3D.size());
    for (size_t iq = 0; iq < xs3D.size(); ++iq) {
        if (xs3D[iq].size() != corr3D[iq].size()) {
            std::cerr << "[UseAccCorrection] Dimension mismatch (level-1)\n";
            return {};
        }
        xsCorr3D[iq].resize(xs3D[iq].size());
        for (size_t it = 0; it < xs3D[iq].size(); ++it) {
            if (xs3D[iq][it].size() != corr3D[iq][it].size()) {
                std::cerr << "[UseAccCorrection] Dimension mismatch (level-2)\n";
                return {};
            }
            xsCorr3D[iq][it].resize(xs3D[iq][it].size());
            for (size_t ix = 0; ix < xs3D[iq][it].size(); ++ix) {
                TH1D* hXS   = xs3D  [iq][it][ix];   // σ_uncorr
                TH1D* hCorr = corr3D[iq][it][ix];   //     c
              
                if (!hXS) { xsCorr3D[iq][it][ix] = nullptr; continue; }
                
                TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

                const int nb = hXS->GetNbinsX();
                for (int b = 1; b <= nb; ++b) {
                    const double xs_val = hXS->GetBinContent(b);
                    const double xs_err = hXS->GetBinError  (b);
                    const double c_val  = (hCorr||hCorr->GetBinContent(b)==0) ? hCorr->GetBinContent(b) : 1.0;
                    const double c_err  = (hCorr||hCorr->GetBinContent(b)==0) ? hCorr->GetBinError  (b) : 0.0;

                    double val_corr = xs_val/(c_val);
                    double err_corr = std::sqrt(
                        std::pow(xs_err / c_val, 2) +
                        std::pow(xs_val * c_err / (c_val * c_val), 2));

                    if (c_val == 0.0||c_err == 0.0) val_corr = xs_val;
                    if (c_val == 0.0||c_err == 0.0) err_corr = xs_err;
                    //std::cout <<"xs_val " << xs_val << " xs_err " << xs_err << " c_val " << c_val << " c_err " << c_err << " val_corr " << val_corr << " err_corr " << err_corr << std::endl;

                    hNew->SetBinContent(b, val_corr);
                    hNew->SetBinError  (b, err_corr);
                }

                xsCorr3D[iq][it][ix] = hNew;
            }
        }
    }
    return xsCorr3D;
  }

  std::vector<std::vector<std::vector<TH1D*>>> UsePi0Correction(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D,
                                                                const std::vector<std::vector<std::vector<TH1D*>>>& xsPi03D,
                                                                  const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
    if (xs3D.size() != corr3D.size()) {
        std::cerr << "[UsePi0Correction] Dimension mismatch (level-0)\n";
        return {};
    }
    std::vector<std::vector<std::vector<TH1D*>>> xsCorr3D(xs3D.size());
    for (size_t iq = 0; iq < xs3D.size(); ++iq) {
        if (xs3D[iq].size() != corr3D[iq].size()) {
            std::cerr << "[UsePi0Correction] Dimension mismatch (level-1)\n";
            return {};
        }
        xsCorr3D[iq].resize(xs3D[iq].size());
        for (size_t it = 0; it < xs3D[iq].size(); ++it) {
            if (xs3D[iq][it].size() != corr3D[iq][it].size()) {
                std::cerr << "[UsePi0Correction] Dimension mismatch (level-2)\n";
                return {};
            }
            xsCorr3D[iq][it].resize(xs3D[iq][it].size());
            for (size_t ix = 0; ix < xs3D[iq][it].size(); ++ix) {
                TH1D* hXS   = xs3D  [iq][it][ix];   // σ_uncorr
                TH1D* hPi0XS = xsPi03D[iq][it][ix]; // σ_pi0
                TH1D* hCorr = corr3D[iq][it][ix];   //     c
              
                if (!hXS) { xsCorr3D[iq][it][ix] = nullptr; continue; }
                
                TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

                const int nb = hXS->GetNbinsX();
                for (int b = 1; b <= nb; ++b) {
                    const double xs_val = hXS->GetBinContent(b);
                    const double xs_err = hXS->GetBinError  (b);
                    const double pi0_val = hPi0XS ? hPi0XS->GetBinContent(b) : 0.0;
                    const double pi0_err = hPi0XS ? hPi0XS->GetBinError  (b) : 0.0;
                    const double c_val  = hCorr ? hCorr->GetBinContent(b) : 0.0;
                    const double c_err  = hCorr ? hCorr->GetBinError  (b) : 0.0;

                    const double val_corr = xs_val*(1.0 - c_val);
                    const double err_corr = std::sqrt(
                        std::pow(xs_err * (1.0 - c_val), 2) +
                        std::pow(xs_val * c_err, 2));

                    hNew->SetBinContent(b, val_corr);
                    hNew->SetBinError  (b, err_corr);
                }

                xsCorr3D[iq][it][ix] = hNew;
            }
        }
    }
    return xsCorr3D;
  }


  std::vector<std::vector<std::vector<TH1D*>>> UsePi0CorrectionForBSA(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D,
    const std::vector<std::vector<std::vector<TH1D*>>>& xsPi03D,
      const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
    if (xs3D.size() != corr3D.size()) {
      std::cerr << "[UsePi0Correction] Dimension mismatch (level-0)\n";
      return {};
    }
    std::vector<std::vector<std::vector<TH1D*>>> xsCorr3D(xs3D.size());
    for (size_t iq = 0; iq < xs3D.size(); ++iq) {
      if (xs3D[iq].size() != corr3D[iq].size()) {
        std::cerr << "[UsePi0Correction] Dimension mismatch (level-1)\n";
        return {};
      }
      xsCorr3D[iq].resize(xs3D[iq].size());
      for (size_t it = 0; it < xs3D[iq].size(); ++it) {
        if (xs3D[iq][it].size() != corr3D[iq][it].size()) {
          std::cerr << "[UsePi0Correction] Dimension mismatch (level-2)\n";
          return {};
        }
        xsCorr3D[iq][it].resize(xs3D[iq][it].size());
        for (size_t ix = 0; ix < xs3D[iq][it].size(); ++ix) {
          TH1D* hXS   = xs3D  [iq][it][ix];   // σ_uncorr
          TH1D* hPi0XS = xsPi03D[iq][it][ix]; // σ_pi0
          TH1D* hCorr = corr3D[iq][it][ix];   //     c

          if (!hXS) { xsCorr3D[iq][it][ix] = nullptr; continue; }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError  (b);
            const double pi0_val = hPi0XS ? hPi0XS->GetBinContent(b) : 0.0;
            const double pi0_err = hPi0XS ? hPi0XS->GetBinError  (b) : 0.0;
            const double c_val  = hCorr ? hCorr->GetBinContent(b) : 0.0;
            const double c_err  = hCorr ? hCorr->GetBinError  (b) : 0.0;

            const double val_corr = xs_val/(1.0 - c_val)-
                                    (pi0_val*c_val/(1.0 - c_val));
            const double err_corr = std::sqrt(
            std::pow(xs_err * (1.0 - c_val), 2) +
            std::pow(xs_val * c_err, 2));

            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError  (b, err_corr);
          }
          xsCorr3D[iq][it][ix] = hNew;
        }
      }
    }
    return xsCorr3D;
  }

  /*
  void ApplyPi0BkgCorr(THnSparseD *correctionHist ){
      kinCalc.SetApplyCorrPi0BKG(true);
      if (correctionHist) {
        kinCalc.SetCorrHist(correctionHist);
      } else {
        std::cerr << "Error: Correction histogram is not set.\n";
      }
  }
  */

 private:

  std::vector<std::string> disvars = {"Q2", "xB", "t", "W", "phi"};
  std::map<std::string, std::pair<double, double>> axisRanges = {{"Q2", {0.0, 15.0}}, {"xB", {0.0, 1.0}}, {"W", {1.0, 10.0}}, {"t", {0.0, 10.0}}, {"phi", {-180.0, 180.0}}};
  std::map<std::string, std::pair<double, double>> kinematicAxisRanges = {
      {"recel_p", {-0.05, 13.0}},    {"recel_theta", {-0.01, 1.0}}, {"recel_phi", {-0.01, 6.1}}, {"recpho_p", {-0.01, 10.0}},
      {"recpho_theta", {-0.01, 1.0}}, {"recpho_phi", {-0.01, 6.1}},   {"recpro_p", {-0.01, 2.0}},  {"recpro_theta", {-0.01, 2.0}},
      {"recpro_phi", {-0.01, 6.1}}
      // Add more as needed
  };

  DISANAMath kinCalc;
  std::unique_ptr<TFile> predFile;
  std::string predFileName = ".";
  double beam_energy;
  bool dopi0corr = false;
  bool doacceptcorr = false;
  std::string ttreeName;
  std::vector<ROOT::RDF::RResultPtr<TH1>> kinematicHistos, disHistos;
  std::vector<std::shared_ptr<TH1>> acceptHistos;
  ROOT::RDF::RNode rdf;
  //ROOT::RDF::RNode rdf_dvcs_data;
  std::optional<ROOT::RDF::RNode> rdf_pi0_data;
  std::optional<ROOT::RDF::RNode> rdf_dvcs_pi0mc;
  std::optional<ROOT::RDF::RNode> rdf_pi0_pi0mc;
  std::optional<ROOT::RDF::RNode> rdf_gen_dvcsmc;
  std::optional<ROOT::RDF::RNode> rdf_accept_dvcsmc;
  
};

#endif  // DISANA_PLOTTER_H
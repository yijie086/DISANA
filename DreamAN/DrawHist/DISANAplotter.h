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

class BinManager;

class DISANAplotter {
 public:
DISANAplotter(
        ROOT::RDF::RNode                 df_dvcs_data,
        double                           beamEnergy,
        std::optional<ROOT::RDF::RNode>  df_pi0_data   = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcs_pi0mc    = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_pi0_pi0mc     = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_gen_dvcsmc = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_accept_dvcsmc = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_bkg = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_nobkg = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_rad = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_norad = std::nullopt
        )
        :rdf(std::move(df_dvcs_data)),
        rdf_pi0_data(std::move(df_pi0_data)),
        rdf_dvcs_pi0mc(std::move(df_dvcs_pi0mc)),
        rdf_pi0_pi0mc(std::move(df_pi0_pi0mc)),
        rdf_gen_dvcsmc(std::move(df_gen_dvcsmc)),
        rdf_accept_dvcsmc(std::move(df_accept_dvcsmc)),
        rdf_dvcsmc_bkg(std::move(df_dvcsmc_bkg)),
        rdf_dvcsmc_nobkg(std::move(df_dvcsmc_nobkg)),
        rdf_dvcsmc_rad(std::move(df_dvcsmc_rad)),
        rdf_dvcsmc_norad(std::move(df_dvcsmc_norad)),
        beam_energy(beamEnergy)
    {}

  ~DISANAplotter() {
    disHistos.clear();
    kinematicHistos.clear();
      // delete cached cross-section histos
    for (auto& row : phi_dsdt_QW_) {
      for (auto* h : row) { delete h; }
    }
    phi_dsdt_QW_.clear();

  }
   
ROOT::RDF::RNode GetRDF() { return rdf; }
// ---------------- NEW: phi mass-fit results & acceptance provider ---------------
  struct PhiMassFitResult {
    double mu{};     // peak position
    double sigma{};  // gaussian width
    double mLo{}, mHi{};
    ULong64_t Nwin{};  // raw counts in ±3σ
    double Nbkg{};     // fitted background under window
    double Nsig{};     // signal = Nwin - Nbkg
    double dNsig{};    // simple stat ⊕ bkg-fit uncertainty
  };

struct PhiMassBin {
    TH1D* hM = nullptr;     // invariant-mass histogram (owned, SetDirectory(0))
    TF1* fTotal = nullptr;  // total fit
    TF1* fSig = nullptr;    // Gaussian signal
    TF1* fBkg = nullptr;    // background
    double mu = 0, sigma = 0, mLo = 0, mHi = 0;
    double Nwin = 0, Nbkg = 0, Nsig = 0, dNsig = 0;
  };

using AccEffProvider = std::function<double(double q2lo, double q2hi, double wlo, double whi, double tlo, double thi)>;
// Default A×ε = 1.0 everywhere
void SetAccEffProvider(AccEffProvider f) { accEffProvider_ = std::move(f); }
void SetPlotApplyCorrection(bool apply) {dopi0corr = apply;}
void SetPlotApplyAcceptanceCorrection(bool apply) {doacceptcorr = apply;}
void SetPlotApplyEfficiencyCorrection(bool apply) {doefficiencycorr = apply;}
void SetPlotApplyRadiativeCorrection(bool apply) {doradiativecorr = apply;}
bool getDoPi0Corr() const { return dopi0corr; }
  // --- in public: add a const accessor ---
const std::vector<std::vector<TH1D*>>& GetPhiDSigmaDt3D() const { return phi_dsdt_QW_; }
std::vector<std::vector<TH1D*>>& GetPhiDSigmaDt3D() { return phi_dsdt_QW_; } // (optional mutable)


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
    if (doefficiencycorr) {
      auto eff3D = ComputeEffCorr(bins);
      result = UseEffCorrection(result, eff3D);
    }
    if (doradiativecorr) {
      auto rad3D = ComputeRadCorr(bins);
      result = UseRadCorrection(result, rad3D);
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

  std::vector<std::vector<std::vector<TH1D*>>> ComputeEffCorr(const BinManager& bins) {
    if (!rdf_dvcsmc_bkg || !rdf_dvcsmc_nobkg) {
      std::cerr << "[ComputeEffCorr] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcEfficiencyCorr(*rdf_dvcsmc_bkg, *rdf_dvcsmc_nobkg, bins);
  }

  std::vector<std::vector<std::vector<TH1D*>>> UseEffCorrection(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D,
                                                                  const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
    if (xs3D.size() != corr3D.size()) {
        std::cerr << "[UseEffCorrection] Dimension mismatch (level-0)\n";
        return {};
    }
    std::vector<std::vector<std::vector<TH1D*>>> xsCorr3D(xs3D.size());
    for (size_t iq = 0; iq < xs3D.size(); ++iq) {
        if (xs3D[iq].size() != corr3D[iq].size()) {
            std::cerr << "[UseEffCorrection] Dimension mismatch (level-1)\n";
            return {};
        }
        xsCorr3D[iq].resize(xs3D[iq].size());
        for (size_t it = 0; it < xs3D[iq].size(); ++it) {
            if (xs3D[iq][it].size() != corr3D[iq][it].size()) {
                std::cerr << "[UseEffCorrection] Dimension mismatch (level-2)\n";
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

  std::vector<std::vector<std::vector<TH1D*>>> ComputeRadCorr(const BinManager& bins) {
    if (!rdf_dvcsmc_rad || !rdf_dvcsmc_norad) {
      std::cerr << "[ComputeRadCorr] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcRadiativeCorr(*rdf_dvcsmc_rad, *rdf_dvcsmc_norad, bins);
  }
  std::vector<std::vector<std::vector<TH1D*>>> UseRadCorrection(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D,
                                                                  const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
    if (xs3D.size() != corr3D.size()) {
        std::cerr << "[UseRadCorrection] Dimension mismatch (level-0)\n";
        return {};
    }
    std::vector<std::vector<std::vector<TH1D*>>> xsCorr3D(xs3D.size());
    for (size_t iq = 0; iq < xs3D.size(); ++iq) {
        if (xs3D[iq].size() != corr3D[iq].size()) {
            std::cerr << "[UseRadCorrection] Dimension mismatch (level-1)\n";
            return {};
        }
        xsCorr3D[iq].resize(xs3D[iq].size());
        for (size_t it = 0; it < xs3D[iq].size(); ++it) {
            if (xs3D[iq][it].size() != corr3D[iq][it].size()) {
                std::cerr << "[UseRadCorrection] Dimension mismatch (level-2)\n";
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


//// phi analysis functions 
  // Header-only helper: fit mass in a filtered node and count ±3σ
  inline PhiMassFitResult FitPhiMassInNode(ROOT::RDF::RNode dfBin, const char* massCol = "invMass_KpKm", int nBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                           bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30) {
    static std::atomic<unsigned long> uid{0};
    PhiMassFitResult out;

    // 1) Make & materialize the histogram (unique name under MT)
    auto hname = Form("hM_tmp_%lu", uid.fetch_add(1));
    auto hR = dfBin.Histo1D(ROOT::RDF::TH1DModel(hname, ";M_{K^{+}K^{-}} [GeV];Counts", (unsigned)nBins, mMin, mMax), massCol);
    hR.GetValue();  // force event loop
    TH1D* h = hR.GetPtr();
    if (!h || h->GetEntries() == 0) return out;

    // 2) Total fit: thresholded background + Gaussian signal
    //    background threshold at ~ 2*m_K ≈ 0.9874
    TF1 fTot("fTot",
             "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874)) + "
             "[3]*TMath::Gaus(x,[4],[5])",
             mMin, mMax);
    fTot.SetParNames("A", "alpha", "lambda", "N", "mu", "sigma");

    // Initial guesses
    const double peakX = 1.0195;
    const double peakY = h->GetBinContent(h->FindBin(peakX));
    fTot.SetParameters(10.0, 0.9, 2.0, std::max(100.0, peakY), peakX, sigmaRef);
    fTot.SetParLimits(4, 0.990, 1.040);  // mu
    const double sLo = std::max(0.001, sigmaRef * (1.0 - sigmaFrac));
    const double sHi = sigmaRef * (1.0 + sigmaFrac);
    fTot.SetParLimits(5, constrainSigma ? sLo : 0.001, constrainSigma ? sHi : 0.030);

    // Quiet, range-only fit (no draw). No "I" option → function value ~ counts/bin
    h->Fit(&fTot, "QR0");  // Q=quiet, R=range, 0=no draw

    // 3) Extract mu/sigma and define ±3σ window
    out.mu = fTot.GetParameter(4);
    out.sigma = std::abs(fTot.GetParameter(5));
    out.mLo = out.mu - 3.0 * out.sigma;
    out.mHi = out.mu + 3.0 * out.sigma;

    // 4) Count data in ±3σ directly from the dataframe (keeps cuts consistent)
    auto dfWin = dfBin.Filter([&](float m) { return m > out.mLo && m < out.mHi; }, {massCol});
    out.Nwin = *dfWin.Count();

    // 5) Background under the window:
    //    Use the background component evaluated at bin centers and summed over bins in [mLo,mHi]
    TF1 fBkg("fBkg", "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))", mMin, mMax);
    fBkg.SetParameters(fTot.GetParameter(0), fTot.GetParameter(1), fTot.GetParameter(2));

    const int ibLo = h->FindBin(out.mLo);
    const int ibHi = h->FindBin(out.mHi);
    double Nbkg = 0.0;
    for (int b = ibLo; b <= ibHi; ++b) {
      const double xc = h->GetXaxis()->GetBinCenter(b);
      Nbkg += fBkg.Eval(xc);  // counts per bin (consistent with Fit without "I")
    }
    out.Nbkg = Nbkg;

    // 6) Signal and uncertainty (simple)
    out.Nsig = std::max(0.0, static_cast<double>(out.Nwin) - out.Nbkg);
    const double dNwin = std::sqrt(std::max(0.0, static_cast<double>(out.Nwin)));
    // crude bkg uncertainty: 10% of Nbkg (tune if you propagate TF1 errors)
    const double dNbkg = 0.10 * out.Nbkg;
    out.dNsig = std::hypot(dNwin, dNbkg);

    return out;
  }

  // Requires BinManager to provide GetQ2Bins(), GetTBins(), and (optionally) GetWBins()
  inline std::vector<std::vector<TH1D*>> ComputePhiDSigmaDt_Fitted(const BinManager& bins, double luminosity_nb_inv, double branching = 0.492, int nMassBins = 120,
                                                                   double mMin = 0.98, double mMax = 1.08, bool constrainSigma = false, double sigmaRef = 0.004,
                                                                   double sigmaFrac = 0.25, std::string modelName = "PhiModel") {
    const auto& q2Edges = bins.GetQ2Bins();
    const auto& tEdges = bins.GetTBins();

    // If not, leave hasW=false and wEdges empty.
    std::vector<double> wEdges;
    bool hasW = false;  

    std::vector<std::vector<TH1D*>> out;
    if (q2Edges.size() < 2 || tEdges.size() < 2) return out;

    out.resize(q2Edges.size() - 1);

    for (size_t iq = 0; iq + 1 < q2Edges.size(); ++iq) {
      const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];
      auto df_q = rdf.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      const size_t nW = hasW ? (wEdges.size() - 1) : 1;
      out[iq].resize(nW);

      for (size_t iw = 0; iw < nW; ++iw) {
        double wLo = 0.0, wHi = 1e9;
        if (hasW) {
          wLo = wEdges[iw];
          wHi = wEdges[iw + 1];
        }

        // Single filter that’s a no-op when hasW==false. No reassignments of RDF types.
        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        auto hName = Form("dsdt_%s_Q%zu_W%zu", modelName.c_str(), iq, iw);
        TH1D* h = new TH1D(hName, "d#sigma/dt; -t [GeV^{2}]; nb/GeV^{2}", static_cast<int>(tEdges.size() - 1), tEdges.data());
        h->SetDirectory(nullptr);  // independent of current file
        for (size_t it = 0; it + 1 < tEdges.size(); ++it) {
          const double tLo = tEdges[it], tHi = tEdges[it + 1];
          const double dT = tHi - tLo;

          auto df_bin = df_qw.Filter([=](double t) { return t > tLo && t <= tHi; }, {"t"});

          auto R = FitPhiMassInNode(df_bin, "invMass_KpKm", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);

          const double Aε = accEffProvider_(qLo, qHi, wLo, wHi, tLo, tHi);

          double val = 0.0, err = 0.0;
          if (luminosity_nb_inv > 0.0 && Aε > 0.0 && branching > 0.0 && dT > 0.0) {
            val = R.Nsig / (luminosity_nb_inv * Aε * branching * dT);
            err = R.dNsig / (luminosity_nb_inv * Aε * branching * dT);
          }

          const int b = static_cast<int>(it + 1);
          h->SetBinContent(b, val);
          h->SetBinError(b, err);
        }

        out[iq][iw] = h;
      }
    }

    return out;
  }
   // Your original wrapper (kept as-is)
  std::vector<TH1D*> ComputePhiCrossSection(const BinManager& bins, double luminosity) { return kinCalc.ComputePhi_CrossSection(rdf, bins, luminosity); }

  // === in public: ==============================================
  struct PhiMassDraw {
    TH1D* h{nullptr};
    TF1* fTot{nullptr};
    TF1* fSig{nullptr};
    TF1* fBkg{nullptr};
    double mu{0}, sigma{0}, mLo{0}, mHi{0};
    TCanvas* c{nullptr};
    std::string name;
  };

inline std::vector<std::vector<std::vector<PhiMassDraw>>> MakePhiMassFitCanvases3D(
    const BinManager& bins,
    const std::string& outDirPerModel,
    int nMassBins = 200,
    double mMin = 0.9874,
    double mMax = 1.120,
    bool constrainSigma = true,
    double sigmaRef = 0.004,
    double sigmaFrac = 0.30,
    // --- NEW (optional): compute dσ/dt when luminosity > 0 ---
    double luminosity_nb_inv = -1.0,
    double branching = 0.492)
{ 
  if (luminosity_nb_inv <= 0.0) {
  std::cerr << "[MakePhiMassFitCanvases3D] luminosity<=0 → will NOT compute dσ/dt cache.\n";
}
  // ensure directory exists
  gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

  const auto& q2Edges = bins.GetQ2Bins();
  const auto& tEdges  = bins.GetTBins();
  const auto& wEdges  = bins.GetWBins();
  const bool  hasW    = !wEdges.empty();

  const size_t nQ = q2Edges.size() > 1 ? q2Edges.size() - 1 : 0;
  const size_t nW = hasW ? (wEdges.size() - 1) : 1;
  const size_t nT = tEdges.size() > 1 ? tEdges.size() - 1 : 0;

  std::vector<std::vector<std::vector<PhiMassDraw>>> out;
  if (!nQ || !nT) return out;
  out.resize(nQ);

  // ---- NEW: prepare the 3D TH1 cache for dσ/dt (Q,W -> TH1 over t) ----
  phi_dsdt_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));
  static std::atomic<unsigned long> uid_dsdt{0};

  for (size_t iq = 0; iq < nQ; ++iq) {
    out[iq].resize(nW);

    const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];
    auto df_q = rdf.Filter([=](double Q2){ return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

    for (size_t iw = 0; iw < nW; ++iw) {
      const double wLo = hasW ? wEdges[iw]     : 0.0;
      const double wHi = hasW ? wEdges[iw + 1] : 1e9;
      auto df_qw = df_q.Filter([=](double W){ return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

      out[iq][iw].resize(nT);

      // --- NEW: create the dσ/dt histogram for this (Q2,W) slice (if we will compute) ---
      if (luminosity_nb_inv > 0.0) {
        auto hname = Form("phi_dsdt_Q%zu_W%zu_%lu", iq, iw, uid_dsdt.fetch_add(1));
        auto htitle = Form("d#sigma/dt; -t [GeV^{2}]; nb/GeV^{2}");
        TH1D* hDsDt = new TH1D(hname, htitle, static_cast<int>(tEdges.size() - 1), tEdges.data());
        hDsDt->SetDirectory(nullptr);
        phi_dsdt_QW_[iq][iw] = hDsDt;
      }

      for (size_t it = 0; it < nT; ++it) {
        const double tLo = tEdges[it], tHi = tEdges[it + 1];
        const double dT  = tHi - tLo;

        auto df_bin = df_qw.Filter([=](double t){ return t > tLo && t <= tHi; }, {"t"});

        static std::atomic<unsigned long> uid{0};
        auto hname = Form("hM_%zu_%zu_%zu_%lu", iq, iw, it, uid.fetch_add(1));
        auto hR = df_bin.Histo1D(ROOT::RDF::TH1DModel(hname, ";M_{K^{+}K^{-}} [GeV];Counts",
                                                      (unsigned)nMassBins, mMin, mMax),
                                 "invMass_KpKm");
        hR.GetValue();
        TH1D* h = (TH1D*)hR.GetPtr();
        if (!h || h->GetEntries() == 0) {
          // still record an empty point in dσ/dt if we're computing it
          if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
            const int b = static_cast<int>(it + 1);
            phi_dsdt_QW_[iq][iw]->SetBinContent(b, 0.0);
            phi_dsdt_QW_[iq][iw]->SetBinError  (b, 0.0);
          }
          continue;
        }

        // ---- (existing) pretty clone for drawing ----
        TH1D* hDraw = (TH1D*)h->Clone(Form("%s_draw", h->GetName()));
        hDraw->SetDirectory(0);
        hDraw->SetTitle("");
        hDraw->SetLineColor(kBlue + 1);
        hDraw->SetLineWidth(2);
        hDraw->SetFillColorAlpha(kBlue - 9, 0.3);
        hDraw->SetMarkerStyle(20);
        hDraw->SetMarkerSize(1.2);
        hDraw->SetMarkerColor(kBlue + 2);
        hDraw->GetXaxis()->SetTitle("M(K^{+}K^{-}) [GeV]");
        hDraw->GetYaxis()->SetTitle("Counts");
        hDraw->GetXaxis()->SetRangeUser(mMin - 0.02, mMax);

        // ---- (existing) total fit: thresholded bkg + gaussian signal ----
        TF1* fTot = new TF1(Form("fTot_%s", hname),
          "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874)) + [3]*TMath::Gaus(x,[4],[5])",
          mMin, mMax);
        fTot->SetParNames("A","alpha","lambda","N","mu","sigma");
        const double peakX = 1.0195;
        const double peakY = 1.25 * (h->GetBinContent(h->FindBin(peakX)));
        fTot->SetParameters(10.0, 0.9, 2.0, std::max(1.0, peakY), peakX, sigmaRef);
        fTot->SetParLimits(4, 0.990, 1.045);
        const double sLo = std::max(0.001, sigmaRef * (1.0 - sigmaFrac));
        const double sHi = sigmaRef * (1.0 + sigmaFrac);
        (void)sLo; (void)sHi; // we leave a wide free sigma in this drawer
        fTot->SetParLimits(5, 0.0001, 0.040);
        hDraw->Fit(fTot, "QR0");

        TF1* fSig = new TF1(Form("fSig_%s", hname), "[0]*TMath::Gaus(x,[1],[2])", mMin, mMax);
        fSig->SetParameters(fTot->GetParameter(3), fTot->GetParameter(4), fTot->GetParameter(5));
        fSig->SetLineColor(kOrange + 1);
        fSig->SetLineStyle(2);
        fSig->SetLineWidth(2);
        fSig->SetFillColorAlpha(kOrange - 3, 0.3);
        fSig->SetFillStyle(1001);

        TF1* fBkg = new TF1(Form("fBkg_%s", hname),
          "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))", mMin, mMax);
        fBkg->SetParameters(fTot->GetParameter(0), fTot->GetParameter(1), fTot->GetParameter(2));
        fBkg->SetLineColor(kGreen + 2);
        fBkg->SetLineStyle(3);
        fBkg->SetLineWidth(2);
        fBkg->SetFillColorAlpha(kGreen - 7, 0.3);
        fBkg->SetFillStyle(1001);

        const double mu    = fTot->GetParameter(4);
        const double sigma = std::abs(fTot->GetParameter(5));
        const double mLo   = mu - 3.0 * sigma;
        const double mHi   = mu + 3.0 * sigma;

        // background under ±3σ (bin-center sum; consistent with fit mode)
        double Nbkg = 0.0;
        const int ibLo = hDraw->FindBin(mLo);
        const int ibHi = hDraw->FindBin(mHi);
        for (int b = ibLo; b <= ibHi; ++b) {
          Nbkg += fBkg->Eval(hDraw->GetXaxis()->GetBinCenter(b));
        }

        // window counts & signal
        auto dfWin     = df_bin.Filter([=](float m){ return m > mLo && m < mHi; }, {"invMass_KpKm"});
        const ULong64_t Nwin = *dfWin.Count();
        const double Nsig    = std::max(0.0, static_cast<double>(Nwin) - Nbkg);
        const double dNsig   = std::hypot(std::sqrt(std::max(0.0, (double)Nwin)), 0.10 * Nbkg);

        // ---- NEW: fill dσ/dt if requested ----
        if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
          const double Aeps = accEffProvider_(qLo, qHi, wLo, wHi, tLo, tHi); // user-supplied or default 1
          const int    b    = static_cast<int>(it + 1);
          double val = 0.0, err = 0.0;
          if (Aeps > 0.0 && branching > 0.0 && dT > 0.0) {
            val = Nsig  / (luminosity_nb_inv * Aeps * branching * dT);
            err = dNsig / (luminosity_nb_inv * Aeps * branching * dT);
          }
          phi_dsdt_QW_[iq][iw]->SetBinContent(b, val);
          phi_dsdt_QW_[iq][iw]->SetBinError(b, err);
        }

        // ---- (existing) canvas beautification & save ----
        TCanvas* c = new TCanvas(Form("c_%s", hname), "K^{+}K^{-} mass", 1200, 1000);
        gStyle->Reset();
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleSize(0.05, "XYZ");
        gStyle->SetLabelSize(0.04, "XYZ");
        c->SetTicks(1,1);
        c->SetTopMargin(0.04);
        c->SetRightMargin(0.03);
        c->SetBottomMargin(0.11);
        c->SetLeftMargin(0.12);
        gStyle->SetCanvasColor(0);
        gStyle->SetPadColor(0);
        c->SetFillColor(0);
        if (c->GetPad(0)) c->GetPad(0)->SetFillColor(0);
        c->cd();
        hDraw->SetMinimum(0.0);  
        hDraw->Draw("PE");
        fTot->SetLineColor(kRed + 1);
        fBkg->Draw("SAME FC");
        fSig->Draw("SAME FC");
        fTot->Draw("SAME C");

        double yMax = hDraw->GetMaximum();
        TLine* L1 = new TLine(mLo, 0.0, mLo, yMax * 0.75);
        TLine* L2 = new TLine(mHi, 0.0, mHi, yMax * 0.75);
        L1->SetLineColor(kMagenta + 2);
        L2->SetLineColor(kMagenta + 2);
        L1->SetLineStyle(2); L2->SetLineStyle(2);
        L1->SetLineWidth(2); L2->SetLineWidth(2);
        L1->Draw("SAME");    L2->Draw("SAME");

        TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hDraw, "K^{+}K^{-} Inv. Mass", "lep");
        leg->AddEntry(fTot,  "Total Fit: Exp + Gauss", "l");
        leg->AddEntry(fSig,  "Signal (Gauss)", "f");
        leg->AddEntry(fBkg,  "Background (Exp)", "f");
        leg->AddEntry(L1,    "#pm 3#sigma", "l");
        leg->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.15, 0.87, Form("#mu = %.3f GeV", mu));
        latex.DrawLatex(0.15, 0.82, Form("#sigma = %.4f GeV", sigma));
        const double chi2ndf = fTot->GetNDF() > 0 ? fTot->GetChisquare()/fTot->GetNDF() : 0.0;
        latex.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf = %.2f", chi2ndf));
        latex.DrawLatex(0.15, 0.72, Form("Fit range: %.3f-%.3f GeV", mMin, mMax));
        latex.DrawLatex(0.15, 0.67, Form("Q^{2}[%.2f,%.2f]  %s  -t[%.2f,%.2f]",
                               qLo, qHi, hasW ? Form("W[%.1f,%.1f]", wLo, wHi) : "", tLo, tHi));
        latex.DrawLatex(0.15, 0.62, Form("N_{win}= %llu, N_{bkg}= %.1f, N_{sig}= %.1f",
                               (unsigned long long)Nwin, Nbkg, Nsig));

        auto tagRaw = hasW ? Form("Q2_%.3f_%.3f__W_%.3f_%.3f__t_%.3f_%.3f", qLo,qHi,wLo,wHi,tLo,tHi)
                           : Form("Q2_%.3f_%.3f__t_%.3f_%.3f", qLo,qHi,tLo,tHi);
        std::string tag = tagRaw; std::replace(tag.begin(), tag.end(), '.', '_');
        c->SaveAs(Form("%s/KKmass_%s.pdf", outDirPerModel.c_str(), tag.c_str()));

        // pack draw products
        PhiMassDraw pack;
        pack.h = hDraw; pack.fTot = fTot; pack.fSig = fSig; pack.fBkg = fBkg;
        pack.mu = mu; pack.sigma = sigma; pack.mLo = mLo; pack.mHi = mHi;
        pack.c = c; pack.name = tag;
        out[iq][iw][it] = pack;

        // clean transient objects we don't cache (canvas snapshot saved already)
        delete c;
        //delete hDraw;
        if (gPad) { gPad->Clear(); gPad->Update(); gPad->SetSelected(nullptr); }
        delete leg;
        delete L1;
        delete L2;
      } // it (t bins)
    } // iw (W bins)
  } // iq (Q2 bins)
  return out;
}

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
  bool doefficiencycorr = false;
  bool doradiativecorr = false;
  std::string ttreeName;
  std::vector<ROOT::RDF::RResultPtr<TH1>> kinematicHistos, disHistos;
  std::vector<std::vector<TH1D*>> phi_dsdt_QW_;
  std::vector<std::shared_ptr<TH1>> acceptHistos;
  ROOT::RDF::RNode rdf;
  
  //ROOT::RDF::RNode rdf_dvcs_data;
  std::optional<ROOT::RDF::RNode> rdf_pi0_data;
  std::optional<ROOT::RDF::RNode> rdf_dvcs_pi0mc;
  std::optional<ROOT::RDF::RNode> rdf_pi0_pi0mc;
  std::optional<ROOT::RDF::RNode> rdf_gen_dvcsmc;
  std::optional<ROOT::RDF::RNode> rdf_accept_dvcsmc;
  std::optional<ROOT::RDF::RNode> rdf_dvcsmc_bkg;
  std::optional<ROOT::RDF::RNode> rdf_dvcsmc_nobkg;
  std::optional<ROOT::RDF::RNode> rdf_dvcsmc_rad;
  std::optional<ROOT::RDF::RNode> rdf_dvcsmc_norad;

  AccEffProvider accEffProvider_ = [](double, double, double, double, double, double) { return 1.0; };
};

#endif  // DISANA_PLOTTER_H
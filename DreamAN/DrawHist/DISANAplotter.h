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
        double                          luminosity,
        std::optional<ROOT::RDF::RNode>  df_pi0_data   = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcs_pi0mc    = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_pi0_pi0mc     = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_gen_dvcsmc = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_accept_dvcsmc = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_bkg = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_nobkg = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_rad = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_norad = std::nullopt,
        std::optional<ROOT::RDF::RNode>  df_dvcsmc_p1cut = std::nullopt
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
        rdf_dvcsmc_p1cut(std::move(df_dvcsmc_p1cut)),
        beam_energy(beamEnergy),
        luminosity_nb_inv(luminosity)
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
void SetLuminosity(double luminosity_nb_inv) { luminosity_nb_inv = luminosity_nb_inv; }
ROOT::RDF::RNode GetRDF() { return rdf; }
ROOT::RDF::RNode GetRDF_Pi0Data() { return rdf_pi0_data.value(); }
ROOT::RDF::RNode GetRDF_DVCSPi0MC() { return rdf_pi0_pi0mc.value(); }
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
void SetPlotApplyP1Cut(bool apply) {dop1cut = apply;}
bool getDoPi0Corr() const { return dopi0corr; }
  // --- in public: add a const accessor ---
const std::vector<std::vector<TH1D*>>& GetPhiDSigmaDt3D() const { return phi_dsdt_QW_; }
std::vector<std::vector<TH1D*>>& GetPhiDSigmaDt3D() { return phi_dsdt_QW_; } // (optional mutable)


void GenerateKinematicHistos(const std::string& type) {
    std::vector<std::string> vars = {"p", "theta", "phi"}; // for DVCS analysis
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

void GeneratePhiKinematicHistos(const std::string& type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi", "vz"}; //for Phi analysis
    //std::vector<std::string> vars = {"p", "theta", "phi"}; // for DVCS analysis
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

  for (const auto& var : Phivars) {
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


void GeneratePi0KinematicHistos(const std::string& type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi", "vz"}; //for Phi analysis
    //std::vector<std::string> vars = {"p", "theta", "phi"}; // for DVCS analysis
    for (const auto& v : vars) {
      std::string base = "rec" + type + "_" + v;
      std::string basename = "pi0_" + base;
      //std::cout << " base is " << base << std::endl;
      double histMin, histMax;
      auto it = kinematicAxisRanges.find(base);
      if (it != kinematicAxisRanges.end()) {
        histMin = it->second.first;
        histMax = it->second.second;
      } else {
        auto df_tmp = (v == "p")
            ? rdf_pi0_data->Filter(base + " > -100")
            : *rdf_pi0_data;
        auto lo = *df_tmp.Min(base);
        auto hi = *df_tmp.Max(base);
        double margin = std::max(1e-3, 0.05 * (hi - lo));
        histMin = lo - margin;
        histMax = hi + margin;
        kinematicAxisRanges[base] = {histMin, histMax};  // store it
      }
    
       // auto h_all = rdf.Histo1D({(base + "_all").c_str(), "", 100, histMin, histMax}, base);
    auto h_acc = rdf_pi0_data->Histo1D({(basename).c_str(), "", 100, histMin, histMax}, base);
    kinematicHistos.push_back(h_acc);
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

  std::vector<std::vector<std::vector<TH1D*>>> ComputeDVCS_CrossSection(const BinManager& bins) {
    auto result = kinCalc.ComputeDVCS_CrossSection(rdf, bins, luminosity_nb_inv);
    if (dopi0corr) {
      auto sigma_pi0_3d = kinCalc.ComputeDVCS_CrossSection(*rdf_pi0_data, bins, luminosity_nb_inv);
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
    if (dop1cut) {
      auto p13D = ComputeP1CutEffect(bins);
      result = UseP1Cut(result, p13D);
    }
    return result;
  }

  /// BSA computations // we need to have refined version of this codes
  std::vector<std::vector<std::vector<TH1D*>>> ComputeBSA(const BinManager& bins,  double pol = 1.0) {
    // 1) Helicity selection
    auto rdf_pos = rdf.Filter("REC_Event_helicity ==  1");
    auto rdf_neg = rdf.Filter("REC_Event_helicity == -1");
    // 2) One-pass cross-section per helicity

    auto sigma_pos_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pos, bins, luminosity_nb_inv);
    auto sigma_neg_3d = kinCalc.ComputeDVCS_CrossSection(rdf_neg, bins, luminosity_nb_inv);

    auto result = kinCalc.ComputeBeamSpinAsymmetry(sigma_pos_3d, sigma_neg_3d, pol);

    if (dopi0corr) {
      auto rdf_pi0_data_pos = rdf_pi0_data->Filter("REC_Event_helicity ==  1");
      auto rdf_pi0_data_neg = rdf_pi0_data->Filter("REC_Event_helicity == -1");
    // 2) One-pass cross-section per helicity

      auto sigma_pi0_pos_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pi0_data_pos, bins, luminosity_nb_inv);
      auto sigma_pi0_neg_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pi0_data_neg, bins, luminosity_nb_inv);
      
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

  std::vector<std::vector<std::vector<TH1D*>>> ComputePi0DVCSdiffmc(const BinManager& bins) {
    if (!rdf_dvcs_pi0mc || !rdf_pi0_pi0mc) {
      std::cerr << "[ComputePi0DVCSdiffmc] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcAcceptanceCorr(*rdf_pi0_pi0mc, *rdf_dvcs_pi0mc, bins);

  }

  std::vector<std::vector<std::vector<TH1D*>>> ComputePi0DVCSdiffexp(const BinManager& bins) {
    if (!rdf_pi0_data) {
      std::cerr << "[ComputePi0DVCSdiffexp] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcAcceptanceCorr(rdf, *rdf_pi0_data, bins);

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

  std::vector<std::vector<std::vector<TH1D*>>> ComputeP1CutEffect(const BinManager& bins) {
    if (!rdf_dvcsmc_p1cut || !rdf_dvcsmc_norad) {
      std::cerr << "[ComputeP1CutEffect] Missing input RDFs.\n";
      return {};
    }
    return kinCalc.CalcP1Cut(*rdf_dvcsmc_p1cut, *rdf_dvcsmc_norad, bins);
  }

  std::vector<std::vector<std::vector<TH1D*>>> UseP1Cut(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D,
                                                                  const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
    if (xs3D.size() != corr3D.size()) {
        std::cerr << "[UseP1CutCorrection] Dimension mismatch (level-0)\n";
        return {};
    }
    std::vector<std::vector<std::vector<TH1D*>>> xsCorr3D(xs3D.size());
    for (size_t iq = 0; iq < xs3D.size(); ++iq) {
        if (xs3D[iq].size() != corr3D[iq].size()) {
            std::cerr << "[UseP1CutCorrection] Dimension mismatch (level-1)\n";
            return {};
        }
        xsCorr3D[iq].resize(xs3D[iq].size());
        for (size_t it = 0; it < xs3D[iq].size(); ++it) {
            if (xs3D[iq][it].size() != corr3D[iq][it].size()) {
                std::cerr << "[UseP1CutCorrection] Dimension mismatch (level-2)\n";
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
                    double val_corr = xs_val;
                    double err_corr = xs_err;
                    if (c_val <= 0.9) val_corr = -1;
                    if (c_val <= 0.9) err_corr = -1;
                    if (c_val <= 0.9) std::cout << " P1 cut applied: c_val " << c_val << "at bin " << b << std::endl;
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

 //std::vector<TH1D*> ComputePhiCrossSection(const BinManager& bins) { return kinCalc.ComputePhi_CrossSection(rdf, bins, luminosity_nb_inv); }

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
    double branching = 0.492)
{ 
  if (luminosity_nb_inv <= 0.0) {
  std::cerr << "[MakePhiMassFitCanvases3D] luminosity<=0 → will NOT compute dσ/dt cache.\n";
}
  // ensure directory exists
  gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

  const auto& q2Edges = bins.GetQ2Bins();
  const auto& tPrimeEdges  = bins.GetTprimeBins();
  const auto& wEdges  = bins.GetWBins();
  const bool  hasW    = !wEdges.empty();

  const size_t nQ = q2Edges.size() > 1 ? q2Edges.size() - 1 : 0;
  const size_t nW = hasW ? (wEdges.size() - 1) : 1;
  const size_t nT = tPrimeEdges.size() > 1 ? tPrimeEdges.size() - 1 : 0;

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
        TH1D* hDsDt = new TH1D(hname, htitle, static_cast<int>(tPrimeEdges.size() - 1), tPrimeEdges.data());
        hDsDt->SetDirectory(nullptr);
        phi_dsdt_QW_[iq][iw] = hDsDt;
      }

      for (size_t it = 0; it < nT; ++it) {
        const double tLo = tPrimeEdges[it], tHi = tPrimeEdges[it + 1];
        const double dT  = tHi - tLo;
        const double dQ2 = qHi - qLo;
        const double dW  = wHi - wLo;
        double binVol =  dQ2 * dW * dT;

       // auto df_bin = df_qw.Filter([=](double t){ return t > tLo && t <= tHi; }, {"t"});
        auto df_bin = df_qw.Filter([=](double mtp){ return mtp > tLo && mtp <= tHi; }, {"mtprime"});

        static std::atomic<unsigned long> uid{0};
        auto hname = Form("hM_%zu_%zu_%zu_%lu", iq, iw, it, uid.fetch_add(1));
        auto hR = df_bin.Histo1D(ROOT::RDF::TH1DModel(hname, ";M_{K^{+}K^{-}} [GeV];Counts", 200, 0.8, 1.8),
                                 "invMass_KpKm");
        hR.GetValue();
        TH1D* h = (TH1D*)hR.GetPtr();
        if (!h || h->GetEntries() <=20) {
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
        hDraw->GetYaxis()->SetTitleOffset(1.2);
        hDraw->SetFillColorAlpha(kBlue - 9, 0.3);
        hDraw->SetMarkerStyle(20);
        hDraw->SetMarkerSize(1.2);
        hDraw->SetMarkerColor(kBlue + 2);
        hDraw->GetXaxis()->SetTitle("M(K^{+}K^{-}) [GeV]");
        hDraw->GetYaxis()->SetTitle("Counts");
        hDraw->GetXaxis()->SetRangeUser(mMin - 0.02, mMax+.020);

        const double bw = hDraw->GetXaxis()->GetBinWidth(1);
        std::cout << Form("[MakePhiMassFitCanvases3D] Fitting %s with %.0f bins, M=[%.4f,%.4f], BW=%.4f GeV",
                         hDraw->GetName(), (double)hDraw->GetNbinsX(), mMin, mMax, bw) << std::endl;
        
        std::string formula = Form(
            // Background: [0]*(x-0.9874)^[1] * exp(-[2]*(x-0.9874))
            "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))"
            " + "
            // Signal (counts/bin): [3] (=Nsig) * bw * Gaussian_pdf(x; [4]=mean, [5]=sigma)
            "[3]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[4])/([5]),2))/TMath::Abs([5])",
            bw);

        TF1* fTot = new TF1(Form("fitTotal_%s",hname), formula.c_str(), mMin, mMax);

        //[0]*0.398942*bw*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])
        // params: [0]=Abkg(scale), [1]=alpha, [2]=lambda, [3]=Nsig(yield), [4]=mu, [5]=sigma
        fTot->SetParNames("A","alpha","lambda","N","mu","sigma");
        fTot->SetParameters(4, 0.9, 2, 10, 1.02, 0.010);
        fTot->SetParLimits(4, 1.005, 1.022);
        fTot->SetParLimits(5, 0.0025, 0.025);
        fTot->SetParLimits(3, 0.000001, 100000.0);
        fTot->SetLineColor(kRed + 1);
        fTot->SetLineWidth(3);
        hDraw->Fit(fTot, "R0QL");
        fTot->Draw("SAME C");

        const double A = fTot->GetParameter(0);
        const double alpha = fTot->GetParameter(1);
        const double lambda = fTot->GetParameter(2);
        const double N = fTot->GetParameter(3);
        const double mu = fTot->GetParameter(4);
        const double sigma = fTot->GetParameter(5);
        const double chi2 = fTot->GetChisquare();
        const double ndf = fTot->GetNDF();

        const double mLo   = mu - 3.0 * sigma;
        const double mHi   = mu + 3.0 * sigma;

        TF1* fSig = new TF1(Form("fSig_%s",hname), Form("[0]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[1])/([2]),2))/TMath::Abs([2])",bw), mMin, mMax);
        fSig->SetParameters(N, mu, sigma);
        fSig->SetLineColor(kOrange + 1);
        fSig->SetLineStyle(2);
        fSig->SetLineWidth(2);
        fSig->SetFillColorAlpha(kOrange - 3, 0.3);
        fSig->SetFillStyle(1001);
        fSig->Draw("SAME FC");
      
        TF1* fBkg = new TF1(Form("fBkg_%s",hname),
                              "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))",
                              mMin, mMax);
        fBkg->SetParameters(A, alpha, lambda);
        fBkg->SetLineColor(kGreen + 2);
        fBkg->SetLineStyle(3);
        fBkg->SetLineWidth(2);
        fBkg->SetFillColorAlpha(kGreen - 7, 0.3);
        fBkg->SetFillStyle(1001);
        fBkg->Draw("SAME FC");
        // background under ±3σ (bin-center sum; consistent with fit mode)
        double Nbkg = 0.0;
        const double Nsig = fTot->GetParameter(3);
        const double Nsig_err = fTot->GetParError(3);

       // window counts & signal
        auto dfWin     = df_bin.Filter([=](float m){ return m > mLo && m < mHi; }, {"invMass_KpKm"});
        const ULong64_t Nwin = *dfWin.Count();
      
        // ---- NEW: fill dσ/dt if requested ----
        if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
          const double Aeps = accEffProvider_(qLo, qHi, wLo, wHi, tLo, tHi); // user-supplied or default 1
          const int    b    = static_cast<int>(it + 1);
          double val = 0.0, err = 0.0;
          if (Aeps > 0.0 && branching > 0.0 && dT > 0.0) {
            val = Nsig  / (luminosity_nb_inv * Aeps * branching * binVol);
            err = Nsig_err / (luminosity_nb_inv * Aeps * branching * binVol);
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
        c->SetLeftMargin(0.13);
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

        TLegend* leg = new TLegend(0.55, 0.5, 0.85, 0.75);
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
        latex.DrawLatex(0.45, 0.89, Form("Q^{2}[%.2f,%.2f]  %s  -t'[%.2f,%.2f]", qLo, qHi, hasW ? Form("W[%.1f,%.1f]", wLo, wHi) : "", tLo, tHi));
        latex.DrawLatex(0.55, 0.85, Form("#mu = %.3f GeV, #sigma = %.3f GeV", mu, sigma));
        const double chi2ndf = fTot->GetNDF() > 0 ? fTot->GetChisquare()/fTot->GetNDF() : 0.0;
        latex.DrawLatex(0.55, 0.81, Form("#chi^{2}/ndf = %.2f", chi2ndf));
        latex.DrawLatex(0.55, 0.77, Form("N_{#phi} (total) = %.1f #pm %.1f", Nsig, Nsig_err));

        auto tagRaw = hasW ? Form("Q2_%.1f_%.1f__W_%.1f_%.2f__tprime_%.1f_%.1f", qLo,qHi,wLo,wHi,tLo,tHi)
                           : Form("Q2_%.1f_%.1f__tprime_%.1f_%.1f", qLo,qHi,tLo,tHi);
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
  std::vector<std::string> Phivars = {"Q2", "xB", "t", "W", "phi","mtprime"};
  std::map<std::string, std::pair<double, double>> axisRanges = {{"Q2", {0.0, 15.0}}, {"xB", {0.0, 1.0}}, {"W", {1.0, 10.0}}, {"t", {0.0, 10.0}}, {"mtprime", {0.0, 10.0}}, {"phi", {-180.0, 180.0}}};
  std::map<std::string, std::pair<double, double>> kinematicAxisRanges = {
      {"recel_p", {-0.05, 13.0}},    {"recel_theta", {-0.01, 1.0}}, {"recel_phi", {-0.01, 6.1}}, {"recpho_p", {-0.01, 10.0}}, {"recpho2_p", {-0.01, 4.0}},
      {"recpho_theta", {-0.01, 1.0}}, {"recpho_phi", {-0.01, 6.1}},   {"recpro_p", {-0.01, 8.0}},  {"recpro_theta", {-0.01, 2.0}},
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
  bool dop1cut = false;
  std::string ttreeName;
  double luminosity_nb_inv = -1.0;
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
  std::optional<ROOT::RDF::RNode> rdf_dvcsmc_p1cut;

  AccEffProvider accEffProvider_ = [](double, double, double, double, double, double) { return 1.0; };
};

#endif  // DISANA_PLOTTER_H
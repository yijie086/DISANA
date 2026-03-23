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
struct DVCSModeTag {};
struct PhiModeTag {};

class DISANAplotter {
 public:
  // for DVCS analysis
  DISANAplotter(DVCSModeTag, ROOT::RDF::RNode df_dvcs_data, double beamEnergy, double luminosity, std::optional<ROOT::RDF::RNode> df_pi0_data = std::nullopt,
                std::optional<ROOT::RDF::RNode> df_dvcs_pi0mc = std::nullopt, std::optional<ROOT::RDF::RNode> df_pi0_pi0mc = std::nullopt,
                std::optional<ROOT::RDF::RNode> df_gen_dvcsmc = std::nullopt, std::optional<ROOT::RDF::RNode> df_accept_dvcsmc = std::nullopt,
                std::optional<ROOT::RDF::RNode> df_dvcsmc_bkg = std::nullopt, std::optional<ROOT::RDF::RNode> df_dvcsmc_nobkg = std::nullopt,
                std::optional<ROOT::RDF::RNode> df_dvcsmc_rad = std::nullopt, std::optional<ROOT::RDF::RNode> df_dvcsmc_norad = std::nullopt,
                std::optional<ROOT::RDF::RNode> df_dvcsmc_p1cut = std::nullopt)
      : rdf(std::move(df_dvcs_data)),
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
        luminosity_nb_inv(luminosity) {}

  // for Phi Analysis
  DISANAplotter(PhiModeTag, ROOT::RDF::RNode df_phi_data, double beamEnergy, double luminosity, std::optional<ROOT::RDF::RNode> df_gen_phimc = std::nullopt,
                std::optional<ROOT::RDF::RNode> df_accept_phimc = std::nullopt, std::optional<ROOT::RDF::RNode> df_phimc_radRatio = std::nullopt)
      : rdf(std::move(df_phi_data)),
        rdf_gen_phimc(std::move(df_gen_phimc)),
        rdf_accept_phimc(std::move(df_accept_phimc)),
        rdf_phimc_radRatio(std::move(df_phimc_radRatio)),
        beam_energy(beamEnergy),
        luminosity_nb_inv(luminosity) {}

  ~DISANAplotter() {
    disHistos.clear();
    kinematicHistos.clear();
    // delete cached cross-section histos
    for (auto& row : phi_dsdt_QW_) {
      for (auto* h : row) {
        delete h;
      }
    }
    phi_dsdt_QW_.clear();

    for (auto& row : phi_nsig_QW_) {
      for (auto* h : row) delete h;
    }
    phi_nsig_QW_.clear();

    for (auto& row : phi_alu_costh_QW_) {
      for (auto* h : row) {
        delete h;
      }
    }
    phi_alu_costh_QW_.clear();
    for (auto& row : phi_alu_cos_QW_) {
      for (auto* h : row) {
        delete h;
      }
    }
    phi_alu_cos_QW_.clear();

    for (auto& row : phi_accept_QW_) {
      for (auto* h : row) delete h;
    }
    phi_accept_QW_.clear();
    for (auto& row : phi_rad_QW_)
      for (auto* h : row) delete h;
    phi_rad_QW_.clear();

    for (auto& row : phi_bsa_trentophi_QW_) {
      for (auto* h : row) {
        delete h;
      }
    }
    phi_bsa_trentophi_QW_.clear();
    for (auto& row : phi_alu_zphi_QW_)
      for (auto* h : row) delete h;
    phi_alu_zphi_QW_.clear();
  }
  void SetLuminosity(double L) { luminosity_nb_inv = L; }
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

  using AccEffProvider = std::function<double(double q2lo, double q2hi, double wlo, double whi, double tlo, double thi)>;
  // Default A×ε = 1.0 everywhere
  void SetAccEffProvider(AccEffProvider f) { accEffProvider_ = std::move(f); }
  void SetPlotApplyCorrection(bool apply) { dopi0corr = apply; }
  void SetPlotApplyAcceptanceCorrection(bool apply) { doacceptcorr = apply; }
  void SetPlotApplyEfficiencyCorrection(bool apply) { doefficiencycorr = apply; }
  void SetPlotApplyRadiativeCorrection(bool apply) { doradiativecorr = apply; }
  void SetPlotApplyP1Cut(bool apply) { dop1cut = apply; }
  bool getDoPi0Corr() const { return dopi0corr; }
  // --- in public: add a const accessor ---
  const std::vector<std::vector<TH1D*>>& GetPhiDSigmaDt3D() const { return phi_dsdt_QW_; }
  std::vector<std::vector<TH1D*>>& GetPhiDSigmaDt3D() { return phi_dsdt_QW_; }  // (optional mutable)
  const std::vector<std::vector<TH1D*>>& GetPhiRawCounts3D() const { return phi_nsig_QW_; }
  std::vector<std::vector<TH1D*>>& GetPhiRawCounts3D() { return phi_nsig_QW_; }

  // --- in public: add a const accessor ---
  const std::vector<std::vector<TH1D*>>& GetPhiBSATrentoPhi3D() const { return phi_bsa_trentophi_QW_; }
  std::vector<std::vector<TH1D*>>& GetPhiBSATrentoPhi3D() { return phi_bsa_trentophi_QW_; }

  const std::vector<std::vector<TH1D*>>& GetPhiALUCosTheta3D() const { return phi_alu_cos_QW_; }
  std::vector<std::vector<TH1D*>>& GetPhiALUCosTheta3D() { return phi_alu_cos_QW_; }
  const std::vector<std::vector<TH1D*>>& GetPhiALUZPhi3D() const { return phi_alu_zphi_QW_; }
  const std::vector<std::vector<TH1D*>>& GetPhiAcceptance3D() const { return phi_accept_QW_; }
  std::vector<std::vector<TH1D*>>& GetPhiAcceptance3D() { return phi_accept_QW_; }
  const std::vector<std::vector<TH1D*>>& GetPhiRadCorr3D() const { return phi_rad_QW_; }
  std::vector<std::vector<TH1D*>>& GetPhiRadCorr3D() { return phi_rad_QW_; }

 double GetPhiMeanXB(size_t iq, size_t iw, size_t it) const {
    if (iq >= phi_mean_xB_QW_.size()) return std::numeric_limits<double>::quiet_NaN();
    if (iw >= phi_mean_xB_QW_[iq].size()) return std::numeric_limits<double>::quiet_NaN();
    if (it >= phi_mean_xB_QW_[iq][iw].size()) return std::numeric_limits<double>::quiet_NaN();
    return phi_mean_xB_QW_[iq][iw][it];
  }
  double GetPhiMeanW(size_t iq, size_t iw, size_t it) const {
    if (iq >= phi_mean_W_QW_.size()) return std::numeric_limits<double>::quiet_NaN();
    if (iw >= phi_mean_W_QW_[iq].size()) return std::numeric_limits<double>::quiet_NaN();
    if (it >= phi_mean_W_QW_[iq][iw].size()) return std::numeric_limits<double>::quiet_NaN();
    return phi_mean_W_QW_[iq][iw][it];
  }
  double GetPhiMeanGammaV(size_t iq, size_t iw, size_t it) const {
    if (iq >= phi_mean_GammaV_QW_.size()) return std::numeric_limits<double>::quiet_NaN();
    if (iw >= phi_mean_GammaV_QW_[iq].size()) return std::numeric_limits<double>::quiet_NaN();
    if (it >= phi_mean_GammaV_QW_[iq][iw].size()) return std::numeric_limits<double>::quiet_NaN();
    return phi_mean_GammaV_QW_[iq][iw][it];
  }


  void GenerateKinematicHistos(const std::string& type) {
    std::vector<std::string> vars = {"p", "theta", "phi"};  // for DVCS analysis
    for (const auto& v : vars) {
      std::string base = "rec" + type + "_" + v;
      // std::cout << " base is " << base << std::endl;
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

      // auto h_acceptance = dynamic_cast<TH1*>(h_acc->Clone((var + "_acceptance").c_str()));
      // h_acceptance->SetTitle(("Acceptance for " + var).c_str());
      // acceptHistos.push_back(std::shared_ptr<TH1>(h_acceptance));
    }
  }

  void GeneratePhiKinematicHistos(const std::string& type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi", "vz"};  // for Phi analysis
    // std::vector<std::string> vars = {"p", "theta", "phi"}; // for DVCS analysis
    for (const auto& v : vars) {
      std::string base = "rec" + type + "_" + v;
      // std::cout << " base is " << base << std::endl;
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

      // auto h_acceptance = dynamic_cast<TH1*>(h_acc->Clone((var + "_acceptance").c_str()));
      // h_acceptance->SetTitle(("Acceptance for " + var).c_str());
      // acceptHistos.push_back(std::shared_ptr<TH1>(h_acceptance));
    }
  }

  void GeneratePi0KinematicHistos(const std::string& type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi", "vz"};  // for Phi analysis
    // std::vector<std::string> vars = {"p", "theta", "phi"}; // for DVCS analysis
    for (const auto& v : vars) {
      std::string base = "rec" + type + "_" + v;
      std::string basename = "pi0_" + base;
      // std::cout << " base is " << base << std::endl;
      double histMin, histMax;
      auto it = kinematicAxisRanges.find(base);
      if (it != kinematicAxisRanges.end()) {
        histMin = it->second.first;
        histMax = it->second.second;
      } else {
        auto df_tmp = (v == "p") ? rdf_pi0_data->Filter(base + " > -100") : *rdf_pi0_data;
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
    auto result = kinCalc.ComputeDVCS_CrossSection_Weighted(rdf, bins, luminosity_nb_inv);
    if (dopi0corr) {
      auto sigma_pi0_3d = kinCalc.ComputeDVCS_CrossSection(*rdf_pi0_data, bins, luminosity_nb_inv);
      result = UsePi0Correction(result, sigma_pi0_3d, ComputePi0Corr(bins));
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

  std::vector<std::vector<std::vector<TH1D*>>> ComputePolDVCS_CrossSection(int pol, const BinManager& bins) {
    auto rdf_pol = rdf.Filter(Form("REC_Event_helicity == %d", pol));
    const auto n_all = static_cast<double>(rdf.Count().GetValue());
    const auto n_pol = static_cast<double>(rdf_pol.Count().GetValue());
    auto lumi_pol = luminosity_nb_inv * n_pol / n_all;
    std::cout << " Pol " << pol << " fraction: " << double(n_pol) / n_all << ", effective luminosity: " << lumi_pol << " nb^-1\n";
    
    auto result = kinCalc.ComputeDVCS_CrossSection_Weighted(rdf_pol, bins, lumi_pol);
    if (dopi0corr) {
      auto sigma_pi0_3d = kinCalc.ComputeDVCS_CrossSection(*rdf_pi0_data, bins, luminosity_nb_inv);
      result = UsePi0Correction(result, sigma_pi0_3d, ComputePi0Corr(bins));
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
  std::vector<std::vector<std::vector<TH1D*>>> ComputeBSA(const BinManager& bins, double pol = 1.0) {
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
          TH1D* hXS = xs3D[iq][it][ix];      // σ_uncorr
          TH1D* hCorr = corr3D[iq][it][ix];  //     c

          if (!hXS) {
            xsCorr3D[iq][it][ix] = nullptr;
            continue;
          }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError(b);
            const double c_val = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinContent(b) : 1.0;
            const double c_err = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinError(b) : 0.0;

            double val_corr = xs_val / (c_val);
            double err_corr = std::sqrt(std::pow(xs_err / c_val, 2) + std::pow(xs_val * c_err / (c_val * c_val), 2));

            if (c_val == 0.0 || c_err == 0.0) val_corr = xs_val;
            if (c_val == 0.0 || c_err == 0.0) err_corr = xs_err;
            // std::cout <<"xs_val " << xs_val << " xs_err " << xs_err << " c_val " << c_val << " c_err " << c_err << " val_corr " << val_corr << " err_corr " << err_corr <<
            // std::endl;

            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError(b, err_corr);
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
          TH1D* hXS = xs3D[iq][it][ix];      // σ_uncorr
          TH1D* hCorr = corr3D[iq][it][ix];  //     c

          if (!hXS) {
            xsCorr3D[iq][it][ix] = nullptr;
            continue;
          }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError(b);
            const double c_val = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinContent(b) : 1.0;
            const double c_err = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinError(b) : 0.0;
            double val_corr = xs_val / (c_val);
            double err_corr = std::sqrt(std::pow(xs_err / c_val, 2) + std::pow(xs_val * c_err / (c_val * c_val), 2));
            if (c_val == 0.0 || c_err == 0.0) val_corr = xs_val;
            if (c_val == 0.0 || c_err == 0.0) err_corr = xs_err;
            // std::cout <<"xs_val " << xs_val << " xs_err " << xs_err << " c_val " << c_val << " c_err " << c_err << " val_corr " << val_corr << " err_corr " << err_corr <<
            // std::endl;
            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError(b, err_corr);
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
          TH1D* hXS = xs3D[iq][it][ix];      // σ_uncorr
          TH1D* hCorr = corr3D[iq][it][ix];  //     c

          if (!hXS) {
            xsCorr3D[iq][it][ix] = nullptr;
            continue;
          }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError(b);
            const double c_val = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinContent(b) : 1.0;
            const double c_err = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinError(b) : 0.0;
            double val_corr = xs_val / (c_val);
            double err_corr = std::sqrt(std::pow(xs_err / c_val, 2) + std::pow(xs_val * c_err / (c_val * c_val), 2));
            if (c_val == 0.0 || c_err == 0.0) val_corr = xs_val;
            if (c_val == 0.0 || c_err == 0.0) err_corr = xs_err;
            if ((b == 1 || b == nb) && (c_val <= 0.8)) val_corr = -1;
            if ((b == 1 || b == nb) && (c_val <= 0.8)) err_corr = -1;
            // std::cout <<"xs_val " << xs_val << " xs_err " << xs_err << " c_val " << c_val << " c_err " << c_err << " val_corr " << val_corr << " err_corr " << err_corr <<
            // std::endl;
            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError(b, err_corr);
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

  std::vector<std::vector<std::vector<TH1D*>>> UseP1Cut(const std::vector<std::vector<std::vector<TH1D*>>>& xs3D, const std::vector<std::vector<std::vector<TH1D*>>>& corr3D) {
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
          TH1D* hXS = xs3D[iq][it][ix];      // σ_uncorr
          TH1D* hCorr = corr3D[iq][it][ix];  //     c

          if (!hXS) {
            xsCorr3D[iq][it][ix] = nullptr;
            continue;
          }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError(b);
            const double c_val = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinContent(b) : 1.0;
            const double c_err = (hCorr || hCorr->GetBinContent(b) == 0) ? hCorr->GetBinError(b) : 0.0;
            double val_corr = xs_val;
            double err_corr = xs_err;
            if (c_val <= 0.9) val_corr = -1;
            if (c_val <= 0.9) err_corr = -1;
            if (c_val <= 0.9) std::cout << " P1 cut applied: c_val " << c_val << "at bin " << b << std::endl;
            // std::cout <<"xs_val " << xs_val << " xs_err " << xs_err << " c_val " << c_val << " c_err " << c_err << " val_corr " << val_corr << " err_corr " << err_corr <<
            // std::endl;
            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError(b, err_corr);
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
          TH1D* hXS = xs3D[iq][it][ix];        // σ_uncorr
          TH1D* hPi0XS = xsPi03D[iq][it][ix];  // σ_pi0
          TH1D* hCorr = corr3D[iq][it][ix];    //     c

          if (!hXS) {
            xsCorr3D[iq][it][ix] = nullptr;
            continue;
          }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError(b);
            const double pi0_val = hPi0XS ? hPi0XS->GetBinContent(b) : 0.0;
            const double pi0_err = hPi0XS ? hPi0XS->GetBinError(b) : 0.0;
            const double c_val = hCorr ? hCorr->GetBinContent(b) : 0.0;
            const double c_err = hCorr ? hCorr->GetBinError(b) : 0.0;

            const double val_corr = xs_val * (1.0 - c_val);
            const double err_corr = std::sqrt(std::pow(xs_err * (1.0 - c_val), 2) + std::pow(xs_val * c_err, 2));

            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError(b, err_corr);
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
          TH1D* hXS = xs3D[iq][it][ix];        // σ_uncorr
          TH1D* hPi0XS = xsPi03D[iq][it][ix];  // σ_pi0
          TH1D* hCorr = corr3D[iq][it][ix];    //     c

          if (!hXS) {
            xsCorr3D[iq][it][ix] = nullptr;
            continue;
          }

          TH1D* hNew = dynamic_cast<TH1D*>(hXS->Clone(Form("%s_corr", hXS->GetName())));

          const int nb = hXS->GetNbinsX();
          for (int b = 1; b <= nb; ++b) {
            const double xs_val = hXS->GetBinContent(b);
            const double xs_err = hXS->GetBinError(b);
            const double pi0_val = hPi0XS ? hPi0XS->GetBinContent(b) : 0.0;
            const double pi0_err = hPi0XS ? hPi0XS->GetBinError(b) : 0.0;
            const double c_val = hCorr ? hCorr->GetBinContent(b) : 0.0;
            const double c_err = hCorr ? hCorr->GetBinError(b) : 0.0;

            const double val_corr = xs_val / (1.0 - c_val) - (pi0_val * c_val / (1.0 - c_val));
            const double err_corr = std::sqrt(std::pow(xs_err * (1.0 - c_val), 2) + std::pow(xs_val * c_err, 2));

            hNew->SetBinContent(b, val_corr);
            hNew->SetBinError(b, err_corr);
          }
          xsCorr3D[iq][it][ix] = hNew;
        }
      }
    }
    return xsCorr3D;
  }

  // std::vector<TH1D*> ComputePhiCrossSection(const BinManager& bins) { return kinCalc.ComputePhi_CrossSection(rdf, bins, luminosity_nb_inv); }

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

  inline std::vector<std::vector<std::vector<PhiMassDraw>>> MakePhiMassFitCanvases3D(const BinManager& bins, const std::string& outDirPerModel, int nMassBins = 200,
                                                                                     double mMin = 0.9874, double mMax = 1.120, bool constrainSigma = true, double sigmaRef = 0.004,
                                                                                     double sigmaFrac = 0.30,
                                                                                     // --- NEW (optional): compute dσ/dt when luminosity > 0 ---
                                                                                     double branching = 0.492) {
    if (luminosity_nb_inv <= 0.0) {
      std::cerr << "[MakePhiMassFitCanvases3D] luminosity<=0 → will NOT compute dσ/dt cache.\n";
    }
    // ensure directory exists
    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    const auto& q2Edges = bins.GetQ2Bins();
    const auto& tPrimeEdges = bins.GetTprimeBins();
    const auto& wEdges = bins.GetWBins();
    const bool hasW = !wEdges.empty();

    const size_t nQ = q2Edges.size() > 1 ? q2Edges.size() - 1 : 0;
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;
    const size_t nT = tPrimeEdges.size() > 1 ? tPrimeEdges.size() - 1 : 0;

    std::vector<std::vector<std::vector<PhiMassDraw>>> out;
    if (!nQ || !nT) return out;
    out.resize(nQ);

    // ---- NEW: prepare the 3D TH1 cache for dσ/dt (Q,W -> TH1 over t) ----
    phi_dsdt_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));
    phi_nsig_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));
    phi_mean_xB_QW_.assign(nQ, std::vector<std::vector<double>>(nW, std::vector<double>(nT, std::numeric_limits<double>::quiet_NaN())));
    phi_mean_W_QW_.assign (nQ, std::vector<std::vector<double>>(nW, std::vector<double>(nT, std::numeric_limits<double>::quiet_NaN())));
    phi_mean_GammaV_QW_.assign(nQ, std::vector<std::vector<double>>(nW, std::vector<double>(nT, std::numeric_limits<double>::quiet_NaN())));
    static std::atomic<unsigned long> uid_dsdt{0};

    for (size_t iq = 0; iq < nQ; ++iq) {
      out[iq].resize(nW);

      const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];
      auto df_q = rdf.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;
        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        out[iq][iw].resize(nT);

        // --- NEW: create the dσ/dt histogram for this (Q2,W) slice (if we will compute) ---
        if (luminosity_nb_inv > 0.0) {
          auto hname = Form("phi_dsdt_Q%zu_W%zu_%lu", iq, iw, uid_dsdt.fetch_add(1));
          auto htitle = Form("d#sigma/dt; -t [GeV^{2}]; nb/GeV^{2}");
          TH1D* hDsDt = new TH1D(hname, htitle, static_cast<int>(tPrimeEdges.size() - 1), tPrimeEdges.data());
          hDsDt->SetDirectory(nullptr);
          phi_dsdt_QW_[iq][iw] = hDsDt;

          auto hname_n = Form("phi_nsig_Q%zu_W%zu_%lu", iq, iw, uid_dsdt.load());
          TH1D* hNsig = new TH1D(hname_n, ";-t [GeV^{2}];N_{sig} (raw)", static_cast<int>(tPrimeEdges.size() - 1), tPrimeEdges.data());
          hNsig->SetDirectory(nullptr);
          hNsig->Sumw2();
          phi_nsig_QW_[iq][iw] = hNsig;
        }

        for (size_t it = 0; it < nT; ++it) {
          const double tLo = tPrimeEdges[it], tHi = tPrimeEdges[it + 1];
          const double dT  = tHi - tLo;
          const double dQ2 = qHi - qLo;
          // W is in GeV → W² in GeV²; use d(W²) so every dimension is GeV²
          // dσ/dt   [nb/GeV²]   uses binVol = dT   (Q2,W are bin labels)
          // dσ/dQ²dt [nb/GeV⁴] uses binVol = dQ2 * dT
          // d⁴σ/dQ²dW²dt [nb/GeV⁶] uses binVol = dQ2 * dW2 * dT
          const double dW2 = hasW ? (wHi * wHi - wLo * wLo) : 1.0;  // GeV²
          // Default: dσ/dt — divide only by dT; Q2 and W are bin labels
          double binVol = dT;

          // auto df_bin = df_qw.Filter([=](double t){ return t > tLo && t <= tHi; }, {"t"});
          auto df_bin = df_qw.Filter([=](double mtp) { return mtp > tLo && mtp <= tHi; }, {"mtprime"});
          static std::atomic<unsigned long> uid{0};
          auto hname = Form("hM_%zu_%zu_%zu_%lu", iq, iw, it, uid.fetch_add(1));
          auto hR = df_bin.Histo1D(ROOT::RDF::TH1DModel(hname, ";M_{K^{+}K^{-}} [GeV];Counts", 200, 0.8, 1.8), "invMass_KpKm");
          hR.GetValue();
          TH1D* h = (TH1D*)hR.GetPtr();
          if (!h || h->GetEntries() <= 20) {
            // still record an empty point in dσ/dt if we're computing it
            if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
              const int b = static_cast<int>(it + 1);
              phi_dsdt_QW_[iq][iw]->SetBinContent(b, 0.0);
              phi_dsdt_QW_[iq][iw]->SetBinError(b, 0.0);
            }
            continue;
          }

          const double mean_xB = *df_bin.Mean("xB");
          const double mean_W  = *df_bin.Mean("W");
          const double mean_Gv = *df_bin.Mean("Gamma_v");  // optional but recommended
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
          hDraw->GetXaxis()->SetRangeUser(mMin - 0.02, mMax + .020);

          const double bw = hDraw->GetXaxis()->GetBinWidth(1);
          std::cout << Form("[MakePhiMassFitCanvases3D] Fitting %s with %.0f bins, M=[%.4f,%.4f], BW=%.4f GeV", hDraw->GetName(), (double)hDraw->GetNbinsX(), mMin, mMax, bw)
                    << std::endl;

          std::string formula = Form(
              // Background: [0]*(x-0.9874)^[1] * exp(-[2]*(x-0.9874))
              "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))"
              " + "
              // Signal (counts/bin): [3] (=Nsig) * bw * Gaussian_pdf(x; [4]=mean, [5]=sigma)
              "[3]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[4])/([5]),2))/TMath::Abs([5])",
              bw);

          TF1* fTot = new TF1(Form("fitTotal_%s", hname), formula.c_str(), mMin, mMax);

          //[0]*0.398942*bw*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])
          // params: [0]=Abkg(scale), [1]=alpha, [2]=lambda, [3]=Nsig(yield), [4]=mu, [5]=sigma
          fTot->SetParNames("A", "alpha", "lambda", "N", "mu", "sigma");
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

          const double mLo = mu - 3.0 * sigma;
          const double mHi = mu + 3.0 * sigma;

          TF1* fSig = new TF1(Form("fSig_%s", hname), Form("[0]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[1])/([2]),2))/TMath::Abs([2])", bw), mMin, mMax);
          fSig->SetParameters(N, mu, sigma);
          fSig->SetLineColor(kOrange + 1);
          fSig->SetLineStyle(2);
          fSig->SetLineWidth(2);
          fSig->SetFillColorAlpha(kOrange - 3, 0.3);
          fSig->SetFillStyle(1001);
          fSig->Draw("SAME FC");

          TF1* fBkg = new TF1(Form("fBkg_%s", hname), "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))", mMin, mMax);
          fBkg->SetParameters(A, alpha, lambda);
          fBkg->SetLineColor(kGreen + 2);
          fBkg->SetLineStyle(3);
          fBkg->SetLineWidth(2);
          fBkg->SetFillColorAlpha(kGreen - 7, 0.3);
          fBkg->SetFillStyle(1001);
          fBkg->Draw("SAME FC");

          double Nbkg = 0.0;
          const double Nsig = fTot->GetParameter(3);
          const double Nsig_err = fTot->GetParError(3);

          // window counts & signal
          auto dfWin = df_bin.Filter([=](float m) { return m > mLo && m < mHi; }, {"invMass_KpKm"});
          const ULong64_t Nwin = *dfWin.Count();

          if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
            double Aeps = 1.0;
          
          if (doacceptcorr && iq < phi_accept_QW_.size()
              && phi_accept_QW_[iq][iw] != nullptr) {
              const double a = phi_accept_QW_[iq][iw]->GetBinContent(int(it+1));
              if (a > 0.0) Aeps = a;
          }
            const int b = static_cast<int>(it + 1);
            double val = 0.0, err = 0.0;
            val = Nsig / (luminosity_nb_inv * Aeps * branching * binVol);
            err = Nsig_err / (luminosity_nb_inv * Aeps * branching * binVol);
            

            if (doradiativecorr && iq < phi_rad_QW_.size() && iw < phi_rad_QW_[iq].size() && phi_rad_QW_[iq][iw]) {
              const int b = int(it + 1);
              const double Crad = phi_rad_QW_[iq][iw]->GetBinContent(b);
              if (Crad > 0.0) {
                val /= Crad;
                err /= Crad;
              }
            }
            phi_dsdt_QW_[iq][iw]->SetBinContent(b, val);
            phi_dsdt_QW_[iq][iw]->SetBinError(b, err);

            // Store raw fit yield (Nsig, Nsig_err) regardless of corrections
            if (phi_nsig_QW_[iq][iw]) {
              phi_nsig_QW_[iq][iw]->SetBinContent(b, Nsig);
              phi_nsig_QW_[iq][iw]->SetBinError(b, Nsig_err);
            }
          }
          phi_mean_xB_QW_[iq][iw][it]     = mean_xB;
          phi_mean_W_QW_[iq][iw][it]      = mean_W;
          phi_mean_GammaV_QW_[iq][iw][it] = mean_Gv;
          // ---- (existing) canvas beautification & save ----
          TCanvas* c = new TCanvas(Form("c_%s", hname), "K^{+}K^{-} mass", 1200, 1000);
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
          L1->SetLineStyle(2);
          L2->SetLineStyle(2);
          L1->SetLineWidth(2);
          L2->SetLineWidth(2);
          L1->Draw("SAME");
          L2->Draw("SAME");

          TLegend* leg = new TLegend(0.55, 0.5, 0.85, 0.75);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->AddEntry(hDraw, "K^{+}K^{-} Inv. Mass", "lep");
          leg->AddEntry(fTot, "Total Fit: Exp + Gauss", "l");
          leg->AddEntry(fSig, "Signal (Gauss)", "f");
          leg->AddEntry(fBkg, "Background (Exp)", "f");
          leg->AddEntry(L1, "#pm 3#sigma", "l");
          leg->Draw();

          TLatex latex;
          latex.SetNDC();
          latex.SetTextSize(0.035);
          latex.DrawLatex(0.45, 0.89, Form("Q^{2}[%.2f,%.2f]  %s  -t'[%.2f,%.2f]", qLo, qHi, hasW ? Form("W[%.1f,%.1f]", wLo, wHi) : "", tLo, tHi));
          latex.DrawLatex(0.55, 0.85, Form("#mu = %.3f GeV, #sigma = %.3f GeV", mu, sigma));
          const double chi2ndf = fTot->GetNDF() > 0 ? fTot->GetChisquare() / fTot->GetNDF() : 0.0;
          latex.DrawLatex(0.55, 0.81, Form("#chi^{2}/ndf = %.2f", chi2ndf));
          latex.DrawLatex(0.55, 0.77, Form("N_{#phi} (total) = %.1f #pm %.1f", Nsig, Nsig_err));

          auto tagRaw = hasW ? Form("Q2_%.1f_%.1f__W_%.1f_%.2f__tprime_%.1f_%.1f", qLo, qHi, wLo, wHi, tLo, tHi) : Form("Q2_%.1f_%.1f__tprime_%.1f_%.1f", qLo, qHi, tLo, tHi);
          std::string tag = tagRaw;
          std::replace(tag.begin(), tag.end(), '.', '_');
          c->SaveAs(Form("%s/KKmass_%s.pdf", outDirPerModel.c_str(), tag.c_str()));

          // pack draw products
          PhiMassDraw pack;
          pack.h = hDraw;
          pack.fTot = fTot;
          pack.fSig = fSig;
          pack.fBkg = fBkg;
          pack.mu = mu;
          pack.sigma = sigma;
          pack.mLo = mLo;
          pack.mHi = mHi;
          pack.c = c;
          pack.name = tag;
          out[iq][iw][it] = pack;

          // clean transient objects we don't cache (canvas snapshot saved already)
          delete c;
          // delete hDraw;
          if (gPad) {
            gPad->Clear();
            gPad->Update();
            gPad->SetSelected(nullptr);
          }
          delete leg;
          delete L1;
          delete L2;
        }  // it (t bins)
      }  // iw (W bins)
    }  // iq (Q2 bins)
    return out;
  }

  inline std::vector<std::vector<TH1D*>> MakePhiAcceptanceCorrection3D(const BinManager& bins, const std::string& outDirPerModel, const std::string& genTPrimeVar = "mtprime",
                                                                       const std::string& recTPrimeVar = "mtprime", int minGenEntries = 10) {
    // ── Guard: both MC RDFs must be present ─────────────────────────────────
    if (!rdf_gen_phimc || !rdf_accept_phimc) {
      std::cerr << "[MakePhiAcceptanceCorrection3D] "
                   "rdf_gen_phimc or rdf_accept_phimc is not set. "
                   "Call SetPhiMCDataFrames() before this function.\n";
      return {};
    }

    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    // ── Bin edges ────────────────────────────────────────────────────────────
    const auto& q2Edges = bins.GetQ2Bins();
    const auto& tPrimeEdges = bins.GetTprimeBins();
    const auto& wEdges = bins.GetWBins();
    const bool hasW = !wEdges.empty();

    const size_t nQ = q2Edges.size() > 1 ? q2Edges.size() - 1 : 0;
    const size_t nW = hasW ? wEdges.size() - 1 : 1;
    const size_t nT = tPrimeEdges.size() > 1 ? tPrimeEdges.size() - 1 : 0;

    if (!nQ || !nT) {
      std::cerr << "[MakePhiAcceptanceCorrection3D] Empty Q2 or t' bin list.\n";
      return {};
    }

    // ── Allocate output and cache ────────────────────────────────────────────
    std::vector<std::vector<TH1D*>> out(nQ, std::vector<TH1D*>(nW, nullptr));
    phi_accept_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));

    static std::atomic<unsigned long> uid_acc{0};

    // ── Loop over (Q2, W) slices ─────────────────────────────────────────────
    for (size_t iq = 0; iq < nQ; ++iq) {
      const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];

      // Filter the MC RDFs to this Q2 slice
      auto df_gen_q = rdf_gen_phimc->Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});
      auto df_rec_q = rdf_accept_phimc->Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;

        // Filter the MC RDFs to this W slice
        auto df_gen_qw = df_gen_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});
        auto df_rec_qw = df_rec_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        // ── Build the acceptance histogram for this (Q2, W) cell ─────────────
        const auto hname_acc = Form("phi_accept_Q%zu_W%zu_%lu", iq, iw, uid_acc.fetch_add(1));
        const auto htitle_acc = hasW ? Form("Acceptance; -t' [GeV^{2}]; A(t')") : Form("Acceptance; -t' [GeV^{2}]; A(t')");

        TH1D* hAcc = new TH1D(hname_acc, htitle_acc, static_cast<int>(tPrimeEdges.size() - 1), tPrimeEdges.data());
        hAcc->SetDirectory(nullptr);
        hAcc->Sumw2();

        // ── N_gen: raw count histogram of generated MC events vs t' ────────────
        // Gen events are pure signal (no background), so raw counts are correct.
        const auto hname_gen = Form("hGen_Q%zu_W%zu_%lu", iq, iw, uid_acc.load());
        auto hGenR = df_gen_qw.Histo1D(
            ROOT::RDF::TH1DModel(hname_gen, ";-t' [GeV^{2}];N_{gen}",
                                 static_cast<int>(tPrimeEdges.size() - 1), tPrimeEdges.data()),
            genTPrimeVar);
        TH1D* hGen = (TH1D*)hGenR.GetPtr();

        // ── N_rec: phi signal yield from K+K- mass fit on reconstructed MC ────
        // We MUST fit (not just count) because the reconstructed MC sample
        // contains combinatorial background exactly like the data.  Using raw
        // event counts would overestimate N_rec and therefore overestimate A,
        // which would make the cross-section systematically too low.
        // The fit model is identical to MakePhiMassFitCanvases3D so that the
        // acceptance and the data yield are extracted in exactly the same way.

        // Build a per-t'-bin array of fitted yields (same loop structure as data fit)
        std::vector<double> nRecFit(nT, 0.0);
        std::vector<double> nRecFitErr(nT, 0.0);

        // Also keep a raw-count hRec for the diagnostic canvas (display only)
        const auto hname_rec = Form("hRec_Q%zu_W%zu_%lu", iq, iw, uid_acc.load() + 1);
        auto hRecR = df_rec_qw.Histo1D(
            ROOT::RDF::TH1DModel(hname_rec, ";-t' [GeV^{2}];N_{rec} (raw)",
                                 static_cast<int>(tPrimeEdges.size() - 1), tPrimeEdges.data()),
            recTPrimeVar);
        TH1D* hRec = (TH1D*)hRecR.GetPtr();  // raw counts — used only for display

        for (size_t it = 0; it < nT; ++it) {
          const double tLo = tPrimeEdges[it], tHi = tPrimeEdges[it + 1];

          // Filter reconstructed MC to this t' bin
          auto df_rec_bin = df_rec_qw.Filter(
              [=](double mtp) { return mtp > tLo && mtp <= tHi; }, {"mtprime"});

          static std::atomic<unsigned long> uid_mfit{0};
          const auto hmname = Form("hMrec_acc_Q%zu_W%zu_T%zu_%lu", iq, iw, it, uid_mfit.fetch_add(1));

          // Build K+K- invariant mass histogram for MC rec in this bin
          auto hMrecR = df_rec_bin.Histo1D(
              ROOT::RDF::TH1DModel(hmname, ";M_{K^{+}K^{-}} [GeV];Counts", 200, 0.8, 1.8),
              "invMass_KpKm");
          hMrecR.GetValue();
          TH1D* hMrec = (TH1D*)hMrecR.GetPtr();

          if (!hMrec || hMrec->GetEntries() < 5) {
            std::cout << Form(
                "[MakePhiAcceptanceCorrection3D] "
                "MC rec: too few events in Q2[%.2f,%.2f] %s t'[%.3f,%.3f] "
                "(%.0f entries) — N_rec set to 0.\n",
                qLo, qHi, hasW ? Form("W[%.2f,%.2f]", wLo, wHi) : "",
                tLo, tHi, hMrec ? (double)hMrec->GetEntries() : 0.0);
            continue;
          }

          // ── Same fit model as MakePhiMassFitCanvases3D ──────────────────────
          // Background: A*(x-0.9874)^alpha * exp(-lambda*(x-0.9874))
          // Signal:     Nsig * bw * Gauss(x; mu, sigma) / |sigma|
          const double mMin_fit = 0.9874, mMax_fit = 1.120;
          const double bw = hMrec->GetXaxis()->GetBinWidth(1);

          std::string formula = Form(
              "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))"
              " + "
              "[3]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[4])/([5]),2))/TMath::Abs([5])",
              bw);

          TF1* fMC = new TF1(Form("fMCrec_%s", hmname), formula.c_str(), mMin_fit, mMax_fit);
          fMC->SetParNames("A", "alpha", "lambda", "Nsig", "mu", "sigma");
          fMC->SetParameters(4, 0.9, 2, 10, 1.019, 0.004);
          fMC->SetParLimits(4, 1.005, 1.022);   // mu
          fMC->SetParLimits(5, 0.0025, 0.025);  // sigma
          fMC->SetParLimits(3, 0.0, 1e7);        // Nsig >= 0
          fMC->SetLineColor(kRed + 1);
          fMC->SetLineWidth(2);

          // Save the mass-fit canvas for MC rec alongside the acceptance canvas
          TH1D* hMrecDraw = (TH1D*)hMrec->Clone(Form("%s_draw", hmname));
          hMrecDraw->SetDirectory(nullptr);
          hMrecDraw->Fit(fMC, "R0QL");

          const double Nsig_rec    = fMC->GetParameter(3);
          const double Nsig_rec_err= fMC->GetParError(3);

          // Save mass-fit canvas for MC rec (diagnostic)
          {
            auto tagRaw = hasW
                ? Form("MCrecMassFit_Q2_%.2f_%.2f__W_%.2f_%.2f__T_%.3f_%.3f", qLo, qHi, wLo, wHi, tLo, tHi)
                : Form("MCrecMassFit_Q2_%.2f_%.2f__T_%.3f_%.3f",              qLo, qHi,             tLo, tHi);
            std::string tagStr = tagRaw;
            std::replace(tagStr.begin(), tagStr.end(), '.', '_');

            std::unique_ptr<TCanvas> cMfit(new TCanvas(
                Form("c_mcrec_%s", tagStr.c_str()), "MC Rec K+K- mass", 900, 700));
            gStyle->SetOptStat(0);
            gStyle->SetOptFit(0);
            hMrecDraw->SetLineColor(kBlue + 1);
            hMrecDraw->SetLineWidth(2);
            hMrecDraw->SetMarkerStyle(20);
            hMrecDraw->GetXaxis()->SetTitle("M(K^{+}K^{-}) [GeV]");
            hMrecDraw->GetYaxis()->SetTitle("Counts");
            hMrecDraw->GetXaxis()->SetRangeUser(mMin_fit - 0.02, mMax_fit + 0.02);
            hMrecDraw->Draw("PE");
            fMC->SetLineColor(kRed + 1);
            fMC->Draw("SAME C");

            TLatex lx;
            lx.SetNDC(); lx.SetTextSize(0.035);
            lx.DrawLatex(0.55, 0.85, Form("N_{#phi}^{rec} = %.1f #pm %.1f", Nsig_rec, Nsig_rec_err));
            lx.DrawLatex(0.55, 0.81, Form("#mu = %.4f GeV,  #sigma = %.4f GeV",
                                          fMC->GetParameter(4), fMC->GetParameter(5)));

            cMfit->SaveAs(Form("%s/%s.pdf", outDirPerModel.c_str(), tagStr.c_str()));
          }

          delete hMrecDraw;
          delete fMC;

          if (Nsig_rec > 0.0) {
            nRecFit[it]    = Nsig_rec;
            nRecFitErr[it] = Nsig_rec_err;
          }
        }

        // ── Per-bin acceptance calculation using fitted N_rec ─────────────────
        for (size_t it = 0; it < nT; ++it) {
          const int b = static_cast<int>(it + 1);  // ROOT bin index (1-based)

          const double nGen    = hGen ? hGen->GetBinContent(b) : 0.0;
          const double nRec    = nRecFit[it];     // ← fitted phi yield from MC rec
          const double nRecErr = nRecFitErr[it];

          double acc    = 0.0;
          double accErr = 0.0;

          if (nGen >= minGenEntries && nRec > 0.0) {
            acc = nRec / nGen;
            // Combined uncertainty: σ_A = A * sqrt( (σ_Nrec/Nrec)^2 + (1/Ngen) )
            // (fit uncertainty on Nrec + Poisson on Ngen)
            accErr = acc * std::sqrt(
                (nRecErr > 0.0 ? (nRecErr / nRec) * (nRecErr / nRec) : 0.0)
                + 1.0 / nGen);
          } else if (nGen < minGenEntries && nGen > 0) {
            std::cout << Form(
                "[MakePhiAcceptanceCorrection3D] "
                "Low statistics in Q2[%.2f,%.2f] %s t'[%.3f,%.3f]: "
                "N_gen=%.0f < %d — setting A=0.\n",
                qLo, qHi, hasW ? Form("W[%.2f,%.2f]", wLo, wHi) : "",
                tPrimeEdges[it], tPrimeEdges[it + 1], nGen, minGenEntries);
          }

          hAcc->SetBinContent(b, acc);
          hAcc->SetBinError(b, accErr);
        }

        // ── Diagnostic canvas (mirrors MakePhiMassFitCanvases3D style) ───────
        {
          const auto cname = Form("c_acc_Q%zu_W%zu_%lu", iq, iw, uid_acc.fetch_add(1));
          std::unique_ptr<TCanvas> c(new TCanvas(cname, "Acceptance vs t'", 1200, 1000));
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(0);
          c->SetTicks(1, 1);
          c->SetTopMargin(0.06);
          c->SetRightMargin(0.03);
          c->SetBottomMargin(0.11);
          c->SetLeftMargin(0.13);

          // ── Draw N_gen and N_rec on the same pad (left axis = counts) ──────
          TPad* padTop = new TPad("padTop", "", 0.0, 0.35, 1.0, 1.0);
          padTop->SetBottomMargin(0.02);
          padTop->SetLeftMargin(0.13);
          padTop->SetRightMargin(0.03);
          padTop->Draw();

          TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.35);
          padBot->SetTopMargin(0.02);
          padBot->SetBottomMargin(0.24);
          padBot->SetLeftMargin(0.13);
          padBot->SetRightMargin(0.03);
          padBot->Draw();

          padTop->cd();
          TH1D* hGenDraw = (TH1D*)hGen->Clone(Form("%s_genDraw", hname_gen));
          TH1D* hRecDraw = (TH1D*)hRec->Clone(Form("%s_recDraw", hname_rec));
          hGenDraw->SetDirectory(nullptr);
          hRecDraw->SetDirectory(nullptr);
          hGenDraw->SetLineColor(kBlue + 1);
          hGenDraw->SetFillColorAlpha(kBlue - 9, 0.3);
          hGenDraw->SetLineWidth(2);
          hGenDraw->GetXaxis()->SetLabelSize(0);
          hGenDraw->GetYaxis()->SetTitle("Events");
          hGenDraw->SetTitle("");
          hRecDraw->SetLineColor(kRed + 1);
          hRecDraw->SetFillColorAlpha(kRed - 9, 0.3);
          hRecDraw->SetLineWidth(2);
          hRecDraw->SetTitle("");

          const double yMax = std::max(hGenDraw->GetMaximum(), hRecDraw->GetMaximum()) * 1.25;
          hGenDraw->GetYaxis()->SetRangeUser(0.0, yMax);
          hGenDraw->Draw("HIST");
          hRecDraw->Draw("HIST SAME");

          TLegend legTop(0.55, 0.72, 0.90, 0.88);
          legTop.SetBorderSize(0);
          legTop.SetFillStyle(0);
          legTop.SetTextSize(0.045);
          legTop.AddEntry(hGenDraw, "Generated", "f");
          legTop.AddEntry(hRecDraw, "Reconstructed (raw, display only)", "f");
          legTop.Draw();

          // Bin label (Q2, W)
          TLatex latTop;
          latTop.SetNDC();
          latTop.SetTextFont(42);
          latTop.SetTextSize(0.048);
          latTop.DrawLatex(0.15, 0.92, hasW ? Form("Q^{2}[%.2f,%.2f]  W[%.2f,%.2f]", qLo, qHi, wLo, wHi) : Form("Q^{2}[%.2f,%.2f]", qLo, qHi));

          padBot->cd();
          TH1D* hAccDraw = (TH1D*)hAcc->Clone(Form("%s_draw", hname_acc));
          hAccDraw->SetDirectory(nullptr);
          hAccDraw->SetTitle("");
          hAccDraw->SetLineColor(kBlack);
          hAccDraw->SetMarkerStyle(20);
          hAccDraw->SetMarkerSize(1.1);
          hAccDraw->SetLineWidth(2);
          hAccDraw->GetXaxis()->SetTitle("-t' [GeV^{2}]");
          hAccDraw->GetYaxis()->SetTitle("A(t')");
          hAccDraw->GetYaxis()->SetRangeUser(0.0, 1.05);
          hAccDraw->GetXaxis()->SetTitleSize(0.10);
          hAccDraw->GetXaxis()->SetLabelSize(0.09);
          hAccDraw->GetYaxis()->SetTitleSize(0.09);
          hAccDraw->GetYaxis()->SetLabelSize(0.09);
          hAccDraw->GetYaxis()->SetTitleOffset(0.65);
          hAccDraw->Draw("E1");

          // Draw horizontal lines at 0 and 1 for reference
          TLine lo(tPrimeEdges.front(), 0.0, tPrimeEdges.back(), 0.0);
          TLine hi(tPrimeEdges.front(), 1.0, tPrimeEdges.back(), 1.0);
          lo.SetLineStyle(2);
          lo.SetLineColor(kGray + 1);
          hi.SetLineStyle(2);
          hi.SetLineColor(kGray + 1);
          lo.Draw("SAME");
          hi.Draw("SAME");

          // Build the file-safe tag
          auto tagRaw = hasW ? Form("Q2_%.2f_%.2f__W_%.2f_%.2f", qLo, qHi, wLo, wHi) : Form("Q2_%.2f_%.2f", qLo, qHi);
          std::string tag = tagRaw;
          std::replace(tag.begin(), tag.end(), '.', '_');

          c->SaveAs(Form("%s/Acceptance_%s.pdf", outDirPerModel.c_str(), tag.c_str()));

          delete hGenDraw;
          delete hRecDraw;
          delete hAccDraw;
        }  // end diagnostic canvas scope

        // ── Store in both output vector and the persistent cache ─────────────
        out[iq][iw] = hAcc;
        phi_accept_QW_[iq][iw] = hAcc;

      }  // iw
    }  // iq

    return out;
  }

  inline std::vector<std::vector<TH1D*>> MakePhiRadiativeCorrection3D_FromRatio(
    const BinManager& bins,
    const std::string& outDirPerModel,
    const std::string& ratioVar = "rad_corr",   // column in RDF with radiative ratio
    const std::string& tPrimeVar = "mtprime",   // same t' variable you use elsewhere
    const std::string& weightVar = "",          // optional weight column
    int minEntries = 10)
{
  if (!rdf_phimc_radRatio) {
    std::cerr << "[MakePhiRadiativeCorrection3D_FromRatio] rdf_phimc_radRatio not set.\n";
    return {};
  }

  gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

  const auto& q2Edges     = bins.GetQ2Bins();
  const auto& wEdges      = bins.GetWBins();
  const auto& tPrimeEdges = bins.GetTprimeBins();

  const bool hasW = !wEdges.empty();
  const size_t nQ = (q2Edges.size() > 1) ? (q2Edges.size() - 1) : 0;
  const size_t nW = hasW ? (wEdges.size() - 1) : 1;

  if (!nQ || tPrimeEdges.size() < 2) return {};

  phi_rad_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));

  // Build per-(Q,W) histograms of the mean ratio in each t' bin
  for (size_t iq = 0; iq < nQ; ++iq) {
    const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];

    auto df_q = rdf_phimc_radRatio->Filter(
        [=](double Q2) { return (Q2 > qLo && Q2 <= qHi); }, {"Q2"});

    for (size_t iw = 0; iw < nW; ++iw) {
      const double wLo = hasW ? wEdges[iw] : 0.0;
      const double wHi = hasW ? wEdges[iw + 1] : 1e9;

      auto df_qw = df_q.Filter(
          [=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

      // IMPORTANT: keep the RResultPtr alive until after we read the profile.
      ROOT::RDF::RResultPtr<TProfile> pR;

      if (!weightVar.empty()) {
        pR = df_qw.Profile1D(
            ROOT::RDF::TProfile1DModel(
                Form("pRad_Q%zu_W%zu", iq, iw),
                ";-t' [GeV^{2}];C_{rad}",
                (int)(tPrimeEdges.size() - 1), tPrimeEdges.data()),
            tPrimeVar, ratioVar, weightVar);
      } else {
        pR = df_qw.Profile1D(
            ROOT::RDF::TProfile1DModel(
                Form("pRad_Q%zu_W%zu", iq, iw),
                ";-t' [GeV^{2}];C_{rad}",
                (int)(tPrimeEdges.size() - 1), tPrimeEdges.data()),
            tPrimeVar, ratioVar);
      }

      // Force evaluation so profile is filled before we read bin contents/entries.
      pR.GetValue();
      TProfile* p = pR.GetPtr();

      auto h = new TH1D(
          Form("hRadCorr_Q%zu_W%zu", iq, iw),
          ";-t' [GeV^{2}];C_{rad}",
          (int)(tPrimeEdges.size() - 1), tPrimeEdges.data());
      h->SetDirectory(nullptr);
      h->Sumw2();

      // Convert profile bin means/errors into TH1
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double n = (p ? p->GetBinEntries(b) : 0.0);

        if (p && n >= minEntries && std::isfinite(p->GetBinContent(b))) {
          std::cout << Form(
              "[MakePhiRadiativeCorrection3D_FromRatio] Q[%zu] W[%zu] t' bin %d: n=%.0f entries, C_rad=%.6g ± %.6g",
              iq, iw, b, n, p->GetBinContent(b), p->GetBinError(b))
                    << std::endl;

          h->SetBinContent(b, p->GetBinContent(b));
          h->SetBinError(b, p->GetBinError(b));
        } else {
          // fallback: no correction
          h->SetBinContent(b, 1.0);
          h->SetBinError(b, 0.0);
        }
      }

      phi_rad_QW_[iq][iw] = h;

      // Optional: write out per-(Q,W) correction hist
      {
        TFile fOut((outDirPerModel + "/phi_radCorr_Q" + std::to_string(iq) + "_W" + std::to_string(iw) + ".root").c_str(), "RECREATE");
        h->Write();
        if (p) p->Write(); // useful for debugging
        fOut.Close();
      }
    }
  }

  return phi_rad_QW_;
}


  struct YieldRes {
    double N{0.0}, dN{0.0};
    double mu{0.0}, sigma{0.0};
    double mLo{0.0}, mHi{0.0};
  };

  inline YieldRes FitPhiMassYieldAndSave(ROOT::RDF::RNode df_in, const std::string& outDir, const std::string& tagRaw,
                                         const std::string& helTag,  // "P" or "M"
                                         int nMassBins, double mMin, double mMax, bool constrainSigma, double sigmaRef, double sigmaFrac) {
    static std::atomic<unsigned long> uid{0};
    const auto hname = Form("hM_BSA_%s_%s_%lu", tagRaw.c_str(), helTag.c_str(), uid.fetch_add(1));

    auto hR = df_in.Histo1D(ROOT::RDF::TH1DModel(hname, ";M_{K^{+}K^{-}} [GeV];Counts", nMassBins, mMin, mMax), "invMass_KpKm");
    hR.GetValue();
    TH1D* h0 = (TH1D*)hR.GetPtr();
    if (!h0 || h0->GetEntries() < 20) return {};

    std::unique_ptr<TH1D> hDraw((TH1D*)h0->Clone(Form("%s_draw", hname)));
    hDraw->SetDirectory(nullptr);
    hDraw->SetTitle("");
    hDraw->SetLineColor(kBlue + 1);
    hDraw->SetLineWidth(2);
    hDraw->SetMarkerStyle(20);
    hDraw->SetMarkerSize(1.1);
    hDraw->GetXaxis()->SetTitle("M(K^{+}K^{-}) [GeV]");
    hDraw->GetYaxis()->SetTitle("Counts");

    const double bw = hDraw->GetXaxis()->GetBinWidth(1);

    // bkg + (Nsig * bw * Gaussian_pdf)
    std::string formula = Form(
        "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))"
        " + "
        "[3]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[4])/([5]),2))/TMath::Abs([5])",
        bw);

    std::unique_ptr<TF1> fTot(new TF1(Form("fTot_%s", hname), formula.c_str(), mMin, mMax));
    fTot->SetParNames("A", "alpha", "lambda", "Nsig", "mu", "sigma");

    const double peakX = 1.019;
    const double peakY = hDraw->GetMaximum();
    fTot->SetParameters(std::max(1.0, peakY), 0.9, 2.0, std::max(1.0, peakY * bw * 10.0), peakX, sigmaRef);

    fTot->SetParLimits(4, 1.005, 1.022);
    if (constrainSigma)
      fTot->SetParLimits(5, sigmaRef * (1.0 - sigmaFrac), sigmaRef * (1.0 + sigmaFrac));
    else
      fTot->SetParLimits(5, 0.0025, 0.025);
    fTot->SetParLimits(3, 0.0, 1e9);

    // Fit (no draw)
    hDraw->Fit(fTot.get(), "R0Q");

    const double A = fTot->GetParameter(0);
    const double alpha = fTot->GetParameter(1);
    const double lambda = fTot->GetParameter(2);
    const double Nsig = fTot->GetParameter(3);
    const double dNsig = fTot->GetParError(3);
    const double mu = fTot->GetParameter(4);
    const double sigma = fTot->GetParameter(5);

    const double mLo = mu - 3.0 * sigma;
    const double mHi = mu + 3.0 * sigma;

    std::unique_ptr<TF1> fSig(new TF1(Form("fSig_%s", hname), Form("[0]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[1])/([2]),2))/TMath::Abs([2])", bw), mMin, mMax));
    fSig->SetParameters(Nsig, mu, sigma);
    fSig->SetLineColor(kOrange + 1);
    fSig->SetLineStyle(2);
    fSig->SetLineWidth(2);
    fSig->SetFillColorAlpha(kOrange - 3, 0.30);
    fSig->SetFillStyle(1001);

    std::unique_ptr<TF1> fBkg(new TF1(Form("fBkg_%s", hname), "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))", mMin, mMax));
    fBkg->SetParameters(A, alpha, lambda);
    fBkg->SetLineColor(kGreen + 2);
    fBkg->SetLineStyle(3);
    fBkg->SetLineWidth(2);
    fBkg->SetFillColorAlpha(kGreen - 7, 0.30);
    fBkg->SetFillStyle(1001);

    // Draw & save
    std::unique_ptr<TCanvas> c(new TCanvas(Form("c_%s", hname), "K^{+}K^{-} mass", 1200, 1000));
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    c->SetTicks(1, 1);
    c->SetTopMargin(0.04);
    c->SetRightMargin(0.03);
    c->SetBottomMargin(0.11);
    c->SetLeftMargin(0.13);

    hDraw->SetMinimum(0.0);
    hDraw->Draw("PE");
    fBkg->Draw("SAME FC");
    fSig->Draw("SAME FC");
    fTot->SetLineColor(kRed + 1);
    fTot->SetLineWidth(3);
    fTot->Draw("SAME C");

    const double yMax = hDraw->GetMaximum();
    TLine L1(mLo, 0.0, mLo, yMax * 0.75);
    TLine L2(mHi, 0.0, mHi, yMax * 0.75);
    L1.SetLineColor(kMagenta + 2);
    L2.SetLineColor(kMagenta + 2);
    L1.SetLineStyle(2);
    L2.SetLineStyle(2);
    L1.SetLineWidth(2);
    L2.SetLineWidth(2);
    L1.Draw("SAME");
    L2.Draw("SAME");

    TLegend leg(0.55, 0.5, 0.85, 0.75);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hDraw.get(), "K^{+}K^{-} Inv. Mass", "lep");
    leg.AddEntry(fTot.get(), "Total Fit", "l");
    leg.AddEntry(fSig.get(), "Signal", "f");
    leg.AddEntry(fBkg.get(), "Background", "f");
    leg.AddEntry(&L1, "#pm 3#sigma", "l");
    leg.Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.55, 0.85, Form("#mu = %.3f GeV, #sigma = %.3f GeV", mu, sigma));
    latex.DrawLatex(0.55, 0.81, Form("N_{#phi} = %.1f #pm %.1f", Nsig, dNsig));

    std::string tag = tagRaw;
    std::replace(tag.begin(), tag.end(), '.', '_');
    c->SaveAs(Form("%s/KKmass_%s_hel%s.pdf", outDir.c_str(), tag.c_str(), helTag.c_str()));

    return {Nsig, dNsig, mu, sigma, mLo, mHi};
  }

  // ---------------- NEW: build A_LU(cos#theta_{KK}) cache and save helicity-separated mass-fit canvases ---------------
  inline void MakePhiBSAMassFitCanvases3D(const BinManager& bins, const std::string& outDirPerModel, int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                          bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    const auto& q2Edges = bins.GetQ2Bins();
    const auto& tPrimeEdges = bins.GetTprimeBins();
    const auto& wEdges = bins.GetWBins();
    const bool hasW = !wEdges.empty();
    const auto& cEdges = bins.GetCosThetaKKBins();

    const size_t nQ = (q2Edges.size() > 1) ? (q2Edges.size() - 1) : 0;
    const size_t nT = (tPrimeEdges.size() > 1) ? (tPrimeEdges.size() - 1) : 0;
    // If user provided 1..-1, reverse it
    // if (cEdges.size() >= 2 && cEdges.front() > cEdges.back()) std::reverse(cEdges.begin(), cEdges.end());

    // Verify strictly increasing
    for (size_t i = 1; i < cEdges.size(); ++i) {
      if (!(cEdges[i] > cEdges[i - 1])) {
        std::cout << "Bin edges:\n";
        for (const auto& e : cEdges) std::cout << e << " ";
        std::cout << std::endl;
        std::cerr << "[MakePhiBSAMassFitCanvases3D] cos(theta) bin edges not strictly increasing!\n";
        return;
      }
    }
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;
    const size_t nC = cEdges.size() - 1;
    if (!nQ || !nC) {
      std::cerr << "[MakePhiALUCosThetaFitCanvases3D] missing Q2 or cos(theta_KK) binning.\n";
      return;
    }

    // allocate cache: (Q,W)->TH1(t')
    phi_alu_cos_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));
    static std::atomic<unsigned long> uidH{0};

    for (size_t iq = 0; iq < nQ; ++iq) {
      const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];
      auto df_q = rdf.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;
        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        // cache hist for this (Q,W)
        {
          auto hname = Form("phi_alu_t_Q%zu_W%zu_%lu", iq, iw, uidH.fetch_add(1));
          auto htitle = "A_{LU}(cos#theta_{KK});cos#theta_{KK};A_{LU}";
          TH1D* hALU = new TH1D(hname, htitle, (int)(cEdges.size() - 1), cEdges.data());
          hALU->SetDirectory(nullptr);
          phi_alu_cos_QW_[iq][iw] = hALU;
        }

        for (size_t ic = 0; ic < nC; ++ic) {
          const double cLo = cEdges[ic], cHi = cEdges[ic + 1];

          auto df_bin = df_qw.Filter([=](double mcp) { return mcp > cLo && mcp <= cHi; }, {"cos_thetaKK"});

          // helicity-separated
          auto df_pos = df_bin.Filter("REC_Event_helicity ==  1");
          auto df_neg = df_bin.Filter("REC_Event_helicity == -1");

          // tag for filenames
          auto tagRawC = hasW ? Form("Q2_%.2f_%.2f__W_%.2f_%.2f__tprime_%.3f_%.3f", qLo, qHi, wLo, wHi, cLo, cHi) : Form("Q2_%.2f_%.2f__tprime_%.3f_%.3f", qLo, qHi, cLo, cHi);
          std::string tagRaw = tagRawC;

          const auto yp = FitPhiMassYieldAndSave(df_pos, outDirPerModel, tagRaw, "P", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);
          const auto ym = FitPhiMassYieldAndSave(df_neg, outDirPerModel, tagRaw, "M", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);

          const double Np = yp.N, dNp = yp.dN;
          const double Nm = ym.N, dNm = ym.dN;

          double A = 0.0, dA = 0.0;
          const double denom = (Np + Nm);
          if (denom > 0.0 && beamPol != 0.0) {
            const double Araw = (Np - Nm) / denom;

            // error propagation
            const double dAraw_dNp = (2.0 * Nm) / (denom * denom);
            const double dAraw_dNm = -(2.0 * Np) / (denom * denom);
            const double varAraw = dAraw_dNp * dAraw_dNp * dNp * dNp + dAraw_dNm * dAraw_dNm * dNm * dNm;

            A = Araw / beamPol;
            dA = std::sqrt(std::max(0.0, varAraw)) / std::fabs(beamPol);
          }

          const int b = (int)(ic + 1);
          phi_alu_cos_QW_[iq][iw]->SetBinContent(b, A);
          phi_alu_cos_QW_[iq][iw]->SetBinError(b, dA);
        }  // it
      }  // iw
    }  // iq
  }

  void MakePhiBSATrentoPhiMassFitCanvases3D(const BinManager& bins, const std::string& outDirPerModel, int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                            bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    const auto& q2Edges = bins.GetQ2Bins();
    const auto& wEdges = bins.GetWBins();
    const auto& phiEdges = bins.GetTrentoPhiBins();  // NEW
    const bool hasW = !wEdges.empty();

    const size_t nQ = (q2Edges.size() > 1) ? (q2Edges.size() - 1) : 0;
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;
    const size_t nPhi = (phiEdges.size() > 1) ? (phiEdges.size() - 1) : 0;

    if (!nQ || !nPhi) {
      std::cerr << "[MakePhiBSATrentoPhiMassFitCanvases3D] missing Q2 or Trento-phi binning.\n";
      return;
    }

    // Verify strictly increasing φ edges
    for (size_t i = 1; i < phiEdges.size(); ++i) {
      if (!(phiEdges[i] > phiEdges[i - 1])) {
        std::cerr << "[MakePhiBSATrentoPhiMassFitCanvases3D] Trento-phi edges not strictly increasing!\n";
        return;
      }
    }

    // allocate cache: (Q,W)->TH1(phi_trento)
    if (!phi_bsa_trentophi_QW_.empty()) {
      for (auto& row : phi_bsa_trentophi_QW_)
        for (auto* h : row) delete h;
      phi_bsa_trentophi_QW_.clear();
    }
    phi_bsa_trentophi_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));

    static std::atomic<unsigned long> uidH{0};

    for (size_t iq = 0; iq < nQ; ++iq) {
      const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];
      auto df_q = rdf.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;
        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        // cache hist for this (Q,W)
        {
          auto hname = Form("phi_bsa_phiTrento_Q%zu_W%zu_%lu", iq, iw, uidH.fetch_add(1));
          auto htitle = "A_{LU}(#phi_{Trento});#phi_{Trento} [deg];A_{LU}";
          TH1D* hALU = new TH1D(hname, htitle, (int)nPhi, phiEdges.data());
          hALU->SetDirectory(nullptr);
          phi_bsa_trentophi_QW_[iq][iw] = hALU;
        }

        for (size_t ip = 0; ip < nPhi; ++ip) {
          const double pLo = phiEdges[ip], pHi = phiEdges[ip + 1];

          // NOTE: using column name "phi" (Trento φ in your codebase)
          auto df_bin = df_qw.Filter([=](double ph) { return ph > pLo && ph <= pHi; }, {"phi"});

          auto df_pos = df_bin.Filter("REC_Event_helicity ==  1");
          auto df_neg = df_bin.Filter("REC_Event_helicity == -1");

          // tag for filenames (so you get inv-mass fit PDFs inside this workflow)
          auto tagRawC =
              hasW ? Form("Q2_%.2f_%.2f__W_%.2f_%.2f__phiTrento_%.1f_%.1f", qLo, qHi, wLo, wHi, pLo, pHi) : Form("Q2_%.2f_%.2f__phiTrento_%.1f_%.1f", qLo, qHi, pLo, pHi);
          std::string tagRaw = tagRawC;

          const auto yp = FitPhiMassYieldAndSave(df_pos, outDirPerModel, tagRaw, "P", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);
          const auto ym = FitPhiMassYieldAndSave(df_neg, outDirPerModel, tagRaw, "M", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);

          const double Np = yp.N, dNp = yp.dN;
          const double Nm = ym.N, dNm = ym.dN;

          double A = 0.0, dA = 0.0;
          const double denom = (Np + Nm);
          if (denom > 0.0 && beamPol != 0.0) {
            const double Araw = (Np - Nm) / denom;

            const double dAraw_dNp = (2.0 * Nm) / (denom * denom);
            const double dAraw_dNm = -(2.0 * Np) / (denom * denom);
            const double varAraw = dAraw_dNp * dAraw_dNp * dNp * dNp + dAraw_dNm * dAraw_dNm * dNm * dNm;

            A = Araw / beamPol;
            dA = std::sqrt(std::max(0.0, varAraw)) / std::fabs(beamPol);
          }

          const int b = (int)(ip + 1);
          phi_bsa_trentophi_QW_[iq][iw]->SetBinContent(b, A);
          phi_bsa_trentophi_QW_[iq][iw]->SetBinError(b, dA);
        }
      }
    }
  }

  // Same signature as the cos(theta_KK) version, but binning on z_phi
  inline void MakePhiALUZPhiMassFitCanvases3D(const BinManager& bins, const std::string& outDirPerModel, int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                              bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    gSystem->Exec(Form("mkdir -p \"%s\"", outDirPerModel.c_str()));

    const auto& q2Edges = bins.GetQ2Bins();
    const auto& wEdges = bins.GetWBins();
    const auto& zEdges = bins.GetZPhiBins();  // <– NEW
    const bool hasW = !wEdges.empty();

    const size_t nQ = (q2Edges.size() > 1) ? (q2Edges.size() - 1) : 0;
    const size_t nZ = (zEdges.size() > 1) ? (zEdges.size() - 1) : 0;
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;

    if (!nQ || !nZ) {
      std::cerr << "[MakePhiALUZPhiMassFitCanvases3D] missing Q2 or z_phi binning.\n";
      return;
    }

    // check strictly increasing z_phi bin edges
    for (size_t i = 1; i < zEdges.size(); ++i) {
      if (!(zEdges[i] > zEdges[i - 1])) {
        std::cerr << "[MakePhiALUZPhiMassFitCanvases3D] z_phi bin edges not strictly increasing!\n";
        return;
      }
    }

    // clear old cache if any
    if (!phi_alu_zphi_QW_.empty()) {
      for (auto& row : phi_alu_zphi_QW_)
        for (auto* h : row) delete h;
      phi_alu_zphi_QW_.clear();
    }
    phi_alu_zphi_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));
    static std::atomic<unsigned long> uidH{0};

    for (size_t iq = 0; iq < nQ; ++iq) {
      const double qLo = q2Edges[iq];
      const double qHi = q2Edges[iq + 1];
      auto df_q = rdf.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;
        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        // Cache A_LU(z_phi) hist for this (Q2,W)
        {
          auto hname = Form("phi_alu_zphi_Q%zu_W%zu_%lu", iq, iw, uidH.fetch_add(1));
          auto htitle = "A_{LU}(z_{#phi});z_{#phi};A_{LU}";
          TH1D* hALU = new TH1D(hname, htitle, (int)nZ, zEdges.data());
          hALU->SetDirectory(nullptr);
          phi_alu_zphi_QW_[iq][iw] = hALU;
        }

        // Now loop over z_phi bins, fit φ mass separately in each bin
        for (size_t iz = 0; iz < nZ; ++iz) {
          const double zLo = zEdges[iz];
          const double zHi = zEdges[iz + 1];

          auto df_bin = df_qw.Filter([=](double z) { return z > zLo && z <= zHi; }, {"z_phi"});

          auto df_pos = df_bin.Filter("REC_Event_helicity ==  1");
          auto df_neg = df_bin.Filter("REC_Event_helicity == -1");

          auto tagRawC = hasW ? Form("Q2_%.2f_%.2f__W_%.2f_%.2f__zphi_%.3f_%.3f", qLo, qHi, wLo, wHi, zLo, zHi) : Form("Q2_%.2f_%.2f__zphi_%.3f_%.3f", qLo, qHi, zLo, zHi);
          std::string tagRaw = tagRawC;

          const auto yp = FitPhiMassYieldAndSave(df_pos, outDirPerModel, tagRaw, "P", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);
          const auto ym = FitPhiMassYieldAndSave(df_neg, outDirPerModel, tagRaw, "M", nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac);

          const double Np = yp.N;
          const double dNp = yp.dN;
          const double Nm = ym.N;
          const double dNm = ym.dN;

          double A = 0.0, dA = 0.0;
          const double denom = (Np + Nm);

          if (denom > 0.0 && beamPol != 0.0) {
            const double Araw = (Np - Nm) / denom;
            const double dAraw = std::sqrt((4.0 * (Nm * Nm * dNp * dNp + Np * Np * dNm * dNm)) / std::pow(denom * denom, 2.0));

            A = Araw / beamPol;
            dA = dAraw / std::fabs(beamPol);
          }

          // fill (center of z_phi bin)
          const double zCenter = 0.5 * (zLo + zHi);
          auto* hALU = phi_alu_zphi_QW_[iq][iw];
          if (hALU) {
            const int bin = hALU->FindBin(zCenter);
            hALU->SetBinContent(bin, A);
            hALU->SetBinError(bin, dA);
          }
        }  // iz
      }  // iw
    }  // iq
  }

  // Build A_LU(cos(theta_KK)) using the sin(phi_Trento) "moment" method,
  // restricted to a given K+K- invariant-mass window [mMin, mMax].
  // Formula (per (Q2, W, cosθ_KK) bin):
  //   A_LU = ( Σ_+ sinφ  -  Σ_- sinφ ) / ( 2 Σ sin^2φ )  * 1/beamPol
  // where φ is Trento phi in degrees (converted to radians for sin).
  inline void MakePhiALUCosThetaSinPhiMoment3D(const BinManager& bins, const std::string& outDirPerModel, double mMin, double mMax, double beamPol = 1.0) {
    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    const auto& q2Edges = bins.GetQ2Bins();
    const auto& wEdges = bins.GetWBins();
    const auto& cEdges = bins.GetCosThetaKKBins();

    const bool hasW = !wEdges.empty();
    const size_t nQ = (q2Edges.size() > 1) ? (q2Edges.size() - 1) : 0;
    const size_t nC = (cEdges.size() > 1) ? (cEdges.size() - 1) : 0;
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;

    if (!nQ || !nC) {
      std::cerr << "[MakePhiALUCosThetaSinPhiMoment3D] missing Q2 or cos(theta_KK) binning.\n";
      return;
    }

    // sanity check: cos(theta_KK) bins strictly increasing
    for (size_t i = 1; i < cEdges.size(); ++i) {
      if (!(cEdges[i] > cEdges[i - 1])) {
        std::cerr << "[MakePhiALUCosThetaSinPhiMoment3D] cos(theta_KK) bin edges not strictly increasing!\n";
        return;
      }
    }

    // clear any previous cache
    if (!phi_alu_cos_QW_.empty()) {
      for (auto& row : phi_alu_cos_QW_) {
        for (auto* h : row) {
          delete h;
        }
      }
      phi_alu_cos_QW_.clear();
    }

    phi_alu_cos_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));
    static std::atomic<unsigned long> uidH{0};

    auto df_all = rdf;  // base dataframe in this plotter

    const double deg2rad = 3.14159265358979323846 / 180.0;

    for (size_t iq = 0; iq < nQ; ++iq) {
      const double qLo = q2Edges[iq];
      const double qHi = q2Edges[iq + 1];

      auto df_q = df_all.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;

        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        // Histogram A_LU(cosθ_KK) for this (Q2, W) slice
        {
          auto hname = Form("phi_alu_costh_sinphi_Q%zu_W%zu_%lu", iq, iw, uidH.fetch_add(1));
          auto htitle = "A_{LU}^{sin#phi}(cos#theta_{KK});cos#theta_{KK};A_{LU}^{sin#phi}";
          TH1D* hALU = new TH1D(hname, htitle, static_cast<int>(nC), cEdges.data());
          hALU->SetDirectory(nullptr);
          phi_alu_cos_QW_[iq][iw] = hALU;
        }

        // Loop over cos(theta_KK) bins and compute the sinφ moment
        for (size_t ic = 0; ic < nC; ++ic) {
          const double cLo = cEdges[ic];
          const double cHi = cEdges[ic + 1];

          // Select cos(theta_KK) bin
          auto df_costh = df_qw.Filter([=](double c) { return c > cLo && c <= cHi; }, {"cos_thetaKK"});

          // Restrict to invariant-mass window
          auto df_mass = df_costh.Filter([=](float m) { return m > mMin && m < mMax; }, {"invMass_KpKm"});

          // Define sin(phi_Trento) (phi is in degrees in the tree)
          auto df_sin = df_mass
                            .Define("sinphi", [=](double phi_deg) { return std::sin(phi_deg * deg2rad); }, {"phi"})
                            // helicity-weighted sinφ: +1 * sinφ for +hel, -1 * sinφ for -hel
                            .Define("w_sinphi", [](short hel, double s) { return static_cast<double>(hel) * s; }, {"REC_Event_helicity", "sinphi"})
                            // sin²φ for the denominator
                            .Define("sin2", [](double s) { return s * s; }, {"sinphi"});

          // Numerator: Σ hel * sinφ  over all events in this bin
          auto sum_num = df_sin.Sum("w_sinphi");
          // Denominator piece: Σ sin²φ over all events in this bin
          auto sum_s2 = df_sin.Sum("sin2");

          const double Nnum = static_cast<double>(*sum_num);  // Σ λ_i sinφ_i
          const double S2 = static_cast<double>(*sum_s2);     // Σ sin²φ_i

          double A = 0.0;
          double dA = 0.0;

          if (S2 > 0.0 && beamPol != 0.0) {
            // A_raw = (Σ_+ sinφ - Σ_- sinφ) / (2 Σ sin²φ)
            const double Araw = Nnum / (2.0 * S2);

            // Approximate variance: var(Nnum) ≈ Σ sin²φ = S2
            // => var(Araw) ≈ S2 / (4 S2²) = 1 / (4 S2)
            const double dAraw = 1.0 / (2.0 * std::sqrt(S2));

            A = Araw / beamPol;
            dA = dAraw / std::fabs(beamPol);
          }

          TH1D* hALU = phi_alu_cos_QW_[iq][iw];
          if (hALU) {
            const double cCenter = 0.5 * (cLo + cHi);
            const int b = hALU->FindBin(cCenter);
            hALU->SetBinContent(b, A);
            hALU->SetBinError(b, dA);
          }
        }  // ic (cosθ_KK)
      }  // iw (W)
    }  // iq (Q²)
  }
  // A_LU(cos(theta_KK)) from a fit of A(phi_Trento) = A * sin(phi)/(1 + b cos(phi))
  // in each (Q2, W, cos(theta_KK)) bin, restricted to an invMass_KpKm window [mMin, mMax].
  inline void MakePhiALUCosTheta_SinOver1PlusbCosFit3D(const BinManager& bins, const std::string& outDirPerModel, double mMin, double mMax, double beamPol = 1.0) {
    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    const auto& q2Edges = bins.GetQ2Bins();
    const auto& wEdges = bins.GetWBins();
    const auto& cEdges = bins.GetCosThetaKKBins();
    const auto& phiEdges = bins.GetTrentoPhiBins();  // φ binning

    const bool hasW = !wEdges.empty();
    const size_t nQ = (q2Edges.size() > 1) ? (q2Edges.size() - 1) : 0;
    const size_t nC = (cEdges.size() > 1) ? (cEdges.size() - 1) : 0;
    const size_t nPhi = (phiEdges.size() > 1) ? (phiEdges.size() - 1) : 0;
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;

    if (!nQ || !nC || !nPhi) {
      std::cerr << "[MakePhiALUCosTheta_SinOver1PlusbCosFit3D] missing Q2, cos(theta_KK) or phi binning.\n";
      return;
    }

    // sanity: strictly increasing edges
    for (size_t i = 1; i < cEdges.size(); ++i) {
      if (!(cEdges[i] > cEdges[i - 1])) {
        std::cerr << "[MakePhiALUCosTheta_SinOver1PlusbCosFit3D] cos(theta_KK) edges not strictly increasing!\n";
        return;
      }
    }
    for (size_t i = 1; i < phiEdges.size(); ++i) {
      if (!(phiEdges[i] > phiEdges[i - 1])) {
        std::cerr << "[MakePhiALUCosTheta_SinOver1PlusbCosFit3D] phi edges not strictly increasing!\n";
        return;
      }
    }

    // reset A_LU cache
    if (!phi_alu_cos_QW_.empty()) {
      for (auto& row : phi_alu_cos_QW_)
        for (auto* h : row) delete h;
      phi_alu_cos_QW_.clear();
    }
    phi_alu_cos_QW_.assign(nQ, std::vector<TH1D*>(nW, nullptr));

    static std::atomic<unsigned long> uidH{0};  // for hist names
    static std::atomic<unsigned long> uidC{0};  // for canvas names

    auto df_all = rdf;

    for (size_t iq = 0; iq < nQ; ++iq) {
      const double qLo = q2Edges[iq];
      const double qHi = q2Edges[iq + 1];

      auto df_q = df_all.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;

        auto df_qw = df_q.Filter([=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        // A_LU(cosθ_KK) hist for this (Q2, W)
        {
          auto hname = Form("phi_alu_costh_sinOver1plusbcos_Q%zu_W%zu_%lu", iq, iw, uidH.fetch_add(1));
          auto htitle = "A_{LU}(cos#theta_{KK});cos#theta_{KK};A_{LU}";
          TH1D* hALU = new TH1D(hname, htitle, (int)nC, cEdges.data());
          hALU->SetDirectory(nullptr);
          phi_alu_cos_QW_[iq][iw] = hALU;
        }

        for (size_t ic = 0; ic < nC; ++ic) {
          const double cLo = cEdges[ic];
          const double cHi = cEdges[ic + 1];

          // cos(theta_KK) slice
          auto df_costh = df_qw.Filter([=](double c) { return c > cLo && c <= cHi; }, {"cos_thetaKK"});

          // inv-mass window (branch is float!)
          auto df_mass = df_costh.Filter([=](float m) { return m > mMin && m < mMax; }, {"invMass_KpKm"});

          // helicity separation
          auto df_pos = df_mass.Filter("REC_Event_helicity ==  1");
          auto df_neg = df_mass.Filter("REC_Event_helicity == -1");

          // N^+(phi), N^-(phi)
          auto hPlus_ptr = df_pos.Histo1D({"hPlus_tmp", "hPlus_tmp", (int)nPhi, phiEdges.data()}, "phi");
          auto hMinus_ptr = df_neg.Histo1D({"hMinus_tmp", "hMinus_tmp", (int)nPhi, phiEdges.data()}, "phi");

          std::unique_ptr<TH1D> hPlus(dynamic_cast<TH1D*>(hPlus_ptr->Clone("hPlus_local")));
          std::unique_ptr<TH1D> hMinus(dynamic_cast<TH1D*>(hMinus_ptr->Clone("hMinus_local")));
          if (hPlus) hPlus->SetDirectory(nullptr);
          if (hMinus) hMinus->SetDirectory(nullptr);

          // BSA(phi)
          std::unique_ptr<TH1D> hBSA(new TH1D("hBSA_local", "BSA vs #phi_{Trento};#phi_{Trento} [deg];A_{LU}", (int)nPhi, phiEdges.data()));
          hBSA->SetDirectory(nullptr);
          hBSA->SetTitle("");

          double totalCounts = 0.0;
          for (int ib = 1; ib <= (int)nPhi; ++ib) {
            const double Np = hPlus ? hPlus->GetBinContent(ib) : 0.0;
            const double Nm = hMinus ? hMinus->GetBinContent(ib) : 0.0;
            const double dNp = std::sqrt(Np);
            const double dNm = std::sqrt(Nm);
            totalCounts += Np + Nm;

            double A = 0.0;
            double dA = 0.0;
            const double denom = Np + Nm;

            if (denom > 0.0 && beamPol != 0.0) {
              const double Araw = (Np - Nm) / denom;

              const double dAraw_dNp = (2.0 * Nm) / (denom * denom);
              const double dAraw_dNm = -(2.0 * Np) / (denom * denom);
              const double varAraw = dAraw_dNp * dAraw_dNp * dNp * dNp + dAraw_dNm * dAraw_dNm * dNm * dNm;

              const double dAraw = std::sqrt(std::max(0.0, varAraw));

              A = Araw / beamPol;
              dA = dAraw / std::fabs(beamPol);
            }

            hBSA->SetBinContent(ib, A);
            hBSA->SetBinError(ib, dA);
          }

          double A_fit = 0.0;
          double dA_fit = 0.0;

          // Do the fit and draw/save a canvas, similar style to KK-mass fits
          if (totalCounts > 0.0 && hBSA) {
            const double xMin = phiEdges.front();
            const double xMax = phiEdges.back();

            TF1 fALU("fALU", "[0]*sin(x*TMath::DegToRad())/(1.0 + [1]*cos(x*TMath::DegToRad()))", xMin, xMax);
            fALU.SetParName(0, "A_{LU}");
            fALU.SetParName(1, "b");
            fALU.SetParameters(0.0, 0.0);     // initial guesses
            fALU.SetParLimits(1, -0.9, 0.9);  // avoid denominator singularities

            hBSA->Fit(&fALU, "R0Q");  // quiet, no draw yet

            A_fit = fALU.GetParameter(0);
            dA_fit = fALU.GetParError(0);
            const double b_fit = fALU.GetParameter(1);
            const double db_fit = fALU.GetParError(1);
            const double chi2 = fALU.GetChisquare();
            const int ndf = fALU.GetNDF();

            // --- NOW draw & save the fit canvas (like invMass_KK fits) ---
            std::string tagRaw =
                hasW ? Form("Q2_%.2f_%.2f__W_%.2f_%.2f__costh_%.3f_%.3f", qLo, qHi, wLo, wHi, cLo, cHi) : Form("Q2_%.2f_%.2f__costh_%.3f_%.3f", qLo, qHi, cLo, cHi);
            std::string tag = tagRaw;
            std::replace(tag.begin(), tag.end(), '.', '_');

            const auto cname = Form("c_BSAphi_%s_%lu", tag.c_str(), uidC.fetch_add(1));

            std::unique_ptr<TCanvas> c(new TCanvas(cname, "BSA vs #phi", 1200, 1000));
            gStyle->SetOptStat(0);
            gStyle->SetOptFit(0);
            c->SetTicks(1, 1);
            c->SetTopMargin(0.04);
            c->SetRightMargin(0.03);
            c->SetBottomMargin(0.11);
            c->SetLeftMargin(0.13);

            hBSA->SetMarkerStyle(20);
            hBSA->SetMarkerSize(1.1);
            hBSA->SetLineWidth(2);
            hBSA->Draw("E1");

            fALU.SetLineColor(kRed + 1);
            fALU.SetLineWidth(2);
            fALU.Draw("SAME");

            // horizontal zero line
            double ymin = hBSA->GetMinimum();
            double ymax = hBSA->GetMaximum();
            double absMax = std::max(std::fabs(ymin), std::fabs(ymax));
            absMax = std::max(absMax, 0.10);
            hBSA->GetYaxis()->SetRangeUser(-1.2 * absMax, 1.2 * absMax);

            TLine zline(xMin, 0.0, xMax, 0.0);
            zline.SetLineStyle(2);
            zline.SetLineColor(kGray + 2);
            zline.Draw("SAME");

            // Legend & text
            TLegend leg(0.55, 0.75, 0.90, 0.90);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextSize(0.035);
            leg.AddEntry(hBSA.get(), "data", "lep");
            leg.AddEntry(&fALU, "fit A sin#phi/(1 + b cos#phi)", "l");
            leg.Draw();

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.15, 0.93,
                            hasW ? Form("Q^{2}[%.2f, %.2f]  W[%.2f, %.2f]  cos#theta_{KK}[%.2f, %.2f]", qLo, qHi, wLo, wHi, cLo, cHi)
                                 : Form("Q^{2}[%.2f, %.2f]  cos#theta_{KK}[%.2f, %.2f]", qLo, qHi, cLo, cHi));

            latex.SetTextSize(0.035);
            latex.DrawLatex(0.16, 0.86, Form("A_{LU} = %.3f #pm %.3f", A_fit, dA_fit));
            latex.DrawLatex(0.16, 0.82, Form("b = %.3f #pm %.3f", b_fit, db_fit));
            latex.DrawLatex(0.16, 0.78, Form("#chi^{2}/ndf = %.1f/%d", chi2, ndf));

            c->SaveAs(Form("%s/BSAphi_%s.pdf", outDirPerModel.c_str(), tag.c_str()));
          }

          // finally, save amplitude into A_LU(cosθ_KK) cache
          TH1D* hALU = phi_alu_cos_QW_[iq][iw];
          if (hALU) {
            const double cCenter = 0.5 * (cLo + cHi);
            const int ibin = hALU->FindBin(cCenter);
            hALU->SetBinContent(ibin, A_fit);
            hALU->SetBinError(ibin, dA_fit);
          }
        }  // ic
      }  // iw
    }  // iq
  }

 private:
  std::vector<std::string> disvars = {"Q2", "xB", "t", "W", "phi"};
  std::vector<std::string> Phivars = {"Q2", "xB", "t", "W", "phi", "mtprime"};
  std::map<std::string, std::pair<double, double>> axisRanges = {{"Q2", {0.0, 15.0}}, {"xB", {0.0, 1.0}},       {"W", {1.0, 10.0}},
                                                                 {"t", {0.0, 10.0}},  {"mtprime", {0.0, 10.0}}, {"phi", {-180.0, 180.0}}};
  std::map<std::string, std::pair<double, double>> kinematicAxisRanges = {
      {"recel_p", {-0.05, 13.0}},     {"recel_theta", {-0.01, 1.0}},  {"recel_phi", {-0.01, 6.1}},  {"recpho_p", {-0.01, 10.0}},
      {"recpho2_p", {-0.01, 4.0}},    {"recpho_theta", {-0.01, 1.0}}, {"recpho_phi", {-0.01, 6.1}}, {"recpro_p", {-0.01, 8.0}},
      {"recpro_theta", {-0.01, 2.0}}, {"recpro_phi", {-0.01, 6.1}}
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
  std::vector<std::vector<std::vector<double>>> phi_mean_xB_QW_;
  std::vector<std::vector<std::vector<double>>> phi_mean_W_QW_;
  std::vector<std::vector<std::vector<double>>> phi_mean_GammaV_QW_;
  
  std::vector<ROOT::RDF::RResultPtr<TH1>> kinematicHistos, disHistos;
  std::vector<std::vector<TH1D*>> phi_dsdt_QW_;
  std::vector<std::vector<TH1D*>> phi_nsig_QW_;   // raw fit yields Nsig per (Q2,W,t) bin
  std::vector<std::vector<TH1D*>> phi_alu_costh_QW_;
  std::vector<std::vector<TH1D*>> phi_alu_cos_QW_;
  std::vector<std::vector<TH1D*>> phi_bsa_trentophi_QW_;
  std::vector<std::vector<TH1D*>> phi_alu_zphi_QW_;

  std::vector<std::shared_ptr<TH1>> acceptHistos;
  std::vector<std::vector<TH1D*>> phi_accept_QW_;
  std::vector<std::vector<TH1D*>> phi_rad_QW_;

  //
  ROOT::RDF::RNode rdf;

  // ROOT::RDF::RNode rdf_dvcs_data;
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
  // rdfs for phi analysis
  std::optional<ROOT::RDF::RNode> rdf_gen_phimc;
  std::optional<ROOT::RDF::RNode> rdf_accept_phimc;
  std::optional<ROOT::RDF::RNode> rdf_phimc_radRatio;

  AccEffProvider accEffProvider_ = [](double, double, double, double, double, double) { return 1.0; };
};

#endif  // DISANA_PLOTTER_H
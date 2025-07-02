#ifndef DISANAMATH_H
#define DISANAMATH_H

// ROOT headers and standard libraries
#include <ROOT/RDataFrame.hxx>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "TLorentzVector.h"
#include "TMath.h"

// Constants used throughout the analysis
constexpr double pi = 3.14159265358979323846;
const double m_e = 0.000511;  // Electron mass in GeV
const double m_p = 0.938272;  // Proton mass in GeV

// --- Utility Functions (unnamed namespace) ---
// Converts spherical coordinates (p, θ, φ) to Cartesian 3-vector
namespace {
TVector3 SphericalToCartesian(double p, double theta, double phi) {
  double px = p * std::sin(theta) * std::cos(phi);
  double py = p * std::sin(theta) * std::sin(phi);
  double pz = p * std::cos(theta);
  return TVector3(px, py, pz);
}

// Builds a 4-vector using spherical angles and mass
TLorentzVector Build4Vector(double p, double theta, double phi, double mass) {
  TVector3 vec = SphericalToCartesian(p, theta, phi);
  double E = std::sqrt(vec.Mag2() + mass * mass);  // E² = p² + m²
  return TLorentzVector(vec, E);
}
}  // end anonymous namespace

// --- BinManager class ---
// Manages Q², t, x_B bins for cross-section calculations
class BinManager {
 public:
  BinManager() {
    q2_bins_ = {1.0, 2.0, 4.0, 6.0};
    t_bins_ = {0.1, 0.3, 0.6, 1.0};
    xb_bins_ = {0.1, 0.2, 0.4, 0.6};
  }

  // Getters for bin edges
  const std::vector<double> &GetQ2Bins() const { return q2_bins_; }
  const std::vector<double> &GetTBins() const { return t_bins_; }
  const std::vector<double> &GetXBBins() const { return xb_bins_; }

  // Setters if dynamic binning is needed
  void SetQ2Bins(const std::vector<double> &bins) { q2_bins_ = bins; }
  void SetTBins(const std::vector<double> &bins) { t_bins_ = bins; }
  void SetXBBins(const std::vector<double> &bins) { xb_bins_ = bins; }

 private:
  std::vector<double> q2_bins_, t_bins_, xb_bins_;
};

// --- DISANAMath class ---
// Central class for computing DVCS kinematics, exclusivity variables, and cross-sections
class DISANAMath {
 private:
  // Kinematic variables
  double Q2_, xB_, t_, phi_deg_, W_, nu_, y_;
  bool applyCorrection = false;
  THnSparseD *correctionHist = nullptr;

  // Exclusivity variables
  double mx2_ep_;         // Missing mass² of ep system
  double emiss_;          // Missing energy
  double ptmiss_;         // Transverse missing momentum
  double mx2_epg_;        // Missing mass² of epγ system
  double delta_phi_;      // Coplanarity angle Δφ (proton - q)
  double theta_gg_;       // Angle between γ and missing vector
  double mx2_egamma_;     // Invariant mass of e + γ
  double Theta_e_gamma_;  // Angle between e' and γ
  double DeltaE_;         // Energy balance: (initial - final)

 public:
  DISANAMath() = default;

  // Constructor: Takes measured quantities and builds all required kinematic variables
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double g_p, double g_theta, double g_phi) {
    TLorentzVector electron_in(0, 0, e_in_E, e_in_E);  // beam along z
    TLorentzVector electron_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector proton_in(0, 0, 0, m_p);  // at rest
    TLorentzVector proton_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector photon = Build4Vector(g_p, g_theta, g_phi, 0.0);  // massless

    ComputeKinematics(electron_in, electron_out, proton_in, proton_out, photon);
  }
  void SetApplyCorrPi0BKG(bool enable) { applyCorrection = enable; }
  void SetCorrHist(THnSparseD *hist) { correctionHist = hist; }

  // Accessor methods for computed values
  double GetQ2() const { return Q2_; }
  double GetxB() const { return xB_; }
  double GetT() const { return t_; }
  double GetPhi() const { return phi_deg_; }
  double GetW() const { return W_; }
  double GetNu() const { return nu_; }

  double Gety() const { return y_; }
  static std::vector<TH1D *> FlattenHists(const std::vector<std::vector<std::vector<TH1D *>>> &h3d);

  double GetMx2_ep() const { return mx2_ep_; }
  double GetEmiss() const { return emiss_; }
  double GetPTmiss() const { return ptmiss_; }
  double GetMx2_epg() const { return mx2_epg_; }
  double GetDeltaPhi() const { return delta_phi_; }

  double GetTheta_gamma_gamma() const { return theta_gg_; }
  double GetMx2_egamma() const { return mx2_egamma_; }
  double GetTheta_e_gamma() const { return Theta_e_gamma_; }
  double GetDeltaE() const { return DeltaE_; }

  double GetCorrectionFactor(double Q2, double t, double xB, double phi_deg) const {
    if (!applyCorrection || !correctionHist) return 1.0;

    int bins[4] = {correctionHist->GetAxis(0)->FindBin(Q2), correctionHist->GetAxis(1)->FindBin(t), correctionHist->GetAxis(2)->FindBin(xB),
                   correctionHist->GetAxis(3)->FindBin(phi_deg)};

    return correctionHist->GetBinContent(bins);
  }

  double ComputePhiH(const TVector3 &q1_, const TVector3 &k1_, const TVector3 &q2_) const {
    double t1 = ((q1_.Cross(k1_)).Dot(q2_)) / std::abs((q1_.Cross(k1_)).Dot(q2_));

    TVector3 t2 = q1_.Cross(k1_);
    TVector3 t3 = q1_.Cross(q2_);
    double n2 = t2.Mag();
    double n3 = t3.Mag();
    double t4 = (t2.Dot(t3)) / (n2 * n3);

    return t1 * std::acos(t4) * 180. / pi + 180.;
  }
  // --- Core computation function ---
  void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &photon) {
    TLorentzVector q = electron_in - electron_out;  // virtual photon

    Q2_ = -q.Mag2();
    nu_ = q.E();
    y_ = nu_ / electron_in.E();
    W_ = (proton_in + q).Mag();
    xB_ = Q2_ / (2.0 * proton_in.Dot(q));
    t_ = std::abs((proton_in - proton_out).Mag2());  // Mandelstam t

    // Azimuthal angle φ between lepton and hadron planes
    TVector3 n_L = electron_in.Vect().Cross(electron_out.Vect()).Unit();
    TVector3 n_H = q.Vect().Cross(proton_out.Vect()).Unit();
    double cos_phi = n_L.Dot(n_H);
    double sin_phi = (n_L.Cross(n_H)).Dot(q.Vect().Unit());
    double phi = std::atan2(sin_phi, cos_phi) + pi;  // Ensure φ is in [0, 2π]
    phi_deg_ = phi * 180.0 / pi;

    // Composite 4-vectors
    TLorentzVector ep_combined = electron_out + photon;
    TLorentzVector total_initial = electron_in + proton_in;
    TLorentzVector total_final = electron_out + proton_out + photon;
    TLorentzVector missing = total_initial - total_final;

    // Exclusivity observables
    mx2_ep_ = (total_initial - electron_out - proton_out).Mag2();
    emiss_ = missing.E();
    ptmiss_ = missing.Vect().Perp();
    mx2_epg_ = missing.Mag2();

    // Coplanarity Δφ between outgoing proton and q-vector in transverse plane
    TVector3 q_vec = q.Vect();
    TVector3 electron_vec = electron_in.Vect();
    TVector3 photon_vec = photon.Vect();
    TVector3 p_vec = proton_out.Vect();
    double phi_q = std::atan2(q_vec.Y(), q_vec.X());
    double phi_p = std::atan2(p_vec.Y(), p_vec.X());
    double delta_phi_rad = TVector2::Phi_mpi_pi(phi_p - phi_q);
    // delta_phi_ = std::abs(delta_phi_rad) * 180.0 / pi;
    delta_phi_ = abs(ComputePhiH(q_vec, electron_vec, photon_vec) - ComputePhiH(q_vec, electron_vec, -p_vec));

    // θ(γ, missing): photon direction vs. missing momentum
    theta_gg_ = photon.Angle((total_initial - (electron_out + proton_out)).Vect()) * 180.0 / pi;

    // Invariant mass of e + γ
    // mx2_egamma_ = (proton_in - (electron_out + photon)).Mag2();
    mx2_egamma_ = (electron_in + proton_in - electron_out - photon).Mag2();

    // Angle between proton and photon
    Theta_e_gamma_ = electron_out.Angle(photon.Vect()) * 180.0 / pi;

    // Energy imbalance (should be 0 for exclusive DVCS)
    DeltaE_ = (electron_in.E() + proton_in.E()) - (electron_out.E() + proton_out.E() + photon.E());
  }

  // Integrated luminosity in cm⁻² (example, you can scale it out if not known)
  std::vector<std::vector<std::vector<TH1D *>>> ComputeDVCS_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
    TStopwatch timer;
    timer.Start();
    constexpr double phi_min = 0.0;
    constexpr double phi_max = 360.0;
    constexpr int n_phi_bins = 18;

    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const auto &xb_bins = bins.GetXBBins();

    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t = t_bins.size() - 1;
    const size_t n_xb = xb_bins.size() - 1;

    std::vector<std::vector<std::vector<TH1D *>>> histograms(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
          const double tmin = t_bins[it], tmax = t_bins[it + 1];
          const double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];

          std::string name = Form("hphi_q%.1f_t%.1f_xb%.2f", qmin, tmin, xbmin);
          std::string title = Form(
              "d#sigma/d#phi (Q^{2}=[%.1f,%.1f], "
              "t=[%.1f,%.1f], x_{B}=[%.2f,%.2f])",
              qmin, qmax, tmin, tmax, xbmin, xbmax);

          histograms[ix][iq][it] = new TH1D(name.c_str(), title.c_str(), n_phi_bins, phi_min, phi_max);
          histograms[ix][iq][it]->SetDirectory(nullptr);  // 独立于当前文件
        }

    auto findBin = [](double val, const std::vector<double> &edges) -> int {
      auto it = std::upper_bound(edges.begin(), edges.end(), val);
      if (it == edges.begin() || it == edges.end()) return -1;  // 越界
      return static_cast<int>(it - edges.begin()) - 1;
    };

    auto filler = [&](double Q2, double t, double xB, double phi) {
      const int iq = findBin(Q2, q2_bins);
      const int it = findBin(t, t_bins);
      const int ix = findBin(xB, xb_bins);

      if (iq >= 0 && it >= 0 && ix >= 0) {
        double factor = 1.0;
        // double q2xBtbin_size = (q2_bins[iq+1]-q2_bins[iq])*(t_bins[it+1]-t_bins[it])*(xb_bins[ix+1]-xb_bins[ix]);
        double q2xBtbin_size = 1;
        if (applyCorrection && correctionHist) {
          factor = GetCorrectionFactor(Q2, t, xB, phi);
        }
        histograms[ix][iq][it]->Fill(phi, factor / q2xBtbin_size);
      }
    };

    df.Foreach(filler, {"Q2", "t", "xB", "phi"});
    // Normalize histograms

    const double bin_width = (phi_max - phi_min) / n_phi_bins;
    for (auto &vec_q2 : histograms)
      for (auto &vec_t : vec_q2)
        for (TH1D *h : vec_t)
          if (h) {
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
              const double raw = h->GetBinContent(b);
              const double norm = raw / (luminosity * bin_width);
              const double err = std::sqrt(raw) / (luminosity * bin_width);

              // Set normalized content and error
              h->SetBinContent(b, norm);
              h->SetBinError(b, err);
            }
          }

    std::cout << "DVCS cross-sections computed in a single pass.\n";
    timer.Stop();
    std::cout << "Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return histograms;
  }

  // === NEW: 3D Beam-Spin Asymmetry ============================================Add commentMore actions
  std::vector<std::vector<std::vector<TH1D *>>> ComputeBeamSpinAsymmetry(const std::vector<std::vector<std::vector<TH1D *>>> &sigma_pos,
                                                                         const std::vector<std::vector<std::vector<TH1D *>>> &sigma_neg, double pol = 1.0) {
    if (sigma_pos.size() != sigma_neg.size()) {
      std::cerr << "ERROR: xB dim mismatch!\n";
      return {};
    }
    const size_t n_xb = sigma_pos.size();
    std::vector<std::vector<std::vector<TH1D *>>> asym(n_xb);
    for (size_t ix = 0; ix < n_xb; ++ix) {
      if (sigma_pos[ix].size() != sigma_neg[ix].size()) {
        std::cerr << "ERROR: Q² dim mismatch at xB " << ix << "!\n";
        return {};
      }

      const size_t n_q2 = sigma_pos[ix].size();
      asym[ix].resize(n_q2);
      for (size_t iq = 0; iq < n_q2; ++iq) {
        if (sigma_pos[ix][iq].size() != sigma_neg[ix][iq].size()) {
          std::cerr << "ERROR: t dim mismatch at (xB,q2)=(" << ix << ',' << iq << ")!\n";
          return {};
        }

        const size_t n_t = sigma_pos[ix][iq].size();
        asym[ix][iq].resize(n_t, nullptr);
        for (size_t it = 0; it < n_t; ++it) {
          TH1D *hp = sigma_pos[ix][iq][it];
          TH1D *hm = sigma_neg[ix][iq][it];
          if (!hp || !hm) {
            std::cerr << "WARNING: null hist @(" << ix << ',' << iq << ',' << it << ")\n";
            continue;
          }

          std::string nm = std::string(hp->GetName()) + "_BSA";
          TH1D *ha = dynamic_cast<TH1D *>(hp->Clone(nm.c_str()));
          ha->Reset();
          ha->SetTitle(("Beam Spin Asymmetry of " + std::string(hp->GetTitle())).c_str());
          for (int b = 1; b <= hp->GetNbinsX(); ++b) {
            const double Np = hp->GetBinContent(b), Nm = hm->GetBinContent(b);
            const double Ep = hp->GetBinError(b), Em = hm->GetBinError(b);
            const double den = Np + Nm, num = Np - Nm;
            double A = (den != 0) ? num / den : 0.0;
            double E = (den != 0) ? 2.0 / (den * den) * std::sqrt(std::pow(Nm * Ep, 2) + std::pow(Np * Em, 2)) : 0.0;
            ha->SetBinContent(b, A / pol);
            ha->SetBinError(b, E / pol);
          }
          asym[ix][iq][it] = ha;
        }
      }
    }
    std::cout << "Beam-Spin Asymmetries (3D) computed.\n";
    return asym;
  }

  /// pi0 correction histo for plotting
  std::vector<std::vector<std::vector<TH1D *>>> computePi0Corr(const BinManager &bins) {
    constexpr double phi_min = 0.0;
    constexpr double phi_max = 360.0;
    constexpr int n_phi_bins = 18;

    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const auto &xb_bins = bins.GetXBBins();

    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t = t_bins.size() - 1;
    const size_t n_xb = xb_bins.size() - 1;

    std::vector<std::vector<std::vector<TH1D *>>> histograms(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          std::string name = Form("pi0corr_q%d_t%d_xb%d", (int)iq, (int)it, (int)ix);
          std::string title = Form("π⁰ Corr (Q² [%g,%g], t [%g,%g], xB [%g,%g])", q2_bins[iq], q2_bins[iq + 1], t_bins[it], t_bins[it + 1], xb_bins[ix], xb_bins[ix + 1]);

          histograms[ix][iq][it] = new TH1D(name.c_str(), title.c_str(), n_phi_bins, phi_min, phi_max);
          histograms[ix][iq][it]->SetDirectory(nullptr);
        }

    // Utility lambda to find bin index
    auto findBin = [](double val, const std::vector<double> &edges) -> int {
      auto it = std::upper_bound(edges.begin(), edges.end(), val);
      if (it == edges.begin() || it == edges.end()) return -1;
      return static_cast<int>(it - edges.begin()) - 1;
    };

    Long64_t nentries = correctionHist->GetNbins();
    const Int_t ndim = correctionHist->GetNdimensions();

    for (Long64_t i = 0; i < nentries; ++i) {
      if (correctionHist->GetBinContent(i) == 0) continue;
    auto* axQ2 = correctionHist->GetAxis(0);
    auto* axT = correctionHist->GetAxis(1);
    auto* axXB = correctionHist->GetAxis(2);
    auto* axPhi = correctionHist->GetAxis(3);

    const int nQ2 = axQ2->GetNbins();
    const int nT = axT->GetNbins();
    const int nXB = axXB->GetNbins();
    const int nPhi = axPhi->GetNbins();
      std::vector<Int_t> bins = {nQ2, nT, nXB, nPhi};

      double Q2 = correctionHist->GetAxis(0)->GetBinCenter(bins[0]);
      double t = correctionHist->GetAxis(1)->GetBinCenter(bins[1]);
      double xB = correctionHist->GetAxis(2)->GetBinCenter(bins[2]);
      double phi = correctionHist->GetAxis(3)->GetBinCenter(bins[3]);

      double val = correctionHist->GetBinContent(i);

      int iq = findBin(Q2, q2_bins);
      int it = findBin(t, t_bins);
      int ix = findBin(xB, xb_bins);
      if (iq >= 0 && it >= 0 && ix >= 0) {
        histograms[ix][iq][it]->Fill(phi, val);
      }
    }

    std::cout << "π⁰ correction histograms created from THnSparse.\n";
    return histograms;
  }
};

#endif  // DISANAMATH_H

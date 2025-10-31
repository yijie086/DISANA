#ifndef DISANAMATH_H
#define DISANAMATH_H

// ROOT + std
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TVector2.h>
#include <TVector3.h>

#include <ROOT/RDataFrame.hxx>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
constexpr double pi = 3.14159265358979323846;
const double m_e = 0.000511;       // GeV
//const double m_e = 0.0;       // GeV
const double m_p = 0.938272;       // GeV
//const double m_p = 0.938;       // GeV
const double m_kMinus = 0.493677;  // GeV
const double m_kPlus = 0.493677;   // GeV

// -----------------------------------------------------------------------------
// (p,theta,phi) helpers
// -----------------------------------------------------------------------------
namespace {
TVector3 SphericalToCartesian(double p, double theta, double phi) {
  const double px = p * std::sin(theta) * std::cos(phi);
  const double py = p * std::sin(theta) * std::sin(phi);
  const double pz = p * std::cos(theta);
  return TVector3(px, py, pz);
}
TLorentzVector Build4Vector(double p, double theta, double phi, double mass) {
  TVector3 v = SphericalToCartesian(p, theta, phi);
  const double E = std::sqrt(v.Mag2() + mass * mass);
  return TLorentzVector(v, E);
}
}  // namespace

// -----------------------------------------------------------------------------
// Bin manager 
// -----------------------------------------------------------------------------
class BinManager {
 public:
  BinManager() {
    q2_bins_ = {1.0, 2.0, 4.0, 6.0};
    t_bins_ = {0.1, 0.3, 0.6, 1.0};
    xb_bins_ = {0.1, 0.2, 0.4, 0.6};
    W_bins_ = {0.1, 10.0};
  }
  const std::vector<double> &GetQ2Bins() const { return q2_bins_; }
  const std::vector<double> &GetTBins() const { return t_bins_; }
  const std::vector<double> &GetXBBins() const { return xb_bins_; }
  const std::vector<double> &GetWBins() const { return W_bins_; }
  void SetQ2Bins(const std::vector<double> &v) { q2_bins_ = v; }
  void SetTBins(const std::vector<double> &v) { t_bins_ = v; }
  void SetXBBins(const std::vector<double> &v) { xb_bins_ = v; }
  void SetWBins(const std::vector<double> &v) { W_bins_ = v; }

 private:
  std::vector<double> q2_bins_, t_bins_, xb_bins_, W_bins_;
};

// -----------------------------------------------------------------------------
// DISANAMath
// -----------------------------------------------------------------------------
struct Pi0Tag {};
class DISANAMath {
 private:
  // Kinematics
  double Q2_{}, xB_{}, t_{}, phi_deg_{}, W_{}, nu_{}, y_{};

  // Exclusivity
  double mx2_ep_{};
  double emiss_{};
  double ptmiss_{};
  double mx2_epg_{};
  double mx2_epKpKm_{};
  double delta_phi_{};
  double theta_gg_{};
  double theta_gphi_{};
  double mx2_egamma_{};
  double mx2_eKpKm_{};
  double mx2_epKp_{-1.};
  double mx2_epKm_{};
  double Theta_e_gamma_{};
  double Theta_e_phimeson_{};
  double Cone_Kp_{};
  double Cone_Km_{};
  double Cone_p_{};
  double coplanarity_had_normals_deg_{};

  double Mass_pi0_{};
  double Mx2_eppi0_{};
  double emiss_pi0_{};
  double Mx2_ep_pi0_{};
  double Mx2_epi0_{};
  double ptmiss_pi0_{};
  double theta_pi0pi0_{};
  double delta_phi_pi0_{};
  double theta_epho1_{};
  double theta_epho2_{};

  double DeltaE_{};

 public:
  DISANAMath() = default;

  // DVCS-style (photon provided as p,theta,phi)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double g_p, double g_theta, double g_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector photon = Build4Vector(g_p, g_theta, g_phi, 0.0);
    ComputeKinematics(e_in, e_out, p_in, p_out, photon);
  }

  // DVPi0-style (photon provided as p,theta,phi)
  DISANAMath(Pi0Tag, double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double g_p, double g_theta, double g_phi,
             double g2_p, double g2_theta, double g2_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector photon = Build4Vector(g_p, g_theta, g_phi, 0.0);
    TLorentzVector photon2 = Build4Vector(g2_p, g2_theta, g2_phi, 0.0);
    ComputePi0Kinematics(e_in, e_out, p_in, p_out, photon, photon2);
  }

  // φ→K+K− analysis (both kaons given)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double kMinus_p, double kMinus_theta,
             double kMinus_phi, double kPlus_p, double kPlus_theta, double kPlus_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector kMinus = Build4Vector(kMinus_p, kMinus_theta, kMinus_phi, m_kMinus);
    TLorentzVector kPlus = Build4Vector(kPlus_p, kPlus_theta, kPlus_phi, m_kPlus);
    ComputeKinematics(e_in, e_out, p_in, p_out, kPlus, kMinus);
  }

  // Mixed: one kaon measured, the other missing (set IsMissing=true if k is the missing one)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double k_p, double k_theta, double k_phi,
             bool IsMissing) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector k_rec = Build4Vector(k_p, k_theta, k_phi, m_kMinus);  // treat as kaon mass
    TLorentzVector k_miss = e_in + p_in - e_out - p_out - k_rec;
    TLorentzVector kMinus, kPlus;
    if (IsMissing) {
      kMinus = k_miss;
      kPlus = k_rec;
    } else {
      kMinus = k_rec;
      kPlus = k_miss;
    }
    ComputeKinematics(e_in, e_out, p_in, p_out, kPlus, kMinus);
  }

  // Accessors
  double GetQ2() const { return Q2_; }
  double GetxB() const { return xB_; }
  double GetT() const { return t_; }
  double GetPhi() const { return phi_deg_; }
  double GetW() const { return W_; }
  double GetNu() const { return nu_; }
  double Gety() const { return y_; }

  double GetMx2_ep() const { return mx2_ep_; }
  double GetEmiss() const { return emiss_; }
  double GetPTmiss() const { return ptmiss_; }
  double GetMx2_epg() const { return mx2_epg_; }
  double GetMx2_epKpKm() const { return mx2_epKpKm_; }
  double GetDeltaPhi() const { return delta_phi_; }
  double GetTheta_gamma_gamma() const { return theta_gg_; }
  double GetMx2_egamma() const { return mx2_egamma_; }
  double GetMx2_eKpKm() const { return mx2_eKpKm_; }
  double GetMx2_epKp() const { return mx2_epKp_; }
  double GetMx2_epKm() const { return mx2_epKm_; }
  double GetTheta_e_gamma() const { return Theta_e_gamma_; }
  double GetTheta_g_phimeson() const { return theta_gphi_; }
  double GetTheta_e_phimeson() const { return Theta_e_phimeson_; }
  double GetDeltaE() const { return DeltaE_; }
  double GetMx_epKp() const { return (mx2_epKp_ > 0) ? std::sqrt(mx2_epKp_) : -999.0; }
  double GetCone_Kp() const { return Cone_Kp_; }
  double GetCone_Km() const { return Cone_Km_; }
  double GetCone_p() const { return Cone_p_; }
  

  double GetMass_pi0() const { return Mass_pi0_; }
  double GetMx2_eppi0() const { return Mx2_eppi0_; }
  double GetEmiss_pi0() const { return emiss_pi0_; }
  double GetMx2_ep_pi0() const { return Mx2_ep_pi0_; }
  double GetMx2_epi0() const { return Mx2_epi0_; }
  double GetPTmiss_pi0() const { return ptmiss_pi0_; }
  double GetTheta_pi0pi0() const { return theta_pi0pi0_; }
  double GetTheta_epho1() const { return theta_epho1_; }
  double GetTheta_epho2() const { return theta_epho2_; }
  double GetDeltaPhi_pi0() const { return delta_phi_pi0_; }


  double GetCoplanarity_had_normals_deg() const { return coplanarity_had_normals_deg_; }

  // Helpers
  double ComputePhiH(const TVector3 &q1_, const TVector3 &k1_, const TVector3 &q2_) const {
    const double t1 = ((q1_.Cross(k1_)).Dot(q2_)) / std::abs((q1_.Cross(k1_)).Dot(q2_));
    const TVector3 t2 = q1_.Cross(k1_);
    const TVector3 t3 = q1_.Cross(q2_);
    const double n2 = t2.Mag(), n3 = t3.Mag();
    const double t4 = (t2.Dot(t3)) / (n2 * n3);
    return t1 * std::acos(t4) * 180. / pi + 180.;
  }

  void ComputePi0Kinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &photon1, const TLorentzVector &photon2) {
    TLorentzVector q = electron_in - electron_out;
    TLorentzVector pi0 = photon1 + photon2;
    Mass_pi0_ = pi0.Mag();
    Mx2_eppi0_ = (electron_in + proton_in - electron_out - proton_out - pi0).Mag2();
    emiss_pi0_ = (electron_in + proton_in - electron_out - proton_out - pi0).E();
    Mx2_ep_pi0_ = (electron_in + proton_in - electron_out - proton_out).Mag2();
    Mx2_epi0_ = (electron_in + proton_in - electron_out - pi0).Mag2();
    ptmiss_pi0_ = (electron_in + proton_in - electron_out - proton_out - pi0).Vect().Perp();
    theta_pi0pi0_ = pi0.Angle((electron_in + proton_in - (electron_out + proton_out)).Vect()) * 180. / pi;
    delta_phi_pi0_ = std::abs(ComputePhiH(q.Vect(), electron_in.Vect(), pi0.Vect()) - ComputePhiH(q.Vect(), electron_in.Vect(), -proton_out.Vect()));
    theta_epho1_ = electron_out.Angle(photon1.Vect()) * 180. / pi;
    theta_epho2_ = electron_out.Angle(photon2.Vect()) * 180. / pi;
  }

  // Core (DVCS-like)
  void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &photon) {
    TLorentzVector q = electron_in - electron_out;

    Q2_ = -q.Mag2();
    nu_ = q.E();
    y_ = nu_ / electron_in.E();
    W_ = (proton_in + q).Mag();
    xB_ = Q2_ / (2.0 * proton_in.Dot(q));
    t_ = std::abs((proton_in - proton_out).Mag2());

    TVector3 nL = electron_in.Vect().Cross(electron_out.Vect()).Unit();
    TVector3 nH = photon.Vect().Cross(q.Vect()).Unit();
    const double cos_phi = nL.Dot(nH);
    const double sin_phi = (nL.Cross(nH)).Dot(q.Vect().Unit());
    phi_deg_ = (std::atan2(sin_phi, cos_phi) + pi) * 180. / pi;

    const TLorentzVector totI = electron_in + proton_in;
    const TLorentzVector totF = electron_out + proton_out + photon;
    const TLorentzVector miss = totI - totF;

    mx2_ep_ = (totI - electron_out - proton_out).Mag2();
    emiss_ = miss.E();
    ptmiss_ = miss.Vect().Perp();
    mx2_epg_ = miss.Mag2();

    // Coplanarity (hadron vs q)
    TVector3 qv = q.Vect();
    TVector3 ev = electron_in.Vect();
    TVector3 gv = photon.Vect();
    TVector3 pv = proton_out.Vect();
    delta_phi_ = std::abs(ComputePhiH(qv, ev, gv) - ComputePhiH(qv, ev, -pv));
    theta_gg_ = photon.Angle((totI - (electron_out + proton_out)).Vect()) * 180. / pi;

    mx2_egamma_ = (electron_in + proton_in - electron_out - photon).Mag2();
    Theta_e_gamma_ = electron_out.Angle(photon.Vect()) * 180. / pi;
    DeltaE_ = (electron_in.E() + proton_in.E()) - (electron_out.E() + proton_out.E() + photon.E());
  }

  // Core (phi)
  void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &kPlus, const TLorentzVector &kMinus) {
    TLorentzVector q = electron_in - electron_out;
    TLorentzVector phi = kPlus + kMinus;

    Q2_ = -q.Mag2();
    nu_ = q.E();
    y_ = nu_ / electron_in.E();
    W_ = (proton_in + q).Mag();
    xB_ = Q2_ / (2.0 * proton_in.Dot(q));
    t_ = std::abs((proton_in - proton_out).Mag2());

    TVector3 nL = electron_in.Vect().Cross(electron_out.Vect()).Unit();
    TVector3 nH = q.Vect().Cross(proton_out.Vect()).Unit();
    const double cos_phi = nL.Dot(nH);
    const double sin_phi = (nL.Cross(nH)).Dot(q.Vect().Unit());
    phi_deg_ = (std::atan2(sin_phi, cos_phi) + pi) * 180. / pi;

    const TLorentzVector totI = electron_in + proton_in;
    const TLorentzVector totF = electron_out + proton_out + phi;
    const TLorentzVector miss = totI - totF;

    mx2_ep_ = (totI - electron_out - proton_out).Mag2();
    emiss_ = miss.E();
    ptmiss_ = miss.Vect().Perp();
    mx2_epKpKm_ = miss.Mag2();

    TVector3 qv = q.Vect();
    TVector3 ev = electron_in.Vect();
    TVector3 phiv = phi.Vect();
    TVector3 pv = proton_out.Vect();
    TVector3 kpv = kPlus.Vect();
    TVector3 kmv = kMinus.Vect();
    TVector3 n_qphi = qv.Cross(phiv);
    TVector3 n_pphi = pv.Cross(phiv);

    delta_phi_ = std::abs(ComputePhiH(qv, ev, phiv) - ComputePhiH(qv, ev, -pv));
    theta_gphi_ = phiv.Angle((totI - (electron_out + proton_out)).Vect()) * 180. / pi;
    Cone_p_ = (electron_in + proton_in - electron_out - phi).Angle(pv) * 180. / pi;  // p vs K+
    Cone_Km_ = (electron_in + proton_in - electron_out - kPlus - proton_out).Angle(kmv) * 180. / pi;  // p vs K−
    Cone_Kp_ = (electron_in + proton_in - electron_out - kMinus - proton_out).Angle(kpv) * 180. / pi;  // p vs K−

    mx2_eKpKm_ = (electron_in + proton_in - electron_out - phi).Mag2();
    mx2_epKp_ = (electron_in + proton_in - electron_out - kPlus - proton_out).Mag2();
    mx2_epKm_ = (electron_in + proton_in - electron_out - kMinus - proton_out).Mag2();

    Theta_e_phimeson_ = electron_out.Angle(phi.Vect()) * 180. / pi;
    DeltaE_ = (electron_in.E() + proton_in.E()) - (electron_out.E() + proton_out.E() + phi.E());
    double coplanarity_had_normals_deg_ = std::numeric_limits<double>::quiet_NaN();
    if (n_qphi.Mag() > 0 && n_pphi.Mag() > 0) {
      n_qphi = n_qphi.Unit();
      n_pphi = n_pphi.Unit();

      double c = n_qphi.Dot(n_pphi);
      // numerical safety
      c = std::max(-1.0, std::min(1.0, c));
      double ang = std::acos(c) * 180. / pi; // in [0, 180]

      // treat ±normals as equivalent (optional but common)
      coplanarity_had_normals_deg_ = std::min(ang, 180.0 - ang);
    }
  }

  // ---------------------------------------------------------------------------
  // dσ/dφ (DVCS) and helpers (unchanged except includes)
  // ---------------------------------------------------------------------------
  std::vector<std::vector<std::vector<TH1D *>>> ComputeDVCS_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
    TStopwatch timer;
    timer.Start();
    constexpr double phi_min = 0.0, phi_max = 360.0;
    constexpr int n_phi_bins = 18;

    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const auto &xb_bins = bins.GetXBBins();

    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t = t_bins.size() - 1;
    const size_t n_xb = xb_bins.size() - 1;

    std::vector<std::vector<std::vector<TH1D *>>> hist(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    std::vector<std::vector<std::vector<double>>> q2xBtbins(n_xb, std::vector<std::vector<double>>(n_q2, std::vector<double>(n_t, 0.0)));

    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
          const double tmin = t_bins[it], tmax = t_bins[it + 1];
          const double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];
          auto name = Form("hphi_q%.1f_t%.1f_xb%.2f", qmin, tmin, xbmin);
          auto title = Form("d#sigma/d#phi (Q^{2}=[%.1f,%.1f], t=[%.1f,%.1f], x_{B}=[%.2f,%.2f])", qmin, qmax, tmin, tmax, xbmin, xbmax);
          hist[ix][iq][it] = new TH1D(name, title, n_phi_bins, phi_min, phi_max);
          hist[ix][iq][it]->SetDirectory(nullptr);
          q2xBtbins[ix][iq][it] = (qmax - qmin) * (tmax - tmin) * (xbmax - xbmin);
        }

    auto findBin = [](double v, const std::vector<double> &e) -> int {
      auto it = std::upper_bound(e.begin(), e.end(), v);
      if (it == e.begin() || it == e.end()) return -1;
      return int(it - e.begin()) - 1;
    };

    auto fill = [&](double Q2, double t, double xB, double phi) {
      int iq = findBin(Q2, q2_bins), it = findBin(t, t_bins), ix = findBin(xB, xb_bins);
      if (iq >= 0 && it >= 0 && ix >= 0) hist[ix][iq][it]->Fill(phi);
    };

    df.Foreach(fill, {"Q2", "t", "xB", "phi"});

    const double bin_width = (pi / 180.) * (phi_max - phi_min) / n_phi_bins;
    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          TH1D *h = hist[ix][iq][it];
          if (!h) continue;
          const double vol = q2xBtbins[ix][iq][it];
          for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double N = h->GetBinContent(b);
            const double ds = N / (luminosity * bin_width * vol);
            const double es = std::sqrt(std::max(0.0, N)) / (luminosity * bin_width * vol);
            h->SetBinContent(b, ds);
            h->SetBinError(b, es);
          }
        }

    std::cout << "DVCS cross-sections computed in a single pass.\n";
    timer.Stop();
    std::cout << "Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return hist;
  }

  // φ dσ/dt (unchanged)
  std::vector<TH1D *> ComputePhi_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
    TStopwatch timer;
    timer.Start();
    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const size_t n_q2 = q2_bins.size() - 1, n_t = t_bins.size() - 1;
    std::vector<TH1D *> hist(n_q2, nullptr);

    for (size_t iq = 0; iq < n_q2; ++iq) {
      const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
      auto hname = Form("phi_cs_q2bin%zu", iq);
      auto htitle = Form("d#sigma/dt for Q^{2}=[%.2f, %.2f]", qmin, qmax);
      hist[iq] = new TH1D(hname, htitle, n_t, &t_bins[0]);
      hist[iq]->SetDirectory(nullptr);
      hist[iq]->GetXaxis()->SetTitle("-t [GeV^{2}]");
      hist[iq]->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    }

    auto findBin = [](double v, const std::vector<double> &e) -> int {
      auto it = std::upper_bound(e.begin(), e.end(), v);
      if (it == e.begin() || it == e.end()) return -1;
      return int(it - e.begin()) - 1;
    };
    auto fill = [&](double Q2, double t) {
      int iq = findBin(Q2, q2_bins), it = findBin(t, t_bins);
      if (iq >= 0 && it >= 0) hist[iq]->Fill(t);
    };
    df.Foreach(fill, {"Q2", "t"});

    for (auto *h : hist) {
      if (!h) continue;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double w = h->GetBinWidth(b), N = h->GetBinContent(b);
        const double ds = N / (luminosity * w), es = std::sqrt(std::max(0.0, N)) / (luminosity * w);
        h->SetBinContent(b, ds);
        h->SetBinError(b, es);
      }
    }
    timer.Stop();
    std::cout << "[ComputePhi_CrossSection] Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return hist;
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

  std::vector<std::vector<std::vector<TH1D *>>> CalcPi0Corr(ROOT::RDF::RNode df_dvcs_pi0mc, ROOT::RDF::RNode df_pi0_pi0mc, ROOT::RDF::RNode df_dvcs_data,
                                                            ROOT::RDF::RNode df_pi0_data, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath pi0Corr;
    auto df_dvcs_pi0mc_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_dvcs_pi0mc, xBins, 1);
    auto df_pi0_pi0mc_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_pi0_pi0mc, xBins, 1);
    auto df_dvcs_data_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_dvcs_data, xBins, 1);
    auto df_pi0_data_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_pi0_data, xBins, 1);

    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcs_pi0mc = df_dvcs_pi0mc_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_pi0_pi0mc = df_pi0_pi0mc_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcs_data = df_dvcs_data_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_pi0_data = df_pi0_data_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_dvcs_pi0mc || !h_pi0_pi0mc || !h_dvcs_data || !h_pi0_data) {
            std::cerr << "Pi0 corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_dvcs_pi0mc->Clone(Form("hPi0Corr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_dvcs_pi0mc, h_pi0_pi0mc);
          hRatio->Multiply(hRatio, h_pi0_data);
          hRatio->Divide(hRatio, h_dvcs_data);
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcAcceptanceCorr(ROOT::RDF::RNode df_gen_dvcsmc, ROOT::RDF::RNode df_accept_dvcsmc, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath AccCorr;
    auto df_gen_dvcsmc_CrossSection = AccCorr.ComputeDVCS_CrossSection(df_gen_dvcsmc, xBins, 1);
    auto df_accept_dvcsmc_CrossSection = AccCorr.ComputeDVCS_CrossSection(df_accept_dvcsmc, xBins, 1);

    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_gen_dvcsmc = df_gen_dvcsmc_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_accept_dvcsmc = df_accept_dvcsmc_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_gen_dvcsmc || !h_accept_dvcsmc) {
            std::cerr << "Acceptance Corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_gen_dvcsmc->Clone(Form("hAccCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_accept_dvcsmc, h_gen_dvcsmc);
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcEfficiencyCorr(ROOT::RDF::RNode df_dvcsmc_bkg, ROOT::RDF::RNode df_dvcsmc_nobkg, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath EffCorr;
    auto df_dvcsmc_bkg_CrossSection = EffCorr.ComputeDVCS_CrossSection(df_dvcsmc_bkg, xBins, 1);
    auto df_dvcsmc_nobkg_CrossSection = EffCorr.ComputeDVCS_CrossSection(df_dvcsmc_nobkg, xBins, 1);
    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcsmc_bkg = df_dvcsmc_bkg_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcsmc_nobkg = df_dvcsmc_nobkg_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_dvcsmc_bkg || !h_dvcsmc_nobkg) {
            std::cerr << "Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_dvcsmc_bkg->Clone(Form("hEffCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_dvcsmc_bkg, h_dvcsmc_nobkg);
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcRadiativeCorr(ROOT::RDF::RNode df_dvcs_rad, ROOT::RDF::RNode df_dvcs_norad, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath RadCorr;
    auto df_dvcs_rad_CrossSection = RadCorr.ComputeDVCS_CrossSection(df_dvcs_rad, xBins, 1);
    auto df_dvcs_norad_CrossSection = RadCorr.ComputeDVCS_CrossSection(df_dvcs_norad, xBins, 1);
    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcs_rad = df_dvcs_rad_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcs_norad = df_dvcs_norad_CrossSection[xb_bin][q2_bin][t_bin];
          if (!h_dvcs_rad || !h_dvcs_norad) {
            std::cerr << " Rad Corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
         TH1D *hRatio = static_cast<TH1D *>(h_dvcs_rad->Clone(Form("hRadCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
         hRatio->SetDirectory(nullptr); // avoid ROOT ownership surprises
         hRatio->Divide(h_dvcs_norad); // uses default error propagation
         hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }
};

#endif  // DISANAMATH_H

#ifndef DISANAMATH_H
#define DISANAMATH_H
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"

constexpr double pi = 3.14159265358979323846;
const double m_e = 0.000511;  // GeV
const double m_p = 0.938272;  // GeV

namespace {
TVector3 SphericalToCartesian(double p, double theta, double phi) {
  double px = p * std::sin(theta) * std::cos(phi);
  double py = p * std::sin(theta) * std::sin(phi);
  double pz = p * std::cos(theta);
  return TVector3(px, py, pz);
}

TLorentzVector Build4Vector(double p, double theta, double phi, double mass) {
  TVector3 vec = SphericalToCartesian(p, theta, phi);
  double E = std::sqrt(vec.Mag2() + mass * mass);
  return TLorentzVector(vec, E);
}
} // anonymous namespace
// require for managing bins in the analysis
class BinManager {
public:
    BinManager(){
    q2_bins_ = {1.0, 2.0, 4.0, 6.0};
    t_bins_  = {0.1, 0.3, 0.6, 1.0};
    xb_bins_ = {0.1, 0.2, 0.4, 0.6};
 }

    // Getters
    const std::vector<double>& GetQ2Bins() const { return q2_bins_; }
    const std::vector<double>& GetTBins() const { return t_bins_; }
    const std::vector<double>& GetXBBins() const { return xb_bins_; }

    // Setters
    void SetQ2Bins(const std::vector<double>& bins) { q2_bins_ = bins; }
    void SetTBins(const std::vector<double>& bins)  { t_bins_ = bins; }
    void SetXBBins(const std::vector<double>& bins) { xb_bins_ = bins; }

private:
    std::vector<double> q2_bins_;
    std::vector<double> t_bins_;
    std::vector<double> xb_bins_;
};
/// 
class DISANAMath {
  private:
  // Internal kinematic variables
  double Q2_;
  double xB_;
  double t_;
  double phi_deg_;
  double W_;
  double nu_;
  double y_;

  // Internal compute function
 // void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out, const TLorentzVector &photon);

 public:
  DISANAMath() = default;
  // Constructor: takes momenta and computes kinematics
  //DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double g_p, double g_theta, double g_phi);
      // Getter functions
  
  double GetQ2() const {
    return Q2_;
  }
  double GetxB() const { return xB_; }
  double GetT() const { return t_; }
  double GetPhi() const { return phi_deg_; }
  double GetW() const { return W_; }
  double GetNu() const { return nu_; }
  double Gety() const { return y_; }

  //std::vector<TH1D *> ComputeDVCS_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity);
 // std::vector<std::vector<std::vector<TH1D*>>> ComputeDVCS_CrossSection2(ROOT::RDF::RNode df, const BinManager &bins, double luminosity);
 // std::vector<TH1D *> ComputeBeamSpinAsymmetry(const std::vector<TH1D *> &sigma_pos,const std::vector<TH1D *> &sigma_neg);


DISANAMath(double e_in_E,
                       double e_out_p, double e_out_theta, double e_out_phi,
                       double p_out_p, double p_out_theta, double p_out_phi,
                       double g_p, double g_theta, double g_phi) {
  TLorentzVector electron_in(0, 0, e_in_E, e_in_E);
  TLorentzVector electron_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
  TLorentzVector proton_in(0, 0, 0, m_p);
  TLorentzVector proton_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
  TLorentzVector photon = Build4Vector(g_p, g_theta, g_phi, 0.0);
  ComputeKinematics(electron_in, electron_out, proton_in, proton_out, photon);
}

void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out,
                                   const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                                   const TLorentzVector &photon) {
  TLorentzVector q = electron_in - electron_out;
  Q2_ = -q.Mag2();
  nu_ = q.E();
  y_ = nu_ / electron_in.E();
  W_ = (proton_in + q).Mag();
  xB_ = Q2_ / (2.0 * proton_in.Dot(q));
  t_ = std::abs((proton_in - proton_out).Mag2());

  TVector3 n_L = electron_in.Vect().Cross(electron_out.Vect()).Unit();
  TVector3 n_H = q.Vect().Cross(proton_out.Vect()).Unit();

  double cos_phi = n_L.Dot(n_H);
  double sin_phi = (n_L.Cross(n_H)).Dot(q.Vect().Unit());

  double phi = std::atan2(sin_phi, cos_phi);
  if (phi < 0.0) {
    phi += 2.0*TMath::Pi(); // Ensure phi is in [0, 360] degrees
  }
  phi = phi * 180.0 / TMath::Pi();

  phi_deg_ = phi;
}

// Integrated luminosity in cm⁻² (example, you can scale it out if not known)

std::vector<TH1D *> ComputeDVCS_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
  constexpr double phi_min = 0.0;
  constexpr double phi_max = 360.0;
  constexpr int n_phi_bins = 12;

  const auto &q2_bins = bins.GetQ2Bins();
  const auto &t_bins = bins.GetTBins();
  const auto &xb_bins = bins.GetXBBins();

  std::vector<TH1D *> histograms;

  for (size_t iq = 0; iq < q2_bins.size() - 1; ++iq) {
    for (size_t it = 0; it < t_bins.size() - 1; ++it) {
      for (size_t ix = 0; ix < xb_bins.size() - 1; ++ix) {
        double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
        double tmin = t_bins[it], tmax = t_bins[it + 1];
        double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];

        std::string cut = Form(
            "(Q2 > %.3f && Q2 <= %.3f) && (t > %.3f && t <= "
            "%.3f) && (xB > %.3f && xB <= %.3f)",
            qmin, qmax, tmin, tmax, xbmin, xbmax);

       // df.Display({"Q2", "t", "xB"})->Print();
        auto df_bin = df.Filter(cut);
        auto nEntries = df_bin.Count();
        //std::cout << "Entries in bin (Q2=[" << qmin << "," << qmax << "], t=[" << tmin << "," << tmax << "], xB=[" << xbmin << "," << xbmax << "]): " << *nEntries << std::endl;

        std::string histName = Form("hphi_q%.1f_t%.1f_xb%.2f", qmin, tmin, xbmin);
        std::string histTitle = Form("d#sigma/d#phi (Q^{2}=[%.1f,%.1f], t=[%.1f,%.1f], xB=[%.2f,%.2f])", qmin, qmax, tmin, tmax, xbmin, xbmax);

        auto hphi = df_bin.Histo1D({histName.c_str(), histTitle.c_str(), n_phi_bins, phi_min, phi_max}, "phi");
        hphi->Draw(); // materialize
        //std::cout << " hphi values are obtained here" << std::endl;
        // Clone histogram from RDataFrame-managed memory
        TH1D *hclone = dynamic_cast<TH1D *>(hphi->Clone(histName.c_str()));
        double bin_width = (phi_max - phi_min) / n_phi_bins;

        for (int b = 1; b <= hclone->GetNbinsX(); ++b) {
          double raw = hclone->GetBinContent(b);
          double norm = raw / (luminosity * bin_width);
          double err = std::sqrt(raw) / (luminosity * bin_width);
          hclone->SetBinContent(b, norm);
          hclone->SetBinError(b, err);
        }

        histograms.push_back(hclone);
      }
    }
  }

  std::cout << "DVCS cross-sections computed and histograms returned.\n";
  return histograms;
}

std::vector<std::vector<std::vector<TH1D*>>> ComputeDVCS_CrossSection2(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
  constexpr double phi_min = 0.0;
  constexpr double phi_max = 360.0;
  constexpr int n_phi_bins = 12;

  const auto &q2_bins = bins.GetQ2Bins();
  const auto &t_bins = bins.GetTBins();
  const auto &xb_bins = bins.GetXBBins();

  size_t n_q2 = q2_bins.size() - 1;
  size_t n_t  = t_bins.size() - 1;
  size_t n_xb = xb_bins.size() - 1;

  // Allocate 3D vector: x_B → Q² → t → histogram
  std::vector<std::vector<std::vector<TH1D*>>> histograms(
      n_xb,
      std::vector<std::vector<TH1D*>>(n_q2, std::vector<TH1D*>(n_t, nullptr))
  );

  for (size_t iq = 0; iq < n_q2; ++iq) {
      for (size_t it = 0; it < n_t; ++it) {
          for (size_t ix = 0; ix < n_xb; ++ix) {
              double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
              double tmin = t_bins[it], tmax = t_bins[it + 1];
              double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];

              std::string cut = Form(
                  "(Q2 > %.3f && Q2 <= %.3f) && (t > %.3f && t <= %.3f) && (xB > %.3f && xB <= %.3f)",
                  qmin, qmax, tmin, tmax, xbmin, xbmax);

              auto df_bin = df.Filter(cut);

              std::string histName = Form("hphi_q%.1f_t%.1f_xb%.2f", qmin, tmin, xbmin);
              std::string histTitle = Form("d#sigma/d#phi (Q^{2}=[%.1f,%.1f], t=[%.1f,%.1f], xB=[%.2f,%.2f])",
                                           qmin, qmax, tmin, tmax, xbmin, xbmax);

              auto hphi = df_bin.Histo1D({histName.c_str(), " ", n_phi_bins, phi_min, phi_max}, "phi");
              hphi->Draw(); // Force evaluation

              TH1D *hclone = dynamic_cast<TH1D *>(hphi->Clone(histName.c_str()));
              hclone->SetDirectory(0);

              double bin_width = (phi_max - phi_min) / n_phi_bins;

              for (int b = 1; b <= hclone->GetNbinsX(); ++b) {
                  double raw = hclone->GetBinContent(b);
                  double norm = raw / (luminosity * bin_width);
                  double err = std::sqrt(raw) / (luminosity * bin_width);
                  hclone->SetBinContent(b, norm);
                  hclone->SetBinError(b, err);
              }

              // Store in 3D structure
              histograms[ix][iq][it] = hclone;
          }
      }
  }

  std::cout << "DVCS cross-sections computed in 3D bin structure.\n";
  return histograms;
}
/// BSA studeis 
std::vector<TH1D *> ComputeBeamSpinAsymmetry(const std::vector<TH1D *> &sigma_pos,
  const std::vector<TH1D *> &sigma_neg) {
std::vector<TH1D *> asymmetry_histograms;

if (sigma_pos.size() != sigma_neg.size()) {
std::cerr << "ERROR: Helicity + and - histogram vectors must have the same size!\n";
return asymmetry_histograms;
}

for (size_t i = 0; i < sigma_pos.size(); ++i) {
TH1D *h_plus = sigma_pos[i];
TH1D *h_minus = sigma_neg[i];

if (!h_plus || !h_minus) {
std::cerr << "WARNING: Null histogram at index " << i << "\n";
continue;
}

std::string asym_name = std::string(h_plus->GetName()) + "_BSA";
std::string asym_title = std::string("Beam Spin Asymmetry of ") + h_plus->GetTitle();

TH1D *h_asym = (TH1D *)h_plus->Clone(asym_name.c_str());
h_asym->SetTitle(" ");
h_asym->Reset();

for (int bin = 1; bin <= h_plus->GetNbinsX(); ++bin) {
double N_plus = h_plus->GetBinContent(bin);
double N_minus = h_minus->GetBinContent(bin);
double E_plus = h_plus->GetBinError(bin);
double E_minus = h_minus->GetBinError(bin);

double denom = N_plus + N_minus;
double num = N_plus - N_minus;

double asym = (denom != 0) ? num / denom : 0.0;
double err = (denom != 0) ?
2.0 / (denom * denom) *
std::sqrt(std::pow(N_minus * E_plus, 2) + std::pow(N_plus * E_minus, 2)) : 0.0;

h_asym->SetBinContent(bin, asym);
h_asym->SetBinError(bin, err);
}

asymmetry_histograms.push_back(h_asym);
}

std::cout << "Beam Spin Asymmetries computed.\n";
return asymmetry_histograms;
}
};
#endif  // DISANAMATH_H

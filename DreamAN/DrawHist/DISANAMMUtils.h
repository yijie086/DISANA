#include <cmath>
#include <memory>

#include "DISANAMath.h"  // <-- REQUIRED: provides DISANAMath class & methods
#include "DISANAMathFitUtils.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLorentzVector.h"

using ROOT::VecOps::RVec;
// --- constants for masses
static constexpr double kMe = 0.000511;
static constexpr double kMp = 0.938272;
static constexpr double kMK = 0.493677;
static constexpr double alpha = 1.0 / 137.035999084;

// convenience 3-vector helpers
static double MomentumFunc(float px, float py, float pz) { return std::sqrt(px * px + py * py + pz * pz); }
static double ThetaFunc(float px, float py, float pz) { return std::acos(pz / std::sqrt(px * px + py * py + pz * pz)); }
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}

inline double HandGammaV(double E, double Q2, double W) {

  if (!(E > 0.0) || !(Q2 > 0.0) || !(W > 0.0)) return std::numeric_limits<double>::quiet_NaN();

  const double W2 = W * W;

  // nu = (W^2 + Q^2 - Mp^2) / (2 Mp)
  const double nu = (W2 + Q2 - kMp*kMp) / (2.0 * kMp);

  const double y  = nu / E;
  const double Ep = E - nu;
  if (!(Ep > 0.0)) return std::numeric_limits<double>::quiet_NaN();

  // K = (W^2 - Mp^2) / (2 Mp)
  const double K = (W2 - kMp*kMp) / (2.0 * kMp);

  const double eps_num = 1.0 - y - Q2 / (4.0 * E * E);
  const double eps_den = 1.0 - y + 0.5 * y * y + Q2 / (4.0 * E * E);
  if (eps_den == 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double eps = eps_num / eps_den;

  if ((1.0 - eps) == 0.0) return std::numeric_limits<double>::quiet_NaN();

  // Hand convention
  const double Gamma =
      (alpha / (2.0 * pi * pi)) * (Ep / E) * (K / Q2) * (1.0 / (1.0 - eps));

  return Gamma;
}

inline int ChooseBestElectron(const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz,
                              const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<int>& pass) {
  int best = -1;
  float bestP = -1.f;
  // Prefer electrons in the target window; choose highest |p|
  for (size_t i = 0; i < pid.size(); ++i) {
    if (pid[i] != 11 || pass[i] == 0) continue;
    if (vz[i] < -10.f || vz[i] > 3.f) continue;
    float P = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
    if (P > bestP) {
      bestP = P;
      best = (int)i;
    }
  }
  if (best >= 0) return best;
  // Fallback: highest |p| among passing tracks
  for (size_t i = 0; i < pid.size(); ++i) {
    if (pid[i] != 11 || pass[i] == 0) continue;
    float P = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
    if (P > bestP) {
      bestP = P;
      best = (int)i;
    }
  }
  return best;  // -1 if no usable e⁻
}
inline int DetRegionFromStatus(short s) {
  int a = std::abs(s);
  if (a >= 1000 && a < 2000) return 0;  // FT
  if (a >= 2000 && a < 3000) return 1;  // FD
  if (a >= 4000 && a < 5000) return 2;  // CD
  return -1;
}
// phi event selection cuts for single photon contaminations
ROOT::RDF::RNode SelectExclusivePhiEvent(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass, const ROOT::VecOps::RVec<int>& daughterPass) {
        int e = 0, km = 0, kp = 0, p = 0;
        bool hasPhiDaughter = true;

        for (size_t i = 0; i < pid.size(); ++i) {
          if (pass[i] == 0) continue;

          if (pid[i] == 11) {
            e++;
          } else if (pid[i] == 321) {
            kp++;
            hasPhiDaughter = hasPhiDaughter /*|| daughterPass[i]/*/;  // check if kaon is from phi
          } else if (pid[i] == -321) {
            km++;
            hasPhiDaughter = hasPhiDaughter /*|| daughterPass[i]*/;  // check if kaon is from phi
          } else if (pid[i] == 2212) {
            p++;
          }
        }

        return (e == 1 && kp == 1 && km == 1 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: 1 e⁻, ≥1 K⁺, ≥1 K⁻, 1 proton, ≥1 kaon from phi");
}

ROOT::RDF::RNode SelectPhiEvent_MissingKm(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
        int e = 0, kp = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (pass[i] == 0) continue;
          if (pid[i] == 11)
            ++e;
          else if (pid[i] == 321)
            ++kp;
          else if (pid[i] == 2212)
            ++p;
        }
        return (e >= 1 && kp >= 1 && p >= 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass"}, "Cut: 1 e⁻, ≥1 K⁺, 1 p (Missing K⁻ workflow)");
}

ROOT::RDF::RNode SelectPhiEvent_MissingKp(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
        int e = 0, km = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (pass[i] == 0) continue;
          if (pid[i] == 11)
            ++e;
          else if (pid[i] == -321)
            ++km;
          else if (pid[i] == 2212)
            ++p;
        }
        return (e >= 1 && km >= 1 && p >= 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass"}, "Cut: 1 e⁻, ≥1 K⁺, 1 p (Missing K+ workflow)");
}
//
ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass, const ROOT::VecOps::RVec<int>& daughterPass) {
        int e = 0, g = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (pass[i] == 0) continue;
          if (daughterPass[i]) return false;  // reject if any daughter particle is a pi0
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22)
            g++;  // photon must NOT be from pi0
          else if (pid[i] == 2212)
            p++;
        }
        return (e == 1 && g == 1 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p");
}
//
template <typename Method>
ROOT::RDF::RNode define_DISCAT(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Define(name,
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi, double reckMinus_p,
                                           double reckMinus_theta, double reckMinus_phi, double reckPlus_p, double reckPlus_theta, double reckPlus_phi) {
                       return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, reckMinus_p, reckMinus_theta, reckMinus_phi, reckPlus_p,
                                          reckPlus_theta, reckPlus_phi).*
                               method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "reckMinus_p", "reckMinus_theta", "reckMinus_phi", "reckPlus_p",
                      "reckPlus_theta", "reckPlus_phi"});
}
// -----------------------------------------------------------------------------
// InitKinematics_MissingKm : K⁻ omitted (exclusive K⁺ channel)
// -----------------------------------------------------------------------------
ROOT::RDF::RNode InitKinematics_MissingKm(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  // pick best e, p, K⁺ (your style)
  *df_ = df_->Define("ele_px_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recel_vz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             ////
             .Define("reckMinus_vz",  // dummy for K minus vz plotting
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("reckPlus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kPlus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kPlus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kPlus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})

             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})

             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
             .Define("RunNumber", "RUN_config_run")
             .Define("EventNumber", "RUN_config_event")
             .Define("nElectrons",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
                       int n = 0;
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) ++n;
                       return n;
                     },
                     {"REC_Particle_pid", "REC_Particle_pass"})
             // Choose ONE best electron per event
             .Define("bestEle_idx", ChooseBestElectron, {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz", "REC_Particle_vz", "REC_Particle_pass"})

             // If you want to **keep** multi-e events but only use the best e⁻, leave the next line.
             // If you want to **drop** multi-e events (strict DVEP), add: .Filter("nElectrons==1","Cut: exactly 1 e⁻")
             .Filter("bestEle_idx >= 0", "At least one usable e⁻")
             .Define("REC_pass_bestE",
                     [](int best, const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
                       ROOT::VecOps::RVec<int> keep(pass.size(), 0);
                       if (best >= 0) keep[best] = 1;        // only best e⁻
                       for (size_t i = 0; i < pid.size(); ++i)  // keep non-e tracks as they were
                         if (pid[i] != 11) keep[i] = pass[i];
                       return keep;
                     },
                     {"bestEle_idx", "REC_Particle_pid", "REC_Particle_pass"})

             .Define("nElectrons_best",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& keep) {
                       int n = 0;
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && keep[i]) ++n;
                       return n;
                     },
                     {"REC_Particle_pid", "REC_pass_bestE"})
             // Project the chosen e⁻ to the scalar columns used everywhere else
             .Define("ele_px", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_px"})
             .Define("ele_py", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_py"})
             .Define("ele_pz", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_pz"})
             .Define("recel_vz", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_vz"})
             .Define("ele_det_region", [](int i, const ROOT::VecOps::RVec<short>& st) { return i >= 0 ? DetRegionFromStatus(st[i]) : -1; }, {"bestEle_idx", "REC_Particle_status"});

  // missing K⁻ 4-vector (components); keep both px/py/pz and derived p,θ,φ
  *df_ = df_->Define("kMinus_miss_px",
                     [beam_energy](float epx, float epy, float epz, float ppx, float ppy, float ppz, float kpx, float kpy, float kpz) -> float {
                       TLorentzVector pBeam(0, 0, beam_energy, beam_energy), pTarg(0, 0, 0, kMp);
                       TLorentzVector e(epx, epy, epz, std::sqrt(epx * epx + epy * epy + epz * epz + kMe * kMe));
                       TLorentzVector p(ppx, ppy, ppz, std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz + kMp * kMp));
                       TLorentzVector kp(kpx, kpy, kpz, std::sqrt(kpx * kpx + kpy * kpy + kpz * kpz + kMK * kMK));
                       return float((pBeam + pTarg - e - p - kp).Px());
                     },
                     {"ele_px", "ele_py", "ele_pz", "pro_px", "pro_py", "pro_pz", "kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("kMinus_miss_py",
                     [beam_energy](float epx, float epy, float epz, float ppx, float ppy, float ppz, float kpx, float kpy, float kpz) -> float {
                       TLorentzVector pBeam(0, 0, beam_energy, beam_energy), pTarg(0, 0, 0, kMp);
                       TLorentzVector e(epx, epy, epz, std::sqrt(epx * epx + epy * epy + epz * epz + kMe * kMe));
                       TLorentzVector p(ppx, ppy, ppz, std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz + kMp * kMp));
                       TLorentzVector kp(kpx, kpy, kpz, std::sqrt(kpx * kpx + kpy * kpy + kpz * kpz + kMK * kMK));
                       return float((pBeam + pTarg - e - p - kp).Py());
                     },
                     {"ele_px", "ele_py", "ele_pz", "pro_px", "pro_py", "pro_pz", "kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("kMinus_miss_pz",
                     [beam_energy](float epx, float epy, float epz, float ppx, float ppy, float ppz, float kpx, float kpy, float kpz) -> float {
                       TLorentzVector pBeam(0, 0, beam_energy, beam_energy), pTarg(0, 0, 0, kMp);
                       TLorentzVector e(epx, epy, epz, std::sqrt(epx * epx + epy * epy + epz * epz + kMe * kMe));
                       TLorentzVector p(ppx, ppy, ppz, std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz + kMp * kMp));
                       TLorentzVector kp(kpx, kpy, kpz, std::sqrt(kpx * kpx + kpy * kpy + kpz * kpz + kMK * kMK));
                       return float((pBeam + pTarg - e - p - kp).Pz());
                     },
                     {"ele_px", "ele_py", "ele_pz", "pro_px", "pro_py", "pro_pz", "kPlus_px", "kPlus_py", "kPlus_pz"})

             .Filter([](float ex, float kPlusx, float px) { return ex != -999 && kPlusx != -999 && px != -999; }, {"ele_px", "kPlus_px", "pro_px"})
             // derived angles/magnitudes
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("reckPlus_p", MomentumFunc, {"kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("reckPlus_theta", ThetaFunc, {"kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("reckPlus_phi", PhiFunc, {"kPlus_px", "kPlus_py"})
             .Define("reckMinus_p", MomentumFunc, {"kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz"})
             .Define("reckMinus_theta", ThetaFunc, {"kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz"})
             .Define("reckMinus_phi", PhiFunc, {"kMinus_miss_px", "kMinus_miss_py"})
             .Define("kMinus_det_region",  // dummy for K minus det region plotting
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 321 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("kPlus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 321 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("ele_det_region_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("invMass_KpKm",
                     [](float px1, float py1, float pz1, float px2, float py2, float pz2) -> float {
                       constexpr float mK = 0.493677;  // Kaon mass in GeV/c²
                       float E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mK * mK);
                       float E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mK * mK);
                       float px = px1 + px2;
                       float py = py1 + py2;
                       float pz = pz1 + pz2;
                       float E = E1 + E2;
                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"kPlus_px", "kPlus_py", "kPlus_pz", "kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz"})
             // p–K⁻ invariant mass
             .Define("invMass_pKminus",
                     [](float p_px, float p_py, float p_pz, float k_px, float k_py, float k_pz) -> float {
                       constexpr float mP = 0.938272f;  // Proton mass in GeV/c²
                       constexpr float mK = 0.493677f;  // Kaon mass in GeV/c²

                       float E_p = std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz + mP * mP);
                       float E_k = std::sqrt(k_px * k_px + k_py * k_py + k_pz * k_pz + mK * mK);

                       float px = p_px + k_px;
                       float py = p_py + k_py;
                       float pz = p_pz + k_pz;
                       float E = E_p + E_k;

                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"pro_px", "pro_py", "pro_pz", "kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz"})

             // p–K⁺ invariant mass
             .Define("invMass_pKplus",
                     [](float p_px, float p_py, float p_pz, float k_px, float k_py, float k_pz) -> float {
                       constexpr float mP = 0.938272f;  // Proton mass in GeV/c²
                       constexpr float mK = 0.493677f;  // Kaon mass in GeV/c²

                       float E_p = std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz + mP * mP);
                       float E_k = std::sqrt(k_px * k_px + k_py * k_py + k_pz * k_pz + mK * mK);

                       float px = p_px + k_px;
                       float py = p_py + k_py;
                       float pz = p_pz + k_pz;
                       float E = E_p + E_k;

                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"pro_px", "pro_py", "pro_pz", "kPlus_px", "kPlus_py", "kPlus_pz"});
  // φ mass built from missing K+ and measured K-
  // DISANAMath-driven observables
  *df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "tmin", &DISANAMath::GetTmin, beam_energy);
  *df_ = df_->Define("mtprime",  // non-negative: this is what you’ll plot
                     [](double mt, double tmin) { return std::abs(mt + tmin); }, {"t", "tmin"})
             .Define("tprime",  // optional signed t' ≤ 0
                     [](double mtp) { return -mtp; }, {"mtprime"});
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);
  *df_ = define_DISCAT(*df_, "cos_thetaKK", &DISANAMath::GetCosTheta_KK, beam_energy);
  *df_ = define_DISCAT(*df_, "cos_phiKK", &DISANAMath::GetCosPhi_KK, beam_energy);

  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKp", &DISANAMath::GetMx2_eKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx_eKp", &DISANAMath::GetMx_eKp, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_p", &DISANAMath::GetCone_p, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Kp", &DISANAMath::GetCone_Kp, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Km", &DISANAMath::GetCone_Km, beam_energy);
  *df_ = define_DISCAT(*df_, "z_phi", &DISANAMath::GetZ_phi, beam_energy);
  *df_ = define_DISCAT(*df_, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  *df_ = df_->Define("Gamma_v", [beam_energy](double Q2, double W) {
                 return HandGammaV(beam_energy, Q2, W);
               }, {"Q2", "W"});
  // --- new: vm_cut and Egamma_star -----------------------------------------
  *df_ = df_->Define("vm_cut",
    [beam_energy](double Q2, double xB, double t) -> double {
      return DISANAMath::ComputeVmApprox(Q2, xB, t, m_phi, beam_energy);
    }, {"Q2", "xB", "t"})
  .Define("Egamma_star",
    [](double Mx2) -> double {
      if (Mx2 <= 0.0) return -999.0;
      const double Mx = std::sqrt(Mx2);
      return (Mx2 - kMp * kMp) / (2.0 * kMp);
    }, {"Mx2_eKpKm"});
  return *df_;
}

// -----------------------------------------------------------------------------
// InitKinematics_MissingKp : K⁺ omitted (exclusive K⁻ channel)
// -----------------------------------------------------------------------------
ROOT::RDF::RNode InitKinematics_MissingKp(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);

  *df_ = df_->Define("ele_px_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recel_vz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             ////
             .Define("reckMinus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kMinus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kMinus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kMinus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("reckPlus_vz",  // dummy for K plus vz plotting
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})

             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
             .Define("RunNumber", "RUN_config_run")
             .Define("EventNumber", "RUN_config_event")
             .Define("nElectrons",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
                       int n = 0;
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) ++n;
                       return n;
                     },
                     {"REC_Particle_pid", "REC_Particle_pass"})
             // Choose ONE best electron per event
             .Define("bestEle_idx", ChooseBestElectron, {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz", "REC_Particle_vz", "REC_Particle_pass"})

             // If you want to **keep** multi-e events but only use the best e⁻, leave the next line.
             // If you want to **drop** multi-e events (strict DVEP), add: .Filter("nElectrons==1","Cut: exactly 1 e⁻")
             .Filter("bestEle_idx >= 0", "At least one usable e⁻")
             .Define("REC_pass_bestE",
                     [](int best, const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
                       ROOT::VecOps::RVec<int> keep(pass.size(), 0);
                       if (best >= 0) keep[best] = 1;        // only best e⁻
                       for (size_t i = 0; i < pid.size(); ++i)  // keep non-e tracks as they were
                         if (pid[i] != 11) keep[i] = pass[i];
                       return keep;
                     },
                     {"bestEle_idx", "REC_Particle_pid", "REC_Particle_pass"})

             .Define("nElectrons_best",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& keep) {
                       int n = 0;
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && keep[i]) ++n;
                       return n;
                     },
                     {"REC_Particle_pid", "REC_pass_bestE"})
             // Project the chosen e⁻ to the scalar columns used everywhere else
             .Define("ele_px", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_px"})
             .Define("ele_py", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_py"})
             .Define("ele_pz", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_pz"})
             .Define("recel_vz", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_vz"})
             .Define("ele_det_region", [](int i, const ROOT::VecOps::RVec<short>& st) { return i >= 0 ? DetRegionFromStatus(st[i]) : -1; }, {"bestEle_idx", "REC_Particle_status"});

  // missing K⁺ 4-vector (components)

  *df_ = df_->Define("kPlus_miss_px",
                     [beam_energy](float epx, float epy, float epz, float ppx, float ppy, float ppz, float kmx, float kmy, float kmz) -> float {
                       TLorentzVector pBeam(0, 0, beam_energy, beam_energy), pTarg(0, 0, 0, kMp);
                       TLorentzVector e(epx, epy, epz, std::sqrt(epx * epx + epy * epy + epz * epz + kMe * kMe));
                       TLorentzVector p(ppx, ppy, ppz, std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz + kMp * kMp));
                       TLorentzVector km(kmx, kmy, kmz, std::sqrt(kmx * kmx + kmy * kmy + kmz * kmz + kMK * kMK));
                       return float((pBeam + pTarg - e - p - km).Px());
                     },
                     {"ele_px", "ele_py", "ele_pz", "pro_px", "pro_py", "pro_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("kPlus_miss_py",
                     [beam_energy](float epx, float epy, float epz, float ppx, float ppy, float ppz, float kmx, float kmy, float kmz) -> float {
                       TLorentzVector pBeam(0, 0, beam_energy, beam_energy), pTarg(0, 0, 0, kMp);
                       TLorentzVector e(epx, epy, epz, std::sqrt(epx * epx + epy * epy + epz * epz + kMe * kMe));
                       TLorentzVector p(ppx, ppy, ppz, std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz + kMp * kMp));
                       TLorentzVector km(kmx, kmy, kmz, std::sqrt(kmx * kmx + kmy * kmy + kmz * kmz + kMK * kMK));
                       return float((pBeam + pTarg - e - p - km).Py());
                     },
                     {"ele_px", "ele_py", "ele_pz", "pro_px", "pro_py", "pro_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("kPlus_miss_pz",
                     [beam_energy](float epx, float epy, float epz, float ppx, float ppy, float ppz, float kmx, float kmy, float kmz) -> float {
                       TLorentzVector pBeam(0, 0, beam_energy, beam_energy), pTarg(0, 0, 0, kMp);
                       TLorentzVector e(epx, epy, epz, std::sqrt(epx * epx + epy * epy + epz * epz + kMe * kMe));
                       TLorentzVector p(ppx, ppy, ppz, std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz + kMp * kMp));
                       TLorentzVector km(kmx, kmy, kmz, std::sqrt(kmx * kmx + kmy * kmy + kmz * kmz + kMK * kMK));
                       return float((pBeam + pTarg - e - p - km).Pz());
                     },
                     {"ele_px", "ele_py", "ele_pz", "pro_px", "pro_py", "pro_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})

             .Filter([](float ex, float kMinusx, float px) { return ex != -999 && kMinusx != -999 && px != -999; }, {"ele_px", "kMinus_px", "pro_px"})
             // derived angles/magnitudes
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("reckMinus_p", MomentumFunc, {"kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("reckMinus_theta", ThetaFunc, {"kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("reckMinus_phi", PhiFunc, {"kMinus_px", "kMinus_py"})
             .Define("reckPlus_p", MomentumFunc, {"kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz"})
             .Define("reckPlus_theta", ThetaFunc, {"kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz"})
             .Define("reckPlus_phi", PhiFunc, {"kPlus_miss_px", "kPlus_miss_py"})
             .Define("kMinus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == -321 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
              .Define("kPlus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == -321 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("ele_det_region_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("invMass_KpKm",
                     [](float px1, float py1, float pz1, float px2, float py2, float pz2) -> float {
                       constexpr float mK = 0.493677;  // Kaon mass in GeV/c²
                       float E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mK * mK);
                       float E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mK * mK);
                       float px = px1 + px2;
                       float py = py1 + py2;
                       float pz = pz1 + pz2;
                       float E = E1 + E2;
                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})
             // p–K⁻ invariant mass
             .Define("invMass_pKminus",
                     [](float p_px, float p_py, float p_pz, float k_px, float k_py, float k_pz) -> float {
                       constexpr float mP = 0.938272f;  // Proton mass in GeV/c²
                       constexpr float mK = 0.493677f;  // Kaon mass in GeV/c²

                       float E_p = std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz + mP * mP);
                       float E_k = std::sqrt(k_px * k_px + k_py * k_py + k_pz * k_pz + mK * mK);

                       float px = p_px + k_px;
                       float py = p_py + k_py;
                       float pz = p_pz + k_pz;
                       float E = E_p + E_k;

                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"pro_px", "pro_py", "pro_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})

             // p–K⁺ invariant mass
             .Define("invMass_pKplus",
                     [](float p_px, float p_py, float p_pz, float k_px, float k_py, float k_pz) -> float {
                       constexpr float mP = 0.938272f;  // Proton mass in GeV/c²
                       constexpr float mK = 0.493677f;  // Kaon mass in GeV/c²

                       float E_p = std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz + mP * mP);
                       float E_k = std::sqrt(k_px * k_px + k_py * k_py + k_pz * k_pz + mK * mK);

                       float px = p_px + k_px;
                       float py = p_py + k_py;
                       float pz = p_pz + k_pz;
                       float E = E_p + E_k;

                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"pro_px", "pro_py", "pro_pz", "kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz"});
  // φ mass built from missing K+ and measured K-
  // DISANAMath-driven observables
  *df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "tmin", &DISANAMath::GetTmin, beam_energy);
  *df_ = df_->Define("mtprime",  // non-negative: this is what you’ll plot
                     [](double mt, double tmin) { return std::abs(mt + tmin); }, {"t", "tmin"})
             .Define("tprime",  // optional signed t' ≤ 0
                     [](double mtp) { return -mtp; }, {"mtprime"});
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);
  *df_ = define_DISCAT(*df_, "cos_thetaKK", &DISANAMath::GetCosTheta_KK, beam_energy);
  *df_ = define_DISCAT(*df_, "cos_phiKK", &DISANAMath::GetCosPhi_KK, beam_energy);

  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKp", &DISANAMath::GetMx2_eKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx_eKp", &DISANAMath::GetMx_eKp, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_p", &DISANAMath::GetCone_p, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Kp", &DISANAMath::GetCone_Kp, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Km", &DISANAMath::GetCone_Km, beam_energy);
  *df_ = define_DISCAT(*df_, "z_phi", &DISANAMath::GetZ_phi, beam_energy);
  *df_ = define_DISCAT(*df_, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  *df_ = df_->Define("Gamma_v", [beam_energy](double Q2, double W) {
                 return HandGammaV(beam_energy, Q2, W);
               }, {"Q2", "W"});
  // --- new: vm_cut and Egamma_star -----------------------------------------
  *df_ = df_->Define("vm_cut",
    [beam_energy](double Q2, double xB, double t) -> double {
      return DISANAMath::ComputeVmApprox(Q2, xB, t, m_phi, beam_energy);
    }, {"Q2", "xB", "t"})
  .Define("Egamma_star",
    [](double Mx2) -> double {
      if (Mx2 <= 0.0) return -999.0;
      const double Mx = std::sqrt(Mx2);
      return (Mx2 - kMp * kMp) / (2.0 * kMp);
    }, {"Mx2_eKpKm"});
  return *df_;
}

// =============================================================================
// ---------------------------------------------------------------------------
// Q2BinMeans  — holds the per-Q^2-bin mean kinematics computed from data
// ---------------------------------------------------------------------------
struct Q2BinMeans {
  std::vector<double> xB;  // mean x_B  in each Q^2 bin
  std::vector<double> Q2;  // mean Q^2  in each Q^2 bin
};

// ---------------------------------------------------------------------------
// ComputeMeanKinPerQ2Bin
//
// Books Mean("xB") and Mean("Q2") for every Q^2 bin in a single pass over
// the RDataFrame (all lazy actions are registered before any result is read).
// Bins with no entries receive fallback values (xB=0.1, Q2=bin centre).
// ---------------------------------------------------------------------------
inline Q2BinMeans ComputeMeanKinPerQ2Bin(ROOT::RDF::RNode df,
                                          const std::vector<double>& q2Bins)
{
  const int nBins = static_cast<int>(q2Bins.size()) - 1;
  Q2BinMeans out;
  if (nBins <= 0) return out;

  // Book all lazy actions first — single event loop
  std::vector<ROOT::RDF::RResultPtr<double>> meanXb, meanQ2;
  meanXb.reserve(nBins);
  meanQ2.reserve(nBins);
  for (int iq = 0; iq < nBins; ++iq) {
    const double lo = q2Bins[iq];
    const double hi = q2Bins[iq + 1];
    auto df_bin = df.Filter(
      [lo, hi](double Q2){ return Q2 >= lo && Q2 < hi; }, {"Q2"},
      Form("vmkin_q2bin_%d", iq));
    meanXb.push_back(df_bin.Mean("xB"));
    meanQ2.push_back(df_bin.Mean("Q2"));
  }

  // Retrieve (triggers the event loop once for all bins and both quantities)
  out.xB.resize(nBins);
  out.Q2.resize(nBins);
  for (int iq = 0; iq < nBins; ++iq) {
    const double vxb = *meanXb[iq];
    const double vq2 = *meanQ2[iq];
    const double q2_centre = 0.5 * (q2Bins[iq] + q2Bins[iq + 1]);
    out.xB[iq] = (std::isfinite(vxb) && vxb > 0.0) ? vxb : 0.1;
    out.Q2[iq] = (std::isfinite(vq2) && vq2 > 0.0) ? vq2 : q2_centre;
    std::cout << Form("[ComputeMeanKinPerQ2Bin] Q2 bin [%.3f, %.3f]:"
                      "  mean Q2 = %.4f  mean xB = %.4f\n",
                      q2Bins[iq], q2Bins[iq+1], out.Q2[iq], out.xB[iq]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// PlotVmCut
//
// For each Q^2 bin draws vm(t) = f(mean_Q2, mean_xB, |t|, mv, beam_energy)
// using the per-bin mean kinematics from ComputeMeanKinPerQ2Bin.
// One colour-coded TGraph per Q^2 bin on a log-x canvas, saved to outFile.
// The legend shows both the bin edges and the actual mean Q^2 / xB used.
//
// Arguments:
//   q2Bins      - Q^2 bin-edge vector, e.g. {1.0, 2.0, 4.0, 6.0}
//   means       - Q2BinMeans struct from ComputeMeanKinPerQ2Bin
//   beam_energy - beam energy [GeV]
//   mv          - vector meson mass [GeV]  (default: phi mass 1.019461)
//   outFile     - output PDF / PNG file path
//   nPoints     - |t| scan points per curve (default 300, log-spaced)
// ---------------------------------------------------------------------------
inline void PlotVmCut(const std::vector<double>& q2Bins,
                      const Q2BinMeans& means,
                      double beam_energy,
                      double mv = m_phi,
                      const std::string& outFile = "vm_cut.pdf",
                      int nPoints = 300)
{
  const int nBins = static_cast<int>(q2Bins.size()) - 1;
  if (nBins <= 0) {
    std::cerr << "[PlotVmCut] ERROR: q2Bins must have at least 2 entries.\n";
    return;
  }

  const bool xb_ok = (static_cast<int>(means.xB.size()) == nBins);
  const bool q2_ok = (static_cast<int>(means.Q2.size()) == nBins);
  if (!xb_ok || !q2_ok)
    std::cerr << "[PlotVmCut] WARNING: means size mismatch, using fallbacks.\n";

  const std::vector<int> colours = {
    kBlue+1, kRed+1, kGreen+2, kOrange+1, kMagenta+1,
    kCyan+2, kViolet+1, kTeal+2, kPink+1, kAzure+2
  };

  auto* mg  = new TMultiGraph();
  auto* leg = new TLegend(0.2, 0.35, 0.6, 0.65);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.028);
  leg->SetHeader("Q^{2} bin [lo, hi]  (#bar{Q}^{2}, #bar{x}_{B})", "C");

  double vm_max_global = 0.0;

  for (int iq = 0; iq < nBins; ++iq) {
    const double q2_lo = q2Bins[iq];
    const double q2_hi = q2Bins[iq + 1];

    // Use data mean Q^2 and mean x_B for this bin
    const double q2_c  = (q2_ok  && means.Q2[iq] > 0.0) ? means.Q2[iq]
                                  : 0.5 * (q2_lo + q2_hi);
    const double xB_c  = (xb_ok  && means.xB[iq] > 0.0) ? means.xB[iq] : 0.1;

    // Kinematic endpoint: vm = 0 when |t| = (W^2 - mv^2) / 2
    const double S        = 2.0 * m_p * beam_energy;
    const double y_c      = (xB_c * S > 0.0) ? q2_c / (xB_c * S) : 0.0;
    const double Sx_c     = y_c * S;
    const double W2_c     = m_p * m_p + Sx_c - q2_c;
    const double tmax_kin = (W2_c > mv * mv) ? (W2_c - mv * mv) / 2.0 : 1.0;
    const double t_hi     = std::min(tmax_kin * 0.99, 5.0);

    if (t_hi <= 0.0) continue;

    auto* gr = new TGraph(nPoints);
    double vm_bin_max = 0.0;

    for (int ip = 0; ip < nPoints; ++ip) {
      // Log-spaced |t| for better resolution at small-t
      const double log_lo = std::log10(1e-3);
      const double log_hi = std::log10(t_hi);
      const double t_abs  = std::pow(10.0,
                              log_lo + ip * (log_hi - log_lo) / (nPoints - 1));
      const double vm = DISANAMath::ComputeVmApprox(q2_c, xB_c, t_abs,
                                                     mv, beam_energy);
      gr->SetPoint(ip, t_abs, vm);
      if (vm > vm_bin_max) vm_bin_max = vm;
    }
    vm_max_global = std::max(vm_max_global, vm_bin_max);

    const int col = colours[iq % static_cast<int>(colours.size())];
    gr->SetLineColor(col);
    gr->SetLineWidth(2);
    mg->Add(gr, "L");

    // Legend: bin edges + actual mean Q^2 and x_B used for the curve
    leg->AddEntry(gr,
      Form("[%.2f, %.2f]  (#bar{Q}^{2}=%.3f, #bar{x}_{B}=%.3f)",
           q2_lo, q2_hi, q2_c, xB_c),
      "L");
  }

  TCanvas* c = new TCanvas("cVmCut", "vm cut vs -t", 1000, 680);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.04);
  c->SetLogx();
  c->SetGrid();

  mg->SetTitle(";-t [GeV^{2}];v_{m} [GeV^{2}]");
  mg->Draw("A");
  mg->GetXaxis()->SetTitle("-t [GeV^{2}]");
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetTitle("v_{m} [GeV^{2}]");
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetLabelSize(0.04);
  if (vm_max_global > 0.0)
    mg->GetYaxis()->SetRangeUser(0.0, vm_max_global * 1.18);

  TLatex lbl;
  lbl.SetNDC();
  lbl.SetTextFont(42);
  lbl.SetTextSize(0.031);
  lbl.DrawLatex(0.2, 0.78,
    Form("E_{beam} = %.2f GeV,  m_{V} = %.4f GeV", beam_energy, mv));

  leg->Draw();
  c->Update();
  c->SaveAs(outFile.c_str());
  std::cout << "[PlotVmCut] Saved -> " << outFile << "\n";
}


// =============================================================================
// PlotRadiativeKinematics  (v2 – beautified, unified x-axis, difference panel)
//
// Each cell contains TWO vertically-stacked panels:
//   TOP    (65 %)  – E*_gamma and E_miss overlaid, normalised to unit area
//   BOTTOM (35 %)  – Difference  ΔN = E*_gamma − E_miss
//
// Layout:
//   Cells 1..nBins  – per-Q² bin
//   Last cell        – all-Q² integrated
//   Top legend strip – variable identification
//
// Variable styling:
//   E*_gamma  deep blue  kAzure-4,  solid,  width 2, light fill
//   E_miss    warm red   kOrange+7, dashed, width 2, light fill
//   Difference  steel grey kGray+1, filled
//
// Axis range (unified for both variables):
//   Default xLo = -0.4,  xHi = 1.0   [GeV]
//
// Arguments:
//   df         – RDataFrame with "Q2", "Egamma_star" (or "Mx2_eKpKm"), "Emiss"
//   q2Bins     – Q² bin-edge vector
//   means      – Q2BinMeans from ComputeMeanKinPerQ2Bin
//   outFile    – output PDF / PNG path
//   nBinsHist  – histogram bins per variable (default 80)
//   xLo / xHi – shared x-axis range [GeV]   (default -0.4 / 1.0)
// =============================================================================
inline void PlotRadiativeKinematics(ROOT::RDF::RNode df,
                                     const std::vector<double>& q2Bins,
                                     const Q2BinMeans& means,
                                     const std::string& outFile = "rad_kinematics.pdf",
                                     int    nBinsHist = 80,
                                     double xLo      = -0.4,
                                     double xHi      =  1.0)
{
  const int nBins = static_cast<int>(q2Bins.size()) - 1;
  if (nBins <= 0) {
    std::cerr << "[PlotRadiativeKinematics] ERROR: q2Bins must have >= 2 entries.\n";
    return;
  }
  const bool q2_ok = (static_cast<int>(means.Q2.size()) == nBins);
  if (!q2_ok)
    std::cerr << "[PlotRadiativeKinematics] WARNING: means.Q2 size mismatch.\n";
 
  // ── Ensure Egamma_star column exists ───────────────────────────────────
  auto colNames = df.GetColumnNames();
  auto hasCol   = [&colNames](const std::string& n) {
    return std::find(colNames.begin(), colNames.end(), n) != colNames.end();
  };
  if (!hasCol("Egamma_star")) {
    std::cout << "[PlotRadiativeKinematics] INFO: computing Egamma_star"
                 " on-the-fly from Mx2_eKpKm.\n";
    df = df.Define("Egamma_star",
      [](double Mx2) -> double {
        if (Mx2 <= 0.0) return -999.0;
        const double Mx = std::sqrt(Mx2);
        return (Mx2 - kMp * kMp) / (2.0 * kMp);
      }, {"Mx2_eKpKm"});
  }
 
  // ── Global style tweaks (scoped; restored after) ────────────────────────
  const int    savedOptStat  = gStyle->GetOptStat();
  const int    savedOptTitle = gStyle->GetOptTitle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);               // titles drawn manually as TPaveText
  gStyle->SetFrameLineWidth(1);
  gStyle->SetTickLength(0.03, "X");
  gStyle->SetTickLength(0.03, "Y");
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(504, "Y");
  gStyle->SetEndErrorSize(4);
 
  // ── Palette ─────────────────────────────────────────────────────────────
  const int kColEgam   = kAzure  - 4;   // deep blue
  const int kColEmiss  = kOrange + 7;   // warm red
  const int kColDiff   = kGray   + 1;   // steel grey for difference
  const int kLsDashed  = 7;             // kDashed style
 
  // ── Book Histo1D lazily (single event loop) ─────────────────────────────
  using H1Ptr = ROOT::RDF::RResultPtr<TH1D>;
  std::vector<H1Ptr> hbEgam(nBins), hbEmiss(nBins);
  H1Ptr hIntEgam, hIntEmiss;
 
  for (int iq = 0; iq < nBins; ++iq) {
    const double lo = q2Bins[iq], hi = q2Bins[iq + 1];
    auto df_bin = df.Filter(
      [lo, hi](double Q2){ return Q2 >= lo && Q2 < hi; }, {"Q2"},
      Form("radkin_bin_%d", iq));
    hbEgam[iq]  = df_bin.Histo1D(
      {Form("hbEgam%d",  iq), "", nBinsHist, xLo, xHi}, "Egamma_star");
    hbEmiss[iq] = df_bin.Histo1D(
      {Form("hbEmiss%d", iq), "", nBinsHist, xLo, xHi}, "Emiss");
  }
  hIntEgam  = df.Histo1D({"hIntEgam",  "", nBinsHist, xLo, xHi}, "Egamma_star");
  hIntEmiss = df.Histo1D({"hIntEmiss", "", nBinsHist, xLo, xHi}, "Emiss");
 
  // ── Helper: clone, normalize, style ─────────────────────────────────────
  auto prepHist = [](H1Ptr& src, const char* newName,
                     int col, int lstyle, int lwidth,
                     float alpha = 0.10f) -> TH1D* {
    auto* h = static_cast<TH1D*>(src->Clone(newName));
    h->SetDirectory(nullptr);
    const double intg = h->Integral();
    if (intg > 0) h->Scale(1.0 / intg);
    h->SetLineColor(col);
    h->SetLineStyle(lstyle);
    h->SetLineWidth(lwidth);
    h->SetFillColorAlpha(col, alpha);
    h->SetFillStyle(1001);
    h->SetStats(0);
    return h;
  };
 
  // ── Retrieve and prepare all histograms ─────────────────────────────────
  std::vector<TH1D*> vEgam(nBins), vEmiss(nBins), vDiff(nBins);
  for (int iq = 0; iq < nBins; ++iq) {
    vEgam[iq]  = prepHist(hbEgam[iq],  Form("cEgam%d",  iq), kColEgam,  1,       2, 0.10f);
    vEmiss[iq] = prepHist(hbEmiss[iq], Form("cEmiss%d", iq), kColEmiss, kLsDashed, 2, 0.10f);
    // Difference histogram: reuse Egamma clone, then subtract Emiss
    vDiff[iq] = static_cast<TH1D*>(vEgam[iq]->Clone(Form("cDiff%d", iq)));
    vDiff[iq]->SetDirectory(nullptr);
    vDiff[iq]->Add(vEmiss[iq], -1.0);
    vDiff[iq]->SetLineColor(kColDiff);
    vDiff[iq]->SetLineStyle(1);
    vDiff[iq]->SetLineWidth(1);
    vDiff[iq]->SetFillColorAlpha(kColDiff, 0.40f);
    vDiff[iq]->SetFillStyle(1001);
  }
  TH1D* hIE = prepHist(hIntEgam,  "cIntEgam",  kColEgam,  1,       2, 0.10f);
  TH1D* hIM = prepHist(hIntEmiss, "cIntEmiss", kColEmiss, kLsDashed, 2, 0.10f);
  TH1D* hID = static_cast<TH1D*>(hIE->Clone("cIntDiff"));
  hID->SetDirectory(nullptr);
  hID->Add(hIM, -1.0);
  hID->SetLineColor(kColDiff); hID->SetLineStyle(1); hID->SetLineWidth(1);
  hID->SetFillColorAlpha(kColDiff, 0.40f); hID->SetFillStyle(1001);
 
  // ── Canvas geometry ──────────────────────────────────────────────────────
  const int nPads = nBins + 1;
  const int nCols = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(nPads))));
  const int nRows = static_cast<int>(std::ceil(static_cast<double>(nPads) / nCols));
 
  const double legFrac  = 0.07;          // top legend strip height fraction
  const double padBot   = 0.0;
  const double padTop   = 1.0 - legFrac;
 
  const int cW = std::max(600, 380 * nCols);
  const int cH = static_cast<int>(std::max(500, 380 * nRows) / (1.0 - legFrac));
 
  TCanvas* c = new TCanvas("cRadKin", "Radiative kinematics", cW, cH);
  c->SetFillColorAlpha(0, 0);  // transparent background
 
  // ── Legend strip ─────────────────────────────────────────────────────────
  c->cd();
  TPad* pLeg = new TPad("pLeg", "", 0.0, padTop, 1.0, 1.0);
  pLeg->SetFillColorAlpha(0, 0); pLeg->SetBorderSize(0);
  pLeg->Draw(); pLeg->cd();
 
  auto* leg = new TLegend(0.04, 0.05, 0.96, 0.95);
  leg->SetNColumns(2);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(43); leg->SetTextSize(20);
  // Dummy objects carrying the visual style
  auto* dEgam  = new TH1F("_dleg_egam",  "", 1, 0, 1);
  dEgam->SetLineColor(kColEgam);  dEgam->SetLineWidth(3); dEgam->SetLineStyle(1);
  dEgam->SetFillColorAlpha(kColEgam,  0.20f); dEgam->SetFillStyle(1001);
  auto* dEmiss = new TH1F("_dleg_emiss", "", 1, 0, 1);
  dEmiss->SetLineColor(kColEmiss); dEmiss->SetLineWidth(3); dEmiss->SetLineStyle(kLsDashed);
  dEmiss->SetFillColorAlpha(kColEmiss, 0.20f); dEmiss->SetFillStyle(1001);
  auto* dDiff  = new TH1F("_dleg_diff",  "", 1, 0, 1);
  dDiff->SetLineColor(kColDiff);  dDiff->SetLineWidth(2); dDiff->SetLineStyle(1);
  dDiff->SetFillColorAlpha(kColDiff,  0.45f); dDiff->SetFillStyle(1001);
  leg->AddEntry(dEgam,  "E*_{#gamma} = (M_{X}^{2} #minus M_{p}^{2}) / 2M_{p}", "LF");
  leg->AddEntry(dEmiss, "E_{miss} = E_{beam} + M_{p} #minus #sum E_{final}",    "LF");
  leg->AddEntry(dDiff,  "#DeltaN = E*_{#gamma} #minus E_{miss}",                "LF");
  leg->Draw();
 
  // ── Helper: draw one cell (outer pad split into top/bottom sub-pads) ────
  // splitFrac = fraction of the outer pad height used by the top (main) panel
  const double splitFrac = 0.63;
 
  auto drawCell = [&](int iPad,
                      TH1D* hEgam, TH1D* hMiss, TH1D* hDiff,
                      const std::string& title)
  {
    const int iRow = iPad / nCols;
    const int iCol = iPad % nCols;
    const double x0 = static_cast<double>(iCol)     / nCols;
    const double x1 = static_cast<double>(iCol + 1) / nCols;
    const double y1 = padTop - static_cast<double>(iRow)     / nRows * padTop;
    const double y0 = padTop - static_cast<double>(iRow + 1) / nRows * padTop;
 
    // Outer pad (no border, transparent)
    c->cd();
    TPad* pOuter = new TPad(Form("pOuter%d", iPad), "", x0, y0, x1, y1);
    pOuter->SetFillColorAlpha(0, 0); pOuter->SetBorderSize(0);
    pOuter->SetLeftMargin(0); pOuter->SetRightMargin(0);
    pOuter->SetTopMargin(0);  pOuter->SetBottomMargin(0);
    pOuter->Draw(); pOuter->cd();
 
    // ── Top sub-pad (main distributions) ───────────────────────────────
    TPad* pTop = new TPad(Form("pTop%d", iPad), "",
                          0.0, 1.0 - splitFrac, 1.0, 1.0);
    pTop->SetFillColorAlpha(0, 0); pTop->SetBorderSize(0);
    pTop->SetLeftMargin(0.20); pTop->SetRightMargin(0.04);
    pTop->SetTopMargin(0.18);  pTop->SetBottomMargin(0.02);
    pTop->SetTickx(1); pTop->SetTicky(1);
    pTop->SetGridx(); pTop->SetGridy();
    pTop->Draw(); pTop->cd();
 
    const double ymax = std::max(hEgam->GetMaximum(), hMiss->GetMaximum()) * 1.32;
 
    // --- Egamma_star ---
    hEgam->GetXaxis()->SetLabelSize(0);        // no x labels on top panel
    hEgam->GetXaxis()->SetTitleSize(0);
    hEgam->GetXaxis()->SetTickLength(0.04);
    hEgam->GetYaxis()->SetTitle("Norm. counts");
    hEgam->GetYaxis()->SetTitleFont(43);
    hEgam->GetYaxis()->SetTitleSize(14);
    hEgam->GetYaxis()->SetTitleOffset(1.6);
    hEgam->GetYaxis()->SetLabelFont(43);
    hEgam->GetYaxis()->SetLabelSize(12);
    hEgam->GetYaxis()->SetNdivisions(504);
    hEgam->GetYaxis()->SetRangeUser(0.0, ymax);
    hEgam->Draw("HIST");
    hMiss->GetYaxis()->SetRangeUser(0.0, ymax);
    hMiss->Draw("HIST SAME");
 
    // Vertical reference line at x = 0
    TLine* l0 = new TLine(0.0, 0.0, 0.0, ymax);
    l0->SetLineColor(kGray + 2); l0->SetLineWidth(1); l0->SetLineStyle(3);
    l0->Draw();
 
    // Pad title as TPaveText (top-right corner, inside pad)
    double ptx0 = 0.22, ptx1 = 0.97, pty0 = 0.76, pty1 = 0.97;
    TPaveText* pt = new TPaveText(ptx0, pty0, ptx1, pty1, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.65); pt->SetBorderSize(0);
    pt->SetTextFont(43); pt->SetTextSize(13);
    pt->SetTextAlign(22);
    pt->AddText(title.c_str());
    pt->Draw();
 
    // ── Bottom sub-pad (difference) ─────────────────────────────────────
    pOuter->cd();
    TPad* pBot = new TPad(Form("pBot%d", iPad), "",
                          0.0, 0.0, 1.0, 1.0 - splitFrac);
    pBot->SetFillColorAlpha(0, 0); pBot->SetBorderSize(0);
    pBot->SetLeftMargin(0.20); pBot->SetRightMargin(0.04);
    pBot->SetTopMargin(0.02);  pBot->SetBottomMargin(0.34);
    pBot->SetTickx(1); pBot->SetTicky(1);
    pBot->SetGridx(); pBot->SetGridy();
    pBot->Draw(); pBot->cd();
 
    // Symmetrise the y-axis of the difference around zero
    const double diffMax = std::max(std::abs(hDiff->GetMinimum()),
                                     std::abs(hDiff->GetMaximum())) * 1.45;
    const double diffRange = (diffMax > 0.0) ? diffMax : 0.05;
 
    hDiff->GetXaxis()->SetTitle("Energy  [GeV]");
    hDiff->GetXaxis()->SetTitleFont(43); hDiff->GetXaxis()->SetTitleSize(14);
    hDiff->GetXaxis()->SetTitleOffset(2.4);
    hDiff->GetXaxis()->SetLabelFont(43); hDiff->GetXaxis()->SetLabelSize(12);
    hDiff->GetXaxis()->SetTickLength(0.07);
    hDiff->GetYaxis()->SetTitle("#DeltaN");
    hDiff->GetYaxis()->SetTitleFont(43); hDiff->GetYaxis()->SetTitleSize(13);
    hDiff->GetYaxis()->SetTitleOffset(1.6);
    hDiff->GetYaxis()->SetLabelFont(43); hDiff->GetYaxis()->SetLabelSize(11);
    hDiff->GetYaxis()->SetNdivisions(503);
    hDiff->GetYaxis()->SetRangeUser(-diffRange, +diffRange);
    hDiff->Draw("HIST");
 
    // Zero reference line
    TLine* lz = new TLine(xLo, 0.0, xHi, 0.0);
    lz->SetLineColor(kBlack); lz->SetLineWidth(1); lz->SetLineStyle(2);
    lz->Draw();
  };
 
  // ── Draw per-Q² cells ───────────────────────────────────────────────────
  for (int iq = 0; iq < nBins; ++iq) {
    const double q2_mean = (q2_ok && means.Q2[iq] > 0.0)
                           ? means.Q2[iq]
                           : 0.5 * (q2Bins[iq] + q2Bins[iq + 1]);
    const std::string title = Form(
      "Q^{2} #in [%.2f, %.2f],  #bar{Q}^{2} = %.3f",
      q2Bins[iq], q2Bins[iq + 1], q2_mean);
    drawCell(iq, vEgam[iq], vEmiss[iq], vDiff[iq], title);
  }
 
  // ── Integrated cell (last) ───────────────────────────────────────────────
  drawCell(nBins, hIE, hIM, hID,
    Form("All Q^{2} #in [%.2f, %.2f]", q2Bins.front(), q2Bins.back()));
 
  // ── Save and restore style ───────────────────────────────────────────────
  c->cd();
  c->Update();
  c->SaveAs(outFile.c_str());
  std::cout << "[PlotRadiativeKinematics] Saved -> " << outFile << "\n";
 
  gStyle->SetOptStat(savedOptStat);
  gStyle->SetOptTitle(savedOptTitle);
}




// =============================================================================
// PlotInvMassVsMissingMass
//
// Single 2D COLZ heat-map of  M(K+K-)  vs  MM(ep),  integrated over all Q².
// Reference lines mark phi(1020) and M_p.
//
// Arguments:
//   df         – RDataFrame with "invMass_KpKm" and "Mx2_ep" (or "MM_ep")
//   outFile    – output PDF path
//   nBinsX     – bins along x = M(K+K-)  (default 100)
//   nBinsY     – bins along y = MM(ep)   (default 100)
//   xLo/xHi   – M(K+K-)  range [GeV]    (default 0.98 / 1.10)
//   yLo/yHi   – MM(ep)   range [GeV]    (default 0.70 / 1.20)
// =============================================================================
inline void PlotInvMassVsMissingMass(
    ROOT::RDF::RNode    df,
    const std::string&  outFile = "invmass_vs_missmass.pdf",
    int    nBinsX = 100,
    int    nBinsY = 100,
    double xLo    = 0.98,
    double xHi    = 1.10,
    double yLo    = 0.70,
    double yHi    = 1.20)
{
  // ── Ensure MM_ep column exists ────────────────────────────────────────────
  {
    auto cols = df.GetColumnNames();
    bool hasMM  = std::find(cols.begin(), cols.end(), "MM_ep")  != cols.end();
    bool hasMx2 = std::find(cols.begin(), cols.end(), "Mx2_ep") != cols.end();
    if (!hasMM) {
      if (!hasMx2) {
        std::cerr << "[PlotInvMassVsMissingMass] ERROR: neither MM_ep nor Mx2_ep found.\n";
        return;
      }
      df = df.Define("MM_ep", [](double Mx2) -> double {
        return (Mx2 > 0.0) ? std::sqrt(Mx2) : -999.0;
      }, {"Mx2_ep"});
    }
  }

  static constexpr double kMphi = 1.01946;

  // ── Global style (scoped) ────────────────────────────────────────────────
  const int savedOptStat  = gStyle->GetOptStat();
  const int savedOptTitle = gStyle->GetOptTitle();
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetPalette(kViridis);

  // ── Fill 2D histogram (single event loop) ────────────────────────────────
  // Unique names prevent ROOT global-directory collisions across calls.
  static std::atomic<unsigned long> uid2d{0};
  const auto hname = Form("h2_imm_%lu", uid2d.fetch_add(1));
  auto h2R = df.Histo2D(
      {hname, "",
       nBinsX, xLo, xHi,
       nBinsY, yLo, yHi},
      "invMass_KpKm", "MM_ep");
  TH2D* h2 = (TH2D*)h2R.GetPtr()->Clone(Form("%s_draw", hname));
  h2->SetDirectory(nullptr);

  // ── Canvas ────────────────────────────────────────────────────────────────
  TCanvas* c = new TCanvas("cInvMassVsMM", "M(K^{+}K^{-}) vs MM(ep)", 800, 700);
  c->SetLeftMargin(0.14); c->SetRightMargin(0.16);
  c->SetTopMargin(0.08);  c->SetBottomMargin(0.12);
  c->SetLogz();
  c->SetTickx(1); c->SetTicky(1);

  // ── Axes ─────────────────────────────────────────────────────────────────
  h2->GetXaxis()->SetTitle("M(K^{+}K^{-})  [GeV]");
  h2->GetXaxis()->SetTitleFont(43); h2->GetXaxis()->SetTitleSize(18);
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetXaxis()->SetLabelFont(43); h2->GetXaxis()->SetLabelSize(15);
  h2->GetXaxis()->SetNdivisions(505);

  h2->GetYaxis()->SetTitle("MM(ep) = #sqrt{M_{X}^{2}(ep)}  [GeV]");
  h2->GetYaxis()->SetTitleFont(43); h2->GetYaxis()->SetTitleSize(16);
  h2->GetYaxis()->SetTitleOffset(1.6);
  h2->GetYaxis()->SetLabelFont(43); h2->GetYaxis()->SetLabelSize(15);
  h2->GetYaxis()->SetNdivisions(504);

  h2->GetZaxis()->SetTitle("Counts / bin");
  h2->GetZaxis()->SetTitleFont(43); h2->GetZaxis()->SetTitleSize(14);
  h2->GetZaxis()->SetLabelFont(43); h2->GetZaxis()->SetLabelSize(12);

  h2->Draw("COLZ");

  // ── Reference lines ───────────────────────────────────────────────────────
  // x: phi(1020) mass  — signal peak in M(K+K-)
  TLine lPhiX(kMphi, yLo, kMphi, yHi);
  lPhiX.SetLineColor(kRed+1); lPhiX.SetLineWidth(2); lPhiX.SetLineStyle(7);
  lPhiX.Draw();

  // y: phi(1020) mass  — MM(ep) = missing mass of (e,p) system = M(KK) for
  //    exclusive phi production, so it also peaks at M_phi, NOT at M_p.
  TLine lPhiY(xLo, kMphi, xHi, kMphi);
  lPhiY.SetLineColor(kAzure-4); lPhiY.SetLineWidth(2); lPhiY.SetLineStyle(7);
  lPhiY.Draw();

  // ── Labels for reference lines ────────────────────────────────────────────
  TLatex ltx;
  ltx.SetTextFont(43); ltx.SetTextSize(14); ltx.SetNDC(false);
  ltx.SetTextColor(kRed+1); ltx.SetTextAngle(90);
  ltx.DrawLatex(kMphi + 0.002, yLo + 0.05*(yHi-yLo), "#phi(1020)");
  ltx.SetTextColor(kAzure-4); ltx.SetTextAngle(0);
  ltx.DrawLatex(xLo + 0.01*(xHi-xLo), kMphi + 0.012*(yHi-yLo), "#phi(1020)");

  // ── Title ─────────────────────────────────────────────────────────────────
  TLatex title;
  title.SetNDC(); title.SetTextFont(43); title.SetTextSize(16);
  title.SetTextAlign(21);
  title.DrawLatex(0.50, 0.94,
    "M(K^{+}K^{-}) vs MM(ep)  [all Q^{2} integrated]");

  c->Update();
  c->SaveAs(outFile.c_str());
  std::cout << "[PlotInvMassVsMissingMass] Saved -> " << outFile << "\n";

  gStyle->SetOptStat(savedOptStat);
  gStyle->SetOptTitle(savedOptTitle);
}


// =============================================================================
// PlotInvMassVsMissingMassProjections
//
// 1D projections of the 2D M(K+K-) vs MM(ep) distribution, side-by-side:
//   Left pad  – X projection: M(K+K-)  (full y range)
//   Right pad – Y projection: MM(ep)   projected from the phi SIGNAL WINDOW
//               in x only (invMass_KpKm in [kMphi-sigWin, kMphi+sigWin]).
//               This suppresses combinatorial background so the phi peak in
//               MM(ep) is clearly visible.  The background under the window
//               is shown as a hatched histogram for comparison.
//
// Arguments:
//   df       – RDataFrame with "invMass_KpKm" and "Mx2_ep" (or "MM_ep")
//   outFile  – output PDF path
//   nBinsX   – bins along x = M(K+K-)              (default 100)
//   nBinsY   – bins along y = MM(ep)               (default 80)
//   xLo/xHi – M(K+K-)  range [GeV]                (default 0.98 / 1.10)
//   yLo/yHi – MM(ep)   range [GeV]                (default 0.90 / 1.15)
//   sigWin   – half-width of phi signal window [GeV] (default 0.010 = ~2.5σ)
// =============================================================================
inline void PlotInvMassVsMissingMassProjections(
    ROOT::RDF::RNode    df,
    const std::string&  outFile = "invmass_missmass_projections.pdf",
    int    nBinsX  = 100,
    int    nBinsY  = 80,
    double xLo     = 0.98,
    double xHi     = 1.10,
    double yLo     = 0.90,
    double yHi     = 1.15,
    double sigWin  = 0.010)
{
  // ── Ensure MM_ep column exists ────────────────────────────────────────────
  {
    auto cols = df.GetColumnNames();
    bool hasMM  = std::find(cols.begin(), cols.end(), "MM_ep")  != cols.end();
    bool hasMx2 = std::find(cols.begin(), cols.end(), "Mx2_ep") != cols.end();
    if (!hasMM) {
      if (!hasMx2) {
        std::cerr << "[PlotInvMassVsMissingMassProjections] ERROR: neither MM_ep nor Mx2_ep found.\n";
        return;
      }
      df = df.Define("MM_ep", [](double Mx2) -> double {
        return (Mx2 > 0.0) ? std::sqrt(Mx2) : -999.0;
      }, {"Mx2_ep"});
    }
  }

  static constexpr double kMphi = 1.01946;

  // Signal and sideband windows in M(K+K-)
  // Signal:   [kMphi - sigWin, kMphi + sigWin]
  // Sideband: [xLo, kMphi - 3*sigWin] ∪ [kMphi + 3*sigWin, xHi]  (background estimate)
  const double sigLo  = kMphi - sigWin;
  const double sigHi  = kMphi + sigWin;
  const double bkgLo1 = xLo;
  const double bkgHi1 = kMphi - 3.0 * sigWin;
  const double bkgLo2 = kMphi + 3.0 * sigWin;
  const double bkgHi2 = xHi;

  // ── Global style (scoped) ────────────────────────────────────────────────
  const int savedOptStat  = gStyle->GetOptStat();
  const int savedOptTitle = gStyle->GetOptTitle();
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(1);

  // ── Book all histograms lazily (single event loop) ────────────────────────
  static std::atomic<unsigned long> uidP{0};
  const unsigned long uid = uidP.fetch_add(1);

  // X projection: full y range
  auto hXR = df.Histo1D(
      {Form("hProjKK_%lu", uid), "", nBinsX, xLo, xHi},
      "invMass_KpKm");

  // Y projection — signal window in M(K+K-)
  auto df_sig = df.Filter(
      [sigLo, sigHi](float m){ return m >= sigLo && m <= sigHi; },
      {"invMass_KpKm"});
  auto hYsigR = df_sig.Histo1D(
      {Form("hProjMM_sig_%lu", uid), "", nBinsY, yLo, yHi},
      "MM_ep");

  // Y projection — sidebands (background control)
  auto df_bkg = df.Filter(
      [bkgLo1, bkgHi1, bkgLo2, bkgHi2](float m){
        return (m >= bkgLo1 && m <= bkgHi1) ||
               (m >= bkgLo2 && m <= bkgHi2);
      }, {"invMass_KpKm"});
  auto hYbkgR = df_bkg.Histo1D(
      {Form("hProjMM_bkg_%lu", uid), "", nBinsY, yLo, yHi},
      "MM_ep");

  // Force evaluation
  hXR.GetValue();

  // Clone and detach from ROOT directory
  TH1D* hX    = (TH1D*)hXR.GetPtr()   ->Clone(Form("cKK_%lu",     uid));
  TH1D* hYsig = (TH1D*)hYsigR.GetPtr()->Clone(Form("cMMsig_%lu",  uid));
  TH1D* hYbkg = (TH1D*)hYbkgR.GetPtr()->Clone(Form("cMMbkg_%lu",  uid));
  hX->SetDirectory(nullptr);
  hYsig->SetDirectory(nullptr);
  hYbkg->SetDirectory(nullptr);

  // Normalise to unit area
  if (hX   ->Integral() > 0) hX   ->Scale(1.0 / hX   ->Integral());
  if (hYsig->Integral() > 0) hYsig->Scale(1.0 / hYsig->Integral());
  if (hYbkg->Integral() > 0) hYbkg->Scale(1.0 / hYbkg->Integral());

  // ── Style ─────────────────────────────────────────────────────────────────
  // M(K+K-) projection
  hX->SetLineColor(kAzure-4); hX->SetLineWidth(3);
  hX->SetFillColorAlpha(kAzure-4, 0.15f); hX->SetFillStyle(1001);
  hX->GetXaxis()->SetTitle("M(K^{+}K^{-})  [GeV]");
  hX->GetXaxis()->SetTitleFont(43); hX->GetXaxis()->SetTitleSize(18);
  hX->GetXaxis()->SetTitleOffset(1.2);
  hX->GetXaxis()->SetLabelFont(43); hX->GetXaxis()->SetLabelSize(15);
  hX->GetYaxis()->SetTitle("Norm. counts");
  hX->GetYaxis()->SetTitleFont(43); hX->GetYaxis()->SetTitleSize(18);
  hX->GetYaxis()->SetTitleOffset(1.5);
  hX->GetYaxis()->SetLabelFont(43); hX->GetYaxis()->SetLabelSize(15);
  hX->GetYaxis()->SetRangeUser(0.0, hX->GetMaximum() * 1.35);
  hX->GetYaxis()->SetNdivisions(504);

  // MM(ep) signal
  hYsig->SetLineColor(kRed-4); hYsig->SetLineWidth(3);
  hYsig->SetFillColorAlpha(kRed-4, 0.20f); hYsig->SetFillStyle(1001);
  hYsig->GetXaxis()->SetTitle("MM(ep) = #sqrt{M_{X}^{2}(ep)}  [GeV]");
  hYsig->GetXaxis()->SetTitleFont(43); hYsig->GetXaxis()->SetTitleSize(16);
  hYsig->GetXaxis()->SetTitleOffset(1.2);
  hYsig->GetXaxis()->SetLabelFont(43); hYsig->GetXaxis()->SetLabelSize(15);
  hYsig->GetYaxis()->SetTitle("Norm. counts");
  hYsig->GetYaxis()->SetTitleFont(43); hYsig->GetYaxis()->SetTitleSize(18);
  hYsig->GetYaxis()->SetTitleOffset(1.5);
  hYsig->GetYaxis()->SetLabelFont(43); hYsig->GetYaxis()->SetLabelSize(15);
  const double ymaxMM = std::max(hYsig->GetMaximum(), hYbkg->GetMaximum()) * 1.40;
  hYsig->GetYaxis()->SetRangeUser(0.0, ymaxMM);
  hYsig->GetYaxis()->SetNdivisions(504);

  // MM(ep) background sideband
  hYbkg->SetLineColor(kGray+1); hYbkg->SetLineWidth(2); hYbkg->SetLineStyle(7);
  hYbkg->SetFillColorAlpha(kGray+1, 0.20f); hYbkg->SetFillStyle(3354);

  // ── Canvas: two side-by-side pads ─────────────────────────────────────────
  TCanvas* c = new TCanvas("cProjComp", "Projections", 1400, 600);
  c->SetFillColorAlpha(0, 0);

  TPad* pL = new TPad("pLeft",  "", 0.0,  0.0, 0.50, 1.0);
  pL->SetLeftMargin(0.16); pL->SetRightMargin(0.05);
  pL->SetTopMargin(0.12);  pL->SetBottomMargin(0.14);
  pL->SetTickx(1); pL->SetTicky(1); pL->SetGridx(); pL->SetGridy();
  pL->Draw();

  TPad* pR = new TPad("pRight", "", 0.50, 0.0, 1.00, 1.0);
  pR->SetLeftMargin(0.16); pR->SetRightMargin(0.05);
  pR->SetTopMargin(0.12);  pR->SetBottomMargin(0.14);
  pR->SetTickx(1); pR->SetTicky(1); pR->SetGridx(); pR->SetGridy();
  pR->Draw();

  // ── Left: M(K+K-) ─────────────────────────────────────────────────────────
  pL->cd();
  hX->Draw("HIST");

  // Signal window box
  TBox* sigBox = new TBox(sigLo, 0.0, sigHi, hX->GetMaximum()*1.30);
  sigBox->SetFillColorAlpha(kRed-4, 0.10f); sigBox->SetFillStyle(1001);
  sigBox->SetLineColor(kRed-4); sigBox->SetLineWidth(1); sigBox->SetLineStyle(3);
  sigBox->Draw("SAME L");
  hX->Draw("HIST SAME");  // redraw on top of box

  TLine lPhiKK(kMphi, 0.0, kMphi, hX->GetMaximum()*1.25);
  lPhiKK.SetLineColor(kRed+1); lPhiKK.SetLineWidth(2); lPhiKK.SetLineStyle(7);
  lPhiKK.Draw();

  TLatex lxKK; lxKK.SetNDC(false);
  lxKK.SetTextFont(43); lxKK.SetTextSize(13);
  lxKK.SetTextColor(kRed+1); lxKK.SetTextAngle(90);
  lxKK.DrawLatex(kMphi + 0.001, hX->GetMaximum()*0.08, "#phi(1020)");

  TLegend legL(0.40, 0.65, 0.93, 0.86);
  legL.SetBorderSize(0); legL.SetFillStyle(0);
  legL.SetTextFont(43);  legL.SetTextSize(14);
  legL.AddEntry(hX,     "M(K^{+}K^{-})",                          "LF");
  legL.AddEntry(&lPhiKK,"#phi(1020) = 1.0195 GeV",                "L");
  legL.AddEntry(sigBox, Form("Signal window #pm%.0f MeV",sigWin*1000), "F");
  legL.Draw();

  TLatex titL; titL.SetNDC();
  titL.SetTextFont(43); titL.SetTextSize(16); titL.SetTextAlign(21);
  titL.DrawLatex(0.57, 0.93, "M(K^{+}K^{-})  projection");

  // ── Right: MM(ep) — signal window only ───────────────────────────────────
  pR->cd();
  hYbkg->GetXaxis()->SetTitle(hYsig->GetXaxis()->GetTitle());
  hYbkg->GetXaxis()->SetTitleFont(43); hYbkg->GetXaxis()->SetTitleSize(16);
  hYbkg->GetXaxis()->SetTitleOffset(1.2);
  hYbkg->GetXaxis()->SetLabelFont(43); hYbkg->GetXaxis()->SetLabelSize(15);
  hYbkg->GetYaxis()->SetTitle("Norm. counts");
  hYbkg->GetYaxis()->SetTitleFont(43); hYbkg->GetYaxis()->SetTitleSize(18);
  hYbkg->GetYaxis()->SetTitleOffset(1.5);
  hYbkg->GetYaxis()->SetLabelFont(43); hYbkg->GetYaxis()->SetLabelSize(15);
  hYbkg->GetYaxis()->SetRangeUser(0.0, ymaxMM);
  hYbkg->GetYaxis()->SetNdivisions(504);
  hYbkg->Draw("HIST");
  hYsig->Draw("HIST SAME");

  TLine lPhiMM(kMphi, 0.0, kMphi, ymaxMM * 0.92);
  lPhiMM.SetLineColor(kRed+1); lPhiMM.SetLineWidth(2); lPhiMM.SetLineStyle(7);
  lPhiMM.Draw();

  TLatex lxMM; lxMM.SetNDC(false);
  lxMM.SetTextFont(43); lxMM.SetTextSize(13);
  lxMM.SetTextColor(kRed+1); lxMM.SetTextAngle(90);
  lxMM.DrawLatex(kMphi + 0.002, ymaxMM*0.05, "#phi(1020)");

  TLegend legR(0.40, 0.65, 0.93, 0.86);
  legR.SetBorderSize(0); legR.SetFillStyle(0);
  legR.SetTextFont(43);  legR.SetTextSize(14);
  legR.AddEntry(hYsig, Form("Signal: M(KK)#in[%.3f,%.3f]", sigLo, sigHi), "LF");
  legR.AddEntry(hYbkg, "Sideband (background)",                            "LF");
  legR.AddEntry(&lPhiMM, "#phi(1020) = 1.0195 GeV",                       "L");
  legR.Draw();

  TLatex titR; titR.SetNDC();
  titR.SetTextFont(43); titR.SetTextSize(14); titR.SetTextAlign(21);
  titR.DrawLatex(0.57, 0.93, "MM(ep)  [#phi signal window vs sideband]");

  c->Update();
  c->SaveAs(outFile.c_str());
  std::cout << "[PlotInvMassVsMissingMassProjections] Saved -> " << outFile << "\n";

  gStyle->SetOptStat(savedOptStat);
  gStyle->SetOptTitle(savedOptTitle);
}


// =============================================================================
// Internal colour palette matching DISANAcomparer::modelShades
// =============================================================================
namespace detail {
static const std::vector<std::tuple<float,float,float>> kModelRGB = {
    {0.20f, 0.30f, 0.85f},   // Blue
    {0.90f, 0.45f, 0.10f},   // Orange
    {0.00f, 0.60f, 0.60f},   // Teal
    {0.00f, 0.70f, 0.00f},   // Green
    {0.60f, 0.30f, 0.80f},   // Purple
    {0.85f, 0.10f, 0.25f},   // Red
    {0.40f, 0.40f, 0.40f},   // Gray
};
// Allocate or reuse a ROOT TColor and return its index.
inline int ModelColor(int im) {
    static const int base = 5200;
    const int idx = base + im * 17;   // 17-step stride avoids collisions
    auto [r, g, b] = kModelRGB[im % kModelRGB.size()];
    if (!gROOT->GetColor(idx)) new TColor(idx, r, g, b);
    return idx;
}
} // namespace detail


// =============================================================================
// PlotInvMassKKComparison
//
// 1D M(K+K-) per Q² bin, multiple datasets overlaid on the same canvas.
// All distributions normalised to unit area for shape comparison.
// A vertical dashed red line marks the phi(1020) nominal mass.
//
// Arguments:
//   dfs        – one RDataFrame per model/dataset (needs "Q2", "invMass_KpKm")
//   labels     – model labels, same order as dfs
//   q2Bins     – Q² bin-edge vector
//   means      – Q2BinMeans from ComputeMeanKinPerQ2Bin (first dataset)
//   outFile    – output PDF path
//   nBinsHist  – histogram bins (default 80)
//   xLo / xHi – M(K+K-) range [GeV]  (default 0.98 / 1.10)
// =============================================================================
inline void PlotInvMassKKComparison(
    std::vector<ROOT::RDF::RNode> dfs,
    const std::vector<std::string>&      labels,
    const std::vector<double>&           q2Bins,
    const Q2BinMeans&                    means,
    const std::string&                   outFile   = "invmassKK_comparison.pdf",
    int    nBinsHist = 80,
    double xLo       = 0.98,
    double xHi       = 1.10)
{
  if (dfs.empty() || dfs.size() != labels.size()) {
    std::cerr << "[PlotInvMassKKComparison] ERROR: dfs and labels must be non-empty and same size.\n";
    return;
  }
  const int nBins = static_cast<int>(q2Bins.size()) - 1;
  if (nBins <= 0) {
    std::cerr << "[PlotInvMassKKComparison] ERROR: q2Bins must have >= 2 entries.\n";
    return;
  }
  const int  nModels = static_cast<int>(dfs.size());
  const bool q2_ok   = (static_cast<int>(means.Q2.size()) == nBins);
  static constexpr double kMphi_nom = 1.01946;

  // ── Global style (scoped) ────────────────────────────────────────────────
  const int savedOptStat  = gStyle->GetOptStat();
  const int savedOptTitle = gStyle->GetOptTitle();
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetNdivisions(505, "X"); gStyle->SetNdivisions(504, "Y");

  // ── Book all histograms lazily ────────────────────────────────────────────
  using H1Ptr = ROOT::RDF::RResultPtr<TH1D>;
  std::vector<std::vector<H1Ptr>> hb(nModels, std::vector<H1Ptr>(nBins));
  std::vector<H1Ptr>              hInt(nModels);

  for (int im = 0; im < nModels; ++im) {
    for (int iq = 0; iq < nBins; ++iq) {
      const double lo = q2Bins[iq], hi = q2Bins[iq + 1];
      auto df_bin = dfs[im].Filter(
          [lo, hi](double Q2){ return Q2 >= lo && Q2 < hi; }, {"Q2"},
          Form("invkk_m%d_b%d", im, iq));
      hb[im][iq] = df_bin.Histo1D(
          {Form("hKK_m%d_b%d", im, iq), "", nBinsHist, xLo, xHi},
          "invMass_KpKm");
    }
    hInt[im] = dfs[im].Histo1D(
        {Form("hKK_m%d_int", im), "", nBinsHist, xLo, xHi},
        "invMass_KpKm");
    hInt[im].GetValue();   // trigger lazy evaluation
  }

  // ── Canvas geometry ───────────────────────────────────────────────────────
  const int    nPads   = nBins + 1;
  const int    nCols   = (int)std::ceil(std::sqrt((double)nPads));
  const int    nRows   = (int)std::ceil((double)nPads / nCols);
  const double legFrac = 0.08, padTop = 1.0 - legFrac;
  const int cW = std::max(700, 400 * nCols);
  const int cH = (int)(std::max(550, 370 * nRows) / (1.0 - legFrac));

  TCanvas* c = new TCanvas("cInvMassKKComp", "M(K^{+}K^{-}) comparison", cW, cH);
  c->SetFillColorAlpha(0, 0);

  // ── Legend strip ─────────────────────────────────────────────────────────
  c->cd();
  TPad* pLeg = new TPad("pLegKK", "", 0.0, padTop, 1.0, 1.0);
  pLeg->SetFillColorAlpha(0, 0); pLeg->SetBorderSize(0);
  pLeg->Draw(); pLeg->cd();

  auto* leg = new TLegend(0.02, 0.05, 0.98, 0.95);
  leg->SetNColumns(std::min(nModels + 1, 4));
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(43);  leg->SetTextSize(18);
  for (int im = 0; im < nModels; ++im) {
    const int col = detail::ModelColor(im);
    auto* d = new TH1F(Form("_dlegKK_%d", im), "", 1, 0, 1);
    d->SetLineColor(col); d->SetLineWidth(3); d->SetLineStyle(1 + im % 5);
    d->SetFillColorAlpha(col, 0.0f);
    leg->AddEntry(d, labels[im].c_str(), "L");
  }
  auto* dPhi = new TH1F("_dlegKK_phi", "", 1, 0, 1);
  dPhi->SetLineColor(kRed + 1); dPhi->SetLineWidth(2); dPhi->SetLineStyle(7);
  leg->AddEntry(dPhi, "#phi(1020) = 1.0195 GeV", "L");
  leg->Draw();

  // ── drawCell lambda ───────────────────────────────────────────────────────
  auto drawCell = [&](int iPad, std::vector<TH1D*>& vH, const std::string& title)
  {
    const double x0 = (double)(iPad % nCols)       / nCols;
    const double x1 = (double)(iPad % nCols + 1)   / nCols;
    const double y1 = padTop - (double)(iPad / nCols)     / nRows * padTop;
    const double y0 = padTop - (double)(iPad / nCols + 1) / nRows * padTop;

    c->cd();
    TPad* pO = new TPad(Form("pOKK%d", iPad), "", x0, y0, x1, y1);
    pO->SetFillColorAlpha(0,0); pO->SetBorderSize(0);
    pO->SetLeftMargin(0); pO->SetRightMargin(0);
    pO->SetTopMargin(0);  pO->SetBottomMargin(0);
    pO->Draw(); pO->cd();

    TPad* pM = new TPad(Form("pMKK%d", iPad), "", 0.0, 0.0, 1.0, 1.0);
    pM->SetFillColorAlpha(0,0); pM->SetBorderSize(0);
    pM->SetLeftMargin(0.20); pM->SetRightMargin(0.04);
    pM->SetTopMargin(0.18);  pM->SetBottomMargin(0.20);
    pM->SetTickx(1); pM->SetTicky(1); pM->SetGridx(); pM->SetGridy();
    pM->Draw(); pM->cd();

    double ymax = 0.0;
    for (auto* h : vH) if (h) ymax = std::max(ymax, h->GetMaximum());
    ymax *= 1.35;

    bool first = true;
    for (int im = 0; im < nModels; ++im) {
      TH1D* h = vH[im]; if (!h) continue;
      const int col = detail::ModelColor(im);
      h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(1 + im % 5);
      h->SetFillStyle(0);
      h->GetXaxis()->SetTitle("M(K^{+}K^{-})  [GeV]");
      h->GetXaxis()->SetTitleFont(43); h->GetXaxis()->SetTitleSize(13);
      h->GetXaxis()->SetTitleOffset(1.7);
      h->GetXaxis()->SetLabelFont(43); h->GetXaxis()->SetLabelSize(11);
      h->GetXaxis()->SetTickLength(0.05);
      h->GetYaxis()->SetTitle("Norm. counts");
      h->GetYaxis()->SetTitleFont(43); h->GetYaxis()->SetTitleSize(13);
      h->GetYaxis()->SetTitleOffset(1.9);
      h->GetYaxis()->SetLabelFont(43); h->GetYaxis()->SetLabelSize(11);
      h->GetYaxis()->SetRangeUser(0.0, ymax);
      h->GetYaxis()->SetNdivisions(504);
      h->Draw(first ? "HIST" : "HIST SAME");
      first = false;
    }
    TLine* lPhi = new TLine(kMphi_nom, 0.0, kMphi_nom, ymax);
    lPhi->SetLineColor(kRed+1); lPhi->SetLineWidth(2); lPhi->SetLineStyle(7);
    lPhi->Draw();

    TPaveText* pt = new TPaveText(0.21, 0.82, 0.97, 0.97, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.65); pt->SetBorderSize(0);
    pt->SetTextFont(43); pt->SetTextSize(12); pt->SetTextAlign(22);
    pt->AddText(title.c_str()); pt->Draw();
  };

  // ── Per-Q² cells ─────────────────────────────────────────────────────────
  for (int iq = 0; iq < nBins; ++iq) {
    std::vector<TH1D*> vH(nModels, nullptr);
    for (int im = 0; im < nModels; ++im) {
      TH1D* h = (TH1D*)hb[im][iq].GetPtr()->Clone(Form("cKK_m%d_b%d", im, iq));
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      vH[im] = h;
    }
    const double q2m = (q2_ok && means.Q2[iq] > 0.0)
                       ? means.Q2[iq] : 0.5*(q2Bins[iq]+q2Bins[iq+1]);
    drawCell(iq, vH, Form("Q^{2} #in [%.2f, %.2f],  #bar{Q}^{2} = %.3f",
                          q2Bins[iq], q2Bins[iq+1], q2m));
  }
  // ── Integrated cell ───────────────────────────────────────────────────────
  {
    std::vector<TH1D*> vH(nModels, nullptr);
    for (int im = 0; im < nModels; ++im) {
      TH1D* h = (TH1D*)hInt[im].GetPtr()->Clone(Form("cKK_m%d_int", im));
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      vH[im] = h;
    }
    drawCell(nBins, vH, Form("All Q^{2} #in [%.2f, %.2f]",
                              q2Bins.front(), q2Bins.back()));
  }

  c->cd(); c->Update();
  c->SaveAs(outFile.c_str());
  std::cout << "[PlotInvMassKKComparison] Saved -> " << outFile << "\n";
  gStyle->SetOptStat(savedOptStat); gStyle->SetOptTitle(savedOptTitle);
}


// =============================================================================
// PlotMissingMassEpComparison
//
// 1D MM(ep) = sqrt(Mx2_ep) per Q² bin, multiple datasets overlaid.
// Same canvas layout and normalisation as PlotInvMassKKComparison.
// A vertical dashed blue line marks the proton mass M_p = 0.938 GeV.
//
// Arguments:
//   dfs_in     – one RDataFrame per model (needs "Q2" and "Mx2_ep" or "MM_ep")
//   labels     – model labels, same order as dfs_in
//   q2Bins     – Q² bin-edge vector
//   means      – Q2BinMeans from ComputeMeanKinPerQ2Bin (first dataset)
//   outFile    – output PDF path
//   nBinsHist  – histogram bins (default 80)
//   xLo / xHi – MM(ep) range [GeV]  (default 0.70 / 1.20)
// =============================================================================
inline void PlotMissingMassEpComparison(
    std::vector<ROOT::RDF::RNode> dfs_in,
    const std::vector<std::string>&      labels,
    const std::vector<double>&           q2Bins,
    const Q2BinMeans&                    means,
    const std::string&                   outFile   = "missmassEp_comparison.pdf",
    int    nBinsHist = 80,
    double xLo       = 0.70,
    double xHi       = 1.20)
{
  if (dfs_in.empty() || dfs_in.size() != labels.size()) {
    std::cerr << "[PlotMissingMassEpComparison] ERROR: dfs and labels must be non-empty and same size.\n";
    return;
  }
  const int nBins   = static_cast<int>(q2Bins.size()) - 1;
  if (nBins <= 0) {
    std::cerr << "[PlotMissingMassEpComparison] ERROR: q2Bins must have >= 2 entries.\n";
    return;
  }
  const int  nModels = static_cast<int>(dfs_in.size());
  const bool q2_ok   = (static_cast<int>(means.Q2.size()) == nBins);

  // ── Ensure MM_ep column exists on every RDF ───────────────────────────────
  std::vector<ROOT::RDF::RNode> dfs;
  dfs.reserve(nModels);
  for (int im = 0; im < nModels; ++im) {
    auto cols  = dfs_in[im].GetColumnNames();
    bool hasMM = std::find(cols.begin(), cols.end(), "MM_ep")  != cols.end();
    bool hasMx = std::find(cols.begin(), cols.end(), "Mx2_ep") != cols.end();
    if (!hasMM && !hasMx) {
      std::cerr << "[PlotMissingMassEpComparison] ERROR model " << im
                << ": neither MM_ep nor Mx2_ep found.\n";
      return;
    }
    dfs.push_back(hasMM ? dfs_in[im]
                        : dfs_in[im].Define("MM_ep", [](double Mx2) -> double {
                              return (Mx2 > 0.0) ? std::sqrt(Mx2) : -999.0;
                          }, {"Mx2_ep"}));
  }

  // ── Global style (scoped) ────────────────────────────────────────────────
  const int savedOptStat  = gStyle->GetOptStat();
  const int savedOptTitle = gStyle->GetOptTitle();
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetNdivisions(505, "X"); gStyle->SetNdivisions(504, "Y");

  // ── Book histograms lazily ────────────────────────────────────────────────
  using H1Ptr = ROOT::RDF::RResultPtr<TH1D>;
  std::vector<std::vector<H1Ptr>> hb(nModels, std::vector<H1Ptr>(nBins));
  std::vector<H1Ptr>              hInt(nModels);

  for (int im = 0; im < nModels; ++im) {
    for (int iq = 0; iq < nBins; ++iq) {
      const double lo = q2Bins[iq], hi = q2Bins[iq + 1];
      auto df_bin = dfs[im].Filter(
          [lo, hi](double Q2){ return Q2 >= lo && Q2 < hi; }, {"Q2"},
          Form("mmep_m%d_b%d", im, iq));
      hb[im][iq] = df_bin.Histo1D(
          {Form("hMMep_m%d_b%d", im, iq), "", nBinsHist, xLo, xHi},
          "MM_ep");
    }
    hInt[im] = dfs[im].Histo1D(
        {Form("hMMep_m%d_int", im), "", nBinsHist, xLo, xHi},
        "MM_ep");
    hInt[im].GetValue();
  }

  // ── Canvas geometry ───────────────────────────────────────────────────────
  const int    nPads   = nBins + 1;
  const int    nCols   = (int)std::ceil(std::sqrt((double)nPads));
  const int    nRows   = (int)std::ceil((double)nPads / nCols);
  const double legFrac = 0.08, padTop = 1.0 - legFrac;
  const int cW = std::max(700, 400 * nCols);
  const int cH = (int)(std::max(550, 370 * nRows) / (1.0 - legFrac));

  TCanvas* c = new TCanvas("cMMepComp", "MM(ep) comparison", cW, cH);
  c->SetFillColorAlpha(0, 0);

  // ── Legend strip ─────────────────────────────────────────────────────────
  c->cd();
  TPad* pLeg = new TPad("pLegMMep", "", 0.0, padTop, 1.0, 1.0);
  pLeg->SetFillColorAlpha(0, 0); pLeg->SetBorderSize(0);
  pLeg->Draw(); pLeg->cd();

  auto* leg = new TLegend(0.02, 0.05, 0.98, 0.95);
  leg->SetNColumns(std::min(nModels + 1, 4));
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(43);  leg->SetTextSize(18);
  for (int im = 0; im < nModels; ++im) {
    const int col = detail::ModelColor(im);
    auto* d = new TH1F(Form("_dlegMM_%d", im), "", 1, 0, 1);
    d->SetLineColor(col); d->SetLineWidth(3); d->SetLineStyle(1 + im % 5);
    d->SetFillColorAlpha(col, 0.0f);
    leg->AddEntry(d, labels[im].c_str(), "L");
  }
  auto* dPro = new TH1F("_dlegMM_pro", "", 1, 0, 1);
  dPro->SetLineColor(kAzure - 4); dPro->SetLineWidth(2); dPro->SetLineStyle(7);
  leg->AddEntry(dPro, Form("M_{p} = %.3f GeV", kMp), "L");
  leg->Draw();

  // ── drawCell lambda ───────────────────────────────────────────────────────
  auto drawCell = [&](int iPad, std::vector<TH1D*>& vH, const std::string& title)
  {
    const double x0 = (double)(iPad % nCols)       / nCols;
    const double x1 = (double)(iPad % nCols + 1)   / nCols;
    const double y1 = padTop - (double)(iPad / nCols)     / nRows * padTop;
    const double y0 = padTop - (double)(iPad / nCols + 1) / nRows * padTop;

    c->cd();
    TPad* pO = new TPad(Form("pOMM%d", iPad), "", x0, y0, x1, y1);
    pO->SetFillColorAlpha(0,0); pO->SetBorderSize(0);
    pO->SetLeftMargin(0); pO->SetRightMargin(0);
    pO->SetTopMargin(0);  pO->SetBottomMargin(0);
    pO->Draw(); pO->cd();

    TPad* pM = new TPad(Form("pMMM%d", iPad), "", 0.0, 0.0, 1.0, 1.0);
    pM->SetFillColorAlpha(0,0); pM->SetBorderSize(0);
    pM->SetLeftMargin(0.20); pM->SetRightMargin(0.04);
    pM->SetTopMargin(0.18);  pM->SetBottomMargin(0.20);
    pM->SetTickx(1); pM->SetTicky(1); pM->SetGridx(); pM->SetGridy();
    pM->Draw(); pM->cd();

    double ymax = 0.0;
    for (auto* h : vH) if (h) ymax = std::max(ymax, h->GetMaximum());
    ymax *= 1.35;

    bool first = true;
    for (int im = 0; im < nModels; ++im) {
      TH1D* h = vH[im]; if (!h) continue;
      const int col = detail::ModelColor(im);
      h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(1 + im % 5);
      h->SetFillStyle(0);
      h->GetXaxis()->SetTitle("MM(ep)  [GeV]");
      h->GetXaxis()->SetTitleFont(43); h->GetXaxis()->SetTitleSize(13);
      h->GetXaxis()->SetTitleOffset(1.7);
      h->GetXaxis()->SetLabelFont(43); h->GetXaxis()->SetLabelSize(11);
      h->GetXaxis()->SetTickLength(0.05);
      h->GetYaxis()->SetTitle("Norm. counts");
      h->GetYaxis()->SetTitleFont(43); h->GetYaxis()->SetTitleSize(13);
      h->GetYaxis()->SetTitleOffset(1.9);
      h->GetYaxis()->SetLabelFont(43); h->GetYaxis()->SetLabelSize(11);
      h->GetYaxis()->SetRangeUser(0.0, ymax);
      h->GetYaxis()->SetNdivisions(504);
      h->Draw(first ? "HIST" : "HIST SAME");
      first = false;
    }
    TLine* lPro = new TLine(kMp, 0.0, kMp, ymax);
    lPro->SetLineColor(kAzure-4); lPro->SetLineWidth(2); lPro->SetLineStyle(7);
    lPro->Draw();

    TPaveText* pt = new TPaveText(0.21, 0.82, 0.97, 0.97, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.65); pt->SetBorderSize(0);
    pt->SetTextFont(43); pt->SetTextSize(12); pt->SetTextAlign(22);
    pt->AddText(title.c_str()); pt->Draw();
  };

  // ── Per-Q² cells ─────────────────────────────────────────────────────────
  for (int iq = 0; iq < nBins; ++iq) {
    std::vector<TH1D*> vH(nModels, nullptr);
    for (int im = 0; im < nModels; ++im) {
      TH1D* h = (TH1D*)hb[im][iq].GetPtr()->Clone(Form("cMM_m%d_b%d", im, iq));
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      vH[im] = h;
    }
    const double q2m = (q2_ok && means.Q2[iq] > 0.0)
                       ? means.Q2[iq] : 0.5*(q2Bins[iq]+q2Bins[iq+1]);
    drawCell(iq, vH, Form("Q^{2} #in [%.2f, %.2f],  #bar{Q}^{2} = %.3f",
                          q2Bins[iq], q2Bins[iq+1], q2m));
  }
  // ── Integrated cell ───────────────────────────────────────────────────────
  {
    std::vector<TH1D*> vH(nModels, nullptr);
    for (int im = 0; im < nModels; ++im) {
      TH1D* h = (TH1D*)hInt[im].GetPtr()->Clone(Form("cMM_m%d_int", im));
      h->SetDirectory(nullptr);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      vH[im] = h;
    }
    drawCell(nBins, vH, Form("All Q^{2} #in [%.2f, %.2f]",
                              q2Bins.front(), q2Bins.back()));
  }

  c->cd(); c->Update();
  c->SaveAs(outFile.c_str());
  std::cout << "[PlotMissingMassEpComparison] Saved -> " << outFile << "\n";
  gStyle->SetOptStat(savedOptStat); gStyle->SetOptTitle(savedOptTitle);
}


inline ROOT::RDF::RNode InitKinematics_ExclusiveKp(const std::string& f, const std::string& t, float E) {
  return InitKinematics_MissingKm(f, t, E);  // exclusive K+ == K- omitted
}
inline ROOT::RDF::RNode InitKinematics_ExclusiveKm(const std::string& f, const std::string& t, float E) {
  return InitKinematics_MissingKp(f, t, E);  // exclusive K- == K+ omitted
}

/// Data
ROOT::RDF::RNode InitKinematics(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  *df_ = df_->Define("ele_px_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recel_vz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("reckMinus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kMinus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kMinus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kMinus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("reckPlus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kPlus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kPlus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kPlus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})

             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
             .Define("RunNumber", "RUN_config_run")
             .Define("EventNumber", "RUN_config_event")
             .Define("nElectrons",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
                       int n = 0;
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) ++n;
                       return n;
                     },
                     {"REC_Particle_pid", "REC_Particle_pass"})
             // Choose ONE best electron per event
             .Define("bestEle_idx", ChooseBestElectron, {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz", "REC_Particle_vz", "REC_Particle_pass"})

             // If you want to **keep** multi-e events but only use the best e⁻, leave the next line.
             // If you want to **drop** multi-e events (strict DVEP), add: .Filter("nElectrons==1","Cut: exactly 1 e⁻")
             .Filter("bestEle_idx >= 0", "At least one usable e⁻")
             .Define("REC_pass_bestE",
                     [](int best, const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& pass) {
                       ROOT::VecOps::RVec<int> keep(pass.size(), 0);
                       if (best >= 0) keep[best] = 1;        // only best e⁻
                       for (size_t i = 0; i < pid.size(); ++i)  // keep non-e tracks as they were
                         if (pid[i] != 11) keep[i] = pass[i];
                       return keep;
                     },
                     {"bestEle_idx", "REC_Particle_pid", "REC_Particle_pass"})

             .Define("nElectrons_best",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<int>& keep) {
                       int n = 0;
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && keep[i]) ++n;
                       return n;
                     },
                     {"REC_Particle_pid", "REC_pass_bestE"})
             // Project the chosen e⁻ to the scalar columns used everywhere else
             .Define("ele_px", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_px"})
             .Define("ele_py", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_py"})
             .Define("ele_pz", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_pz"})
             .Define("recel_vz", [](int i, const ROOT::VecOps::RVec<float>& v) { return i >= 0 ? v[i] : -999.f; }, {"bestEle_idx", "REC_Particle_vz"})
             .Define("ele_det_region", [](int i, const ROOT::VecOps::RVec<short>& st) { return i >= 0 ? DetRegionFromStatus(st[i]) : -1; }, {"bestEle_idx", "REC_Particle_status"})

             .Filter([](float ex, float kMinusx, float kPlusx, float px) { return ex != -999 && kMinusx != -999 && kPlusx != -999 && px != -999; },
                     {"ele_px", "kMinus_px", "kPlus_px", "pro_px"})
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("reckMinus_p", MomentumFunc, {"kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("reckMinus_theta", ThetaFunc, {"kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("reckMinus_phi", PhiFunc, {"kMinus_px", "kMinus_py"})
             .Define("reckPlus_p", MomentumFunc, {"kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("reckPlus_theta", ThetaFunc, {"kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("reckPlus_phi", PhiFunc, {"kPlus_px", "kPlus_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("kMinus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == -321 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("kPlus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 321 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})

             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("ele_det_region_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<int>& pass) -> int {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})

             .Define("invMass_KpKm",
                     [](float px1, float py1, float pz1, float px2, float py2, float pz2) -> float {
                       constexpr float mK = 0.493677;  // Kaon mass in GeV/c²
                       float E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mK * mK);
                       float E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mK * mK);
                       float px = px1 + px2;
                       float py = py1 + py2;
                       float pz = pz1 + pz2;
                       float E = E1 + E2;
                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"kPlus_px", "kPlus_py", "kPlus_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})
             // p–K⁻ invariant mass
             .Define("invMass_pKminus",
                     [](float p_px, float p_py, float p_pz, float k_px, float k_py, float k_pz) -> float {
                       constexpr float mP = 0.938272f;  // Proton mass in GeV/c²
                       constexpr float mK = 0.493677f;  // Kaon mass in GeV/c²

                       float E_p = std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz + mP * mP);
                       float E_k = std::sqrt(k_px * k_px + k_py * k_py + k_pz * k_pz + mK * mK);

                       float px = p_px + k_px;
                       float py = p_py + k_py;
                       float pz = p_pz + k_pz;
                       float E = E_p + E_k;

                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"pro_px", "pro_py", "pro_pz", "kMinus_px", "kMinus_py", "kMinus_pz"})

             // p–K⁺ invariant mass
             .Define("invMass_pKplus",
                     [](float p_px, float p_py, float p_pz, float k_px, float k_py, float k_pz) -> float {
                       constexpr float mP = 0.938272f;  // Proton mass in GeV/c²
                       constexpr float mK = 0.493677f;  // Kaon mass in GeV/c²

                       float E_p = std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz + mP * mP);
                       float E_k = std::sqrt(k_px * k_px + k_py * k_py + k_pz * k_pz + mK * mK);

                       float px = p_px + k_px;
                       float py = p_py + k_py;
                       float pz = p_pz + k_pz;
                       float E = E_p + E_k;

                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"pro_px", "pro_py", "pro_pz", "kPlus_px", "kPlus_py", "kPlus_pz"});

  // DISANAMath-driven observables
  *df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "tmin", &DISANAMath::GetTmin, beam_energy);

  *df_ = df_->Define("mtprime",  // non-negative: this is what you’ll plot
                     [](double mt, double tmin) { return std::abs(mt + tmin); }, {"t", "tmin"})
             .Define("tprime",  // optional signed t' ≤ 0
                     [](double mtp) { return -mtp; }, {"mtprime"});

  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);
  *df_ = define_DISCAT(*df_, "cos_thetaKK", &DISANAMath::GetCosTheta_KK, beam_energy);
  *df_ = define_DISCAT(*df_, "cos_phiKK", &DISANAMath::GetCosPhi_KK, beam_energy);

  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKp", &DISANAMath::GetMx2_eKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx_eKp", &DISANAMath::GetMx_eKp, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_p", &DISANAMath::GetCone_p, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Kp", &DISANAMath::GetCone_Kp, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Km", &DISANAMath::GetCone_Km, beam_energy);
  *df_ = define_DISCAT(*df_, "z_phi", &DISANAMath::GetZ_phi, beam_energy);
  *df_ = define_DISCAT(*df_, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  *df_ = df_->Define("Gamma_v", [beam_energy](double Q2, double W) {
                 return HandGammaV(beam_energy, Q2, W);
               }, {"Q2", "W"});
  // --- new: vm_cut and Egamma_star -----------------------------------------
  *df_ = df_->Define("vm_cut",
    [beam_energy](double Q2, double xB, double t) -> double {
      return DISANAMath::ComputeVmApprox(Q2, xB, t, m_phi, beam_energy);
    }, {"Q2", "xB", "t"})
  .Define("Egamma_star",
    [](double Mx2) -> double {
      if (Mx2 <= 0.0) return -999.0;
      const double Mx = std::sqrt(Mx2);
      return (Mx2 - kMp * kMp) / (2.0 * kMp);
    }, {"Mx2_eKpKm"});

  return *df_;
}

// MC
ROOT::RDF::RNode InitGenKinematics(const std::string& filename_,
                                  const std::string& treename_,
                                  float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  ROOT::RDF::RNode df = rdf;

  // pick "best" particle of a given PID = highest momentum
  auto bestIdxByPid = [](int targetPid,
                         const ROOT::VecOps::RVec<int>& pid,
                         const ROOT::VecOps::RVec<float>& px,
                         const ROOT::VecOps::RVec<float>& py,
                         const ROOT::VecOps::RVec<float>& pz) -> int {
    int best = -1;
    float bestP2 = -1.f;
    for (size_t i = 0; i < pid.size(); ++i) {
      if (pid[i] != targetPid) continue;
      float p2 = px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
      if (p2 > bestP2) { bestP2 = p2; best = (int)i; }
    }
    return best;
  };

  // safe component getter
  auto atOr = [](int idx, const ROOT::VecOps::RVec<float>& v) -> float {
    return (idx >= 0 && (size_t)idx < v.size()) ? v[idx] : -999.0f;
  };

df = df
  // indices
  .Define("gen_e_idx",
          [bestIdxByPid](const RVec<int>& pid,
                         const RVec<float>& px,
                         const RVec<float>& py,
                         const RVec<float>& pz) {
            return bestIdxByPid(11, pid, px, py, pz);
          },
          {"MC_Particle_pid","MC_Particle_px","MC_Particle_py","MC_Particle_pz"})

  .Define("gen_p_idx",
          [bestIdxByPid](const RVec<int>& pid,
                         const RVec<float>& px,
                         const RVec<float>& py,
                         const RVec<float>& pz) {
            return bestIdxByPid(2212, pid, px, py, pz);
          },
          {"MC_Particle_pid","MC_Particle_px","MC_Particle_py","MC_Particle_pz"})

  .Define("gen_kp_idx",
          [bestIdxByPid](const RVec<int>& pid,
                         const RVec<float>& px,
                         const RVec<float>& py,
                         const RVec<float>& pz) {
            return bestIdxByPid(321, pid, px, py, pz);
          },
          {"MC_Particle_pid","MC_Particle_px","MC_Particle_py","MC_Particle_pz"})

  .Define("gen_km_idx",
          [bestIdxByPid](const RVec<int>& pid,
                         const RVec<float>& px,
                         const RVec<float>& py,
                         const RVec<float>& pz) {
            return bestIdxByPid(-321, pid, px, py, pz);
          },
          {"MC_Particle_pid","MC_Particle_px","MC_Particle_py","MC_Particle_pz"})

  // components (also avoid `auto` here)
  .Define("ele_px",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_e_idx","MC_Particle_px"})
  .Define("ele_py",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_e_idx","MC_Particle_py"})
  .Define("ele_pz",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_e_idx","MC_Particle_pz"})

  .Define("pro_px",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_p_idx","MC_Particle_px"})
  .Define("pro_py",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_p_idx","MC_Particle_py"})
  .Define("pro_pz",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_p_idx","MC_Particle_pz"})

  .Define("kPlus_px",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_kp_idx","MC_Particle_px"})
  .Define("kPlus_py",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_kp_idx","MC_Particle_py"})
  .Define("kPlus_pz",  [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_kp_idx","MC_Particle_pz"})

  .Define("kMinus_px", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_km_idx","MC_Particle_px"})
  .Define("kMinus_py", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_km_idx","MC_Particle_py"})
  .Define("kMinus_pz", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_km_idx","MC_Particle_pz"})
  .Define("recel_vz", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_e_idx","MC_Particle_pz"})
  .Define("recpro_vz", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_p_idx","MC_Particle_pz"})
  .Define("reckPlus_vz", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_kp_idx","MC_Particle_pz"})
  .Define("reckMinus_vz", [atOr](int i, const RVec<float>& v){ return atOr(i, v); }, {"gen_km_idx","MC_Particle_pz"})

    // require all 4 exist (keep if you want "all-particle" GEN only)
    .Filter([](int ei, int pi, int kpi, int kmi){
        return (ei >= 0 && pi >= 0 && kpi >= 0 && kmi >= 0);
      }, {"gen_e_idx","gen_p_idx","gen_kp_idx","gen_km_idx"})

    // additionally guard against -999 components
    .Filter([](float ex, float kMx, float kPx, float px){
        return ex != -999.f && kMx != -999.f && kPx != -999.f && px != -999.f;
      }, {"ele_px","kMinus_px","kPlus_px","pro_px"})

    // rec* kinematics
    .Define("recel_p",     MomentumFunc, {"ele_px","ele_py","ele_pz"})
    .Define("recel_theta", ThetaFunc,    {"ele_px","ele_py","ele_pz"})
    .Define("recel_phi",   PhiFunc,      {"ele_px","ele_py"})

    .Define("recpro_p",     MomentumFunc, {"pro_px","pro_py","pro_pz"})
    .Define("recpro_theta", ThetaFunc,    {"pro_px","pro_py","pro_pz"})
    .Define("recpro_phi",   PhiFunc,      {"pro_px","pro_py"})

    .Define("reckPlus_p",     MomentumFunc, {"kPlus_px","kPlus_py","kPlus_pz"})
    .Define("reckPlus_theta", ThetaFunc,    {"kPlus_px","kPlus_py","kPlus_pz"})
    .Define("reckPlus_phi",   PhiFunc,      {"kPlus_px","kPlus_py"})

    .Define("reckMinus_p",     MomentumFunc, {"kMinus_px","kMinus_py","kMinus_pz"})
    .Define("reckMinus_theta", ThetaFunc,    {"kMinus_px","kMinus_py","kMinus_pz"})
    .Define("reckMinus_phi",   PhiFunc,      {"kMinus_px","kMinus_py"})

    // "regions" placeholders (truth has no region unless you map it)
    .Define("kPlus_det_region",
            [](const ROOT::VecOps::RVec<int>& pid){
              for (size_t i=0;i<pid.size();++i) if (pid[i]==321) return 1;
              return -1;
            }, {"MC_Particle_pid"})
    .Define("kMinus_det_region",
            [](const ROOT::VecOps::RVec<int>& pid){
              for (size_t i=0;i<pid.size();++i) if (pid[i]==-321) return 1; // <-- FIXED SIGN
              return -1;
            }, {"MC_Particle_pid"})
    .Define("pro_det_region",
            [](const ROOT::VecOps::RVec<int>& pid){
              for (size_t i=0;i<pid.size();++i) if (pid[i]==2212) return 1;
              return -1;
            }, {"MC_Particle_pid"})
    .Define("ele_det_region",
            [](const ROOT::VecOps::RVec<int>& pid){
              for (size_t i=0;i<pid.size();++i) if (pid[i]==11) return 1;
              return -1;
            }, {"MC_Particle_pid"})

    // invariant masses (with numerical guard)
    .Define("invMass_KpKm",
            [](float px1,float py1,float pz1,float px2,float py2,float pz2){
              constexpr float mK = 0.493677f;
              float E1 = std::sqrt(px1*px1 + py1*py1 + pz1*pz1 + mK*mK);
              float E2 = std::sqrt(px2*px2 + py2*py2 + pz2*pz2 + mK*mK);
              float px = px1 + px2, py = py1 + py2, pz = pz1 + pz2, E = E1 + E2;
              float m2 = E*E - (px*px + py*py + pz*pz);
              return (m2 > 0.f) ? std::sqrt(m2) : 0.f;
            }, {"kPlus_px","kPlus_py","kPlus_pz","kMinus_px","kMinus_py","kMinus_pz"})

    .Define("invMass_pKminus",
            [](float ppx,float ppy,float ppz,float kpx,float kpy,float kpz){
              constexpr float mP = 0.938272f, mK = 0.493677f;
              float Ep = std::sqrt(ppx*ppx + ppy*ppy + ppz*ppz + mP*mP);
              float Ek = std::sqrt(kpx*kpx + kpy*kpy + kpz*kpz + mK*mK);
              float px = ppx + kpx, py = ppy + kpy, pz = ppz + kpz, E = Ep + Ek;
              float m2 = E*E - (px*px + py*py + pz*pz);
              return (m2 > 0.f) ? std::sqrt(m2) : 0.f;
            }, {"pro_px","pro_py","pro_pz","kMinus_px","kMinus_py","kMinus_pz"})

    .Define("invMass_pKplus",
            [](float ppx,float ppy,float ppz,float kpx,float kpy,float kpz){
              constexpr float mP = 0.938272f, mK = 0.493677f;
              float Ep = std::sqrt(ppx*ppx + ppy*ppy + ppz*ppz + mP*mP);
              float Ek = std::sqrt(kpx*kpx + kpy*kpy + kpz*kpz + mK*mK);
              float px = ppx + kpx, py = ppy + kpy, pz = ppz + kpz, E = Ep + Ek;
              float m2 = E*E - (px*px + py*py + pz*pz);
              return (m2 > 0.f) ? std::sqrt(m2) : 0.f;
            }, {"pro_px","pro_py","pro_pz","kPlus_px","kPlus_py","kPlus_pz"});

  // DISANAMath-driven observables
  df = define_DISCAT(df, "Q2",  &DISANAMath::GetQ2,  beam_energy);
  df = define_DISCAT(df, "xB",  &DISANAMath::GetxB,  beam_energy);
  df = define_DISCAT(df, "t",   &DISANAMath::GetT,   beam_energy);
  df = define_DISCAT(df, "tmin",&DISANAMath::GetTmin,beam_energy);

  df = df.Define("mtprime", [](double t, double tmin){ return std::abs(t + tmin); }, {"t","tmin"})
         .Define("tprime",  [](double mtp){ return -mtp; }, {"mtprime"});

  df = define_DISCAT(df, "phi", &DISANAMath::GetPhi, beam_energy);
  df = define_DISCAT(df, "W",   &DISANAMath::GetW,   beam_energy);
  df = define_DISCAT(df, "nu",  &DISANAMath::GetNu,  beam_energy);
  df = define_DISCAT(df, "y",   &DISANAMath::Gety,   beam_energy);

  df = define_DISCAT(df, "cos_thetaKK", &DISANAMath::GetCosTheta_KK, beam_energy);
  df = define_DISCAT(df, "cos_phiKK",   &DISANAMath::GetCosPhi_KK,   beam_energy);

  df = define_DISCAT(df, "Mx2_ep",       &DISANAMath::GetMx2_ep,       beam_energy);
  df = define_DISCAT(df, "Emiss",        &DISANAMath::GetEmiss,        beam_energy);
  df = define_DISCAT(df, "PTmiss",       &DISANAMath::GetPTmiss,       beam_energy);
  df = define_DISCAT(df, "Mx2_epKpKm",   &DISANAMath::GetMx2_epKpKm,   beam_energy);
  df = define_DISCAT(df, "Mx2_eKpKm",    &DISANAMath::GetMx2_eKpKm,    beam_energy);
  df = define_DISCAT(df, "Mx2_epKm",     &DISANAMath::GetMx2_epKm,     beam_energy);
  df = define_DISCAT(df, "Mx2_epKp",     &DISANAMath::GetMx2_epKp,     beam_energy);
  df = define_DISCAT(df, "Mx2_eKp",      &DISANAMath::GetMx2_eKp,      beam_energy);
  df = define_DISCAT(df, "Mx_eKp",       &DISANAMath::GetMx_eKp,       beam_energy);

  df = define_DISCAT(df, "DeltaPhi",                 &DISANAMath::GetDeltaPhi,                 beam_energy);
  df = define_DISCAT(df, "Theta_g_phimeson",         &DISANAMath::GetTheta_g_phimeson,         beam_energy);
  df = define_DISCAT(df, "Theta_e_phimeson",         &DISANAMath::GetTheta_e_phimeson,         beam_energy);
  df = define_DISCAT(df, "DeltaE",                   &DISANAMath::GetDeltaE,                   beam_energy);
  df = define_DISCAT(df, "Cone_p",                   &DISANAMath::GetCone_p,                   beam_energy);
  df = define_DISCAT(df, "Cone_Kp",                  &DISANAMath::GetCone_Kp,                  beam_energy);
  df = define_DISCAT(df, "Cone_Km",                  &DISANAMath::GetCone_Km,                  beam_energy);
  df = define_DISCAT(df, "z_phi",                    &DISANAMath::GetZ_phi,                    beam_energy);
  df = define_DISCAT(df, "Coplanarity_had_normals_deg",
                     &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  df = df.Define("Gamma_v", [beam_energy](double Q2, double W) {
                 return HandGammaV(beam_energy, Q2, W);
               }, {"Q2", "W"});
  // --- new: vm_cut and Egamma_star -----------------------------------------
  df = df.Define("vm_cut",
    [beam_energy](double Q2, double xB, double t) -> double {
      return DISANAMath::ComputeVmApprox(Q2, xB, t, m_phi, beam_energy);
    }, {"Q2", "xB", "t"})
  .Define("Egamma_star",
    [](double Mx2) -> double {
      if (Mx2 <= 0.0) return -999.0;
      const double Mx = std::sqrt(Mx2);
      return (Mx2 - kMp * kMp) / (2.0 * kMp);
    }, {"Mx2_eKpKm"});

  return df;
}

ROOT::RDF::RNode WriteSlimAndReload_exclusive(ROOT::RDF::RNode df, const std::string& outFile, const std::string& outTree) {
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep = {
      // Original single-particle projections
      "REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass", "ele_px_org", "ele_py_org", "ele_pz_org", "recel_vz_org", "reckMinus_vz", "kMinus_px", "kMinus_py",
      "kMinus_pz", "reckPlus_vz", "kPlus_px", "kPlus_py", "kPlus_pz", "pro_px", "pro_py", "pro_pz", "recpro_vz", "REC_Event_helicity",

      // Run/event and counting
      "RunNumber", "EventNumber", "nElectrons", "bestEle_idx", "REC_pass_bestE", "nElectrons_best",

      // Best-e scalar copies
      "ele_px", "ele_py", "ele_pz", "recel_vz", "ele_det_region",

      // Basic kinematics (magnitudes/angles)
      "recel_p", "recel_theta", "recel_phi", "reckMinus_p", "reckMinus_theta", "reckMinus_phi", "reckPlus_p", "reckPlus_theta", "reckPlus_phi", "recpro_p", "recpro_theta",
      "recpro_phi",

      // Det region tags
      "kMinus_det_region", "kPlus_det_region", "pro_det_region", "ele_det_region_org",

      // Simple composites
      "invMass_KpKm", "invMass_pKminus", "invMass_pKplus",

      // DISANAMath-derived
      "Q2", "xB", "t", "cos_thetaKK", "cos_phiKK", "tmin", "mtprime", "tprime", "phi", "W", "nu", "y", "Gamma_v", "z_phi", "Mx2_ep", "Emiss", "PTmiss", "Mx2_epKpKm", "Mx2_eKpKm", "Mx2_eKp",
      "Mx_eKp", "Mx2_epKm", "Mx2_epKp", "DeltaPhi", "Theta_g_phimeson", "Theta_e_phimeson", "DeltaE", "Cone_p", "Cone_Kp", "Cone_Km", "Coplanarity_had_normals_deg",
                                         "vm_cut", "Egamma_star"};

  // Write the slim tree (this triggers the event loop)
  df.Snapshot(outTree, outFile, keep);

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim;  // implicitly converts to RNode
}

ROOT::RDF::RNode WriteSlimAndReload_missingKm(ROOT::RDF::RNode df, const std::string& outFile, const std::string& outTree) {
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep = {// Original single-particle projections
                                         "REC_Particle_pid", "REC_Particle_pass", "ele_px_org", "ele_py_org", "ele_pz_org", "recel_vz_org", "reckMinus_vz", "reckPlus_vz",
                                         "kPlus_px", "kPlus_py", "kPlus_pz", "pro_px", "pro_py", "pro_pz", "recpro_vz", "REC_Event_helicity",

                                         // Run/event and counting
                                         "RunNumber", "EventNumber", "nElectrons", "bestEle_idx", "REC_pass_bestE", "nElectrons_best",

                                         // Best-e scalar copies
                                         "ele_px", "ele_py", "ele_pz", "recel_vz", "ele_det_region",

                                         // Basic kinematics (magnitudes/angles)
                                         "recel_p", "recel_theta", "recel_phi", "reckMinus_p", "reckMinus_theta", "reckMinus_phi", "reckPlus_p", "reckPlus_theta", "reckPlus_phi",
                                         "recpro_p", "recpro_theta", "recpro_phi", "kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz",

                                         // Det region tags
                                         "kMinus_det_region", "kPlus_det_region", "pro_det_region", "ele_det_region_org",

                                         // Simple composites
                                         "invMass_KpKm", "invMass_pKminus", "invMass_pKplus",

                                         // DISANAMath-derived
                                         "Q2", "xB", "t", "tmin", "cos_thetaKK", "cos_phiKK", "mtprime", "tprime", "phi", "W", "nu", "y", "Gamma_v", "z_phi", "Mx2_ep", "Emiss", "PTmiss",
                                         "Mx2_epKpKm", "Mx2_eKpKm", "Mx2_eKp", "Mx_eKp", "Mx2_epKm", "Mx2_epKp", "DeltaPhi", "Theta_g_phimeson", "Theta_e_phimeson", "DeltaE",
                                         "Cone_p", "Cone_Kp", "Cone_Km", "Coplanarity_had_normals_deg",
                                         "vm_cut", "Egamma_star"};

  // Write the slim tree (this triggers the event loop)
  df.Snapshot(outTree, outFile, keep);

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim;  // implicitly converts to RNode
}
ROOT::RDF::RNode WriteSlimAndReload_missingKp(ROOT::RDF::RNode df, const std::string& outFile, const std::string& outTree) {
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep = {// Original single-particle projections
                                         "REC_Particle_pid", "REC_Particle_pass", "ele_px_org", "ele_py_org", "ele_pz_org", "recel_vz_org", "reckMinus_vz", "kMinus_px",
                                         "kMinus_py", "kMinus_pz", "reckPlus_vz", "pro_px", "pro_py", "pro_pz", "recpro_vz", "REC_Event_helicity",

                                         // Run/event and counting
                                         "RunNumber", "EventNumber", "nElectrons", "bestEle_idx", "REC_pass_bestE", "nElectrons_best",

                                         // Best-e scalar copies
                                         "ele_px", "ele_py", "ele_pz", "recel_vz", "ele_det_region",

                                         // Basic kinematics (magnitudes/angles)
                                         "recel_p", "recel_theta", "recel_phi", "reckMinus_p", "reckMinus_theta", "reckMinus_phi", "reckPlus_p", "reckPlus_theta", "reckPlus_phi",
                                         "recpro_p", "recpro_theta", "recpro_phi", "kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz",

                                         // Det region tags
                                         "kMinus_det_region", "kPlus_det_region", "pro_det_region", "ele_det_region_org",

                                         // Simple composites
                                         "invMass_KpKm", "invMass_pKminus", "invMass_pKplus",

                                         // DISANAMath-derived
                                         "Q2", "xB", "t", "tmin", "cos_thetaKK", "cos_phiKK", "mtprime", "tprime", "phi", "W", "nu", "y", "Gamma_v","z_phi", "Mx2_ep", "Emiss", "PTmiss",
                                         "Mx2_epKpKm", "Mx2_eKpKm", "Mx2_eKp", "Mx_eKp", "Mx2_epKm", "Mx2_epKp", "DeltaPhi", "Theta_g_phimeson", "Theta_e_phimeson", "DeltaE",
                                         "Cone_p", "Cone_Kp", "Cone_Km", "Coplanarity_had_normals_deg",
                                         "vm_cut", "Egamma_star"};
  // Write the slim tree (this triggers the event loop)
  df.Snapshot(outTree, outFile, keep);

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim;  // implicitly converts to RNode
}
ROOT::RDF::RNode GetSlim_missingKm(ROOT::RDF::RNode src, const std::string& f, const std::string& t) {
  const bool fileExists = !gSystem->AccessPathName(f.c_str());  // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
    std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_missingKm(src, f, t);
  }
}
ROOT::RDF::RNode GetSlim_missingKp(ROOT::RDF::RNode src, const std::string& f, const std::string& t) {
  const bool fileExists = !gSystem->AccessPathName(f.c_str());  // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
    std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_missingKp(src, f, t);
  }
}
ROOT::RDF::RNode GetSlim_exclusive(ROOT::RDF::RNode src, const std::string& f, const std::string& t) {
  const bool fileExists = !gSystem->AccessPathName(f.c_str());  // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
    std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_exclusive(src, f, t);
  }
}
//  Dump one line per event with 26 columns:
//  1)  RunNumber
//  2)  EventNumber
//  3)  helicity state
//  4–6)  electron   P, theta, phi (Lab)
//  7–9)  proton     P, theta, phi (Lab)
// 10–12) K+        P, theta, phi (Lab)
// 13–15) K-        P, theta, phi (Lab)
// 16–18) x_B, Q^2, W
// 19)  phi* (Trento phi of the Phi)
// 20)  cos(theta)  (K+ angle in Phi rest frame)
// 21)  cos(varphi) (angle of decay plane)
// 22)  cone angle between calculated and measured proton
// 23)  missing mass of e' p X
// 24)  missing mass of e' K+ K- X
// 25)  t
// 26)  t'
ROOT::RDF::RNode DumpExclusiveTxt(ROOT::RDF::RNode df, const std::string& outTxt) {
  std::ofstream out(outTxt);
  if (!out.is_open()) {
    throw std::runtime_error("DumpExclusiveTxt: cannot open file " + outTxt);
  }

  // Header documenting all 26 columns
  out << "# 1:RunNumber 2:EventNumber 3:Helicity "
      << "4:ele_p 5:ele_theta 6:ele_phi "
      << "7:pro_p 8:pro_theta 9:pro_phi "
      << "10:Kp_p 11:Kp_theta 12:Kp_phi "
      << "13:Km_p 14:Km_theta 15:Km_phi "
      << "16:xB 17:Q2 18:W "
      << "19:phi_star 20:cos_theta 21:cos_varphi "
      << "22:cone_p 23:MM_ep 24:MM_eKpKm "
      << "25:t 26:tprime"
      << "MM_eKp\n";

  df.Foreach(
      [&](int run,      // RunNumber (int)
          int evt,      // EventNumber (int)
          Short_t hel,  // REC_Event_helicity (short)
          double ele_p, double ele_theta, double ele_phi, double pro_p, double pro_theta, double pro_phi, double kp_p, double kp_theta, double kp_phi, double km_p, double km_theta,
          double km_phi, double xB, double Q2, double W, double phi_trento, double cos_theta_Kp, double cos_phi_decay, double cone_p, double Mx2_ep, double Mx2_eKpKm, double t,
          double tprime, double Mx2_eKp) {
        // 23–24: missing masses from squared values
        const double MM_ep = std::sqrt(std::max(0.0, Mx2_ep));
        const double MM_eKpKm = std::sqrt(std::max(0.0, Mx2_eKpKm));
        const double MM_eKp = std::sqrt(std::max(0.0, Mx2_eKp));

        out << run << ' '            // 1
            << evt << ' '            // 2
            << hel << ' '            // 3
            << ele_p << ' '          // 4
            << ele_theta << ' '      // 5
            << ele_phi << ' '        // 6
            << pro_p << ' '          // 7
            << pro_theta << ' '      // 8
            << pro_phi << ' '        // 9
            << kp_p << ' '           // 10
            << kp_theta << ' '       // 11
            << kp_phi << ' '         // 12
            << km_p << ' '           // 13
            << km_theta << ' '       // 14
            << km_phi << ' '         // 15
            << xB << ' '             // 16
            << Q2 << ' '             // 17
            << W << ' '              // 18
            << phi_trento << ' '     // 19
            << cos_theta_Kp << ' '   // 20
            << cos_phi_decay << ' '  // 21
            << cone_p << ' '         // 22
            << MM_ep << ' '          // 23
            << MM_eKpKm << ' '       // 24
            << t << ' '              // 25
            << tprime << ' '         // 26
            << MM_eKp << ' '         // 27
            << '\n';
      },
      {"RunNumber",           // 1
       "EventNumber",         // 2
       "REC_Event_helicity",  // 3
       "recel_p",
       "recel_theta",
       "recel_phi",  // 4–6
       "recpro_p",
       "recpro_theta",
       "recpro_phi",  // 7–9
       "reckPlus_p",
       "reckPlus_theta",
       "reckPlus_phi",  // 10–12
       "reckMinus_p",
       "reckMinus_theta",
       "reckMinus_phi",  // 13–15
       "xB",
       "Q2",
       "W",            // 16–18
       "phi",          // 19: phi* (Trento)
       "cos_thetaKK",  // 20: cos(theta)
       "cos_phiKK",    // 21: cos(varphi)
       "Cone_p",       // 22: cone angle (p)
       "Mx2_ep",       // 23: MM^2(e' p X)
       "Mx2_eKpKm",    // 24: MM^2(e' K+K- X)
       "t",            // 25
       "tprime",       // 26
       "Mx2_eKp"});

  return df;
}
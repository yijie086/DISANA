#include <cmath>
#include <memory>

#include "DISANAMath.h"  // <-- REQUIRED: provides DISANAMath class & methods
#include "DISANAMathFitUtils.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLorentzVector.h"

// --- constants for masses
static constexpr double kMe = 0.000511;
static constexpr double kMp = 0.938272;
static constexpr double kMK = 0.493677;

// convenience 3-vector helpers
static double MomentumFunc(float px, float py, float pz) { return std::sqrt(px * px + py * py + pz * pz); }
static double ThetaFunc(float px, float py, float pz) { return std::acos(pz / std::sqrt(px * px + py * py + pz * pz)); }
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}

inline int ChooseBestElectron(const ROOT::VecOps::RVec<int>& pid,
                              const ROOT::VecOps::RVec<float>& px,
                              const ROOT::VecOps::RVec<float>& py,
                              const ROOT::VecOps::RVec<float>& pz,
                              const ROOT::VecOps::RVec<float>& vz,
                              const ROOT::VecOps::RVec<bool>& pass) {
  int best = -1; float bestP = -1.f;
  // Prefer electrons in the target window; choose highest |p|
  for (size_t i=0;i<pid.size();++i) {
    if (pid[i]!=11 || !pass[i]) continue;
    if (vz[i] < -10.f || vz[i] > 3.f) continue;
    float P = std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
    if (P>bestP){ bestP=P; best=(int)i; }
  }
  if (best>=0) return best;
  // Fallback: highest |p| among passing tracks
  for (size_t i=0;i<pid.size();++i) {
    if (pid[i]!=11 || !pass[i]) continue;
    float P = std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
    if (P>bestP){ bestP=P; best=(int)i; }
  }
  return best; // -1 if no usable e⁻
}
inline int DetRegionFromStatus(short s){
  int a = std::abs(s);
  if (a>=1000 && a<2000) return 0; // FT
  if (a>=2000 && a<3000) return 1; // FD
  if (a>=4000 && a<5000) return 2; // CD
  return -1;
}
// phi event selection cuts for single photon contaminations
ROOT::RDF::RNode SelectExclusivePhiEvent(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& daughterPass) {
        int e = 0, km = 0, kp = 0, p = 0;
        bool hasPhiDaughter = true;

        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;

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
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass) {
        int e = 0, kp = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
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
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass) {
        int e = 0, km = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
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
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& daughterPass) {
        int e = 0, g = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recel_vz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

                     ////
              .Define("reckMinus_vz", // dummy for K minus vz plotting
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("reckPlus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kPlus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kPlus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kPlus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})

             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})

             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
              .Define("RunNumber",   "RUN_config_run")
              .Define("EventNumber", "RUN_config_event")
              .Define("nElectrons",
                   [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass){
                     int n=0; for (size_t i=0;i<pid.size();++i) if (pid[i]==11 && pass[i]) ++n; return n;
                   },
                   {"REC_Particle_pid","REC_Particle_pass"})
            // Choose ONE best electron per event
             .Define("bestEle_idx", ChooseBestElectron,
                    {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz",
                    "REC_Particle_vz","REC_Particle_pass"})

            // If you want to **keep** multi-e events but only use the best e⁻, leave the next line.
            // If you want to **drop** multi-e events (strict DVEP), add: .Filter("nElectrons==1","Cut: exactly 1 e⁻")
            .Filter("bestEle_idx >= 0", "At least one usable e⁻")
              .Define("REC_pass_bestE",
    [](int best, const ROOT::VecOps::RVec<int>& pid,
       const ROOT::VecOps::RVec<bool>& pass){
      ROOT::VecOps::RVec<bool> keep(pass.size(), false);
      if (best >= 0) keep[best] = true;        // only best e⁻
      for (size_t i=0;i<pid.size();++i)        // keep non-e tracks as they were
        if (pid[i]!=11) keep[i] = pass[i];
      return keep;
    }, {"bestEle_idx","REC_Particle_pid","REC_Particle_pass"})
            
    .Define("nElectrons_best", [](const ROOT::VecOps::RVec<int>& pid,
       const ROOT::VecOps::RVec<bool>& keep){    int n=0; for (size_t i=0;i<pid.size();++i) if (pid[i]==11 && keep[i]) ++n; return n;
    }, {"REC_Particle_pid","REC_pass_bestE"})
     // Project the chosen e⁻ to the scalar columns used everywhere else
            .Define("ele_px", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_px"})
            .Define("ele_py", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_py"})
            .Define("ele_pz", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_pz"})
            .Define("recel_vz",[](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_vz"})
            .Define("ele_det_region",[](int i,const ROOT::VecOps::RVec<short>& st){ return i>=0 ? DetRegionFromStatus(st[i]) : -1; },
                    {"bestEle_idx","REC_Particle_status"});

    // missing K⁻ 4-vector (components); keep both px/py/pz and derived p,θ,φ
  constexpr double kMe = 0.000511, kMp = 0.938272, kMK = 0.493677;
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

            .Filter([](float ex, float kPlusx, float px) { return ex != -999 && kPlusx != -999 && px != -999; },
                     {"ele_px", "kPlus_px", "pro_px"})
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
             .Define("kMinus_det_region", // dummy for K minus det region plotting
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     {"kPlus_px", "kPlus_py", "kPlus_pz", "kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz"});
  // φ mass built from missing K+ and measured K-
  // DISANAMath-driven observables
*df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "tmin", &DISANAMath::GetTmin, beam_energy);
  *df_=df_->Define("mtprime",        // non-negative: this is what you’ll plot
                 [](double mt, double tmin){ return std::abs(-mt - std::abs(tmin)); },
                 {"t", "tmin"})
      .Define("tprime",         // optional signed t' ≤ 0
                 [](double mtp){ return -mtp; },
                 {"mtprime"});
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);

  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_p", &DISANAMath::GetCone_p, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Kp", &DISANAMath::GetCone_Kp, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Km", &DISANAMath::GetCone_Km, beam_energy);
  *df_ = define_DISCAT(*df_, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  return *df_;
}

// -----------------------------------------------------------------------------
// InitKinematics_MissingKp : K⁺ omitted (exclusive K⁻ channel)
// -----------------------------------------------------------------------------
ROOT::RDF::RNode InitKinematics_MissingKp(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);

   *df_ = df_->Define("ele_px_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recel_vz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

                     ////
              .Define("reckMinus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kMinus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kMinus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kMinus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("reckPlus_vz", // dummy for K plus vz plotting
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

              .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})

             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
              .Define("RunNumber",   "RUN_config_run")
              .Define("EventNumber", "RUN_config_event")
              .Define("nElectrons",
                   [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass){
                     int n=0; for (size_t i=0;i<pid.size();++i) if (pid[i]==11 && pass[i]) ++n; return n;
                   },
                   {"REC_Particle_pid","REC_Particle_pass"})
            // Choose ONE best electron per event
             .Define("bestEle_idx", ChooseBestElectron,
                    {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz",
                    "REC_Particle_vz","REC_Particle_pass"})

            // If you want to **keep** multi-e events but only use the best e⁻, leave the next line.
            // If you want to **drop** multi-e events (strict DVEP), add: .Filter("nElectrons==1","Cut: exactly 1 e⁻")
            .Filter("bestEle_idx >= 0", "At least one usable e⁻")
              .Define("REC_pass_bestE",
    [](int best, const ROOT::VecOps::RVec<int>& pid,
       const ROOT::VecOps::RVec<bool>& pass){
      ROOT::VecOps::RVec<bool> keep(pass.size(), false);
      if (best >= 0) keep[best] = true;        // only best e⁻
      for (size_t i=0;i<pid.size();++i)        // keep non-e tracks as they were
        if (pid[i]!=11) keep[i] = pass[i];
      return keep;
    }, {"bestEle_idx","REC_Particle_pid","REC_Particle_pass"})
            
    .Define("nElectrons_best", [](const ROOT::VecOps::RVec<int>& pid,
       const ROOT::VecOps::RVec<bool>& keep){    int n=0; for (size_t i=0;i<pid.size();++i) if (pid[i]==11 && keep[i]) ++n; return n;
    }, {"REC_Particle_pid","REC_pass_bestE"})
     // Project the chosen e⁻ to the scalar columns used everywhere else
            .Define("ele_px", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_px"})
            .Define("ele_py", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_py"})
            .Define("ele_pz", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_pz"})
            .Define("recel_vz",[](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_vz"})
            .Define("ele_det_region",[](int i,const ROOT::VecOps::RVec<short>& st){ return i>=0 ? DetRegionFromStatus(st[i]) : -1; },
                    {"bestEle_idx","REC_Particle_status"});

          
  // missing K⁺ 4-vector (components)
  constexpr double kMe = 0.000511, kMp = 0.938272, kMK = 0.493677;
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

            .Filter([](float ex, float kMinusx, float px) { return ex != -999 && kMinusx != -999 && px != -999; },
                     {"ele_px", "kMinus_px", "pro_px"})
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     {"kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz", "kMinus_px", "kMinus_py", "kMinus_pz"});
  // φ mass built from missing K+ and measured K-
  // DISANAMath-driven observables
*df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "tmin", &DISANAMath::GetTmin, beam_energy);
  *df_=df_->Define("mtprime",        // non-negative: this is what you’ll plot
                 [](double mt, double tmin){ return std::abs(-mt - std::abs(tmin)); },
                 {"t", "tmin"})
      .Define("tprime",         // optional signed t' ≤ 0
                 [](double mtp){ return -mtp; },
                 {"mtprime"});
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);

  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_p", &DISANAMath::GetCone_p, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Kp", &DISANAMath::GetCone_Kp, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Km", &DISANAMath::GetCone_Km, beam_energy);
  *df_ = define_DISCAT(*df_, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  return *df_;
}

// Thin aliases per your request (semantic names for “exclusive” channels)
inline ROOT::RDF::RNode InitKinematics_ExclusiveKp(const std::string& f, const std::string& t, float E) {
  return InitKinematics_MissingKm(f, t, E);  // exclusive K⁺ == K⁻ omitted
}
inline ROOT::RDF::RNode InitKinematics_ExclusiveKm(const std::string& f, const std::string& t, float E) {
  return InitKinematics_MissingKp(f, t, E);  // exclusive K⁻ == K⁺ omitted
}
ROOT::RDF::RNode InitKinematics(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  *df_ = df_->Define("ele_px_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recel_vz_org",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("reckMinus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kMinus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kMinus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kMinus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("reckPlus_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})

             .Define("kPlus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kPlus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kPlus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})

             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
              .Define("RunNumber",   "RUN_config_run")
              .Define("EventNumber", "RUN_config_event")
              .Define("nElectrons",
                   [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass){
                     int n=0; for (size_t i=0;i<pid.size();++i) if (pid[i]==11 && pass[i]) ++n; return n;
                   },
                   {"REC_Particle_pid","REC_Particle_pass"})
            // Choose ONE best electron per event
             .Define("bestEle_idx", ChooseBestElectron,
                    {"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz",
                    "REC_Particle_vz","REC_Particle_pass"})

            // If you want to **keep** multi-e events but only use the best e⁻, leave the next line.
            // If you want to **drop** multi-e events (strict DVEP), add: .Filter("nElectrons==1","Cut: exactly 1 e⁻")
            .Filter("bestEle_idx >= 0", "At least one usable e⁻")
              .Define("REC_pass_bestE",
    [](int best, const ROOT::VecOps::RVec<int>& pid,
       const ROOT::VecOps::RVec<bool>& pass){
      ROOT::VecOps::RVec<bool> keep(pass.size(), false);
      if (best >= 0) keep[best] = true;        // only best e⁻
      for (size_t i=0;i<pid.size();++i)        // keep non-e tracks as they were
        if (pid[i]!=11) keep[i] = pass[i];
      return keep;
    }, {"bestEle_idx","REC_Particle_pid","REC_Particle_pass"})
            
    .Define("nElectrons_best", [](const ROOT::VecOps::RVec<int>& pid,
       const ROOT::VecOps::RVec<bool>& keep){    int n=0; for (size_t i=0;i<pid.size();++i) if (pid[i]==11 && keep[i]) ++n; return n;
    }, {"REC_Particle_pid","REC_pass_bestE"})
            // Project the chosen e⁻ to the scalar columns used everywhere else
            .Define("ele_px", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_px"})
            .Define("ele_py", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_py"})
            .Define("ele_pz", [](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_pz"})
            .Define("recel_vz",[](int i,const ROOT::VecOps::RVec<float>& v){ return i>=0 ? v[i] : -999.f; }, {"bestEle_idx","REC_Particle_vz"})
            .Define("ele_det_region",[](int i,const ROOT::VecOps::RVec<short>& st){ return i>=0 ? DetRegionFromStatus(st[i]) : -1; },
                    {"bestEle_idx","REC_Particle_status"})

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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) -> int {
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
                     {"kPlus_px", "kPlus_py", "kPlus_pz", "kMinus_px", "kMinus_py", "kMinus_pz"});


  // DISANAMath-driven observables
  *df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "tmin", &DISANAMath::GetTmin, beam_energy);

  *df_=df_->Define("mtprime",        // non-negative: this is what you’ll plot
                 [](double mt, double tmin){ return std::abs(-mt - std::abs(tmin)); },
                 {"t", "tmin"})
      .Define("tprime",         // optional signed t' ≤ 0
                 [](double mtp){ return -mtp; },
                 {"mtprime"});

  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);

  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_p", &DISANAMath::GetCone_p, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Kp", &DISANAMath::GetCone_Kp, beam_energy);
  *df_ = define_DISCAT(*df_, "Cone_Km", &DISANAMath::GetCone_Km, beam_energy);
  *df_ = define_DISCAT(*df_, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, beam_energy);
  
  return *df_;
}

ROOT::RDF::RNode WriteSlimAndReload_exclusive(ROOT::RDF::RNode df,
                                    const std::string& outFile,
                                    const std::string& outTree)
{
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep = {
    // Original single-particle projections
    "REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass",
    "ele_px_org","ele_py_org","ele_pz_org","recel_vz_org",
    "reckMinus_vz","kMinus_px","kMinus_py","kMinus_pz",
    "reckPlus_vz","kPlus_px","kPlus_py","kPlus_pz",
    "pro_px","pro_py","pro_pz","recpro_vz",

    // Run/event and counting
    "RunNumber","EventNumber","nElectrons","bestEle_idx",
    "REC_pass_bestE","nElectrons_best",

    // Best-e scalar copies
    "ele_px","ele_py","ele_pz","recel_vz","ele_det_region",

    // Basic kinematics (magnitudes/angles)
    "recel_p","recel_theta","recel_phi",
    "reckMinus_p","reckMinus_theta","reckMinus_phi",
    "reckPlus_p","reckPlus_theta","reckPlus_phi",
    "recpro_p","recpro_theta","recpro_phi",

    // Det region tags
    "kMinus_det_region","kPlus_det_region","pro_det_region","ele_det_region_org",

    // Simple composites
    "invMass_KpKm",

    // DISANAMath-derived
    "Q2","xB","t","tmin","mtprime","tprime","phi","W","nu","y",
    "Mx2_ep","Emiss","PTmiss","Mx2_epKpKm","Mx2_eKpKm",
    "Mx2_epKm","Mx2_epKp","DeltaPhi","Theta_g_phimeson",
    "Theta_e_phimeson","DeltaE","Cone_p","Cone_Kp","Cone_Km",
    "Coplanarity_had_normals_deg"
  };

  // Write the slim tree (this triggers the event loop)
  df.Snapshot(outTree, outFile, keep);

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim; // implicitly converts to RNode
}

ROOT::RDF::RNode WriteSlimAndReload_missingKm(ROOT::RDF::RNode df,
                                    const std::string& outFile,
                                    const std::string& outTree)
{
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep = {
    // Original single-particle projections
    "REC_Particle_pid", "REC_Particle_pass",
    "ele_px_org","ele_py_org","ele_pz_org","recel_vz_org",
    "reckMinus_vz","reckPlus_vz","kPlus_px","kPlus_py","kPlus_pz",
    "pro_px","pro_py","pro_pz","recpro_vz",

    // Run/event and counting
    "RunNumber","EventNumber","nElectrons","bestEle_idx",
    "REC_pass_bestE","nElectrons_best",

    // Best-e scalar copies
    "ele_px","ele_py","ele_pz","recel_vz","ele_det_region",

    // Basic kinematics (magnitudes/angles)
    "recel_p","recel_theta","recel_phi",
    "reckMinus_p","reckMinus_theta","reckMinus_phi",
    "reckPlus_p","reckPlus_theta","reckPlus_phi",
    "recpro_p","recpro_theta","recpro_phi","kMinus_miss_px", "kMinus_miss_py", "kMinus_miss_pz",

    // Det region tags
    "kMinus_det_region","kPlus_det_region","pro_det_region","ele_det_region_org",

    // Simple composites
    "invMass_KpKm",

    // DISANAMath-derived
    "Q2","xB","t","tmin","mtprime", "tprime","phi","W","nu","y",
    "Mx2_ep","Emiss","PTmiss","Mx2_epKpKm","Mx2_eKpKm",
    "Mx2_epKm","Mx2_epKp","DeltaPhi","Theta_g_phimeson",
    "Theta_e_phimeson","DeltaE","Cone_p","Cone_Kp","Cone_Km",
    "Coplanarity_had_normals_deg"
  };

  // Write the slim tree (this triggers the event loop)
  df.Snapshot(outTree, outFile, keep);

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim; // implicitly converts to RNode
}
ROOT::RDF::RNode WriteSlimAndReload_missingKp(ROOT::RDF::RNode df,
                                    const std::string& outFile,
                                    const std::string& outTree)
{
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep = {
    // Original single-particle projections
    "REC_Particle_pid", "REC_Particle_pass",
    "ele_px_org","ele_py_org","ele_pz_org","recel_vz_org",
    "reckMinus_vz","kMinus_px","kMinus_py","kMinus_pz",
    "reckPlus_vz", "pro_px","pro_py","pro_pz","recpro_vz",

    // Run/event and counting
    "RunNumber","EventNumber","nElectrons","bestEle_idx",
    "REC_pass_bestE","nElectrons_best",

    // Best-e scalar copies
    "ele_px","ele_py","ele_pz","recel_vz","ele_det_region",

    // Basic kinematics (magnitudes/angles)
    "recel_p","recel_theta","recel_phi",
    "reckMinus_p","reckMinus_theta","reckMinus_phi",
    "reckPlus_p","reckPlus_theta","reckPlus_phi",
    "recpro_p","recpro_theta","recpro_phi","kPlus_miss_px", "kPlus_miss_py", "kPlus_miss_pz",

    // Det region tags
    "kMinus_det_region","kPlus_det_region","pro_det_region","ele_det_region_org",

    // Simple composites
    "invMass_KpKm",

    // DISANAMath-derived
    "Q2","xB","t","tmin","mtprime","tprime","phi","W","nu","y",
    "Mx2_ep","Emiss","PTmiss","Mx2_epKpKm","Mx2_eKpKm",
    "Mx2_epKm","Mx2_epKp","DeltaPhi","Theta_g_phimeson",
    "Theta_e_phimeson","DeltaE","Cone_p","Cone_Kp","Cone_Km",
    "Coplanarity_had_normals_deg"
  };
  // Write the slim tree (this triggers the event loop)
  df.Snapshot(outTree, outFile, keep);

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim; // implicitly converts to RNode
}
ROOT::RDF::RNode GetSlim_missingKm(ROOT::RDF::RNode src, const std::string& f, const std::string& t)
{
  const bool fileExists = !gSystem->AccessPathName(f.c_str()); // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
    std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_missingKm(src, f, t);
  }
}
ROOT::RDF::RNode GetSlim_missingKp(ROOT::RDF::RNode src, const std::string& f, const std::string& t)
{
  const bool fileExists = !gSystem->AccessPathName(f.c_str()); // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
     std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_missingKp(src, f, t);
  }
}
ROOT::RDF::RNode GetSlim_exclusive(ROOT::RDF::RNode src, const std::string& f, const std::string& t)
{
  const bool fileExists = !gSystem->AccessPathName(f.c_str()); // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
    std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_exclusive(src, f, t);
  }
}
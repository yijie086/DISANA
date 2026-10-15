// Coherent He3 DVCS GEN/REC comparison using the He3-specific comparer path.
// Run from DISANA/rootmacros: clas12root DISANA_XplotterHe3DVCS.cpp

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace He3DVCSANA {
using ROOT::VecOps::RVec;

constexpr int kElectronPid = 11;
constexpr int kPhotonPid = 22;
constexpr int kGenHe3Pid = 49;//1000020030;
constexpr int kRecHe3Pid = 49;
constexpr double kElectronMass = 0.00051099895;
constexpr double kNucleonMass = 0.938;//2720813;
constexpr double kHe3Mass = 2.80839160743;

// Fast ALERT response for an MC-truth He3 embedded in the reconstructed RDF.
// Values marked "proposal" come from the ALERT Tagged-DVCS proposal.  Values
// marked "toy" are explicit working assumptions and must eventually be
// replaced by matched ALERT GEN/REC response maps.
struct FakeHe3Config {
  bool enabled = false;             // Master switch: false preserves e'gamma REC.
  double charge = 2.0;              // He3 is doubly charged; curvature measures pT/Z.
  double magneticFieldT = 5.0;      // proposal: CLAS12 solenoid field [tesla].
  double scintRadiusM = 0.085;       // proposal: AHDC outer radius, scintillator just outside [m].
  double detectorHalfLengthM = 0.15;// proposal: about 300 mm active length [m].
  double phiGapCenterRad = 0.0;      // toy: set from the actual ALERT geometry when known.
  double phiGapWidthRad = 20.0*M_PI/180.0; // proposal: 340 deg azimuthal coverage.
  double rigidityTurnOnGeV = 0.060; // toy: pT/Z turn-on near 50--64 MeV/c scale.
  double rigidityWidthGeV = 0.015;  // toy: smooth threshold width, avoids a hard step.
  double pidEfficiency = 1.00;      // Default: ideal conditional He3 PID efficiency.
  double sigmaRigidityFraction = 0.10; // toy: proposal's ~10% transverse-momentum scale.
  double sigmaThetaRad = 0.005;     // toy: 5 mrad inferred from 200 um / 53 mm lever arm.
  double sigmaPhiRad = 0.005;       // toy: same angular-resolution scale [rad].
  double sigmaVzCm = 0.30;          // proposal: about 3 mm z resolution [cm].
  std::uint64_t randomSeed = 0xA13E3ULL; // Change for an independent toy realization.
};

// Public configuration object: after loading the macro, parameters can be
// changed from the ROOT prompt before calling DISANA_XplotterHe3DVCS().
FakeHe3Config gFakeHe3Config;

// Final GEN/REC selection hook.  The REC fake-He3 requirement is the current
// event definition.  Future exclusivity cuts can be enabled independently for
// GEN and REC by filling the expressions below; no plotting code must change.
struct FinalSelectionConfig {
  bool requireFakeHe3ForRec = true;
  bool applyGenExclusivity = false;
  bool applyRecExclusivity = true;
  std::string genExclusivityExpression;
  // REC windows are fixed at the measured pre-cut mean +/- 2 standard
  // deviations.  They are deliberately explicit so every boundary remains
  // visible and independently adjustable from the ROOT prompt.
  std::string recExclusivityExpression =
      "mm2_ehe3 > -0.2248886 && mm2_ehe3 < 0.2450634 && "
      "mm2_ehe3gamma > -0.03205732 && mm2_ehe3gamma < 0.03000788 && "
      "mm2_egamma > 6.915168 && mm2_egamma < 9.350772 && "
      "missing_energy > -0.2940161 && missing_energy < 0.4578559 && "
      "missing_pt > -0.0597752 && missing_pt < 0.1378676 && "
      "theta_gamma_gamma_pred > -0.657282 && theta_gamma_gamma_pred < 1.232862 && "
      "theta_e_gamma > 9.84466 && theta_e_gamma < 23.44614 && "
      "coplanarity > -7.13105 && coplanarity < 14.10843";
};

FinalSelectionConfig gFinalSelectionConfig;

struct FakeHe3Result {
  bool hasTruth = false, radiusPass = false, zPass = false, phiPass = false;
  bool pass = false;
  double acceptance = 0.0;
  float px = -999.f, py = -999.f, pz = -999.f, vz = -999.f;
};

// Stateless random numbers make the result reproducible with ROOT implicit MT.
std::uint64_t Mix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
double Uniform01(std::uint64_t seed) {
  return (Mix64(seed) >> 11) * (1.0 / 9007199254740992.0);
}
double Gaussian01(std::uint64_t seed) {
  const double u1 = std::max(Uniform01(seed), 1.0e-15);
  const double u2 = Uniform01(seed ^ 0xD1B54A32D192ED03ULL);
  return std::sqrt(-2.0*std::log(u1))*std::cos(2.0*M_PI*u2);
}
double WrapToPi(double x) {
  while (x > M_PI) x -= 2.0*M_PI;
  while (x <= -M_PI) x += 2.0*M_PI;
  return x;
}
double WrapTo2Pi(double x) {
  while (x >= 2.0*M_PI) x -= 2.0*M_PI;
  while (x < 0.0) x += 2.0*M_PI;
  return x;
}

FakeHe3Result MakeFakeHe3(const RVec<int>& pid,
                          const RVec<float>& px,
                          const RVec<float>& py,
                          const RVec<float>& pz,
                          const RVec<float>& vz,
                          ULong64_t entry,
                          const FakeHe3Config& cfg) {
  FakeHe3Result out;
  const size_t n = std::min({pid.size(),px.size(),py.size(),pz.size()});
  size_t index = n;
  for (size_t i=0; i<n; ++i) {
    // Some production files use the PDG nuclear code, others the GEANT3 code 49.
    if (pid[i] == 1000020030 || pid[i] == 49) { index=i; break; }
  }
  if (index == n) return out;
  out.hasTruth = true;

  const double pt = std::hypot(double(px[index]),double(py[index]));
  const double p = std::hypot(pt,double(pz[index]));
  if (!(pt > 0.0) || !(p > 0.0)) return out;
  const double theta = std::atan2(pt,double(pz[index]));
  const double phi = WrapTo2Pi(std::atan2(double(py[index]),double(px[index])));
  const double rigidity = pt/cfg.charge; // ALERT bending observable pT/Z.

  // In a uniform solenoid the transverse orbit radius is pT/(0.3 Z B).
  // A track born near the beam line reaches at most 2*orbitRadius radially.
  const double orbitRadius = pt/(0.3*cfg.charge*cfg.magneticFieldT);
  out.radiusPass = 2.0*orbitRadius >= cfg.scintRadiusM;
  if (out.radiusPass) {
    const double arg = std::clamp(cfg.scintRadiusM/(2.0*orbitRadius),0.0,1.0);
    const double transverseArc = orbitRadius*2.0*std::asin(arg);
    const double z0m = index < vz.size() ? 0.01*double(vz[index]) : 0.0;
    const double zHit = z0m + transverseArc*double(pz[index])/pt;
    out.zPass = std::abs(zHit) < cfg.detectorHalfLengthM;
  }
  out.phiPass = std::abs(WrapToPi(phi-cfg.phiGapCenterRad)) >
                0.5*cfg.phiGapWidthRad;

  const double turnOn = 0.5*(1.0+std::erf(
      (rigidity-cfg.rigidityTurnOnGeV)/
      (std::sqrt(2.0)*cfg.rigidityWidthGeV)));
  out.acceptance = out.radiusPass && out.zPass && out.phiPass
                     ? std::clamp(turnOn*cfg.pidEfficiency,0.0,1.0) : 0.0;

  const std::uint64_t base = Mix64(cfg.randomSeed ^ std::uint64_t(entry) ^
                                   (std::uint64_t(index)<<32));
  out.pass = Uniform01(base ^ 0x01ULL) < out.acceptance;
  if (!out.pass) return out;

  // Smear the measured quantities, not Cartesian components: ALERT measures
  // pT/Z from curvature and angles from the stereo-wire helix fit.
  const double rigidityRec = std::max(1.0e-6, rigidity*(1.0+
      cfg.sigmaRigidityFraction*Gaussian01(base ^ 0x11ULL)));
  const double ptRec = cfg.charge*rigidityRec;
  const double thetaRec = std::clamp(theta+
      cfg.sigmaThetaRad*Gaussian01(base ^ 0x22ULL),1.0e-4,M_PI-1.0e-4);
  const double phiRec = WrapTo2Pi(phi+
      cfg.sigmaPhiRad*Gaussian01(base ^ 0x33ULL));
  out.px = float(ptRec*std::cos(phiRec));
  out.py = float(ptRec*std::sin(phiRec));
  out.pz = float(ptRec/std::tan(thetaRec));
  const double truthVz = index < vz.size() ? vz[index] : 0.0;
  out.vz = float(truthVz+cfg.sigmaVzCm*Gaussian01(base ^ 0x44ULL));
  return out;
}

void RequireColumns(ROOT::RDF::RNode df,
                    const std::vector<std::string>& required) {
  std::string missing;
  for (const auto& name : required) {
    if (!df.HasColumn(name)) missing += (missing.empty() ? "" : ", ") + name;
  }
  if (!missing.empty())
    throw std::runtime_error("Missing required branch(es): " + missing);
}

// Convert an entry number plus a fixed seed into a reproducible uniform number
// in [0,1). This stateless hash is safe with RDataFrame implicit MT, and every
// graph evaluation makes the same keep/reject decision for a given tree entry.
std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

double EntryHashToUnit(ULong64_t entry, std::uint64_t seed) {
  const std::uint64_t value = SplitMix64(
      static_cast<std::uint64_t>(entry) ^ seed);
  return static_cast<double>(value >> 11) *
         (1.0 / 9007199254740992.0); // 2^53.
}

ROOT::RDF::RNode ApplyIndependentSampling(ROOT::RDF::RNode df,
                                         double fraction,
                                         std::uint64_t seed,
                                         const std::string& label) {
  if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
    throw std::runtime_error("Sampling fraction must be finite and in [0,1]");
  if (fraction >= 1.0) return df;

  return df.Filter(
      [fraction, seed](ULong64_t entry) {
        return EntryHashToUnit(entry, seed) < fraction;
      },
      {"rdfentry_"}, label);
}

// Event-level pseudo-experiment for finite longitudinal polarizations.  A
// true spin is retained with probability (1+P)/2 and flipped with probability
// (1-P)/2.  Stateless entry hashes make the reassignment reproducible and MT
// safe. GEN and REC receive different base seeds because no common event ID is
// persisted in the two trees.
ROOT::RDF::RNode ApplyPolarizationResampling(ROOT::RDF::RNode df,
                                            double beamPolarization,
                                            double targetPolarization,
                                            std::uint64_t seed) {
  if (!std::isfinite(beamPolarization) || !std::isfinite(targetPolarization) ||
      beamPolarization <= 0.0 || beamPolarization > 1.0 ||
      targetPolarization <= 0.0 || targetPolarization > 1.0)
    throw std::runtime_error("He3 beam/target polarizations must be in (0,1]");

  return df.Define("real_beam_spin",
      [beamPolarization,seed](int spin,ULong64_t entry) {
        if (spin==0) return 0;
        const bool keep=EntryHashToUnit(entry,seed^0xb34d5eedULL)<
                        0.5*(1.0+beamPolarization);
        return keep?spin:-spin;
      },{"beam_spin","rdfentry_"})
    .Define("real_target_spin",
      [targetPolarization,seed](int spin,ULong64_t entry) {
        if (spin==0) return 0;
        const bool keep=EntryHashToUnit(entry,seed^0x7a296e7ULL)<
                        0.5*(1.0+targetPolarization);
        return keep?spin:-spin;
      },{"target_spin","rdfentry_"});
}

ROOT::RDF::RNode ApplyFinalSelection(ROOT::RDF::RNode df, bool isRec,
                                    const FinalSelectionConfig& cfg) {
  if (isRec && cfg.requireFakeHe3ForRec)
    df = df.Filter("fake_he3_pass", "REC final state: e' He3' gamma");

  const bool applyExclusivity = isRec ? cfg.applyRecExclusivity
                                      : cfg.applyGenExclusivity;
  const std::string& expression = isRec ? cfg.recExclusivityExpression
                                        : cfg.genExclusivityExpression;
  if (applyExclusivity) {
    if (expression.empty())
      throw std::runtime_error(std::string(isRec ? "REC" : "GEN") +
                               " exclusivity cut is enabled but empty");
    std::cout << "[He3DVCS final selection] Applying "
              << (isRec ? "REC" : "GEN")
              << " exclusivity cut: " << expression << '\n';
    df = df.Filter(expression, isRec ? "REC exclusivity selection"
                                     : "GEN exclusivity selection");
  }
  return df;
}

bool BranchIsBool(const std::string& file, const std::string& tree,
                  const std::string& branch) {
  TFile input(file.c_str(), "READ");
  auto* t = dynamic_cast<TTree*>(input.Get(tree.c_str()));
  if (!t) throw std::runtime_error("Tree not found: " + tree);
  auto* leaf = t->GetLeaf(branch.c_str());
  if (!leaf) throw std::runtime_error("Branch not found: " + branch);
  const std::string type = leaf->GetTypeName();
  return type == "Bool_t" || type == "bool" || type == "vector<bool>";
}

ROOT::RDF::RNode NormalizePassColumns(const std::string& file,
                                      const std::string& tree) {
  ROOT::RDataFrame input(tree, file);
  ROOT::RDF::RNode df(input);
  if (BranchIsBool(file, tree, "REC_Particle_pass")) {
    df = df.Define("REC_Particle_pass_std",
      [](const RVec<bool>& v) { return v; }, {"REC_Particle_pass"});
  } else {
    df = df.Define("REC_Particle_pass_std", [](const RVec<int>& v) {
      RVec<bool> out(v.size());
      for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] != 0;
      return out;
    }, {"REC_Particle_pass"});
  }
  if (BranchIsBool(file, tree, "REC_Photon_MaxE")) {
    df = df.Define("REC_Photon_MaxE_std",
      [](const RVec<bool>& v) { return v; }, {"REC_Photon_MaxE"});
  } else {
    df = df.Define("REC_Photon_MaxE_std", [](const RVec<int>& v) {
      RVec<bool> out(v.size());
      for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] != 0;
      return out;
    }, {"REC_Photon_MaxE"});
  }
  return df;
}

float FindGen(const RVec<int>& pid, const RVec<float>& value, int wanted) {
  for (size_t i = 0, n = std::min(pid.size(), value.size()); i < n; ++i)
    if (pid[i] == wanted) return value[i];
  return -999.0f;
}

float FindRec(const RVec<int>& pid, const RVec<float>& value,
              const RVec<bool>& pass, const RVec<bool>& maxe,
              int wanted, bool require_maxe) {
  const size_t n = std::min({pid.size(), value.size(), pass.size()});
  for (size_t i = 0; i < n; ++i) {
    if (pid[i] == wanted && pass[i] &&
        (!require_maxe || (i < maxe.size() && maxe[i]))) return value[i];
  }
  return -999.0f;
}

int DetectorRegion(const RVec<int>& pid, const RVec<short>& status,
                   const RVec<bool>& pass, const RVec<bool>& maxe,
                   int wanted, bool require_maxe) {
  const size_t n = std::min({pid.size(), status.size(), pass.size()});
  for (size_t i = 0; i < n; ++i) {
    if (pid[i] != wanted || !pass[i] ||
        (require_maxe && (i >= maxe.size() || !maxe[i]))) continue;
    const int s = std::abs(static_cast<int>(status[i]));
    if (s >= 1000 && s < 2000) return 0;
    if (s >= 2000 && s < 3000) return 1;
    if (s >= 4000 && s < 5000) return 2;
    return -1;
  }
  return -1;
}

double Momentum(float x, float y, float z) {
  return std::sqrt(double(x)*x + double(y)*y + double(z)*z);
}
double ThetaDeg(float x, float y, float z) {
  const double p = Momentum(x, y, z);
  return p > 0 ? std::acos(std::clamp(double(z)/p, -1.0, 1.0))*180.0/M_PI : -999.0;
}
double PhiDeg(float x, float y) {
  double phi = std::atan2(y, x)*180.0/M_PI;
  return phi < 0 ? phi + 360.0 : phi;
}
TLorentzVector P4(float x, float y, float z, double mass) {
  TLorentzVector p; p.SetXYZM(x, y, z, mass); return p;
}
TLorentzVector GammaP4(float x, float y, float z) {
  TLorentzVector p; p.SetPxPyPzE(x, y, z, Momentum(x,y,z)); return p;
}
double TrentoPhi(const TLorentzVector& beam, const TLorentzVector& electron,
                 const TLorentzVector& photon) {
  const TVector3 q = (beam-electron).Vect();
  const TVector3 nl = beam.Vect().Cross(electron.Vect());
  const TVector3 nh = q.Cross(photon.Vect());
  if (!q.Mag2() || !nl.Mag2() || !nh.Mag2()) return -999.0;
  double phi = std::atan2(q.Unit().Dot(nl.Unit().Cross(nh.Unit())),
                          nl.Unit().Dot(nh.Unit()))*180.0/M_PI;
  return phi < 0 ? phi + 360.0 : phi;
}

ROOT::RDF::RNode DefineKinematics(ROOT::RDF::RNode df, double ebeam,
                                  bool require_he3) {
  df = require_he3
      ? df.Filter("ele_px != -999.f && pho_px != -999.f && he3_px != -999.f",
                  "GEN final state: e' He3' gamma")
      : df.Filter("ele_px != -999.f && pho_px != -999.f",
                  "REC final state: e' gamma");
  return df
    .Define("recel_p",Momentum,{"ele_px","ele_py","ele_pz"})
    .Define("recel_theta",ThetaDeg,{"ele_px","ele_py","ele_pz"})
    .Define("recel_phi",PhiDeg,{"ele_px","ele_py"})
    .Define("recpho_p",Momentum,{"pho_px","pho_py","pho_pz"})
    .Define("recpho_theta",ThetaDeg,{"pho_px","pho_py","pho_pz"})
    .Define("recpho_phi",PhiDeg,{"pho_px","pho_py"})
    .Filter("recpho_p >= 2.0", "photon energy >= 2 GeV")
    .Define("reche3_p",Momentum,{"he3_px","he3_py","he3_pz"})
    .Define("reche3_theta",ThetaDeg,{"he3_px","he3_py","he3_pz"})
    .Define("reche3_phi",PhiDeg,{"he3_px","he3_py"})
    .Define("beam_p4",[ebeam](){return TLorentzVector(0,0,ebeam,ebeam);})
    .Define("target_p4",[](){return TLorentzVector(0,0,0,kHe3Mass);})
    .Define("electron_p4",[](float x,float y,float z){return P4(x,y,z,kElectronMass);},{"ele_px","ele_py","ele_pz"})
    .Define("photon_p4",GammaP4,{"pho_px","pho_py","pho_pz"})
    .Define("he3_p4",[](float x,float y,float z){return P4(x,y,z,kHe3Mass);},{"he3_px","he3_py","he3_pz"})
    .Define("q_p4",[](const TLorentzVector& b,const TLorentzVector& e){return b-e;},{"beam_p4","electron_p4"})
    .Define("missing_p4",[](const TLorentzVector& b,const TLorentzVector& t,const TLorentzVector& e,const TLorentzVector& g,const TLorentzVector& h){return b+t-e-g-h;},{"beam_p4","target_p4","electron_p4","photon_p4","he3_p4"})
    .Define("missing_ehe3_p4",[](const TLorentzVector& b,const TLorentzVector& t,const TLorentzVector& e,const TLorentzVector& h){return b+t-e-h;},{"beam_p4","target_p4","electron_p4","he3_p4"})
    .Define("missing_egamma_p4",[](const TLorentzVector& b,const TLorentzVector& t,const TLorentzVector& e,const TLorentzVector& g){return b+t-e-g;},{"beam_p4","target_p4","electron_p4","photon_p4"})
    .Define("Q2",[](const TLorentzVector& q){return -q.M2();},{"q_p4"})
    .Define("nu",[](const TLorentzVector& q){return q.E();},{"q_p4"})
    .Define("xB",[](double q2,double nu){return nu>0?q2/(2*kNucleonMass*nu):-999.;},{"Q2","nu"})
    .Define("WN",[](const TLorentzVector& q){double w2=(TLorentzVector(0,0,0,kNucleonMass)+q).M2();return w2>=0?std::sqrt(w2):-999.;},{"q_p4"})
    .Define("WHe3",[](const TLorentzVector& q){double w2=(TLorentzVector(0,0,0,kHe3Mass)+q).M2();return w2>=0?std::sqrt(w2):-999.;},{"q_p4"})
    // Backward-compatible name used by existing cuts and plotting code.
    .Define("W",[](double wn){return wn;},{"WN"})
    .Define("y",[ebeam](double nu){return nu/ebeam;},{"nu"})
    .Define("t",[](const TLorentzVector& q,const TLorentzVector& g){return (q-g).M2();},{"q_p4","photon_p4"})
    .Define("minus_t",[](double t){return -t;},{"t"})
    .Define("phi",TrentoPhi,{"beam_p4","electron_p4","photon_p4"})
    .Define("missing_energy",[](const TLorentzVector& p){return p.E();},{"missing_p4"})
    .Define("missing_pt",[](const TLorentzVector& p){return p.Pt();},{"missing_p4"})
    .Define("mm2_ehe3",[](const TLorentzVector& p){return p.M2();},{"missing_ehe3_p4"})
    .Define("mm2_ehe3gamma",[](const TLorentzVector& p){return p.M2();},{"missing_p4"})
    .Define("mm2_egamma",[](const TLorentzVector& p){return p.M2();},{"missing_egamma_p4"})
    .Define("theta_gamma_gamma_pred",[](const TLorentzVector& g,const TLorentzVector& pred){return g.P()>0&&pred.P()>0?g.Vect().Angle(pred.Vect())*180/M_PI:-999.;},{"photon_p4","missing_ehe3_p4"})
    .Define("theta_e_gamma",[](const TLorentzVector& e,const TLorentzVector& g){return e.Vect().Angle(g.Vect())*180/M_PI;},{"electron_p4","photon_p4"})
    .Define("coplanarity",[](const TLorentzVector& q,const TLorentzVector& g,const TLorentzVector& h){double d=q.Vect().Cross(g.Vect()).Angle(q.Vect().Cross(-h.Vect()))*180/M_PI;return d;},{"q_p4","photon_p4","he3_p4"})
    .Define("event_weight",[](){return 1.0;});
}
} // namespace He3DVCSANA

ROOT::RDF::RNode InitGenKinematicsHe3DVCS(const std::string& file,
    const std::string& tree="dfSelectedMC", double ebeam=10.6,
    const He3DVCSANA::FakeHe3Config& fake={}) {
  using namespace He3DVCSANA;
  ROOT::RDataFrame input(tree,file); ROOT::RDF::RNode df(input);
  RequireColumns(df,{"MC_Particle_pid","MC_Particle_px","MC_Particle_py","MC_Particle_pz","MC_Particle_vz",
                     "MC_Event_pbeam","MC_Event_ptarget"});
  auto pick=[](int wanted){return [wanted](const RVec<int>& p,const RVec<float>& v){return FindGen(p,v,wanted);};};
  df=df.Define("ele_px",pick(kElectronPid),{"MC_Particle_pid","MC_Particle_px"}).Define("ele_py",pick(kElectronPid),{"MC_Particle_pid","MC_Particle_py"}).Define("ele_pz",pick(kElectronPid),{"MC_Particle_pid","MC_Particle_pz"})
       .Define("pho_px",pick(kPhotonPid),{"MC_Particle_pid","MC_Particle_px"}).Define("pho_py",pick(kPhotonPid),{"MC_Particle_pid","MC_Particle_py"}).Define("pho_pz",pick(kPhotonPid),{"MC_Particle_pid","MC_Particle_pz"})
       .Define("he3_px",pick(kGenHe3Pid),{"MC_Particle_pid","MC_Particle_px"}).Define("he3_py",pick(kGenHe3Pid),{"MC_Particle_pid","MC_Particle_py"}).Define("he3_pz",pick(kGenHe3Pid),{"MC_Particle_pid","MC_Particle_pz"})
       // pbeam/ptarget are floating-point polarization values in these trees.
       // JIT expressions let ROOT use the exact on-disk scalar type.
       .Define("beam_spin","MC_Event_pbeam > 0 ? 1 : (MC_Event_pbeam < 0 ? -1 : 0)")
       .Define("target_spin","MC_Event_ptarget > 0 ? 1 : (MC_Event_ptarget < 0 ? -1 : 0)")
       .Define("recel_vz",[](){return 0.;}).Define("recpho_vz",[](){return 0.;}).Define("reche3_vz",[](){return 0.;})
       .Define("ele_det_region",[](){return 1;}).Define("pho_det_region",[](){return 1;}).Define("he3_det_region",[](){return 1;})
       .Define("fake_he3_gen_result",[fake](const RVec<int>& pid,
          const RVec<float>& px,const RVec<float>& py,const RVec<float>& pz,
          const RVec<float>& vz,ULong64_t entry) {
            return MakeFakeHe3(pid,px,py,pz,vz,entry,fake);
          },{"MC_Particle_pid","MC_Particle_px","MC_Particle_py",
             "MC_Particle_pz","MC_Particle_vz","rdfentry_"})
       .Define("fake_he3_acceptance",[](const FakeHe3Result& r){return r.acceptance;},
               {"fake_he3_gen_result"});
  return DefineKinematics(df,ebeam,true);
}

ROOT::RDF::RNode InitKinematicsHe3DVCS(const std::string& file,
    const std::string& tree="dfSelected", double ebeam=10.6,
    const He3DVCSANA::FakeHe3Config& fake={}) {
  using namespace He3DVCSANA;
  auto df=NormalizePassColumns(file,tree);
  RequireColumns(df,{"REC_Particle_pid","REC_Particle_px","REC_Particle_py","REC_Particle_pz","REC_Particle_vz","REC_Particle_status",
                     "MC_Event_pbeam","MC_Event_ptarget"});
  if (fake.enabled) RequireColumns(df,{"MC_Particle_pid","MC_Particle_px",
      "MC_Particle_py","MC_Particle_pz","MC_Particle_vz"});
  auto pick=[](int wanted,bool maxe){return [wanted,maxe](const RVec<int>& p,const RVec<float>& v,const RVec<bool>& pass,const RVec<bool>& m){return FindRec(p,v,pass,m,wanted,maxe);};};
  auto reg=[](int wanted,bool maxe){return [wanted,maxe](const RVec<int>& p,const RVec<short>& s,const RVec<bool>& pass,const RVec<bool>& m){return DetectorRegion(p,s,pass,m,wanted,maxe);};};
  auto cols=[&](const char* v){return std::vector<std::string>{"REC_Particle_pid",v,"REC_Particle_pass_std","REC_Photon_MaxE_std"};};
  df=df.Define("ele_px",pick(kElectronPid,false),cols("REC_Particle_px")).Define("ele_py",pick(kElectronPid,false),cols("REC_Particle_py")).Define("ele_pz",pick(kElectronPid,false),cols("REC_Particle_pz"))
       .Define("pho_px",pick(kPhotonPid,true),cols("REC_Particle_px")).Define("pho_py",pick(kPhotonPid,true),cols("REC_Particle_py")).Define("pho_pz",pick(kPhotonPid,true),cols("REC_Particle_pz"))
       .Define("recel_vz",pick(kElectronPid,false),cols("REC_Particle_vz")).Define("recpho_vz",pick(kPhotonPid,true),cols("REC_Particle_vz"))
       .Define("ele_det_region",reg(kElectronPid,false),cols("REC_Particle_status")).Define("pho_det_region",reg(kPhotonPid,true),cols("REC_Particle_status"))
       .Define("beam_spin","MC_Event_pbeam > 0 ? 1 : (MC_Event_pbeam < 0 ? -1 : 0)")
       .Define("target_spin","MC_Event_ptarget > 0 ? 1 : (MC_Event_ptarget < 0 ? -1 : 0)");
  // Standard CLAS12 status convention: 1=FD, 0=FT, 2=CD.
  // This is a REC-only detector acceptance requirement; GEN has no REC status.
  df=df.Filter("ele_det_region == 1","electron in Forward Detector");

  if (fake.enabled) {
    // One correlated result object guarantees that acceptance and all smeared
    // components use the same MC He3 and the same event-level random draw.
    df=df.Define("fake_he3_result",[fake](const RVec<int>& pid,
        const RVec<float>& px,const RVec<float>& py,const RVec<float>& pz,
        const RVec<float>& vz,ULong64_t entry) {
          return MakeFakeHe3(pid,px,py,pz,vz,entry,fake);
        },{"MC_Particle_pid","MC_Particle_px","MC_Particle_py",
           "MC_Particle_pz","MC_Particle_vz","rdfentry_"})
       .Define("fake_he3_has_truth",[](const FakeHe3Result& r){return r.hasTruth;},{"fake_he3_result"})
       .Define("fake_he3_radius_pass",[](const FakeHe3Result& r){return r.radiusPass;},{"fake_he3_result"})
       .Define("fake_he3_z_pass",[](const FakeHe3Result& r){return r.zPass;},{"fake_he3_result"})
       .Define("fake_he3_phi_pass",[](const FakeHe3Result& r){return r.phiPass;},{"fake_he3_result"})
       .Define("fake_he3_acceptance",[](const FakeHe3Result& r){return r.acceptance;},{"fake_he3_result"})
       .Define("fake_he3_pass",[](const FakeHe3Result& r){return r.pass;},{"fake_he3_result"})
       .Define("he3_px",[](const FakeHe3Result& r){return r.px;},{"fake_he3_result"})
       .Define("he3_py",[](const FakeHe3Result& r){return r.py;},{"fake_he3_result"})
       .Define("he3_pz",[](const FakeHe3Result& r){return r.pz;},{"fake_he3_result"})
       .Define("reche3_vz",[](const FakeHe3Result& r){return double(r.vz);},{"fake_he3_result"})
       .Define("he3_det_region",[](){return 3;}); // 3 = fake ALERT, not FT/FD/CD.
  } else {
    df=df.Define("he3_px",pick(kRecHe3Pid,false),cols("REC_Particle_px"))
       .Define("he3_py",pick(kRecHe3Pid,false),cols("REC_Particle_py"))
       .Define("he3_pz",pick(kRecHe3Pid,false),cols("REC_Particle_pz"))
       .Define("reche3_vz",pick(kRecHe3Pid,false),cols("REC_Particle_vz"))
       .Define("he3_det_region",reg(kRecHe3Pid,false),cols("REC_Particle_status"))
       .Define("fake_he3_pass",[](){return false;})
       .Define("fake_he3_acceptance",[](){return 0.0;});
  }
  return DefineKinematics(df,ebeam,false);
}

void DISANA_XplotterHe3DVCS(const std::string& directory="../build/He3DVCSrate",
                            bool enable_fake_he3=true,
                            double sampling_fraction=0.05276*15,
                            double beam_polarization=0.80,
                            double target_polarization=0.60,
                            bool draw_corrected_asymmetry=true) {
  ROOT::EnableImplicitMT(40); constexpr double ebeam=10.6;
  He3DVCSANA::FakeHe3Config fake = He3DVCSANA::gFakeHe3Config;
  fake.enabled = enable_fake_he3; // All ALERT fast-simulation parameters are above.
  auto gen=InitGenKinematicsHe3DVCS(directory+"/dfSelectedMC.root","dfSelectedMC",ebeam,fake);
  auto rec=InitKinematicsHe3DVCS(directory+"/dfSelected_afterFid.root","dfSelected_afterFid",ebeam,fake);
  // No common event identifier is persisted in the two trees. Therefore GEN
  // and REC are thinned independently. Distinct fixed seeds keep the samples
  // independent, deterministic, reproducible, and safe with implicit MT.
  gen=He3DVCSANA::ApplyIndependentSampling(
      gen,sampling_fraction,0x47554e5f484533ULL,
      "GEN independent luminosity sampling");
  rec=He3DVCSANA::ApplyIndependentSampling(
      rec,sampling_fraction,0x5245435f484533ULL,
      "REC independent luminosity sampling");
  gen=He3DVCSANA::ApplyFinalSelection(
      gen,false,He3DVCSANA::gFinalSelectionConfig);
  rec=He3DVCSANA::ApplyFinalSelection(
      rec,true,He3DVCSANA::gFinalSelectionConfig);
  gen=He3DVCSANA::ApplyPolarizationResampling(
      gen,beam_polarization,target_polarization,0x47454e5f504f4cULL);
  rec=He3DVCSANA::ApplyPolarizationResampling(
      rec,beam_polarization,target_polarization,0x5245435f504f4cULL);
  DISANAcomparer comparer;
  comparer.SetOutputDir("He3DVCS_gen_rec_plots");
  comparer.AddModelHe3DVCS(gen);  // model 0: GEN
  comparer.AddModelHe3DVCS(rec);  // model 1: REC
  comparer.PlotNormalizedComparisonHe3DVCS(enable_fake_he3);
  comparer.PlotHe3DVCSAcc();
  comparer.PlotAsymmetryHe3(beam_polarization,target_polarization,
                            draw_corrected_asymmetry);
  auto ngen=gen.Count(), nrec=rec.Count();
  auto nfake=rec.Filter("fake_he3_pass","fake ALERT He3 accepted").Count();
  std::cout<<"He3DVCS, Ebeam="<<ebeam<<" GeV\nGEN: "<<ngen.GetValue()
           <<"\nREC e'gamma: "<<nrec.GetValue()<<"\nFake He3: "
           <<(enable_fake_he3?"ON":"OFF")
           <<"\nREC fake-He3 accepted: "<<nfake.GetValue()<<"\nWeights: 1\n";
  std::cout<<"Independent luminosity sampling: "
           <<(sampling_fraction<1.0?"ON":"OFF")
           <<"\nSampling fraction: "<<sampling_fraction
           <<"\nGEN and REC use independent deterministic seeds.\n";
  std::cout<<"Beam polarization Pb: "<<beam_polarization
           <<"\nTarget polarization Pt: "<<target_polarization
           <<"\nPolarization simulation: event-level deterministic resampling"
           <<"\nAsymmetry output: "
           <<(draw_corrected_asymmetry?"reconstructed A":"raw A")<<'\n';
}

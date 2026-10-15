// He3 DIS GEN/REC comparison for epp, epn, and ed final states.
// Run from DISANA/rootmacros:
//   clas12root DISANA_XplotterHe3DIS.cpp
// The default fake spectator response is He3Detector with one 30-degree
// module.  The final arguments select detector, module count/width/offset;
// choosing "ALERT" reproduces the retained legacy fake response.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <TVector3.h>

#include <dlfcn.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../../Simulation/He3DIS/He3DIS/src/spin_model_api.h"
#include "../../../../Simulation/He3DIS/He3DIS/src/nucleus/polarization.h"

namespace He3DISPlot {
using ROOT::VecOps::RVec;

constexpr int kElectronPid = 11;
constexpr int kProtonPid = 2212;
constexpr int kNeutronPid = 2112;
constexpr int kDeuteronRecPid = 45;
constexpr int kDeuteronPdgPid = 1000010020;
constexpr double kNucleonMass = 0.93892;
constexpr double kHe3Mass = 2.80839160743;

// Shared x_active binning for the He3DIS count projection, all three g1
// extraction methods, and the polarized-EMC ratio plots.  Edit only here.
constexpr int kHe3DISXActiveBins = 6;
constexpr double kHe3DISXActiveMin = 0.0;
constexpr double kHe3DISXActiveMax = 0.96;
constexpr double kHe3DISModelXMax = 0.85;

// Double-spectator epp selection from the tagged-neutron study. In the
// stationary-He3 frame this is also the magnitude of the struck neutron's
// initial momentum: p_n = -(p_s1+p_s2).
constexpr double kEppSpectatorPSumMaxGeV = 0.20;
// For ed two-body breakup, the spectator system is the deuteron itself, so
// |p_d| equals the magnitude of the struck proton's initial momentum in the
// stationary-He3 frame.
constexpr double kEdSpectatorMomentumMaxGeV = 0.20;

bool IsRequestedPid(int pid, int requested);

// Fast ALERT response applied to MC-truth charged spectators embedded in the
// reconstructed RDF.  Geometry is shared, while each species has independent
// efficiency, resolution, turn-on, charge, and random seed.  pidEfficiency is
// conditional: total efficiency = geometry * rigidity turn-on * PID efficiency.
struct FakeChargedConfig {
  bool enabled = true;
  double charge = 1.0;
  double magneticFieldT = 5.0;
  double scintRadiusM = 0.085;
  double detectorHalfLengthM = 0.15;
  double phiGapCenterRad = 0.0;
  double phiGapWidthRad = 20.0*M_PI/180.0;
  double rigidityTurnOnGeV = 0.060;
  double rigidityWidthGeV = 0.015;
  double pidEfficiency = 0.50;      // Default: ideal conditional p/d PID efficiency.
  double sigmaRigidityFraction = 0.10;
  double sigmaThetaRad = 0.005;
  double sigmaPhiRad = 0.005;
  double sigmaVzCm = 0.30;
  std::uint64_t randomSeed = 0xA13D15ULL;
};

struct FakeSpectatorConfig {
  bool enabled = true;
  FakeChargedConfig proton;
  FakeChargedConfig deuteron;
  bool neutronEnabled = false; // Current assumption: ALERT cannot detect n.
  FakeSpectatorConfig() {
    proton.charge = 1.0;
    proton.randomSeed = 0xA13A2212ULL;
    deuteron.charge = 1.0;
    deuteron.randomSeed = 0xA13A0045ULL;
  }
};

FakeSpectatorConfig gFakeSpectatorConfig;

// Select which charged-spectator response supplies the reconstructed p/d
// columns.  ALERT remains available unchanged; He3Detector is the new
// three-layer silicon telescope toy model; GEMC uses the stored REC particles.
enum class FakeDetectorMode { ALERT, He3Detector, GEMC };

// Local NIST proton CSDA tables. Tables are loaded once before RDF execution;
// all event threads share immutable data. No network access occurs at runtime.
struct PStarRangeTable {
  std::vector<double> energyMeV, rangeGcm2;
  std::string source;

  static PStarRangeTable Load(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open PSTAR table: "+path.string());
    PStarRangeTable table;
    table.source=path.string();
    std::string line;
    while (std::getline(input,line)) {
      const auto first=line.find_first_not_of(" \t\r");
      if (first==std::string::npos || line[first]=='#') continue;
      std::istringstream row(line);
      double t,se,sn,st,rc,rp,detour;
      if (!(row>>t>>se>>sn>>st>>rc>>rp>>detour) ||
          !std::isfinite(t) || !std::isfinite(rc) || !(t>0.0) || !(rc>0.0))
        throw std::runtime_error("Invalid PSTAR numeric row in "+path.string());
      if (!table.energyMeV.empty() &&
          (t<=table.energyMeV.back() || rc<=table.rangeGcm2.back()))
        throw std::runtime_error("PSTAR energy/range must increase in "+path.string());
      table.energyMeV.push_back(t);
      table.rangeGcm2.push_back(rc); // Fifth column: CSDA, not projected range.
    }
    if (table.energyMeV.size()<2)
      throw std::runtime_error("PSTAR table has fewer than two rows: "+path.string());
    return table;
  }

  double Interpolate(double value, const std::vector<double>& x,
                     const std::vector<double>& y) const {
    if (!std::isfinite(value) || value<0.0)
      throw std::invalid_argument("Invalid PSTAR query in "+source);
    if (value==0.0) return 0.0;
    // Explicit approximate tail below the first tabulated point; log(0) is
    // never evaluated. The same tail is used for the inverse query.
    if (value<=x.front()) return y.front()*value/x.front();
    if (value>x.back())
      throw std::out_of_range("PSTAR query exceeds table domain in "+source);
    const auto upper=std::lower_bound(x.begin(),x.end(),value);
    const size_t hi=static_cast<size_t>(upper-x.begin());
    if (*upper==value) return y[hi];
    const size_t lo=hi-1;
    const double f=std::log(value/x[lo])/std::log(x[hi]/x[lo]);
    return y[lo]*std::exp(f*std::log(y[hi]/y[lo]));
  }
  double Range(double tMeV) const {
    return Interpolate(tMeV,energyMeV,rangeGcm2);
  }
  double Energy(double rGcm2) const {
    return Interpolate(rGcm2,rangeGcm2,energyMeV);
  }
};

struct He3PStarTables {
  PStarRangeTable aluminum, silicon;
};

const He3PStarTables& GetHe3PStarTables() {
  static const He3PStarTables tables=[] {
    const auto directory=(std::filesystem::absolute(__FILE__).parent_path()/
                          "data/pstar").lexically_normal();
    return He3PStarTables{
        PStarRangeTable::Load(directory/"aluminum.dat"),
        PStarRangeTable::Load(directory/"silicon.dat")};
  }();
  return tables;
}

constexpr double kProtonMassGeV=0.9382720813;
constexpr double kDeuteronMassGeV=1.87561294257;
constexpr double kHelixGeVPerTeslaM=0.299792458;

double He3ChargedMassGeV(int requested) {
  if (requested==kProtonPid) return kProtonMassGeV;
  if (requested==kDeuteronRecPid || requested==kDeuteronPdgPid)
    return kDeuteronMassGeV;
  throw std::invalid_argument("PSTAR He3 response supports only proton/deuteron");
}

struct He3TrackState {
  TVector3 positionM;
  TVector3 momentumGeV;
  double massGeV=kProtonMassGeV;
  double charge=1.0;
  bool stopped=false;
};

double He3KineticEnergyMeV(const He3TrackState& state) {
  const double p2=state.momentumGeV.Mag2();
  // Stable at low momentum: sqrt(p^2+m^2)-m = p^2/(sqrt(p^2+m^2)+m).
  return 1000.0*p2/(std::sqrt(p2+state.massGeV*state.massGeV)+state.massGeV);
}

struct FakeHe3DetectorLayerConfig {
  double radiusM = 0.0;
  double thicknessM = 0.0; // Instantaneous PSTAR loss at the nominal radius.
  double efficiency = 0.50;
  int strips = 0;
};

struct FakeHe3DetectorSpeciesConfig {
  double charge = 1.0;
  // Overall conditional efficiency, applied once after transport/geometry.
  double pidEfficiency = 0.50;
  double sigmaMomentumFraction = 0.10;
  double sigmaThetaRad = 0.005;
  double sigmaPhiRad = 0.005;
  double sigmaVzCm = 0.30;
  std::uint64_t randomSeed = 0xA13D15ULL;
};

struct FakeHe3DetectorConfig {
  bool enabled = true;
  double magneticFieldT = 5.0;
  int fieldPolarity = +1;

  // Thin-face target model: local Al window or non-window Al wall.
  // Gas, dead layers, scattering and straggling are not modeled.
  double targetRadiusM = 0.010;
  double targetHalfLengthM = 0.100;
  double windowThicknessM = 2.5e-6;
  double cellWallThicknessM = 1.0e-4;

  // Identical, uniformly spaced azimuthal telescope modules.  moduleOffsetRad
  // is the center of module zero.  All layers share this layout in v1.
  int moduleCount = 1;
  double moduleWidthRad = 30.0*M_PI/180.0;
  double moduleOffsetRad = 0.0;
  // Minimum required silicon arrival: 1=L1, 2=L2, 3=L3.
  // The particle may stop in that layer; later layers are optional.
  int requiredLayer = 2;
  bool requireSameModule = true;
  bool windowFollowsModules = true;
  double detectorHalfLengthM = 0.050;

  std::array<FakeHe3DetectorLayerConfig,3> layers{{
      {0.0133, 60.e-6, 1.0, 0},
      {0.0200, 100.e-6, 1.0, 32},
      {0.0333, 1.0e-3, 1.0, 200}}};
  FakeHe3DetectorSpeciesConfig proton;
  FakeHe3DetectorSpeciesConfig deuteron;

  FakeHe3DetectorConfig() {
    proton.charge = 1.0;
    proton.randomSeed = 0x3E3A2212ULL;
    deuteron.charge = 1.0;
    deuteron.randomSeed = 0x3E3A0045ULL;
  }
};

FakeHe3DetectorConfig gFakeHe3DetectorConfig;

std::uint64_t Mix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
double Uniform01(std::uint64_t seed) {
  return (Mix64(seed) >> 11) * (1.0/9007199254740992.0);
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

FakeDetectorMode ParseFakeDetectorMode(const std::string& name) {
  if (name=="ALERT" || name=="Alert" || name=="alert")
    return FakeDetectorMode::ALERT;
  if (name=="He3Detector" || name=="HE3DETECTOR" ||
      name=="he3detector" || name=="He3")
    return FakeDetectorMode::He3Detector;
  if (name=="GEMC" || name=="Gemc" || name=="gemc")
    return FakeDetectorMode::GEMC;
  throw std::invalid_argument(
      "fake_detector must be ALERT, He3Detector, or GEMC (received '"+
      name+"')");
}

const char* FakeDetectorModeName(FakeDetectorMode mode) {
  if (mode==FakeDetectorMode::ALERT) return "ALERT";
  if (mode==FakeDetectorMode::He3Detector) return "He3Detector";
  return "GEMC";
}

void ValidateFakeHe3DetectorConfig(const FakeHe3DetectorConfig& cfg) {
  if (!cfg.enabled) return;
  if (cfg.requiredLayer<1 || cfg.requiredLayer>3)
    throw std::invalid_argument("Fake He3Detector requiredLayer must be 1, 2, or 3");
  if (cfg.moduleCount<1)
    throw std::invalid_argument("Fake He3Detector moduleCount must be >= 1");
  const double spacing=2.0*M_PI/cfg.moduleCount;
  if (!(cfg.moduleWidthRad>0.0) || cfg.moduleWidthRad>spacing+1.e-12)
    throw std::invalid_argument(
        "Fake He3Detector module width must be in (0, 360/moduleCount] deg");
  if (!(cfg.magneticFieldT>0.0) || cfg.fieldPolarity==0)
    throw std::invalid_argument(
        "Fake He3Detector requires nonzero field polarity and B > 0");
  if (!(cfg.targetRadiusM>0.0) || !(cfg.targetHalfLengthM>0.0) ||
      !(cfg.detectorHalfLengthM>0.0) || cfg.windowThicknessM<0.0 ||
      cfg.cellWallThicknessM<0.0)
    throw std::invalid_argument("Invalid Fake He3Detector target dimensions");
  double previous=cfg.targetRadiusM;
  for (const auto& layer:cfg.layers) {
    if (!(layer.radiusM>previous) || !(layer.thicknessM>0.0) ||
        !(layer.efficiency>=0.0 && layer.efficiency<=1.0))
      throw std::invalid_argument(
          "Fake He3Detector layer radii must increase outside the target, "
          "with positive thickness and efficiency in [0,1]");
    previous=layer.radiusM;
  }
  for (const auto* species:{&cfg.proton,&cfg.deuteron}) {
    if (!(std::abs(species->charge)>0.0) ||
        !(species->pidEfficiency>=0.0 && species->pidEfficiency<=1.0) ||
        species->sigmaMomentumFraction<0.0 || species->sigmaThetaRad<0.0 ||
        species->sigmaPhiRad<0.0 || species->sigmaVzCm<0.0)
      throw std::invalid_argument(
          "Invalid Fake He3Detector charge, PID efficiency, or resolution");
  }
}

int FindHe3DetectorModule(double phi,
                          const FakeHe3DetectorConfig& cfg) {
  for (int module=0;module<cfg.moduleCount;++module) {
    const double center=cfg.moduleOffsetRad+
                        2.0*M_PI*module/cfg.moduleCount;
    if (std::abs(WrapToPi(phi-center))<=0.5*cfg.moduleWidthRad)
      return module;
  }
  return -1;
}

struct FakeHe3DetectorHit {
  bool reachable=false;
  bool zPass=false;
  int module=-1;
  double phi=0.0; // POSITION azimuth; distinct from the momentum azimuth.
  double zM=0.0;
  double transverseArcM=0.0;
  He3TrackState state; // Full position and tangent momentum at this crossing.
};

FakeHe3DetectorHit PropagateToHe3DetectorRadius(
    const He3TrackState& start, double radiusM, double halfLengthM,
    const FakeHe3DetectorConfig& cfg, bool requireModule) {
  FakeHe3DetectorHit hit;
  hit.state=start;
  const double pt=start.momentumGeV.Perp();
  if (start.stopped || !(pt>0.0) || !(radiusM>0.0)) return hit;
  const double x0=start.positionM.X(), y0=start.positionM.Y();
  const double phi0=start.momentumGeV.Phi();
  const double bz=cfg.magneticFieldT*(cfg.fieldPolarity>=0 ? 1.0 : -1.0);
  const double qb=kHelixGeVPerTeslaM*start.charge*bz;
  double arcM=std::numeric_limits<double>::infinity();
  double xHit=0.0,yHit=0.0,phiMomentum=phi0;
  if (std::abs(qb)<1.e-14) {
    const double ux=start.momentumGeV.X()/pt, uy=start.momentumGeV.Y()/pt;
    const double b=x0*ux+y0*uy;
    const double disc=b*b+radiusM*radiusM-x0*x0-y0*y0;
    if (disc<0.0) return hit;
    for (double s:{-b-std::sqrt(disc),-b+std::sqrt(disc)})
      if (s>1.e-12 && s<arcM) arcM=s;
    if (!std::isfinite(arcM)) return hit;
    xHit=x0+arcM*ux; yHit=y0+arcM*uy;
  } else {
    // Lorentz-force convention: positive q*Bz rotates momentum clockwise.
    // Recompute the orbit center from the CURRENT crossing and momentum.
    const double cx=x0+start.momentumGeV.Y()/qb;
    const double cy=y0-start.momentumGeV.X()/qb;
    const double d=std::hypot(cx,cy);
    if (!(d>1.e-14)) return hit; // Concentric orbit cannot reach a new radius.
    const double signedCurvature=qb/pt;
    const double orbitRadius=1.0/std::abs(signedCurvature);
    // d^2-r_orbit^2 without subtracting two large squared orbit radii.
    const double offset=x0*x0+y0*y0+
        2.0*(x0*start.momentumGeV.Y()-y0*start.momentumGeV.X())/qb;
    const double a=(radiusM*radiusM+offset)/(2.0*d);
    double h2=radiusM*radiusM-a*a;
    if (h2 < -1.e-12) return hit;
    h2=std::max(0.0,h2);
    const double ux=cx/d, uy=cy/d, h=std::sqrt(h2);
    const double bendSign=qb>0.0 ? 1.0 : -1.0;
    for (double side:{-1.0,1.0}) {
      const double x=a*ux-side*h*uy, y=a*uy+side*h*ux;
      const double tangentPhi=std::atan2(-signedCurvature*(x-cx),
                                         signedCurvature*(y-cy));
      double angle=std::fmod(bendSign*(phi0-tangentPhi),2.0*M_PI);
      if (angle<0.0) angle+=2.0*M_PI;
      const double s=orbitRadius*angle;
      if (s>1.e-12 && s<arcM) {
        arcM=s; xHit=x; yHit=y; phiMomentum=tangentPhi;
      }
    }
    if (!std::isfinite(arcM)) return hit;
  }
  const double zHit=start.positionM.Z()+arcM*start.momentumGeV.Z()/pt;
  hit.reachable=true;
  hit.transverseArcM=arcM;
  hit.zM=zHit;
  hit.phi=std::atan2(yHit,xHit);
  hit.zPass=std::isfinite(zHit) && std::abs(zHit)<halfLengthM;
  hit.module=requireModule ? FindHe3DetectorModule(hit.phi,cfg) : 0;
  hit.state.positionM.SetXYZ(xHit,yHit,zHit);
  hit.state.momentumGeV.SetXYZ(pt*std::cos(phiMomentum),
                              pt*std::sin(phiMomentum),start.momentumGeV.Z());
  return hit;
}

struct FakeHe3MaterialCrossing {
  FakeHe3DetectorHit hit;
  He3TrackState outgoing;
  bool active=false; // Thin window / active silicon actually present here.
  bool thinWindow=false; // Region zero otherwise corresponds to the Al wall.
  double thicknessM=0.0;
  double incidenceCos=0.0;
  double effectivePathM=0.0; // Full thin-face path t/|u.n|.
  double pathM=0.0;          // Material path before stopping, capped by range.
  double tinMeV=std::numeric_limits<double>::quiet_NaN();
  double toutMeV=std::numeric_limits<double>::quiet_NaN();
  double depositMeV=0.0;
  double pinGeV=std::numeric_limits<double>::quiet_NaN();
  double poutGeV=std::numeric_limits<double>::quiet_NaN();
};

void ApplyInstantHe3EnergyLoss(FakeHe3MaterialCrossing& crossing,
                               const PStarRangeTable& table,
                               double densityGcm3, double thicknessM) {
  const auto& incoming=crossing.hit.state;
  crossing.outgoing=incoming;
  crossing.thicknessM=thicknessM;
  crossing.pinGeV=incoming.momentumGeV.Mag();
  crossing.poutGeV=crossing.pinGeV;
  crossing.tinMeV=He3KineticEnergyMeV(incoming);
  crossing.toutMeV=crossing.tinMeV;
  if (!(thicknessM>0.0) || !(crossing.pinGeV>0.0)) return;
  const TVector3 normal(incoming.positionM.X(),incoming.positionM.Y(),0.0);
  crossing.incidenceCos=std::clamp(std::abs(
      incoming.momentumGeV.Unit().Dot(normal.Unit())),0.0,1.0);
  crossing.effectivePathM=crossing.incidenceCos>0.0
      ? thicknessM/crossing.incidenceCos
      : std::numeric_limits<double>::infinity();
  // Equal-velocity mass scaling: exact proton table, approximate deuteron
  // total CSDA scaling (including the approximated nuclear stopping tail).
  const double massRatio=incoming.massGeV/kProtonMassGeV;
  const double rIn=massRatio*table.Range(crossing.tinMeV/massRatio);
  const double densityPerM=100.0*densityGcm3; // m -> cm for g/cm^2.
  crossing.pathM=std::min(crossing.effectivePathM,rIn/densityPerM);
  const double fullArealDensity=densityPerM*crossing.effectivePathM;
  if (fullArealDensity>=rIn) {
    crossing.toutMeV=0.0;
    crossing.poutGeV=0.0;
    crossing.outgoing.momentumGeV.SetXYZ(0.0,0.0,0.0);
    crossing.outgoing.stopped=true;
  } else {
    crossing.toutMeV=std::clamp(
        massRatio*table.Energy((rIn-fullArealDensity)/massRatio),
        0.0,crossing.tinMeV);
    const double tGeV=crossing.toutMeV/1000.0;
    crossing.poutGeV=std::sqrt(tGeV*(tGeV+2.0*incoming.massGeV));
    // Scale ALL three components: energy loss itself does not change direction.
    crossing.outgoing.momentumGeV*=crossing.poutGeV/crossing.pinGeV;
    crossing.outgoing.stopped=!(crossing.poutGeV>0.0);
  }
  crossing.depositMeV=crossing.tinMeV-crossing.toutMeV;
  // Thin-face approximation: no position/z displacement inside the material.
}

struct FakeHe3TransportRecord {
  // Region 0=target Al window/wall, 1=L1, 2=L2, 3=L3.
  std::array<FakeHe3MaterialCrossing,4> crossings;
  int stopRegion=-1; // -1 means not stopped in the modeled material.
  int unreachableRegion=-1;
  int failedGeometryRegion=-1;
  int siliconLayersReached=0;
  bool allLayersPass=false; // Strict L3 diagnostic, independent of config.
  bool acceptedGeometry=false; // Window + prefix through requiredLayer.
};

FakeHe3TransportRecord TransportHe3Charged(
    He3TrackState state, const FakeHe3DetectorConfig& cfg,
    const He3PStarTables& tables) {
  FakeHe3TransportRecord record;
  if (!std::isfinite(state.positionM.Z()) ||
      std::abs(state.positionM.Z())>=cfg.targetHalfLengthM) {
    record.failedGeometryRegion=0; // End-window transport is not modeled.
    return record;
  }
  for (size_t region=0;region<record.crossings.size();++region) {
    auto& crossing=record.crossings[region];
    const bool cell=region==0;
    const double radius=cell ? cfg.targetRadiusM : cfg.layers[region-1].radiusM;
    const double halfLength=cell ? cfg.targetHalfLengthM : cfg.detectorHalfLengthM;
    crossing.hit=PropagateToHe3DetectorRadius(
        state,radius,halfLength,cfg,cell ? cfg.windowFollowsModules : true);
    if (!crossing.hit.reachable) {
      record.unreachableRegion=static_cast<int>(region);
      break;
    }
    state=crossing.hit.state;
    crossing.outgoing=state;
    crossing.pinGeV=crossing.poutGeV=state.momentumGeV.Mag();
    crossing.tinMeV=crossing.toutMeV=He3KineticEnergyMeV(state);
    crossing.active=crossing.hit.zPass &&
        (cell ? (!cfg.windowFollowsModules || crossing.hit.module>=0)
              : crossing.hit.module>=0);
    if (!crossing.active && record.failedGeometryRegion<0)
      record.failedGeometryRegion=static_cast<int>(region);
    if (cell && !crossing.hit.zPass) break;
    if (cell) {
      crossing.thinWindow=crossing.active;
      ApplyInstantHe3EnergyLoss(crossing,tables.aluminum,2.70,
          crossing.thinWindow ? cfg.windowThicknessM : cfg.cellWallThicknessM);
    } else if (crossing.active) {
      ++record.siliconLayersReached;
      ApplyInstantHe3EnergyLoss(crossing,tables.silicon,2.329,
                               cfg.layers[region-1].thicknessM);
    }
    // Outside an active Si sector there is no modeled silicon sheet. Such
    // trajectories fail if that missing layer belongs to the required prefix.
    state=crossing.outgoing;
    if (state.stopped) {
      record.stopRegion=static_cast<int>(region);
      break;
    }
  }
  auto passesRequiredPrefix=[&](int lastLayer) {
    bool pass=true;
    for (int region=0;region<=lastLayer;++region) {
      const auto& crossing=record.crossings[region];
      pass=pass && crossing.hit.reachable && crossing.active;
    }
    if (pass && cfg.requireSameModule) {
      const int module=record.crossings[1].hit.module;
      for (int layer=2;layer<=lastLayer;++layer)
        pass=pass && record.crossings[layer].hit.module==module;
      if (cfg.windowFollowsModules)
        pass=pass && record.crossings[0].hit.module==module;
    }
    return pass;
  };
  // Continue transport for diagnostics, but only the configured prefix gates
  // detection. Arrival is sufficient, including stopping in the required layer.
  record.allLayersPass=passesRequiredPrefix(3);
  record.acceptedGeometry=passesRequiredPrefix(cfg.requiredLayer);
  return record;
}

struct FakeChargedParticles {
  RVec<int> pid;
  RVec<float> px, py, pz, vz;
  RVec<bool> pass;
  // One entry for every generated particle of this species, including those
  // outside acceptance.  These arrays are used for probability-weighted
  // acceptance maps and are independent of the random fake-pass decision.
  RVec<float> truthP, truthTheta, truthPhi, acceptance;
  // He3-only transport records, aligned with the truth arrays above.
  RVec<int> transportMCIndex;
  RVec<FakeHe3TransportRecord> transport;
};

FakeChargedParticles MakeFakeHe3DetectorCharged(
    const RVec<int>& pid, const RVec<float>& px, const RVec<float>& py,
    const RVec<float>& pz, const RVec<float>& vz, ULong64_t entry,
    int requested, const FakeHe3DetectorSpeciesConfig& species,
    const FakeHe3DetectorConfig& cfg) {
  FakeChargedParticles out;
  const size_t n=std::min({pid.size(),px.size(),py.size(),pz.size()});
  for (size_t i=0;i<n;++i) {
    if (!IsRequestedPid(pid[i],requested)) continue;
    const double pt=std::hypot(double(px[i]),double(py[i]));
    const double p=std::hypot(pt,double(pz[i]));
    if (!(pt>0.0) || !(p>0.0) || !cfg.enabled) continue;
    const double theta=std::atan2(pt,double(pz[i]));
    const double phi=std::atan2(double(py[i]),double(px[i]));
    const double z0M=i<vz.size() ? 0.01*double(vz[i]) : 0.0;

    He3TrackState state;
    state.positionM.SetXYZ(0.0,0.0,z0M);
    state.momentumGeV.SetXYZ(px[i],py[i],pz[i]);
    state.massGeV=He3ChargedMassGeV(requested);
    state.charge=species.charge;
    const auto transport=TransportHe3Charged(state,cfg,GetHe3PStarTables());
    const bool geometryPass=transport.acceptedGeometry;
    out.transportMCIndex.push_back(static_cast<int>(i));
    out.transport.push_back(transport);

    double efficiency=geometryPass ? species.pidEfficiency : 0.0;
    for (int layer=0;layer<cfg.requiredLayer;++layer)
      efficiency*=cfg.layers[layer].efficiency;
    efficiency=std::clamp(efficiency,0.0,1.0);
    double phiDeg=phi*180.0/M_PI;
    if (phiDeg<0.0) phiDeg+=360.0;
    out.truthP.push_back(float(p));
    out.truthTheta.push_back(float(theta*180.0/M_PI));
    out.truthPhi.push_back(float(phiDeg));
    out.acceptance.push_back(float(efficiency));

    const std::uint64_t base=Mix64(species.randomSeed^
        std::uint64_t(entry)^(std::uint64_t(i)<<32));
    if (Uniform01(base^0x01ULL)>=efficiency) continue;

    // Retain the vertex-level truth-smearing reconstruction assumption.
    // Local post-material momenta live in the transport diagnostics;
    // a stopped particle's zero momentum is not a vertex observable.
    const double pRec=std::max(1.e-6,p*(1.0+
        species.sigmaMomentumFraction*Gaussian01(base^0x11ULL)));
    const double thetaRec=std::clamp(theta+
        species.sigmaThetaRad*Gaussian01(base^0x22ULL),1.e-4,M_PI-1.e-4);
    const double phiRec=phi+species.sigmaPhiRad*Gaussian01(base^0x33ULL);
    const double ptRec=pRec*std::sin(thetaRec);
    out.pid.push_back(requested);
    out.px.push_back(float(ptRec*std::cos(phiRec)));
    out.py.push_back(float(ptRec*std::sin(phiRec)));
    out.pz.push_back(float(pRec*std::cos(thetaRec)));
    const double truthVz=i<vz.size() ? vz[i] : 0.0;
    out.vz.push_back(float(truthVz+
        species.sigmaVzCm*Gaussian01(base^0x44ULL)));
    out.pass.push_back(true);
  }
  return out;
}

// These columns contain one entry per transported MC particle, including
// particles rejected by geometry/stopping or by the final 50% random decision.
ROOT::RDF::RNode DefineHe3TransportColumns(
    ROOT::RDF::RNode df, const std::string& result) {
  const std::string prefix=result+"_transport";
  df=df.Define(prefix+"_mc_index",[](const FakeChargedParticles& r) {
        return r.transportMCIndex;
      },{result});
  const std::array<const char*,6> eventFields{{
      "stop_region","unreachable_region","failed_geometry_region",
      "silicon_layers_reached","all_layers_pass","required_layer_pass"}};
  for (size_t field=0;field<eventFields.size();++field)
    df=df.Define(prefix+"_"+eventFields[field],
      [field](const FakeChargedParticles& r) {
        RVec<int> values;
        for (const auto& record:r.transport) {
          int value=0;
          switch (field) {
            case 0: value=record.stopRegion; break;
            case 1: value=record.unreachableRegion; break;
            case 2: value=record.failedGeometryRegion; break;
            case 3: value=record.siliconLayersReached; break;
            case 4: value=record.allLayersPass ? 1 : 0; break;
            case 5: value=record.acceptedGeometry ? 1 : 0; break;
          }
          values.push_back(value);
        }
        return values;
      },{result});
  const std::array<const char*,4> regions{{"window","L1","L2","L3"}};
  const std::array<const char*,3> hitFields{{"reached","active","module"}};
  const std::array<const char*,13> materialFields{{
      "x_m","y_m","z_m","position_phi_deg","momentum_phi_deg",
      "path_m","pin_GeV","pout_GeV","Tin_MeV","Tout_MeV","dE_MeV",
      "incidence_cos","effective_path_m"}};
  for (size_t region=0;region<regions.size();++region) {
    const std::string surface=prefix+"_"+regions[region];
    for (size_t field=0;field<hitFields.size();++field)
      df=df.Define(surface+"_"+hitFields[field],
        [region,field](const FakeChargedParticles& r) {
          RVec<int> values;
          for (const auto& record:r.transport) {
            const auto& c=record.crossings[region];
            values.push_back(field==0 ? int(c.hit.reachable) :
                             field==1 ? int(c.active) : c.hit.module);
          }
          return values;
        },{result});
    for (size_t field=0;field<materialFields.size();++field)
      df=df.Define(surface+"_"+materialFields[field],
        [region,field](const FakeChargedParticles& r) {
          RVec<float> values;
          for (const auto& record:r.transport) {
            const auto& c=record.crossings[region];
            double value=std::numeric_limits<double>::quiet_NaN();
            if (c.hit.reachable) {
              switch (field) {
                case 0: value=c.hit.state.positionM.X(); break;
                case 1: value=c.hit.state.positionM.Y(); break;
                case 2: value=c.hit.state.positionM.Z(); break;
                case 3: value=c.hit.phi*180.0/M_PI; break;
                case 4: value=c.hit.state.momentumGeV.Phi()*180.0/M_PI; break;
                case 5: value=c.pathM; break;
                case 6: value=c.pinGeV; break;
                case 7: value=c.poutGeV; break;
                case 8: value=c.tinMeV; break;
                case 9: value=c.toutMeV; break;
                case 10: value=c.depositMeV; break;
                case 11: value=c.incidenceCos; break;
                case 12: value=c.effectivePathM; break;
              }
            }
            values.push_back(float(value));
          }
          return values;
        },{result});
  }
  return df;
}

FakeChargedParticles MakeFakeCharged(
    const RVec<int>& pid, const RVec<float>& px, const RVec<float>& py,
    const RVec<float>& pz, const RVec<float>& vz, ULong64_t entry,
    int requested, const FakeChargedConfig& cfg) {
  FakeChargedParticles out;
  const size_t n = std::min({pid.size(),px.size(),py.size(),pz.size()});
  for (size_t i=0; i<n; ++i) {
    if (!IsRequestedPid(pid[i], requested)) continue;
    const double pt = std::hypot(double(px[i]),double(py[i]));
    const double p = std::hypot(pt,double(pz[i]));
    if (!(pt > 0.0) || !(p > 0.0) || !cfg.enabled) continue;
    const double theta = std::atan2(pt,double(pz[i]));
    const double phi = std::atan2(double(py[i]),double(px[i]));
    const double rigidity = pt/std::abs(cfg.charge);
    const double orbitRadius = pt/(0.3*std::abs(cfg.charge)*cfg.magneticFieldT);
    const bool radiusPass = 2.0*orbitRadius >= cfg.scintRadiusM;
    bool zPass = false;
    if (radiusPass) {
      const double arg = std::clamp(cfg.scintRadiusM/(2.0*orbitRadius),0.0,1.0);
      const double transverseArc = 2.0*orbitRadius*std::asin(arg);
      const double z0m = i < vz.size() ? 0.01*double(vz[i]) : 0.0;
      const double zHit = z0m + transverseArc*double(pz[i])/pt;
      zPass = std::abs(zHit) < cfg.detectorHalfLengthM;
    }
    const bool phiPass = std::abs(WrapToPi(phi-cfg.phiGapCenterRad)) >
                         0.5*cfg.phiGapWidthRad;
    const double turnOn = 0.5*(1.0+std::erf(
        (rigidity-cfg.rigidityTurnOnGeV)/
        (std::sqrt(2.0)*cfg.rigidityWidthGeV)));
    const double efficiency = radiusPass && zPass && phiPass
        ? std::clamp(turnOn*cfg.pidEfficiency,0.0,1.0) : 0.0;
    double phiDeg=phi*180.0/M_PI;
    if (phiDeg<0.0) phiDeg+=360.0;
    out.truthP.push_back(float(p));
    out.truthTheta.push_back(float(theta*180.0/M_PI));
    out.truthPhi.push_back(float(phiDeg));
    out.acceptance.push_back(float(efficiency));
    const std::uint64_t base = Mix64(cfg.randomSeed ^ std::uint64_t(entry) ^
                                     (std::uint64_t(i)<<32));
    if (Uniform01(base ^ 0x01ULL) >= efficiency) continue;

    const double rigidityRec = std::max(1.0e-6, rigidity*(1.0+
        cfg.sigmaRigidityFraction*Gaussian01(base ^ 0x11ULL)));
    const double ptRec = std::abs(cfg.charge)*rigidityRec;
    const double thetaRec = std::clamp(theta+
        cfg.sigmaThetaRad*Gaussian01(base ^ 0x22ULL),1.0e-4,M_PI-1.0e-4);
    const double phiRec = phi+cfg.sigmaPhiRad*Gaussian01(base ^ 0x33ULL);
    out.pid.push_back(requested);
    out.px.push_back(float(ptRec*std::cos(phiRec)));
    out.py.push_back(float(ptRec*std::sin(phiRec)));
    out.pz.push_back(float(ptRec/std::tan(thetaRec)));
    const double truthVz = i < vz.size() ? vz[i] : 0.0;
    out.vz.push_back(float(truthVz+
        cfg.sigmaVzCm*Gaussian01(base ^ 0x44ULL)));
    out.pass.push_back(true);
  }
  return out;
}

struct ParticleValues {
  RVec<float> p, theta, phi, vz;
};

void RequireColumns(ROOT::RDF::RNode df,
                    const std::vector<std::string>& columns) {
  std::string missing;
  for (const auto& column : columns) {
    if (!df.HasColumn(column))
      missing += (missing.empty() ? "" : ", ") + column;
  }
  if (!missing.empty())
    throw std::runtime_error("Missing required branch(es): " + missing);
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

ROOT::RDF::RNode NormalizeRecPass(const std::string& file,
                                  const std::string& tree) {
  ROOT::RDataFrame input(tree, file);
  ROOT::RDF::RNode df(input);
  if (BranchIsBool(file, tree, "REC_Particle_pass")) {
    return df.Define("REC_Particle_pass_std",
                     [](const RVec<bool>& pass) { return pass; },
                     {"REC_Particle_pass"});
  }
  return df.Define("REC_Particle_pass_std",
                   [](const RVec<int>& pass) {
                     RVec<bool> result(pass.size());
                     for (size_t i = 0; i < pass.size(); ++i)
                       result[i] = pass[i] != 0;
                     return result;
                   },
                   {"REC_Particle_pass"});
}

bool IsRequestedPid(int pid, int requested) {
  if (requested == kDeuteronRecPid)
    return pid == kDeuteronRecPid || pid == kDeuteronPdgPid;
  return pid == requested;
}

ParticleValues SelectParticles(const RVec<int>& pid,
                               const RVec<float>& px,
                               const RVec<float>& py,
                               const RVec<float>& pz,
                               const RVec<float>& vz,
                               const RVec<bool>& pass,
                               int requested) {
  const size_t n = std::min(
      {pid.size(), px.size(), py.size(), pz.size(), pass.size()});
  std::vector<size_t> indices;
  for (size_t i = 0; i < n; ++i)
    if (pass[i] && IsRequestedPid(pid[i], requested)) indices.push_back(i);
  std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
    const double pa2 = double(px[a])*px[a] + double(py[a])*py[a] +
                       double(pz[a])*pz[a];
    const double pb2 = double(px[b])*px[b] + double(py[b])*py[b] +
                       double(pz[b])*pz[b];
    return pa2 > pb2;
  });

  ParticleValues result;
  result.p.reserve(indices.size());
  result.theta.reserve(indices.size());
  result.phi.reserve(indices.size());
  result.vz.reserve(indices.size());
  for (const size_t i : indices) {
    const double pt = std::hypot(double(px[i]), double(py[i]));
    const double momentum = std::hypot(pt, double(pz[i]));
    double phi = std::atan2(double(py[i]), double(px[i]))*180.0/M_PI;
    if (phi < 0.0) phi += 360.0;
    result.p.push_back(momentum);
    result.theta.push_back(
        momentum > 0.0 ? std::atan2(pt, double(pz[i]))*180.0/M_PI : -999.0);
    result.phi.push_back(phi);
    result.vz.push_back(i < vz.size() ? vz[i] : 0.0f);
  }
  return result;
}

float AtOrMissing(const RVec<float>& values, size_t index) {
  return index < values.size() ? values[index] : -999.0f;
}

float AtLargestMomentum(const RVec<float>& momentum,
                        const RVec<float>& values) {
  const size_t n=std::min(momentum.size(),values.size());
  if (n==0) return -999.0f;
  size_t best=0;
  for (size_t i=1;i<n;++i) if (momentum[i]>momentum[best]) best=i;
  return values[best];
}

int CountPid(const RVec<int>& pid, int requested) {
  return std::count_if(pid.begin(), pid.end(),
      [requested](int value) { return IsRequestedPid(value, requested); });
}

ROOT::RDF::RNode DefineParticleColumns(
    ROOT::RDF::RNode df, const std::string& prefix, int pid,
    const std::vector<std::string>& inputs) {
  df = df.Define(prefix + "_set",
      [pid](const RVec<int>& pids, const RVec<float>& px,
            const RVec<float>& py, const RVec<float>& pz,
            const RVec<float>& vz, const RVec<bool>& pass) {
        return SelectParticles(pids, px, py, pz, vz, pass, pid);
      }, inputs);
  return df
    .Define(prefix + "_p_all",
      [](const ParticleValues& values) { return values.p; },
      {prefix + "_set"})
    .Define(prefix + "_theta_all",
      [](const ParticleValues& values) { return values.theta; },
      {prefix + "_set"})
    .Define(prefix + "_phi_all",
      [](const ParticleValues& values) { return values.phi; },
      {prefix + "_set"})
    .Define(prefix + "_vz_all",
      [](const ParticleValues& values) { return values.vz; },
      {prefix + "_set"});
}

ROOT::RDF::RNode DefineScalar(ROOT::RDF::RNode df,
                              const std::string& output,
                              const std::string& vectorColumn,
                              size_t index) {
  return df.Define(output, [index](const RVec<float>& values) {
    return AtOrMissing(values, index);
  }, {vectorColumn});
}

// Only event spin labels come from MC metadata. No truth kinematics or MC
// asymmetry enters the REC estimator. JIT casts preserve the stored scalar
// type (GEMC often writes float; other converters may write short/double).
ROOT::RDF::RNode DefineProjectionSpin(ROOT::RDF::RNode df,
                                      double beamPolarization,
                                      double targetPolarization) {
  if (!(beamPolarization > 0.0 && beamPolarization <= 1.0) ||
      !(targetPolarization > 0.0 && targetPolarization <= 1.0))
    throw std::invalid_argument(
        "He3DIS projection polarizations must be in (0,1]");
  const auto names = df.GetColumnNames();
  for (const std::string column : {"MC_Event_pbeam", "MC_Event_ptarget"}) {
    if (std::find(names.begin(), names.end(), column) == names.end()) {
      std::cerr << "[He3DIS projection] Missing " << column
                << ": yields available, asymmetry disabled for this model.\n";
      return df.Define("projection_spin", [](){ return 0; });
    }
    const std::string type = df.GetColumnType(column);
    if (type.find("Vec") != std::string::npos ||
        type.find("vector") != std::string::npos)
      throw std::runtime_error("Projection requires scalar spin label: "+column);
  }
  return df.Define("projection_beam", "static_cast<double>(MC_Event_pbeam)")
    .Define("projection_target", "static_cast<double>(MC_Event_ptarget)")
    // A partially polarized ensemble is a mixture of the two pure spin
    // states.  Reassign each full-polarization MC label with probabilities
    // (1+P)/2 and (1-P)/2.  The entry-based hashes make the result exactly
    // reproducible and thread-order independent.
    .Define("projection_beam_mixed", [beamPolarization](double b,
        ULong64_t entry) {
      if (!std::isfinite(b) || std::abs(std::abs(b)-1.) > 1.e-6) return 0.0;
      const bool keep=Uniform01(0xB3A0C80ULL ^ std::uint64_t(entry)) <
                      0.5*(1.0+beamPolarization);
      return keep ? b : -b;
    }, {"projection_beam","rdfentry_"})
    .Define("projection_target_mixed", [targetPolarization](double t,
        ULong64_t entry) {
      if (!std::isfinite(t) || std::abs(std::abs(t)-1.) > 1.e-6) return 0.0;
      const bool keep=Uniform01(0x7A26E60ULL ^ std::uint64_t(entry)) <
                      0.5*(1.0+targetPolarization);
      return keep ? t : -t;
    }, {"projection_target","rdfentry_"})
    .Define("projection_spin", [](double b, double t) {
      if (!std::isfinite(b) || !std::isfinite(t) ||
          std::abs(std::abs(b)-1.) > 1.e-6 ||
          std::abs(std::abs(t)-1.) > 1.e-6) return 0;
      return b*t > 0 ? 1 : -1; // parallel / antiparallel to beam axis
    }, {"projection_beam_mixed", "projection_target_mixed"});
}

// Spectator observables in the target rest frame. The longitudinal axis is
// q = k - k', NOT the beam/magnetic-field axis used by the fake detector.
// Columns are topology-prefixed so that rec.all remains electron-inclusive.
ROOT::RDF::RNode DefineSpectatorKinematics(ROOT::RDF::RNode df,
                                         double beamEnergy) {
  for (const std::string topo : {"epp", "epn", "ed"}) {
    const bool single = topo == "ed";
    const std::string first = single ? "deu" : "pro1";
    const std::string second = topo == "epn" ? "neu" : "pro2";
    const double m1 = single ? 1.87561294257 : 0.9382720813;
    const double m2 = topo == "epn" ? 0.9395654133 : 0.9382720813;
    const std::string other = single ? first : second;
    df = df.Define(topo + "_spectator_kinematics",
      [beamEnergy, single, m1, m2](float pe, float te, float fe,
          float p1, float t1, float f1, float p2, float t2, float f2) {
        // Missing candidates must not turn into physical momenta. No event
        // filter is introduced here; invalid values go outside plot ranges.
        if (!(pe > 0) || p1 < 0 || (!single && p2 < 0))
          return RVec<double>{-999., -999., -999., -999., -999.};
        auto momentum = [](double p, double theta, double phi) {
          const double t = theta*M_PI/180., f = phi*M_PI/180.;
          return TVector3(p*std::sin(t)*std::cos(f),
                          p*std::sin(t)*std::sin(f), p*std::cos(t));
        };
        const TVector3 q = TVector3(0,0,beamEnergy)-momentum(pe,te,fe);
        if (!(q.Mag2() > 0))
          return RVec<double>{-999., -999., -999., -999., -999.};
        TVector3 ps = momentum(p1,t1,f1);
        double energy = std::sqrt(m1*m1+double(p1)*p1);
        if (!single) {
          ps += momentum(p2,t2,f2);
          energy += std::sqrt(m2*m2+double(p2)*p2);
        }
        const double longitudinal = ps.Dot(q.Unit());
        // p_sum is |p1+p2|, not |p1|+|p2|. alpha_S is the spectator
        // minus light-cone fraction: 3(E_S-p_L,S)/M_He3 (near 2 at rest).
        return RVec<double>{ps.Mag(), 3.*(energy-longitudinal)/kHe3Mass,
                            (ps-longitudinal*q.Unit()).Mag(),energy,
                            longitudinal};
      }, {"recel_p", "recel_theta", "recel_phi", first+"_p",
          first+"_theta", first+"_phi", other+"_p", other+"_theta", other+"_phi"});
    int index = 0;
    for (const std::string variable :
         {"p_sum", "alpha", "pT", "spectator_energy", "spectator_pL"}) {
      df = df.Define(topo+"_"+variable,
          [index](const RVec<double>& v) { return v[index]; },
          {topo+"_spectator_kinematics"});
      ++index;
    }
  }
  return df;
}

// Coordinates used to query the generator spin model. They are computed from
// the current node's kinematics: truth for GEN, reconstructed/fake-smeared for
// REC. The active nucleon carries the remaining He3 minus momentum fraction.
ROOT::RDF::RNode DefineSpinModelCoordinates(ROOT::RDF::RNode df,
                                            double beamEnergy) {
  df = df.Define("model_cos_theta_q", [beamEnergy](double q2, double nu) {
    const double qmag2=nu*nu+q2;
    return qmag2>0.0 ? std::clamp((nu+q2/(2.0*beamEnergy))/std::sqrt(qmag2),-1.0,1.0)
                     : -999.0;
  }, {"Q2","nu"});
  for (const std::string topology : {"epp","epn","ed"}) {
    df = df.Define(topology+"_alpha_active",
        [](double spectatorAlpha) { return 3.0-spectatorAlpha; },
        {topology+"_alpha"});
    // Exact invariant Bjorken x seen by the struck nucleon.  In the stationary
    // He3 frame p_active=P_He3-p_s, hence
    //   p_active.q=(M_He3-E_s)*nu+p_L,s*|q|.
    // Keep the electron-only xB column unchanged and separate.
    df = df.Define(topology+"_x_active",
        [](double q2, double nu, double spectatorEnergy,
           double spectatorLongitudinal) {
          if(!(q2>0.0) || !(nu>0.0) || !(spectatorEnergy>0.0) ||
             !std::isfinite(spectatorLongitudinal)) return -999.0;
          const double qMagnitude=std::sqrt(nu*nu+q2);
          const double activeEnergy=kHe3Mass-spectatorEnergy;
          const double activeDotQ=activeEnergy*nu+
                                  spectatorLongitudinal*qMagnitude;
          return activeDotQ>0.0 ? q2/(2.0*activeDotQ) : -999.0;
        }, {"Q2","nu",topology+"_spectator_energy",
                         topology+"_spectator_pL"});
  }
  return df;
}

// Momentum-dependent effective polarization of the struck nucleon in each
// tagged breakup channel. PN is a generator-level nuclear-structure quantity:
// both GEN and REC therefore use the unsmeared MC spectator momenta. REC/fake
// spectator kinematics remain the plotting axes and event-selection inputs,
// but detector smearing must not redefine PN for the same physical event.
ROOT::RDF::RNode DefineEffectiveNucleonPolarization(ROOT::RDF::RNode df) {
  df = df.Define("MC_spectator_p_sums",
    [](const RVec<int>& pid, const RVec<float>& px,
       const RVec<float>& py, const RVec<float>& pz) {
      const size_t n=std::min({pid.size(),px.size(),py.size(),pz.size()});
      TVector3 protonSum, neutronSum, deuteronSum;
      bool hasProton=false,hasNeutron=false,hasDeuteron=false;
      for(size_t i=0;i<n;++i) {
        const TVector3 momentum(px[i],py[i],pz[i]);
        if(IsRequestedPid(pid[i],kProtonPid)) {
          protonSum+=momentum;hasProton=true;
        } else if(IsRequestedPid(pid[i],kNeutronPid)) {
          neutronSum+=momentum;hasNeutron=true;
        } else if(IsRequestedPid(pid[i],kDeuteronRecPid)) {
          deuteronSum+=momentum;hasDeuteron=true;
        }
      }
      return RVec<double>{hasProton ? protonSum.Mag() : -999.0,
                          hasProton && hasNeutron
                              ? (protonSum+neutronSum).Mag() : -999.0,
                          hasDeuteron ? deuteronSum.Mag() : -999.0};
    }, {"MC_Particle_pid","MC_Particle_px","MC_Particle_py",
        "MC_Particle_pz"})
    .Define("epp_PN", [](const RVec<double>& p) {
      return p[0]>=0.0 ? polarization::GetPn(p[0])
                       : std::numeric_limits<double>::quiet_NaN();
    }, {"MC_spectator_p_sums"})
    .Define("epn_PN", [](const RVec<double>& p) {
      return p[1]>=0.0 ? polarization::GetPp3BBU(p[1])
                       : std::numeric_limits<double>::quiet_NaN();
    }, {"MC_spectator_p_sums"})
    .Define("ed_PN", [](const RVec<double>& p) {
      return p[2]>=0.0 ? polarization::GetPp2BBU(p[2])
                       : std::numeric_limits<double>::quiet_NaN();
    }, {"MC_spectator_p_sums"});
  return df;
}

ROOT::RDF::RNode DefineCommonKinematics(ROOT::RDF::RNode df,
                                        double beamEnergy) {
  return df
    .Define("Q2", [beamEnergy](float p, float thetaDeg) {
      const double theta = thetaDeg*M_PI/180.0;
      return p > 0.0f ? 2.0*beamEnergy*p*(1.0-std::cos(theta)) : -999.0;
    }, {"recel_p", "recel_theta"})
    .Define("nu", [beamEnergy](float p) {
      return p > 0.0f ? beamEnergy-p : -999.0;
    }, {"recel_p"})
    .Define("xB", [](double q2, double nu) {
      return nu > 0.0 ? q2/(2.0*kNucleonMass*nu) : -999.0;
    }, {"Q2", "nu"})
    .Define("WN", [](double q2, double nu) {
      const double w2 = kNucleonMass*kNucleonMass +
                        2.0*kNucleonMass*nu-q2;
      return w2 >= 0.0 ? std::sqrt(w2) : -999.0;
    }, {"Q2", "nu"})
    .Define("WHe3", [](double q2, double nu) {
      const double w2 = kHe3Mass*kHe3Mass +
                        2.0*kHe3Mass*nu-q2;
      return w2 >= 0.0 ? std::sqrt(w2) : -999.0;
    }, {"Q2", "nu"})
    // Backward-compatible name used by existing cuts and plotting code.
    .Define("W", [](double wn) { return wn; }, {"WN"})
    .Define("y", [beamEnergy](double nu) {
      return beamEnergy > 0.0 ? nu/beamEnergy : -999.0;
    }, {"nu"})
    .Define("event_weight", []() { return 1.0; });
}

ROOT::RDF::RNode DefineAllParticles(
    ROOT::RDF::RNode df, const std::vector<std::string>& inputs,
    double beamEnergy) {
  df = DefineParticleColumns(df, "ele", kElectronPid, inputs);
  df = DefineParticleColumns(df, "pro", kProtonPid, inputs);
  df = DefineParticleColumns(df, "neu", kNeutronPid, inputs);
  df = DefineParticleColumns(df, "deu", kDeuteronRecPid, inputs);

  df = DefineScalar(df, "recel_p", "ele_p_all", 0);
  df = DefineScalar(df, "recel_theta", "ele_theta_all", 0);
  df = DefineScalar(df, "recel_phi", "ele_phi_all", 0);
  df = DefineScalar(df, "recel_vz", "ele_vz_all", 0);
  for (size_t i = 0; i < 2; ++i) {
    const std::string number = std::to_string(i+1);
    df = DefineScalar(df, "pro" + number + "_p", "pro_p_all", i);
    df = DefineScalar(df, "pro" + number + "_theta", "pro_theta_all", i);
    df = DefineScalar(df, "pro" + number + "_phi", "pro_phi_all", i);
    df = DefineScalar(df, "pro" + number + "_vz", "pro_vz_all", i);
  }
  for (const auto& prefix : {std::string("neu"), std::string("deu")}) {
    df = DefineScalar(df, prefix + "_p", prefix + "_p_all", 0);
    df = DefineScalar(df, prefix + "_theta", prefix + "_theta_all", 0);
    df = DefineScalar(df, prefix + "_phi", prefix + "_phi_all", 0);
    df = DefineScalar(df, prefix + "_vz", prefix + "_vz_all", 0);
  }
  return DefineCommonKinematics(df, beamEnergy);
}

ROOT::RDF::RNode DefineFakeRecParticles(
    ROOT::RDF::RNode df, double beamEnergy,
    const FakeSpectatorConfig& fake) {
  df = DefineParticleColumns(df, "ele", kElectronPid,
      {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_py",
       "REC_Particle_pz", "REC_Particle_vz", "REC_Particle_pass_std"});

  auto defineSpecies = [&df](const std::string& prefix, int requested,
                             const FakeChargedConfig& cfg) {
    const std::string result = "fake_" + prefix + "_result";
    df = df.Define(result,
        [requested,cfg](const RVec<int>& pid, const RVec<float>& px,
                        const RVec<float>& py, const RVec<float>& pz,
                        const RVec<float>& vz, ULong64_t entry) {
          return MakeFakeCharged(pid,px,py,pz,vz,entry,requested,cfg);
        }, {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
            "MC_Particle_pz", "MC_Particle_vz", "rdfentry_"})
      .Define(result+"_pid", [](const FakeChargedParticles& r){return r.pid;}, {result})
      .Define(result+"_px", [](const FakeChargedParticles& r){return r.px;}, {result})
      .Define(result+"_py", [](const FakeChargedParticles& r){return r.py;}, {result})
      .Define(result+"_pz", [](const FakeChargedParticles& r){return r.pz;}, {result})
      .Define(result+"_vz", [](const FakeChargedParticles& r){return r.vz;}, {result})
      .Define(result+"_pass", [](const FakeChargedParticles& r){return r.pass;}, {result});
    df = DefineParticleColumns(df, prefix, requested,
        {result+"_pid",result+"_px",result+"_py",result+"_pz",
         result+"_vz",result+"_pass"});
  };
  defineSpecies("pro", kProtonPid, fake.proton);
  defineSpecies("deu", kDeuteronRecPid, fake.deuteron);

  df = df.Define("neu_p_all", [](){return RVec<float>{};})
         .Define("neu_theta_all", [](){return RVec<float>{};})
         .Define("neu_phi_all", [](){return RVec<float>{};})
         .Define("neu_vz_all", [](){return RVec<float>{};});
  df = DefineScalar(df, "recel_p", "ele_p_all", 0);
  df = DefineScalar(df, "recel_theta", "ele_theta_all", 0);
  df = DefineScalar(df, "recel_phi", "ele_phi_all", 0);
  df = DefineScalar(df, "recel_vz", "ele_vz_all", 0);
  for (size_t i=0; i<2; ++i) {
    const std::string number=std::to_string(i+1);
    df=DefineScalar(df,"pro"+number+"_p","pro_p_all",i);
    df=DefineScalar(df,"pro"+number+"_theta","pro_theta_all",i);
    df=DefineScalar(df,"pro"+number+"_phi","pro_phi_all",i);
    df=DefineScalar(df,"pro"+number+"_vz","pro_vz_all",i);
  }
  for (const auto& prefix : {std::string("neu"),std::string("deu")}) {
    df=DefineScalar(df,prefix+"_p",prefix+"_p_all",0);
    df=DefineScalar(df,prefix+"_theta",prefix+"_theta_all",0);
    df=DefineScalar(df,prefix+"_phi",prefix+"_phi_all",0);
    df=DefineScalar(df,prefix+"_vz",prefix+"_vz_all",0);
  }
  return DefineCommonKinematics(df,beamEnergy);
}

ROOT::RDF::RNode DefineFakeHe3DetectorRecParticles(
    ROOT::RDF::RNode df, double beamEnergy,
    const FakeHe3DetectorConfig& fake) {
  // The scattered electron remains the stored GEMC reconstruction.  Only the
  // charged spectator p/d collections are replaced by the He3Detector model.
  df=DefineParticleColumns(df,"ele",kElectronPid,
      {"REC_Particle_pid","REC_Particle_px","REC_Particle_py",
       "REC_Particle_pz","REC_Particle_vz","REC_Particle_pass_std"});

  auto defineSpecies=[&df,&fake](const std::string& prefix,int requested,
                                 const FakeHe3DetectorSpeciesConfig& species) {
    const std::string result="fake_he3det_"+prefix+"_result";
    df=df.Define(result,
        [requested,species,fake](const RVec<int>& pid,const RVec<float>& px,
                                 const RVec<float>& py,const RVec<float>& pz,
                                 const RVec<float>& vz,ULong64_t entry) {
          return MakeFakeHe3DetectorCharged(
              pid,px,py,pz,vz,entry,requested,species,fake);
        }, {"MC_Particle_pid","MC_Particle_px","MC_Particle_py",
            "MC_Particle_pz","MC_Particle_vz","rdfentry_"})
      .Define(result+"_pid",[](const FakeChargedParticles& r){return r.pid;},{result})
      .Define(result+"_px",[](const FakeChargedParticles& r){return r.px;},{result})
      .Define(result+"_py",[](const FakeChargedParticles& r){return r.py;},{result})
      .Define(result+"_pz",[](const FakeChargedParticles& r){return r.pz;},{result})
      .Define(result+"_vz",[](const FakeChargedParticles& r){return r.vz;},{result})
      .Define(result+"_pass",[](const FakeChargedParticles& r){return r.pass;},{result});
    df=DefineHe3TransportColumns(df,result);
    df=DefineParticleColumns(df,prefix,requested,
        {result+"_pid",result+"_px",result+"_py",result+"_pz",
         result+"_vz",result+"_pass"});
  };
  defineSpecies("pro",kProtonPid,fake.proton);
  defineSpecies("deu",kDeuteronRecPid,fake.deuteron);

  df=df.Define("neu_p_all",[](){return RVec<float>{};})
       .Define("neu_theta_all",[](){return RVec<float>{};})
       .Define("neu_phi_all",[](){return RVec<float>{};})
       .Define("neu_vz_all",[](){return RVec<float>{};});
  df=DefineScalar(df,"recel_p","ele_p_all",0);
  df=DefineScalar(df,"recel_theta","ele_theta_all",0);
  df=DefineScalar(df,"recel_phi","ele_phi_all",0);
  df=DefineScalar(df,"recel_vz","ele_vz_all",0);
  for (size_t i=0;i<2;++i) {
    const std::string number=std::to_string(i+1);
    df=DefineScalar(df,"pro"+number+"_p","pro_p_all",i);
    df=DefineScalar(df,"pro"+number+"_theta","pro_theta_all",i);
    df=DefineScalar(df,"pro"+number+"_phi","pro_phi_all",i);
    df=DefineScalar(df,"pro"+number+"_vz","pro_vz_all",i);
  }
  for (const auto& prefix:{std::string("neu"),std::string("deu")}) {
    df=DefineScalar(df,prefix+"_p",prefix+"_p_all",0);
    df=DefineScalar(df,prefix+"_theta",prefix+"_theta_all",0);
    df=DefineScalar(df,prefix+"_phi",prefix+"_phi_all",0);
    df=DefineScalar(df,prefix+"_vz",prefix+"_vz_all",0);
  }
  return DefineCommonKinematics(df,beamEnergy);
}

ROOT::RDF::RNode DefineGenAcceptanceColumns(
    ROOT::RDF::RNode df, const FakeSpectatorConfig& fake) {
  auto defineSpecies = [&df](const std::string& name, int requested,
                             const FakeChargedConfig& cfg) {
    const std::string result="acc_"+name+"_result";
    df=df.Define(result,
        [requested,cfg](const RVec<int>& pid,const RVec<float>& px,
                        const RVec<float>& py,const RVec<float>& pz,
                        const RVec<float>& vz,ULong64_t entry) {
          return MakeFakeCharged(pid,px,py,pz,vz,entry,requested,cfg);
        }, {"MC_Particle_pid","MC_Particle_px","MC_Particle_py",
            "MC_Particle_pz","MC_Particle_vz","rdfentry_"})
      .Define("acc_"+name+"_p",[](const FakeChargedParticles& r){return r.truthP;},{result})
      .Define("acc_"+name+"_theta",[](const FakeChargedParticles& r){return r.truthTheta;},{result})
      .Define("acc_"+name+"_phi",[](const FakeChargedParticles& r){return r.truthPhi;},{result})
      .Define("acc_"+name+"_weight",[](const FakeChargedParticles& r){return r.acceptance;},{result});
  };
  defineSpecies("p",kProtonPid,fake.proton);
  defineSpecies("d",kDeuteronRecPid,fake.deuteron);

  // Probability that at least two generated protons are detected.  This also
  // handles events containing more than two protons without choosing a pair.
  df=df.Define("acc_event_epp",[](const RVec<float>& a) {
        double p0=1.0;
        for (float x:a) p0*=1.0-std::clamp(double(x),0.0,1.0);
        double p1=0.0;
        for (size_t i=0;i<a.size();++i) {
          double term=std::clamp(double(a[i]),0.0,1.0);
          for (size_t j=0;j<a.size();++j)
            if (j!=i) term*=1.0-std::clamp(double(a[j]),0.0,1.0);
          p1+=term;
        }
        return std::clamp(1.0-p0-p1,0.0,1.0);
      },{"acc_p_weight"})
    .Define("acc_event_ed",[](const RVec<float>& a) {
        double none=1.0;
        for (float x:a) none*=1.0-std::clamp(double(x),0.0,1.0);
        return std::clamp(1.0-none,0.0,1.0);
      },{"acc_d_weight"})
    .Define("acc_epp_p",[](const RVec<float>& p){return AtLargestMomentum(p,p);},{"acc_p_p"})
    .Define("acc_epp_theta",[](const RVec<float>& p,const RVec<float>& v){return AtLargestMomentum(p,v);},{"acc_p_p","acc_p_theta"})
    .Define("acc_epp_phi",[](const RVec<float>& p,const RVec<float>& v){return AtLargestMomentum(p,v);},{"acc_p_p","acc_p_phi"})
    .Define("acc_ed_p",[](const RVec<float>& v){return AtOrMissing(v,0);},{"acc_d_p"})
    .Define("acc_ed_theta",[](const RVec<float>& v){return AtOrMissing(v,0);},{"acc_d_theta"})
    .Define("acc_ed_phi",[](const RVec<float>& v){return AtOrMissing(v,0);},{"acc_d_phi"});
  return df;
}

ROOT::RDF::RNode DefineGenHe3DetectorAcceptanceColumns(
    ROOT::RDF::RNode df,const FakeHe3DetectorConfig& fake) {
  auto defineSpecies=[&df,&fake](const std::string& name,int requested,
                                 const FakeHe3DetectorSpeciesConfig& species) {
    const std::string result="acc_he3det_"+name+"_result";
    df=df.Define(result,
        [requested,species,fake](const RVec<int>& pid,const RVec<float>& px,
                                 const RVec<float>& py,const RVec<float>& pz,
                                 const RVec<float>& vz,ULong64_t entry) {
          return MakeFakeHe3DetectorCharged(
              pid,px,py,pz,vz,entry,requested,species,fake);
        }, {"MC_Particle_pid","MC_Particle_px","MC_Particle_py",
            "MC_Particle_pz","MC_Particle_vz","rdfentry_"})
      .Define("acc_"+name+"_p",[](const FakeChargedParticles& r){return r.truthP;},{result})
      .Define("acc_"+name+"_theta",[](const FakeChargedParticles& r){return r.truthTheta;},{result})
      .Define("acc_"+name+"_phi",[](const FakeChargedParticles& r){return r.truthPhi;},{result})
      .Define("acc_"+name+"_weight",[](const FakeChargedParticles& r){return r.acceptance;},{result});
    df=DefineHe3TransportColumns(df,result);
  };
  defineSpecies("p",kProtonPid,fake.proton);
  defineSpecies("d",kDeuteronRecPid,fake.deuteron);

  df=df.Define("acc_event_epp",[](const RVec<float>& a) {
        double p0=1.0;
        for (float x:a) p0*=1.0-std::clamp(double(x),0.0,1.0);
        double p1=0.0;
        for (size_t i=0;i<a.size();++i) {
          double term=std::clamp(double(a[i]),0.0,1.0);
          for (size_t j=0;j<a.size();++j)
            if (j!=i) term*=1.0-std::clamp(double(a[j]),0.0,1.0);
          p1+=term;
        }
        return std::clamp(1.0-p0-p1,0.0,1.0);
      },{"acc_p_weight"})
    .Define("acc_event_ed",[](const RVec<float>& a) {
        double none=1.0;
        for (float x:a) none*=1.0-std::clamp(double(x),0.0,1.0);
        return std::clamp(1.0-none,0.0,1.0);
      },{"acc_d_weight"})
    .Define("acc_epp_p",[](const RVec<float>& p){return AtLargestMomentum(p,p);},{"acc_p_p"})
    .Define("acc_epp_theta",[](const RVec<float>& p,const RVec<float>& v){return AtLargestMomentum(p,v);},{"acc_p_p","acc_p_theta"})
    .Define("acc_epp_phi",[](const RVec<float>& p,const RVec<float>& v){return AtLargestMomentum(p,v);},{"acc_p_p","acc_p_phi"})
    .Define("acc_ed_p",[](const RVec<float>& v){return AtOrMissing(v,0);},{"acc_d_p"})
    .Define("acc_ed_theta",[](const RVec<float>& v){return AtOrMissing(v,0);},{"acc_d_theta"})
    .Define("acc_ed_phi",[](const RVec<float>& v){return AtOrMissing(v,0);},{"acc_d_phi"});
  return df;
}
}  // namespace He3DISPlot

struct He3DISDataFrames {
  ROOT::RDF::RNode all, epp, epn, ed;
};

He3DISDataFrames InitGenKinematicsHe3DIS(
    const std::string& file, const std::string& tree = "dfSelectedMC",
    double beamEnergy = 10.6,
    const He3DISPlot::FakeSpectatorConfig& fakeAlert = {},
    const He3DISPlot::FakeHe3DetectorConfig& fakeHe3Detector = {},
    He3DISPlot::FakeDetectorMode detectorMode =
        He3DISPlot::FakeDetectorMode::ALERT) {
  using namespace He3DISPlot;
  ROOT::RDataFrame input(tree, file);
  ROOT::RDF::RNode df(input);
  RequireColumns(df, {"MC_Particle_pid", "MC_Particle_px",
      "MC_Particle_py", "MC_Particle_pz", "MC_Particle_vz"});
  // Use explicit lambdas because each count has a different requested PID.
  df = df.Define("GEN_nProton",
      [](const RVec<int>& pid) { return CountPid(pid, kProtonPid); },
      {"MC_Particle_pid"})
    .Define("GEN_nNeutron",
      [](const RVec<int>& pid) { return CountPid(pid, kNeutronPid); },
      {"MC_Particle_pid"})
    .Define("GEN_nDeuteron",
      [](const RVec<int>& pid) { return CountPid(pid, kDeuteronRecPid); },
      {"MC_Particle_pid"})
    .Define("GEN_pass_epp", "GEN_nProton >= 2")
    .Define("GEN_pass_epn", "GEN_nProton >= 1 && GEN_nNeutron >= 1")
    .Define("GEN_pass_ed", "GEN_nDeuteron >= 1")
    .Define("truth_pass",
      [](const RVec<int>& pid) { return RVec<bool>(pid.size(), true); },
      {"MC_Particle_pid"});
  df = DefineAllParticles(df, {"MC_Particle_pid", "MC_Particle_px",
      "MC_Particle_py", "MC_Particle_pz", "MC_Particle_vz", "truth_pass"},
      beamEnergy);
  if (detectorMode==FakeDetectorMode::He3Detector)
    df=DefineGenHe3DetectorAcceptanceColumns(df,fakeHe3Detector);
  else
    df=DefineGenAcceptanceColumns(df,fakeAlert);
  // GEN remains the undiluted, 100%-polarization physics reference. Partial
  // beam/target polarization is an experimental REC effect only.
  df = DefineProjectionSpin(df,1.0,1.0);
  df = DefineSpectatorKinematics(df, beamEnergy);
  df = DefineEffectiveNucleonPolarization(df);
  df = DefineSpinModelCoordinates(df, beamEnergy);
  df = df.Filter(
      "recel_p > 2.0 && (GEN_pass_epp || GEN_pass_epn || GEN_pass_ed)",
      "GEN He3DIS: electron p > 2 GeV and epp OR epn OR ed");
  df = df.Filter("Q2 > 1.0", "GEN DIS Q2 > 1 GeV^2");
  // W is WN (nucleon mass convention), not the whole-He3 invariant mass.
  df = df.Filter("W > 2.0", "GEN DIS W > 2 GeV");
  auto epp=df.Filter("GEN_pass_epp")
      .Filter([](double pSum) {
          return pSum>=0.0 && pSum<kEppSpectatorPSumMaxGeV;
        }, {"epp_p_sum"}, "GEN epp |p_s1+p_s2| < 0.5 GeV/c");
  auto ed=df.Filter("GEN_pass_ed")
      .Filter([](double pd) {
          return pd>=0.0 && pd<kEdSpectatorMomentumMaxGeV;
        }, {"ed_p_sum"}, "GEN ed p_d < 0.5 GeV/c");
  return {df, epp, df.Filter("GEN_pass_epn"), ed};
}

He3DISDataFrames InitRecKinematicsHe3DIS(
    const std::string& file, const std::string& tree = "dfSelected_afterFid",
    double beamEnergy = 10.6,
    const He3DISPlot::FakeSpectatorConfig& fakeAlert = {},
    const He3DISPlot::FakeHe3DetectorConfig& fakeHe3Detector = {},
    He3DISPlot::FakeDetectorMode detectorMode =
        He3DISPlot::FakeDetectorMode::ALERT,
    double beamPolarization = 0.80,
    double targetPolarization = 0.50) {
  using namespace He3DISPlot;
  auto df = NormalizeRecPass(file, tree);
  RequireColumns(df, {"REC_Particle_pid", "REC_Particle_px",
      "REC_Particle_py", "REC_Particle_pz", "REC_Particle_vz",
      "REC_Particle_status", "REC_Particle_pass", "REC_Event_pass_epp",
      "REC_Event_pass_epn", "REC_Event_pass_ed",
      "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
      "MC_Particle_pz"});
  // Explicit plotting-level detector requirement.  Under the standard CLAS12
  // status convention, abs(status) in [2000,3000) identifies the FD.
  df = df.Define("REC_electron_in_FD",
      [](const RVec<int>& pid, const RVec<short>& status,
         const RVec<bool>& pass) {
        const size_t n = std::min({pid.size(),status.size(),pass.size()});
        for (size_t i=0; i<n; ++i) {
          const int s=std::abs(static_cast<int>(status[i]));
          if (pass[i] && pid[i]==kElectronPid && s>=2000 && s<3000)
            return true;
        }
        return false;
      }, {"REC_Particle_pid","REC_Particle_status","REC_Particle_pass_std"})
    .Filter("REC_electron_in_FD", "REC electron in Forward Detector");
  if (detectorMode==FakeDetectorMode::ALERT) {
    RequireColumns(df, {"MC_Particle_px", "MC_Particle_py",
        "MC_Particle_pz", "MC_Particle_vz"});
    df=DefineFakeRecParticles(df,beamEnergy,fakeAlert);
  } else if (detectorMode==FakeDetectorMode::He3Detector) {
    RequireColumns(df,{"MC_Particle_px","MC_Particle_py",
        "MC_Particle_pz","MC_Particle_vz"});
    df=DefineFakeHe3DetectorRecParticles(df,beamEnergy,fakeHe3Detector);
  } else {
    df = DefineAllParticles(df, {"REC_Particle_pid", "REC_Particle_px",
        "REC_Particle_py", "REC_Particle_pz", "REC_Particle_vz",
        "REC_Particle_pass_std"}, beamEnergy);
  }
  df = df
    .Define("GEN_nProton",
      [](const RVec<int>& pid) { return CountPid(pid, kProtonPid); },
      {"MC_Particle_pid"})
    .Define("GEN_nNeutron",
      [](const RVec<int>& pid) { return CountPid(pid, kNeutronPid); },
      {"MC_Particle_pid"})
    .Define("GEN_nDeuteron",
      [](const RVec<int>& pid) { return CountPid(pid, kDeuteronRecPid); },
      {"MC_Particle_pid"})
    .Define("GEN_pass_epp", "GEN_nProton >= 2")
    .Define("GEN_pass_epn", "GEN_nProton >= 1 && GEN_nNeutron >= 1")
    .Define("GEN_pass_ed", "GEN_nDeuteron >= 1")
    .Filter("recel_p > 2", "REC He3DIS electron exists");
  if (detectorMode!=FakeDetectorMode::GEMC) {
    df = df.Define("fake_REC_pass_epp",
                   [](const RVec<float>& p){return p.size()>=2;}, {"pro_p_all"})
           .Define("fake_REC_pass_epn", [](){return false;})
           .Define("fake_REC_pass_ed",
                   [](const RVec<float>& p){return !p.empty();}, {"deu_p_all"})
           .Redefine("REC_Event_pass_epp", "fake_REC_pass_epp")
           .Redefine("REC_Event_pass_epn", "fake_REC_pass_epn")
           .Redefine("REC_Event_pass_ed", "fake_REC_pass_ed");
  }
  df = DefineSpectatorKinematics(df, beamEnergy);
  df = DefineEffectiveNucleonPolarization(df);
  df = DefineSpinModelCoordinates(df, beamEnergy);
  df = DefineProjectionSpin(df,beamPolarization,targetPolarization);
  df = df.Filter("Q2 > 1.0", "REC DIS Q2 > 1 GeV^2");
  // Apply the same DIS threshold using the reconstructed electron's W.
  df = df.Filter("W > 2.0", "REC DIS W > 2 GeV");
  auto epp=df.Filter("REC_Event_pass_epp")
      .Filter([](double pSum) {
          return pSum>=0.0 && pSum<kEppSpectatorPSumMaxGeV;
        }, {"epp_p_sum"}, "REC epp |p_s1+p_s2| < 0.5 GeV/c");
  auto ed=df.Filter("REC_Event_pass_ed")
      .Filter([](double pd) {
          return pd>=0.0 && pd<kEdSpectatorMomentumMaxGeV;
        }, {"ed_p_sum"}, "REC ed p_d < 0.5 GeV/c");
  return {df, epp, df.Filter("REC_Event_pass_epn"),
          ed};
}

using He3DISModelEvaluator =
    std::function<double(bool,double,double,double,double,double,double)>;
using He3DISNucleonEvaluator =
    std::function<std::array<double,7>(bool,double,double,double)>;
using He3DISG1UncertaintyEvaluator =
    std::function<std::array<double,4>(bool,double,double,double)>;

struct He3DISGeneratorModel {
  He3DISModelEvaluator aparallel;
  He3DISNucleonEvaluator nucleon;
  He3DISG1UncertaintyEvaluator g1Uncertainty;
  explicit operator bool() const { return bool(aparallel); }
};

He3DISGeneratorModel LoadHe3DISGeneratorModel() {
  using ModelFunction = int (*)(int,double,double,double,double,double,double,
                                double*,double*);
  using NucleonFunction = int (*)(int,double,double,double,double*,double*,
                                  double*,double*,double*,double*,double*);
  using G1UncertaintyFunction = int (*)(int,double,double,double,double*,
                                        double*,double*,double*);
  std::vector<std::string> candidates;
  if (const char* configured=std::getenv("HE3DIS_SPIN_MODEL_LIB"))
    candidates.emplace_back(configured);
  const auto generatorRoot=(std::filesystem::absolute(__FILE__).parent_path()/
      "../../../../Simulation/He3DIS/He3DIS").lexically_normal();
  candidates.push_back((generatorRoot/"build-ifarm/libhe3dis_spin_model.so").string());
  candidates.push_back((generatorRoot/"build/libhe3dis_spin_model.dylib").string());
  candidates.push_back((generatorRoot/"build/libhe3dis_spin_model.so").string());

  void* handle=nullptr;
  std::string loadedPath;
  for (const auto& path:candidates) {
    handle=dlopen(path.c_str(),RTLD_NOW|RTLD_LOCAL);
    if(handle) {loadedPath=path;break;}
  }
  if(!handle) {
    std::cerr<<"[He3DIS generator model] shared library not found. Build it with:\n"
      <<"  cd "<<generatorRoot<<"\n"
      <<"  source setup_ifarm.csh build\n"
      <<"or set HE3DIS_SPIN_MODEL_LIB. Existing count projection is retained.\n";
    return {};
  }
  dlerror();
  auto* function=reinterpret_cast<ModelFunction>(dlsym(handle,"he3dis_spin_response"));
  if(const char* error=dlerror();error || !function) {
    std::cerr<<"[He3DIS generator model] missing he3dis_spin_response in "
             <<loadedPath<<": "<<(error?error:"unknown error")<<'\n';
    dlclose(handle);return {};
  }
  dlerror();
  auto* nucleonFunction=reinterpret_cast<NucleonFunction>(
      dlsym(handle,"he3dis_nucleon_observables"));
  if(const char* error=dlerror();error || !nucleonFunction) {
    std::cerr<<"[He3DIS generator model] missing "
                "he3dis_nucleon_observables in "<<loadedPath<<": "
             <<(error?error:"unknown error")
             <<". Rebuild libhe3dis_spin_model before drawing A1.\n";
    nucleonFunction=nullptr;
  }
  dlerror();
  auto* g1UncertaintyFunction=reinterpret_cast<G1UncertaintyFunction>(
      dlsym(handle,"he3dis_g1_uncertainty"));
  if(const char* error=dlerror();error || !g1UncertaintyFunction) {
    std::cerr<<"[He3DIS generator model] missing he3dis_g1_uncertainty in "
             <<loadedPath<<": "<<(error?error:"unknown error")
             <<". Rebuild libhe3dis_spin_model before drawing polarized EMC.\n";
    g1UncertaintyFunction=nullptr;
  }
  // Deliberately retain the handle: the callback is used later by comparer.
  std::cout<<"[He3DIS generator model] loaded "<<loadedPath<<'\n';
  He3DISGeneratorModel model;
  model.aparallel=[function,loadedPath](bool neutron,double xb,double alpha,
                                       double q2,double y,double spar,
                                       double sperpCos) {
    double unpolarized=0.,asymmetry=std::numeric_limits<double>::quiet_NaN();
    const int status=function(neutron?1:0,xb,alpha,q2,y,spar,sperpCos,
                              &unpolarized,&asymmetry);
    if(status!=0) {
      // Print only a few examples so an empty graph is diagnosable without
      // flooding the terminal once for every plot and every bin.
      static std::array<unsigned int,9> reported{};
      const size_t index=status>=0 && status<static_cast<int>(reported.size())
                           ? static_cast<size_t>(status) : 0U;
      if(reported[index]++<3U)
        std::cerr<<"[He3DIS generator model] status "<<status
                 <<" for neutron="<<neutron<<", xB="<<xb
                 <<", alpha="<<alpha<<", Q2="<<q2<<", y="<<y
                 <<" ("<<loadedPath<<")\n";
      return std::numeric_limits<double>::quiet_NaN();
    }
    return asymmetry;
  };
  if(nucleonFunction) model.nucleon=[nucleonFunction,loadedPath](bool neutron,
      double xb,double alpha,double q2) {
    std::array<double,7> values;
    values.fill(std::numeric_limits<double>::quiet_NaN());
    const int status=nucleonFunction(neutron?1:0,xb,alpha,q2,
        &values[0],&values[1],&values[2],&values[3],&values[4],
        &values[5],&values[6]);
    if(status!=0) {
      static std::array<unsigned int,7> reported{};
      const size_t index=status>=0 && status<static_cast<int>(reported.size())
                           ? static_cast<size_t>(status) : 0U;
      if(reported[index]++<3U)
        std::cerr<<"[He3DIS nucleon model] status "<<status
                 <<" for neutron="<<neutron<<", xB="<<xb
                 <<", alpha="<<alpha<<", Q2="<<q2
                 <<" ("<<loadedPath<<")\n";
      values.fill(std::numeric_limits<double>::quiet_NaN());
    }
    return values;
  };
  if(g1UncertaintyFunction)
    model.g1Uncertainty=[g1UncertaintyFunction,loadedPath](bool neutron,
        double xb,double alpha,double q2) {
      std::array<double,4> values;
      values.fill(std::numeric_limits<double>::quiet_NaN());
      const int status=g1UncertaintyFunction(neutron?1:0,xb,alpha,q2,
          &values[0],&values[1],&values[2],&values[3]);
      if(status!=0) {
        static std::array<unsigned int,10> reported{};
        const size_t index=status>=0 && status<static_cast<int>(reported.size())
                             ? static_cast<size_t>(status) : 0U;
        if(reported[index]++<3U)
          std::cerr<<"[He3DIS JAM25 uncertainty] status "<<status
                   <<" for neutron="<<neutron<<", xB="<<xb
                   <<", alpha="<<alpha<<", Q2="<<q2
                   <<" ("<<loadedPath<<")\n";
        values.fill(std::numeric_limits<double>::quiet_NaN());
      }
      return values;
    };
  return model;
}

void DISANA_XplotterHe3DIS(
    const std::string& directory = "../build/He3DISrate",
    bool enable_fake_spectators = true,
    double projection_to_count_ratio = 34.46*4.0,
    double beam_polarization = 0.80,
    double target_polarization = 0.50,
    bool apply_pn_to_legacy_projection = true,
    const std::string& fake_detector = "He3Detector",
    int he3det_module_count = 1,
    double he3det_module_width_deg = 30.0,
    double he3det_module_offset_deg = 0.0) {
  ROOT::EnableImplicitMT(40);
  constexpr double beamEnergy = 10.6;
  He3DISPlot::FakeDetectorMode detectorMode=enable_fake_spectators
      ? He3DISPlot::ParseFakeDetectorMode(fake_detector)
      : He3DISPlot::FakeDetectorMode::GEMC;
  He3DISPlot::FakeSpectatorConfig fakeAlert=
      He3DISPlot::gFakeSpectatorConfig;
  fakeAlert.enabled=detectorMode==He3DISPlot::FakeDetectorMode::ALERT;
  He3DISPlot::FakeHe3DetectorConfig fakeHe3Detector=
      He3DISPlot::gFakeHe3DetectorConfig;
  fakeHe3Detector.enabled=
      detectorMode==He3DISPlot::FakeDetectorMode::He3Detector;
  fakeHe3Detector.moduleCount=he3det_module_count;
  fakeHe3Detector.moduleWidthRad=he3det_module_width_deg*M_PI/180.0;
  fakeHe3Detector.moduleOffsetRad=he3det_module_offset_deg*M_PI/180.0;
  if (fakeHe3Detector.enabled) {
    He3DISPlot::ValidateFakeHe3DetectorConfig(fakeHe3Detector);
    (void)He3DISPlot::GetHe3PStarTables(); // Load once before event loops.
  }
  auto gen = InitGenKinematicsHe3DIS(
      directory+"/dfSelectedMC.root","dfSelectedMC",beamEnergy,
      fakeAlert,fakeHe3Detector,detectorMode);
  auto rec = InitRecKinematicsHe3DIS(
      directory+"/dfSelected_afterFid.root","dfSelected_afterFid",
      beamEnergy,fakeAlert,fakeHe3Detector,detectorMode,
      beam_polarization,target_polarization);

  DISANAcomparer comparer;
  comparer.SetHe3DISXActiveBinning(
      He3DISPlot::kHe3DISXActiveBins,
      He3DISPlot::kHe3DISXActiveMin,
      He3DISPlot::kHe3DISXActiveMax,
      He3DISPlot::kHe3DISModelXMax);
  const std::string outputDirectory=
      detectorMode==He3DISPlot::FakeDetectorMode::He3Detector
          ? "He3DIS_gen_rec_plots_He3Detector"
          : detectorMode==He3DISPlot::FakeDetectorMode::GEMC
              ? "He3DIS_gen_rec_plots_GEMC"
              : "He3DIS_gen_rec_plots";
  comparer.SetOutputDir(outputDirectory);
  comparer.AddModelHe3DIS(gen.all, gen.epp, gen.epn, gen.ed, "GEN");
  comparer.AddModelHe3DIS(rec.all,rec.epp,rec.epn,rec.ed,
      std::string("REC ")+He3DISPlot::FakeDetectorModeName(detectorMode));
  //comparer.PlotNormalizedComparisonHe3DIS();
  //comparer.PlotHe3DISAcc(gen.all, rec.all);
  //comparer.PlotTopologyMigrationHe3DIS(rec.all);
  // Keep the finite-count projection and write the generator-model centers
  // to a separate directory/PDF set. Run the count projection first because
  // it defines the y-axis ranges reused by the generator-model plots.
  comparer.PlotProjection(projection_to_count_ratio,beam_polarization,
                          target_polarization,
                          apply_pn_to_legacy_projection);
  if(auto model=LoadHe3DISGeneratorModel()) {
    comparer.PlotGeneratorModelProjection(model.aparallel,projection_to_count_ratio,
        beam_polarization,target_polarization,
        apply_pn_to_legacy_projection);
    if(model.nucleon)
      comparer.PlotA1ProjectionHe3DIS(model.aparallel,model.nucleon,
                                      projection_to_count_ratio);
    if(model.nucleon)
      comparer.PlotG1ProjectionHe3DIS(model.aparallel,model.nucleon,
                                      projection_to_count_ratio);
    if(model.nucleon && model.g1Uncertainty) {
      const auto emcRatioTable=(std::filesystem::absolute(__FILE__).parent_path()/
          "g1p_g1n_in_He3.dat").lexically_normal();
      comparer.PlotPolarizedEMCEd(model.nucleon,model.g1Uncertainty,
                                  emcRatioTable.string(),
                                  projection_to_count_ratio);
    }
  }

  auto ngen = gen.all.Count(), ngenEpp = gen.epp.Count();
  auto ngenEpn = gen.epn.Count(), ngenEd = gen.ed.Count();
  // Report the union of detectable tagged topologies, counting each event
  // once. Keep rec.all electron-inclusive for the existing plots/acceptance.
  auto nrec = rec.all.Filter(
      "REC_Event_pass_epp || REC_Event_pass_epn || REC_Event_pass_ed",
      "REC detectable tagged topology (reported count only)").Count();
  auto nrecEpp = rec.epp.Count();
  auto nrecEpn = rec.epn.Count(), nrecEd = rec.ed.Count();
  std::cout << "He3DIS, Ebeam=" << beamEnergy << " GeV\n"
            << "GEN all/epp/epn/ed: " << ngen.GetValue() << " / "
            << ngenEpp.GetValue() << " / " << ngenEpn.GetValue() << " / "
            << ngenEd.GetValue() << "\nREC all (tagged)/epp/epn/ed: "
            << nrec.GetValue() << " / " << nrecEpp.GetValue() << " / "
            << nrecEpn.GetValue() << " / " << nrecEd.GetValue()
            << "\nSpectator detector model: "
            << He3DISPlot::FakeDetectorModeName(detectorMode)
            << "; neutron detection: OFF\n"
            << "Projection polarization: Pb=" << beam_polarization
            << ", Pt=" << target_polarization
            << ", Pb*Pt=" << beam_polarization*target_polarization
            << "; plots show corrected A_parallel only\n"
            << "PN projection for legacy counts: "
            << (apply_pn_to_legacy_projection?"ON":"OFF (input already includes PN)")
            << "; plotting A_parallel(tagged) and A_parallel/PN\n";
  if (fakeHe3Detector.enabled) {
    const double spacingDeg=360.0/fakeHe3Detector.moduleCount;
    const double gapDeg=spacingDeg-he3det_module_width_deg;
    std::cout
        << "Fake He3Detector: Bz=" << fakeHe3Detector.magneticFieldT
        << " T, modules=" << fakeHe3Detector.moduleCount
        << ", width=" << he3det_module_width_deg
        << " deg, offset=" << he3det_module_offset_deg << " deg\n"
        << "  active phi="
        << fakeHe3Detector.moduleCount*he3det_module_width_deg
        << " deg, gap/module=" << gapDeg << " deg\n"
        << "  detection requires window and Si layers through L"
        << fakeHe3Detector.requiredLayer << " (stopping in that layer allowed)\n"
        << "  required window/layers in same module: "
        << (fakeHe3Detector.requireSameModule?"YES":"NO") << "\n"
        << "  radii [cm]: target=" << 100.0*fakeHe3Detector.targetRadiusM
        << ", L1=" << 100.0*fakeHe3Detector.layers[0].radiusM
        << ", L2=" << 100.0*fakeHe3Detector.layers[1].radiusM
        << ", L3=" << 100.0*fakeHe3Detector.layers[2].radiusM << "\n"
        << "  instantaneous PSTAR energy loss: Al window="
        << 1.e6*fakeHe3Detector.windowThicknessM
        << " um, non-window cell wall="
        << 1.e6*fakeHe3Detector.cellWallThicknessM << " um\n"
        << "  Si thicknesses [um]: "
        << 1.e6*fakeHe3Detector.layers[0].thicknessM << " / "
        << 1.e6*fakeHe3Detector.layers[1].thicknessM << " / "
        << 1.e6*fakeHe3Detector.layers[2].thicknessM << "\n"
        << "  sequential helix from each crossing; gas/dead layers/scattering OFF\n"
        << "  REC vertex momentum: retained truth-smearing assumption\n"
        << "  PSTAR tables: " << He3DISPlot::GetHe3PStarTables().aluminum.source
        << " ; " << He3DISPlot::GetHe3PStarTables().silicon.source << "\n"
        << "  output directory: " << outputDirectory << "\n";
  }
}

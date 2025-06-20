#ifndef MOMENTUM_CORRECTION_H
#define MOMENTUM_CORRECTION_H

#include <functional>
#include <map>
#include <vector>
#include <memory>

class MomentumCorrection {
public:
  enum class DetectorRegion { ANY, FT, FD, CD };

  static constexpr DetectorRegion FT = DetectorRegion::FT;
  static constexpr DetectorRegion FD = DetectorRegion::FD;
  static constexpr DetectorRegion CD = DetectorRegion::CD;
  static constexpr DetectorRegion ANY = DetectorRegion::ANY;


  struct RegionWithDetector {
    double pmin, pmax;
    double thetamin, thetamax;
    double phimin, phimax;
    DetectorRegion detector;
  };

  using CorrectionFunction = std::function<double(double p, double theta, double phi)>;
  using RECExtendStoreType = std::function<std::vector<float>(
    const std::vector<int>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<short>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<short>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&
  )>;

  struct RegionCorrection {
    RegionWithDetector region;
    CorrectionFunction func;
  };

  void AddPiecewiseCorrection(int pid, const RegionWithDetector& region, CorrectionFunction func);

  RECExtendStoreType RECParticlePxCorrected() const;
  RECExtendStoreType RECParticlePyCorrected() const;
  RECExtendStoreType RECParticlePzCorrected() const;

private:
  std::map<int, std::vector<RegionCorrection>> p_corrections_;

  double GetCorrectedP(int pid, double p, double theta, double phi, short status) const;
  static bool InRegion(const RegionWithDetector& region, double p, double theta, double phi, short status);
};

#endif  // MOMENTUM_CORRECTION_H

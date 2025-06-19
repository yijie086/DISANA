#ifndef MOMENTUM_CORRECTION_H
#define MOMENTUM_CORRECTION_H

#include <functional>
#include <map>
#include <tuple>
#include <vector>
#include <memory>

class MomentumCorrection {
public:
  using Region = std::tuple<double, double, double, double, double, double>;
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
    Region region;
    CorrectionFunction func;
  };

  void AddPiecewiseCorrection(int pid, const Region& region, CorrectionFunction func);

  RECExtendStoreType RECParticlePxCorrected() const;
  RECExtendStoreType RECParticlePyCorrected() const;
  RECExtendStoreType RECParticlePzCorrected() const;

private:
  std::map<int, std::vector<RegionCorrection>> p_corrections_;

  double GetCorrectedP(int pid, double p, double theta, double phi) const;
  static bool InRegion(const Region& region, double p, double theta, double phi);
};

#endif  // MOMENTUM_CORRECTION_H

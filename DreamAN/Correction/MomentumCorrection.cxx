#include "MomentumCorrection.h"
#include <cmath>
#include <iostream>

void MomentumCorrection::AddPiecewiseCorrection(int pid, const Region& region, CorrectionFunction func) {
  p_corrections_[pid].emplace_back(RegionCorrection{region, func});
}

double MomentumCorrection::GetCorrectedP(int pid, double p, double theta, double phi) const {
  auto it = p_corrections_.find(pid);
  if (it == p_corrections_.end()) return p;

  for (const auto& rc : it->second) {
    if (InRegion(rc.region, p, theta, phi)) {
      return rc.func(p, theta, phi);
    }
  }
  return p;
}

bool MomentumCorrection::InRegion(const Region& region, double p, double theta, double phi) {
  auto [pmin, pmax, tmin, tmax, phimin, phimax] = region;
  return (p >= pmin && p < pmax &&
          theta >= tmin && theta < tmax &&
          phi >= phimin && phi < phimax);
}

MomentumCorrection::RECExtendStoreType MomentumCorrection::RECParticlePxCorrected() const {
  return [this](const std::vector<int>& pid,
                const std::vector<float>& /*px*/,
                const std::vector<float>& /*py*/,
                const std::vector<float>& /*pz*/,
                const std::vector<float>& /*vx*/,
                const std::vector<float>& /*vy*/,
                const std::vector<float>& /*vz*/,
                const std::vector<float>& /*vt*/,
                const std::vector<short>& /*charge*/,
                const std::vector<float>& /*beta*/,
                const std::vector<float>& /*chi2pid*/,
                const std::vector<short>& /*status*/,
                const std::vector<float>& phi,
                const std::vector<float>& theta,
                const std::vector<float>& p) -> std::vector<float> {
    std::vector<float> result(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      float p_corr = GetCorrectedP(pid[i], p[i], theta[i], phi[i]);
      result[i] = p_corr * std::sin(theta[i]) * std::cos(phi[i]);
    }
    return result;
  };
}

MomentumCorrection::RECExtendStoreType MomentumCorrection::RECParticlePyCorrected() const {
  return [this](const std::vector<int>& pid,
                const std::vector<float>& /*px*/,
                const std::vector<float>& /*py*/,
                const std::vector<float>& /*pz*/,
                const std::vector<float>& /*vx*/,
                const std::vector<float>& /*vy*/,
                const std::vector<float>& /*vz*/,
                const std::vector<float>& /*vt*/,
                const std::vector<short>& /*charge*/,
                const std::vector<float>& /*beta*/,
                const std::vector<float>& /*chi2pid*/,
                const std::vector<short>& /*status*/,
                const std::vector<float>& phi,
                const std::vector<float>& theta,
                const std::vector<float>& p) -> std::vector<float> {
    std::vector<float> result(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      float p_corr = GetCorrectedP(pid[i], p[i], theta[i], phi[i]);
      result[i] = p_corr * std::sin(theta[i]) * std::sin(phi[i]);
    }
    return result;
  };
}

MomentumCorrection::RECExtendStoreType MomentumCorrection::RECParticlePzCorrected() const {
  return [this](const std::vector<int>& pid,
                const std::vector<float>& /*px*/,
                const std::vector<float>& /*py*/,
                const std::vector<float>& /*pz*/,
                const std::vector<float>& /*vx*/,
                const std::vector<float>& /*vy*/,
                const std::vector<float>& /*vz*/,
                const std::vector<float>& /*vt*/,
                const std::vector<short>& /*charge*/,
                const std::vector<float>& /*beta*/,
                const std::vector<float>& /*chi2pid*/,
                const std::vector<short>& /*status*/,
                const std::vector<float>& phi,
                const std::vector<float>& theta,
                const std::vector<float>& p) -> std::vector<float> {
    std::vector<float> result(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      float p_corr = GetCorrectedP(pid[i], p[i], theta[i], phi[i]);
      result[i] = p_corr * std::cos(theta[i]);
    }
    return result;
  };
}

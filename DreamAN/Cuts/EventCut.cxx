
#include "EventCut.h"
#include "./../Math/ParticleMassTable.h"

#include <cmath>
#include <iostream>

EventCut::EventCut() = default;
EventCut::~EventCut() = default;

void EventCut::AddParticleCut(const std::string& name, const ParticleCut& userCut) {
  ParticleCut cut = userCut;

  if (name == "proton") {
    cut.pid = 2212;
    cut.charge = 1;
    cut.minVz = -8;
    cut.maxVz = 2;

  } else if (name == "electron") {
    cut.pid = 11;
    cut.charge = -1;
    cut.minVz = -8;
    cut.maxVz = 2;
  } else if (name == "photon") {
    cut.pid = 22;
    cut.minVz = -8;
    cut.maxVz = 2;
    cut.minBeta = 0.9;  // Photon beta cut
    cut.maxBeta = 1.1;  // Photon beta cut
  }else if (name == "Neg Kaon") {
    cut.pid = -321;
    cut.minVz = -8;
    cut.maxVz = 2;
  } else if (name == "Pos Kaon") {
    cut.pid = 321;
    cut.minVz = -8;
    cut.maxVz = 2;
  }

  fParticleCuts[name] = cut;
}

void EventCut::AddParticleMotherCut(const std::string& name, const TwoBodyMotherCut& userCut) {
  TwoBodyMotherCut cut = userCut;
  if (name == "pi0") {
    cut.pidDaug1 = 22;
    cut.pidDaug2 = 22;
  } else if (name == "phi") {
    cut.pidDaug1 = 321;
    cut.pidDaug2 = -321;
  } else if (name == "rho") {
    cut.pidDaug1 = 211;
    cut.pidDaug2 = -211;
  }

  fTwoBodyMotherCuts[name] = cut;
}

const ParticleCut* EventCut::GetParticleCut(const std::string& name) const {
  auto it = fParticleCuts.find(name);
  return (it != fParticleCuts.end()) ? &it->second : nullptr;
}

EventCut* EventCut::ProtonCuts() {
  EventCut* cuts = new EventCut();
  cuts->AddParticleCut("proton", ParticleCut());
  return cuts;
}

EventCut* EventCut::ElectronCuts() {
  EventCut* cuts = new EventCut();
  cuts->AddParticleCut("electron", ParticleCut());
  return cuts;
}

EventCut* EventCut::PhotonCuts() {
  EventCut* cuts = new EventCut();
  cuts->AddParticleCut("photon", ParticleCut());
  return cuts;
}

EventCutResult EventCut::operator()(const std::vector<int>& pid, const std::vector<float>& px, const std::vector<float>& py, const std::vector<float>& pz,
                                    const std::vector<float>& vx, const std::vector<float>& vy, const std::vector<float>& vz, const std::vector<float>& vt,
                                    const std::vector<short>& charge, const std::vector<float>& beta, const std::vector<float>& chi2pid, const std::vector<short>& status,
                                    const std::vector<int>& REC_Track_pass_fid) const {
  EventCutResult result;
  result.particlePass.resize(pid.size(), false);
  result.particleDaughterPass.resize(pid.size(), false);
  result.MaxPhotonEnergyPass.resize(pid.size(), false);
  result.MotherMass.resize(pid.size(), -999);

  bool allCutsPassed = true;
  float MaxEphotonEnergy = 0.0f;
  float MaxPhotonEnergyIndex = 0;

  for (const auto& [name, cut] : fParticleCuts) {
    int count = 0;
    for (size_t i = 0; i < pid.size(); ++i) {
      const float p2 = px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i];
      if (p2 < 1e-4f) continue;

      if (pid[i] != cut.pid || charge[i] != cut.charge || REC_Track_pass_fid[i] != 1) continue;
      if (!IsInRange(chi2pid[i], cut.minChi2PID, cut.maxChi2PID)) continue;

      const float momentum = std::sqrt(p2);
      const float theta = std::atan2(std::sqrt(px[i] * px[i] + py[i] * py[i]), pz[i]);
      float phi = std::atan2(py[i], px[i]);
      if (phi < 0) phi += 2 * M_PI;
      int statusAbs = std::abs(status[i]);
      // cut based on momentum and status of th detector of the particle
      bool momentumFTCut = IsInRange(momentum, cut.minFTMomentum, cut.maxFTMomentum) && (statusAbs >= 1000 && statusAbs < 2000);
      bool momentumFDCut = IsInRange(momentum, cut.minFDMomentum, cut.maxFDMomentum) && (statusAbs >= 2000 && statusAbs < 3000);
      bool momentumCDCut = IsInRange(momentum, cut.minCDMomentum, cut.maxCDMomentum) && (statusAbs >= 4000 && statusAbs < 5000);
      bool momentumCut = momentumFTCut || momentumFDCut || momentumCDCut;

      bool betaCut = IsInRange(beta[i], cut.minBeta, cut.maxBeta);
      bool thetaCut = IsInRange(theta, cut.minTheta, cut.maxTheta);
      bool phiCut = IsInRange(phi, cut.minPhi, cut.maxPhi);
      bool vzCut = IsInRange(vz[i], cut.minVz, cut.maxVz);
      if (momentumCut && betaCut && thetaCut && phiCut && vzCut) {
        result.particlePass[i] = true;
        if (pid[i] == 22 && momentum > MaxEphotonEnergy) {
            MaxEphotonEnergy = momentum;
            result.MaxPhotonEnergyPass[MaxPhotonEnergyIndex] = false;
            result.MaxPhotonEnergyPass[i] = true;
            MaxPhotonEnergyIndex = i;
        }
        ++count;
      }
    }


    if (!IsInRange(count, cut.minCount, cut.maxCount)) {
      allCutsPassed = false;
    }
  }

  if (fAcceptEverything) {
    allCutsPassed = true;
  }

  if (fCutTwoBodyMotherDecay) {
    for (const auto& [name, cut] : fTwoBodyMotherCuts) {
      for (size_t i = 0; i < pid.size(); ++i) {
        if (pid[i] != cut.pidDaug1) continue;
        for (size_t j = i + 1; j < pid.size(); ++j) {
          if (pid[j] != cut.pidDaug2) continue;

          float minMass = cut.expectedMotherMass - cut.massSigma * cut.nSigmaMass;
          float maxMass = cut.expectedMotherMass + cut.massSigma * cut.nSigmaMass;

          float massDaug1 = getParticleMass( pid[i]);
          float massDaug2 = getParticleMass( pid[j]);

          float E1 = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i] + massDaug1 * massDaug1);
          float E2 = std::sqrt(px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j] + massDaug2 * massDaug2);
          float px_sum = px[i] + px[j];
          float py_sum = py[i] + py[j];
          float pz_sum = pz[i] + pz[j];
          float E_sum = E1 + E2;

          float invMass2 = E_sum * E_sum - (px_sum * px_sum + py_sum * py_sum + pz_sum * pz_sum);
          float invMass = (invMass2 > 0) ? std::sqrt(invMass2) : 0;
          // Store invariant mass regardless
          result.MotherMass[i] = invMass;
          result.MotherMass[j] = invMass;
          // Flag if passes the mass window
          if (invMass >= minMass && invMass <= maxMass) {
            result.particleDaughterPass[i] = true;
            result.particleDaughterPass[j] = true;
          }
        }
      }
    }
  }

  result.eventPass = allCutsPassed;
  
  return result;
}

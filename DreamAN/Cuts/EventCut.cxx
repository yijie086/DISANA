#include "EventCut.h"

#include <cmath>
#include <iostream>

EventCut::EventCut() = default;
EventCut::~EventCut() = default;

void EventCut::AddParticleCut(const std::string& name, const ParticleCut& userCut) {
  ParticleCut cut = userCut;

  if (name == "proton") {
     cut.pid = 2212;
     cut.charge = 1;
     /*cut.minMomentum = 0.0f;
     cut.maxMomentum = 500.0f;
      cut.minTheta = 0.0f;
      cut.maxTheta = 100.7f;
      cut.minVz = -100.0f;
      cut.maxVz = 100.0f;
      cut.minChi2PID = 0.0f;
      cut.maxChi2PID = 1000000.0f;
      cut.minCount = 1;
      cut.maxCount = 100;*/
    }
  else if (name == "electron") {
    cut.pid = 11;
    cut.charge = -1;
    /*cut.minMomentum = 0.0f;
      cut.maxMomentum = 500.0f;
      cut.minTheta = 0.0f;
      cut.maxTheta = 100.7f;
      cut.minVz = -100.0f;
      cut.maxVz = 100.0f;
      cut.minChi2PID = 0.0f;
      cut.maxChi2PID = 1000000.0f;
      cut.minCount = 1;
      cut.maxCount = 100;*/
  }else if (name == "photon") {
     cut.pid = 22;
     /*cut.charge = 0;
     cut.minMomentum = 0.0f;
     cut.maxMomentum = 500.0f;
     cut.minTheta = 0.0f;
     cut.maxTheta = 100.7f; 
     cut.minChi2PID = 0.0f;
      cut.maxChi2PID = 1000000.0f;    
      cut.minCount = 1;
      cut.maxCount = 999;*/
  }
  fParticleCuts[name] = cut;
}


const ParticleCut* EventCut::GetParticleCut(const std::string& name) const {
  auto it = fParticleCuts.find(name);
  return (it != fParticleCuts.end()) ? &it->second : nullptr;
}

EventCut* EventCut::ProtonCuts() {
  EventCut* cuts = new EventCut();
  ParticleCut proton;
  cuts->AddParticleCut("proton", proton);
  return cuts;
}

EventCut* EventCut::ElectronCuts() {
  EventCut* cuts = new EventCut();
  ParticleCut electron;
  cuts->AddParticleCut("electron", electron);
  return cuts;
}

EventCut* EventCut::PhotonCuts() {
  EventCut* cuts = new EventCut();
  ParticleCut photon;
  cuts->AddParticleCut("photon", photon);
  return cuts;
}

bool EventCut::operator()(const std::vector<int>& pid,
                          const std::vector<float>& px, const std::vector<float>& py, const std::vector<float>& pz,
                          const std::vector<float>& vx, const std::vector<float>& vy, const std::vector<float>& vz,
                          const std::vector<float>& vt,
                          const std::vector<short>& charge,
                          const std::vector<float>& beta,
                          const std::vector<float>& chi2pid,
                          const std::vector<short>& status,
                          const std::vector<int>& REC_Track_pass_fid) const {
  for (const auto& [name, cut] : fParticleCuts) {
    int pidCount = 0;

    for (size_t i = 0; i < pid.size(); ++i) {
      const float px_i = px[i], py_i = py[i], pz_i = pz[i];
      const float px2py2 = px_i * px_i + py_i * py_i;
      const float p2 = px2py2 + pz_i * pz_i;

      if (p2 < 1e-4f) continue;

      if (pid[i] != cut.pid || charge[i] != cut.charge || REC_Track_pass_fid[i] != 1 ||
          !IsInRange(chi2pid[i], cut.minChi2PID, cut.maxChi2PID)) {
        continue;
      }

      const float momentum = std::sqrt(p2);
      const float theta = std::atan2(std::sqrt(px2py2), pz_i);
      float phi = std::atan2(py_i, px_i);
      if (phi < 0) phi += 2 * M_PI;

      if (!(IsInRange(momentum, cut.minMomentum, cut.maxMomentum) &&
            IsInRange(theta, cut.minTheta, cut.maxTheta) &&
            IsInRange(phi, cut.minPhi, cut.maxPhi) &&
            IsInRange(vz[i], cut.minVz, cut.maxVz))) {
        continue;
      }

      ++pidCount;
    }

    if (!IsInRange(pidCount, cut.minCount, cut.maxCount)) {
      return false;
    }
  }

  return true;
}

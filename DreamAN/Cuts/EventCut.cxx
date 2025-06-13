#include "EventCut.h"

#include <cmath>
#include <iomanip>
#include <iostream>

// Constructor
EventCut::EventCut() = default;

// Destructor
EventCut::~EventCut() = default;

// Set charge range
void EventCut::SetChargeCut(int Charge) { fCharge = Charge; }

// Set momentum range
void EventCut::SetMomentumCut(float minMomentum, float maxMomentum) {
  fMinMomentum = minMomentum;
  fMaxMomentum = maxMomentum;
}

void EventCut::SetthetaCut(float minTheta, float maxTheta) {
  fMinTheta = minTheta;
  fMaxTheta = maxTheta;
}

void EventCut::SetPhiCut(float minPhi, float maxPhi) {
  fMinPhi = minPhi;
  fMaxPhi = maxPhi;
}

// Set vertex position range
void EventCut::SetVzCut(float minVz, float maxVz) {
  fMinVz = minVz;
  fMaxVz = maxVz;
}

// Set Chi2PID range
void EventCut::SetChi2PIDCut(float minChi2PID, float maxChi2PID) {
  fMinChi2PID = minChi2PID;
  fMaxChi2PID = maxChi2PID;
}

// Set the count limit for a specific PID
void EventCut::SetPIDCountCut(int selectedPID, int minCount, int maxCount) {
  fSelectedPID = selectedPID;
  fMinPIDCount = minCount;
  fMaxPIDCount = maxCount;
}

EventCut* EventCut::ProtonCuts() {
  EventCut* cuts = new EventCut();
  cuts->SetChargeCut(1);
  cuts->SetPIDCountCut(2212, 1, 10); // at this place you may take just one electron which is first in the list
  return cuts;
}

EventCut* EventCut::ElectronCuts() {
  EventCut* cuts = new EventCut();
  cuts->SetChargeCut(-1);
  cuts->SetPIDCountCut(11, 1, 99); // at this place you may choose all protons
  return cuts;
}

EventCut* EventCut::PhotonCuts() {
  EventCut* cuts = new EventCut();
  cuts->SetChargeCut(0); 
  cuts->SetPIDCountCut(22, 1, 999); // at this place you may choose all protons
  return cuts;
}

// Filtering logic with Fiducial cuts
bool EventCut::operator()(const std::vector<int>& pid, const std::vector<float>& px, const std::vector<float>& py, const std::vector<float>& pz, const std::vector<float>& vx,
                          const std::vector<float>& vy, const std::vector<float>& vz, const std::vector<float>& vt, const std::vector<int>& charge, const std::vector<float>& beta,
                          const std::vector<float>& chi2pid, const std::vector<int>& status, const std::vector<int>& REC_Track_pass_fid) const {
  int pidCount = 0;
  for (size_t i = 0; i < pid.size(); ++i) {
    // Skip trivial momentum
    const float px_i = px[i], py_i = py[i], pz_i = pz[i];
    const float px2py2 = px_i * px_i + py_i * py_i;
    const float p2 = px2py2 + pz_i * pz_i;

    if (p2 < 1e-4f) continue; // Equivalent to momentum < 0.01

    if (pid[i] != fSelectedPID || static_cast<int8_t>(charge[i]) != fCharge || REC_Track_pass_fid[i] != 1 || !IsInRange(chi2pid[i], fMinChi2PID, fMaxChi2PID))
      continue;

    const float momentum = std::sqrt(p2);
    const float theta = std::atan2(std::sqrt(px2py2), pz_i);
    float phi = std::atan2(py_i, px_i);
    if (phi < 0) phi += 2 * M_PI;

    if (!(IsInRange(momentum, fMinMomentum, fMaxMomentum) &&
          IsInRange(theta, fMinTheta, fMaxTheta) &&
          IsInRange(phi, fMinPhi, fMaxPhi) &&
          IsInRange(vz[i], fMinVz, fMaxVz))) {
      return false;  // Reject immediately on any kinematic failure
    }

    ++pidCount;
  }

  return IsInRange(pidCount, fMinPIDCount, fMaxPIDCount);
}
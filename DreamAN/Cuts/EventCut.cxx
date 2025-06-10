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
  cuts->SetPIDCountCut(2212, 1, 1);
  return cuts;
}

EventCut* EventCut::ElectronCuts() {
  EventCut* cuts = new EventCut();
  cuts->SetChargeCut(-1);
  cuts->SetPIDCountCut(11, 1, 1);
  return cuts;
}

EventCut* EventCut::PhotonCuts() {
  EventCut* cuts = new EventCut();
  cuts->SetChargeCut(0);
  cuts->SetPIDCountCut(22, 1, 1);
  return cuts;
}

// Filtering logic with Fiducial cuts
bool EventCut::operator()(const std::vector<int>& pid, const std::vector<float>& px, const std::vector<float>& py, const std::vector<float>& pz, const std::vector<float>& vx,
                          const std::vector<float>& vy, const std::vector<float>& vz, const std::vector<float>& vt, const std::vector<int>& charge, const std::vector<float>& beta,
                          const std::vector<float>& chi2pid, const std::vector<int>& status, const std::vector<int>& REC_Traj_pass,
                          const std::vector<int>& REC_Calorimeter_pass) const {
  int pidCount = 0;      // Count the number of target PIDs
  bool selected = true;  // Whether the target PID is selected

  for (size_t i = 0; i < pid.size(); ++i) {
    // Count the number of target PIDs
    if (pid[i] == fSelectedPID && std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]) > 0.01 && (static_cast<int8_t>(charge[i]) == fCharge) && REC_Traj_pass[i] == 1 &&
        REC_Calorimeter_pass[i] == 1 && IsInRange(chi2pid[i], fMinChi2PID, fMaxChi2PID)) {
      pidCount++;
      float momentum = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
      float theta = std::atan2(std::sqrt(px[i] * px[i] + py[i] * py[i]), pz[i]);
      float phi = std::atan2(py[i], px[i]);
      if (phi < 0) {
        phi += 2 * M_PI;
      }

      bool passMomentum = IsInRange(momentum, fMinMomentum, fMaxMomentum);
      bool passTheta = IsInRange(theta, fMinTheta, fMaxTheta);
      bool passPhi = IsInRange(phi, fMinPhi, fMaxPhi);
      bool passVz = IsInRange(vz[i], fMinVz, fMaxVz);

      if (!(passMomentum && passTheta && passPhi && passVz)) {
        selected = false;
        // std::cout << "Event rejected " << std::endl;
        //  std::cout << "edgeOk" << edgeOk << "EdgeCut: " << edgeCut << " EdgeVal: " << edgeVal <<"region number: "<<region<< std::endl;
      }
    }
  }
  // Check if the number of target PIDs is within the range
  if (IsInRange(pidCount, fMinPIDCount, fMaxPIDCount) && selected) {
    // std::cout << "Event selected!" << std::endl;
    return true;
  }
  return false;
}
#include "EventCut.h"
#include <iostream>
#include <cmath>
#include <iomanip>

// Constructor
EventCut::EventCut() = default;

// Destructor
EventCut::~EventCut() = default;

// Set charge range
void EventCut::SetChargeCut(int Charge) {
    fCharge = Charge;
}

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

// Filtering logic
bool EventCut::operator()(const std::vector<int>& pid,
                          const std::vector<float>& px,
                          const std::vector<float>& py,
                          const std::vector<float>& pz,
                          const std::vector<float>& vx,
                          const std::vector<float>& vy,
                          const std::vector<float>& vz,
                          const std::vector<float>& vt,
                          const std::vector<int>& charge,
                          const std::vector<float>& beta,
                          const std::vector<float>& chi2pid,
                          const std::vector<int>& status,
                          const std::vector<int>& REC_Traj_pass
                          ) const {
    int pidCount = 0; // Count the number of target PIDs
    bool selected = true; // Whether the target PID is selected
    
    for (size_t i = 0; i < pid.size(); ++i) {
        // Count the number of target PIDs
        if (pid[i] == fSelectedPID && 
            std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]) > 0.01 &&
            (static_cast<int8_t>(charge[i]) == fCharge) &&
            REC_Traj_pass[i] == 1 &&
            IsInRange(chi2pid[i], fMinChi2PID, fMaxChi2PID)) {
            pidCount++;
            float momentum = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            float theta = std::atan2(std::sqrt(px[i] * px[i] + py[i] * py[i]), pz[i]);
            float phi = std::atan2(py[i], px[i]);
            if (phi < 0) {
                phi += 2 * M_PI;
            }
            if (!((static_cast<int8_t>(charge[i]) == fCharge) &&
                IsInRange(momentum, fMinMomentum, fMaxMomentum) &&
                IsInRange(theta, fMinTheta, fMaxTheta) &&
                IsInRange(phi, fMinPhi, fMaxPhi) &&
                IsInRange(vz[i], fMinVz, fMaxVz))) {
                selected = false;
                //std::cout << pid[i] << " " << actual_charge << " " << momentum << " " << vz[i] << " " << chi2pid[i] << std::endl;
            }
        }
    }
    
    // Check if the number of target PIDs is within the range
    if (IsInRange(pidCount, fMinPIDCount, fMaxPIDCount) && selected) {
        //std::cout << "Event selected!" << std::endl;
        return true;
    }
    // Debug output in table format
    /*
    if (pidCount > 0) {
    std::cout << "Debug Information:" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "| Index | PID    | Charge | Momentum  | Vz       | Chi2PID | Selected |";
    std::cout << "| charge| momentum|  vz   | chi2pid   |" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;

    for (size_t i = 0; i < pid.size(); ++i) {
        float momentum = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
        bool isSelected = (pid[i] == fSelectedPID);
        std::cout << "| " << std::setw(5) << i 
                  << " | " << std::setw(6) << pid[i] 
                  << " | " << std::setw(6) << charge[i]
                  << " | " << std::setw(9) << std::fixed << std::setprecision(3) << momentum 
                  << " | " << std::setw(8) << vz[i] 
                  << " | " << std::setw(7) << chi2pid[i] 
                  << " | " << std::setw(8) << (isSelected ? "Yes" : "No") 
                  << " |" << std::setw(9) << (charge[i] == fCharge)
                  << " |" << std::setw(10) << IsInRange(momentum, fMinMomentum, fMaxMomentum)
                  << " |" << std::setw(11) << IsInRange(vz[i], fMinVz, fMaxVz)
                  << " |" << std::setw(12) << IsInRange(chi2pid[i], fMinChi2PID, fMaxChi2PID)
                  << " |" << std::setw(10) << selected
                  << " |" << std::endl;
    }

    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "Total particles: " << pid.size() << std::endl;
    std::cout << "Selected PID: " << fSelectedPID << std::endl;
    std::cout << "PID Count: " << pidCount << std::endl;

    std::cout << "Event rejected!" << std::endl;
    }
    */
    
    return false;
}
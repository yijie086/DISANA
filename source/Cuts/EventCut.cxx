#include "EventCut.h"
#include <iostream>
#include <cmath>
#include <iomanip>

// 构造函数
EventCut::EventCut() = default;

// 析构函数
EventCut::~EventCut() = default;

// 设置电荷范围
void EventCut::SetChargeCut(int Charge, bool ChargeSelection) {
    fCharge = Charge;
    fChargeSelection = ChargeSelection;
}

// 设置动量范围
void EventCut::SetMomentumCut(float minMomentum, float maxMomentum) {
    fMinMomentum = minMomentum;
    fMaxMomentum = maxMomentum;
}

// 设置顶点位置范围
void EventCut::SetVzCut(float minVz, float maxVz) {
    fMinVz = minVz;
    fMaxVz = maxVz;
}

// 设置 Chi2PID 范围
void EventCut::SetChi2PIDCut(float minChi2PID, float maxChi2PID) {
    fMinChi2PID = minChi2PID;
    fMaxChi2PID = maxChi2PID;
}

// 设置特定 PID 的数量限制
void EventCut::SetPIDCountCut(int selectedPID, int minCount, int maxCount) {
    fSelectedPID = selectedPID;
    fMinPIDCount = minCount;
    fMaxPIDCount = maxCount;
}

// 筛选逻辑
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
                          const std::vector<int>& status) const {
    int pidCount = 0; // 统计目标 PID 的数量
    bool selected = true; // 是否选择了目标 PID
    
    for (size_t i = 0; i < pid.size(); ++i) {
        // 统计目标 PID 的数量
        if (pid[i] == fSelectedPID && 
            std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i])>0.01 &&
            (static_cast<int8_t>(charge[i])==fCharge)) {
            pidCount++;
            float momentum = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (!((static_cast<int8_t>(charge[i])==fCharge) &&
                IsInRange(momentum, fMinMomentum, fMaxMomentum) &&
                IsInRange(vz[i], fMinVz, fMaxVz) &&
                IsInRange(chi2pid[i], fMinChi2PID, fMaxChi2PID))) {
                selected = false;
                //std::cout <<pid[i] << " " << actual_charge << " " << momentum << " " << vz[i] << " " << chi2pid[i] << std::endl;
            }
        }
    }
    
    // 检查目标 PID 的数量是否在范围内
    if (IsInRange(pidCount, fMinPIDCount, fMaxPIDCount) && selected) {
        //std::cout << "Event selected!" << std::endl;
        return true;
    }
    // Debug 输出为表的形式
    /*
    if (pidCount>0){
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
                  << " |" << std::setw(9) << (charge[i]==fCharge)
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
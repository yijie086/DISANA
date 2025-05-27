#ifndef EVENTCUT_H_
#define EVENTCUT_H_

#include <vector>
#include <string>
#include <tuple>

class EventCut {
public:
    EventCut();
    virtual ~EventCut();

    // Set filtering conditions
    void SetChargeCut(int Charge, bool ChargeSelection);
    void SetMomentumCut(float minMomentum, float maxMomentum);
    void SetVzCut(float minVz, float maxVz);
    void SetChi2PIDCut(float minChi2PID, float maxChi2PID);
    void SetPIDCountCut(int selectedPID, int minCount, int maxCount);

    // Filter function for RDataFrame
    bool operator()(const std::vector<int>& pid,
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
                    ) const;

private:
    // Filtering condition ranges
    int fCharge = -1; bool fChargeSelection = false;
    float fMinMomentum = -99999, fMaxMomentum = 99999;
    float fMinVz = -999999, fMaxVz = 999999;
    float fMinChi2PID = -999999, fMaxChi2PID = 999999;

    int fSelectedPID = 0; // Target PID
    int fMinPIDCount = 0, fMaxPIDCount = 9999;

    // Utility function: check if a value is within a range
    template <typename T>
    bool IsInRange(T value, T min, T max) const {
        return value >= min && value <= max;
    }
};

#endif // EVENTCUT_H_
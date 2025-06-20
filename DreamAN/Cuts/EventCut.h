#ifndef EVENTCUT_H_
#define EVENTCUT_H_

#include <string>
#include <map>
#include <vector>
#include <cfloat>
#include <cmath>

struct ParticleCut {
  int charge = 0;
  int pid = 0;
  int minCount = 1;
  int maxCount = 999;
  float minMomentum = 0;
  float maxMomentum = 20;
  float minTheta = -999;
  float maxTheta = M_PI;
  float minPhi = 0;
  float maxPhi = 2*M_PI;
  float minVz = -999;
  float maxVz = 999;
  float minChi2PID = -9999;
  float maxChi2PID = 999999;
};

class EventCut {
 public:
  EventCut();
  virtual ~EventCut();

  // Add particle-specific cuts; default values auto-filled for known names
  void AddParticleCut(const std::string& name, const ParticleCut& cut);
  const ParticleCut* GetParticleCut(const std::string& name) const;

  // Predefined sets
  static EventCut* ProtonCuts();
  static EventCut* ElectronCuts();
  static EventCut* PhotonCuts();

  // Apply cuts to an event
  bool operator()(const std::vector<int>& pid,
                  const std::vector<float>& px, const std::vector<float>& py, const std::vector<float>& pz,
                  const std::vector<float>& vx, const std::vector<float>& vy, const std::vector<float>& vz,
                  const std::vector<float>& vt,
                  const std::vector<short>& charge,
                  const std::vector<float>& beta,
                  const std::vector<float>& chi2pid,
                  const std::vector<short>& status,
                  const std::vector<int>& REC_Track_pass_fid) const;

 private:
  std::map<std::string, ParticleCut> fParticleCuts;

  template <typename T>
  bool IsInRange(T value, T min, T max) const {
    return value >= min && value <= max;
  }
};

#endif  // EVENTCUT_H_

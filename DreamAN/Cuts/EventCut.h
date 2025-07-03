
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
  /// detector dependent momentum cuts
  float minCDMomentum = 0;
  float minFDMomentum = 0;
  float minFTMomentum = 0;
  float maxCDMomentum = 999;
  float maxFDMomentum = 999;
  float maxFTMomentum = 999;
  //
  float minBeta = -999;
  float maxBeta = 999;
  float minTheta = -999;
  float maxTheta = M_PI;
  float minPhi = 0;
  float maxPhi = 2*M_PI;
  float minVz = -999;
  float maxVz = 999;
  float minChi2PID = -999999;
  float maxChi2PID = 999999;
};

struct TwoBodyMotherCut {
  int charge = 0;
  int pidDaug1 = 0;
  int pidDaug2 = 0;
  float expectedMotherMass = -999.0f;
  float massSigma = 0.005f;
  float nSigmaMass = 3.0;
};

struct EventCutResult {
  bool eventPass = false;
  std::vector<bool> particlePass;
  std::vector<bool> particleDaughterPass;
  std::vector<bool> MaxPhotonEnergyPass;
  std::vector<float> MotherMass; // corresponding
};

class EventCut {
 public:
  EventCut();
  virtual ~EventCut();
  void SetDoCutMotherInvMass(bool doCut) { fCutTwoBodyMotherDecay = doCut; }
  void AddParticleCut(const std::string& name, const ParticleCut& cut);
  void AddParticleMotherCut(const std::string& name, const TwoBodyMotherCut& cut);

  const ParticleCut* GetParticleCut(const std::string& name) const;

  static EventCut* ProtonCuts();
  static EventCut* ElectronCuts();
  static EventCut* PhotonCuts();

  EventCutResult operator()(const std::vector<int>& pid,
                            const std::vector<float>& px, const std::vector<float>& py, const std::vector<float>& pz,
                            const std::vector<float>& vx, const std::vector<float>& vy, const std::vector<float>& vz,
                            const std::vector<float>& vt,
                            const std::vector<short>& charge,
                            const std::vector<float>& beta,
                            const std::vector<float>& chi2pid,
                            const std::vector<short>& status,
                            const std::vector<int>& REC_Track_pass_fid) const;

 private:
  bool fCutTwoBodyMotherDecay = false;
  std::map<std::string, ParticleCut> fParticleCuts;
  std::map<std::string, TwoBodyMotherCut> fTwoBodyMotherCuts;

  template <typename T>
  bool IsInRange(T value, T min, T max) const {
    return value >= min && value <= max;
  }
};

#endif  // EVENTCUT_H_

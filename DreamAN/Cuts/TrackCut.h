#ifndef TRACKCUT_H_
#define TRACKCUT_H_

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

struct FiducialAxisCut {
  std::vector<std::pair<float, float>> excludedRanges;  // e.g., {{100, 120}, {240, 260}}
  std::set<float> excludedStrips;                       // e.g., {128.0, 256.0}
};

struct FiducialCut3D {
  FiducialAxisCut luCut;
  FiducialAxisCut lvCut;
  FiducialAxisCut lwCut;
};
class TrackCut {
 public:
  TrackCut();
  virtual ~TrackCut();
  // it is requited when copying the object itself
  TrackCut(const TrackCut& other);

  void SetSectorCut(int SSector, int selectpid, int selectdetector, bool selectSector);
  void SetPositionCut(float minX, float maxX, float minY, float maxY, float minZ, float maxZ);
  void SetDirectionCut(float minCX, float maxCX, float minCY, float maxCY, float minCZ, float maxCZ);
  void SetPathLengthCut(float minPath, float maxPath);
  void SetDCEdgeCut(float minEdge, float maxEdge);
  void SetECALEdgeCut(float minEdge, float maxEdge);
  void SetSectorCut_Bhawani(const std::vector<int>& sectors, int selectpid, int selectdetector, bool selectSector);
  void SetDoFiducialCut(bool isFiducial) { IsDoFiducial = isFiducial; }
  void SetDCEdgeCut(float minTheta, float maxTheta, int dcRegion, float edgeCut);
  void SetThetaBins(const std::vector<std::pair<float, float>>& thetaBins);
  void SetEdgeCuts(const std::vector<float>& edgeCutsPerRegion);
  float GetEdgeCut(int region) const;
  const std::vector<float>& GetEdgeCuts() const;

  bool operator()(const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int>& detector, const std::vector<int>& layer,
                  const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx, const std::vector<float>& cy,
                  const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge) const;

  // DC filter function
  std::function<std::vector<int>(const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int16_t>& detector, const std::vector<int16_t>& layer,
                                 const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx, const std::vector<float>& cy,
                                 const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge, const std::vector<int>& pid,
                                 const int& REC_Particle_num)>
  RECTrajPass() const;

  // Calorimeter filter function
  std::function<std::vector<int>(const std::vector<int16_t>&,  // index
                                 const std::vector<int16_t>&,  // pindex
                                 const std::vector<int16_t>&,  // detector
                                 const std::vector<int16_t>&,  // sector
                                 const std::vector<int16_t>&,  // layer
                                 const std::vector<float>&,    // energy
                                 const std::vector<float>&,    // time
                                 const std::vector<float>&,    // path
                                 const std::vector<float>&,    // chi2
                                 const std::vector<float>&,    // x
                                 const std::vector<float>&,    // y
                                 const std::vector<float>&,    // z
                                 const std::vector<float>&,    // hx
                                 const std::vector<float>&,    // hy
                                 const std::vector<float>&,    // hz
                                 const std::vector<float>&,    // lu
                                 const std::vector<float>&,    // lv
                                 const std::vector<float>&,    // lw
                                 const std::vector<float>&,    // du
                                 const std::vector<float>&,    // dv
                                 const std::vector<float>&,    // dw
                                 const std::vector<float>&,    // m2u
                                 const std::vector<float>&,    // m2v
                                 const std::vector<float>&,    // m2w
                                 const std::vector<float>&,    // m3u
                                 const std::vector<float>&,    // m3v
                                 const std::vector<float>&,    // m3w
                                 const std::vector<int>&,      // status
                                 const std::vector<int>&,      // pid
                                 const int& REC_Particle_num)>
  RECCalorimeterPass() const;

  // ===== New sector-specific fiducial setters =====
  void AddPCalFiducialStrip(int sector, const std::string& axis, float value);
  void AddPCalFiducialRange(int sector, const std::string& axis, float min, float max);
  void AddECinFiducialStrip(int sector, const std::string& axis, float value);
  void AddECinFiducialRange(int sector, const std::string& axis, float min, float max);
  void AddECoutFiducialStrip(int sector, const std::string& axis, float value);
  void AddECoutFiducialRange(int sector, const std::string& axis, float min, float max);

 private:
  bool fselectSector = false;
  bool IsDoFiducial = false;

  struct FiducialAxisCut {
    std::vector<std::pair<float, float>> excludedRanges;
    std::set<float> excludedStrips;
  };

  struct FiducialCut3D {
    FiducialAxisCut luCut;
    FiducialAxisCut lvCut;
    FiducialAxisCut lwCut;
  };

  float fSector = -1;
  int fselectPID = -1;
  int fselectdetector = 1000;
  std::set<int> fSectors;  // allowed sectors for DC
  float fMinX = -999999, fMaxX = 999999;
  float fMinY = -999999, fMaxY = 999999;
  float fMinZ = -999999, fMaxZ = 999999;
  float fMinCX = -999999, fMaxCX = 999999;
  float fMinCY = -999999, fMaxCY = 999999;
  float fMinCZ = -999999, fMaxCZ = 999999;
  float fMinPath = -999999, fMaxPath = 999999;
  float fDCMinEdge = -999999, fDCMaxEdge = 999999;
  float fECALMinEdge = -999999, fECALMaxEdge = 999999;
  std::vector<std::pair<float, float>> fThetaBins;                   // Still used for reference or binning
  std::vector<float> fEdgeCutsPerRegion{9999.0f, 9999.0f, 9999.0f};  // size = 3 → region 1, 2, 3

  /// ECin ECout and PCal Fiducial cuts
  std::map<int, FiducialCut3D> fFiducialCutsPCal;
  std::map<int, FiducialCut3D> fFiducialCutsECin;
  std::map<int, FiducialCut3D> fFiducialCutsECout;

  template <typename T>
  bool IsInRange(T value, T min, T max) const {
    return value >= min && value <= max;
  }
};

#endif  // TRACKCUT_H_
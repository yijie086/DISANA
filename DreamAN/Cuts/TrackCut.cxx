#include "TrackCut.h"

#include <cmath>
#include <iostream>
#include <map>


TrackCut::TrackCut() = default;
TrackCut::~TrackCut() = default;
TrackCut::TrackCut(const TrackCut& other) {
  this->fEdgeCutsPerRegion = other.fEdgeCutsPerRegion;
  this->fThetaBins = other.fThetaBins;
  this->fselectPID = other.fselectPID;
  this->fselectdetector = other.fselectdetector;
  this->fselectSector = other.fselectSector;
  this->fSector = other.fSector;
  this->fSectors = other.fSectors;
  this->IsDoFiducial = other.IsDoFiducial;

  // Also copy other necessary cuts if needed
  this->fMinX = other.fMinX;
  this->fMaxX = other.fMaxX;
  this->fMinY = other.fMinY;
  this->fMaxY = other.fMaxY;
  this->fMinZ = other.fMinZ;
  this->fMaxZ = other.fMaxZ;

  this->fMinCX = other.fMinCX;
  this->fMaxCX = other.fMaxCX;
  this->fMinCY = other.fMinCY;
  this->fMaxCY = other.fMaxCY;
  this->fMinCZ = other.fMinCZ;
  this->fMaxCZ = other.fMaxCZ;

  this->fMinPath = other.fMinPath;
  this->fMaxPath = other.fMaxPath;

  this->fDCMinEdge = other.fDCMinEdge;
  this->fDCMaxEdge = other.fDCMaxEdge;
  this->fECALMinEdge = other.fECALMinEdge;
  this->fECALMaxEdge = other.fECALMaxEdge;

  this->fFiducialCutsPCal = other.fFiducialCutsPCal;
  this->fFiducialCutsECin = other.fFiducialCutsECin;
  this->fFiducialCutsECout = other.fFiducialCutsECout;
}

void TrackCut::SetSectorCut(int SSector, int selectpid, int selectdetector, bool selectSector) {
  fSector = SSector;
  fselectSector = selectSector;
  fselectPID = selectpid;  
  fselectdetector = selectdetector;
}

void TrackCut::SetSectorCut_Bhawani(const std::vector<int>& sectors, int selectpid, int selectdetector, bool selectSector) {
  fSectors = std::set<int>(sectors.begin(), sectors.end());
  fselectSector = selectSector;
  fselectPID = selectpid;
  fselectdetector = selectdetector;
}

void TrackCut::SetPositionCut(float minX, float maxX, float minY, float maxY, float minZ, float maxZ) {
  fMinX = minX;
  fMaxX = maxX;
  fMinY = minY;
  fMaxY = maxY;
  fMinZ = minZ;
  fMaxZ = maxZ;
}

void TrackCut::SetDirectionCut(float minCX, float maxCX, float minCY, float maxCY, float minCZ, float maxCZ) {
  fMinCX = minCX;
  fMaxCX = maxCX;
  fMinCY = minCY;
  fMaxCY = maxCY;
  fMinCZ = minCZ;
  fMaxCZ = maxCZ;
}

void TrackCut::SetPathLengthCut(float minPath, float maxPath) {
  fMinPath = minPath;
  fMaxPath = maxPath;
}


void TrackCut::SetDCEdgeCut(float minEdge, float maxEdge) {
  fDCMinEdge = minEdge;
  fDCMaxEdge = maxEdge;
}

void TrackCut::SetECALEdgeCut(float minEdge, float maxEdge) {
  fECALMinEdge = minEdge;
  fECALMaxEdge = maxEdge;
}

void TrackCut::SetThetaBins(const std::vector<std::pair<float, float>>& thetaBins) { fThetaBins = thetaBins; }

void TrackCut::SetEdgeCuts(const std::vector<float>& edgeCutsPerRegion) {
  fEdgeCutsPerRegion = edgeCutsPerRegion;
  std::cout << "[Info] DC edge cuts are set for regions (R1, R2, R3): ";
  for (auto e : fEdgeCutsPerRegion) std::cout << e << " ";
  std::cout << std::endl;
}

float TrackCut::GetEdgeCut(int region) const {
  if (fEdgeCutsPerRegion.size() < 3) {
    throw std::runtime_error("Error: Edge cuts not initialized properly!");
  }
  return fEdgeCutsPerRegion[region - 1];
}

void TrackCut::AddPCalFiducialStrip(int sector, const std::string& axis, float value) {
  if (axis == "lu")
    fFiducialCutsPCal[sector].luCut.excludedStrips.insert(value);
  else if (axis == "lv")
    fFiducialCutsPCal[sector].lvCut.excludedStrips.insert(value);
  else if (axis == "lw")
    fFiducialCutsPCal[sector].lwCut.excludedStrips.insert(value);
  std::cout << "[Info] Added PCal fiducial strip: sector " << sector << ", axis " << axis << ", value " << value << std::endl;
}

void TrackCut::AddPCalFiducialRange(int sector, const std::string& axis, float min, float max) {
  if (axis == "lu")
    fFiducialCutsPCal[sector].luCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lv")
    fFiducialCutsPCal[sector].lvCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lw")
    fFiducialCutsPCal[sector].lwCut.excludedRanges.emplace_back(min, max);
  std::cout << "[Info] Added PCal fiducial range: sector " << sector << ", axis " << axis << ", range (" << min << ", " << max << ")" << std::endl;
}

void TrackCut::AddECinFiducialStrip(int sector, const std::string& axis, float value) {
  if (axis == "lu")
    fFiducialCutsECin[sector].luCut.excludedStrips.insert(value);
  else if (axis == "lv")
    fFiducialCutsECin[sector].lvCut.excludedStrips.insert(value);
  else if (axis == "lw")
    fFiducialCutsECin[sector].lwCut.excludedStrips.insert(value);
  std::cout << "[Info] Added ECin fiducial strip: sector " << sector << ", axis " << axis << ", value " << value << std::endl;
}

void TrackCut::AddECinFiducialRange(int sector, const std::string& axis, float min, float max) {
  if (axis == "lu")
    fFiducialCutsECin[sector].luCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lv")
    fFiducialCutsECin[sector].lvCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lw")
    fFiducialCutsECin[sector].lwCut.excludedRanges.emplace_back(min, max);
  std::cout << "[Info] Added ECin fiducial range: sector " << sector << ", axis " << axis << ", range (" << min << ", " << max << ")" << std::endl;
}

void TrackCut::AddECoutFiducialStrip(int sector, const std::string& axis, float value) {
  if (axis == "lu")
    fFiducialCutsECout[sector].luCut.excludedStrips.insert(value);
  else if (axis == "lv")
    fFiducialCutsECout[sector].lvCut.excludedStrips.insert(value);
  else if (axis == "lw")
    fFiducialCutsECout[sector].lwCut.excludedStrips.insert(value);
  std::cout << "[Info] Added ECout fiducial strip: sector " << sector << ", axis " << axis << ", value " << value << std::endl;
}

void TrackCut::AddECoutFiducialRange(int sector, const std::string& axis, float min, float max) {
  if (axis == "lu")
    fFiducialCutsECout[sector].luCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lv")
    fFiducialCutsECout[sector].lvCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lw")
    fFiducialCutsECout[sector].lwCut.excludedRanges.emplace_back(min, max);
  std::cout << "[Info] Added ECout fiducial range: sector " << sector << ", axis " << axis << ", range (" << min << ", " << max << ")" << std::endl;
}

const std::vector<float>& TrackCut::GetEdgeCuts() const { return fEdgeCutsPerRegion; }

bool TrackCut::operator()(const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int>& detector, const std::vector<int>& layer, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx,
                          const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge) const {
  return true;  //
}

/// copy the fiducical cuts here and should be used it from the EventFilte
std::function<std::vector<int>(const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int16_t>& detector, const std::vector<int16_t>& layer, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx,
                               const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge, const std::vector<int>& pid, const int& REC_Particle_num)>
TrackCut::RECTrajPass() const {
  return [this](const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int16_t>& detector, const std::vector<int16_t>& layer, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx,
                const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge, const std::vector<int>& pid, const int& REC_Particle_num) -> std::vector<int> {
    std::vector<int> pass_values(REC_Particle_num, 1);
    for (size_t i = 0; i < pindex.size(); ++i) {
      if (detector[i] == 6) {  // DC
        if (IsDoFiducial) {
          int region = 0;
          int absLayer = std::abs(layer[i]);
          if (absLayer == 6) {
            region = 1;
          } else if (absLayer == 18) {
            region = 2;
          } else if (absLayer == 36) {
            region = 3;
          }
          float edgeCut = this->GetEdgeCut(region);
          if (edge[i] <= edgeCut) {
            pass_values[pindex[i]] = 0;
            continue;
          }
        }
      }
    }

    return pass_values;
  };
}

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
TrackCut::RECCalorimeterPass() const {
  return [this](const std::vector<int16_t>& index, const std::vector<int16_t>& pindex, const std::vector<int16_t>& detector, const std::vector<int16_t>& sector, const std::vector<int16_t>& layer, const std::vector<float>& energy, const std::vector<float>& time, const std::vector<float>& path,
                const std::vector<float>& chi2, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& hx, const std::vector<float>& hy, const std::vector<float>& hz, const std::vector<float>& lu, const std::vector<float>& lv,
                const std::vector<float>& lw, const std::vector<float>& du, const std::vector<float>& dv, const std::vector<float>& dw, const std::vector<float>& m2u, const std::vector<float>& m2v, const std::vector<float>& m2w, const std::vector<float>& m3u, const std::vector<float>& m3v,
                const std::vector<float>& m3w, const std::vector<int>& status, const std::vector<int>& pid, const int& REC_Particle_num) -> std::vector<int> {
    // Initialize return_values with size REC_Particle_num and default value 9999.0
    std::vector<int> return_values(REC_Particle_num, 1);
    auto isExcluded = [](float value, const FiducialAxisCut& cut) -> bool {
      for (const auto& range : cut.excludedRanges) {
        if (value >= range.first && value <= range.second) return true;
      }
      return cut.excludedStrips.find(value) != cut.excludedStrips.end();
    };

    for (size_t i = 0; i < pindex.size(); ++i) {
      if (detector[i] == 7) {
        if(IsDoFiducial){
        const std::map<int, FiducialCut3D>* cutMap = nullptr;
        if (layer[i] == 1)
          cutMap = &fFiducialCutsPCal;
        else if (layer[i] == 4)
          cutMap = &fFiducialCutsECin;
        else if (layer[i] == 7)
          cutMap = &fFiducialCutsECout;
        if (cutMap) {
          auto it = cutMap->find(sector[i]);
          if (it != cutMap->end()) {
            const FiducialCut3D& cut = it->second;
            if (isExcluded(lu[i], cut.luCut) || isExcluded(lv[i], cut.lvCut) || isExcluded(lw[i], cut.lwCut)) {
              return_values[pindex[i]] = 0;
              continue;
            }
          }
        }
      }
    }
    }
        return return_values;
  };
}
#include "TrackCut.h"
#include <cmath>
#include <iostream>

// 构造函数
TrackCut::TrackCut() = default;

// 析构函数
TrackCut::~TrackCut() = default;

// 设置位置范围筛选条件

void TrackCut::SetSectorCut(int SSector, int selectpid, int selectdetector, bool selectSector) {
    fSector = SSector;
    fselectSector = selectSector;
    fselectPID = selectpid; // 设置选择的 PID
    fselectdetector = selectdetector;
}

void TrackCut::SetPositionCut(float minX, float maxX, float minY, float maxY, float minZ, float maxZ) {
    fMinX = minX; fMaxX = maxX;
    fMinY = minY; fMaxY = maxY;
    fMinZ = minZ; fMaxZ = maxZ;
}

// 设置方向余弦范围筛选条件
void TrackCut::SetDirectionCut(float minCX, float maxCX, float minCY, float maxCY, float minCZ, float maxCZ) {
    fMinCX = minCX; fMaxCX = maxCX;
    fMinCY = minCY; fMaxCY = maxCY;
    fMinCZ = minCZ; fMaxCZ = maxCZ;
}

// 设置路径长度范围筛选条件
void TrackCut::SetPathLengthCut(float minPath, float maxPath) {
    fMinPath = minPath; fMaxPath = maxPath;
}

// 设置边缘距离范围筛选条件
void TrackCut::SetDCEdgeCut(float minEdge, float maxEdge) {
    fDCMinEdge = minEdge; fDCMaxEdge = maxEdge;
}

void TrackCut::SetECALEdgeCut(float minEdge, float maxEdge) {
    fECALMinEdge = minEdge; fECALMaxEdge = maxEdge;
}

// 过滤逻辑
bool TrackCut::operator()(const std::vector<int16_t>& pindex,
                          const std::vector<int16_t>& index,
                          const std::vector<int>& detector,
                          const std::vector<int>& layer,
                          const std::vector<float>& x,
                          const std::vector<float>& y,
                          const std::vector<float>& z,
                          const std::vector<float>& cx,
                          const std::vector<float>& cy,
                          const std::vector<float>& cz,
                          const std::vector<float>& path,
                          const std::vector<float>& edge) const {
/*    for (size_t i = 0; i < pindex.size(); ++i) {
        if (IsInRange(x[i], fMinX, fMaxX) &&
            IsInRange(y[i], fMinY, fMaxY) &&
            IsInRange(z[i], fMinZ, fMaxZ) &&
            IsInRange(cx[i], fMinCX, fMaxCX) &&
            IsInRange(cy[i], fMinCY, fMaxCY) &&
            IsInRange(cz[i], fMinCZ, fMaxCZ) &&
            IsInRange(path[i], fMinPath, fMaxPath) &&
            IsInRange(edge[i], fMinEdge, fMaxEdge)) {
            return true; // 如果满足所有条件，则选择该轨迹
        }
    }
*/
    return true; // 如果没有轨迹满足条件，则拒绝事件
}

// RECTrajPass 实现
std::function<std::vector<int>(const std::vector<int16_t>& pindex,
                               const std::vector<int16_t>& index,
                               const std::vector<int16_t>& detector,
                               const std::vector<int16_t>& layer,
                               const std::vector<float>& x,
                               const std::vector<float>& y,
                               const std::vector<float>& z,
                               const std::vector<float>& cx,
                               const std::vector<float>& cy,
                               const std::vector<float>& cz,
                               const std::vector<float>& path,
                               const std::vector<float>& edge,
                               const std::vector<int>& pid,
                               const int& REC_Particle_num)> TrackCut::RECTrajPass() const {
    return [this](const std::vector<int16_t>& pindex,
                  const std::vector<int16_t>& index,
                  const std::vector<int16_t>& detector,
                  const std::vector<int16_t>& layer,
                  const std::vector<float>& x,
                  const std::vector<float>& y,
                  const std::vector<float>& z,
                  const std::vector<float>& cx,
                  const std::vector<float>& cy,
                  const std::vector<float>& cz,
                  const std::vector<float>& path,
                  const std::vector<float>& edge,
                  const std::vector<int>& pid,
                  const int& REC_Particle_num) -> std::vector<int> {
        // 初始化 pass_values，大小为 REC_Particle_num，默认所有粒子通过
        //std::cout << "TrackCut::RECTrajPass called with REC_Particle_num: " << REC_Particle_num << std::endl;
        std::vector<int> pass_values(REC_Particle_num, 1);

        // 遍历 RECTraj 的所有行
        for (size_t i = 0; i < pindex.size(); ++i) {
            //std::cout << pid[pindex[i]] << " " << detector[i] << " " << layer[i] << " " << edge[i] << std::endl;
            if (detector[i] == 6){//DC
                if (!IsInRange(edge[i], fDCMinEdge, fDCMaxEdge)) {
                    pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                }
                float temp_phi = std::atan2(y[i], x[i]);
                if (temp_phi < 0) {
                    temp_phi += 2 * M_PI; // 确保 phi 在 [0, 2π] 范围内
                }
                if (fselectSector && layer[i] == 6 && pid[pindex[i]]==fselectPID && detector[i]==fselectdetector) { // R1
                    if ((temp_phi>=0 && temp_phi<30*M_PI/180)||(temp_phi>=330*M_PI/180 && temp_phi<360*M_PI/180)){
                        if (fSector != 1) {
                            pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                        }
                        //std::cout << "PCAL sector 1 " << temp_phi*180/M_PI <<std::endl;
                    }
                    if (temp_phi>=30*M_PI/180 && temp_phi<90*M_PI/180) {
                        if (fSector != 2) {
                            pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                        }
                        //std::cout << "PCAL sector 2 " << temp_phi*180/M_PI<<std::endl;
                    }
                    if (temp_phi>=90*M_PI/180 && temp_phi<150*M_PI/180) {
                        if (fSector != 3) {
                            pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                        }
                        //std::cout << "PCAL sector 3 " << temp_phi*180/M_PI<<std::endl;
                    }
                    if (temp_phi>=150*M_PI/180 && temp_phi<210*M_PI/180) {
                        if (fSector != 4) {
                            pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                        }
                        //std::cout << "PCAL sector 4 " << temp_phi*180/M_PI<<std::endl;
                    }
                    if (temp_phi>=210*M_PI/180 && temp_phi<270*M_PI/180) {
                        if (fSector != 5) {
                            pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                        }
                        //std::cout << "PCAL sector 5 " << temp_phi*180/M_PI<<std::endl;
                    }
                    if (temp_phi>=270*M_PI/180 && temp_phi<330*M_PI/180) {
                        if (fSector != 6) {
                            pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                        }
                        //std::cout << "PCAL sector 6 " << temp_phi*180/M_PI<<std::endl;
                    }
                }
            }
        }

        return pass_values;
    };
}

std::function<std::vector<int>(const std::vector<int16_t>&,      // index
                                 const std::vector<int16_t>&,      // pindex
                                 const std::vector<int16_t>&,      // detector
                                 const std::vector<int16_t>&,      // sector
                                 const std::vector<int16_t>&,      // layer
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
                                 const std::vector<int>&, //pid
                                 const int& REC_Particle_num)> TrackCut::RECCalorimeterPass() const{
    return [this](const std::vector<int16_t>& index,
                  const std::vector<int16_t>& pindex,
                  const std::vector<int16_t>& detector,
                  const std::vector<int16_t>& sector,
                  const std::vector<int16_t>& layer,
                  const std::vector<float>& energy,
                  const std::vector<float>& time,
                  const std::vector<float>& path,
                  const std::vector<float>& chi2,
                  const std::vector<float>& x,
                  const std::vector<float>& y,
                  const std::vector<float>& z,
                  const std::vector<float>& hx,
                  const std::vector<float>& hy,
                  const std::vector<float>& hz,
                  const std::vector<float>& lu,
                  const std::vector<float>& lv,
                  const std::vector<float>& lw,
                  const std::vector<float>& du,
                  const std::vector<float>& dv,
                  const std::vector<float>& dw,
                  const std::vector<float>& m2u,
                  const std::vector<float>& m2v,
                  const std::vector<float>& m2w,
                  const std::vector<float>& m3u,
                  const std::vector<float>& m3v,
                  const std::vector<float>& m3w,
                  const std::vector<int>& status,
                  const std::vector<int>& pid,
                  const int& REC_Particle_num) -> std::vector<int> {
        // Initialize return_values with size REC_Particle_num and default value 9999.0
        std::vector<int> return_values(REC_Particle_num, 1);
        
        for (size_t i = 0; i < pindex.size(); ++i) {
            if (detector[i] == 7) {
                if (fselectSector && layer[i] == 1 && pid[pindex[i]]==fselectPID && detector[i]==fselectdetector){
                    if (sector[i] == 1 && fSector != 1) {
                        return_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                    } else if (sector[i] == 2 && fSector != 2) {
                        return_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                    } else if (sector[i] == 3 && fSector != 3) {
                        return_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                    } else if (sector[i] == 4 && fSector != 4) {
                        return_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                    } else if (sector[i] == 5 && fSector != 5) {
                        return_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                    } else if (sector[i] == 6 && fSector != 6) {
                        return_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                    }
                    //std::cout << "sector: " << sector[i] << " layer: " << layer[i] << " pid: " << pid[pindex[i]] << " detector: " << detector[i] << std::endl;
                    //std::cout << "selected: " << return_values[pindex[i]] << std::endl;
                }
                
            }            
        }

        return return_values;
    };
}
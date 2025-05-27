#include "TrackCut.h"
#include <cmath>
#include <iostream>

// 构造函数
TrackCut::TrackCut() = default;

// 析构函数
TrackCut::~TrackCut() = default;

// 设置位置范围筛选条件
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
                  const int& REC_Particle_num) -> std::vector<int> {
        // 初始化 pass_values，大小为 REC_Particle_num，默认所有粒子通过
        //std::cout << "TrackCut::RECTrajPass called with REC_Particle_num: " << REC_Particle_num << std::endl;
        std::vector<int> pass_values(REC_Particle_num, 1);

        // 遍历 RECTraj 的所有行
        for (size_t i = 0; i < pindex.size(); ++i) {
            if (detector[i] == 7){//ECAL
                if (!IsInRange(edge[i], fECALMinEdge, fECALMaxEdge)) {
                    pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                }
                
            }
            if (detector[i] == 6){//DC
                if (!IsInRange(edge[i], fDCMinEdge, fDCMaxEdge)) {
                    pass_values[pindex[i]] = 0; // 如果不满足条件，则标记为 0
                }
            }
        }

        return pass_values;
    };
}

std::function<std::vector<float>(const std::vector<int16_t>& pindex,
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
                                 const int& REC_Particle_num)> TrackCut::RECTrajedge(int target_detector, int target_layer) const {
    return [this,target_detector, target_layer](
                  const std::vector<int16_t>& pindex,
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
                  const int& REC_Particle_num) -> std::vector<float> {
        // 初始化 pass_values，大小为 REC_Particle_num，默认所有粒子通过
        std::vector<float> edge_values(REC_Particle_num, 9999);
        for (size_t i = 0; i < pindex.size(); ++i) {
            if (detector[i] == target_detector && layer[i]==target_layer){
                edge_values[pindex[i]] = edge[i];
            }
        }

        return edge_values;
    };
}
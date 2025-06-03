#ifndef TRACKCUT_H_
#define TRACKCUT_H_

#include <vector>
#include <string>
#include <functional>

class TrackCut {
public:
    TrackCut();
    virtual ~TrackCut();

    // 设置筛选条件
    void SetSectorCut(int SSector, int selectpid, bool selectSector);
    void SetPositionCut(float minX, float maxX, float minY, float maxY, float minZ, float maxZ);
    void SetDirectionCut(float minCX, float maxCX, float minCY, float maxCY, float minCZ, float maxCZ);
    void SetPathLengthCut(float minPath, float maxPath);
    void SetDCEdgeCut(float minEdge, float maxEdge);
    void SetECALEdgeCut(float minEdge, float maxEdge);

    // 用于 RDataFrame 的过滤函数
    bool operator()(const std::vector<int16_t>& pindex,
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
                    const std::vector<float>& edge) const;

    // 用于生成 REC_Traj_pass 列
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
                                   const int& REC_Particle_num)> RECTrajPass() const;
    

private:
    // 筛选条件
    float fSector = -1; bool fselectSector = false; int fselectPID = -1;
    float fMinX = -999999, fMaxX = 999999;
    float fMinY = -999999, fMaxY = 999999;
    float fMinZ = -999999, fMaxZ = 999999;
    float fMinCX = -999999, fMaxCX = 999999;
    float fMinCY = -999999, fMaxCY = 999999;
    float fMinCZ = -999999, fMaxCZ = 999999;
    float fMinPath = -999999, fMaxPath = 999999;
    float fDCMinEdge = -999999, fDCMaxEdge = 999999;
    float fECALMinEdge = -999999, fECALMaxEdge = 999999;

    // 工具函数：检查值是否在范围内
    template <typename T>
    bool IsInRange(T value, T min, T max) const {
        return value >= min && value <= max;
    }
};

#endif // TRACKCUT_H_
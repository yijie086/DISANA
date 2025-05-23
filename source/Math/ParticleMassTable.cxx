#include "ParticleMassTable.h"
#include <stdexcept>

// 定义 getParticleMass 函数
double getParticleMass(int pid) {
    // 定义常见粒子的 pid 和质量（单位：GeV/c^2）
    static const std::unordered_map<int, double> particleMasses = {
        {11, 0.000511},    // 电子 (Electron)
        {-11, 0.000511},   // 正电子 (Positron)
        {2212, 0.938272},  // 质子 (Proton)
        {-2212, 0.938272}, // 反质子 (Anti-Proton)
        {211, 0.13957},    // 正π介子 (Pi+)
        {-211, 0.13957},   // 负π介子 (Pi-)
        {22, 0.0},         // 光子 (Photon)
        {13, 0.10566},     // μ子 (Muon)
        {-13, 0.10566},    // 反μ子 (Anti-Muon)
        {111, 0.13498},    // 中性π介子 (Pi0)
        {2112, 0.939565},  // 中子 (Neutron)
        {-2112, 0.939565}  // 反中子 (Anti-Neutron)
    };

    // 查找 pid 对应的质量
    auto it = particleMasses.find(pid);
    if (it != particleMasses.end()) {
        return it->second; // 返回质量
    } else {
        throw std::invalid_argument("Unknown pid: " + std::to_string(pid)); // 未知 pid 抛出异常
    }
}
#include "RECParticleKinematic.h"
#include <cmath>

// 返回计算 theta 的 lambda 函数
RECStoreType RECParticletheta() {
    return [](const std::vector<int>& pid,
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
              const std::vector<int>& status) -> std::vector<float> {
        
        std::vector<float> theta_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float theta = std::atan2(std::sqrt(px[i] * px[i] + py[i] * py[i]), pz[i]);
            theta_values.push_back(theta);
        }
        return theta_values;
    };
}

RECStoreType RECParticlephi() {
    return [](const std::vector<int>& pid,
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
              const std::vector<int>& status) -> std::vector<float> {
        
        std::vector<float> phi_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float phi = std::atan2(py[i], px[i]);
            if (phi < 0) {
                phi += 2 * M_PI;
            }
            phi_values.push_back(phi);
        }
        return phi_values;
    };
}

RECStoreType RECParticleP(){
    return [](const std::vector<int>& pid,
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
              const std::vector<int>& status) -> std::vector<float> {
        
        std::vector<float> p_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float p = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            p_values.push_back(p);
        }
        return p_values;
    };
}
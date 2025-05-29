#include "RECParticle.h"
#include <iostream>
#include <cmath>
#include <functional>

const std::vector<std::string>& RECParticle::All() {
    static const std::vector<std::string> names = {
        "REC_Particle_pid",
        "REC_Particle_px",
        "REC_Particle_py",
        "REC_Particle_pz",
        "REC_Particle_vx",
        "REC_Particle_vy",
        "REC_Particle_vz",
        "REC_Particle_vt",
        "REC_Particle_charge",
        "REC_Particle_beta",
        "REC_Particle_chi2pid",
        "REC_Particle_status"
    };
    return names;
}

const std::vector<std::string>& RECParticle::Extend() {
    static const std::vector<std::string> extended = [] {
        std::vector<std::string> base = RECParticle::All();  
        base.push_back("REC_Particle_phi");
        base.push_back("REC_Particle_theta");
        base.push_back("REC_Particle_p");
        return base;
    }();
    return extended;
}


std::function<std::vector<float>(const std::vector<float>&, const std::vector<int>&, const std::vector<int>&, const std::vector<float>&, const std::vector<int>&)> 
get_RECParticle_float_var(int target_pid, int target_charge) {
    return [target_pid, target_charge](const std::vector<float>& var, const std::vector<int>& pid, const std::vector<int>& charge, const std::vector<float>& p, const std::vector<int>& trackpass) {
        std::vector<float> out;
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && p[i] > 0.02 && trackpass[i] == 1) {
                out.push_back(var[i]);
            }
        }
        return out;
    };
}

std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&, const std::vector<int>&, const std::vector<float>&, const std::vector<int>&)> 
get_RECParticle_int_var(int target_pid, int target_charge) {
    return [target_pid, target_charge](const std::vector<int>& var, const std::vector<int>& pid, const std::vector<int>& charge, const std::vector<float>& p, const std::vector<int>& trackpass) {
        std::vector<int> out;
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && p[i] > 0.02 && trackpass[i] == 1) {
                out.push_back(var[i]);
            }
        }
        return out;
    };
}


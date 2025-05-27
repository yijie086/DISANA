#include "RECParticleKinematic.h"
#include <cmath>
#include "MathKinematicVariable.h"
#include "ParticleMassTable.h"

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

RECStoreType EventQ2(float E, int target_pid, int target_charge) {
    return [E,target_pid,target_charge](const std::vector<int>& pid,
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
        
        std::vector<float> Q2_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float Q2 = 0.0;
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && pz[i] > 0.02) {
                TLorentzVector p4_beam(0.0,0.0,E,E); // px, py, pz, E
                TLorentzVector p4_electron(px[i],py[i],pz[i],std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+getParticleMass(pid[i])*getParticleMass(pid[i]))); // px, py, pz, E
                Q2 = math_Q2(p4_beam, p4_electron);
            }
            Q2_values.push_back(Q2);
        }
        return Q2_values;
    };
}

RECStoreType EventxB(float E, int target_pid, int target_charge, float target_mass) {
    return [E,target_pid,target_charge,target_mass](const std::vector<int>& pid,
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
        
        std::vector<float> xB_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float xB = 0.0;
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && pz[i] > 0.02) {
                TLorentzVector p4_beam(0.0,0.0,E,E); // px, py, pz, E
                TLorentzVector p4_electron(px[i],py[i],pz[i],std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+getParticleMass(pid[i])*getParticleMass(pid[i]))); // px, py, pz, E
                xB = math_xB(p4_beam, p4_electron, target_mass);
            }
            xB_values.push_back(xB);
        }
        return xB_values;
    };
}

RECStoreType EventNu(float E, int target_pid, int target_charge) {
    return [E,target_pid,target_charge](const std::vector<int>& pid,
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
        
        std::vector<float> nu_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float nu = 0.0;
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && pz[i] > 0.02) {
                TLorentzVector p4_beam(0.0,0.0,E,E); // px, py, pz, E
                TLorentzVector p4_electron(px[i],py[i],pz[i],std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+getParticleMass(pid[i])*getParticleMass(pid[i]))); // px, py, pz, E
                nu = math_Nu(p4_beam, p4_electron);
            }
            nu_values.push_back(nu);
        }
        return nu_values;
    };
}

RECStoreType EventW(float E, int target_pid, int target_charge, float target_mass) {
    return [E,target_pid,target_charge,target_mass](const std::vector<int>& pid,
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
        
        std::vector<float> W_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float W = 0.0;
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && pz[i] > 0.02) {
                TLorentzVector p4_beam(0.0,0.0,E,E); // px, py, pz, E
                TLorentzVector p4_electron(px[i],py[i],pz[i],std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+getParticleMass(pid[i])*getParticleMass(pid[i]))); // px, py, pz, E
                W = math_W(p4_beam, p4_electron, target_mass);
            }
            W_values.push_back(W);
        }
        return W_values;
    };
}

RECStoreType Eventmt(float E, int target_pid, int target_charge, float target_mass) {
    return [E,target_pid,target_charge,target_mass](const std::vector<int>& pid,
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
        
        std::vector<float> t_values;
        for (size_t i = 0; i < pid.size(); ++i) {
            float t = 0.0;
            if (pid[i] == target_pid && static_cast<int8_t>(charge[i])==target_charge && pz[i] > 0.02) {
                TLorentzVector p4_recoil(px[i],py[i],pz[i],std::sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+getParticleMass(pid[i])*getParticleMass(pid[i]))); // px, py, pz, E
                t = -math_t(p4_recoil, target_mass);
            }
            t_values.push_back(t);
        }
        return t_values;
    };
}
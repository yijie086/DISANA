#ifndef RECPARTICLEKINEMATIC_H
#define RECPARTICLEKINEMATIC_H

#include <vector>
#include <functional>

// 定义返回类型
using RECStoreType = std::function<std::vector<float>(
    const std::vector<int>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<int>&,
    const std::vector<float>&,
    const std::vector<float>&,
    const std::vector<int>&
)>;

// 声明函数：返回一个 lambda
RECStoreType RECParticletheta();
RECStoreType RECParticlephi();
RECStoreType RECParticleP();

RECStoreType EventQ2(float E, int target_pid, int target_charge);
RECStoreType EventxB(float E, int target_pid, int target_charge, float target_mass);
RECStoreType EventNu(float E, int target_pid, int target_charge);
RECStoreType EventW(float E, int target_pid, int target_charge, float target_mass);
RECStoreType Eventmt(float E, int target_pid, int target_charge, float target_mass);


#endif

#ifndef ELECTRONCUT_H
#define ELECTRONCUT_H

#include <vector>

bool ElectronCut(
    const std::vector<int>& pid,
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
    const std::vector<int>& status
);

#endif
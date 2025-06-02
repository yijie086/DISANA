#include "RECTraj.h"
#include <iostream>
#include <cmath>
#include <functional>

const std::vector<std::string>& RECTraj::All() {
    static const std::vector<std::string> names = {
        "REC_Traj_pindex",
        "REC_Traj_index",
        "REC_Traj_detector",
        "REC_Traj_layer",
        "REC_Traj_x",
        "REC_Traj_y",
        "REC_Traj_z",
        "REC_Traj_cx",
        "REC_Traj_cy",
        "REC_Traj_cz",
        "REC_Traj_path",
        "REC_Traj_edge"
    };
    return names;
}

const std::vector<std::string>& RECTraj::Extend() {
    static const std::vector<std::string> extended = [] {
        std::vector<std::string> base = RECTraj::All();
        // 如果需要扩展更多变量，可以在这里添加
        return base;
    }();
    return extended;
}

// RECTrajX 函数实现
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
                                 const int& REC_Particle_num)> RECTrajXYZ(int target_detector, int target_layer, int xyz) {
    return [target_detector, target_layer, xyz](
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
        // 初始化 pass_values，大小为 REC_Particle_num，默认所有粒子9999为没有hit
        std::vector<float> return_values(REC_Particle_num, 9999);
        if (xyz == 1){
            for (size_t i = 0; i < pindex.size(); ++i) {
                if (detector[i] == target_detector && layer[i] == target_layer) {
                    return_values[pindex[i]] = x[i];
                }
            }
        } else if (xyz == 2) {
            for (size_t i = 0; i < pindex.size(); ++i) {
                if (detector[i] == target_detector && layer[i] == target_layer) {
                    return_values[pindex[i]] = y[i];
                }
            }
        } else if (xyz == 3) {
            for (size_t i = 0; i < pindex.size(); ++i) {
                if (detector[i] == target_detector && layer[i] == target_layer) {
                    return_values[pindex[i]] = z[i];
                }
            }
        }
        return return_values;
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
                                 const int& REC_Particle_num)> RECTrajedge(int target_detector, int target_layer) {
    return [target_detector, target_layer](
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


std::function<std::vector<float>(const std::vector<float>&,
                                 const std::vector<int16_t>&,
                                 const std::vector<int16_t>&,
                                 const std::vector<int16_t>&,
                                 const std::vector<int>&,
                                 const std::vector<int>&,
                                 const std::vector<float>&,
                                 const std::vector<int>&)>
get_RECTraj_float_var(int target_detector, int target_layer, int target_pid, int target_charge) {
    return [target_detector, target_layer, target_pid, target_charge](
               const std::vector<float>& var,
               const std::vector<int16_t>& detector,
               const std::vector<int16_t>& layer,
               const std::vector<int16_t>& pindex,
               const std::vector<int>& pid,
               const std::vector<int>& charge,
               const std::vector<float>& p,
               const std::vector<int>& trackpass) -> std::vector<float> {
        std::vector<float> out;
        for (size_t i = 0; i < pid.size(); ++i) {
            for (size_t j = 0; j < detector.size(); ++j) {
                if (detector[j] == target_detector && layer[j] == target_layer &&
                    pindex[j] == i &&
                    pid[i] == target_pid && static_cast<int8_t>(charge[i]) == target_charge &&
                    p[i] > 0.02 && trackpass[i] == 1) {
                    out.push_back(var[i]);
                }
            }
        }
        return out;
    };
}

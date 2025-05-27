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

std::function<std::vector<float>(const std::vector<float>&, const std::vector<int>&, const std::vector<int>&)> 
get_RECTraj_float_var(int target_detector, int target_layer) {
    return [target_detector, target_layer](const std::vector<float>& var, const std::vector<int>& detector, const std::vector<int>& layer) {
        std::vector<float> out;
        for (size_t i = 0; i < detector.size(); ++i) {
            if (detector[i] == target_detector && layer[i] == target_layer) {
                out.push_back(var[i]);
            }
        }
        return out;
    };
}

std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&, const std::vector<int>&)> 
get_RECTraj_int_var(int target_detector, int target_layer) {
    return [target_detector, target_layer](const std::vector<int>& var, const std::vector<int>& detector, const std::vector<int>& layer) {
        std::vector<int> out;
        for (size_t i = 0; i < detector.size(); ++i) {
            if (detector[i] == target_detector && layer[i] == target_layer) {
                out.push_back(var[i]);
            }
        }
        return out;
    };
}
#ifndef RECTRAJ_H
#define RECTRAJ_H

#include <vector>
#include <string>
#include <functional>

struct RECTraj {
    static const std::vector<std::string>& All();
    static const std::vector<std::string>& Extend();

    using AllTypes = std::tuple<
        const std::vector<int>&,      // pindex
        const std::vector<int>&,      // index
        const std::vector<int>&,      // detector
        const std::vector<int>&,      // layer
        const std::vector<float>&,    // x
        const std::vector<float>&,    // y
        const std::vector<float>&,    // z
        const std::vector<float>&,    // cx
        const std::vector<float>&,    // cy
        const std::vector<float>&,    // cz
        const std::vector<float>&,    // path
        const std::vector<float>&     // edge
    >;
};

std::function<std::vector<float>(const std::vector<float>&, const std::vector<int>&, const std::vector<int>&)> 
get_RECTraj_float_var(int target_detector, int target_layer);

std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&, const std::vector<int>&)> 
get_RECTraj_int_var(int target_detector, int target_layer);

#endif
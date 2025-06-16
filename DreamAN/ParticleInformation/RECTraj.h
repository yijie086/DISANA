#ifndef RECTRAJ_H
#define RECTRAJ_H

#include <vector>
#include <string>
#include <functional>

struct RECTraj {
    static const std::vector<std::string>& All();
    static const std::vector<std::string>& Extend();
        // Minimal set needed for matching to REC::Particle
    static const std::vector<std::string>& ForFiducialCut();

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
                                 const int& REC_Particle_num)> RECTrajXYZ(int target_detector, int target_layer, int xyz);

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
                                const int& REC_Particle_num)> RECTrajedge(int target_detector, int target_layer);

std::function<std::vector<float>(const std::vector<float>&,
                                 const std::vector<int16_t>&,
                                 const std::vector<int16_t>&,
                                 const std::vector<int16_t>&,
                                 const std::vector<int>&,
                                 const std::vector<int>&,
                                 const std::vector<float>&,
                                 const std::vector<int>&)>
get_RECTraj_float_var(int target_detector, int target_layer, int target_pid, int target_charge);



#endif
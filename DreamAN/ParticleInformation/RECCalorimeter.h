#ifndef RECCALORIMETER_H
#define RECCALORIMETER_H

#include <vector>
#include <string>
#include <functional>

struct RECCalorimeter {
    static const std::vector<std::string>& All();
    static const std::vector<std::string>& Extend();

    using AllTypes = std::tuple<
        const std::vector<int16_t>&,      // index
        const std::vector<int16_t>&,      // pindex
        const std::vector<int16_t>&,      // detector
        const std::vector<int16_t>&,      // sector
        const std::vector<int16_t>&,      // layer
        const std::vector<float>&,    // energy
        const std::vector<float>&,    // time
        const std::vector<float>&,    // path
        const std::vector<float>&,    // chi2
        const std::vector<float>&,    // x
        const std::vector<float>&,    // y
        const std::vector<float>&,    // z
        const std::vector<float>&,    // hx
        const std::vector<float>&,    // hy
        const std::vector<float>&,    // hz
        const std::vector<float>&,    // lu
        const std::vector<float>&,    // lv
        const std::vector<float>&,    // lw
        const std::vector<float>&,    // du
        const std::vector<float>&,    // dv
        const std::vector<float>&,    // dw
        const std::vector<float>&,    // m2u
        const std::vector<float>&,    // m2v
        const std::vector<float>&,    // m2w
        const std::vector<float>&,    // m3u
        const std::vector<float>&,    // m3v
        const std::vector<float>&,    // m3w
        const std::vector<int>&       // status
    >;
};

std::function<std::vector<float>(const std::vector<int16_t>&,      // index
                                 const std::vector<int16_t>&,      // pindex
                                 const std::vector<int16_t>&,      // detector
                                 const std::vector<int16_t>&,      // sector
                                 const std::vector<int16_t>&,      // layer
                                 const std::vector<float>&,    // energy
                                 const std::vector<float>&,    // time
                                 const std::vector<float>&,    // path
                                 const std::vector<float>&,    // chi2
                                 const std::vector<float>&,    // x
                                 const std::vector<float>&,    // y
                                 const std::vector<float>&,    // z
                                 const std::vector<float>&,    // hx
                                 const std::vector<float>&,    // hy
                                 const std::vector<float>&,    // hz
                                 const std::vector<float>&,    // lu
                                 const std::vector<float>&,    // lv
                                 const std::vector<float>&,    // lw
                                 const std::vector<float>&,    // du
                                 const std::vector<float>&,    // dv
                                 const std::vector<float>&,    // dw
                                 const std::vector<float>&,    // m2u
                                 const std::vector<float>&,    // m2v
                                 const std::vector<float>&,    // m2w
                                 const std::vector<float>&,    // m3u
                                 const std::vector<float>&,    // m3v
                                 const std::vector<float>&,    // m3w
                                 const std::vector<int>&,      // status
                                 const int& REC_Particle_num)>
RECCalorimeterluvw(int target_detector, int target_layer, int uvw);

#endif // RECCALORIMETER_H
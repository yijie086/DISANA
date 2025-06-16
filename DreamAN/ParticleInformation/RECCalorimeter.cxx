#include "RECCalorimeter.h"
#include <iostream>
#include <cmath>
#include <functional>

const std::vector<std::string>& RECCalorimeter::All() {
    static const std::vector<std::string> names = {
        "REC_Calorimeter_index",
        "REC_Calorimeter_pindex",
        "REC_Calorimeter_detector",
        "REC_Calorimeter_sector",
        "REC_Calorimeter_layer",
        "REC_Calorimeter_energy",
        "REC_Calorimeter_time",
        "REC_Calorimeter_path",
        "REC_Calorimeter_chi2",
        "REC_Calorimeter_x",
        "REC_Calorimeter_y",
        "REC_Calorimeter_z",
        "REC_Calorimeter_hx",
        "REC_Calorimeter_hy",
        "REC_Calorimeter_hz",
        "REC_Calorimeter_lu",
        "REC_Calorimeter_lv",
        "REC_Calorimeter_lw",
        "REC_Calorimeter_du",
        "REC_Calorimeter_dv",
        "REC_Calorimeter_dw",
        "REC_Calorimeter_m2u",
        "REC_Calorimeter_m2v",
        "REC_Calorimeter_m2w",
        "REC_Calorimeter_m3u",
        "REC_Calorimeter_m3v",
        "REC_Calorimeter_m3w",
        "REC_Calorimeter_status"
    };
    return names;
}
const std::vector<std::string>& RECCalorimeter::ForFiducialCut() {
    static const std::vector<std::string> minimal = {
        "REC_Calorimeter_sector",  // sector
        "REC_Calorimeter_energy",
        "REC_Calorimeter_time",
        "REC_Calorimeter_path",
        "REC_Calorimeter_chi2",
        "REC_Calorimeter_x",
        "REC_Calorimeter_y",
        "REC_Calorimeter_z",
        "REC_Calorimeter_hx",
        "REC_Calorimeter_hy",
        "REC_Calorimeter_hz",
        "REC_Calorimeter_lu",
        "REC_Calorimeter_lv",
        "REC_Calorimeter_lw",
        "REC_Calorimeter_du",
        "REC_Calorimeter_dv",
        "REC_Calorimeter_dw",
        "REC_Calorimeter_m2u",
        "REC_Calorimeter_m2v",
        "REC_Calorimeter_m2w",
        "REC_Calorimeter_m3u",
        "REC_Calorimeter_m3v",
        "REC_Calorimeter_m3w",
        "REC_Calorimeter_status"
    };
    return minimal;
}


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
RECCalorimeterluvw(int target_detector, int target_layer, int uvw) {
    return [target_detector, target_layer, uvw](
                  const std::vector<int16_t>& index,
                  const std::vector<int16_t>& pindex,
                  const std::vector<int16_t>& detector,
                  const std::vector<int16_t>& sector,
                  const std::vector<int16_t>& layer,
                  const std::vector<float>& energy,
                  const std::vector<float>& time,
                  const std::vector<float>& path,
                  const std::vector<float>& chi2,
                  const std::vector<float>& x,
                  const std::vector<float>& y,
                  const std::vector<float>& z,
                  const std::vector<float>& hx,
                  const std::vector<float>& hy,
                  const std::vector<float>& hz,
                  const std::vector<float>& lu,
                  const std::vector<float>& lv,
                  const std::vector<float>& lw,
                  const std::vector<float>& du,
                  const std::vector<float>& dv,
                  const std::vector<float>& dw,
                  const std::vector<float>& m2u,
                  const std::vector<float>& m2v,
                  const std::vector<float>& m2w,
                  const std::vector<float>& m3u,
                  const std::vector<float>& m3v,
                  const std::vector<float>& m3w,
                  const std::vector<int>& status,
                  const int& REC_Particle_num) -> std::vector<float> {
        // Initialize return_values with size REC_Particle_num and default value 9999.0
        std::vector<float> return_values(REC_Particle_num, 9999.0);
        
        for (size_t i = 0; i < pindex.size(); ++i) {
            //std::cout << "index: "<< index[i]<< " pindex: " << pindex[i] << " detector: " << detector[i] << " layer: " << layer[i] << std::endl;
            //std::cout << "lu: " << lu[i] << " lv: " << lv[i] << " lw: " << lw[i] << std::endl;
            if (detector[i] == target_detector && layer[i] == target_layer && uvw == 1) {
                return_values[pindex[i]] = lu[i];
            } else if (detector[i] == target_detector && layer[i] == target_layer && uvw == 2) {
                return_values[pindex[i]] = lv[i];
            } else if (detector[i] == target_detector && layer[i] == target_layer && uvw == 3) {
                return_values[pindex[i]] = lw[i];
            }
        }

        return return_values;
    };
}

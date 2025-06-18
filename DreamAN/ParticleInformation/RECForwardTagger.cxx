#include "RECForwardTagger.h"
#include <iostream>
#include <cmath>
#include <functional>

const std::vector<std::string>& RECForwardTagger::All() {
    static const std::vector<std::string> names = {
        "REC_ForwardTagger_index",
        "REC_ForwardTagger_pindex",
        "REC_ForwardTagger_detector",
        "REC_ForwardTagger_layer",
        "REC_ForwardTagger_energy",
        "REC_ForwardTagger_time",
        "REC_ForwardTagger_path",
        "REC_ForwardTagger_chi2",
        "REC_ForwardTagger_x",
        "REC_ForwardTagger_y",
        "REC_ForwardTagger_z",
        "REC_ForwardTagger_dx",
        "REC_ForwardTagger_dy",
        "REC_ForwardTagger_radius",
        "REC_ForwardTagger_size",
        "REC_ForwardTagger_status",
    };
    return names;
}


const std::vector<std::string>& RECForwardTagger::Extend() {
    static const std::vector<std::string> extended = [] {
        std::vector<std::string> base = RECForwardTagger::All();
        // 如果需要扩展更多变量，可以在这里添加
        return base;
    }();
    return extended;
}


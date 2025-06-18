#ifndef RECFORWARDTAGGER_H
#define RECFORWARDTAGGER_H

#include <vector>
#include <string>
#include <functional>

struct RECForwardTagger {
    static const std::vector<std::string>& All();
    static const std::vector<std::string>& Extend();
        // Minimal set needed for matching to REC::Particle

    using AllTypes = std::tuple<
        const std::vector<int>&,      // index
        const std::vector<int>&,      // pindex
        const std::vector<int>&,      // detector
        const std::vector<int>&,      // layer
        const std::vector<float>&,    // energy
        const std::vector<float>&,    // time
        const std::vector<float>&,    // path
        const std::vector<float>&,    // chi2
        const std::vector<float>&,    // x
        const std::vector<float>&,    // y
        const std::vector<float>&,    // z
        const std::vector<float>&,    // dx
        const std::vector<float>&,    // dy
        const std::vector<float>&,    // radius
        const std::vector<int>&,    // size
        const std::vector<int>&     // edge
    >;
};



#endif
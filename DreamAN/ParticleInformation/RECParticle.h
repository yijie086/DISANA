#ifndef RECPARTICLE_H
#define RECPARTICLE_H

#include <vector>
#include <string>
#include <functional>

struct RECParticle {
    static const std::vector<std::string>& All();
    static const std::vector<std::string>& Extend();


    using AllTypes = std::tuple<
        const std::vector<int>&,      // pid
        const std::vector<float>&,    // px
        const std::vector<float>&,    // py
        const std::vector<float>&,    // pz
        const std::vector<float>&,    // vx
        const std::vector<float>&,    // vy
        const std::vector<float>&,    // vz
        const std::vector<float>&,    // vt
        const std::vector<short>&,      // charge
        const std::vector<float>&,    // beta
        const std::vector<float>&,    // chi2pid
        const std::vector<short>&       // status
    >;
};



std::function<std::vector<float>(const std::vector<float>&, const std::vector<int>&, const std::vector<short>&, const std::vector<float>&, const std::vector<int>&)> 
get_RECParticle_float_var(int target_pid, int target_charge);


std::function<std::vector<int>(const std::vector<int>&, const std::vector<int>&, const std::vector<short>&, const std::vector<float>&, const std::vector<int>&)> 
get_RECParticle_int_var(int target_pid, int target_charge);

#endif
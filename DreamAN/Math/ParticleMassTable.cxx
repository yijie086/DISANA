#include "ParticleMassTable.h"
#include <stdexcept>

// Define the getParticleMass function
double getParticleMass(int pid) {
    // Define the pid and mass (unit: GeV/c^2) of common particles
    static const std::unordered_map<int, double> particleMasses = {
        {11, 0.000511},    // Electron
        {-11, 0.000511},   // Positron
        {2212, 0.938272},  // Proton
        {-2212, 0.938272}, // Anti-Proton
        {211, 0.13957},    // Pi+
        {-211, 0.13957},   // Pi-
        {22, 0.0},         // Photon
        {13, 0.10566},     // Muon
        {-13, 0.10566},    // Anti-Muon
        {111, 0.13498},    // Pi0
        {2112, 0.939565},  // Neutron
        {-2112, 0.939565}  // Anti-Neutron
    };

    // Look up the mass corresponding to the pid
    auto it = particleMasses.find(pid);
    if (it != particleMasses.end()) {
        return it->second; // Return the mass
    } else {
        throw std::invalid_argument("Unknown pid: " + std::to_string(pid)); // Throw an exception for unknown pid
    }
}
// Initialize He3DIS generator-level kinematics in the DISANA RDataFrame style.
// This file intentionally does not make plots yet.
//
// Run from DISANA/rootmacros with:
//   clas12root DISANA_XplotterHe3.cpp

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TString.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace He3DISANA {

using ROOT::VecOps::RVec;

constexpr double kNucleonMass = 0.93892;

void RequireColumns(ROOT::RDF::RNode& dataframe,
                    const std::vector<std::string>& columns) {
    std::vector<std::string> missing;
    for (const auto& column : columns) {
        if (!dataframe.HasColumn(column)) missing.push_back(column);
    }

    if (!missing.empty()) {
        std::string message = "input tree is missing required branch(es): ";
        for (std::size_t i = 0; i < missing.size(); ++i) {
            if (i != 0) message += ", ";
            message += missing[i];
        }
        throw std::runtime_error(message);
    }
}

float FindParticleComponent(const RVec<int>& pid,
                            const RVec<float>& component,
                            int requested_pid) {
    const std::size_t size = std::min(pid.size(), component.size());
    for (std::size_t i = 0; i < size; ++i) {
        if (pid[i] == requested_pid) return component[i];
    }
    return -999.0F;
}

double Momentum(float px, float py, float pz) {
    return std::sqrt(static_cast<double>(px) * px +
                     static_cast<double>(py) * py +
                     static_cast<double>(pz) * pz);
}

double Theta(float px, float py, float pz) {
    const double momentum = Momentum(px, py, pz);
    if (momentum == 0.0) return 0.0;
    return std::acos(std::clamp(static_cast<double>(pz) / momentum,
                               -1.0, 1.0));
}

double Phi(float px, float py) {
    double phi = std::atan2(py, px);
    return phi < 0.0 ? phi + 2.0 * M_PI : phi;
}

RVec<int> SpectatorPid(const RVec<int>& pid) {
    RVec<int> result;
    result.reserve(pid.size());
    for (const int value : pid) {
        if (value == 2212 || value == 2112 || value == 1000010020) {
            result.push_back(value);
        }
    }
    return result;
}

RVec<float> SpectatorComponent(const RVec<int>& pid,
                               const RVec<float>& component) {
    RVec<float> result;
    const std::size_t size = std::min(pid.size(), component.size());
    result.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        if (pid[i] == 2212 || pid[i] == 2112 || pid[i] == 1000010020) {
            result.push_back(component[i]);
        }
    }
    return result;
}

RVec<double> SpectatorMomentum(const RVec<float>& px,
                               const RVec<float>& py,
                               const RVec<float>& pz) {
    RVec<double> result;
    const std::size_t size = std::min({px.size(), py.size(), pz.size()});
    result.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        result.push_back(Momentum(px[i], py[i], pz[i]));
    }
    return result;
}

}  // namespace He3DISANA

struct He3FinalStateDataFrames {
    ROOT::RDF::RNode all;
    ROOT::RDF::RNode epp;
    ROOT::RDF::RNode ed;
    ROOT::RDF::RNode epn;
};

He3FinalStateDataFrames InitGenKinematicsHe3(
    const std::string& filename,
    const std::string& treename = "MC",
    double beam_energy = 10.6) {
    using namespace He3DISANA;

    ROOT::RDataFrame rdf(treename, filename);
    ROOT::RDF::RNode df(rdf);

    RequireColumns(
        df,
        {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
         "MC_Particle_pz", "MC_Event_helicity",
         "MC_Event_targetPolarization", "MC_Event_xB", "MC_Event_y",
         "MC_Event_W2", "MC_Event_Q2", "MC_Event_nu"});

    // Follow DISANA_Xplotter2csv.cpp::InitGenKinematics: first select the
    // generated particles, then define their kinematics, then define DIS
    // variables. He3DIS has an inclusive scattered electron plus spectators,
    // rather than the exclusive e' p' gamma final state used by the DVCS code.
    df = df
        .Define("ele_px",
                [](const RVec<int>& pid, const RVec<float>& value) {
                    return FindParticleComponent(pid, value, 11);
                },
                {"MC_Particle_pid", "MC_Particle_px"})
        .Define("ele_py",
                [](const RVec<int>& pid, const RVec<float>& value) {
                    return FindParticleComponent(pid, value, 11);
                },
                {"MC_Particle_pid", "MC_Particle_py"})
        .Define("ele_pz",
                [](const RVec<int>& pid, const RVec<float>& value) {
                    return FindParticleComponent(pid, value, 11);
                },
                {"MC_Particle_pid", "MC_Particle_pz"})
        .Filter("ele_px != -999.0f", "He3DIS: generated electron exists")
        .Define("recel_p", Momentum, {"ele_px", "ele_py", "ele_pz"})
        .Define("recel_theta", Theta, {"ele_px", "ele_py", "ele_pz"})
        .Define("recel_phi", Phi, {"ele_px", "ele_py"})
        .Define("recel_vz", []() { return 0.0; })
        .Define("ele_det_region", []() { return 1; })
        .Define("spectator_pid", SpectatorPid, {"MC_Particle_pid"})
        .Define("spectator_px", SpectatorComponent,
                {"MC_Particle_pid", "MC_Particle_px"})
        .Define("spectator_py", SpectatorComponent,
                {"MC_Particle_pid", "MC_Particle_py"})
        .Define("spectator_pz", SpectatorComponent,
                {"MC_Particle_pid", "MC_Particle_pz"})
        .Define("spectator_p", SpectatorMomentum,
                {"spectator_px", "spectator_py", "spectator_pz"})
        .Define("nSpectators",
                [](const RVec<int>& pid) {
                    return static_cast<int>(pid.size());
                },
                {"spectator_pid"})
        .Define("he3_final_state",
                [](const RVec<int>& spectator_pid) {
                    int n_proton = 0;
                    int n_neutron = 0;
                    int n_deuteron = 0;
                    for (const int pid : spectator_pid) {
                        if (pid == 2212) ++n_proton;
                        if (pid == 2112) ++n_neutron;
                        if (pid == 1000010020) ++n_deuteron;
                    }
                    if (n_proton == 2 && n_neutron == 0 &&
                        n_deuteron == 0) return 1;  // e' + p + p
                    if (n_proton == 0 && n_neutron == 0 &&
                        n_deuteron == 1) return 2;  // e' + d
                    if (n_proton == 1 && n_neutron == 1 &&
                        n_deuteron == 0) return 3;  // e' + p + n
                    return 0;
                },
                {"spectator_pid"})
        // Preserve the generator-header values under explicit gen_* names.
        .Define("gen_xB", [](double value) { return value; },
                {"MC_Event_xB"})
        .Define("gen_Q2", [](double value) { return value; },
                {"MC_Event_Q2"})
        .Define("gen_W", [](double w2) {
                    return std::sqrt(std::max(0.0, w2));
                },
                {"MC_Event_W2"})
        .Define("gen_nu", [](double value) { return value; },
                {"MC_Event_nu"})
        .Define("gen_y", [](double value) { return value; },
                {"MC_Event_y"})
        // Standard DISANA column names, recalculated from the scattered
        // electron for the supplied beam energy.
        .Define("Q2",
                [beam_energy](double scattered_energy, double theta) {
                    return 2.0 * beam_energy * scattered_energy *
                           (1.0 - std::cos(theta));
                },
                {"recel_p", "recel_theta"})
        .Define("nu",
                [beam_energy](double scattered_energy) {
                    return beam_energy - scattered_energy;
                },
                {"recel_p"})
        .Define("xB",
                [](double q2, double nu) {
                    return nu > 0.0
                               ? q2 / (2.0 * kNucleonMass * nu)
                               : -999.0;
                },
                {"Q2", "nu"})
        .Define("W",
                [](double q2, double nu) {
                    const double w2 = kNucleonMass * kNucleonMass +
                                      2.0 * kNucleonMass * nu - q2;
                    return w2 >= 0.0 ? std::sqrt(w2) : -999.0;
                },
                {"Q2", "nu"})
        .Define("y",
                [beam_energy](double nu) {
                    return beam_energy > 0.0 ? nu / beam_energy : -999.0;
                },
                {"nu"})
        .Define("beam_target_spin",
                [](Short_t beam, Short_t target) {
                    return static_cast<int>(beam * target);
                },
                {"MC_Event_helicity", "MC_Event_targetPolarization"});

    // These filtered nodes are lazy views of the same event graph. No event
    // data are copied, and downstream actions can process them concurrently.
    return {
        df,
        df.Filter("he3_final_state == 1", "He3 final state: e' p p"),
        df.Filter("he3_final_state == 2", "He3 final state: e' d"),
        df.Filter("he3_final_state == 3", "He3 final state: e' p n")
    };
}

void DISANA_XplotterHe3() {
    ROOT::EnableImplicitMT(40);

    std::string input_path_from_He3DIS_mc ="/work/clas12/yijie/Simulation/He3DIS/He3DIS";
    std::string filename_He3DIS_mc = Form("%s/CLAS12_3He_10p6_rdf.root",input_path_from_He3DIS_mc.c_str());

    std::string treename_He3DIS_mc = "MC";
    float beam_energy = 10.6;

    He3FinalStateDataFrames df_He3DIS_mc_init = InitGenKinematicsHe3(filename_He3DIS_mc, treename_He3DIS_mc, beam_energy);

    DISANAcomparer comparer;
    comparer.AddModelHe3(df_He3DIS_mc_init.all);
    comparer.PlotParticleKinematicHe3();
    comparer.PlotDISKinematicHe3();

    auto all_count = df_He3DIS_mc_init.all.Count();
    auto epp_count = df_He3DIS_mc_init.epp.Count();
    auto ed_count = df_He3DIS_mc_init.ed.Count();
    auto epn_count = df_He3DIS_mc_init.epn.Count();

    std::cout << "Initialized " << all_count.GetValue()
              << " He3DIS events with " << beam_energy
              << " GeV beam energy using 40 ROOT threads.\n";
    std::cout << "  e' + p + p : " << epp_count.GetValue() << '\n'
              << "  e' + d     : " << ed_count.GetValue() << '\n'
              << "  e' + p + n : " << epn_count.GetValue() << '\n';
    std::cout << "The returned RDataFrame nodes have Q2, xB, W, nu, y, "
                 "electron, spectator, and polarization columns ready for "
                 "analysis.\n";
}

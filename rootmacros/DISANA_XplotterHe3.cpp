// He3DIS generator-level analysis for the DISANA environment.
//
// Example:
//   root -l -q 'rootmacros/DISANA_XplotterHe3.cpp(
//       "he3dis_rdf.root","He3DIS_plots")'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace He3Plotter {

using ROOT::VecOps::RVec;

double momentum(float px, float py, float pz) {
    return std::sqrt(static_cast<double>(px) * px +
                     static_cast<double>(py) * py +
                     static_cast<double>(pz) * pz);
}

double electron_value(const RVec<int>& pid, const RVec<float>& px,
                      const RVec<float>& py, const RVec<float>& pz,
                      int requested_value) {
    const auto size = std::min({pid.size(), px.size(), py.size(), pz.size()});
    for (std::size_t i = 0; i < size; ++i) {
        if (pid[i] != 11) continue;
        const double p = momentum(px[i], py[i], pz[i]);
        if (requested_value == 0) return p;
        if (p == 0.0) return 0.0;
        const double cosine = std::clamp(static_cast<double>(pz[i]) / p,
                                         -1.0, 1.0);
        return std::acos(cosine) * 180.0 / M_PI;
    }
    return -999.0;
}

RVec<double> spectator_momenta(const RVec<int>& pid, const RVec<float>& px,
                               const RVec<float>& py,
                               const RVec<float>& pz) {
    RVec<double> result;
    const auto size = std::min({pid.size(), px.size(), py.size(), pz.size()});
    result.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        // He3DIS spectators are protons, neutrons, or deuterons.
        if (pid[i] == 2212 || pid[i] == 2112 || pid[i] == 1000010020) {
            result.push_back(momentum(px[i], py[i], pz[i]));
        }
    }
    return result;
}

RVec<int> spectator_pids(const RVec<int>& pid) {
    RVec<int> result;
    result.reserve(pid.size());
    for (const int value : pid) {
        if (value == 2212) {
            result.push_back(1);  // proton
        } else if (value == 2112) {
            result.push_back(2);  // neutron
        } else if (value == 1000010020) {
            result.push_back(3);  // deuteron
        }
    }
    return result;
}

void require_columns(const ROOT::RDF::RNode& dataframe,
                     const std::vector<std::string>& columns) {
    std::vector<std::string> missing;
    for (const auto& column : columns) {
        if (!dataframe.HasColumn(column)) missing.push_back(column);
    }
    if (!missing.empty()) {
        std::string message = "input MC tree is missing required branch(es): ";
        for (std::size_t i = 0; i < missing.size(); ++i) {
            if (i != 0) message += ", ";
            message += missing[i];
        }
        throw std::runtime_error(message);
    }
}

}  // namespace He3Plotter

void DISANA_XplotterHe3(const char* input_file = "he3dis_rdf.root",
                        const char* output_directory = "He3DIS_plots") {
    using namespace He3Plotter;

    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);
    gSystem->mkdir(output_directory, true);

    ROOT::RDataFrame input("MC", input_file);
    ROOT::RDF::RNode dataframe(input);

    require_columns(
        dataframe,
        {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
         "MC_Particle_pz", "MC_Event_helicity",
         "MC_Event_targetPolarization", "MC_Event_xB", "MC_Event_y",
         "MC_Event_W2", "MC_Event_Q2", "MC_Event_nu"});

    auto df = dataframe
        .Define("he3_electron_p",
                [](const RVec<int>& pid, const RVec<float>& px,
                   const RVec<float>& py, const RVec<float>& pz) {
                    return electron_value(pid, px, py, pz, 0);
                },
                {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
                 "MC_Particle_pz"})
        .Define("he3_electron_theta",
                [](const RVec<int>& pid, const RVec<float>& px,
                   const RVec<float>& py, const RVec<float>& pz) {
                    return electron_value(pid, px, py, pz, 1);
                },
                {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
                 "MC_Particle_pz"})
        .Define("he3_spectator_p", spectator_momenta,
                {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py",
                 "MC_Particle_pz"})
        .Define("he3_spectator_type", spectator_pids,
                {"MC_Particle_pid"})
        .Define("he3_spin_product",
                [](Short_t beam, Short_t target) {
                    return static_cast<double>(beam * target);
                },
                {"MC_Event_helicity", "MC_Event_targetPolarization"})
        .Define("he3_spin_state",
                [](Short_t beam, Short_t target) {
                    if (beam > 0 && target > 0) return 1;
                    if (beam > 0 && target < 0) return 2;
                    if (beam < 0 && target > 0) return 3;
                    if (beam < 0 && target < 0) return 4;
                    return 0;
                },
                {"MC_Event_helicity", "MC_Event_targetPolarization"});

    const auto event_count = df.Count();

    auto h_xb = df.Histo1D(
        {"h_xB", "He3DIS;x_{B};Events", 100, 0.0, 1.0}, "MC_Event_xB");
    auto h_q2 = df.Histo1D(
        {"h_Q2", "He3DIS;Q^{2} [GeV^{2}];Events", 100, 0.0, 12.0},
        "MC_Event_Q2");
    auto h_w2 = df.Histo1D(
        {"h_W2", "He3DIS;W^{2} [GeV^{2}];Events", 100, 0.0, 25.0},
        "MC_Event_W2");
    auto h_y = df.Histo1D(
        {"h_y", "He3DIS;y;Events", 100, 0.0, 1.0}, "MC_Event_y");
    auto h_nu = df.Histo1D(
        {"h_nu", "He3DIS;#nu [GeV];Events", 100, 0.0, 11.0},
        "MC_Event_nu");
    auto h_q2_xb = df.Histo2D(
        {"h_Q2_vs_xB", "He3DIS;x_{B};Q^{2} [GeV^{2}]", 100, 0.0,
         1.0, 100, 0.0, 12.0},
        "MC_Event_xB", "MC_Event_Q2");

    auto h_electron_p = df.Histo1D(
        {"h_electron_p", "Scattered electron;p_{e'} [GeV];Events", 100,
         0.0, 11.0},
        "he3_electron_p");
    auto h_electron_theta = df.Histo1D(
        {"h_electron_theta", "Scattered electron;#theta_{e'} [deg];Events",
         100, 0.0, 50.0},
        "he3_electron_theta");
    auto h_spectator_p = df.Histo1D(
        {"h_spectator_p", "Spectators;p [GeV];Particles", 120, 0.0, 1.2},
        "he3_spectator_p");
    auto h_spectator_type = df.Histo1D(
        {"h_spectator_type", "Spectator composition;particle;Particles", 3,
         0.5, 3.5},
        "he3_spectator_type");
    auto h_spin_state = df.Histo1D(
        {"h_spin_state", "Beam-target polarization;state;Events", 4, 0.5,
         4.5},
        "he3_spin_state");
    auto p_double_spin_xb = df.Profile1D(
        {"p_double_spin_xB",
         "Raw double-spin observable;x_{B};<#it{h}_{e}#it{h}_{T}>", 20,
         0.0, 1.0, -1.1, 1.1},
        "MC_Event_xB", "he3_spin_product");

    // Trigger the event loop before drawing and report the exact entry count.
    std::cout << "Read " << *event_count << " He3DIS events from "
              << input_file << '\n';

    h_spectator_type->GetXaxis()->SetBinLabel(1, "p");
    h_spectator_type->GetXaxis()->SetBinLabel(2, "n");
    h_spectator_type->GetXaxis()->SetBinLabel(3, "d");
    h_spin_state->GetXaxis()->SetBinLabel(1, "beam+, target+");
    h_spin_state->GetXaxis()->SetBinLabel(2, "beam+, target-");
    h_spin_state->GetXaxis()->SetBinLabel(3, "beam-, target+");
    h_spin_state->GetXaxis()->SetBinLabel(4, "beam-, target-");

    const std::string directory(output_directory);

    TCanvas kinematics("c_he3_kinematics", "He3DIS kinematics", 1500, 900);
    kinematics.Divide(3, 2);
    kinematics.cd(1); h_xb->Draw("hist");
    kinematics.cd(2); h_q2->Draw("hist");
    kinematics.cd(3); h_w2->Draw("hist");
    kinematics.cd(4); h_y->Draw("hist");
    kinematics.cd(5); h_nu->Draw("hist");
    kinematics.cd(6); h_q2_xb->Draw("colz");
    kinematics.SaveAs((directory + "/He3DIS_kinematics.pdf").c_str());

    TCanvas particles("c_he3_particles", "He3DIS particles", 1500, 900);
    particles.Divide(3, 2);
    particles.cd(1); h_electron_p->Draw("hist");
    particles.cd(2); h_electron_theta->Draw("hist");
    particles.cd(3); h_spectator_p->Draw("hist");
    particles.cd(4); h_spectator_type->Draw("hist");
    particles.cd(5); h_spin_state->Draw("hist");
    particles.cd(6); p_double_spin_xb->SetMinimum(-1.1);
    p_double_spin_xb->SetMaximum(1.1);
    p_double_spin_xb->Draw();
    particles.SaveAs((directory + "/He3DIS_particles_polarization.pdf").c_str());

    TFile histogram_file((directory + "/He3DIS_histograms.root").c_str(),
                         "RECREATE");
    h_xb->Write();
    h_q2->Write();
    h_w2->Write();
    h_y->Write();
    h_nu->Write();
    h_q2_xb->Write();
    h_electron_p->Write();
    h_electron_theta->Write();
    h_spectator_p->Write();
    h_spectator_type->Write();
    h_spin_state->Write();
    p_double_spin_xb->Write();
    histogram_file.Close();

    std::cout << "Wrote plots and histograms to " << output_directory << '\n';
}

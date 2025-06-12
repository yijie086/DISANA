#include <string>
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>
#include <vector>
#include <chrono>

#include "../DreamAN/ParticleInformation/RECParticle.h"
#include "../DreamAN/ParticleInformation/RECTraj.h"
#include "../DreamAN/ParticleInformation/RECTrack.h"
#include "../DreamAN/ParticleInformation/RECCalorimeter.h"
#include "../DreamAN/Cuts/ElectronCut.h"
#include "../DreamAN/DrawHist/DrawAndSave.h"
#include "../DreamAN/core/FilesInPath.h"
#include "../DreamAN/core/Columns.h"
#include "../DreamAN/Cuts/EventCut.h"
#include "../DreamAN/Cuts/TrackCut.h"
#include "../DreamAN/Math/RECParticleKinematic.h"
#include "../DreamAN/Math/ParticleMassTable.h"

void Run102DCfiducialcut_2(const float beam_energy, const std::string& FilePath) {
    ROOT::EnableImplicitMT();

    std::vector<std::string> inputFiles = GetHipoFilesInPath(FilePath);
    if(inputFiles.empty()) {
        std::cout << "No .hipo files found in directory!\n";
        return;
    }

    auto ds = std::make_unique<RHipoDS>(inputFiles);
    ROOT::RDF::RNode df = ROOT::RDataFrame(std::move(ds));

    TFile* fout = new TFile("test.root", "RECREATE");

    // Cuts
    EventCut Electron_cut;
    Electron_cut.SetChargeCut(-1);
    Electron_cut.SetPIDCountCut(11, 1, 1);

    int detector_investigate = 6;
    int layer_investigate1 = 6;
    int layer_investigate2 = 18;
    int layer_investigate3 = 36;

    // Apply all invariant Defines once
    auto dfSelected = df
        .Define("REC_Particle_num", [](const std::vector<int>& pid) {
            return static_cast<int>(pid.size());
        }, {"REC_Particle_pid"})
        .Define("REC_Traj_pedgeR1", RECTrajedge(detector_investigate,layer_investigate1), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}))
        .Define("REC_Traj_pedgeR2", RECTrajedge(detector_investigate,layer_investigate2), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}))
        .Define("REC_Traj_pedgeR3", RECTrajedge(detector_investigate,layer_investigate3), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}))
        .Define("REC_Track_chi2perndf", RECTrackchi2perndf(detector_investigate), CombineColumns(RECTrack::All(), std::vector<std::string>{"REC_Particle_num"}))
        .Define("REC_Particle_theta", RECParticletheta(), RECParticle::All())
        .Define("REC_Particle_phi", RECParticlephi(), RECParticle::All())
        .Define("REC_Particle_p", RECParticleP(), RECParticle::All())
        .Cache();

   auto filter = [&](const std::vector<int>& pid,
                  const std::vector<float>& px,
                  const std::vector<float>& py,
                  const std::vector<float>& pz,
                  const std::vector<float>& vx,
                  const std::vector<float>& vy,
                  const std::vector<float>& vz,
                  const std::vector<float>& vt,
                  const std::vector<int>& charge,
                  const std::vector<float>& beta,
                  const std::vector<float>& chi2pid,
                  const std::vector<int>& status,
                  const std::vector<int>& REC_Traj_pass,
                  const std::vector<int>& REC_Calorimeter_pass) {
    return Electron_cut(pid, px, py, pz, vx, vy, vz, vt, charge, beta, chi2pid, status, REC_Traj_pass, REC_Calorimeter_pass);
};

for (int i = 0; i < 5; i++) {
    for (int j = 1; j <= 6; j++) {
        int theta_min = (i == 4) ? 40 : (5 + i * 5);
        int theta_max = 10 + i * 5;

        std::string base = "theta" + std::to_string(theta_min) + "to" + std::to_string(theta_max) + "_sector" + std::to_string(j);

        TrackCut trackCut;
        trackCut.SetSectorCut(j, 11, 6, true);

        // Set theta cut before the lambda
        Electron_cut.SetthetaCut(theta_min * M_PI / 180.0, theta_max * M_PI / 180.0);

        auto dfLoop = dfSelected
            .Define("REC_Traj_pass", trackCut.RECTrajPass(), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"},  std::vector<std::string>{"REC_Particle_num"}))
            .Define("REC_Calorimeter_pass", trackCut.RECCalorimeterPass(), CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}))
            .Filter(filter, {"pid", "px", "py", "pz", "vx", "vy", "vz", "vt", "charge", "beta", "chi2pid", "status", "REC_Traj_pass", "REC_Calorimeter_pass"});

        std::cout << "df_selected count: " << *dfLoop.Count() << std::endl;

        DrawAndSavechi2perndfvsedge(dfLoop, detector_investigate, layer_investigate1, 11, -1, fout,
            500, 0, 250, 25, 0, 25, {"REC_Traj_pedgeR1"}, base + "_R1");

        DrawAndSavechi2perndfvsedge(dfLoop, detector_investigate, layer_investigate2, 11, -1, fout,
            500, 0, 250, 25, 0, 25, {"REC_Traj_pedgeR2"}, base + "_R2");

        DrawAndSavechi2perndfvsedge(dfLoop, detector_investigate, layer_investigate3, 11, -1, fout,
            500, 0, 250, 25, 0, 25, {"REC_Traj_pedgeR3"}, base + "_R3");
    }
}

    fout->Close();

   // auto t_end = std::chrono::high_resolution_clock::now();
    //std::cout << "Total execution time: " << std::chrono::duration<double>(t_end - t_start).count() << " seconds.\n";
}

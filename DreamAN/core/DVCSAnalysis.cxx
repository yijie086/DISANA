#include "DVCSAnalysis.h"

#include <cmath>
#include <stdexcept>

#include "AnalysisTaskManager.h"

DVCSAnalysis::DVCSAnalysis() : fHistPhotonP(nullptr), fOutFile(nullptr) {}
DVCSAnalysis::~DVCSAnalysis() {}

void DVCSAnalysis::UserCreateOutputObjects() {
}

void DVCSAnalysis::UserExec(ROOT::RDF::RNode& df) {
  using namespace std;

  if (!fTrackCuts || !fTrackCutsElectron || !fTrackCutsProton || !fTrackCutsPhoton) throw std::runtime_error("DVCSAnalysis: One or more cut pointers not set.");
  
  dfSelected.emplace(df.Define("REC_Particle_num",[](const std::vector<int>& pid) { return static_cast<int>(pid.size()); }, {"REC_Particle_pid"}));
  
  dfSelected = dfSelected->Define("REC_Traj_pass", fTrackCuts->RECTrajPass(), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}));

  *dfSelected = dfSelected->Filter(*fTrackCutsElectron, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));
  *dfSelected = dfSelected->Filter(*fTrackCutsProton, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));
  *dfSelected = dfSelected->Filter(*fTrackCutsPhoton, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));

  *dfSelected = dfSelected->Define("REC_Particle_theta", RECParticletheta(), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Particle_phi", RECParticlephi(), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Particle_p", RECParticleP(), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Event_Q2", EventQ2(fbeam_energy, 11, -1), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Event_xB", EventxB(fbeam_energy, 11, -1, getParticleMass(2212)), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Event_Nu", EventNu(fbeam_energy, 11, -1), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Event_W", EventW(fbeam_energy, 11, -1, getParticleMass(2212)), RECParticle::All());
  *dfSelected = dfSelected->Define("REC_Event_mt", Eventmt(fbeam_energy, 2212, 1, getParticleMass(2212)), RECParticle::All());

}

void DVCSAnalysis::SaveOutput() {
  if (!fOutFile || fOutFile->IsZombie()) {
    std::cerr << "DVCSAnalysis::SaveOutput: No valid output file!" << std::endl;
    return;
  }

  if (!dfSelected.has_value()) {
    std::cerr << "DVCSAnalysis::SaveOutput: dfSelected not set!" << std::endl;
    return;
  }
  
  dfSelected->Snapshot("dfSelected", "snapshot.root", 
    CombineColumns(RECParticle::Extend(), RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));

  fOutFile->cd();

  // Save plots (safe if fOutFile is open and valid)
  DrawAndSaveEventQ2(*dfSelected, 11, -1, fOutFile, 500, 0, 5, "");
  DrawAndSaveEventxB(*dfSelected, 11, -1, fOutFile, 500, 0, 1.5, "");
  DrawAndSaveEventNu(*dfSelected, 11, -1, fOutFile, 500, 0, 6, "");
  DrawAndSaveEventW(*dfSelected, 11, -1, fOutFile, 500, 0, 4, "");
  DrawAndSaveEventmt(*dfSelected, 2212, 1, fOutFile, 500, 0, 6, "");
  DrawAndSaveQ2vsxB(*dfSelected, 11, -1, fOutFile, 500, 0, 5, 500, 0, 1.5, "");

  DrawAndSaveParticleHistograms(*dfSelected, 11, -1, fOutFile, 500, 0, 1, 500, 0, 2 * M_PI, 500, 0, 9, "");
  DrawAndSaveParticleHistograms(*dfSelected, 2212, 1, fOutFile, 500, 0, 2, 500, 0, 2 * M_PI, 500, 0, 5, "");
  DrawAndSaveParticleHistograms(*dfSelected, 22, 0, fOutFile, 500, 0, 1, 500, 0, 2 * M_PI, 500, 0, 9, "");
}

///
void DVCSAnalysis::SetOutputFile(TFile* file) { fOutFile = file; }
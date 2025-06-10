#include "DVCSAnalysis.h"

#include <cmath>
#include <stdexcept>

#include "AnalysisTaskManager.h"

DVCSAnalysis::DVCSAnalysis() : fHistPhotonP(nullptr), fOutFile(nullptr) {}
DVCSAnalysis::~DVCSAnalysis() {}

void DVCSAnalysis::UserCreateOutputObjects() {}

void DVCSAnalysis::UserExec(ROOT::RDF::RNode& df) {
  using namespace std;

  int detector_investigate = 6;  // DC

  if (!fTrackCuts || !fTrackCutsElectron || !fTrackCutsProton || !fTrackCutsPhoton) throw std::runtime_error("DVCSAnalysis: One or more cut not set.");

  fTrackCutsNoFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid->SetDoFiducialCut(true);

  // Debug check
  std::cout << "Using edge cuts (R1,R2,R3): ";
  for (auto e : fTrackCutsWithFid->GetEdgeCuts()) std::cout << e << " ";
  std::cout << std::endl;

  // ----------------------------
  // Define basic new columns
  // ----------------------------

  auto dfDefs = df.Define(
                      "REC_Particle_num", [](const std::vector<int>& pid) { return static_cast<int>(pid.size()); }, std::vector<std::string>{"REC_Particle_pid"})
                    .Define("REC_Particle_theta", RECParticletheta(), RECParticle::All())
                    .Define("REC_Particle_phi", RECParticlephi(), RECParticle::All())
                    .Define("REC_Particle_p", RECParticleP(), RECParticle::All())
                    .Define("REC_Event_Q2", EventQ2(fbeam_energy, 11, -1), RECParticle::All())
                    .Define("REC_Event_xB", EventxB(fbeam_energy, 11, -1, getParticleMass(2212)), RECParticle::All())
                    .Define("REC_Event_Nu", EventNu(fbeam_energy, 11, -1), RECParticle::All())
                    .Define("REC_Event_W", EventW(fbeam_energy, 11, -1, getParticleMass(2212)), RECParticle::All())
                    .Define("REC_Event_mt", Eventmt(fbeam_energy, 2212, 1, getParticleMass(2212)), RECParticle::All());

  // ----------------------------
  // Define both pass columns for the fiducial cuts
  // ----------------------------
  auto dfDefsWithTraj = dfDefs
                            .Define("REC_Traj_pass_nofid", fTrackCutsNoFid->RECTrajPass(),
                                    CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}))
                            .Define("REC_Traj_pass_fid", fTrackCutsWithFid->RECTrajPass(),
                                    CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}))
                            .Define("REC_Calorimeter_pass_nofid", fTrackCutsNoFid->RECCalorimeterPass(),
                                    CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}))
                            .Define("REC_Calorimeter_pass_fid", fTrackCutsWithFid->RECCalorimeterPass(),
                                    CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}));

  auto cols_nofid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass_nofid"}, std::vector<std::string>{"REC_Calorimeter_pass_nofid"});
  auto cols_fid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass_fid"}, std::vector<std::string>{"REC_Calorimeter_pass_fid"});

  // ----------------------------
  // Before fiducial cut
  // ----------------------------
  auto dfBefore = dfDefsWithTraj.Filter(*fTrackCutsElectron, cols_nofid).Filter(*fTrackCutsProton, cols_nofid).Filter(*fTrackCutsPhoton, cols_nofid);

  dfSelected = dfBefore;

  // ----------------------------
  // After fiducial cut
  // ----------------------------
  if (fFiducialCut) {
    auto dfAfter = dfDefsWithTraj.Filter(*fTrackCutsElectron, cols_fid).Filter(*fTrackCutsProton, cols_fid).Filter(*fTrackCutsPhoton, cols_fid);

    dfSelected_after = dfAfter;
  }
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

  dfSelected->Snapshot("dfSelected_before", Form("%s/%s", fOutputDir.c_str(), "dfSelected_before_fiducialCuts.root"));

  if (fFiducialCut && dfSelected_after.has_value()) {
    std::cout << "output directory is : " << fOutputDir.c_str() << std::endl;
    std::cout << "Events before fiducial: " << dfSelected->Count().GetValue() << std::endl;
    std::cout << "Events after fiducial: " << dfSelected_after->Count().GetValue() << std::endl;
    dfSelected_after->Snapshot("dfSelected_after", Form("%s/%s", fOutputDir.c_str(), "dfSelected_after_fiducialCuts.root"));
  }

  fOutFile->cd();

  /*DrawAndSaveEventQ2(*dfSelected, 11, -1, fOutFile, 500, 0, 5, "");
  DrawAndSaveEventxB(*dfSelected, 11, -1, fOutFile, 500, 0, 1.5, "");
  DrawAndSaveEventNu(*dfSelected, 11, -1, fOutFile, 500, 0, 6, "");
  DrawAndSaveEventW(*dfSelected, 11, -1, fOutFile, 500, 0, 4, "");
  DrawAndSaveEventmt(*dfSelected, 2212, 1, fOutFile, 500, 0, 6, "");
  DrawAndSaveQ2vsxB(*dfSelected, 11, -1, fOutFile, 500, 0, 5, 500, 0, 1.5, "");

  DrawAndSaveParticleHistograms(*dfSelected, 11, -1, fOutFile, 500, 0, 1, 500, 0, 2 * M_PI, 500, 0, 9, "");
  DrawAndSaveParticleHistograms(*dfSelected, 2212, 1, fOutFile, 500, 0, 2, 500, 0, 2 * M_PI, 500, 0, 5, "");
  DrawAndSaveParticleHistograms(*dfSelected, 22, 0, fOutFile, 500, 0, 1, 500, 0, 2 * M_PI, 500, 0, 9, "");*/
}

void DVCSAnalysis::SetOutputFile(TFile* file) { fOutFile = file; }
void DVCSAnalysis::SetOutputDir(const std::string& dir) { fOutputDir = dir; }

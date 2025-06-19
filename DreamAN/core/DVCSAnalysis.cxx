#include "DVCSAnalysis.h"

#include <cmath>
#include <stdexcept>

#include "AnalysisTaskManager.h"

DVCSAnalysis::DVCSAnalysis(bool IsMC, bool IsReproc) : IsMC(IsMC), IsReproc(IsReproc), fHistPhotonP(nullptr), fOutFile(nullptr) {}
DVCSAnalysis::~DVCSAnalysis() {}

void DVCSAnalysis::UserCreateOutputObjects() {}

void DVCSAnalysis::UserExec(ROOT::RDF::RNode& df) {
  using namespace std;

  if (!fTrackCuts || !fTrackCutsElectron || !fTrackCutsProton || !fTrackCutsPhoton) throw std::runtime_error("DVCSAnalysis: One or more cut not set.");

  fTrackCutsNoFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid->SetDoFiducialCut(true);
  fTrackCutsWithFid->SetFiducialCutOptions(true, true);  // apply both DC and ECAL cuts

  // Debug check
  std::cout << "Using PID-specific DC edge cuts (R1, R2, R3):" << std::endl;
  for (const auto& [pid, edgeCuts] : fTrackCutsWithFid->GetEdgeCuts()) {
    std::cout << "  PID " << pid << ": ";
    for (auto e : edgeCuts) std::cout << e << " ";
    std::cout << std::endl;
  }

  std::cout << "Using PID-specific CVT edge cuts (l1, l3, l5, l7, l12):" << std::endl;
  for (const auto& [pid, edgeCuts] : fTrackCutsWithFid->GetCVTEdgeCuts()) {
    std::cout << "  PID " << pid << ": ";
    for (auto e : edgeCuts) std::cout << e << " ";
    std::cout << std::endl;
  }

  // Cache column names
  //auto colnames = df.GetColumnNames();
  auto dfDefs = df;
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_num", [](const std::vector<int>& pid) { return static_cast<int>(pid.size()); }, {"REC_Particle_pid"});
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_theta", RECParticletheta(), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_phi", RECParticlephi(), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_p", RECParticleP(), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Event_Q2", EventQ2(fbeam_energy, 11, -1), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Event_xB", EventxB(fbeam_energy, 11, -1, getParticleMass(2212)), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Event_Nu", EventNu(fbeam_energy, 11, -1), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Event_W", EventW(fbeam_energy, 11, -1, getParticleMass(2212)), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Event_mt", Eventmt(fbeam_energy, 2212, 1, getParticleMass(2212)), RECParticle::All());

  if (IsMC) {
    dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_phi_1", RECParticlephi(), RECParticle::All());
  }

  // Fiducial cuts
  auto dfDefsWithTraj = dfDefs;
  auto trajCols = CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});
  auto caloCols = CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});
  auto fwdtagCols = CombineColumns(RECForwardTagger::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});

  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_nofid", fTrackCutsNoFid->RECTrajPass(), trajCols);
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Traj_pass_fid", fTrackCutsWithFid->RECTrajPass(), trajCols);
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Calorimeter_pass_fid", fTrackCutsWithFid->RECCalorimeterPass(), caloCols);
  if (fFTonConfig) {
    dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_ForwardTagger_pass_fid", fTrackCutsWithFid->RECForwardTaggerPass(), fwdtagCols);
  }
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_fid", Columns::LogicalAND2(),
                                    CombineColumns(std::vector<std::string>{"REC_Traj_pass_fid"}, std::vector<std::string>{"REC_Calorimeter_pass_fid"}));
  if (fFTonConfig) {
    dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_fid", Columns::LogicalAND2(),
                                    CombineColumns(std::vector<std::string>{"REC_Track_pass_fid"}, std::vector<std::string>{"REC_ForwardTagger_pass_fid"}));
  }
  auto AllCols = CombineColumns(trajCols, caloCols);
  
  /*
  auto AllCols = CombineColumns(RECTraj::All(), RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_fid", fTrackCutsWithFid->RECFiducialPass(), AllCols);
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_nofid", fTrackCutsNoFid->RECFiducialPass(), AllCols);
  */

  auto cols_track_fid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass_fid"});
  auto cols_track_nofid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass_nofid"});

  // Before fiducial cut
  auto dfBefore = dfDefsWithTraj.Filter(*fTrackCutsElectron, cols_track_nofid).Filter(*fTrackCutsPhoton, cols_track_nofid).Filter(*fTrackCutsProton, cols_track_nofid);

  dfSelected = dfBefore;

  // After fiducial cut
  if (fFiducialCut) {
    auto dfAfter = dfDefsWithTraj.Filter(*fTrackCutsElectron, cols_track_fid).Filter(*fTrackCutsPhoton, cols_track_fid).Filter(*fTrackCutsProton, cols_track_fid);
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

  if (!IsReproc) dfSelected->Snapshot("dfSelected_before", Form("%s/%s", fOutputDir.c_str(), "dfSelected_before_fiducialCuts.root"));
  if (fFiducialCut && dfSelected_after.has_value()) {
    std::cout << "output directory is : " << fOutputDir.c_str() << std::endl;
    std::cout << "Events before fiducial: " << dfSelected->Count().GetValue() << std::endl;
    std::cout << "Events after fiducial: " << dfSelected_after->Count().GetValue() << std::endl;
    if (IsReproc && dfSelected_after.has_value()) {
      dfSelected_after->Snapshot("dfSelected_after_reprocessed", Form("%s/%s", fOutputDir.c_str(), "dfSelected_after_fiducialCuts_reprocessed.root"));
    } else {
      dfSelected_after->Snapshot("dfSelected_after", Form("%s/%s", fOutputDir.c_str(), "dfSelected_after_fiducialCuts.root"));
    }
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

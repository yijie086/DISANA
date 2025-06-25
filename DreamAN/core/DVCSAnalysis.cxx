#include "DVCSAnalysis.h"

#include <cmath>
#include <stdexcept>

#include "AnalysisTaskManager.h"

DVCSAnalysis::DVCSAnalysis(bool IsMC, bool IsReproc) : IsMC(IsMC), IsReproc(IsReproc), fHistPhotonP(nullptr), fOutFile(nullptr) {}
DVCSAnalysis::~DVCSAnalysis() {}

void DVCSAnalysis::UserCreateOutputObjects() {}

void DVCSAnalysis::UserExec(ROOT::RDF::RNode& df) {
  using namespace std;

  if (fMaxEvents > 0) {
        df = df.Range(0, fMaxEvents);   // only process the first fMaxEvents
  }

  if (!fTrackCuts || !fEventCuts) throw std::runtime_error("DVCSAnalysis: One or more cut not set.");

  fTrackCutsNoFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid->SetDoFiducialCut(true);
  fTrackCutsWithFid->SetFiducialCutOptions(true, true);  // apply both DC and ECAL cuts

  // Cache column names
  // auto colnames = df.GetColumnNames();
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
  // dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Event_pass","REC_Particle_pass", *fEventCuts, cols_track_fid);

  dfSelected = dfDefsWithTraj;
  dfSelected = DefineOrRedefine(*dfSelected, "EventCutResult", *fEventCuts, cols_track_nofid);
  dfSelected = DefineOrRedefine(*dfSelected, "REC_Event_pass", [](const EventCutResult& result) { return result.eventPass; }, {"EventCutResult"});
  dfSelected = DefineOrRedefine(*dfSelected, "REC_Particle_pass", [](const EventCutResult& result) { return result.particlePass; }, {"EventCutResult"});

  if (fDoInvMassCut) {
    fEventCuts->SetDoCutMotherInvMass(true);
    dfSelected = DefineOrRedefine(*dfSelected, "REC_DaughterParticle_pass", [](const EventCutResult& result) { return result.particleDaughterPass; }, {"EventCutResult"});
    dfSelected = DefineOrRedefine(*dfSelected, "REC_MotherMass", [](const EventCutResult& result) { return result.MotherMass; }, {"EventCutResult"});
  }
  dfSelected = dfSelected->Filter("REC_Event_pass");

  // After fiducial cut
  if (fFiducialCut) {
    dfSelected_afterFid = dfDefsWithTraj;
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "EventCutResult", *fEventCuts, cols_track_fid);
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_Event_pass", [](const EventCutResult& result) { return result.eventPass; }, {"EventCutResult"});
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_Particle_pass", [](const EventCutResult& result) { return result.particlePass; }, {"EventCutResult"});
    if (fDoInvMassCut) {
      fEventCuts->SetDoCutMotherInvMass(true);
      dfSelected_afterFid =
          DefineOrRedefine(*dfSelected_afterFid, "REC_DaughterParticle_pass", [](const EventCutResult& result) { return result.particleDaughterPass; }, {"EventCutResult"});
      dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_MotherMass", [](const EventCutResult& result) { return result.MotherMass; }, {"EventCutResult"});
    }
    dfSelected_afterFid = dfSelected_afterFid->Filter("REC_Event_pass");
  }

  dfSelected_afterFid_afterCorr = dfSelected_afterFid;

  if (fMomCorr && fDoMomentumCorrection) {
    std::cout << "Applying momentum correction..." << std::endl;
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_px", fMomCorr->RECParticlePxCorrected(), RECParticle::Extend());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_py", fMomCorr->RECParticlePyCorrected(), RECParticle::Extend());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_pz", fMomCorr->RECParticlePzCorrected(), RECParticle::Extend());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_theta", RECParticletheta(), RECParticle::All());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_phi", RECParticlephi(), RECParticle::All());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_p", RECParticleP(), RECParticle::All());
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

  if (!IsReproc) SafeSnapshot(*dfSelected, "dfSelected", Form("%s/%s", fOutputDir.c_str(), "dfSelected.root"));
  if (fFiducialCut && dfSelected_afterFid.has_value()) {
    std::cout << "output directory is : " << fOutputDir.c_str() << std::endl;
    std::cout << "Events selected: " << dfSelected->Count().GetValue() << std::endl;
    std::cout << "Events selected after fiducial: " << dfSelected_afterFid->Count().GetValue() << std::endl;
    if (IsReproc && dfSelected_afterFid.has_value()) {
      SafeSnapshot(*dfSelected_afterFid, "dfSelected_afterFid_reprocessed", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_reprocessed.root"));
    } else {
      SafeSnapshot(*dfSelected_afterFid, "dfSelected_afterFid", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid.root"));
    }
  }
  if (fDoMomentumCorrection && dfSelected_afterFid_afterCorr.has_value()) {
    std::cout << "Events selected after fiducial and momentum correction: " << dfSelected_afterFid_afterCorr->Count().GetValue() << std::endl;
    SafeSnapshot(*dfSelected_afterFid_afterCorr, "dfSelected_afterFid_afterCorr", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_afterCorr.root"));
  }

  fOutFile->cd();
}

void DVCSAnalysis::SetOutputFile(TFile* file) { fOutFile = file; }
void DVCSAnalysis::SetOutputDir(const std::string& dir) { fOutputDir = dir; }

#include "PhiAnalysis.h"

#include <cmath>
#include <stdexcept>

#include "AnalysisTaskManager.h"
#include "PerRunCounter.h"
#include <algorithm>
#include "ROOT/RVec.hxx"

// Returns the run and event column names from an RHipoDS-backed node.
// RHipoDS translates "RUN::config.run" -> "RUN_config_run" (:: -> _, . -> _).
static inline std::pair<std::string,std::string>
PickRunEventCols(ROOT::RDF::RNode df) {
  auto cols = df.GetColumnNames();
  auto has = [&](const std::string& n){
    return std::find(cols.begin(), cols.end(), n) != cols.end();
  };
  if (has("RUN_config_run") && has("RUN_config_event"))
    return {"RUN_config_run", "RUN_config_event"};
  if (has("RUN::config.run") && has("RUN::config.event"))
    return {"RUN::config.run", "RUN::config.event"};
  throw std::runtime_error("QADB: cannot find run/event columns");
}

PhiAnalysis::PhiAnalysis(bool IsMC, bool IsReproc, bool IsMinBook) : IsMC(IsMC), IsReproc(IsReproc), IsMinBooking(IsMinBook), fHistPhotonP(nullptr) {}
PhiAnalysis::~PhiAnalysis() {}

void PhiAnalysis::UserCreateOutputObjects() {}

void PhiAnalysis::UserExec(ROOT::RDF::RNode& df) {
  using namespace std;

  if (fMaxEvents > 0) {
    df = df.Range(0, fMaxEvents);  // only process the first fMaxEvents
  }
  if (!fTrackCuts || !fEventCuts) throw std::runtime_error("PhiAnalysis: One or more cut not set.");

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
  dfDefs = DefineOrRedefine(dfDefs, "num_events", [](ULong64_t e) { return e; }, {"rdfentry_"});


  dforginal = dfDefs;
  // Fiducial cuts
  auto dfDefsWithTraj = dfDefs;
  auto trajCols = CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});
  auto caloCols =
      CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_p"}, std::vector<std::string>{"REC_Particle_num"});
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

  auto cols_track_fid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass_fid"});
  auto cols_track_nofid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass_nofid"});

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
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_Photon_MaxE", [](const EventCutResult& result) { return result.MaxPhotonEnergyPass; }, {"EventCutResult"});

    if (fDoInvMassCut) {
      fEventCuts->SetDoCutMotherInvMass(true);
      dfSelected_afterFid =
          DefineOrRedefine(*dfSelected_afterFid, "REC_DaughterParticle_pass", [](const EventCutResult& result) { return result.particleDaughterPass; }, {"EventCutResult"});
      dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_MotherMass", [](const EventCutResult& result) { return result.MotherMass; }, {"EventCutResult"});
    }

    dfSelected_afterFid = dfSelected_afterFid->Filter("REC_Event_pass");
  }

  dfSelected_afterFid_afterCorr = dfSelected_afterFid;
  /// For Phi Analysis since the QADB may changes over time, its better to apply together with momentum correction
    // QADB cuts should be place in the first to reduce the computation load
if (fIsQADBCut && fQADBCuts) {
  std::cout << "Applying QADB cut..." << std::endl;

  auto node = *dfSelected_afterFid_afterCorr;
  auto [runCol, evCol] = PickRunEventCols(node);

  // Use a typed lambda wrapping the functor so ROOT can deduce the return type.
  // QADBCuts has multiple operator() overloads which prevents ROOT's CallableTraits
  // from deducing ret_type if the functor is passed directly.
  // RHipoDS exposes RUN::config.run as a scalar int (nrows=1 bank).
  // We pass the source columns directly — no intermediate Define nodes — so that
  // GetColumnReadersImpl for these columns is called in the single-threaded init
  // phase together with all other source columns, not lazily from worker threads.
  auto qadb = *fQADBCuts;
  node = node.Define("REC_QADB_pass",
                     [qadb](int run, int ev) mutable -> bool {
                       return qadb(run, ev);
                     },
                     {runCol, evCol})
             .Filter("REC_QADB_pass", "QADB pass");

  dfSelected_afterFid_afterCorr = node;
}


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
void PhiAnalysis::SaveOutput() {
  if (IsMC) {
    // snapshot of the MC bank for efficiency and other studies
    dforginal->Snapshot(
        "dfSelectedMC", Form("%s/%s", fOutputDir.c_str(), "dfSelectedMC.root"),
        {"num_events",          // event identity
         "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz",
         "MC_Particle_vx",  "MC_Particle_vy", "MC_Particle_vz", "MC_Particle_vt",
         "MC_Lund_pid", "MC_Lund_px", "MC_Lund_py", "MC_Lund_pz",  // full generator record
         "MC_Lund_parent", "MC_Lund_daughter",
         "MC_RecMatch_pindex", "MC_RecMatch_mcindex",               // reco<->truth links
         "MC_GenMatch_pindex", "MC_GenMatch_mcindex", "MC_GenMatch_quality",
         "MC_Event_weight", "MC_Event_pbeam", "MC_Event_ptarget", "MC_Event_ebeam"});
  }

  if (!dfSelected.has_value()) {
    std::cerr << "PhiAnalysis::SaveOutput: dfSelected not set!" << std::endl;
    return;
  }

  if (!IsReproc && !IsMinBooking) SafeSnapshot(*dfSelected, "dfSelected", Form("%s/%s", fOutputDir.c_str(), "dfSelected.root"));
  if (fFiducialCut && dfSelected_afterFid.has_value()) {
    std::cout << "output directory is : " << fOutputDir.c_str() << std::endl;
    std::cout << "Events selected: " << dfSelected->Count().GetValue() << std::endl;
    std::cout << "Events selected after fiducial: " << dfSelected_afterFid->Count().GetValue() << std::endl;
    if (IsReproc && dfSelected_afterFid.has_value()) {
      SafeSnapshot(*dfSelected_afterFid, "dfSelected_afterFid_reprocessed", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_reprocessed.root"));
    } else {
      if (!IsMinBooking) SafeSnapshot(*dfSelected_afterFid, "dfSelected_afterFid", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid.root"));
    }
  }
  if (fDoMomentumCorrection && dfSelected_afterFid_afterCorr.has_value()) {
    std::cout << "Events selected after fiducial and momentum correction: " << dfSelected_afterFid_afterCorr->Count().GetValue() << std::endl;
    if (!IsMinBooking) SafeSnapshot(*dfSelected_afterFid_afterCorr, "dfSelected_afterFid_afterCorr", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_afterCorr.root"));
  }
  if (fIsQADBCut) {
    std::cout << "\n[QADB] total accumulated charge analyzed: " << fQADBCuts->GetAccumulatedCharge() / 1e6 << " mC (Do NOT use this number if you enable MT)\n";
  }
}
void PhiAnalysis::SetOutputDir(const std::string& dir) { fOutputDir = dir; }
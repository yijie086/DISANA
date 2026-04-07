#include "PhiAnalysis.h"

#include <cmath>
#include <stdexcept>
#include <unordered_set>

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

  // Guard: verify that the essential REC::Particle columns are present.
  // With very few or sparse hipo files the bank may be entirely absent,
  // which makes RHipoDS not register the columns at all.  Crashing with an
  // opaque "Unknown column" error is hard to debug, so we emit a clear
  // warning and return early — the job will produce no output ROOT file,
  // job_ok() will mark it as failed, and it can be relaunched or skipped.
  {
    const std::vector<std::string> required = {
        "REC_Particle_pid", "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz",
        "REC_Particle_vx",  "REC_Particle_vy", "REC_Particle_vz",
        "REC_Particle_charge", "REC_Particle_beta", "REC_Particle_status"
    };
    auto existing = df.GetColumnNames();
    std::unordered_set<std::string> colSet(existing.begin(), existing.end());
    std::vector<std::string> missing;
    for (const auto& col : required)
      if (!colSet.count(col)) missing.push_back(col);

    if (!missing.empty()) {
      std::cerr << "[PhiAnalysis] WARNING: required column(s) missing from dataframe "
                   "(hipo files may have empty/absent REC::Particle bank):\n";
      for (const auto& col : missing)
        std::cerr << "    " << col << "\n";
      std::cerr << "[PhiAnalysis] Skipping analysis for this job — no output will be written.\n";
      return;   // leave all optional<RNode> members unset → SaveOutput() is a no-op
    }
  }

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
  auto trajCols    = CombineColumns(RECTraj::ForFiducialCut(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});
  auto caloCols    = CombineColumns(RECCalorimeter::ForFiducialCut(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_p"}, std::vector<std::string>{"REC_Particle_num"});
  auto fwdtagCols  = CombineColumns(RECForwardTagger::ForFiducialCut(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});

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
// ---------------------------------------------------------------------------
// MinimalColumns – the columns that every snapshot must carry so that the
// analysis can be re-run on the resulting .root file without a HIPO source.
//
// Rule of thumb: include every raw bank column that UserExec reads as an
// *input* to a Define/Redefine lambda, plus the key derived columns that
// downstream plotters need.  Intermediate helper columns (pass-booleans, the
// EventCutResult struct, etc.) that are re-derived each time are NOT needed.
// ---------------------------------------------------------------------------
std::vector<std::string> PhiAnalysis::MinimalColumns() const {
  using V = std::vector<std::string>;

  // -----------------------------------------------------------------------
  // These lists mirror exactly what the analysis lambdas actually read.
  // They must stay in sync with TrackCut::RECTrajPass(),
  // TrackCut::RECCalorimeterPass(), TrackCut::RECForwardTaggerPass(),
  // EventCut::operator(), and the downstream plotters (DISANAMMUtils.h,
  // DISANAplotter.h).
  // -----------------------------------------------------------------------

  // RECParticle::All() — all 12 columns kept because EventCut reads 9 of
  // them and the 3 unused ones (vx, vy, vt) are trivially small.
  auto cols = CombineColumns(RECParticle::All());

  // Fiducial bank columns — only what ForFiducialCut() declares, which
  // exactly matches the slimmed TrackCut lambda signatures:
  //   RECTraj::ForFiducialCut()        → 7 cols  (was ::All() = 12)
  //   RECCalorimeter::ForFiducialCut() → 8 cols  (was ::All() = 28)
  //   RECForwardTagger::ForFiducialCut()→ 5 cols  (was ::All() = 16)
  cols = CombineColumns(cols, RECTraj::ForFiducialCut());
  cols = CombineColumns(cols, RECCalorimeter::ForFiducialCut());
  if (fFTonConfig)
    cols = CombineColumns(cols, RECForwardTagger::ForFiducialCut());

  // REC::Event helicity — raw bank scalar, never Defined, used by
  // DISANAplotter.h (BSA) and DISANAMMUtils.h (helicity-weighted fills).
  cols.push_back("REC_Event_helicity");

  // Run / event IDs — both naming conventions for HIPO vs ROOT source.
  // ResolveSnapshotColumns silently drops whichever variant is absent.
  for (const auto& c : V{"RUN_config_run", "RUN_config_event",
                          "RUN::config.run", "RUN::config.event"})
    cols.push_back(c);

  // Pre-computed kinematics and analysis decision columns.
  for (const auto& c : V{
      "num_events",
      "REC_Particle_num",
      "REC_Particle_theta",
      "REC_Particle_phi",
      "REC_Particle_p",
      "REC_Particle_pass",
      "REC_Event_pass",
      "REC_MotherMass",
      "REC_DaughterParticle_pass",
  })
    cols.push_back(c);

  return cols;
}

void PhiAnalysis::SaveOutput() {
  // If UserExec returned early (missing bank columns), dforginal and dfSelected
  // are both unset. Nothing to write — exit silently.
  if (!dforginal.has_value()) {
    std::cerr << "[PhiAnalysis] SaveOutput: skipped (UserExec did not populate dataframes).\n";
    return;
  }

  if (IsMC) {
    // snapshot of the MC bank for efficiency and other studies
    dforginal->Snapshot(
        "dfSelectedMC", Form("%s/%s", fOutputDir.c_str(), "dfSelectedMC.root"),
        {"num_events",
         "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz",
         "MC_Particle_vx",  "MC_Particle_vy", "MC_Particle_vz", "MC_Particle_vt",
         "MC_Lund_pid", "MC_Lund_px", "MC_Lund_py", "MC_Lund_pz",
         "MC_Lund_parent", "MC_Lund_daughter",
         "MC_RecMatch_pindex", "MC_RecMatch_mcindex",
         "MC_GenMatch_pindex", "MC_GenMatch_mcindex", "MC_GenMatch_quality",
         "MC_Event_weight", "MC_Event_pbeam", "MC_Event_ptarget", "MC_Event_ebeam"});
  }

  if (!dfSelected.has_value()) {
    std::cerr << "PhiAnalysis::SaveOutput: dfSelected not set!" << std::endl;
    return;
  }

  if (fOptimizeColumns) {
    std::cout << "[SaveOutput] Column optimisation ON — writing only analysis-used columns.\n";
  } else {
    std::cout << "[SaveOutput] Column optimisation OFF — writing all columns.\n";
  }

  // -----------------------------------------------------------------------
  // Helper: resolve the output column list for a given dataframe node.
  //   fOptimizeColumns = true  → MinimalColumns() filtered to what exists
  //   fOptimizeColumns = false → every column except EventCutResult (the
  //                              internal struct ROOT cannot serialise)
  // -----------------------------------------------------------------------
  auto resolveColumns = [this](ROOT::RDF::RNode& node) -> std::vector<std::string> {
    if (fOptimizeColumns) {
      return ResolveSnapshotColumns(node, MinimalColumns());
    } else {
      return SafeSnapshotColumns(node, {"EventCutResult"});
    }
  };

  // -----------------------------------------------------------------------
  // Count() in RDF is always lazy — it books an action but waits for a loop
  // trigger.  Snapshot() by default is eager — it triggers the loop immediately.
  // Because all actions booked on the same graph node share one loop, booking
  // Count() BEFORE calling Snapshot() means both are computed in a single pass.
  // We avoid lazy Snapshot (fLazy=true) because *snap; as a statement does not
  // reliably trigger the loop in all ROOT versions on the farm.
  // -----------------------------------------------------------------------

  // --- dfSelected ---------------------------------------------------------
  if (!IsReproc && !IsMinBooking) {
    auto cols = resolveColumns(*dfSelected);
    auto cnt  = dfSelected->Count();           // book (lazy) — shares the loop below
    dfSelected->Snapshot("dfSelected",
                    Form("%s/%s", fOutputDir.c_str(), "dfSelected.root"), cols);  // triggers loop
    if (fFiducialCut)
      std::cout << "Events selected: " << *cnt << std::endl;  // already computed
  }

  // --- dfSelected_afterFid ------------------------------------------------
  // NOTE: dfSelected_afterFid is intentionally NOT written to disk.
  // dfSelected_afterFid_afterCorr already contains the full fiducial + QADB +
  // momentum-correction selection and is the only downstream output needed.
  // Skipping this snapshot saves one full event-loop pass and significant disk.
  if (fFiducialCut && dfSelected_afterFid.has_value() && IsReproc) {
    // In reproc mode we still print the event count for bookkeeping,
    // but do not snapshot — afterCorr is the authoritative output.
    auto cnt_fid = dfSelected_afterFid->Count();
    auto cnt_sel = dfSelected->Count();
    // Touch both counts in the same loop triggered by afterCorr snapshot below
    std::cout << "Events selected (reproc): " << *cnt_sel << std::endl;
    std::cout << "Events after fiducial   : " << *cnt_fid << std::endl;
  }

  // --- dfSelected_afterFid_afterCorr -------------------------------------
  if (fDoMomentumCorrection && dfSelected_afterFid_afterCorr.has_value()) {
    if (!IsMinBooking) {
      auto cols_corr = resolveColumns(*dfSelected_afterFid_afterCorr);
      auto cnt = dfSelected_afterFid_afterCorr->Count();   // book before triggering
      dfSelected_afterFid_afterCorr->Snapshot(
          "dfSelected_afterFid_afterCorr",
          Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_afterCorr.root"),
          cols_corr);                                       // triggers loop → computes cnt
      std::cout << "Events selected after fiducial and momentum correction: " << *cnt << std::endl;
    } else {
      std::cout << "Events selected after fiducial and momentum correction: "
                << dfSelected_afterFid_afterCorr->Count().GetValue() << std::endl;
    }
  }

  if (fIsQADBCut) {
    std::cout << "\n[QADB] total accumulated charge analyzed: "
              << fQADBCuts->GetAccumulatedCharge() / 1e6
              << " mC (Do NOT use this number if you enable MT)\n";
  }
}
void PhiAnalysis::SetOutputDir(const std::string& dir) { fOutputDir = dir; }
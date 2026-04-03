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

  // --- raw bank columns consumed by UserExec lambdas ---
  auto cols = CombineColumns(
      RECParticle::All(),        // pid, px, py, pz, vx, vy, vz, vt, charge, beta, chi2pid, status
      RECTraj::All(),            // needed by RECTrajPass (DC fiducial)
      RECCalorimeter::All()     // needed by RECCalorimeterPass (ECAL fiducial)
      //RECForwardTagger::All()    // needed by RECForwardTaggerPass (FT fiducial) No FT in phi analysis
  );

  // --- REC::Event bank scalars (raw, never computed by a Define) ---
  // REC_Event_helicity is the beam helicity (+1/-1/0) used by DISANAplotter.h
  // and DISANAMMUtils.h for BSA / helicity-weighted tree fills. It comes
  // directly from the HIPO source and is absent from every ::All() list.
  for (const auto& c : V{"REC_Event_helicity"}) {
    cols.push_back(c);
  }

  // --- run / event identifiers (both naming conventions; SelectiveSnapshot
  //     will silently drop whichever one is absent) ---
  for (const auto& c : V{"RUN_config_run",   "RUN_config_event",
                          "RUN::config.run",  "RUN::config.event"}) {
    cols.push_back(c);
  }

  // --- derived columns kept for convenience / downstream analysis ---
  for (const auto& c : V{
      "num_events",              // sequential event index
      "REC_Particle_num",        // particle multiplicity
      "REC_Particle_theta",      // derived kinematics (already recomputed if
      "REC_Particle_phi",        //   missing, but nice to have pre-computed)
      "REC_Particle_p",
      // analysis decision columns
      "REC_Particle_pass",
      "REC_Event_pass",
      "REC_Photon_MaxE",
      "REC_MotherMass",
      "REC_DaughterParticle_pass",
  }) {
    cols.push_back(c);
  }

  return cols;
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

  // -----------------------------------------------------------------------
  // Lazy snapshot options: booking Snapshot() and Count() on the same node
  // with fLazy=true lets ROOT execute BOTH in a single event-loop pass
  // instead of running the full graph twice (once for the file write, once
  // for the counter). With three snapshot targets that halves the total
  // number of event-loop runs from 6 down to 3.
  // -----------------------------------------------------------------------
  ROOT::RDF::RSnapshotOptions lazyOpts;
  lazyOpts.fLazy = true;

  // --- dfSelected ---------------------------------------------------------
  if (!IsReproc && !IsMinBooking) {
    auto cols = ResolveSnapshotColumns(*dfSelected, MinimalColumns());
    auto snap = dfSelected->Snapshot("dfSelected",
                    Form("%s/%s", fOutputDir.c_str(), "dfSelected.root"), cols, lazyOpts);
    auto cnt  = dfSelected->Count();  // booked on same node — shares the loop
    *snap;                            // ONE event loop: writes file + counts
    if (fFiducialCut)
      std::cout << "Events selected: " << *cnt << std::endl;
  }

  // --- dfSelected_afterFid ------------------------------------------------
  if (fFiducialCut && dfSelected_afterFid.has_value()) {
    std::cout << "output directory is : " << fOutputDir.c_str() << std::endl;

    auto cols_fid = ResolveSnapshotColumns(*dfSelected_afterFid, MinimalColumns());

    if (IsReproc) {
      auto snap = dfSelected_afterFid->Snapshot(
                      "dfSelected_afterFid_reprocessed",
                      Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_reprocessed.root"),
                      cols_fid, lazyOpts);
      auto cnt_fid = dfSelected_afterFid->Count();
      auto cnt_sel = dfSelected->Count();
      *snap;                          // ONE loop for dfSelected_afterFid
      *cnt_sel;                       // ONE loop for dfSelected
      std::cout << "Events selected: " << *cnt_sel << std::endl;
      std::cout << "Events selected after fiducial: " << *cnt_fid << std::endl;
    } else {
      if (!IsMinBooking) {
        auto snap = dfSelected_afterFid->Snapshot(
                        "dfSelected_afterFid",
                        Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid.root"),
                        cols_fid, lazyOpts);
        auto cnt_fid = dfSelected_afterFid->Count();
        auto cnt_sel = dfSelected->Count();
        *snap;                        // ONE loop for dfSelected_afterFid
        *cnt_sel;                     // ONE loop for dfSelected (if not run above)
        std::cout << "Events selected: " << *cnt_sel << std::endl;
        std::cout << "Events selected after fiducial: " << *cnt_fid << std::endl;
      }
    }
  }

  // --- dfSelected_afterFid_afterCorr -------------------------------------
  if (fDoMomentumCorrection && dfSelected_afterFid_afterCorr.has_value()) {
    if (!IsMinBooking) {
      auto cols_corr = ResolveSnapshotColumns(*dfSelected_afterFid_afterCorr, MinimalColumns());
      auto snap = dfSelected_afterFid_afterCorr->Snapshot(
                      "dfSelected_afterFid_afterCorr",
                      Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_afterCorr.root"),
                      cols_corr, lazyOpts);
      auto cnt = dfSelected_afterFid_afterCorr->Count();
      *snap;                          // ONE loop: writes file + counts
      std::cout << "Events selected after fiducial and momentum correction: " << *cnt << std::endl;
    } else {
      std::cout << "Events selected after fiducial and momentum correction: "
                << dfSelected_afterFid_afterCorr->Count().GetValue() << std::endl;
    }
  }

  if (fIsQADBCut) {
    std::cout << "\n[QADB] total accumulated charge analyzed: " << fQADBCuts->GetAccumulatedCharge() / 1e6 << " mC (Do NOT use this number if you enable MT)\n";
  }
}
void PhiAnalysis::SetOutputDir(const std::string& dir) { fOutputDir = dir; }
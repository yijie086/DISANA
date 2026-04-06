#!/bin/zsh
# run_disana_parallel.zsh
# Run K DISANA jobs (AnalysisPhi), optionally relaunch failed jobs, optionally run in batches, then merge.
# All outputs (WORKDIR + merged ROOT outputs) go under OUTPUT_BASE so you can run from anywhere.
#
# This version supports merging multiple possible per-job ROOT outputs:
#   dfSelected_afterFid_afterCorr.root
#   dfSelected_afterFid.root
#   dfSelectedMC.root
#   dfSelected.root
#
# It checks per job which files exist (non-empty) and merges each file type separately at the end.

set -u
set -o pipefail
setopt NO_NOMATCH   # avoid zsh "no matches found" on empty globs


usage() {
  cat <<'USAGE'
Usage:
  run_disana_parallel_legacy_defaults_cli.zsh [options]

This version keeps the script's old hardcoded defaults, but adds CLI parsing
on top so flags like -c, -i, -o, -t, --inbending, --hphm, etc. actually work.

Options:
  -i, --input-dir DIR         Override INPUT_DIR
  -o, --output-base DIR       Override OUTPUT_BASE
  -c, --config NAME           Override CFG
  -t, --threads N             Override THREADS_PER_JOB
  -n, --nfile N               Override NEVT
  -k, --jobs N                Override K
      --exe PATH              Override EXE
      --workdir NAME|PATH     Override WORKDIR
      --relaunch-only         Set RELAUNCH_ONLY=1
      --merge-partial         Set MERGE_PARTIAL=1
      --no-merge-partial      Set MERGE_PARTIAL=0
      --run-in-batches        Set RUN_IN_BATCHES=1
      --no-run-in-batches     Set RUN_IN_BATCHES=0
      --jobs-per-batch N      Override JOBS_PER_BATCH
      --two-pass              Set TWO_PASS_MODE=1 (Pass 1: dfSelected only; Pass 2: reproc for fid/QADB/corr)
      --no-two-pass           Set TWO_PASS_MODE=0 (original behaviour)
      --pass2-threads N       Threads for the Pass 2 reproc job (default: 8)
      --success-grep TEXT     Override SUCCESS_GREP
      --mc / --no-mc
      --reproc / --no-reproc
      --inbending / --outbending
      --minimal / --no-minimal
      --missingKm / --no-missingKm
      --hphm / --no-hphm
      --reproc-file NAME
      --reproc-tree NAME
  -h, --help                  Show this help
USAGE
}

# =============================================================================
# USER SETTINGS (TOP)
# =============================================================================

# ---- Parallelization / splitting
K=120

# ---- Executable for skimming
EXE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/build/AnalysisPhi

# ---- Input data
#INPUT_DIR=/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/
#INPUT_DIR=/lustre24/expphy/volatile/clas12/singh/hipo2root/osg/sims_lager/fall2018_outb_50nA_phi_lager/
#INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/avakian/10388/
#INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/yijie/10386/
#INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/yijie/10395/
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/kneupane/10424/ #sp18_outb_50nA
#INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/singh/10423/ #fall_outb_50nA
#INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/singh/10422/ #fall_inb_55nA

#data for hplushminus (Bhawani)
INPUT_DIR=/lustre24/expphy/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/





### efficiencies for phi analysis (Bhawani)
INPUT_DIR=/lustre24/expphy/volatile/clas12/singh/hipo2root/osg/sims_lager/fall2018_outb_50nA_phi_lager/
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/singh/10527/
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/singh/10527/
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/kneupane/10646/ #sp19_inb_50nA
#spring 2018 efficincies 
INPUT_DIR=//lustre24/expphy/volatile/clas12/osg/yijie/10503/ #sp18_inb_50nA
INPUT_DIR=//lustre24/expphy/volatile/clas12/osg/yijie/10504/ #sp18_inb2_50nA
INPUT_DIR=//lustre24/expphy/volatile/clas12/osg/yijie/10505/ #sp18_inb_no_bkg
INPUT_DIR=//lustre24/expphy/volatile/clas12/osg/yijie/10506/ #sp18_inb_no_bkg

INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/yijie/10507/ #sp18_outb_45nA
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/yijie/10508/ #sp18_outb2_45nA
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/yijie/10509/ #sp18_outb_no_bkg
INPUT_DIR=/lustre24/expphy/volatile/clas12/osg/yijie/10510/ #sp18_outb_no_bkg

#data for hplushminus (Bhawani)
INPUT_DIR=/lustre24/expphy/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/


# ---- Output base (ABSOLUTE OR RELATIVE OK)
# Everything will be created under this directory:
#   OUTPUT_BASE/<WORKDIR>/job_###/{input,out,log}
#   OUTPUT_BASE/<merged outputs...>
#OUTPUT_BASE="${OUTPUT_BASE:-$PWD}"     # override with env var if you want
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/clasdis/sp2019_50nA_inb/
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2019/50nA/
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2018_outb/45nA/
#OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/fall2018_outb/50nA/
#OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/fall2018_inb/55nA/

#efficiencies for phi analysis (Bhawani)
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/fall2018_inb/55nA/

#Harut MC GEN test:
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/Harut_MC/p39_10p2/
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2019_inb/50nA
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2019_inb/no_bkg/

# spring 2018 efficiencies:
OUTPUT_BASE=//w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_inb/50nA/
OUTPUT_BASE=//w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_inb2/50nA/
OUTPUT_BASE=//w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_inb/no_bkg/
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_inb2/no_bkg/
OUTPUT_BASE=//w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_outb2/45nA/
OUTPUT_BASE=//w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_outb2/no_bkg/


OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DSTs/hplus_hminus_ep/
OUTPUT_BASE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/nSIDIS/hphm/rgasp19_inb/

# ---- AnalysisPhi args (matches AnalysisDVCS.cpp)
NEVT=-1                # AnalysisPhi: -n / --nfile
THREADS_PER_JOB=1      # AnalysisPhi: -t / --threads
CFG=rgafall18_inb     # AnalysisPhi: -c / --config
CFG=rgasp18_outb        # AnalysisPhi: -c / --config
CFG=rgasp19_inb  

# Boolean flags (0/1) -> translated to presence/absence of flags
IS_MC=0              # --mc
IS_REPROC=0            # --reproc
IS_INBENDING=1         # --inbending
IS_MINIMAL=0           # --minimal
IS_MISSINGKM=0         # --missingKm
IS_HPHM=0              # --hphm  (ep -> e'p'h+h- di-charged-hadron mode)

# Reproc options (only used if IS_REPROC=1; passed if non-empty)
REPROC_FILE="merged.root"         # --reproc-file <filename>
REPROC_TREE="hipo"                # --reproc-tree <treename>

# ---- Expected outputs per job (we will merge each of these independently)
typeset -a OUTPUT_FILES
OUTPUT_FILES=(
  dfSelected_afterFid_afterCorr.root
  dfSelected_afterFid.root
  dfSelectedMC.root
  dfSelected.root
)

# Keep OUTROOT for informational prints (per-job default output name)
OUTROOT=dfSelected.root

# RELAUNCH_ONLY=0 : create new WORKDIR, distribute files, run jobs
# RELAUNCH_ONLY=1 : use existing WORKDIR, relaunch failed jobs only
RELAUNCH_ONLY=0

# When RELAUNCH_ONLY=0, WORKDIR auto-created under OUTPUT_BASE:
WORKDIR="disana_jobs_$(date +%Y%m%d_%H%M%S)"
# When RELAUNCH_ONLY=1, set WORKDIR to an existing one (name only or full path):
#WORKDIR="disana_jobs_20260204_194543"

# Merge behavior
MERGE_PARTIAL=1   # 1 merge what succeeded, 0 require all success

# Success detection (log grep)
SUCCESS_GREP='Finished processing all events'

# ---- Manual override (switchable) ----
# If 0: use auto-detection
# If 1: relaunch ONLY the job IDs listed in FAILED_JOBS
USE_MANUAL_FAILED_LIST=0
FAILED_JOBS=(000 008 009 018 024 029)   # must match ID width
# -------------------------------------

# ---- Batching ----
RUN_IN_BATCHES=1
JOBS_PER_BATCH=20
# -------------------------------------

# ---- Two-pass mode ----
# TWO_PASS_MODE=1 : Pass 1 writes only dfSelected.root per job (--pass1 flag,
#                   skips fiducial/QADB/correction).  After merging, a single
#                   Pass 2 reproc job reads the merged dfSelected.root and
#                   produces dfSelected_afterFid_reprocessed.root +
#                   dfSelected_afterFid_afterCorr.root directly in OUTPUT_BASE.
# TWO_PASS_MODE=0 : Original behaviour (all outputs per job, all merged).
TWO_PASS_MODE=0
PASS2_THREADS=16          # threads for the single-job Pass 2 reproc step
PASS2_REPROC_TREE="dfSelected"   # tree name written by Snapshot in Pass 1
# -------------------------------------

# =============================================================================
# ARGUMENT PARSING (overrides old defaults above)
# =============================================================================

while (( $# > 0 )); do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -i|--input-dir)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      INPUT_DIR="$2"
      shift 2
      ;;
    -o|--output-base)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      OUTPUT_BASE="$2"
      shift 2
      ;;
    -c|--config)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      CFG="$2"
      shift 2
      ;;
    -t|--threads)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      THREADS_PER_JOB="$2"
      shift 2
      ;;
    -n|--nfile)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      NEVT="$2"
      shift 2
      ;;
    -k|--jobs)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      K="$2"
      shift 2
      ;;
    --exe)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      EXE="$2"
      shift 2
      ;;
    --workdir)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      WORKDIR="$2"
      shift 2
      ;;
    --relaunch-only)
      RELAUNCH_ONLY=1
      shift
      ;;
    --merge-partial)
      MERGE_PARTIAL=1
      shift
      ;;
    --no-merge-partial)
      MERGE_PARTIAL=0
      shift
      ;;
    --run-in-batches)
      RUN_IN_BATCHES=1
      shift
      ;;
    --no-run-in-batches)
      RUN_IN_BATCHES=0
      shift
      ;;
    --two-pass)
      TWO_PASS_MODE=1
      shift
      ;;
    --no-two-pass)
      TWO_PASS_MODE=0
      shift
      ;;
    --pass2-threads)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      PASS2_THREADS="$2"
      shift 2
      ;;
    --jobs-per-batch)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      JOBS_PER_BATCH="$2"
      shift 2
      ;;
    --success-grep)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      SUCCESS_GREP="$2"
      shift 2
      ;;
    --mc)
      IS_MC=1
      shift
      ;;
    --no-mc)
      IS_MC=0
      shift
      ;;
    --reproc)
      IS_REPROC=1
      shift
      ;;
    --no-reproc)
      IS_REPROC=0
      shift
      ;;
    --inbending)
      IS_INBENDING=1
      shift
      ;;
    --outbending)
      IS_INBENDING=0
      shift
      ;;
    --minimal)
      IS_MINIMAL=1
      shift
      ;;
    --no-minimal)
      IS_MINIMAL=0
      shift
      ;;
    --missingKm)
      IS_MISSINGKM=1
      shift
      ;;
    --no-missingKm)
      IS_MISSINGKM=0
      shift
      ;;
    --hphm)
      IS_HPHM=1
      shift
      ;;
    --no-hphm)
      IS_HPHM=0
      shift
      ;;
    --reproc-file)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      REPROC_FILE="$2"
      shift 2
      ;;
    --reproc-tree)
      [[ $# -ge 2 ]] || { echo "[ERROR] Missing value for $1"; exit 2; }
      REPROC_TREE="$2"
      shift 2
      ;;
    *)
      echo "[ERROR] Unknown argument: $1"
      echo
      usage
      exit 2
      ;;
  esac
done

# =============================================================================
# Derived / normalized paths
# =============================================================================

# Normalize OUTPUT_BASE to absolute path (best effort)
# If OUTPUT_BASE is relative, make it relative to where you run the script from.
mkdir -p "$OUTPUT_BASE" 2>/dev/null || true
OUTPUT_BASE_ABS="$(cd "$OUTPUT_BASE" 2>/dev/null && pwd)"
if [[ -z "${OUTPUT_BASE_ABS}" ]]; then
  echo "[ERROR] OUTPUT_BASE does not exist or not accessible: $OUTPUT_BASE"
  exit 2
fi

# WORKDIR can be provided as name or full path; ensure final is under OUTPUT_BASE_ABS
if [[ "$WORKDIR" = /* ]]; then
  WORKDIR_ABS="$WORKDIR"
else
  WORKDIR_ABS="$OUTPUT_BASE_ABS/$WORKDIR"
fi

# Where "default" merged file would be (not the only one merged!)
FINAL_MERGED="$OUTPUT_BASE_ABS/$OUTROOT"

# Build AnalysisPhi EXTRA_ARGS array based on booleans / strings
typeset -a EXTRA_ARGS
EXTRA_ARGS=()
(( IS_MC ))        && EXTRA_ARGS+=(--mc)
(( IS_REPROC ))    && EXTRA_ARGS+=(--reproc)
(( IS_INBENDING )) && EXTRA_ARGS+=(--inbending)
(( IS_MINIMAL ))   && EXTRA_ARGS+=(--minimal)
(( IS_MISSINGKM )) && EXTRA_ARGS+=(--missingKm)
(( IS_HPHM ))      && EXTRA_ARGS+=(--hphm)
# Only inject reproc-file / reproc-tree when IS_REPROC=1 — these are
# meaningless (and polluting) for normal hipo-reading jobs.
(( IS_REPROC )) && [[ -n "$REPROC_FILE" ]] && EXTRA_ARGS+=(--reproc-file "$REPROC_FILE")
(( IS_REPROC )) && [[ -n "$REPROC_TREE" ]] && EXTRA_ARGS+=(--reproc-tree "$REPROC_TREE")

# =============================================================================
# Print configuration
# =============================================================================
echo "[INFO] K                     = $K"
echo "[INFO] EXE                   = $EXE"
echo "[INFO] INPUT_DIR             = $INPUT_DIR"
echo "[INFO] OUTPUT_BASE_ABS       = $OUTPUT_BASE_ABS"
echo "[INFO] WORKDIR_ABS           = $WORKDIR_ABS"
echo "[INFO] CFG                   = $CFG"
echo "[INFO] NEVT                  = $NEVT"
echo "[INFO] THREADS_PER_JOB       = $THREADS_PER_JOB"
echo "[INFO] EXTRA_ARGS            = ${EXTRA_ARGS[*]}"
echo "[INFO] Expected outputs/job  = ${OUTPUT_FILES[*]}"
echo "[INFO] RELAUNCH_ONLY         = $RELAUNCH_ONLY"
echo "[INFO] MERGE_PARTIAL         = $MERGE_PARTIAL"
echo "[INFO] USE_MANUAL_FAILED_LIST= $USE_MANUAL_FAILED_LIST"
echo "[INFO] FAILED_JOBS           = ${FAILED_JOBS[*]}"
echo "[INFO] RUN_IN_BATCHES        = $RUN_IN_BATCHES"
echo "[INFO] JOBS_PER_BATCH (J)    = $JOBS_PER_BATCH"
echo "[INFO] TWO_PASS_MODE         = $TWO_PASS_MODE"
(( TWO_PASS_MODE )) && echo "[INFO] PASS2_THREADS         = $PASS2_THREADS"
echo ""

# =============================================================================
# Sanity checks
# =============================================================================
if [[ ! -x "$EXE" ]]; then
  echo "[ERROR] Executable not found or not executable: $EXE"
  exit 2
fi

if (( RELAUNCH_ONLY == 0 )); then
  if [[ ! -d "$INPUT_DIR" ]]; then
    echo "[ERROR] Input directory not found: $INPUT_DIR"
    exit 2
  fi
else
  if [[ ! -d "$WORKDIR_ABS" ]]; then
    echo "[ERROR] WORKDIR not found for relaunch: $WORKDIR_ABS"
    exit 2
  fi
fi

if ! command -v hadd >/dev/null 2>&1; then
  echo "[ERROR] 'hadd' not found in PATH. Source your ROOT environment."
  exit 2
fi

# =============================================================================
# Helpers
# =============================================================================
fmt_id () { printf "%03d" "$1"; }

job_paths () {
  local kid="$1"
  local JOBDIR="$WORKDIR_ABS/job_${kid}"
  local INDIR="$JOBDIR/input"
  local OUTDIR="$JOBDIR/out"
  local LOG="$JOBDIR/job_${kid}.log"
  echo "$JOBDIR" "$INDIR" "$OUTDIR" "$LOG"
}

# A job is OK if:
#  - log exists and contains SUCCESS_GREP
#  - and at least ONE expected output ROOT exists and is non-empty
job_ok () {
  local kid="$1"
  local JOBDIR INDIR OUTDIR LOG
  read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$kid")"

  [[ -f "$LOG" ]] || return 1
  grep -Fq "$SUCCESS_GREP" "$LOG" 2>/dev/null || return 1

  local f
  for f in "${OUTPUT_FILES[@]}"; do
    [[ -s "$OUTDIR/$f" ]] && return 0
  done
  return 1
}

run_job_bg () {
  local kid="$1"
  local JOBDIR INDIR OUTDIR LOG
  read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$kid")"
  mkdir -p "$INDIR" "$OUTDIR"

  # Remove old outputs for this job (all possible)
  local f
  for f in "${OUTPUT_FILES[@]}"; do
    rm -f "$OUTDIR/$f" 2>/dev/null
  done

  "$EXE" \
    -i "$INDIR" \
    -n "$NEVT" \
    -t "$THREADS_PER_JOB" \
    -o "$OUTDIR" \
    -c "$CFG" \
    "${EXTRA_ARGS[@]}" \
    &> "$LOG" &
}

launch_jobs_in_batches () {
  local -a JOBLIST
  JOBLIST=("$@")
  local total=${#JOBLIST[@]}

  if (( total == 0 )); then
    echo "[INFO] No jobs to launch."
    return 0
  fi

  if (( RUN_IN_BATCHES == 0 )); then
    echo "[INFO] Launching $total job(s) all at once..."
    local kid
    for kid in "${JOBLIST[@]}"; do
      run_job_bg "$kid"
    done
    wait
    return 0
  fi

  # Sliding-window pool: keep at most JOBS_PER_BATCH jobs running at once.
  # As soon as any job finishes a slot opens and the next job is dispatched
  # immediately — no waiting for a whole batch to drain.
  echo "[INFO] Launching $total job(s) with sliding window (max $JOBS_PER_BATCH concurrent)..."

  local -a pool_pids=()     # PIDs currently running
  local -a still_running=() # reused each poll cycle — declared here to avoid
  local _pid                #   zsh re-declaring (and echoing) inside the loop
  local reaped=0
  local idx=1

  while (( idx <= total || ${#pool_pids[@]} > 0 )); do

    # ---- Fill empty slots ------------------------------------------------
    while (( idx <= total && ${#pool_pids[@]} < JOBS_PER_BATCH )); do
      run_job_bg "${JOBLIST[idx]}"
      pool_pids+=($!)
      echo "[INFO] Started job ${JOBLIST[idx]} ($idx/$total), active=${#pool_pids[@]}"
      (( idx++ ))
    done

    # ---- Reap any finished PIDs -----------------------------------------
    # Poll each PID: kill -0 returns non-zero once the process has exited.
    still_running=()
    reaped=0
    for _pid in "${pool_pids[@]}"; do
      if kill -0 "$_pid" 2>/dev/null; then
        still_running+=("$_pid")
      else
        wait "$_pid" 2>/dev/null   # collect exit status / avoid zombie
        (( reaped++ ))
      fi
    done
    pool_pids=("${still_running[@]}")

    # If nothing finished yet and we are at capacity, pause briefly before
    # polling again so we do not spin the CPU.
    if (( reaped == 0 && ${#pool_pids[@]} >= JOBS_PER_BATCH )); then
      sleep 0.5
    fi

  done

  echo "[INFO] All $total job(s) complete."
}

# =============================================================================
# Mode A: Full run
# =============================================================================
if (( RELAUNCH_ONLY == 0 )); then
  mkdir -p "$WORKDIR_ABS/lists"

  ALLLIST="$WORKDIR_ABS/lists/all_files.txt"
  find "$INPUT_DIR" -type f -name "*.hipo" | sort > "$ALLLIST"

  NF=$(wc -l < "$ALLLIST")
  if (( NF == 0 )); then
    echo "[ERROR] No *.hipo files found under: $INPUT_DIR"
    exit 3
  fi
  echo "[INFO] Found $NF input files."
  echo ""

  # Create empty lists
  integer i
  for (( i=0; i<K; i++ )); do
    kid="$(fmt_id $i)"
    : > "$WORKDIR_ABS/lists/files_${kid}.txt"
  done

  # Round-robin distribute
  integer j=0
  while IFS= read -r f; do
    kid="$(fmt_id $j)"
    echo "$f" >> "$WORKDIR_ABS/lists/files_${kid}.txt"
    (( j++ ))
    (( j >= K )) && j=0
  done < "$ALLLIST"

  # Prepare job dirs + symlinks
  for (( i=0; i<K; i++ )); do
    kid="$(fmt_id $i)"
    LIST="$WORKDIR_ABS/lists/files_${kid}.txt"
    [[ -s "$LIST" ]] || continue

    read -r _JOBDIR _INDIR _OUTDIR _LOG <<<"$(job_paths "$kid")"
    mkdir -p "$_INDIR" "$_OUTDIR"
    rm -f "$_INDIR"/*.hipo 2>/dev/null
    while IFS= read -r ff; do
      ln -sf "$ff" "$_INDIR/"
    done < "$LIST"
  done

  # Build job list
  typeset -a JOBS_TO_RUN
  JOBS_TO_RUN=()
  for (( i=0; i<K; i++ )); do
    kid="$(fmt_id $i)"
    [[ -s "$WORKDIR_ABS/lists/files_${kid}.txt" ]] && JOBS_TO_RUN+=("$kid")
  done

  echo "[INFO] Launching ${#JOBS_TO_RUN[@]} job(s) total."

  # In TWO_PASS_MODE Pass 1: inject --pass1 so each job skips fid/QADB/correction
  # and only writes dfSelected.root.  This saves both per-job CPU and per-job disk.
  if (( TWO_PASS_MODE )); then
    echo "[INFO] TWO_PASS_MODE: injecting --pass1 into all jobs (fid/QADB/correction skipped per job)."
    EXTRA_ARGS+=(--pass1)
  fi

  launch_jobs_in_batches "${JOBS_TO_RUN[@]}"
fi

# =============================================================================
# Mode B: Relaunch only (existing WORKDIR)
# =============================================================================
typeset -a JOBIDS FAILED
JOBIDS=()
FAILED=()

if (( RELAUNCH_ONLY == 1 )); then
  echo "[INFO] Relaunch mode: scanning $WORKDIR_ABS"
  echo ""

  _d=""
  for _d in "$WORKDIR_ABS"/job_[0-9][0-9][0-9]; do
    [[ -d "$_d" ]] || continue
    kid="${_d:t}"
    kid="${kid#job_}"
    JOBIDS+=("$kid")
  done

  if (( ${#JOBIDS[@]} == 0 )); then
    echo "[ERROR] No job_### directories found under $WORKDIR_ABS"
    exit 3
  fi

  if (( USE_MANUAL_FAILED_LIST == 1 )); then
    FAILED=("${FAILED_JOBS[@]}")
    echo "[INFO] Using MANUAL failed list: ${FAILED[*]}"
  else
    for kid in "${JOBIDS[@]}"; do
      if job_ok "$kid"; then
        echo "[INFO] Job $kid: OK (skip)"
      else
        echo "[WARN] Job $kid: FAILED (relaunch)"
        FAILED+=("$kid")
      fi
    done
  fi

  echo ""
  echo "[INFO] Relaunching ${#FAILED[@]} job(s)."
  launch_jobs_in_batches "${FAILED[@]}"
fi

# =============================================================================
# Summary + merge (multi-output)
# =============================================================================
# In TWO_PASS_MODE Pass 1: only dfSelected.root matters per job.
# Override the file list used for success-checking and merging.
typeset -a MERGE_OUTPUT_FILES
if (( TWO_PASS_MODE )); then
  MERGE_OUTPUT_FILES=(dfSelected.root)
  echo "[INFO] TWO_PASS_MODE: Pass 1 merge restricted to: ${MERGE_OUTPUT_FILES[*]}"
else
  MERGE_OUTPUT_FILES=("${OUTPUT_FILES[@]}")
fi

typeset -a CHECKSET OKJ BADJ
OKJ=()
BADJ=()

# For collecting per-output-file lists and availability stats
typeset -A MERGE_LISTS   # MERGE_LISTS["dfSelected.root"] -> " file1 file2 ..."
typeset -A MISS_COUNTS   # missing per output among OK jobs
typeset -A HAVE_COUNTS   # present per output among OK jobs

# init — plain variables, no local (we are at script top level, not in a function)
_f=""
for _f in "${MERGE_OUTPUT_FILES[@]}"; do
  MERGE_LISTS[$_f]=""
  MISS_COUNTS[$_f]=0
  HAVE_COUNTS[$_f]=0
done

if (( RELAUNCH_ONLY == 1 )); then
  CHECKSET=("${JOBIDS[@]}")
else
  CHECKSET=()
  for (( i=0; i<K; i++ )); do
    kid="$(fmt_id $i)"
    [[ -s "$WORKDIR_ABS/lists/files_${kid}.txt" ]] && CHECKSET+=("$kid")
  done
fi

for kid in "${CHECKSET[@]}"; do
  if job_ok "$kid"; then
    OKJ+=("$kid")

    read -r _JOBDIR _INDIR _OUTDIR _LOG <<<"$(job_paths "$kid")"

    for _f in "${MERGE_OUTPUT_FILES[@]}"; do
      if [[ -s "$_OUTDIR/$_f" ]]; then
        MERGE_LISTS[$_f]="${MERGE_LISTS[$_f]} $_OUTDIR/$_f"
        HAVE_COUNTS[$_f]=$(( HAVE_COUNTS[$_f] + 1 ))
      else
        MISS_COUNTS[$_f]=$(( MISS_COUNTS[$_f] + 1 ))
      fi
    done
  else
    BADJ+=("$kid")
  fi
done

echo "[INFO] Summary:"
echo "       Success: ${#OKJ[@]} job(s)"
echo "       Failed : ${#BADJ[@]} job(s) : ${BADJ[*]}"
echo ""

if (( ${#BADJ[@]} > 0 && MERGE_PARTIAL == 0 )); then
  echo "[ERROR] Not merging (MERGE_PARTIAL=0 and some jobs failed)"
  exit 6
fi

# Merge each output file type independently if we have at least one input
merged_any=0
_inputs=()
_outmerged=""
for _f in "${MERGE_OUTPUT_FILES[@]}"; do
  _inputs=(${=MERGE_LISTS[$_f]})

  if (( ${#_inputs[@]} > 0 )); then
    _outmerged="$OUTPUT_BASE_ABS/$_f"
    echo "[INFO] Merging ${#_inputs[@]} file(s) into: $_outmerged"
    hadd -f -k "$_outmerged" "${_inputs[@]}"
    HSTAT=$?
    if (( HSTAT != 0 )); then
      echo "[ERROR] hadd failed for $_f with status $HSTAT"
      exit 5
    fi
    merged_any=1
  else
    echo "[WARN] No inputs found for $_f -> skipping merge."
  fi
done

echo ""
echo "[INFO] Per-output availability across OK jobs:"
for _f in "${MERGE_OUTPUT_FILES[@]}"; do
  echo "       $_f : present=${HAVE_COUNTS[$_f]} missing=${MISS_COUNTS[$_f]}"
done
echo ""

if (( merged_any == 0 )); then
  echo "[ERROR] No output files found to merge for any of: ${OUTPUT_FILES[*]}"
  exit 4
fi

echo "[INFO] Done."
echo "[INFO] Merged outputs are under: $OUTPUT_BASE_ABS/"
echo "[INFO] Per-job outputs:         $WORKDIR_ABS/job_*/out/"

# =============================================================================
# TWO-PASS MODE: Pass 2 — reproc the merged dfSelected.root
# =============================================================================
if (( TWO_PASS_MODE && merged_any )); then
  PASS2_INPUT="$OUTPUT_BASE_ABS/dfSelected.root"
  if [[ ! -s "$PASS2_INPUT" ]]; then
    echo "[ERROR] TWO_PASS_MODE=1 but merged dfSelected.root not found: $PASS2_INPUT"
    exit 7
  fi

  PASS2_LOG="$OUTPUT_BASE_ABS/pass2_reproc.log"
  echo ""
  echo "[INFO] ============================================================"
  echo "[INFO] TWO-PASS MODE: starting Pass 2 (single reproc job)"
  echo "[INFO]   Input  : $PASS2_INPUT"
  echo "[INFO]   Tree   : $PASS2_REPROC_TREE"
  echo "[INFO]   Output : $OUTPUT_BASE_ABS"
  echo "[INFO]   Threads: $PASS2_THREADS"
  echo "[INFO]   Log    : $PASS2_LOG"
  echo "[INFO] ============================================================"

  # Remove --pass1 from EXTRA_ARGS if it was added — Pass 2 must run full cuts.
  typeset -a EXTRA_ARGS_PASS2
  EXTRA_ARGS_PASS2=()
  _arg=""
  for _arg in "${EXTRA_ARGS[@]}"; do
    [[ "$_arg" == "--pass1" ]] && continue
    EXTRA_ARGS_PASS2+=("$_arg")
  done

  "$EXE" \
    -i "$OUTPUT_BASE_ABS" \
    -n "$NEVT" \
    -t "$PASS2_THREADS" \
    -o "$OUTPUT_BASE_ABS" \
    -c "$CFG" \
    --reproc \
    --reproc-file "dfSelected.root" \
    --reproc-tree "$PASS2_REPROC_TREE" \
    "${EXTRA_ARGS_PASS2[@]}" \
    > "$PASS2_LOG" 2>&1

  P2STAT=$?
  if (( P2STAT == 0 )); then
    echo "[INFO] Pass 2 complete."
    echo "[INFO] After-fid / after-correction outputs written to: $OUTPUT_BASE_ABS"
  else
    echo "[ERROR] Pass 2 failed (exit $P2STAT). See log: $PASS2_LOG"
    exit 7
  fi
fi

exit 0
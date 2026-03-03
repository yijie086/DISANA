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

# =============================================================================
# USER SETTINGS (TOP)
# =============================================================================

# ---- Parallelization / splitting
K=200

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
# ---- AnalysisPhi args (matches AnalysisDVCS.cpp)
NEVT=-1                # AnalysisPhi: -n / --nfile
THREADS_PER_JOB=1      # AnalysisPhi: -t / --threads
CFG=rgafall18_outb     # AnalysisPhi: -c / --config
CFG=rgasp18_outb        # AnalysisPhi: -c / --config
#CFG=rgafall18_inb  

# Boolean flags (0/1) -> translated to presence/absence of flags
IS_MC=1                # --mc
IS_REPROC=0            # --reproc
IS_INBENDING=1         # --inbending
IS_MINIMAL=0           # --minimal
IS_MISSINGKM=0         # --missingKm

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
JOBS_PER_BATCH=40
# -------------------------------------

# =============================================================================
# Derived / normalized paths
# =============================================================================

# Normalize OUTPUT_BASE to absolute path (best effort)
# If OUTPUT_BASE is relative, make it relative to where you run the script from.
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
[[ -n "$REPROC_FILE" ]] && EXTRA_ARGS+=(--reproc-file "$REPROC_FILE")
[[ -n "$REPROC_TREE" ]] && EXTRA_ARGS+=(--reproc-tree "$REPROC_TREE")

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

  local b=1
  local start=1
  while (( start <= total )); do
    local end=$(( start + JOBS_PER_BATCH - 1 ))
    (( end > total )) && end=$total

    echo "[INFO] Batch $b: launching jobs index $start..$end (count=$(( end-start+1 )))"
    local i
    for (( i=start; i<=end; i++ )); do
      run_job_bg "${JOBLIST[i]}"
    done

    echo "[INFO] Batch $b: waiting..."
    wait
    echo "[INFO] Batch $b: done."
    echo ""
    (( b++ ))
    start=$(( end + 1 ))
  done
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

    local JOBDIR INDIR OUTDIR LOG
    read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$kid")"
    mkdir -p "$INDIR" "$OUTDIR"

    rm -f "$INDIR"/*.hipo 2>/dev/null
    while IFS= read -r ff; do
      ln -sf "$ff" "$INDIR/"
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

  local d
  for d in "$WORKDIR_ABS"/job_[0-9][0-9][0-9]; do
    [[ -d "$d" ]] || continue
    kid="${d:t}"
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
    local kid
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
typeset -a CHECKSET OKJ BADJ
OKJ=()
BADJ=()

# For collecting per-output-file lists and availability stats
typeset -A MERGE_LISTS   # MERGE_LISTS["dfSelected.root"] -> " file1 file2 ..."
typeset -A MISS_COUNTS   # missing per output among OK jobs
typeset -A HAVE_COUNTS   # present per output among OK jobs

# init
local f
for f in "${OUTPUT_FILES[@]}"; do
  MERGE_LISTS[$f]=""
  MISS_COUNTS[$f]=0
  HAVE_COUNTS[$f]=0
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

local kid
for kid in "${CHECKSET[@]}"; do
  if job_ok "$kid"; then
    OKJ+=("$kid")

    local JOBDIR INDIR OUTDIR LOG
    read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$kid")"

    for f in "${OUTPUT_FILES[@]}"; do
      if [[ -s "$OUTDIR/$f" ]]; then
        MERGE_LISTS[$f]="${MERGE_LISTS[$f]} $OUTDIR/$f"
        HAVE_COUNTS[$f]=$(( HAVE_COUNTS[$f] + 1 ))
      else
        MISS_COUNTS[$f]=$(( MISS_COUNTS[$f] + 1 ))
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
local merged_any=0
for f in "${OUTPUT_FILES[@]}"; do
  local inputs=(${=MERGE_LISTS[$f]})

  if (( ${#inputs[@]} > 0 )); then
    local outmerged="$OUTPUT_BASE_ABS/$f"
    echo "[INFO] Merging ${#inputs[@]} file(s) into: $outmerged"
    hadd -f -k "$outmerged" "${inputs[@]}"
    HSTAT=$?
    if (( HSTAT != 0 )); then
      echo "[ERROR] hadd failed for $f with status $HSTAT"
      exit 5
    fi
    merged_any=1
  else
    echo "[WARN] No inputs found for $f -> skipping merge."
  fi
done

echo ""
echo "[INFO] Per-output availability across OK jobs:"
for f in "${OUTPUT_FILES[@]}"; do
  echo "       $f : present=${HAVE_COUNTS[$f]} missing=${MISS_COUNTS[$f]}"
done
echo ""

if (( merged_any == 0 )); then
  echo "[ERROR] No output files found to merge for any of: ${OUTPUT_FILES[*]}"
  exit 4
fi

echo "[INFO] Done."
echo "[INFO] Merged outputs are under: $OUTPUT_BASE_ABS/"
echo "[INFO] Per-job outputs:         $WORKDIR_ABS/job_*/out/"
exit 0

#!/bin/zsh
# run_disana_all_datasets.zsh
#
# Orchestrator for DISANA AnalysisPhi across all RG-A datasets.
#
# Structure:
#   1 main job
#   └── 5 parallel sub-jobs (one per dataset: rgasp18_inb/outb, rgafall18_inb/outb, rgasp19_inb)
#       └── K parallel AnalysisPhi jobs per sub-job, run in sliding-window batches → merge
#
# Single required flag:
#   --DVKpKm    Run the DVKpKmP wagon (ep → e'K+K-p)
#   --nSIDIS    Run the nSidis wagon
#
# All other parameters are defined in the USER CONFIGURATION section below.
#
# Usage:
#   ./run_disana_all_datasets.zsh --DVKpKm
#   ./run_disana_all_datasets.zsh --nSIDIS

set -u
set -o pipefail
setopt NO_NOMATCH   # suppress zsh "no matches found" on empty globs

# =============================================================================
# USER CONFIGURATION — edit this block; no other changes needed
# =============================================================================

# ---- Executable
EXE=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/build/AnalysisPhi

# ---- Top-level output directory; mode sub-dir is appended automatically:
#        $OUTPUT_ROOT/DVKpKm/<dataset>/  or  $OUTPUT_ROOT/nSIDIS/<dataset>/
OUTPUT_ROOT=/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed

# ---- Parallelization
K=60                # parallel AnalysisPhi jobs per dataset
JOBS_PER_BATCH=30   # sliding-window concurrency cap (jobs running at one time)
RUN_IN_BATCHES=1    # 1 = sliding window (recommended); 0 = all K at once

# ---- AnalysisPhi job parameters
NEVT=-1             # -n / --nfile  (-1 = all events)
THREADS_PER_JOB=1   # -t / --threads

# ---- Merge / cleanup behaviour
MERGE_PARTIAL=1         # 1 = merge whatever succeeded; 0 = require all jobs OK
CLEANUP_JOB_OUTPUTS=1   # 1 = delete per-job ROOT files after successful merge

# ---- Success detection (searched in each job's log file)
SUCCESS_GREP='Finished processing all events'

# ---- Two-pass mode
TWO_PASS_MODE=0           # 0 = original (all outputs per job); 1 = pass1 → merge → pass2
PASS2_THREADS=16           # threads for the single pass-2 reproc step
PASS2_REPROC_TREE="dfSelected"

# ---- Expected ROOT outputs per job (each file type is merged independently)
typeset -a OUTPUT_FILES
OUTPUT_FILES=(
  dfSelected_afterFid_afterCorr.root
  dfSelected_afterFid.root
  dfSelectedMC.root
  dfSelected.root
)

# ---- Mode-specific AnalysisPhi boolean flags (0/1)
#      These are set per-mode below; override here only if both modes share a value.

# DVKpKm mode
DVKpKm_IS_MC=0
DVKpKm_IS_REPROC=0
DVKpKm_IS_MINIMAL=0
DVKpKm_IS_MISSINGKM=0
DVKpKm_IS_HPHM=1        # ep → e'K+K- di-hadron mode

# nSIDIS mode
nSIDIS_IS_MC=0
nSIDIS_IS_REPROC=0
nSIDIS_IS_MINIMAL=0
nSIDIS_IS_MISSINGKM=0
nSIDIS_IS_HPHM=0

# =============================================================================
# ARGUMENT PARSING — single flag only
# =============================================================================

usage() {
  echo "Usage: $0 --DVKpKm | --nSIDIS"
  echo "  --DVKpKm   Process the DVKpKmP wagon across all RG-A datasets"
  echo "  --nSIDIS   Process the nSidis wagon across all RG-A datasets"
}

if (( $# != 1 )); then
  usage; exit 2
fi

case "$1" in
  --DVKpKm)
    MODE=DVKpKm
    WAGON=DVKpKmP
    IS_MC=$DVKpKm_IS_MC
    IS_REPROC=$DVKpKm_IS_REPROC
    IS_MINIMAL=$DVKpKm_IS_MINIMAL
    IS_MISSINGKM=$DVKpKm_IS_MISSINGKM
    IS_HPHM=$DVKpKm_IS_HPHM

    IN_rgasp18_inb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/DVKpKmP/"
    IN_rgasp18_outb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/DVKpKmP/"
    IN_rgafall18_inb="/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/DVKpKmP/"
    IN_rgafall18_outb="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/DVKpKmP/"
    IN_rgasp19_inb="/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/"
    ;;
  --nSIDIS)
    MODE=nSIDIS
    WAGON=nSidis
    IS_MC=$nSIDIS_IS_MC
    IS_REPROC=$nSIDIS_IS_REPROC
    IS_MINIMAL=$nSIDIS_IS_MINIMAL
    IS_MISSINGKM=$nSIDIS_IS_MISSINGKM
    IS_HPHM=$nSIDIS_IS_HPHM

    IN_rgasp18_inb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/nSidis/"
    IN_rgasp18_outb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/nSidis/"
    IN_rgafall18_inb="/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/"
    IN_rgafall18_outb="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/"
    IN_rgasp19_inb="/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/nSidis/"
    ;;
  -h|--help)
    usage; exit 0
    ;;
  *)
    echo "[ERROR] Unknown flag: $1"
    usage; exit 2
    ;;
esac

# =============================================================================
# DATASET TABLE
# Each row: label | input_dir | output_subdir | cfg_name | is_inbending
#
# Output layout under $OUTPUT_ROOT/$MODE/:
#   rgasp18/inb/    ← spring 2018 inbending
#   rgasp18/outb/   ← spring 2018 outbending
#   rgafall/inb/    ← fall 2018 inbending
#   rgafall/outb/   ← fall 2018 outbending
#   rgasp19/        ← spring 2019 inbending
# =============================================================================

typeset -a DS_LABELS DS_INPUTS DS_OUTDIRS DS_CFGS DS_INBENDING

DS_LABELS=(   rgasp18_inb         rgasp18_outb         rgafall18_inb         rgafall18_outb         rgasp19_inb    )
DS_INPUTS=(   "$IN_rgasp18_inb"   "$IN_rgasp18_outb"   "$IN_rgafall18_inb"   "$IN_rgafall18_outb"   "$IN_rgasp19_inb" )
DS_OUTDIRS=(  "rgasp18/inb"       "rgasp18/outb"       "rgafall/inb"         "rgafall/outb"         "rgasp19"      )
DS_CFGS=(     rgasp18_inb         rgasp18_outb          rgafall18_inb         rgafall18_outb         rgasp19_inb    )
DS_INBENDING=( 1                   0                    1                     0                      1              )

# =============================================================================
# SANITY CHECKS
# =============================================================================

if [[ ! -x "$EXE" ]]; then
  echo "[ERROR] Executable not found or not executable: $EXE"; exit 2
fi
if ! command -v hadd >/dev/null 2>&1; then
  echo "[ERROR] 'hadd' not found in PATH. Source your ROOT environment."; exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "[ERROR] 'python3' not found in PATH."; exit 2
fi
if ! python3 -c "import ROOT" 2>/dev/null; then
  echo "[ERROR] PyROOT not available. Source ROOT before running."; exit 2
fi

# =============================================================================
# HELPERS (shared across all sub-jobs via the subshell environment)
# =============================================================================

fmt_id() { printf "%03d" "$1"; }

# root_merge <output.root> <input1.root> [input2.root ...]
# Uses PyROOT/TFileMerger to bypass the 100 GB TTree size limit.
root_merge() {
  local _out="$1"; shift
  local -a _ins=("$@")
  local _script _add_lines=""
  _script="$(mktemp /tmp/disana_merge_XXXXXX.py)"
  for _f in "${_ins[@]}"; do
    _add_lines+="m.AddFile('${_f}')"$'\n'
  done
  cat > "$_script" <<PYEOF
import sys, ROOT
ROOT.gROOT.SetBatch(True)
ROOT.TTree.SetMaxTreeSize(int(1e15))
m = ROOT.TFileMerger(False)
m.SetFastMethod(True)
m.OutputFile('${_out}', 'RECREATE')
${_add_lines}
if not m.Merge():
    print('[root_merge] TFileMerger::Merge() failed', file=sys.stderr)
    sys.exit(1)
PYEOF
  python3 "$_script"
  local _stat=$?
  rm -f "$_script"
  return $_stat
}

# job_paths <workdir_abs> <kid>  → prints: JOBDIR INDIR OUTDIR LOG
job_paths() {
  local wdir="$1" kid="$2"
  local jd="$wdir/job_${kid}"
  echo "$jd" "$jd/input" "$jd/out" "$jd/job_${kid}.log"
}

# job_ok <workdir_abs> <kid>  → returns 0 if job succeeded, 1 otherwise
job_ok() {
  local wdir="$1" kid="$2"
  local JOBDIR INDIR OUTDIR LOG
  read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$wdir" "$kid")"
  [[ -f "$LOG" ]]                                     || return 1
  grep -Fq "$SUCCESS_GREP" "$LOG" 2>/dev/null         || return 1
  local f
  for f in "${OUTPUT_FILES[@]}"; do
    [[ -s "$OUTDIR/$f" ]] && return 0
  done
  return 1
}

# run_job_bg <wdir> <kid> <cfg> — starts one AnalysisPhi job in the background.
# Reads EXE, NEVT, THREADS_PER_JOB, EXTRA_ARGS_DS, OUTPUT_FILES from enclosing scope.
run_job_bg() {
  local wdir="$1" kid="$2" cfg="$3"
  local JOBDIR INDIR OUTDIR LOG
  read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$wdir" "$kid")"
  mkdir -p "$INDIR" "$OUTDIR"
  local f
  for f in "${OUTPUT_FILES[@]}"; do rm -f "$OUTDIR/$f" 2>/dev/null; done
  "$EXE" \
    -i "$INDIR" \
    -n "$NEVT" \
    -t "$THREADS_PER_JOB" \
    -o "$OUTDIR" \
    -c "$cfg" \
    "${EXTRA_ARGS_DS[@]}" \
    &> "$LOG" &
}

# launch_jobs — sliding-window dispatcher for one dataset's K jobs.
# Reads JOBS_DS, RUN_IN_BATCHES, JOBS_PER_BATCH, WORKDIR_DS, CFG_DS from enclosing scope.
launch_jobs() {
  local total=${#JOBS_DS[@]}
  (( total == 0 )) && { echo "  [INFO] No jobs to launch."; return 0; }

  if (( RUN_IN_BATCHES == 0 )); then
    echo "  [INFO] Launching $total job(s) all at once..."
    local kid
    for kid in "${JOBS_DS[@]}"; do run_job_bg "$WORKDIR_DS" "$kid" "$CFG_DS"; done
    wait
    return 0
  fi

  echo "  [INFO] Launching $total job(s) with sliding window (max $JOBS_PER_BATCH concurrent)..."
  local -a pool_pids=() still_running=()
  typeset -A pid_to_job
  local _pid _new_pid reaped idx=1 completed=0

  while (( idx <= total || ${#pool_pids[@]} > 0 )); do
    # Fill empty slots
    while (( idx <= total && ${#pool_pids[@]} < JOBS_PER_BATCH )); do
      run_job_bg "$WORKDIR_DS" "${JOBS_DS[$idx]}" "$CFG_DS"
      _new_pid=$!
      pool_pids+=($_new_pid)
      pid_to_job[$_new_pid]="${JOBS_DS[$idx]}"
      echo "  [INFO] Started   job ${JOBS_DS[$idx]} ($idx/$total), active=${#pool_pids[@]}"
      (( idx++ ))
    done
    # Reap finished PIDs
    still_running=(); reaped=0
    for _pid in "${pool_pids[@]}"; do
      if kill -0 "$_pid" 2>/dev/null; then
        still_running+=("$_pid")
      else
        wait "$_pid" 2>/dev/null
        (( completed++ )); (( reaped++ ))
        echo "  [INFO] Completed job ${pid_to_job[$_pid]} ($completed/$total done)"
        unset "pid_to_job[$_pid]"
      fi
    done
    pool_pids=("${still_running[@]}")
    # Avoid busy-wait when pool is full and nothing finished
    (( reaped == 0 && ${#pool_pids[@]} >= JOBS_PER_BATCH )) && sleep 0.5
  done
  echo "  [INFO] All $total job(s) complete."
}

# =============================================================================
# run_dataset — full pipeline for a single dataset (called inside a subshell).
#
# Arguments:
#   $1 label        — human-readable name (e.g. rgasp18_inb)
#   $2 input_dir    — path to hipo files
#   $3 output_base  — absolute output directory for this dataset
#   $4 cfg          — AnalysisPhi -c / --config value
#   $5 is_inbending — 1 = inbending, 0 = outbending
# =============================================================================
run_dataset() {
  local label="$1"
  local input_dir="$2"
  local output_base="$3"
  local CFG_DS="$4"
  local is_inbending="$5"
  local P="[$label]"   # log prefix

  # ---- Resolve absolute output path
  mkdir -p "$output_base"
  local OUTPUT_ABS
  OUTPUT_ABS="$(cd "$output_base" && pwd)"

  # ---- Unique workdir for this dataset's job staging area
  local WORKDIR_DS="$OUTPUT_ABS/disana_jobs_$(date +%Y%m%d_%H%M%S)"

  # ---- Build AnalysisPhi extra flags for this dataset
  local -a EXTRA_ARGS_DS
  EXTRA_ARGS_DS=()
  (( IS_MC ))        && EXTRA_ARGS_DS+=(--mc)
  (( IS_REPROC ))    && EXTRA_ARGS_DS+=(--reproc)
  (( is_inbending )) && EXTRA_ARGS_DS+=(--inbending)
  (( IS_MINIMAL ))   && EXTRA_ARGS_DS+=(--minimal)
  (( IS_MISSINGKM )) && EXTRA_ARGS_DS+=(--missingKm)
  (( IS_HPHM ))      && EXTRA_ARGS_DS+=(--hphm)

  echo "$P ============================================================"
  echo "$P Dataset  : $label"
  echo "$P Input    : $input_dir"
  echo "$P Output   : $OUTPUT_ABS"
  echo "$P Config   : $CFG_DS"
  echo "$P Inbending: $is_inbending"
  echo "$P ExtraArgs: ${EXTRA_ARGS_DS[*]}"
  echo "$P K jobs   : $K  |  batch size: $JOBS_PER_BATCH"
  echo "$P ============================================================"

  # ---- Validate input directory
  if [[ ! -d "$input_dir" ]]; then
    echo "$P [ERROR] Input directory not found: $input_dir" >&2
    return 1
  fi

  # ---- Discover and round-robin distribute hipo files across K job lists
  mkdir -p "$WORKDIR_DS/lists"
  local alllist="$WORKDIR_DS/lists/all_files.txt"
  find "$input_dir" -type f -name "*.hipo" | sort > "$alllist"

  local nf
  nf=$(wc -l < "$alllist")
  if (( nf == 0 )); then
    echo "$P [ERROR] No *.hipo files found under: $input_dir" >&2
    return 1
  fi
  echo "$P Found $nf input file(s)."

  integer i j=0
  for (( i=0; i<K; i++ )); do
    : > "$WORKDIR_DS/lists/files_$(fmt_id $i).txt"
  done
  while IFS= read -r f; do
    echo "$f" >> "$WORKDIR_DS/lists/files_$(fmt_id $j).txt"
    (( j++ )); (( j >= K )) && j=0
  done < "$alllist"

  # ---- Create per-job directories and symlink input files
  for (( i=0; i<K; i++ )); do
    local kid="$(fmt_id $i)"
    local lst="$WORKDIR_DS/lists/files_${kid}.txt"
    [[ -s "$lst" ]] || continue
    local JOBDIR INDIR OUTDIR LOG
    read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$WORKDIR_DS" "$kid")"
    mkdir -p "$INDIR" "$OUTDIR"
    rm -f "$INDIR"/*.hipo 2>/dev/null
    while IFS= read -r ff; do ln -sf "$ff" "$INDIR/"; done < "$lst"
  done

  # ---- Build list of non-empty job IDs
  local -a JOBS_DS
  JOBS_DS=()
  for (( i=0; i<K; i++ )); do
    local kid="$(fmt_id $i)"
    [[ -s "$WORKDIR_DS/lists/files_${kid}.txt" ]] && JOBS_DS+=("$kid")
  done
  echo "$P Distributing across ${#JOBS_DS[@]} job(s)."

  # ---- Two-pass: inject --pass1 so per-job step skips fid/QADB/corrections
  if (( TWO_PASS_MODE )); then
    EXTRA_ARGS_DS+=(--pass1)
    echo "$P TWO_PASS_MODE: injecting --pass1 (fid/QADB/corrections deferred to Pass 2)"
  fi

  # ---- Launch all K jobs (sliding-window batching)
  launch_jobs   # uses JOBS_DS, WORKDIR_DS, CFG_DS, EXTRA_ARGS_DS from this scope

  # ==========================================================================
  # Summarise success/failure and build merge lists
  # ==========================================================================
  local -a OKJ BADJ
  OKJ=(); BADJ=()
  typeset -A MERGE_LISTS HAVE_COUNTS MISS_COUNTS
  local _f
  for _f in "${OUTPUT_FILES[@]}"; do
    MERGE_LISTS[$_f]=""; HAVE_COUNTS[$_f]=0; MISS_COUNTS[$_f]=0
  done

  local -a MERGE_OUTPUT_FILES
  if (( TWO_PASS_MODE )); then
    MERGE_OUTPUT_FILES=(dfSelected.root)
  else
    MERGE_OUTPUT_FILES=("${OUTPUT_FILES[@]}")
  fi

  local kid JOBDIR INDIR OUTDIR LOG
  for kid in "${JOBS_DS[@]}"; do
    if job_ok "$WORKDIR_DS" "$kid"; then
      OKJ+=("$kid")
      read -r JOBDIR INDIR OUTDIR LOG <<<"$(job_paths "$WORKDIR_DS" "$kid")"
      for _f in "${MERGE_OUTPUT_FILES[@]}"; do
        if [[ -s "$OUTDIR/$_f" ]]; then
          MERGE_LISTS[$_f]="${MERGE_LISTS[$_f]} $OUTDIR/$_f"
          HAVE_COUNTS[$_f]=$(( HAVE_COUNTS[$_f] + 1 ))
        else
          MISS_COUNTS[$_f]=$(( MISS_COUNTS[$_f] + 1 ))
        fi
      done
    else
      BADJ+=("$kid")
    fi
  done

  echo "$P Success: ${#OKJ[@]} job(s) | Failed: ${#BADJ[@]} job(s)"
  (( ${#BADJ[@]} > 0 )) && echo "$P   Failed IDs: ${BADJ[*]}"

  if (( ${#BADJ[@]} > 0 && MERGE_PARTIAL == 0 )); then
    echo "$P [ERROR] MERGE_PARTIAL=0 and some jobs failed — aborting merge." >&2
    return 1
  fi

  # ==========================================================================
  # Merge each output file type
  # ==========================================================================
  local merged_any=0
  local -a _inputs
  for _f in "${MERGE_OUTPUT_FILES[@]}"; do
    _inputs=(${=MERGE_LISTS[$_f]})
    if (( ${#_inputs[@]} == 0 )); then
      echo "$P [WARN] No inputs found for $_f — skipping merge."
      continue
    fi
    local _outmerged="$OUTPUT_ABS/$_f"
    echo "$P Merging ${#_inputs[@]} file(s) → $_outmerged"
    echo "$P   (present=${HAVE_COUNTS[$_f]} missing=${MISS_COUNTS[$_f]})"
    root_merge "$_outmerged" "${_inputs[@]}" || {
      echo "$P [ERROR] root_merge failed for $_f" >&2
      return 1
    }
    merged_any=1
  done

  if (( merged_any == 0 )); then
    echo "$P [ERROR] No output files found to merge for any type." >&2
    return 1
  fi

  # ==========================================================================
  # Two-pass: Pass 2 — single reproc job on the merged dfSelected.root
  # ==========================================================================
  if (( TWO_PASS_MODE && merged_any )); then
    local pass2_input="$OUTPUT_ABS/dfSelected.root"
    if [[ ! -s "$pass2_input" ]]; then
      echo "$P [ERROR] Pass 2: dfSelected.root not found at $pass2_input" >&2
      return 1
    fi
    local pass2_log="$OUTPUT_ABS/pass2_reproc.log"
    echo "$P TWO_PASS_MODE: starting Pass 2 reproc..."
    echo "$P   Input : $pass2_input"
    echo "$P   Log   : $pass2_log"

    # Strip --pass1 from args — Pass 2 runs the full cut chain
    local -a EXTRA_ARGS_PASS2
    EXTRA_ARGS_PASS2=()
    local _arg
    for _arg in "${EXTRA_ARGS_DS[@]}"; do
      [[ "$_arg" == "--pass1" ]] && continue
      EXTRA_ARGS_PASS2+=("$_arg")
    done

    "$EXE" \
      -i "$OUTPUT_ABS" \
      -n "$NEVT" \
      -t "$PASS2_THREADS" \
      -o "$OUTPUT_ABS" \
      -c "$CFG_DS" \
      --reproc \
      --reproc-file "dfSelected.root" \
      --reproc-tree "$PASS2_REPROC_TREE" \
      "${EXTRA_ARGS_PASS2[@]}" \
      > "$pass2_log" 2>&1 || {
      echo "$P [ERROR] Pass 2 failed. See log: $pass2_log" >&2
      return 1
    }
    echo "$P Pass 2 complete."
  fi

  # ==========================================================================
  # Cleanup: remove per-job ROOT files once merges have succeeded
  # ==========================================================================
  if (( CLEANUP_JOB_OUTPUTS && merged_any )); then
    echo "$P Cleanup: removing per-job ROOT files..."
    local _deleted=0 _freed_kb=0 _jdir _rf
    for _jdir in "$WORKDIR_DS"/job_[0-9][0-9][0-9] "$WORKDIR_DS"/job_[0-9][0-9][0-9][0-9]; do
      [[ -d "$_jdir/out" ]] || continue
      for _rf in "$_jdir/out"/*.root; do
        [[ -f "$_rf" ]] || continue
        _freed_kb=$(( _freed_kb + $(du -k "$_rf" 2>/dev/null | awk '{print $1}') ))
        rm -f "$_rf"
        (( _deleted++ ))
      done
      rmdir "$_jdir/out" 2>/dev/null
    done
    echo "$P Deleted $_deleted per-job ROOT file(s), freed ~$(( _freed_kb / 1024 )) MB."
    echo "$P Job logs preserved under: $WORKDIR_DS/job_*/"
  fi

  echo "$P ============================================================"
  echo "$P Dataset $label DONE. Outputs: $OUTPUT_ABS/"
  echo "$P ============================================================"
}

# =============================================================================
# MAIN — print summary, create log dir, launch 5 dataset sub-jobs in parallel
# =============================================================================

echo "============================================================"
echo " DISANA ALL-DATASETS ORCHESTRATOR"
echo "============================================================"
echo " Mode        : $MODE  ($WAGON wagon)"
echo " Output root : $OUTPUT_ROOT/$MODE/"
echo " K jobs/ds   : $K  |  Batch size: $JOBS_PER_BATCH"
echo " TWO_PASS    : $TWO_PASS_MODE"
echo " Cleanup     : $CLEANUP_JOB_OUTPUTS"
echo "============================================================"
echo " Datasets:"
integer _ds
for (( _ds=1; _ds<=${#DS_LABELS[@]}; _ds++ )); do
  printf "   [%d] %-18s → %s/%s/%s\n" \
    "$_ds" "${DS_LABELS[$_ds]}" "$OUTPUT_ROOT" "$MODE" "${DS_OUTDIRS[$_ds]}"
done
echo "============================================================"
echo ""

# Ensure top-level log directory exists
mkdir -p "$OUTPUT_ROOT/$MODE"

# ---- Launch each dataset as a background sub-job, redirect its output to a log
typeset -a SUBJOB_PIDS SUBJOB_LABELS SUBJOB_LOGS

for (( _ds=1; _ds<=${#DS_LABELS[@]}; _ds++ )); do
  local_label="${DS_LABELS[$_ds]}"
  local_input="${DS_INPUTS[$_ds]}"
  local_outdir="${DS_OUTDIRS[$_ds]}"
  local_cfg="${DS_CFGS[$_ds]}"
  local_inb="${DS_INBENDING[$_ds]}"
  local_output="$OUTPUT_ROOT/$MODE/$local_outdir"
  local_log="$OUTPUT_ROOT/$MODE/${local_label}.log"

  (
    run_dataset "$local_label" "$local_input" "$local_output" "$local_cfg" "$local_inb"
  ) > "$local_log" 2>&1 &

  _pid=$!
  SUBJOB_PIDS+=($_pid)
  SUBJOB_LABELS+=("$local_label")
  SUBJOB_LOGS+=("$local_log")

  echo "[MAIN] Sub-job launched: $local_label"
  echo "[MAIN]   PID    : $_pid"
  echo "[MAIN]   Output : $local_output"
  echo "[MAIN]   Log    : $local_log"
  echo ""
done

echo "[MAIN] All 5 sub-jobs running. Waiting for completion..."
echo ""

# ---- Wait for all sub-jobs and collect exit statuses
all_ok=1
for (( s=1; s<=${#SUBJOB_PIDS[@]}; s++ )); do
  _pid="${SUBJOB_PIDS[$s]}"
  _lbl="${SUBJOB_LABELS[$s]}"
  _log="${SUBJOB_LOGS[$s]}"
  wait "$_pid"
  _st=$?
  if (( _st == 0 )); then
    echo "[MAIN] ✓  $_lbl — completed successfully."
  else
    echo "[MAIN] ✗  $_lbl — FAILED (exit $_st).  See log: $_log"
    all_ok=0
  fi
done

echo ""
echo "============================================================"
if (( all_ok )); then
  echo " All 5 datasets completed successfully."
  echo " Merged outputs are under: $OUTPUT_ROOT/$MODE/"
else
  echo " One or more datasets FAILED — review the logs above."
fi
echo "============================================================"

(( all_ok )) || exit 1
exit 0
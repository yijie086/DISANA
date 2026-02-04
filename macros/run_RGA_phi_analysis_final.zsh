#!/usr/bin/env zsh
# ------------------------------------------------------------
# Batch runner for AnalysisPhi / RunPhiAnalysis
# Improved version with custom reproc file/tree options
# ------------------------------------------------------------

module purge

# --- make sure 'module' command is available ---
if ! command -v module &>/dev/null; then
  if [[ -r /etc/profile.d/modules.sh ]]; then
    source /etc/profile.d/modules.sh
  fi
fi

# --- load CLAS12 environment ---
module load clas12/5.3
module unload qadb/3.2.0  || true
module load qadb/3.4

# ---- Defaults (edit here) ----
EXECUTABLE="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/build/AnalysisPhi"
NFILE=-1        # -1 => all files
THREADS=1       # single-threaded
IS_MC=0
REPROC=0
MINIMAL=0
MISSING_KM=0    # 0 = exclusive DVKpKmP, 1 = kaon-missing nSidis
INBENDING=0     # explicit inbending flag

# NEW: Custom reproc options
REPROC_FILE=""      # Custom ROOT filename
REPROC_TREE=""      # Custom tree name
CUSTOM_INPUT=""     # Custom input directory (overrides default paths)
CUSTOM_OUTPUT=""    # Custom output directory (overrides default paths)

# ---- Known datasets ----
CONFIGS=(rgasp18_inb rgasp18_outb rgafall18_inb rgafall18_outb rgasp19_inb)

# ---- helpers ----
me="$(basename "$0")"

die() {
  echo "$@" >&2
  exit 1
}

# ---- usage ----
show_usage() {
  cat <<EOF
Usage:
  $me -A | -c <dataconfig>  [-n <nfile>] [-t <threads>] [-x <executable>]
         [--mc] [--reproc] [--inbending] [--minimal] [--missingKm]
         [--reproc-file <filename>] [--reproc-tree <treename>]
         [--custom-input <path>] [--custom-output <path>]
         [--dry-run] [-h]

Options:
  -A                         Run all configurations
  -c <dataconfig>            Run specific configuration (e.g., rgasp18_inb)
  -n <nfile>                 Number of files to process (-1 for all)
  -t <threads>               Number of threads (default: 1)
  -x <executable>            Path to executable
  
  --mc                       Monte Carlo data flag
  --reproc                   Reprocessed ROOT file flag
  --inbending                Force inbending flag (auto-detected from config name)
  --minimal                  Minimal booking flag
  --missingKm                Kaon-missing nSidis mode
  
  --reproc-file <filename>   Custom ROOT filename (e.g., "MyData.root")
  --reproc-tree <treename>   Custom tree name (e.g., "myTree", "dfSelected", "hipo")
  --custom-input <path>      Custom input directory (overrides defaults)
  --custom-output <path>     Custom output directory (overrides defaults)
  
  --dry-run                  Show commands without executing
  -h                         Show this help

Notes:
  - Input/output paths are fixed by default based on dataset
  - Use --custom-input/--custom-output to override defaults
  - --missingKm switches input dirs to nSidis/ and passes --missingKm to AnalysisPhi
  - --inbending is auto-detected from config name (e.g., *_inb) but can be forced

Examples:
  # Run all configs with 200 files
  $me -A -n 200
  
  # Run specific config with minimal booking
  $me -c rgasp18_inb --minimal
  
  # Run with kaon-missing mode
  $me -c rgasp18_inb --missingKm
  
  # Use custom reprocessed file with custom tree
  $me -c rgasp18_outb --reproc --reproc-file MyAnalysis.root --reproc-tree dfSelected
  
  # Use completely custom input/output paths
  $me -c rgasp18_inb --custom-input /my/data/path --custom-output /my/results
EOF
}

if (( $# == 0 )); then
  show_usage
  exit 0
fi

# ---- parse args ----
RUN_ALL=0
ONE_CONFIG=""
DRY=0

while (( $# > 0 )); do
  a="$1"
  case "$a" in
    -h|--help)
      show_usage
      exit 0
      ;;
    -A)
      RUN_ALL=1
      shift
      ;;
    -c)
      (( $# < 2 )) && die "ERROR: -c requires a dataconfig"
      ONE_CONFIG="$2"
      shift 2
      ;;
    -n|--nfile)
      (( $# < 2 )) && die "ERROR: -n requires an integer"
      NFILE="$2"
      shift 2
      ;;
    -t|--threads)
      (( $# < 2 )) && die "ERROR: -t/--threads requires an integer"
      THREADS="$2"
      [[ "$THREADS" == <-> ]] || die "ERROR: threads must be an integer (got: $THREADS)"
      (( THREADS >= 1 )) || die "ERROR: threads must be >= 1 (got: $THREADS)"
      shift 2
      ;;
    -x|--executable)
      (( $# < 2 )) && die "ERROR: -x requires a path"
      EXECUTABLE="$2"
      shift 2
      ;;
    --mc)
      IS_MC=1
      shift
      ;;
    --reproc)
      REPROC=1
      shift
      ;;
    --inbending)
      INBENDING=1
      shift
      ;;
    --minimal)
      MINIMAL=1
      shift
      ;;
    --missingKm)
      MISSING_KM=1
      shift
      ;;
    --reproc-file)
      (( $# < 2 )) && die "ERROR: --reproc-file requires a filename"
      REPROC_FILE="$2"
      shift 2
      ;;
    --reproc-tree)
      (( $# < 2 )) && die "ERROR: --reproc-tree requires a tree name"
      REPROC_TREE="$2"
      shift 2
      ;;
    --custom-input)
      (( $# < 2 )) && die "ERROR: --custom-input requires a path"
      CUSTOM_INPUT="$2"
      shift 2
      ;;
    --custom-output)
      (( $# < 2 )) && die "ERROR: --custom-output requires a path"
      CUSTOM_OUTPUT="$2"
      shift 2
      ;;
    --dry-run)
      DRY=1
      shift
      ;;
    *)
      die "ERROR: unknown option: $a"
      ;;
  esac
done

echo "[$me] Start: $(date)"

# ---- decide OUT_BASE now that MISSING_KM is known ----
if [[ -n "$CUSTOM_OUTPUT" ]]; then
  OUT_BASE="$CUSTOM_OUTPUT"
  echo "[$me] Using custom output base: $OUT_BASE"
elif (( MISSING_KM )); then
  OUT_BASE="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/nSIDIS/"
  echo "[$me] Running for Kaon Missing case (nSidis inputs)"
else
  OUT_BASE="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DVKpKm/"
  echo "[$me] Running for Exclusive Kaon case (DVKpKmP inputs)"
fi
echo "[$me] OUT_BASE = $OUT_BASE"

# ---- decide run list ----
RUN_LIST=()
if (( RUN_ALL )); then
  RUN_LIST=("${CONFIGS[@]}")
elif [[ -n "$ONE_CONFIG" ]]; then
  RUN_LIST=("$ONE_CONFIG")
else
  die "ERROR: choose -A or -c <dataconfig>"
fi

# ---- Per-dataset input directories, depending on MISSING_KM ----
# Only set these if CUSTOM_INPUT is not provided
if [[ -z "$CUSTOM_INPUT" ]]; then
  if (( MISSING_KM )); then
    IN_rgasp18_inb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/nSidis/"
    IN_rgasp18_outb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/nSidis/"
    IN_rgafall18_inb="/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/"
    IN_rgafall18_outb="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/"
    IN_rgasp19_inb="/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/nSidis/"
  else
    IN_rgasp18_inb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/DVKpKmP/"
    IN_rgasp18_outb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/DVKpKmP/"
    IN_rgafall18_inb="/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/DVKpKmP/"
    IN_rgafall18_outb="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/DVKpKmP/"
    IN_rgasp19_inb="/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/"
  fi
fi

# ---- Preflight checks ----
echo "[$me] Preflight: checking EXECUTABLE"
echo "[$me] EXECUTABLE = $EXECUTABLE"

if [[ ! -f "$EXECUTABLE" ]]; then
  die "ERROR: executable not found at: $EXECUTABLE"
fi

if [[ ! -x "$EXECUTABLE" ]]; then
  die "ERROR: file exists but is not executable: $EXECUTABLE
HINT: chmod +x $EXECUTABLE"
fi

echo "[$me] ✓ Executable found and is executable"

echo "[$me] Preflight: creating OUT_BASE = $OUT_BASE"
mkdir -p "$OUT_BASE" || die "ERROR: cannot create OUT_BASE: $OUT_BASE"
echo "[$me] ✓ Output base directory ready"

# ---- Validate reproc options ----
if (( REPROC )); then
  echo "[$me] Reprocessed file mode enabled"
  if [[ -n "$REPROC_FILE" ]]; then
    echo "[$me]   Custom ROOT file: $REPROC_FILE"
  else
    echo "[$me]   Using default ROOT file from RunPhiAnalysis"
  fi
  if [[ -n "$REPROC_TREE" ]]; then
    echo "[$me]   Custom tree name: $REPROC_TREE"
  else
    echo "[$me]   Using default tree name from RunPhiAnalysis"
  fi
fi

# ---- run all configs (can be parallelized) ----
pids=()

for cfg in "${RUN_LIST[@]}"; do
  echo ""
  echo "=========================================="
  echo "[$me] Processing config: $cfg"
  echo "=========================================="

  # Validate config
  ok=0
  for k in "${CONFIGS[@]}"; do
    [[ "$cfg" == "$k" ]] && ok=1
  done
  (( ok )) || die "ERROR: unknown dataconfig $cfg"

  # Determine input path
  if [[ -n "$CUSTOM_INPUT" ]]; then
    inpath="$CUSTOM_INPUT"
    echo "[$me]   Using custom input: $inpath"
  else
    # indirect: IN_rgasp18_inb -> inpath
    varname="IN_${cfg}"
    inpath="${(P)varname}"
    
    [[ -z "$inpath" ]] && die "ERROR: input path for $cfg not set in script"
    echo "[$me]   Using standard input: $inpath"
  fi

  if [[ ! -d "$inpath" ]]; then
    die "ERROR: input path for $cfg not a directory: $inpath"
  fi

  # make per-config outdir
  outdir="${OUT_BASE}/${cfg}/"
  echo "[$me]   Creating output dir: $outdir"
  mkdir -p "$outdir" || die "ERROR: cannot create outdir $outdir"
  echo "[$me]   ✓ Output directory ready"

  # timestamp
  stamp="$(/bin/date +%Y%m%d_%H%M%S)"
  logfile="${outdir}/run_${cfg}_${stamp}.txt"

  # build flags array to match AnalysisPhi CLI
  flags=( -i "$inpath" -n "$NFILE" -t "$THREADS" -o "$outdir" -c "$cfg" )
  
  # Add boolean flags
  (( IS_MC ))      && flags+=( --mc )
  (( REPROC ))     && flags+=( --reproc )
  (( MINIMAL ))    && flags+=( --minimal )
  (( MISSING_KM )) && flags+=( --missingKm )
  
  # Add inbending flag: either explicit or auto-detected from config name
  if (( INBENDING )); then
    flags+=( --inbending )
    echo "[$me]   ✓ Inbending mode: FORCED by --inbending flag"
  elif [[ "$cfg" == *_inb ]]; then
    flags+=( --inbending )
    echo "[$me]   ✓ Inbending mode: AUTO-DETECTED from config name"
  else
    echo "[$me]   ✓ Outbending mode (default)"
  fi
  
  # Add reproc file/tree options if provided
  if [[ -n "$REPROC_FILE" ]]; then
    flags+=( --reproc-file "$REPROC_FILE" )
  fi
  if [[ -n "$REPROC_TREE" ]]; then
    flags+=( --reproc-tree "$REPROC_TREE" )
  fi

  echo "------------------------------------------------------------------"
  echo " Config       : $cfg"
  echo " InputDir     : $inpath"
  echo " OutDir       : $outdir"
  echo " NFile        : $NFILE"
  echo " Threads      : $THREADS"
  echo " Flags        : mc=$IS_MC reproc=$REPROC minimal=$MINIMAL missingKm=$MISSING_KM"
  if (( INBENDING )) || [[ "$cfg" == *_inb ]]; then
    echo "                inbending=1"
  fi
  if [[ -n "$REPROC_FILE" ]]; then
    echo " ReprocFile   : $REPROC_FILE"
  fi
  if [[ -n "$REPROC_TREE" ]]; then
    echo " ReprocTree   : $REPROC_TREE"
  fi
  echo " Executable   : $EXECUTABLE"
  echo " Full Command : $EXECUTABLE ${flags[*]}"
  echo " LogFile      : $logfile"
  echo "------------------------------------------------------------------"

  (( DRY )) && continue

  # Run the command
  {
    echo "[$me] Starting $cfg at $(date)"
    "$EXECUTABLE" "${flags[@]}" 2>&1 | tee "$logfile"
    rc=${pipestatus[1]}
    echo "[$me] ${cfg} finished at $(date) with exit status: $rc" | tee -a "$logfile"
  } &

  pids+=$!
done

# ---- wait for all parallel jobs to finish ----
if (( !DRY )); then
  echo ""
  echo "=========================================="
  echo "[$me] Waiting for all jobs to complete..."
  echo "=========================================="
  
  for pid in "${pids[@]}"; do
    wait "$pid"
    wait_status=$?
    echo "[$me] Job $pid completed with status: $wait_status"
  done
  
  echo ""
  echo "=========================================="
  echo "[$me] All requested runs completed at $(date)"
  echo "=========================================="
else
  echo ""
  echo "[$me] Dry run complete. No jobs were actually executed."
fi
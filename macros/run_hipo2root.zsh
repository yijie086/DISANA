#!/usr/bin/env zsh
# ------------------------------------------------------------
# Batch runner for hipo2root converter
# Converts HIPO files to ROOT format with parallel processing
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

# ---- Defaults (edit here) ----
EXECUTABLE="/w/hallb-scshelf2102/clas12/singh/Softwares/HIPO_test/newhipo2root/hipo-utils/build/hipo2root"
OUTPUT_DIR="./test_outputs/"
MERGED_FILE="DVKpKm.root"
TREE_NAME="hipo"
NTHREADS=10
NEVT_INSPECT=1000
INPUT_PATTERN=""
INCLUDE_BANKS=""
EXCLUDE_BANKS=""

# ---- Known datasets ----
DATASETS=(
  "rgasp18_inb:/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/DVKpKmP/"
  "rgasp18_outb:/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/DVKpKmP/"
  "rgafall18_inb:/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/DVKpKmP/"
  "rgafall18_outb:/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/DVKpKmP/"
  "rgasp19_inb:/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/"
  "rgasp19_inb_nsidis:/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/nSidis/"
)

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
  $me -i <input_pattern> [-o <output_dir>] [-m <merged_file>] [-t <tree_name>]
         [-j <threads>] [-n <nevt_inspect>] [-I <banks>] [-E <banks>]
         [-d <dataset>] [--dry-run] [-h]

Required:
  -i <input_pattern>    Input HIPO files (glob pattern, e.g., "/path/*.hipo")
                        OR use -d <dataset> for predefined paths

Optional:
  -o <output_dir>       Output directory (default: $OUTPUT_DIR)
  -m <merged_file>      Merged output ROOT file name (default: $MERGED_FILE)
  -t <tree_name>        TTree name in output (default: $TREE_NAME)
  -j <threads>          Number of parallel jobs (default: $NTHREADS)
  -n <nevt_inspect>     Events to inspect for schema (default: $NEVT_INSPECT)
  -I <banks>            Include only these banks (comma-separated)
  -E <banks>            Exclude these banks (comma-separated)
  -d <dataset>          Use predefined dataset path (see list below)
  --dry-run             Show commands without executing
  -h                    Show this help

Predefined Datasets:
EOF
  for ds in "${DATASETS[@]}"; do
    name="${ds%%:*}"
    path="${ds#*:}"
    printf "  %-20s %s\n" "$name" "$path"
  done
  cat <<EOF

Examples:
  # Convert specific files
  $me -i "/cache/clas12/rg-a/.../DVKpKmP/*.hipo" -o output/ -j 8

  # Use predefined dataset
  $me -d rgasp19_inb -o output/sp19/ -j 16

  # Convert with bank filtering
  $me -d rgasp18_inb -I "REC::Particle,REC::Calorimeter" -o output/

  # Exclude certain banks
  $me -d rgafall18_inb -E "MC::Particle,MC::Event" -j 12

  # Custom tree name and merged file
  $me -i "/path/*.hipo" -o out/ -m combined.root -t mydata
EOF
}

if (( $# == 0 )); then
  show_usage
  exit 0
fi

# ---- parse args ----
DRY_RUN=0
USE_DATASET=""

while (( $# > 0 )); do
  a="$1"
  case "$a" in
    -h|--help)
      show_usage
      exit 0
      ;;
    -i)
      (( $# < 2 )) && die "ERROR: -i requires an input pattern"
      INPUT_PATTERN="$2"
      shift 2
      ;;
    -o)
      (( $# < 2 )) && die "ERROR: -o requires an output directory"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -m)
      (( $# < 2 )) && die "ERROR: -m requires a merged filename"
      MERGED_FILE="$2"
      shift 2
      ;;
    -t)
      (( $# < 2 )) && die "ERROR: -t requires a tree name"
      TREE_NAME="$2"
      shift 2
      ;;
    -j)
      (( $# < 2 )) && die "ERROR: -j requires number of threads"
      NTHREADS="$2"
      [[ "$NTHREADS" == <-> ]] || die "ERROR: threads must be an integer (got: $NTHREADS)"
      (( NTHREADS >= 1 )) || die "ERROR: threads must be >= 1 (got: $NTHREADS)"
      shift 2
      ;;
    -n)
      (( $# < 2 )) && die "ERROR: -n requires number of events"
      NEVT_INSPECT="$2"
      [[ "$NEVT_INSPECT" == <-> ]] || die "ERROR: nevt must be an integer (got: $NEVT_INSPECT)"
      shift 2
      ;;
    -I)
      (( $# < 2 )) && die "ERROR: -I requires bank list"
      INCLUDE_BANKS="$2"
      shift 2
      ;;
    -E)
      (( $# < 2 )) && die "ERROR: -E requires bank list"
      EXCLUDE_BANKS="$2"
      shift 2
      ;;
    -d)
      (( $# < 2 )) && die "ERROR: -d requires a dataset name"
      USE_DATASET="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    *)
      die "ERROR: unknown option: $a"
      ;;
  esac
done

echo "[$me] Start: $(date)"

# ---- resolve input pattern from dataset if specified ----
if [[ -n "$USE_DATASET" ]]; then
  found=0
  for ds in "${DATASETS[@]}"; do
    name="${ds%%:*}"
    path="${ds#*:}"
    if [[ "$name" == "$USE_DATASET" ]]; then
      INPUT_PATTERN="${path}*"
      echo "[$me] Using dataset '$USE_DATASET': $path"
      found=1
      break
    fi
  done
  (( found )) || die "ERROR: unknown dataset '$USE_DATASET'"
fi

# ---- validate required arguments ----
[[ -z "$INPUT_PATTERN" ]] && die "ERROR: input pattern required (-i or -d)"

# ---- preflight checks ----
echo "[$me] Preflight: checking EXECUTABLE"
echo "[$me] EXECUTABLE = $EXECUTABLE"

if [[ ! -f "$EXECUTABLE" ]]; then
  die "ERROR: executable not found at: $EXECUTABLE"
fi

if [[ ! -x "$EXECUTABLE" ]]; then
  die "ERROR: file exists but is not executable: $EXECUTABLE"
fi

echo "[$me] Preflight: creating OUTPUT_DIR = $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
if (( $? != 0 )); then
  die "ERROR: cannot create output directory: $OUTPUT_DIR"
fi

# ---- build command ----
cmd=( "$EXECUTABLE" )
cmd+=( -o "$OUTPUT_DIR" )
cmd+=( -m "$MERGED_FILE" )
cmd+=( -t "$TREE_NAME" )
cmd+=( -j "$NTHREADS" )
cmd+=( -n "$NEVT_INSPECT" )

if [[ -n "$INCLUDE_BANKS" ]]; then
  cmd+=( -I "$INCLUDE_BANKS" )
fi

if [[ -n "$EXCLUDE_BANKS" ]]; then
  cmd+=( -E "$EXCLUDE_BANKS" )
fi

cmd+=( $INPUT_PATTERN )

# ---- display configuration ----
stamp="$(/bin/date +%Y%m%d_%H%M%S)"
logfile="${OUTPUT_DIR}/hipo2root_${stamp}.log"

echo "=================================================================="
echo " hipo2root Conversion"
echo "=================================================================="
echo " Input Pattern  : $INPUT_PATTERN"
echo " Output Dir     : $OUTPUT_DIR"
echo " Merged File    : $MERGED_FILE"
echo " Tree Name      : $TREE_NAME"
echo " Threads        : $NTHREADS"
echo " Events Inspect : $NEVT_INSPECT"
[[ -n "$INCLUDE_BANKS" ]] && echo " Include Banks  : $INCLUDE_BANKS"
[[ -n "$EXCLUDE_BANKS" ]] && echo " Exclude Banks  : $EXCLUDE_BANKS"
echo " Executable     : $EXECUTABLE"
echo " Command        : ${cmd[*]}"
echo " Log File       : $logfile"
echo "=================================================================="

if (( DRY_RUN )); then
  echo "[$me] DRY RUN - command would be executed as shown above"
  exit 0
fi

# ---- execute conversion ----
echo "[$me] Starting conversion..."
"${cmd[@]}" 2>&1 | tee "$logfile"
rc=${pipestatus[1]}

echo ""
echo "[$me] Conversion finished with exit status: $rc" | tee -a "$logfile"
echo "[$me] End: $(date)" | tee -a "$logfile"

if (( rc == 0 )); then
  echo "[$me] SUCCESS: Output file should be at: ${OUTPUT_DIR}/${MERGED_FILE}"
  ls -lh "${OUTPUT_DIR}/${MERGED_FILE}" 2>/dev/null || true
else
  echo "[$me] ERROR: Conversion failed with status $rc" >&2
  echo "[$me] Check log file: $logfile" >&2
  exit $rc
fi

exit 0
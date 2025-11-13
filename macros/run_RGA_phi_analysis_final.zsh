#!/usr/bin/env zsh
# ------------------------------------------------------------
# Batch runner for AnalysisPhi / RunPhiAnalysis (single-threaded)
# zsh version converted from run_RGA_phi_analysis.csh
# ------------------------------------------------------------
module purge
# --- make sure 'module' command is available ---
if ! command -v module &>/dev/null; then
  # Lmod / environment-modules are usually initialized here on many clusters
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
OUT_BASE="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/source/macros/results"   # outputs go to ${OUT_BASE}/${config}
NFILE=-1        # -1 => all files
THREADS=1       # single-threaded
IS_MC=0
REPROC=0
MINIMAL=0

# ---- Known datasets ----
CONFIGS=(rgasp18_inb rgasp18_outb rgafall18_inb rgafall18_outb rgasp19_inb)

# ---- Per-dataset input directories ----
IN_rgasp18_inb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/DVKpKmP/"
IN_rgasp18_outb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/DVKpKmP/"
IN_rgafall18_inb="/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/DVKpKmP/"
IN_rgafall18_outb="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/DVKpKmP/"
IN_rgasp19_inb="/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/"

# ---- helpers ----
me="$(basename "$0")"

die() {
  echo "$@" >&2
  exit 1
}

# ---- usage ----
if (( $# == 0 )); then
  cat <<EOF
Usage:
  $me -A | -c <dataconfig>  [-n <nfile>] [-x <executable>] [--mc] [--reproc] [--minimal] [--dry-run] [-h]
Notes: input/output paths are fixed in the script; outputs go to \${OUT_BASE}/<dataconfig>/
Examples:
  $me -A -n 200
  $me -c rgasp18_inb --minimal
EOF
  exit 0
fi

# ---- parse args (no --in-* or -o supported) ----
RUN_ALL=0
ONE_CONFIG=""
DRY=0

while (( $# > 0 )); do
  a="$1"
  case "$a" in
    -h)
      cat <<EOF
Usage:
  $me -A | -c <dataconfig>  [-n <nfile>] [-x <executable>] [--mc] [--reproc] [--minimal] [--dry-run] [-h]
EOF
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
    -n)
      (( $# < 2 )) && die "ERROR: -n requires an integer"
      NFILE="$2"
      shift 2
      ;;
    -x)
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
    --minimal)
      MINIMAL=1
      shift
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

# ---- decide run list ----
RUN_LIST=()
if (( RUN_ALL )); then
  RUN_LIST=("${CONFIGS[@]}")
elif [[ -n "$ONE_CONFIG" ]]; then
  RUN_LIST=("$ONE_CONFIG")
else
  die "ERROR: choose -A or -c <dataconfig>"
fi

echo "[$me] Test Debug = 1"

# ---- preflight (no 'die', print everything) ----
echo "[$me] Preflight: checking EXECUTABLE"
echo "[$me] EXECUTABLE = $EXECUTABLE"

ls -l "$EXECUTABLE"
ls_status=$?
echo "[$me] ls status = $ls_status"
if (( ls_status != 0 )); then
  echo "ERROR: executable not found at: $EXECUTABLE" >&2
  exit 2
fi

if [[ ! -x "$EXECUTABLE" ]]; then
  echo "ERROR: file exists but is not executable: $EXECUTABLE" >&2
  echo "HINT: chmod +x $EXECUTABLE" >&2
  exit 2
fi

echo "[$me] Preflight: creating OUT_BASE = $OUT_BASE"
mkdir -p "$OUT_BASE"
mk_status=$?
echo "[$me] mkdir status = $mk_status"
# Uncomment if you want hard failure here:
# if (( mk_status != 0 )); then
#   echo "ERROR: cannot create OUT_BASE $OUT_BASE (status $mk_status)" >&2
#   exit 2
# fi

echo "[$me] Test Debug = 2"

# ---- run ----
for cfg in "${RUN_LIST[@]}"; do
  echo "[$me] running config: $cfg"

  ok=0
  for k in "${CONFIGS[@]}"; do
    [[ "$cfg" == "$k" ]] && ok=1
  done
  (( ok )) || die "ERROR: unknown dataconfig $cfg"

  echo "[$me] Test Debug = 3"

  # indirect: IN_rgasp18_inb -> inpath
  varname="IN_${cfg}"
  # ${(P)varname} is zsh: value of variable whose name is in $varname
  inpath="${(P)varname}"

  [[ -z "$inpath" ]] && die "ERROR: input path for $cfg not set in script"
  [[ -d "$inpath" ]] || die "ERROR: input path for $cfg not a directory: $inpath"

  echo "[$me] Test Debug = 4"

  # make per-config outdir
  outdir="${OUT_BASE}/${cfg}/"
  echo "[$me] creating outdir: $outdir"
  /bin/mkdir -p "$outdir"
  mkcfg_status=$?
  echo "[$me] mkdir(outdir) status = $mkcfg_status"
  if (( mkcfg_status != 0 )); then
    echo "ERROR: cannot create outdir $outdir (status $mkcfg_status)" >&2
    exit 2
  fi

  # timestamp
  stamp="$(/bin/date +%Y%m%d_%H%M%S)"
  echo "[$me] stamp = $stamp"

  # logfile path
  logfile="${outdir}/run_${cfg}_${stamp}.txt"
  echo "[$me] logfile = $logfile"

  echo "[$me] Test Debug = 5"

  # build flags array
  flags=( -i "$inpath" -n "$NFILE" -t "$THREADS" -o "$outdir" -c "$cfg" )
  (( IS_MC ))   && flags+=( --mc )
  (( REPROC ))  && flags+=( --reproc )
  (( MINIMAL )) && flags+=( --minimal )

  echo "[$me] Test Debug = 6"
  echo "------------------------------------------------------------------"
  echo " Config   : $cfg"
  echo " InputDir : $inpath"
  echo " OutDir   : $outdir"
  echo " NFile    : $NFILE"
  echo " Threads  : $THREADS"
  echo " Exec     : $EXECUTABLE ${flags[*]}"
  echo " LogFile  : $logfile"
  echo "------------------------------------------------------------------"

  (( DRY )) && continue

  # run executable; keep its exit status (pipestatus[1] = first cmd in pipeline)
  "$EXECUTABLE" "${flags[@]}" 2>&1 | tee "$logfile"
  rc=${pipestatus[1]}
  echo "[$me] ${cfg} exit status: $rc" | tee -a "$logfile"
done

echo "All requested runs completed."

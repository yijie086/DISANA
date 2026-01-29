#!/usr/bin/env zsh
# ------------------------------------------------------------
# Batch runner for AnalysisPhi / RunPhiAnalysis (single-threaded)
# zsh version converted from run_RGA_phi_analysis.csh
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

# ---- Known datasets ----
CONFIGS=(rgasp18_inb rgasp18_outb rgafall18_inb rgafall18_outb rgasp19_inb)

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
  $me -A | -c <dataconfig>  [-n <nfile>] [-x <executable>]
         [--mc] [--reproc] [--minimal] [--missingKm] [--dry-run] [-h]

Notes:
  - input/output paths are fixed in the script
  - outputs go to \${OUT_BASE}/<dataconfig>/
  - --missingKm switches input dirs to nSidis/ and passes --missingKm to AnalysisPhi

Examples:
  $me -A -n 200
  $me -c rgasp18_inb --minimal
  $me -c rgasp18_inb --missingKm
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
  $me -A | -c <dataconfig>  [-n <nfile>] [-x <executable>]
         [--mc] [--reproc] [--minimal] [--missingKm] [--dry-run] [-h]
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
    -t)
      (( $# < 2 )) && die "ERROR: -t/--threads requires an integer"
      THREADS="$2"
      [[ "$THREADS" == <-> ]] || die "ERROR: threads must be an integer (got: $THREADS)"
      (( THREADS >= 1 )) || die "ERROR: threads must be >= 1 (got: $THREADS)"
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
    --missingKm)
      MISSING_KM=1
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

# ---- decide OUT_BASE now that MISSING_KM is known ----
if (( MISSING_KM )); then
  OUT_BASE="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/source/macros/n_sidis_wagon_missing/"
  echo "[$me] Running for Kaon Missing case (nSidis inputs)"
else
  OUT_BASE="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/source/macros/NewCached_KpKm/"
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

# ---- run all configs in parallel ----
pids=()

for cfg in "${RUN_LIST[@]}"; do
  echo "[$me] running config: $cfg"

  ok=0
  for k in "${CONFIGS[@]}"; do
    [[ "$cfg" == "$k" ]] && ok=1
  done
  (( ok )) || die "ERROR: unknown dataconfig $cfg"

  # indirect: IN_rgasp18_inb -> inpath
  varname="IN_${cfg}"
  inpath="${(P)varname}"

  [[ -z "$inpath" ]] && die "ERROR: input path for $cfg not set in script"
  [[ -d "$inpath" ]] || die "ERROR: input path for $cfg not a directory: $inpath"

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
  logfile="${outdir}/run_${cfg}_${stamp}.txt"

  # build flags array to match AnalysisPhi CLI
  flags=( -i "$inpath" -n "$NFILE" -t "$THREADS" -o "$outdir" -c "$cfg" )
  (( IS_MC ))      && flags+=( --mc )
  (( REPROC ))     && flags+=( --reproc )
  (( MINIMAL ))    && flags+=( --minimal )
  (( MISSING_KM )) && flags+=( --missingKm )
  # Add inbending flag automatically based on config name
  if (("$cfg" == *_inb )); then
   flags+=( --inbending )
  fi

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

  {
    "$EXECUTABLE" "${flags[@]}" 2>&1 | tee "$logfile"
    rc=${pipestatus[1]}
    echo "[$me] ${cfg} exit status: $rc" | tee -a "$logfile"
  } &

  pids+=$!
done

# ---- wait for all parallel jobs to finish ----
for pid in "${pids[@]}"; do
  wait "$pid"
done

echo "All requested runs completed."

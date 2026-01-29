#!/bin/tcsh -f
# ------------------------------------------------------------
# Batch runner for AnalysisPhi / RunPhiAnalysis (single-threaded)
# No CLI args for input or output paths; they are defined below.
# ------------------------------------------------------------

# ---- Defaults (edit here) ----
set EXECUTABLE = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/build/AnalysisPhi"
set OUT_BASE   = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/source/macros/results"   # outputs go to ${OUT_BASE}/${config}
set NFILE      = -1            # -1 => all files
set THREADS    = 1             # single-threaded
set IS_MC      = 0
set REPROC     = 0
set MINIMAL    = 0

# ---- Known datasets ----
set CONFIGS = ( rgasp18_inb rgasp18_outb rgafall18_inb rgafall18_outb rgasp19_inb )

# ---- Per-dataset input directories ----
set IN_rgasp18_inb    = "/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/DVKpKmP/"
set IN_rgasp18_outb   = "/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/DVKpKmP/"
set IN_rgafall18_inb  = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/DVKpKmP/"
set IN_rgafall18_outb = "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/DVKpKmP/"
set IN_rgasp19_inb    = "/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/"

# ---- helpers ----
set me = "$0:t"
alias die 'echo "\!:*" >&2; exit 1'
alias ts  'date "+%Y%m%d_%H%M%S"'

# ---- usage ----
if ( $#argv == 0 ) then
  echo "Usage:"
  echo "  $me -A | -c <dataconfig>  [-n <nfile>] [-x <executable>] [--mc] [--reproc] [--minimal] [--dry-run] [-h]"
  echo "Notes: input/output paths are fixed in the script; outputs go to ${OUT_BASE}/<dataconfig>/"
  echo "Examples:"
  echo "  $me -A -n 200"
  echo "  $me -c rgasp18_inb --minimal"
  exit 0
endif

# ---- parse args (no --in-* or -o supported) ----
set RUN_ALL    = 0
set ONE_CONFIG = ""
set DRY        = 0

while ( $#argv > 0 )
  set a = "$1"
  switch ( "$a" )
    case "-h":
      echo "Usage:"
      echo "  $me -A | -c <dataconfig>  [-n <nfile>] [-x <executable>] [--mc] [--reproc] [--minimal] [--dry-run] [-h]"
      exit 0
      breaksw
    case "-A":
      set RUN_ALL = 1; shift; breaksw
    case "-c":
      if ( $#argv < 2 ) die "ERROR: -c requires a dataconfig"
      set ONE_CONFIG = "$2"; shift; shift; breaksw
    case "-n":
      if ( $#argv < 2 ) die "ERROR: -n requires an integer"
      set NFILE = "$2"; shift; shift; breaksw
    case "-x":
      if ( $#argv < 2 ) die "ERROR: -x requires a path"
      set EXECUTABLE = "$2"; shift; shift; breaksw
    case "--mc":
      set IS_MC = 1; shift; breaksw
    case "--reproc":
      set REPROC = 1; shift; breaksw
    case "--minimal":
      set MINIMAL = 1; shift; breaksw
    case "--dry-run":
      set DRY = 1; shift; breaksw
    default:
      die "ERROR: unknown option: $a"
  endsw
end

echo "[$me] Start: `date`"

# ---- decide run list ----
set RUN_LIST = ()
if ( $RUN_ALL ) then
  set RUN_LIST = ( $CONFIGS )
else if ( "$ONE_CONFIG" != "" ) then
  set RUN_LIST = ( "$ONE_CONFIG" )
else
  die "ERROR: choose -A or -c <dataconfig>"
endif

echo "[$me] Test Debug = 1"

# ---- preflight (no 'die', print everything) ----
echo "[$me] Preflight: checking EXECUTABLE"
echo "[$me] EXECUTABLE = $EXECUTABLE"

# show existence + perms
/bin/ls -l "$EXECUTABLE"
set ls_status = $status
echo "[$me] ls status = $ls_status"
if ( $ls_status != 0 ) then
  echo "ERROR: executable not found at: $EXECUTABLE" >&2
  exit 2
endif

# explicit test for -x
if ( ! -x "$EXECUTABLE" ) then
  echo "ERROR: file exists but is not executable: $EXECUTABLE" >&2
  echo "HINT: chmod +x $EXECUTABLE" >&2
  exit 2
endif

echo "[$me] Preflight: creating OUT_BASE = $OUT_BASE"
mkdir -p "$OUT_BASE"
set mk_status = $status
echo "[$me] mkdir status = $mk_status"
#if ( $mk_status != 0 ) then
 # echo "ERROR: cannot create OUT_BASE $OUT_BASE (status $mk_status)" >&2
  #exit 2
#endif

echo "[$me] Test Debug = 2"

# ---- run ----
foreach cfg ( $RUN_LIST )
  echo "[$me] running config: $cfg"
  set ok = 0
  foreach k ( $CONFIGS )
    if ( "$cfg" == "$k" ) set ok = 1
  end
  if ( $ok == 0 ) die "ERROR: unknown dataconfig $cfg"
echo "[$me] Test Debug = 3"
  set varname = "IN_${cfg}"
  eval set inpath = \"\$$varname\"
  if ( "$inpath" == "" ) die "ERROR: input path for $cfg not set in script"
  if ( ! -d "$inpath" ) die "ERROR: input path for $cfg not a directory: $inpath"
echo "[$me] Test Debug = 4"

# make per-config outdir with explicit /bin/mkdir and status print
set outdir = "${OUT_BASE}/${cfg}/"
echo "[$me] creating outdir: $outdir"
/bin/mkdir -p "$outdir"
set mkcfg_status = $status
echo "[$me] mkdir(outdir) status = $mkcfg_status"
if ( $mkcfg_status != 0 ) then
  echo "ERROR: cannot create outdir $outdir (status $mkcfg_status)" >&2
  exit 2
endif

# compute timestamp without using alias expansion
set stamp = "`/bin/date +%Y%m%d_%H%M%S`"
echo "[$me] stamp = $stamp"

# build logfile path
set logfile = "${outdir}/run_${cfg}_${stamp}.txt"
echo "[$me] logfile = $logfile"

echo "[$me] Test Debug = 5"

  set flags = ( -i "$inpath" -n "$NFILE" -t $THREADS -o "$outdir" -c "$cfg" )
  if ( $IS_MC )   set flags = ( $flags --mc )
  if ( $REPROC )  set flags = ( $flags --reproc )
  if ( $MINIMAL ) set flags = ( $flags --minimal )
 echo "[$me] Test Debug = 6"
  echo "------------------------------------------------------------------"
  echo " Config   : $cfg"
  echo " InputDir : $inpath"
  echo " OutDir   : $outdir"
  echo " NFile    : $NFILE"
  echo " Threads  : $THREADS"
  echo " Exec     : $EXECUTABLE $flags"
  echo " LogFile  : $logfile"
  echo "------------------------------------------------------------------"

  if ( $DRY ) continue

  ( $EXECUTABLE $flags |& tee "$logfile" )
  set status = $status
  echo "[$me] ${cfg} exit status: $status" |& tee -a "$logfile"
end

echo "All requested runs completed."

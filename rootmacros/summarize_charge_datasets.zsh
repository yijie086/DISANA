#!/usr/bin/env zsh
# ------------------------------------------------------------
# Batch runner for chargeSumRuns.groovy over fixed RGA datasets
# Directories copied from run_RGA_phi_analysis_final.zsh
module purge
# --- make sure 'module' command is available ---
if ! command -v module &>/dev/null; then
  # Lmod / environment-modules are usually initialized here on many clusters
  if [[ -r /etc/profile.d/modules.sh ]]; then
    source /etc/profile.d/modules.sh
  fi
fi
# ------------------------------------------------------------
# --- load CLAS12 environment ---
module load coatjava/13.4.0
module load groovy/4.0.20
module unload qadb/3.2.0  || true
module load qadb/3.4

set -euo pipefail 

# ---- Known datasets (same as phi analysis script) ----
CONFIGS=(rgasp18_inb rgasp18_outb rgafall18_inb rgafall18_outb rgasp19_inb)

# ---- Per-dataset input directories (copied verbatim) ----
IN_rgasp18_inb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus-1/pass1/dst/train/DVKpKmP/"
IN_rgasp18_outb="/cache/clas12/rg-a/production/recon/spring2018/10.59gev/torus+1/pass1/dst/train/DVKpKmP/"
IN_rgafall18_inb="/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/DVKpKmP/"
IN_rgafall18_outb="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/DVKpKmP/"
IN_rgasp19_inb="/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/DVKpKmP/"

me="$(basename "$0")"

die() {
  echo "ERROR: $@" >&2
  exit 1
}

# ---- usage ----
if (( $# > 1 )); then
  cat <<EOF
Usage:
  $me [output_summary.txt]

If no output file is given, defaults to: charge_summary_RGA.txt

This script:
  * runs chargeSumRuns.groovy on each fixed RGA dataset directory
  * parses the "total accumulated charge" line
  * writes a small summary text file
EOF
  exit 1
fi

# ---- locate chargeSumRuns.groovy ----
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CHARGE_SUM_GROOVY="${SCRIPT_DIR}/chargeSumRuns.groovy"

[[ -f "$CHARGE_SUM_GROOVY" ]] || die "chargeSumRuns.groovy not found at: $CHARGE_SUM_GROOVY"
[[ -n "${COATJAVA:-}" ]] || echo "WARNING: \$COATJAVA is not set; make sure your environment is initialized."

OUTFILE="${1:-charge_summary_RGA.txt}"

echo "[$me] Writing summary to: $OUTFILE"
: > "$OUTFILE"
echo "# dataset  total_charge_nC" >> "$OUTFILE"

for cfg in "${CONFIGS[@]}"; do
  echo
  echo "[$me] === Processing dataset '$cfg' ==="

  # IN variable name, exactly like in run_RGA_phi_analysis_final.zsh
  varname="IN_${cfg}"
  # ${(P)varname} is zsh indirection: value of variable whose name is in \$varname
  inpath="${(P)varname}"

  [[ -z "$inpath" ]] && { echo "[$me] WARNING: input path for $cfg not set"; continue; }
  [[ -d "$inpath" ]] || { echo "[$me] WARNING: '$inpath' is not a directory, skipping $cfg"; continue; }

  echo "[$me] Input dir: $inpath"

  tmpfile="$(mktemp)"
  # Add --recursive if your directory tree is nested with .hipo files
  "$COATJAVA/bin/run-groovy" "$CHARGE_SUM_GROOVY" "$inpath" --exclude 3262 --recursive > "$tmpfile"

  # Parse the final 'total accumulated charge' line
  ##total_charge="$(grep -i 'total accumulated charge:' "$tmpfile" | tail -n 1 | awk '{print $(NF-1)}')"
  total_charge="$(grep -i 'total no hel accumulated charge:' "$tmpfile" | tail -n 1 | awk '{print $(NF-1)}')"

  if [[ -z "$total_charge" ]]; then
    echo "[$me] WARNING: could not extract total charge for '$cfg'" >&2
  else
    echo "[$me] -> total charge = ${total_charge} nC"
    echo "$cfg  $total_charge" >> "$OUTFILE"
  fi

  rm -f "$tmpfile"
done

echo
echo "[$me] Done. Summary written to $OUTFILE"

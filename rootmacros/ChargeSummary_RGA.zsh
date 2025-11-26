#!/usr/bin/env zsh
# ------------------------------------------------------------
# Batch runner for qadb_charge_analysis.groovy over fixed RGA datasets.
#
# This script:
#   1. Runs 'qadb-info charge' to generate a run list for each dataset.
#   2. Runs qadb_charge_analysis.groovy in 'file-only' mode on the list, 
#      applying dataset-specific exclusions.
#   3. Captures the ENTIRE output (logs and totals) from the Groovy script 
#      and writes it to the output file.
#   4. Also writes a final summary table (dataset, total_charge).
# ------------------------------------------------------------
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

# ---- Known datasets ----
CONFIGS=(rgasp18_outb rgafall18_inb rgafall18_outb rgasp19_inb)
#CONFIGS=(rgasp19_inb)

# ---- Per-dataset QADB names and specific run exclusions ----
# The names below correspond to the QADB server entries.
QADB_NAMES=(rga_sp18_outbending rga_fa18_inbending rga_fa18_outbending rga_sp19)
#QADB_NAMES=(rga_sp19)

# Define exclusions for each dataset
# EXCLUSIONS should be in "run1,run2,runA-runB" format.
EXCLUSIONS_rgasp18_inb="3262" # Example: keep run 3262 exclusion from original script logic
EXCLUSIONS_rgasp18_outb=""
EXCLUSIONS_rgafall18_inb=""
EXCLUSIONS_rgafall18_outb=""
EXCLUSIONS_rgasp19_inb=""

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

This script runs the qadb_charge_analysis.groovy script for RGA datasets
and captures the full output, including per-run logs and final totals, 
into the summary file.
EOF
  exit 1
fi

# ---- locate the groovy script ----
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Renamed to match the latest file in the conversation history:
CHARGE_SUM_GROOVY="${SCRIPT_DIR}/chargeSumRuns.groovy" 

[[ -f "$CHARGE_SUM_GROOVY" ]] || die "chargeSumRuns.groovy not found at: $CHARGE_SUM_GROOVY"
[[ -n "${COATJAVA:-}" ]] || echo "WARNING: \$COATJAVA is not set; make sure your environment is initialized."

OUTFILE="${1:-charge_summary_RGA.txt}"

echo "[$me] Writing full logs and summary to: $OUTFILE"
: > "$OUTFILE"

# Initialize variables for the final summary section
SUMMARY_TABLE="# dataset  total_charge_nC\n"

i=1
for cfg in "${CONFIGS[@]}"; do
  qadb_name="${QADB_NAMES[$i]}"
  i=$((i+1))
  
  echo "\n[$me] === Processing dataset '$cfg' ($qadb_name) ===" | tee -a "$OUTFILE"
  echo "--- START LOG FOR $cfg ---" >> "$OUTFILE"

  # Step 1: Generate the run list file using qadb-info
  RUNLIST_FILE="$(mktemp --suffix=.txt)"
  echo "[$me] Generating run list with: qadb-info charge -d $qadb_name > $RUNLIST_FILE" | tee -a "$OUTFILE"
  # This command requires the 'qadb-info' tool to be in your PATH
  qadb-info charge -d "$qadb_name" > "$RUNLIST_FILE"

  if [[ ! -s "$RUNLIST_FILE" ]]; then
    echo "[$me] WARNING: Run list file '$RUNLIST_FILE' is empty or not generated for '$qadb_name', skipping." >&2 | tee -a "$OUTFILE"
    rm -f "$RUNLIST_FILE"
    echo "--- END LOG FOR $cfg ---" >> "$OUTFILE"
    continue
  fi
  
  # Step 2: Determine exclusion runs for this dataset
  exclusions_varname="EXCLUSIONS_${cfg}"
  # Use ${!exclusions_varname} expansion for POSIX compliance, although zsh supports (P)
  exclusions="${(P)exclusions_varname:-}" 
  
  # Prepare command arguments for the groovy script (file-only mode)
  cmd_args=("$CHARGE_SUM_GROOVY" "$RUNLIST_FILE")
  if [[ -n "$exclusions" ]]; then
    cmd_args+=(--exclude "$exclusions")
    echo "[$me] Applying exclusions: $exclusions" | tee -a "$OUTFILE"
  fi
  
  # Step 3: Run the modified Groovy script in file-only mode
  tmpfile="$(mktemp)"
  
  echo "[$me] Running: $COATJAVA/bin/run-groovy ${cmd_args[@]}" | tee -a "$OUTFILE"
  "$COATJAVA/bin/run-groovy" "${cmd_args[@]}" > "$tmpfile"

  # Write the entire Groovy log to the output file
  cat "$tmpfile" >> "$OUTFILE"
  echo "--- END LOG FOR $cfg ---" >> "$OUTFILE"

  # Parse the final 'total accumulated charge' line from the groovy script output
  total_charge="$(grep -i 'total accumulated charge:' "$tmpfile" | tail -n 1 | awk '{print $(NF-1)}')"

  if [[ -z "$total_charge" ]]; then
    echo "[$me] WARNING: could not extract total charge for '$cfg'" >&2
  else
    echo "[$me] -> total charge = ${total_charge} nC (Extracted for summary)"
    SUMMARY_TABLE+="$cfg  $total_charge\n"
  fi

  rm -f "$tmpfile" "$RUNLIST_FILE" # Clean up temporary files
done

# Print the final summary table
echo "\n[$me] === FINAL CHARGE SUMMARY TABLE ===" >> "$OUTFILE"
echo -e "$SUMMARY_TABLE" >> "$OUTFILE"

echo
echo "[$me] Done. Full logs and summary table written to $OUTFILE"
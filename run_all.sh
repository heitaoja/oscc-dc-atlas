#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- run the whole pipeline without `make`
# -----------------------------------------------------------------------------
# Usage:
#   ./run_all.sh              # run the discovery workflow in order
#   ./run_all.sh 12           # validate after downloading GSE181919
#   ./run_all.sh 05 06        # run only the steps whose name starts with 05/06
#   RSCRIPT=/usr/bin/Rscript ./run_all.sh
#
# Works in Git Bash / WSL / Linux. On Windows, if Rscript is not on PATH the
# script looks for the newest R installation under "C:/Program Files/R".
# =============================================================================
set -euo pipefail

STEPS=(
  01_load_qc
  02_integrate_cluster
  03_annotate_major
  04_dc_subset
  05_dc_annotate_de
  06_enrichment
  07_trajectory
  08_tf_activity
  09_composition
  10_qc_figure
  11_donor_pseudobulk
)
EXTERNAL_STEP=12_external_validation

# ---- locate Rscript --------------------------------------------------------
if [[ -z "${RSCRIPT:-}" ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    RSCRIPT="$(command -v Rscript)"
  else
    RSCRIPT="$(ls -d "/c/Program Files/R"/R-*/bin/Rscript.exe 2>/dev/null | sort -V | tail -1 || true)"
  fi
fi
if [[ -z "${RSCRIPT:-}" || ! -x "$RSCRIPT" ]]; then
  echo "ERROR: Rscript not found. Set RSCRIPT=/path/to/Rscript" >&2
  exit 1
fi
echo "Using R: $RSCRIPT"

mkdir -p results/logs results/objects results/figures results/tables

# ---- select steps ----------------------------------------------------------
if [[ $# -gt 0 ]]; then
  SELECTED=()
  for want in "$@"; do
    matched=false
    for s in "${STEPS[@]}" "$EXTERNAL_STEP"; do
      if [[ "$s" == "$want"* ]]; then
        SELECTED+=("$s")
        matched=true
      fi
    done
    if ! "$matched"; then
      echo "Unknown stage: $want" >&2
      exit 2
    fi
  done
else
  SELECTED=("${STEPS[@]}")
fi

# ---- run -------------------------------------------------------------------
start=$(date +%s)
for s in "${SELECTED[@]}"; do
  echo ""
  echo "=============================================================="
  echo "  $s    $(date '+%H:%M:%S')"
  echo "=============================================================="
  "$RSCRIPT" "R/$s.R" 2>&1 | tee "results/logs/$s.log"
done
end=$(date +%s)
echo ""
echo "All selected steps finished in $(( (end-start)/60 )) min $(( (end-start)%60 )) s."

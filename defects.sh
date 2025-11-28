#!/usr/bin/env bash
# defects.sh — run DefectSweepFixedN and plot with gnuplot
#
# Usage:
#   ./defects.sh a b c Nions sweeps lagtime [phi_min] [phi_step] [phi_max] [outfile] [outdir]
#
# Example:
#   ./defects.sh 30 10 10 150 50000 5 0.0 0.05 0.60 data/defect_sweep.tsv plots

set -euo pipefail

if [[ $# -lt 6 ]]; then
  echo "Usage: $0 a b c Nions sweeps lagtime [phi_min] [phi_step] [phi_max] [outfile] [outdir]"
  exit 1
fi

a="$1"; b="$2"; c="$3"; Nions="$4"; sweeps="$5"; lag="$6"
phi_min="${7:-0.0}"
phi_step="${8:-0.05}"
phi_max="${9:-0.95}"
outfile="${10:-data/defect_sweep.tsv}"
outdir="${11:-plots}"

# Paths assuming this script is at repo root
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
RUN_JL="${ROOT_DIR}/run_defects.jl"
GP_SCRIPT="${ROOT_DIR}/gnuplot_scripts/defect_sweep.gp"

mkdir -p "$(dirname "$outfile")" "$outdir"

# 1) Run simulation/export using repo project
julia --project="$ROOT_DIR" "$RUN_JL" "$a" "$b" "$c" "$Nions" "$sweeps" "$lag" \
      "$phi_min" "$phi_step" "$phi_max" "$outfile"

# 2) Plot (expects 'infile' + 'outdir')
if [[ -f "$GP_SCRIPT" ]]; then
  gnuplot -e "infile='${outfile}'; outdir='${outdir}'" "$GP_SCRIPT"
  echo "Plots written to: $outdir/"
else
  echo "Note: gnuplot script not found at $GP_SCRIPT. Skipping plots."
fi

echo "TSV written to: $outfile"

#!/usr/bin/env bash
# morgan.sh
# Run MorganSweep and produce plots
#
# Usage:
#   ./morgan.sh a b c sweeps lagtime [a_lat] [kB] [T] [q] [outfile] [outdir]
#
# Example:
#   ./morgan.sh 30 30 30 50000 20 1.0 1.0 1.0 1.0 data/Morgan_noninteracting.tsv plots

set -euo pipefail

if [ $# -lt 5 ]; then
  echo "Usage: $0 a b c sweeps lagtime [a_lat] [kB] [T] [q] [outfile] [outdir]"
  exit 1
fi

a="$1"; b="$2"; c="$3"; sweeps="$4"; lag="$5"
a_lat="${6:-1.0}"
kB="${7:-1.0}"
T="${8:-1.0}"
q="${9:-1.0}"
outfile="${10:-data/Morgan_noninteracting.tsv}"
outdir="${11:-plots}"

mkdir -p "$(dirname "$outfile")" "$outdir"

# 1) Run simulation/export
./run_morgan.jl "$a" "$b" "$c" "$sweeps" "$lag" "$a_lat" "$kB" "$T" "$q" "$outfile"

# 2) Plot
gnuplot -e "infile='${outfile}'; outdir='${outdir}'" gnuplot_scripts/plot_morgan.gp

echo "All done."
echo "TSV  : $outfile"
echo "Plots: $outdir/"

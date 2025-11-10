# gnuplot_scripts/plot_morgan.gnuplot_scripts
# Usage:
# gnuplot -e "infile='data/Morgan_noninteracting.tsv'; outdir='plots/'" gnuplot_scripts/plot_morgan.gnuplot

if(!exists("infile")) infile = "data/Morgan_noninteracting.tsv"
if(!exists("outdir")) outdir = "plots"

set datafile separator '\t'
set key autotitle columnhead
set grid
set term pngcairo size 1400, 900 font ",14"

# 1) Dtr vs Dbulk
set output sprintf("%s/Dtr_Dbulk.png", outdir)
set xlabel "Occupancy percentile"
set ylabel "Diffusion coefficient (arb. units)"
plot infile using 1:2 with linespoints lw 2 pt 7 t "D_{tr}", \
     ''     using 1:3 with linespoints lw 2 pt 5 t "D_{bulk}"

# 2) Haven ratio
set output sprintf("%s/Haven.png", outdir)
set xlabel "Occupancy percentile"
set ylabel "Haven ratio (Dtr / Dbulk)"
plot infile using 1:4 with linespoints lw 2 pt 7 t "Haven"

# 3) Correlation factors f_tr and f_col
set output sprintf("%s/CorrelationFactors.png", outdir)
set xlabel "Occupancy percentile"
set ylabel "Correlation factor"
plot infile using 1:5 with linespoints lw 2 pt 7 t "f_{tr}", \
     ''     using 1:6 with linespoints lw 2 pt 5 t "f_{col}"

# 4) Reduced conductivity
set output sprintf("%s/ReducedConductivity.png", outdir)
set xlabel "Occupancy percentile"
set ylabel "Reduced conductivity (C q^2 D_{tr} / kBT)"
plot infile using 1:7 with linespoints lw 2 pt 7 t "σ_{red}"
set datafile separator '\t'           # change to ',' if CSV
set key autotitle columnhead
set term pngcairo size 1200,800
set output 'D_vs_loading.png'
set grid
set xlabel 'Occupancy percentile'
set ylabel 'Diffusion coefficient'
#set logscale y                       # optional if values span orders of magnitude
set xrange [0:1]

plot 'haven_sweep.tsv' using 1:2 with linespoints t 'D_tr', \
     ''                using 1:3 with linespoints t 'D_σ'

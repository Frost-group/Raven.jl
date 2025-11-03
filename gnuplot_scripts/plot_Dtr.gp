set datafile separator '\t'
set key autotitle columnhead
set xlabel 'Percentile'
set ylabel 'D_tr'
set grid
set term pngcairo size 1200,800
set output 'percentile_vs_Dtr.png'
plot 'Dtr_sweep.tsv' using 1:2 with linespoints


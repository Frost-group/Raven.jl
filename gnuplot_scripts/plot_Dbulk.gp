set datafile separator '\t'
set key autotitle columnhead
set xlabel 'Percentile'
set ylabel 'D_bulk'
set grid
set term pngcairo size 1200,800
set output 'percentile_vs_Dbulk.png'
plot 'Dbulk_sweep.tsv' using 1:2 with linespoints
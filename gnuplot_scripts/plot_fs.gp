set datafile separator '\t'           # change to ',' if CSV
set key autotitle columnhead
set term pngcairo size 1200,800
set output 'f_vsfI.png'
set grid
set xlabel 'Occupancy percentile'
set ylabel 'Correlation factors'
#set logscale y                       # optional if values span orders of magnitude
set xrange [0:0.95]

plot 'Morgan_noninteracting.tsv' using 1:5 with linespoints t 'f_{tr}', \
     ''                using 1:6 with linespoints t 'f_{col}'
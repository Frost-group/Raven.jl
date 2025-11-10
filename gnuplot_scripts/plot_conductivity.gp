set datafile separator '\t'           # change to ',' if CSV
set key autotitle columnhead
set term pngcairo size 1200,800
set output 'Sigma_vs_Occupancy.png'
set grid
set xlabel 'Occupancy percentile'
set ylabel 'Reduced conductivity σ*'
set xrange [0:1]
set yrange [0:0.04]

plot 'Morgan_noninteracting.tsv' using 1:7 with linespoints t 'conductivity'
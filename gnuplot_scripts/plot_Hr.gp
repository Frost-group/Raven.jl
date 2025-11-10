set datafile separator '\t'           # change to ',' if CSV
set key autotitle columnhead
set term pngcairo size 1200,800
set output 'Haven_vs_Occupancy.png'
set grid
set xlabel 'Occupancy percentile'
set ylabel 'Haven ratio (D_σ / D_tr)'
set xrange [0:0.95]
set yrange [0:1.5]

plot 'Morgan_noninteracting.tsv' using 1:4 with linespoints t 'Haven Ratio'

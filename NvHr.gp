# NvHr.gp
file = "data/S_scan_N_vs_Hr.tsv"
N = 21   # number of xi blocks (e.g. xi = 0..10)

set datafile separator "\t"
set key outside
set xlabel "N"
set ylabel "Hr"

#set logscale y
set yrange [0:1.3]

set term pngcairo size 900,600
set output "NvsHr_disorder_scan.png"

plot for [i=0:N-1] file index i using 1:4 with lines title sprintf("S=%.2f", i*4.0)
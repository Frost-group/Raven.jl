# NvHr.gp
file = "data/NvsHr_test.tsv"
N = 21   # number of xi blocks (e.g. xi = 0..10)

set datafile separator "\t"
set key outside
set xlabel "N"
set ylabel "Hr"

#set logscale y
set yrange [0.8:1.3]

set term pngcairo size 900,600
set output "NvsHr_scan.png"

plot for [i=0:N-1] file index i using 1:5 with lines title sprintf("beta=%.2f", i*0.05)
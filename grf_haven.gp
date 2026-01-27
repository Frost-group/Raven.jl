# grf_haven.gp
file = "data/DtrDbulk_test.tsv"
N = 11   # number of xi blocks (e.g. xi = 0..10)

set datafile separator "\t"
set key outside
set xlabel "S"
set ylabel "Hr"

#set logscale y
set yrange [0:1.3]

set term pngcairo size 900,600
set output "Hr_scan.png"

#f(x) = exp(-x)

plot for [i=0:N-1] file index i using 2:5 with lines title sprintf("xi=%g", i)

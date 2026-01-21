# xi_scan1.gp
file = "data/xi_scan1.tsv"
N = 11   # number of xi blocks (e.g. xi = 0..10)

set datafile separator "\t"
set key outside
set xlabel "S"
set ylabel "D/D0"

#set logscale y

set term pngcairo size 900,600
set output "xi_scan_raw.png"

f(x) = exp(-x)

plot f(x) w l lw 2 title "exp(-S)", \
    for [i=0:N-1] file index i using 2:4 with lines title sprintf("xi=%g", i)

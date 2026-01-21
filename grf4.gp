set datafile separator "\t"
set datafile missing "NaN"
set grid
set key outside
set xlabel "S"
set ylabel "D/D0"

set term pngcairo size 1000,650 enhanced
set output "SvsD_overlay_raw.png"

set xrange [0:20]
#set logscale y
f(x) = exp(-x)

plot f(x) w l lw 2 title "exp(-S)", \
     for [i=0:7] "data/attempt2.tsv" index i every ::1 using 2:3 w lp title sprintf("β block %d", i)

unset output

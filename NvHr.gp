# NvHr.gp
file = "data/particle_scan_prod_run2.tsv"
N = 6   # number of blocks

set datafile separator "\t"
set key outside
set xlabel "N"
set ylabel "Hr"

#set logscale y
#set yrange [0:0.2]

set term pngcairo size 900,600
set output "NvsHr_over_Beta2.png"

plot for [i=0:N-1] file index i using 1:5 with lines title sprintf("β=%.2f", i*4.0)
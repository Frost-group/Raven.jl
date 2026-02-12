# NvHr.gp
file = "data/particle_scan_prod_run3.tsv"
N = 6   # number of blocks

set datafile separator "\t"
set key outside
set xlabel "N"
set ylabel "Db"

#set logscale y
#set yrange [0:0.2]

set term pngcairo size 900,600
set output "NvsDb_over_Beta3.png"

plot for [i=0:N-1] file index i using 1:4 with points title sprintf("β=%.2f", i*4.0)
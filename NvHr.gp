# NvHr.gp
file = "data/beta_vs_sweep_msd_N400.tsv"
N = 5   # number of blocks

set datafile separator "\t"
set key outside
set xlabel "Sweeps"
set ylabel "Tracer MSD"

#set logscale y
#set yrange [0:0.2]

set term pngcairo size 900,600
set output "sweeps_vs_trMSD.png"

plot for [i=0:N-1] file index i using 2:3 with points title sprintf("β=%.2f", (i-1)*0.25)
# NvHr.gp
file  = "data/beta_vs_sweep_msd_N400.tsv"
Temps = "500"   # IMPORTANT: spaces, no commas

set datafile separator "\t"
set datafile missing "NaN"
set key outside
set xlabel "Sweeps"
set ylabel "tracer MSD"
set xrange [0:*]

set term pngcairo size 900,600
set output "sweeps_vs_trMSD-3.png"

n = words(Temps)

# gnuplot block indices start at 0
plot for [i=0:n-1] file index i using 2:3 with lines \
    title sprintf("T=%s", word(Temps, i+1))

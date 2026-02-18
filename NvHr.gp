# NvHr.gp
file  = "data/T0.5-25K_scan_sigma3.0.tsv"
Temps = "0.5 1 2 5 10 15 20 25"   # IMPORTANT: spaces, no commas

set datafile separator "\t"
set datafile missing "NaN"
set key outside
set xlabel "Sweeps"
set ylabel "tracer MSD"
set xrange [0:*]

set term pngcairo size 900,600
set output "T0.5-25K_scan_sigma3.0.png"

n = words(Temps)

# gnuplot block indices start at 0
plot for [i=0:n-1] file index i using 2:3 with lines \
    title sprintf("T=%s", word(Temps, i+1))

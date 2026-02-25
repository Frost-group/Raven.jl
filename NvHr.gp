# NvHr.gp
file  = "data/T25-500K_scan_sigma3.0-upscaledV.tsv"
Temps = "500 400 200 100 50 25"   # IMPORTANT: spaces, no commas

set datafile separator "\t"
set datafile missing "NaN"
set key outside
set xlabel "Sweeps"
set ylabel "tracer MSD"
set xrange [0:*]

set term pngcairo size 900,600
set output "T25-500K_scan_sigma3.0-upscaledV.png"

n = words(Temps)

# gnuplot block indices start at 0
plot for [i=0:n-1] file index i using 2:3 with lines \
    title sprintf("T=%s", word(Temps, i+1))

set datafile separator "\t"
set grid
set key left top
set xlabel "S = β^2 χ0 / 3"
set ylabel "Diffusion coefficient D"
set title "GRF disorder vs diffusion (ξ=2, β=1)"
set logscale x
# set logscale y   # uncomment if you want log D too

set term pngcairo size 900,600
set output "grf_S_vs_D.png"
plot "data/grf_disorder_vs_diffusion.tsv" using 5:6 every ::1 with linespoints title "S vs D"

set output
set term qt
replot
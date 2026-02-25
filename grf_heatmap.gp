set terminal pngcairo size 900,700 enhanced font 'Arial,11'
set output 'grf_heatmap.png'

set title 'GRF slice'
set xlabel 'x'
set ylabel 'y'
set view map
set size ratio -1
set pm3d map
unset key

splot 'data/grf_slice_zmid.dat' using 1:2:3 with pm3d

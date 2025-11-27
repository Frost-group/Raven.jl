set term pngcairo size 1200,900
set datafile sep "\t"
set view map
set size ratio -1
set key off
set xlabel "x"; set ylabel "y"

set output "plots/pes_heatmap.png"
plot "data/pes_slice_z5.tsv" using 1:2:3 with image

set dgrid3d 60, 60 splines
set output "plots/pes_contour.png"
set contour base; unset surface
splot "data/pes_slice_z5.tsv" using 1:2:3 with lines

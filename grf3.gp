# grf_DoverD0_vs_S_filterbeta.gp

set datafile separator "\t"
set datafile missing "NaN"
set grid
set key left top

#set logscale x
set logscale y
set xlabel "S = β^2 χ0 / 3"
set ylabel "D / D0"
set title "GRF: normalized diffusion vs disorder"
set xrange [0:5]
set yrange [0:1.1]

file  = "data/ugrid.tsv"
betas = "0.0005 0.001 0.025 0.05 0.1 0.25 0.5 1.0"

eps = 1e-6
isbeta(x,B) = (abs(x-(B+0)) < eps)

set term pngcairo size 1000,650 enhanced
set output "grf_DoverD0_vs_S_quench.png"

plot for [B in betas] \
  file every ::1 using \
    ((isbeta($3,B) && $6>=0 && $9>0) ? $6 : 1/0) : \
    ((isbeta($3,B) && $6>=0 && $9>0) ? $9 : 1/0) \
  with linespoints title sprintf("β=%g (T=%g)", (B+0), 1.0/(B+0))

unset output

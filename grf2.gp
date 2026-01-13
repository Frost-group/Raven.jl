# grf_DoverD0_vs_S_filterbeta.gp
set datafile separator "\t"
set datafile missing "NaN"
set grid
set key left top
#set logscale x
set xlabel "S = β^2 χ0 / 3"
set ylabel "D / D0"
set title "GRF: normalized diffusion vs disorder (ξ=2)"
set xrange [0:20]
set yrange [0:1.05]

file  = "data/grf_disorder_vs_diffusion.tsv"
betas = "0.004 0.01 0.02 0.04 0.1 0.5 1.0"

set term pngcairo size 1000,650 enhanced
set output "grf_DoverD0_vs_S.png"

plot for [B in betas] \
  file every ::1 using \
    ( (abs($3-(B+0))<1e-12 && $6>0) ? $6 : 1/0 ) : \
    ( (abs($3-(B+0))<1e-12 && $6>0) ? $9 : 1/0 ) \
  with linespoints title sprintf("β=%s (T=%.5g)", B, 1.0/(B+0))

unset output

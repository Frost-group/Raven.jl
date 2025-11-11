# plot_defects.gp
# expects: -e "infile='data/defect_sweep.tsv'; outdir='plots'"

set datafile separator "\t"
set term pngcairo size 1200,800
set grid back
set key top right
set xlabel "Defect fraction (φ)"
set ylabel "D(φ) / D(0)"
set yrange [0:1.02]

set output sprintf("%s/defects_normalized.png", outdir)
plot infile using 1:9  with linespoints pt 7 lw 2 title "D_tr(φ)/D_tr(0)", \
     infile using 1:10 with linespoints pt 5 lw 2 title "D_bulk(φ)/D_bulk(0)"

# Optional: raw Dtr/Dbulk vs φ
set ylabel "Diffusion coefficient"
set output sprintf("%s/defects_raw.png", outdir)
plot infile using 1:3 with linespoints lw 2 title "D_tr(φ)", \
     infile using 1:4 with linespoints lw 2 title "D_bulk(φ)"

# gnuplot_scripts/defect_sweep.gp
set datafile separator whitespace
set key right top
set grid
set term pngcairo size 1200,800 font ",12"

# Point gnuplot at your TSV
file = "data/defect_sweep.tsv"

# Tell gnuplot the first row is headers (so NaN stays NaN)
set datafile commentschars "#"
set key autotitle columnhead

# 1:defect_frac, 2:M, 3:Dtr, 4:Dbulk, 5:Haven, 6:f_tr, 7:f_col, 8:reduced_conductivity

# 1) Dtr and Dbulk vs defect_frac
set output "plots/defects_Dtr_Dbulk.png"
set xlabel "Defect fraction (ϕ)"
set ylabel "Diffusion coefficient (arb. units)"
plot file using 1:3 with linespoints lw 2 pt 7 title "D_{tr}", \
     file using 1:4 with linespoints lw 2 pt 5 title "D_{bulk}"

# 2) Haven ratio vs defect_frac
set output "plots/defects_Haven.png"
set ylabel "Haven ratio (D_{tr} / D_{bulk})"
set yrange [0:1.5]
plot file using 1:5 with linespoints lw 2 pt 7 title "Haven"

# 3) Correlation factors vs defect_frac
set output "plots/defects_correlation_factors.png"
set ylabel "Correlation factor f"
set yrange [0:*]
plot file using 1:6 with linespoints lw 2 pt 7 title "f_{tr}", \
     file using 1:7 with linespoints lw 2 pt 5 title "f_{col}"

# 4) Reduced conductivity vs defect_frac
set output "plots/defects_sigma_red.png"
set ylabel "Reduced conductivity (C q^2 D_{tr} / k_B T)"
set yrange [0:*]
plot file using 1:8 with linespoints lw 2 pt 7 title "σ_{red}"

# 5) Normalized: Dtr_norm (col 9), Dbulk_norm (col 10)
set output "plots/defects_normalized.png"
set xlabel "Defect fraction (ϕ)"
set ylabel "D(ϕ) / D(0)"
set yrange [0:*]
plot file using 1:9  with linespoints lw 2 pt 7 title "D_{tr} (ϕ)/D_{tr} (0)", \
     file using 1:10 with linespoints lw 2 pt 5 title "D_{bulk} (ϕ)/D_{bulk} (0)"